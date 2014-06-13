;;
;; ParaGauss, a program package for high-performance computations
;; of molecular systems
;; Copyright (C) 2014
;; T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
;; M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
;; A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
;; T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
;; M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
;; M. Roderus, N. Rösch
;;
;; This program is free software; you can redistribute it and/or modify it
;; under the terms of the GNU General Public License version 2 as published
;; by the Free Software Foundation [1].
;;
;; This program is distributed in the hope that it will be useful, but
;; WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
;; General Public License for more details.
;;
;; [1] http://www.gnu.org/licenses/gpl-2.0.html
;;
;; Please see the accompanying LICENSE file for further information.
;;
;;;
;;; Do not confuse with (baslib guile paragauss) in
;;;
;;; ./baslib/guile/paragauss.scm
;;;
;;; Some are home-grown  modules.  This needs %load-path to  be set up
;;; properly.
;;;
(define-module (guile paragauss)
  #:use-module (srfi srfi-1)            ; list functions
  #:use-module (ice-9 pretty-print)
  #:use-module (ice-9 getopt-long)
  ;;
  ;; Custom modules, need properly set %load-path:
  ;;
  #:use-module ((guile comm)
                #:select (comm-size
                          comm-rank
                          comm-bcast
                          comm-barrier
                          critical))
  #:use-module ((guile utils)
                #:select (begin0
                          let-env
                          derivatives
                          list-derivatives))
  ;;
  ;; Import basis set library:
  ;;
;;#:use-module ((baslib guile baslib)
;;              #:select (find-basis))
;;#:use-module ((baslib guile paragauss)
;;              #:select (qm-write-input))
  #:export (qm-main
            ;; qm-init
            ;; qm-run
            ;; qm-finalize
            qm-trace-hook
            qm-flush-trace
            ;; qm-find-basis
            qm-vdw-dft-phi
            qm-vdw-dft-phi/asy
            pople-radius
            slater-radius
            ionic-radius
            with-stdout-to-file
            make-prl-plot
            test-xc))

;;;
;;; This name has to be defined on guile startup, see the C sources of
;;; guile_main () in ./guile-qm.c:
;;;
(define guile-paragauss-module-init
  (@@ (guile-user) guile-paragauss-module-init))


;;;
;;; The Fortran  function (see ../modules/paragauss.f90)  decorated by
;;; this C-wrapper  (./guile-qm.c) defines  some of the  auxilary qm-*
;;; symbols.  FIXME: you better do not call any of those whose Fortran
;;; sources were  compiled with Intel.  Not until Intel learns  how to
;;; return a struct properly (as with GCC -freg-struct-return).
;;;
;;; The list of the procedures defined by this call includes:
;;;
;;;   qm-init
;;;   qm-run
;;;   qm-finalize
;;;
;;; and posissibly more, depending on the compilation options.
;;;
(guile-paragauss-module-init)

;;;
;;; Beware that writing to the  same file from multiple workers is not
;;; going to end good:
;;;
(define (with-stdout-to-file file proc)
  (with-output-to-file file
    (lambda ()
      (let ((dup-1 (dup 1))         ; make a dup of the current stdout
            (fd (fileno (current-output-port))))
        (dup2 fd 1)     ; close stdout and redirect it to fd
        (proc)          ; execute thunk with stdout redirected to file
        (dup2 dup-1 1)  ; close file, stdout to original destination
        (close dup-1)))))               ; dont leak file descriptors

;; (with-stdout-to-file "a1"
;;   (lambda ()
;;     (display "Hi!")
;;     (newline)))

;;;
;;; This is artificial intelligence guessing temp dir:
;;;
(define (guess-temp-dir input)
  (let ((prefix (or (getenv "SCRATCH") ; set on SuperMUC
                    (getenv "OPT_TMP") ; used by LRZ
                    (getenv "TEMP")    ; not sure if ever used
                    (and (file-exists? "/scratch") "/scratch")
                    "/tmp")))           ; if none is set use this
    (string-append prefix "/" (getenv "USER") "-"
                   (input-base-name input) "-"
                   (number->string (getpid)))))

;;;
;;; This emulates the behaviour of  the "runpg" bash script that tells
;;; PG to put output into o.name/ for an input i.name:
;;;
(define (guess-output-dir input)
  (cond
    ((string=? "input" input) ".") ; for input named input use CWD
    (else                          ; by default use o.$name for output
     (string-append "o." (input-base-name input)))))

;;;
;;; Normalize input name by dropping inessential suffix or prefix:
;;;
(define (input-base-name input)
  (cond
    ((string-suffix? ".scm" input) (string-drop-right input 4))
    ((string-suffix? ".nml" input) (string-drop-right input 4))
    ((string-prefix? "i." input) (string-drop input 2))
    (else input)))

;;;
;;; Now that we are  responsible for creating directories ourselves we
;;; need to coordinate the work between process within and between SMP
;;; hosts:
;;;
(define (maybe-mkdir! world dirname) ; dirname is the same for all ranks
  (critical world
            (lambda ()
              (if (not (file-exists? dirname))
                  (mkdir dirname)))))

(define (maybe-rm-rf! world dirname) ; dirname is the same for all ranks
  (critical world
            (lambda ()
              (if (file-exists? dirname)
                  (system* "rm" "-rf" dirname))))) ; DANGEROUS !!!

;;;
;;; Call  proc  without  args  with  a  TTFSTMP  environment  variable
;;; pointing to a fresh temporary directory on each rank in the world.
;;; Set TTFSTMP  to something different  from $PWD to  avoid polluting
;;; working directory with temporary files before running PG:
;;;
(define (with-temp-dir world dirname proc)
  (maybe-mkdir! world dirname)          ; make dir on each rank
  (let-env (("TTFSTMP" dirname))        ; temporarily change env
    (begin0                             ; returns the first expr
      (proc)                            ; <- value of begin0 form
      (maybe-rm-rf! world dirname))))   ; remove dir, DANGEROUS !!!


;;;
;;; This runs an input stored in a namelist format in a file. Input is
;;; the name of that file:
;;;
(define (run-path-nml world input)
  (let ((temp-dir
         (comm-bcast world 0 (guess-temp-dir input))) ; prefer the value at rank 0
        (output-dir
         (comm-bcast world 0 (guess-output-dir input))))
    ;;
    (maybe-mkdir! world output-dir)
    ;;
    ;; Set the environment for the time PG is running, then restore:
    ;;
    (let-env (("TTFSINPUT" input)
              ("TTFSOUTPUTDIR" output-dir))
      ;;
      ;; This creates temp-dir, runs the thunk and removes that dir:
      ;;
      (with-temp-dir world
                     temp-dir
                     (lambda () (qm-run world)))) ; this invokes the program
    ;;
    ;; Dump trace log, if not empty, into a file:
    ;;
    (qm-flush-trace world input output-dir)
    ;;
    ;; Return total enery, set in guile-user (global) namespace by
    ;; qm_run (), see modules/paragauss.f90:
    ;;
    (@@ (guile-user) *energy*)))

;;;
;;; This runs  an input stored in  sexp format in a  file. path-sxp is
;;; the name of that file.
;;;
;;; NOTE: for a path-sexp such as base.scm or i.base the file base.nml
;;; will be overwritten!
;;;
(define (run-path-sxp world path-sxp)
  (let ((path-nml (string-append (input-base-name path-sxp) ".nml"))
        (input (with-input-from-file path-sxp read)))
    ;;
    ;; Hm, here is a place where multiple workers may attempt to write
    ;; to the same file. We  anyway assume that the input directory is
    ;; accessible by any worker, so let only the master do that:
    ;;
    (if (equal? (comm-rank world) 0)
        (with-output-to-file path-nml
          (lambda () (qm-write-input input))))
    ;;
    ;; Other  workers should not  attempt processing the  input before
    ;; master has finished writing. Cross the barrer together.
    ;;
    (comm-barrier world)
    (run-path-nml world path-nml)))

;;;
;;; Invoke  (program world)  with the  world communicator  privided by
;;; qm-init:
;;;
(define (call-with-qm-world program)
  (let ((world (qm-init)))           ; intialize MPI, get communicator
    (begin0                          ; returns first expr
      (program world)                ; run the program, note result
      (qm-finalize world))))         ; finalize MPI


(define (qm-main argv)
  ;;
  ;; No support for command line options yet. First element of argv is
  ;; ignored:
  ;;
  (let* ((args (cdr argv))
         ;;
         ;; A path ending in .scm is assumed to contain the input as
         ;; an s-expression.  Otherwise it is assumed to contain a
         ;; plain-old namelist input:
         ;;
         (run-path
          (lambda (world path)
            (if (string-suffix? ".scm" path)
                (run-path-sxp world path)
                (run-path-nml world path))))
         ;;
         ;; A program that  takes a communicator, processes all inputs
         ;; each specified by a path and returns a list of results:
         ;;
         (program
          (lambda (world)
            (map (lambda (path)
                   (run-path world path))
                 args))))
    ;;
    ;; Call the program with an MPI comm and return result:
    ;;
    (call-with-qm-world program)))

;;;
;;; This code  is for  testing XC functionals.   Note that  qm-xc only
;;; becomes defined after executing guile-paragauss-module-init:
;;;
(define (write-row row)
  (display (string-join (map number->string row) "\t"))
  (newline))

(define (transpose x)
  "(transpose lol) -> lol

Transpose a list of lists."
  (apply map list x))

(define (unpack-upper lst)
  (let loop ((n 0)
             (lst (reverse lst))
             (acc '()))
    (if (null? lst)
        (map reverse acc)
        (let* ((n (+ n 1))
               (row (list-head lst n))
               (lst (list-tail lst n)))
          (loop n
                lst
                (cons row acc))))))

(define (get-upper lol)
  (let loop ((n 0)
             (lol lol)
             (acc '()))
    (if (null? lol)
        (reverse acc)
        (let ((row (list-tail (car lol) n)))
          (loop (+ n 1)
                (cdr lol)
                (cons row acc))))))

(define (test-xc xc)
  (let* ((rhos (list 1.0 2.0 3.0 4.0 5.0)) ; ra rb gaa gbb gab
         (h 0.1)         ; initial scale for numerical differentiation
         (fxc (lambda args (first (apply qm-xc xc args)))) ; a scalar function
         (vxc (lambda args (second (apply qm-xc xc args)))) ; a list-valued function
         (kxc (lambda args (third (apply qm-xc xc args))))) ; a "triangular matrix" as a list
    (format #t "Value:\n")
    (pretty-print (apply fxc rhos))
    (newline)
    (format #t "First derivatives (analytical, then numerical):\n")
    (write-row (apply vxc rhos))
    (write-row (derivatives fxc rhos h))
    (newline)
    (format #t "Second derivatives (analytical, then numerical upper, then lower):\n")
    (for-each write-row (unpack-upper (apply kxc rhos)))
    (let ((numerical (list-derivatives vxc rhos h)))
      (for-each write-row (get-upper numerical))
      (for-each write-row (get-upper (transpose numerical))))))

;; (test-xc "xalpha")
;; (test-xc "vwn")
;; (test-xc "pwlda")
;; (test-xc "becke88")
;; (test-xc "becke88-old")
;; (test-xc "pbec")
;; (test-xc "pbecsol")

(define (make-prl-plot n)
  (let* ((pi (* 4 (atan 1)))       ; is there a PI constant somewhere?
         (phi (lambda (D d)        ; 4 * pi * D**2 * phi(D, d)
                (* 4 pi (* D D) (qm-vdw-dft-phi (* D (+ 1 d))
                                                (* D (- 1 d))))))
         (ds (list 0.0 0.5 1.0))        ; three values for deltas
         (Ds (iota n 0.0 (/ 8.0 n))))   ; n values for Ds
    ;;
    ;; Return list  of lists with rows corresponding to  n values of D
    ;; and  columns corresponding to (three) values  of delta. Prepend
    ;; the value of D to each row for convenience:
    ;;
    (map (lambda (D)
           (cons D (map (lambda (d) (phi D d)) ds)))
         Ds)))

;; (qm-vdw-dft-phi 1.0 1.0)
;; 0.117493999122921
;; (qm-vdw-dft-phi 0.0 0.0)
;; 2.21067277513615
;; (qm-vdw-dft-phi 1.0 0.0)
;; 0.207226645435212

(define (qm-vdw-dft-phi/asy d1 d2)
  "(qm-vdw-dft-phi/asy d1 d2) -> double

Asymptotic form of qm-vdw-dft-phi for large d1 and d2."
  (let* ((%pi (* 4 (atan 1)))
         (C (* 12 (expt (* 4/9 %pi) 3)))
         (dd1 (* d1 d1))
         (dd2 (* d2 d2)))
    (/ (- C)
       (* dd1 dd2 (+ dd1 dd2)))))


;; (for-each write-row (make-prl-plot 10))


;;;
;;; This function  is called  from PG in  some cases (notably  when no
;;; basis is supplied in the input)  and is supposed to return a basis
;;; as  specified  by  atom   and  the  type.   The  specification  is
;;; intentionally left vague at the moment.
;;;
;;(define (qm-find-basis atom type)
;;  (let* ((atom (string->symbol (string-downcase (string-trim-both atom))))
;;         (lib (list 'turbomole 'basen atom))
;;         (basis (find-basis lib type)))
;;    (if #f                              ; debug prints
;;        (begin
;;          (display (list "QM-FIND-BASIS:" atom type))
;;          (newline)
;;          (pretty-print basis)
;;          (newline)
;;          (force-output)))              ; flush stdout
;;    basis))                             ; return basis

;;;
;;; One of  these, qm-trace-hook, if defined/imported,  is called from
;;; scheme_trace_hook() in PG.   So uncomment it if you  want a trace.
;;; The counterpart,  qm-flush-trace, will do nothing if  the trace is
;;; empty.
;;;
;;; In  debug runs  here we  collect a  list of  tracing  entries, see
;;; qm-trace-hook, qm-flush-trace:
;;;
(define *trace* '())

;;;
;;; This one is eventually called from PG if the TRACE() macro exapnds
;;; to something proper:
;;;
(define (qm-trace-hook key file line time)
  "In some debug modes this funciton is called, if defined."
  ;; (display (list key file line time)) (newline) (force-output)
  ;;
  ;; Cons a new entry onto the global list:
  ;;
  (set! *trace* (cons (list key file line time) *trace*))
  0)                                    ; return something

;;;
;;; Dump trace entries collected during the run in *trace*:
;;;
(define (qm-flush-trace world input output-dir)
  (if (not (null? *trace*))             ; dont create empty files
      (let* ((size (comm-size world))
             (rank (comm-rank world))
             (path (string-append output-dir "/"
                                  input "-"
                                  (number->string size) "-"
                                  (number->string rank))))
        ;;
        ;; Write the list of trace entries into a file:
        ;;
        (with-output-to-file path
          (lambda () (pretty-print (reverse *trace*))))
        ;;
        ;; Empty list of trace entries:
        ;;
        (set! *trace* '()))))
