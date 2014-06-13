;;;
;;; Copyright (c) 2011-2013 Alexei Matveev
;;;
(define-module (baslib guile baslib)
  #:export (find-basis
            find-ecp
            read-library))

(use-modules (srfi srfi-1)              ; list processing functions
             (srfi srfi-2)              ; and-let*
             (srfi srfi-11)             ; let-values
             (ice-9 match)              ; pattern matching
             (ice-9 pretty-print)
             ((guile utils)
              #:select (memoize))
             ;; Need to augment %load-path for this to work:
             ((baslib guile turbomole)
              #:select (turbomole-read))
             ((baslib guile nwchem)
              #:select (nwchem-read
                        make-basis      ; to create generalized ones
                        basis?
                        basis-name
                        basis-segments
                        ecp?
                        ecp-name
                        ecp-core
                        ecp-shells)))
;;
;; FIXME: description needs update. NWChem parser returns basis- and
;; ecp records.
;;
;; Basis data as returned by  turbomole-read is a list of entries that
;; may look  like this  (the car position  may eventually be  a longer
;; list of aliases). If the selection criteria where non-ambigous then
;; the list will contain just one entry:
;;
;; (((h "cc-pVDZ-F12")) ; In CAR position a list of one or more aliases
;;  (p (0.2845) (1.0))  ; followed by a list of segments with ...
;;  (p (1.1046) (1.0))  ; a) momentum symbol in CAR
;;  (s (0.1027) (1.0))  ; b) exponents in CADR
;;  (s (0.3258) (1.0))  ; c) coefficients in CADDR, CADDDR, ...
;;  (s (1.159 5.095 33.87)            ; -- these are exponents
;;     (0.202822 0.045308 0.006068))) ; -- a (single) contraction
;;
;; Note  that  in  general  a   segment  may  contain  more  than  one
;; (segmented) contraction.  So that the  length of a segment entry is
;; 1  +  1 +  n  (symbol,  exponents and  n  >  0 contractions).   The
;; Turbomole library  happens to  specify one segment  per contraction
;; even if that means duplicating exponents for several segments.  The
;; NWChem library  has both, segments  with more than  one (segmented)
;; contraction  and generalized  contractions (a  single  segment with
;; more than one contraction).
;;

;;;
;;; Convert symbols to strings and joins them with "/" infix.
;;;
;;; A list  of symbols  is interpreted as  a path relative  to baslib/
;;; directory.  The %load-path is assumed to contain that.
;;;
(define (make-path symbols)
  (string-join (map symbol->string (cons 'baslib symbols)) "/"))

;; (display (make-path '(turbomole basen fe)))
;; (newline)
;; ;; -> "turbomole/basen/fe"

;;
;; Returns an path to the library, if found in %load-path, error
;; otherwise:
;;
(define (find-library symbols)
  (let ((file (make-path symbols)))
    (or (search-path %load-path file)
        (error "baslib/find-library: no such file:"
               file
               %load-path))))  ; or should it be #f?

;; (display (find-library '(turbomole basen he)))
;; (newline)
;; ;; -> "/users/alexei/devel/basis/baslib/turbomole/basen/he"
;;
;; (display (find-library '(turbomole heretic)))
;; (newline)
;; ;; -> ERROR: no such file "turbomole/heretic"


;;;
;;; Returns  a list  of basis-  and ECP  entries like  the  one quoted
;;; above.  To be  used as in (read-library '(turbomole  basen fe)) or
;;; (read-library '(nwchem  6-31g)). At least  in the case  of NWChem,
;;; the library  may contain  auxilary key/value pair  indicating ECPs
;;; associated with these bases.
;;;
(define (read-library symbols)          ; e.g. '(turbomole basen fe)
  "(read-library symbols) => library contents

Examples:

    (read-library '(turbomole basen h))
    (read-library '(nwchem 3-21g))
    (read-library '(nwchem stuttgart_rlc_ecp))
"
  (let ((path (find-library symbols))
        (read (case (car symbols)
                ((turbomole) turbomole-read)
                ((nwchem) nwchem-read)
                (else read))))
    (with-input-from-file path read)))


;;;
;;; NOTE: cache is never  emptied, this basically assumes immutability
;;; of the on-disk library during  a single run.  Total disk usage for
;;; nwchem and turbomole  libs is ~28M. However only  the content that
;;; was read stays in memory.
;;;
(set! read-library (memoize read-library))


;; (define *list* (with-input-from-file (find-library '(turbomole list)) read))
;; (for-each (lambda (a) (read-library `(turbomole basen ,a))) *list*)
;; (for-each (lambda (a) (read-library `(turbomole jbasen ,a))) *list*)

;; (define (pretty-print-bases-for-atom atom)
;;   (begin
;;     (display atom)(newline)
;;     (pretty-print (read-library `(turbomole basen ,atom)))
;;     (newline)))

;; (pretty-print-bases-for-atom 'h)
;; (pretty-print-bases-for-atom 'he)
;; (pretty-print-bases-for-atom 'fe)

;;
;; This returns  a "selector"  --- a function  of a basis  returning a
;; trueish value or an  #f to be e.g. passed as a  first argument to a
;; filter. Strings  and symbols are compared  ignoring the upper/lower
;; case.
;;
(define (match-by-name get-name name)
  (lambda (entry)              ; (get-name entry) is a list of aliases
    (member name (get-name entry) is-equal?))) ; srfi-1 member!

(define (is-equal? a b)
  (or (equal? a b)
      (and (string? a)
           (string? b)
           (string-ci=? a b))
      (and (symbol? a)
           (symbol? b)
           (is-equal? (symbol->string a)
                      (symbol->string b)))))
;;
;; Returns a list (may be empty) of bases in INTERNAL representation:
;;
(define (find-matching bases name) ; e.g. (find-matching (read-library) '(h "SV"))
  (filter (match-by-name basis-name name) bases))

;;
;; Returns a basis record in rep  or fails, to be used as in (basis-ref
;; (read-library ...)  '(fe "SV")):
;;
(define (basis-ref bases alias)
  (match (find-matching bases alias)
         ((one) one)                    ; return the record
         (() (error "baslib/basis-ref: not found:" alias))
         ((one . more) (error "baslib/basis-ref: not unique:" alias))))

;;
;; a) Some entries in the database  refer to other "bases" such as one
;; or  more polarization exponents.  These refs  are represented  by a
;; list with an arrow in CAR position, e.g. (-> (h "PP")).
;;
;; b) In  NWChem libs  some segments  are of SP-type  with one  set of
;; exponents and two sets of coefficients.
;;
(define (expand bases items)
  (define (maybe-expand item)
    ;; An item may be a ref to a list of segments of another basis, or
    ;; an SP segment. Both need expansion:
    (match item
           ;;
           ;; TM-style cross ref, lookup a ref (do we need
           ;; recursion?):
           ;;
           (('-> . alias)
            (basis-segments (basis-ref bases alias)))
           ;;
           ;; SP-segment exponents, s- and p-contractions:
           ;;
           (('SP exps s-cont p-cont)
            (list (cons 'P (list exps p-cont))
                  (cons 'S (list exps s-cont)))) ; expands to two segments
           ;;
           ;; Not a ref, regular case:
           ;;
           (other (list other))))
  (concatenate (map maybe-expand items)))

;;
;; Returns a basis record in GENERALIZED contraction or fails:
;;
(define (find-basis library name)
  "(find-basis library name) => basis

Examples:

    (find-basis '(turbomole basen o) '(o \"SV\"))
    (find-basis '(nwchem 3-21g) 'Cs)
    (find-basis '(nwchem 3-21g) \"Cs_3-21G\")
"
  (let* ((bases (filter basis? (read-library library))) ; strip ECPs
         (basis (basis-ref bases name)) ; record with segments, may fail
         (segments (basis-segments basis)) ; basis payload
         (expanded (expand bases segments)) ; expand (-> ...) and (SP ...)
         (generalized (generalized-contraction expanded)))
    (make-basis (basis-name basis)
                generalized)))

;;;
;;; Returns an ECP record or fails.
;;;
;;; At  the latest  in  NWChem 6.3  the  library file  may include  an
;;; ASSOCIATED_ECP keyword  followed by a  file name (as a  string) of
;;; the same  or another library file  from where the  actual ECP data
;;; should be taken.  Note that several orbital bases  may be combined
;;; with the same ECP (other way around is possible but rarely used).
;;;
(define (find-ecp library name)
  "(find-ecp library name) => ecp record

Example:

    (find-ecp '(nwchem stuttgart_rlc_ecp) 'Lr)
"
  ;;
  ;; Substitutes the  last entry (file name) in  the library "path" by
  ;; another name.
  ;;
  (define (another library sym)
    (reverse (cons sym (cdr (reverse library))))) ; replace the last entry
  ;;
  ;; First, parse the basis set library.  If the library refers to
  ;; another file for ECP data, parse that instead.  Note that
  ;; read-library is memoized so that reading library' has no overhead
  ;; if library' == library.
  ;;
  (let* ((contents (read-library library))
         (contents' (or (and-let* ((str (assoc-ref contents 'ASSOCIATED_ECP))
                                   (sym (string->symbol str))
                                   (library' (another library sym)))
                          (read-library library'))
                        contents))
         (ecp-list-all (filter ecp? contents'))
         (match? (match-by-name ecp-name name))
         (ecp-list (filter match? ecp-list-all)))
    ;;
    ;; Filtered results should  contain exactly one entry. Return that
    ;; entry or raise an error otherwise:
    ;;
    (match ecp-list
      ((one) one)                       ; return the record found
      (()
       (error "baslib/find-ecp: not found:"
              (list library name)
              (map ecp-name ecp-list-all)))
      ((one . more)
       (error "baslib/find-ecp: not unique:"
              (list library name)
              (map ecp-name ecp-list-all))))))

(define (shell-symbol->integer sym)
  (case sym
    ((s) 0)
    ((p) 1)
    ((d) 2)
    ((f) 3)
    ((g) 4)
    ((h) 5)
    ((i) 6)
    (else (error "baslib/shell-symbol->integer: no such case-symbol:" sym))))

;; (define (shell-symbol-less? a b)
;;   (< (shell-symbol->integer a)
;;      (shell-symbol->integer b)))

;; (define shell-symbols (list 's 'p 'd 'f 'g 'h 'i))
;; (display (map shell-symbol->integer shell-symbols))
;; (0 1 2 3 4 5 6)

;;
;; Basis is a list of segmented contractions with symbolic momentum
;; index in the car position.  Returns a list of (cons symbol (cons
;; exponents contractions)):
;;
(define (generalized-contraction shells)
  ;;
  ;; These are internal functions, otherwise I am getting lost in
  ;; those parens.
  ;;
  (define (mk-shell momentum data) ; process shells with the same momentum
    ;; (display "mk-shell: ")
    ;; (display momentum)
    ;; (display " data: ")
    ;; (pretty-print data)
    ;; (newline)
    (let* ((exps (concatenate (map car data))) ; extract all exponents
           (exps (sort exps <))                ; sort them
           (exps (delete-duplicates exps))     ; and delete duplicates
           ;;
           ;; Each segment may give one or more (thus a list of)
           ;; contractions:
           ;;
           (contractions (map (lambda (segment)
                                (mk-contr exps segment)) data)))
      ;;
      ;; Therefore we need to concatenate the resulting contraction,
      ;; generalized == flat:
      ;;
      (cons* momentum exps (concatenate contractions))))

  (define (mk-contr exponents segment) ; make a generalized out of segmented
    ;;
    ;; One  segment may  give to more  than one contraction.  Think of
    ;; generalized  contractions -- they can be  considered an example
    ;; of that.
    ;;
    ;; (display "mk-contr: ")
    ;; (display exponents)
    ;; (display segmented)
    ;; (newline)
    ;;
    (let ((seg-exps (car segment))
          ;;
          ;; This is a list of numeric arrays (ok, lists) each
          ;; representig one segmented contraciton, one for each:
          ;;
          (seg-cont (cdr segment)))
      (let go ((acc '())
               (seg-cont seg-cont))
        (if (null? seg-cont)    ; if the list of segments is empty ...
            acc                 ; return accumulated result
            ;; (concatenate acc)           ; DONT concatenate
            (let* ((alist (map cons seg-exps (car seg-cont))) ; association list (e . c)
                   (c (lambda (e)              ; coeff of an exponents
                        (or (assq-ref alist e) ; is either in the list of pairs
                            0.0)))             ; or zero otherwise
                   (gen-cont (map c exponents))
                   (acc+ (cons gen-cont acc))) ; grow accumulator
              (go acc+ (cdr seg-cont))))))) ; return a list of exponents and contractions

  ;;
  ;; Now  go iteratively over  all shells collecting results  into the
  ;; accumulator:
  ;;
  (let go ((acc '())        ; (accumulator for the) result
           (shells shells)) ; (initial) input
    ;;
    ;; If list of shells is empty ...
    ;;
    (if (null? shells)
        ;;
        ;; ... return accumulated result:
        ;;
        (reverse acc)
        ;;
        ;; Otherwise, pick the angular momentum of the first shell:
        ;;
        (let* ((momentum (caar shells))
               ;;
               ;; This  is a predicate to compare  shell momentum with
               ;; the one extracted above:
               ;;
               (same? (lambda (shell) (eq? momentum (car shell)))))
          ;;
          ;; Partition the shells having the same- and different shell
          ;; momentum.  SRFI-1 partition function returns "multiple
          ;; values" to be used with this somewhat obscure syntax.
          ;;
          (let-values (((same other) (partition same? shells)))
            ;;
            ;; Strip redundant, because identical, shell symbols and
            ;; pass the list of segments list to mk-shell:
            ;;
            (let* ((segments (map cdr same))
                   (generalized-shell (mk-shell momentum segments))
                   ;;
                   ;; Grow accumulator ...
                   ;;
                   (acc+ (cons generalized-shell acc)))
              ;;
              ;; ... and proceed to other shells:
              ;;
              (go acc+ other)))))))
