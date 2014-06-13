;;;
;;; Copyright (c) 2011-2013 Alexei Matveev
;;;
(define-module (baslib guile turbomole)
  #:export (turbomole-read))
;;
;; export GUILE_LOAD_PATH=~/devel/basis/guile
;; for this to work:
;;
(use-modules (ice-9 pretty-print)
             ;; These need properly set %load-path:
             ((baslib guile nwchem)
              #:select (make-basis
                        make-ecp)) ; FIXME: a better place for record decls?
             (baslib guile lalr)) ;; (load "lalr.scm")

; (define (tokenize reader eoi?)
;   (let loop ((token (reader)))
;     (if (not (eoi? token))
;       (begin
;         (write (make-token token))
;         (newline)
;         (loop (reader))))))

(define (tokenizer)
  (make-token (read)))

(define (location)
  (port-line (current-input-port)))

(define (make-token token)
  (cond
   ((eof-object? token)  '*eoi*)  ; end of input convention of lalr.scm
   ((eq? token '$end)    '*eoi*)
   ((eq? token '$ehtdef) '*eoi*)
   ((eq? token '$basis)  (make-lexical-token 'BASIS   (location) token))
   ((eq? token '$jbas)   (make-lexical-token 'BASIS   (location) token)) ; same as $basis
   ((eq? token '$ecp)    (make-lexical-token 'ECP     (location) token)) ; start of a section
   ((eq? token '*)       (make-lexical-token 'STAR    (location) token)) ; separator
   ((eq? token '->)      (make-lexical-token 'ARROW   (location) token)) ; ->, as in -> h "6-31g"
   ((symbol? token)      (make-lexical-token 'SYMBOL  (location) token)) ; s, p, d, ...
   ((string? token)      (make-lexical-token 'STRING  (location) token)) ; e.g. "def-SVP"
   ((exact? token)       (make-lexical-token 'EXACT   (location) token)) ; integers
   ((inexact? token)     (make-lexical-token 'INEXACT (location) token)) ; real numbers
   (#t                   (make-lexical-token token    (location) token)))) ; pass through terminals

; (tokenize read eof-object?)

;;
;; Reverse engineered grammar:
;;
(define (make-turbomole-parser)
  (lalr-parser
   ;;
   ;; Terminal symbols:
   ;;
   (BASIS STAR ARROW SYMBOL STRING EXACT INEXACT ECP)
   ;;
   ;; Productions:
   ;;
   (input
    (BASIS bases)   : $2
    (BASIS bases ECP ecps) : (append $2 $4))

   (bases
    (STAR)          : '()
    (bases STAR)    : $1
    (basis)         : $1
    (bases basis)   : (cons $2 $1))

   (basis
    (names STAR data) : (make-basis $1
                                     (reverse $3)))

   (data
    (shells)        : $1
    (refs)          : $1

    ;;
    ;; Only in "ne" lib do literal exponents follow refs.  There was a
    ;; STAR   between  them,  which  was  removed  to  reduce  grammar
    ;; ambiguity:
    ;;
    (refs shells) : (list 'refs: $1 'literal: $2))

   (shells
    (shell)         : (list $1)
    (shells shell)  : (cons $2 $1))

   (shell
    (EXACT SYMBOL doubles) : (make-shell $1 $2 $3)

    ;;
    ;; There appear some "orphaned" exponents without coeff
    ;; scattered over the basis libraries.
    ;; We assume coeff = 1.0:
    ;;
    (EXACT SYMBOL INEXACT) : (make-shell $1 $2 (list (cons $3 1.0))))

   (doubles
    (INEXACT INEXACT)         : (acons $1 $2 '())
    (doubles INEXACT INEXACT) : (acons $2 $3 $1))

   (names
    (name)          : (list $1)
    (names name)    : (cons $2 $1))

   (name
    (SYMBOL STRING) : (list $1 $2))

   (refs
    (ref)           : (list $1)
    (refs ref)      : (cons $2 $1))

   (ref
    (ARROW name)    : (cons $1 $2))

   (ecps
    (STAR)      : '()
    (ecps STAR) : $1
    (ecp)      : (list $1)
    (ecps ecp) : (cons $2 $1))

   (ecp
    (names
     STAR
     assign assign
     ecp-shells) : (make-ecp $1 $3 (reverse $5))

    ;; Spin-Orbit ECPs, need another record field:
    (names
     STAR
     assign assign assign
     ecp-shells
     STAR
     ecp-shells) : (list 'so-ecp: $1 $3 $6 $8))

   (assign                              ; e.g. ncore = 60, lmax = 3
    (SYMBOL SYMBOL EXACT) : $3)

   (ecp-shells
    (ecp-shell)            : (list $1)
    (ecp-shells ecp-shell) : (cons $2 $1))

   (ecp-shell
    (SYMBOL ecp-rows) : (cons $1 (reverse $2)))

   (ecp-rows
    (ecp-row)          : (list $1)
    (ecp-rows ecp-row) : (cons $2 $1))

   (ecp-row                             ; coeff power exponent
    (INEXACT EXACT INEXACT) : (list $2 $3 $1)))) ; power exponent coeff

(define (make-shell count symbol pairs)
  (if (= count (length pairs))          ; omit redundant count
      (list symbol                      ; shell symbol
            (map car pairs)             ; exponents
            (map cdr pairs))            ; (segmented) contraction
      (error "number and length disagree" (list count symbol pairs)))) ; report inconsistency

; (pretty-print
;   (turbomole-parser tokenizer error))
; (newline)
(define (turbomole-read)
  ((make-turbomole-parser) tokenizer error))

;; options for vim:sw=2:expandtab:smarttab:autoindent:syntax
