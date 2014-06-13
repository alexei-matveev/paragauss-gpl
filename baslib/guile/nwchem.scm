;;;
;;; Copyright (c) 2011-2013 Alexei Matveev
;;;
(define-module (baslib guile nwchem)
  #:export (nwchem-read
            make-basis
            basis?
            basis-name
            basis-segments
            ecp?
            ecp-name
            ecp-core
            ecp-shells))

(use-modules (srfi srfi-1)              ; zip
             (srfi srfi-9)              ; records
             (ice-9 rdelim)
             (ice-9 pretty-print))

;;;
;;; Import  "lalr.scm", Guile  2  has  it built  in.   With Guile  1.8
;;; %load-path needs to be set to contain the repository root.
;;;
(cond-expand
 ((not guile-2)
  (use-modules (baslib guile lalr)))    ; more recent?
 (else
  (use-modules (system base lalr))))    ; built in

;;;
;;; Record  constructors  are  macros  in  Guile 2.   So  that  record
;;; definitions must appear before use of constructors (see production
;;; rules in (lalr-parser ...).
;;;
;;; This declares a record type to represent an orbital basis.  FIXME:
;;; need to put it elsewhere to be used also in turbomole.scm:
;;;
(define-record-type basis
  (make-basis name segments)            ; constructor
  basis?                                ; predicate
  (name basis-name)                     ; fields with "getters"
  (segments basis-segments))

;;;
;;; This  declares  a  record  type  to  represent  an  ECP  (not  the
;;; corresponding orbital basis!).  FIXME: need to put it elsewhere to
;;; be used also in turbomole.scm:
;;;
(define-record-type ecp
  (make-ecp name core shells)           ; constructor
  ecp?                                  ; predicate
  (name ecp-name)                       ; fields with "getters"
  (core ecp-core)                       ; (no setters, please)
  (shells ecp-shells))

;;
;; Slurps the whole input into a list and returns in reverse order it
;; if called as in (tokenize '() read eof-object?):
;;
(define (tokenize acc reader eoi?)
  (let ((token (reader)))
    (if (eoi? token)
        acc
        (tokenize (cons token acc) reader eoi?))))

;;
;; Tokenizes a line, by cons'ing additional entries into acc:
;;
(define (lex-line acc line)
  (with-input-from-string line
    (lambda () (tokenize acc read eof-object?))))

;;
;; New lines seem to need to be treated as tokens in NWChem library,
;; the built in Scheme reader strips them as whitespace, we need a
;; different one. This, when called as (tokenize-eol '()) slurps the
;; whole input and returns all tokens by cons'ing them into a
;; (reversed) list:
;;
(define (tokenize-eol acc)
  (let ((line (read-line)))
    (if (eof-object? line)
        acc
        (let ((fields (lex-line acc line)))
          (tokenize-eol (cons '*eol* fields)))))) ; local convention for end of line
;;        (tokenize-eol fields)))))

(define (slurp)
  (reverse (tokenize-eol '())))

;; (define *x* (slurp))
;; (pretty-print *x*)
;; (exit 1)

;;
;; FIXME: slurps the whole input into a list, yileds that elementwise:
;;
(define (make-greedy-tokenizer)
  (let ((*buf* (slurp)))
    ;; (pretty-print *buf*)
    (lambda ()
      (if (null? *buf*)
          '*eoi*                 ; end of input convention of lalr.scm
          (let ((token (car *buf*)))
            (set! *buf* (cdr *buf*))
            (make-token token))))))

(define (location)
  'undefined)

(define (keyword kw)
  (case kw
    ((basis) 'BASIS)
    ((end)   'END)
    ((ecp)   'ECP)
    ((*eol*) 'EOL)                      ; another local convention
    (else    'SYMBOL)))

(define (make-token token)
  (cond
   ((eq? token '*eoi) '*eoi*)
   ((symbol? token)    (make-lexical-token (keyword token) (location) token)) ; basis, ecp, end, s, p, ...
   ((string? token)    (make-lexical-token 'STRING (location) token)) ; quoted string
   ((exact? token)     (make-lexical-token 'EXACT (location) token)) ; integers
   ((inexact? token)   (make-lexical-token 'INEXACT (location) token)) ; real numbers
   (#t                 (make-lexical-token token (location) token)) ; pass through terminals
   ))

;; (tokenize read eof-object?)
;; (exit 1)

;;
;; Reverse engineered grammar:
;;
(define (make-nwchem-parser)
  (lalr-parser
   ;;
   ;; Terminal symbols:
   ;;
   (EOL BASIS ECP END SYMBOL STRING EXACT INEXACT)

   ;;
   ;; Productions:
   ;;

   ;;
   ;; All files contain the leading comment- and empty lines of which
   ;; only the EOL is what is being parsed. There also seem to be some
   ;; comments and stray empty lines following the ASSOCIATED_ECP
   ;; directive:
   ;;
   (input
    (eol+ item) : (list $2)
    (input item) : (cons $2 $1))

   (item
    (ecp) : $1
    (basis) : $1
    (associated-ecp) : $1)

   (eol+
    (EOL) : 'nothing
    (eol+ EOL) : 'nothing)

   ;;
   ;; A line like
   ;;
   ;;   ASSOCIATED_ECP "stuttgart_rsc_1997_ecp"
   ;;
   ;; in a basis file refers to the associated ECP library (which may
   ;; be the very same file).  This line is followed by comments lexed
   ;; into EOL tokens. FIXME: this entry introduces the third kind of
   ;; entries into the parsed list (in addition to basis- and ecp
   ;; records) --- a key-value pair.
   ;;
   (associated-ecp
    (SYMBOL STRING eol+) : (cons $1 $2))

   ;;
   ;; BASES. Each shell in the library (redundantly) specifies atomic
   ;; symbol. CAAR $4 extracts one of them assuming all are the
   ;; same. MAP CDR $4 strips them as redundant from the real data.
   ;;
   (basis
    (BASIS header EOL shells END EOL) : (make-basis (cons (caar $4) $2)
                                                    (reverse (map cdr $4))))

   (header
    (STRING SYMBOL) : (list $1 $2))

   (shells
    (shell) : (list $1)
    (shells shell) : (cons $2 $1))

   ;;
   ;; The first symbol is the atom name, duplicated in the header
   ;; string:
   ;;
   (shell
    (SYMBOL SYMBOL EOL rows) : (make-shell $1 $2 $4))

   ;;
   ;; Real numbers come either in rows, either as a pair (exponent
   ;; coeff), or in triples, (sp-exponent s-coeff-s p-coeff), or even
   ;; longer rows corresponding to generalized contractions (see
   ;; e.g. WTBS):
   ;;
   (row
    (doubles EOL) : (reverse $1))

   (rows
    (row) : (list $1)
    (rows row) : (cons $2 $1))

   (doubles
    (INEXACT) : (list $1)
    (doubles INEXACT) : (cons $2 $1))

   ;;
   ;; ECP,  return a record  with a name  (Lr "Lr_Stuttgart_RLC_ECP"),
   ;; number  of  core  electrons  and  the  data for  shell  specific
   ;; potentials:
   ;;
   (ecp
    (ECP ecp-header
     ecp-shells
     END EOL) : (let ((core (car $2))
                      (name (cdr $2))
                      (shells (reverse $3)))
                  (if (eq? (caar shells) 'ul) ; some NWChem ECPs omit local part
                      (make-ecp name
                                core
                                shells)
                      (make-ecp name
                                core
                                (cons (list 'ul) ; empty list of triples in CDR
                                      shells)))))

   ;; ... "B_Stuttgart_RLC_ECP"
   ;;      B  nelec  2
   (ecp-header
    (STRING EOL SYMBOL SYMBOL EXACT EOL) : (list $5 ; core nelec
                                                 $3 ; atomic symbol
                                                 $1)) ; type string

   ;;
   ;; The first symbol is the atom name, duplicated in the header
   ;; string:
   ;;
   (ecp-shell
    (SYMBOL SYMBOL EOL ecp-rows) : (cons $2 $4))

   (ecp-shells
    (ecp-shell) : (list $1)
    (ecp-shells ecp-shell) : (cons $2 $1))

   (ecp-row
    (EXACT INEXACT INEXACT EOL) : (list $1 $2 $3))

   (ecp-rows
    (ecp-row) : (list $1)
    (ecp-rows ecp-row) : (cons $2 $1))))

;;
;; Other  libraries (Turbomole) specify  more than  one basis  name to
;; address  a basis,  so keep  basis identifiers  as a  list  (here of
;; length  2). The  first  identifier  in the  list  is the  lowercase
;; version of the basis name  (FIXME: we loose the distinction between
;; SPHERICAl/CARTESIAN):
;;
(define (make-aliases header)
  ;; (write header) (newline)
  (let* ((basis-name-mixed (car header))
         (basis-name-lower (string-downcase basis-name-mixed))
         (basis-name-upper (string-upcase basis-name-mixed)))
    ;;
    ;; Augment the list of aliases by a lowercase version of the basis
    ;; name as an alias, keep original ID too:
    ;;
    (list basis-name-mixed header)))

(define (make-shell atom shell-symbol rows)
  (cons* atom shell-symbol (apply zip rows)))

(define (nwchem-read)
  ((make-nwchem-parser) (make-greedy-tokenizer) error))

;; (let ((data (nwchem-read)))
;;   (begin
;;     (pretty-print (reverse data))
;;     (newline)))