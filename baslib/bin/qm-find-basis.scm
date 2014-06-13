#!/usr/bin/env guile
;;; -*- mode: scheme; -*- vim: set syntax=scheme:
!#
;;;
;;; Copyright (c) 2011-2013 Alexei Matveev
;;;
;;; Make sure  the shebang line  and %load-path are set  correctly. In
;;; particular "pkgdatadir"  between @-markers should  be replaced by,
;;; say, a path to PG repository.
;;;
;;; Load  user-specific config, this  may also  set the  %load-path to
;;; search for files/modules.  This should work with both 1.8 and 2.0.
;;; Use  absolute paths  when  extending %load-path  ---  there is  no
;;; globbing for ~/ at compile time.
;;;
(cond-expand
 ((not guile-2) (use-modules (ice-9 syncase))) ; eval-when for 1.8
 (else))                                       ; else nothing


(eval-when
 (eval load compile)      ; extend %load-path at compile time with 2.0
 (set! %load-path (cons "@pkgdatadir@" %load-path))
 (let ((conf (string-append (getenv "HOME") "/.qmrc")))
   (if (file-exists? conf)
       (load conf))))

(use-modules (ice-9 pretty-print)
             (ice-9 getopt-long)
             ((baslib guile paragauss)
              #:select (qm-write-input
                        qm-find-basis
                        qm-find-ecp)))

(define (help)
  (display
"USAGE:

  qm-find-basis [OPTIONS] ATOM TYPE [AUXILARY]

with ATOM being  the atomic symbol, TYPE being the  name of an orbital
basis or an ECP type.  AUXILARY is  a name of a fit basis, defaults to
TYPE when omitted. For ECPs AUXILARY  is the name of the orbital basis
with  TYPE as  the default  and  the fit  basis is  taken from  NWChem
ahlrichs_coulomb_fitting collection unconditionally.

OPTIONS:

    --library=(turbomole|nwchem)
    --ecp

EXAMPLES:

    qm-find-basis fe def-SVP
    qm-find-basis --ecp au ecp-60-mwb ecp-60-mwb-opt
    qm-find-basis --library=nwchem Fe 6-31g ahlrichs_coulomb_fitting
    qm-find-basis --library=nwchem --ecp Au stuttgart_rsc_1997_ecp
"))
;; (help)

(define *options*
  '((library (value #t))                ; accept --library=name
    (ecp (value #f))))                  ; accept --ecp

(define (run)
  (let* ((opts (getopt-long (command-line) *options*)) ; parse command line
         (args (option-ref opts '() '())) ; extract list of (positional) args
         (library (option-ref opts 'library "turbomole")) ; use turbomole lib by default
         (ecp (option-ref opts 'ecp #f)))
    ;;
    ;; Write the basis in PG namelist format to stdout:
    ;;
    (qm-write-input
     (if ecp
         (apply qm-find-ecp library args)
         (apply qm-find-basis library args))))) ; first arg is the library name

(catch
 #t                                  ; #t is a placeholder
 run                                 ; try executing this ...
 (lambda (key . rest)                ; ... except when an error occurs
   (help)                            ; then print help ...
   (apply throw (cons key rest))))   ; ... and throw again
