#!/usr/bin/guile
!#
;;;
;;; Copyright (c) 2011-2013 Alexei Matveev
;;;

;; FIXME: hardcoded path:
(set! %load-path (cons "/users/alexei/devel/bases" %load-path))

(use-modules ((baslib guile paragauss)
	      #:select (qm-write-input
			qm-make-basis
			*h2o*)))

(use-modules ((baslib guile baslib)
	      #:select (find-basis)))

(use-modules (ice-9 pretty-print))

;; (pretty-print (find-basis '(turbomole basen o) "def-SVP"))
;; (exit 1)

(define *oxygen*
  (qm-make-basis
   (find-basis '(turbomole basen o) '(o "def-SVP"))
   (find-basis '(turbomole jbasen o) '(o "universal"))))

(define *hydrogen*
  (qm-make-basis
   (find-basis '(turbomole basen h) '(h "def-SVP"))
   (find-basis '(turbomole jbasen h) '(h "universal"))))

(define *oxygen*
  (qm-make-basis
   (find-basis '(nwchem 6-31g) "O_6-31G")
   (find-basis '(nwchem ahlrichs_coulomb_fitting) "O_Ahlrichs_Coulomb_Fitting")))

(define *hydrogen*
  (qm-make-basis
   (find-basis '(nwchem 6-31g) "H_6-31G")
   (find-basis '(nwchem ahlrichs_coulomb_fitting) "H_Ahlrichs_Coulomb_Fitting")))

(qm-write-input *h2o*)			; example input, no bases
(qm-write-input *oxygen*)
(qm-write-input *hydrogen*)
