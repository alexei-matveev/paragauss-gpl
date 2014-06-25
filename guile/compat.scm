;;;
;;; ParaGauss, a program  package for high-performance computations of
;;; molecular systems
;;;
;;; Copyright (C) 2014     T. Belling,    T. Grauschopf,    S. Krüger,
;;; F. Nörtemann,      M. Staufer,      M. Mayer,      V. A. Nasluzov,
;;; U. Birkenheuer,       A. Hu,       A. V. Matveev,      A. V. Shor,
;;; M. S. K. Fuchs-Rohr,         K. M. Neyman,        D. I. Ganyushin,
;;; T. Kerdcharoen,    A. Woiterski,   A. B. Gordienko,   S. Majumder,
;;; M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
;;; M. Roderus, N. Rösch
;;;
;;; This  program is  free software;  you can  redistribute  it and/or
;;; modify  it under  the  terms  of the  GNU  General Public  License
;;; version 2 as published by the Free Software Foundation [1].
;;;
;;; This program  is distributed in the  hope that it  will be useful,
;;; but  WITHOUT ANY WARRANTY;  without even  the implied  warranty of
;;; MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;; General Public License for more details.
;;;
;;; [1] http://www.gnu.org/licenses/gpl-2.0.html
;;;
;;; Please see the accompanying LICENSE file for further information.
;;;
;;;
;;; Copyright (c) 2011-2013 Alexei Matveev
;;;
(define-module (guile compat))

(cond-expand
 ((not guile-2)
  ;;
  ;; Code to backport selected features to 1.8 goes here:
  ;;
  (use-syntax (ice-9 syncase))

  ;;
  ;; Thanks to taylanub at #guile:
  ;;
  (define-syntax define-syntax-rule
    (syntax-rules ()
      ((define-syntax-rule (keyword . pattern) template)
       (define-syntax keyword
         (syntax-rules ()
           ((keyword . pattern) template))))))

  ;;
  ;; Not in the define-module form:
  ;;
  (export-syntax define-syntax-rule)

  ;;
  ;; When/uless were added in 2.0:
  ;;
  (define-syntax-rule (when test stmt stmt* ...)
    (if test (begin stmt stmt* ...)))

  (export-syntax when)

  (define-syntax-rule (unless test stmt stmt* ...)
    (if (not test) (begin stmt stmt* ...)))

  (export-syntax unless)

  ;;
  ;; Re-export all symbols from (ice-9 syncase). Credit to mark_weaver
  ;; at #guile. Otherwise the importer of this module, (guile compat),
  ;; still  needs to  explicitly import (ice-9  syncase) in  order for
  ;; define-syntax-rule to work.
  ;;
  (let ((names (module-map (lambda (name var) name)
                           (resolve-interface '(ice-9 syncase)))))
    (module-re-export! (current-module) names)))
 (else)
 ;;
 ;; Nothing here.
 ;;
 )
