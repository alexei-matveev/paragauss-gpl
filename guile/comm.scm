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
(define-module (guile comm)
  ;;
  ;; The module (guile comm internal)  is set up on guile startup, see
  ;; the C sources in guile-qm.c (or guile-comm.c):
  ;;
  #:use-module (guile comm internal)    ; see libguile-comm.c
  #:re-export                           ; from (guile comm internal)
  (comm-init
   comm-finalize
   comm-rank
   comm-size
   comm-barrier
   comm-bcast
   comm-send-recv
   comm-send
   comm-recv
   comm-split
   comm-free
   comm-set-name
   comm-pi)
  #:export
  (critical))

;;
;; Call  prog  with  for  each  rank  in  sequence  with  comm-barrier
;; inbetween:
;;
(define (critical world prog)
  (let loop ((rank 0))
    (if (< rank (comm-size world))
        (begin
          (if (= rank (comm-rank world)) ; My turn to ...
              (prog))                    ; ... execute prog.
          (comm-barrier world)           ; Others wait here.
          (loop (+ rank 1))))))

;;;
;;; Apply selectively func to elements of the list of args if index of
;;; the element i == k mod n. For other elements just return #f:
;;;
;;; (round-robin-map (lambda (i) (* 10 i)) (iota 9) 4 1)
;;; => (#f 10 #f #f #f 50 #f #f #f)
;;;
(define (round-robin-map func args n k)
  (let loop ((args args)
             (i 0)
             (acc '()))
    (if (null? args)
        (reverse acc)
        (let ((res (and (= k (modulo i n))
                        (func (car args)))))
          (loop (cdr args)
                (+ 1 i)
                (cons res acc))))))
