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
;;; Copyright (c) 2011-2013 Alexei Matveev
;;;
;;;
;;; You may need to set GUILE_LOAD_PATH=.. for this to work:
;;;
(use-modules (guile comm))

;(define libguile-comm (dynamic-link "./libguile-comm.so"))
;(dynamic-call "guile_comm_init_module" libguile-comm)

;; returns MPI_COMM_WORLD:
(define world
  (comm-init (command-line)))

;; size and rank are communicator specific:
(define size (comm-size world))
(define rank (comm-rank world))

;;
;; Cyclic shift, send to the right, receive from the left:
;;
(let
    ((right-peer (modulo (+ rank 1) size))
     (left-peer (modulo (- rank 1) size))
     (token (list "token of rank" rank "with payload" (+ rank 0.1))) ; token to pass, rank specific
     (message-tag 999))      ; arbitrary number, same for every worker
  (let loop ((i 0) (token token))
    (if (< i size)
        (let
            ((token (comm-send-recv world right-peer left-peer message-tag token)))
          (display (list rank token))(newline)
          (loop (+ i 1) token)))))

;;
;; Compute PI in parallel (by calling native code with our communicator):
;;
(define pi
  (comm-pi world 10000000))

;;
;; Loop over ranks with communication barriers inbetween
;; for proper output formatting (does not always work):
;;
(define (for-each-rank world proc)
  (let loop ((p 0) (size (comm-size world)))
    (if (< p size)
        (begin
          (if (= p (comm-rank world))  ; then it is my turn to act ...
              (let ((result (proc)))
                (display (list world (comm-rank world) (comm-size world) result))
                (newline)))
          (comm-barrier world)   ; others wait here until I finish ...
          (loop (+ p 1) size)))))

(for-each-rank world (lambda () pi))

; (exit 0)

;; each worker is either even or odd:
(define color (modulo rank 2))

;; group even and odd workers into communication groups:
(define country (comm-split world color))
(comm-set-name country (if (= color 0) "even" "odd"))
;; (display country)(newline)

;; let the two groups compete in computing the pi again:
(define pi-2 (comm-pi country 1000))

;; even output first, odd second:
(comm-barrier world)
(if (= color 0) (for-each-rank country (lambda () pi-2)))
(comm-barrier world)
(if (= color 1) (for-each-rank country (lambda () pi-2)))

;; release communicators:
(comm-free country)

(comm-barrier world)

;;
;; FIXME: ugly abstraction:
;;
;; ping-pong, FIXME: this will dead-lock with odd number of workers:
(let
    ((right (modulo (+ rank 1) size))   ; right peer
     (left (modulo (- rank 1) size))    ; left peer
     (token (+ rank 1000))
     (tag 999))
  (if (= color 0)                       ; even send, odd receive ...
      (comm-send world right tag token) ; on even
      (let ((ping (comm-recv world left tag))) ; on odd
        (display (list rank "recv ping" ping))
        (newline)))
  (if (= color 1)                       ; odd send, even receive ...
      (comm-send world left tag token)  ; on odd
      (let ((pong (comm-recv world right tag))) ; on even
        (display (list rank "recv pong" pong))
        (newline))))

;; required by MPI:
(comm-finalize)
