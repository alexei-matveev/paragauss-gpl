#!/usr/bin/env guile-qm
!#
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
;;
;; The shebang line above assumes the interpreter is in the $PATH:
;;
;;      mpirun -np 8 ./test-guile-qm.scm
;;
;; Alternatively run the script specifying the interpeter explicitly:
;;
;;      mpirun -np 8 ../guile-qm -s test-guile-qm.scm
;;

;;
;; This runs ParaGauss in a temp dir.
;;
;; FIXME: every process creates a directory named differently
;;
;; Note that PG assumes the same directory was created on all hosts,
;; so that the separate Scheme process need to agree on directory
;; naming when creating one:
;;
;;
(define (run world)
    (let ((name (tmpnam)))              ; abuse tmpnam for generating directory names
        (begin
            (mkdir name)                ; create the directory
            (setenv "TTFSTMP" name)     ; make sure PG uses it
            (setenv "TTFSOUTPUTDIR" name)
            (qm-run world))))           ; and execute the program

;;
;; Not idempotent, calls MPI_Init. Returns MPI_COMM_WORLD:
;;
(define world (qm-init))

;;
;; By default input is taken from the file named "input".
;; I assume it is not a good idea when several processes
;; try to copy the file in parallel:
;;
(if (= 0 (comm-rank world))
    (copy-file "input.1" "input"))

;;
;; Runs the same input several time, almost idempotent:
;;
(run world)
(run world)

;;
;; Split the "world" in even and odd "countries":
;;
(define rank (comm-rank world))
(define color (modulo rank 2))
(define country (comm-split world color))

; (if (= color 1)
(run country)

;;
;; Release the communicator:
;;
(comm-free country)
; (exit 0)

(run world)
(run world)
(run world)
;; (exit 0)

(if (= 0 (comm-rank world))
    (copy-file "input.2" "input"))
(run world)
(run world)

(if (= 0 (comm-rank world))
    (copy-file "input.1" "input"))
(run world)
(run world)

(if (= 0 (comm-rank world))
    (copy-file "input.2" "input"))
(run world)
(run world)

;;
;; Not idempotent, calls MPI_Finalize
;;
(qm-finalize world)
