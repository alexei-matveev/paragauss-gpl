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
;;; Module/package systems differ between implementaitons:
;;;
(cond-expand
 (guile
  ;;
  ;; Primary implementation:
  ;;
  (define-module (guile scheduling)
    #:use-module (ice-9 pretty-print)
    #:use-module (srfi srfi-1)          ; list manipulations
    #:use-module ((guile utils) #:select (isqrt))
    #:export (qm-mpts->npts
              grid-shape))
  ;; To get rid of deprecation warning make import of syncase
  ;; conditional:
  (cond-expand
   ((not guile-2) (use-modules (ice-9 syncase))) ; define-syntax
   (else)))                                      ; else nothing
 (else
  ;;
  ;; MzScheme, aka PLT Scheme,  aka Racket (needs cond-expand macro in
  ;; ~/.mzschmerc):
  ;;
  (require (lib "1.ss" "srfi"))
  (define (1+ x) (+ 1 x))
  (define (sorted? lst pred?) #t)))     ; FIXME: lies!

;; Use this for debug prints.  Unquote the pretty-print form to enable
;; prints   and   evaluate   (potentially  expensive)   arguments   of
;; debug-print macro:
(define-syntax debug-print
  (syntax-rules ()
    ((_ obj)
     '(pretty-print obj))))

;;;
;;; See "Scheduling malleable and nonmalleable parallel tasks", Walter
;;; Ludwig, Prasoon  Tiwari, SODA '94 Proceedings of  the fifth annual
;;; ACM-SIAM symposium on Discrete algorithms.
;;;

;;
;; Task  is  an "abstract"  data  type  (currently  represented by  an
;; integer size parameter  of the task).  And "partition"  is a subset
;; of workers  (currently represented by an integer  worker count). It
;; is assumed the two are  sufficient to determine the cumulative work
;; function:
;;
(define (make-world size) size)
(define (make-task size) size)
(define (world-size world) world)
(define (task-size task) task)
;; (define (make-world size) (cons 'world size))
;; (define (make-task size) (cons 'task size))
;; (define (world-size world) (cdr world))
;; (define (task-size task) (cdr task))

(define (work-function task world)
  ;;
  ;; Work function should better be a non-decreasing function of
  ;; number of workers:
  ;;
  (let ((n (task-size task))
        (p (world-size world)))
    (+ (* n n n)                        ; O(n**3) work
       (* n n (- p 1)))))               ; O(pn**2) overhead

(debug-print
 (work-function (make-task 100) (make-world 16)))

;;
;; We may use these properties for optimizations, like binary search
;; on bounds:
;;
(define (non-decreasing? lst) (sorted? lst <))
(define (non-increasing? lst) (sorted? lst >))
(define (regular-case? partitions time-matrix time-series)
  (let* ((sizes
          (map world-size partitions))
         (load-matrix
          (map (lambda (times) (map * times sizes))
               time-matrix)))
    (and (non-decreasing? sizes)
         (non-decreasing? time-series)
         (every non-increasing? time-matrix)
         (every non-decreasing? load-matrix))))

(define (feasible-partitions world)
  ;;
  ;; Returns a list of feasible subsets of workers to execute the
  ;; task. Here all feasible partitions of the size ranging from one
  ;; to the number of workers (FIXME: partition as a function of size
  ;; only here, feasible partitions do not depend on task):
  ;;
  (let* ((size (world-size world))
         (sizes (map 1+ (iota size))))
    (map make-world sizes)))

(debug-print
 (feasible-partitions (make-world 10)))

;;
;; Computes f(a, b) for all pairs from the list of a's and the list of
;; b's. Returns a nested list [[f(a, b) | b <- bs] | a <- as]:
;;
(define (outer f as bs)
  (map (lambda (a) (map (lambda (b) (f a b)) bs)) as))

(debug-print
 (outer cons '(1 2 3) '(a b c d)))

;;
;; Returns #f  if any of evaluations  of maybe-proc do  so.  Note that
;; this is not an "and-map" which would return the last result or #f.
;;
(define (maybe-map maybe-proc args)
  (let loop ((args args)
             (acc '()))
    (if (null? args)
        (reverse acc)
        (let ((res (maybe-proc (car args))))
          (and res
               (loop (cdr args) (cons res acc)))))))

(define (test n) (if (= 0 n) #f (/ 1 n)))
(debug-print (maybe-map test '(1 2 3 4)))
(debug-print (maybe-map test '()))
(debug-print (maybe-map test '(1 2 0 4)))

(define (maybe-find threshold candidates cost size)
  ;;
  ;; May return  #f if the list of candidates is  empty or rather when
  ;; there   is  no  candidate   with   a  cost  below   or  equal  to
  ;; threshold. Otherwise return one of the candidates of the smallest
  ;; size with the cost not above threshold.
  ;;
  ;; O(P) complexity, with P being the number of worker candidates.
  ;;
  (let ((choose
         (lambda (new old)
           ;;
           ;; Check  if the  new is "more  appropriate" than  the old,
           ;; return one or another depending on that:
           ;;
           (if (> (cost new) threshold) ; filter them out
               old
               (if (or (not old)    ; not #f, only then compare sizes:
                       (< (size new) (size old)))
                   new
                   old)))))
    (fold choose #f candidates)))

(define (mpts->npts work-function tasks world)
  ;;
  ;; Suggest a partition (number  of workers) for each of three tasks.
  ;; Expects  a list of task  descriptors and a desriptor  of a worker
  ;; set.
  ;;
  (let* ((np
          (world-size world))

         (mean-work
          (lambda (tasks partitions)
            (let* ((loads (map work-function tasks partitions))
                   (total (apply + loads)))
              (/ total np))))

         (cost-function
          (lambda (task world)
            ;;
            ;; Time span is a total work divided by number of workers
            ;; (and conversly). Should better be a non-increasing
            ;; function of number of workers:
            ;;
            (let ((w (work-function task world))
                  (p (world-size world)))
              (/ w p))))

         (partitions
          (feasible-partitions world))

         ;;
         ;; Costs of every task on every partition. O(MP) size.
         ;;
         (time-matrix
          (outer cost-function tasks partitions))

         (time-series
          (delete-duplicates (sort (concatenate time-matrix) <)))

         ;;
         ;; This function finds one subset of the world of the
         ;; smallest size that allows execution of the task in time
         ;; below tau.  In case it does not, it returns #f.  O(P)
         ;; complexity inherited from maybe-find.
         ;;
         (maybe-find-partition
          (lambda (task tau)
            (let ((cost (lambda (world) (cost-function task world)))) ; cost as func of partition
              (maybe-find tau partitions cost world-size))))

         ;;
         ;; This function finds one subset of the world of the
         ;; smallest size that allows execution of the task in time
         ;; below tau. O(MP) complexity with M being the number of
         ;; tasks.
         ;;
         (maybe-find-partitions
          (lambda (tau)
            (maybe-map (lambda (task)
                         (maybe-find-partition task tau))
                       tasks)))

         ;;
         ;; O(MP) complexity, with M being the number of tasks and
         ;; P being the number of workers, inherited from
         ;; maybe-find-partitions.
         ;;
         (maybe-work
          (lambda (tau)
            (let ((partitions (maybe-find-partitions tau)))
              (and partitions (mean-work tasks partitions)))))

         ;;
         ;; O(MP), inherited from maybe-work.
         ;;
         (maybe-bound
          (lambda (tau)
            (let ((w (maybe-work tau)))
              (and w (max tau w)))))

         (regular
          (regular-case? partitions time-matrix time-series)))

    (let* ((series
            (if (not regular)
                ;;
                ;; Brute force minimization. O((MP)**2) complexity,
                ;; the product of that for maybe-bound and the size of
                ;; time-series:
                ;;
                (error "brute force?" time-series)
                ;;
                ;; Otherwise use bisection method to bracket the
                ;; interval. This should return a list of length one
                ;; or two:
                ;;
                (bisect-series maybe-work time-series)))

           (omega
            (apply min (filter-map maybe-bound series)))

           (optimal
            (maybe-find-partitions omega))

           (times
            (map cost-function tasks optimal)))

      ;; (list tasks '@ world '-> 'omega: omega 'partitions: optimal
      ;; 'times: (map exact->inexact times))
      (map cons optimal times))))

(define (bisect-series maybe-work time-series)
  (let* ((less?
          (lambda (tau)
            ;;
            ;; Function, (less?  tau), on a time series to answer if
            ;; the (increasing) tau is still "less" than the
            ;; corresponding (decreasing) work:
            ;;
            (let ((work (maybe-work tau))) ; may return #f for small tau
              (if work                     ; if work is a number ...
                  (< tau work)             ; max(tau, work) == work
                  #t)))) ; otherwise tau too small for work to be meaningfull

         ;;
         ;; Here we bracket the location of the minima to the interval
         ;; [i, j], usually j = i + 1.  This is position where T(j) >=
         ;; W(T(j)), found by bisection:
         ;;
         (j (list-bisect time-series less?))

         ;;
         ;; This is the previous position where T(i) < W(T(i)):
         ;;
         (i (max (- j 1) 0))

         ;;
         ;; Somewhere inbetween is the minimum. A shorter time series,
         ;; ideally a list of length one or two:
         ;;
         (series
          (list-range time-series i j)))
    series))

(define (list-bisect lst less?)
  (let ((n (length lst))
        (less? (lambda (i) (less? (list-ref lst i)))))
    (bisect 0 (- n 1) less?)))

(define (bisect a b less?)
  ;;
  ;; For a non-decreasing function  f(x) on an integer interval [a, b]
  ;; find  (the smallest?)   integer x such  that the  predicate less?
  ;; does not hold for f(x).  Think of finding a zero of monotonically
  ;; increasing  f(x) where  the predicate compares  the value  of the
  ;; function to zero. Note that the function "crosses zero" somewhere
  ;; between  x - 1  and x, so  that either of  two might serve  as an
  ;; approximation for  the root. Since the values  of a function f(x)
  ;; are used only as an  input to the predicate less? we redefine the
  ;; predicate to be a function of x instead.
  ;;
  (if (not (< a b))                     ; FIXME: sanity checks?
      a
      (let ((c (quotient (+ a b) 2)))
        (if (less? c)
            (bisect (+ 1 c) b less?)
            (bisect a c less?)))))

(debug-print
 (let* ((f (lambda (x) (* x x)))
        (less? (lambda (x) (< (f x) 100))))
   (bisect 0 100 less?)))

(define (list-range lst a b)
  (drop (take lst (+ 1 b)) a))

(debug-print
 (list-range '(0 1 2 3 4 5) 2 4))

(debug-print
 (mpts->npts work-function
             (map make-task '(50 50 50 50 50 60 60 60 60 60))
             (make-world 1000)))

;;;
;;; This one is called from  PG, expects integer problem sizes and the
;;; world size:
;;;
(define (qm-mpts->npts sizes workers)
  (mpts->npts work-function
              (map make-task sizes)
              (make-world workers)))

;;;
;;; This one  is called from PG, see  se_eigen_module.f90. Returns the
;;; shape of  the BALACS  grid to be  used by ScaLAPACK  in eigenvalue
;;; problems:
;;;
(define (grid-shape n)
  (if (< n 9)
      (cons 1 n)
      (let ((m (isqrt n)))
        (cons m m))))