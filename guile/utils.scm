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
;;; Utilities that do not require any Guile extension.
;;;
(define-module (guile utils)
  #:use-module (guile compat)           ; for 1.8
  #:export (derivatives
            list-derivatives
            dfridr
            fmap
            ddd
            qtrap
            qsimp
            memoize
            bohr->angstrom
            angstrom->bohr
            isqrt
            ;; Syntax/macros:
            begin0
            let-env
            ))

(use-modules (ice-9 pretty-print))


;;;
;;; Some  codes  operate in  angstroms,  others  use  atomic units  by
;;; convention:
;;;
(define (bohr->angstrom x) (* x 0.52917706))
(define (angstrom->bohr x) (/ x 0.52917706))


;;;
;;; Let-over-lambda  here.   Given a  function  (f  x)  The result  of
;;; (memoize f)  is a  function (f' x)  == (f  x) that will  cache all
;;; results.  FIXME: make it work for arbitrary number of arguments.
;;;
(define (memoize f)
  (let ((*cache* '()))                  ; empty cache
    (lambda (x)
      (cdr (or (assoc x *cache*)        ; first check the cache
               (let* ((y (f x))         ; otherwise invoke f
                      (p (cons x y)))   ; make dictionary pair
                 (set! *cache* (cons p *cache*)) ; cache new pair
                 p))))))                         ; and return it too

;;;
;;; See  qtrap(),  trapzd()  in  Numerical  Recepies,  Integration  of
;;; functions. Almost literal translation:
;;;
(define (trapzd f a b s n)
  (if (= n 1)
      (* 0.5 (- b a) (+ (f a) (f b)))   ; (b - a) * (f(a) + f(b)) / 2
      (let* ((it (expt 2 (- n 2)))      ; it = 2^(n-2)
             (del (/ (- b a) it)))      ; del = (b - a) / it
        (let loop ((sum 0.0)            ; sum f(a + del * [i + 1/2]), 0 <= i < it
                   (i 0))
          (if (< i it)
              (loop (+ sum (f (+ a (* del (+ i 0.5)))))
                    (+ 1 i))
              (* 0.5 (+ s (* del sum)))))))) ; (s + del * sum) / 2

;;;
;;; This function, (qtrap f a b), integrates f from a to b:
;;;
(define (qtrap f a b)
  "(qtrap f a b) integrates f from a to b."
  (let ((JMAX 20)                       ; max iteration
        (EPS 1.0e-8))                   ; relative accuracy
    (let loop ((ost 0.0)                ; old trapzd result
               (j 1))                   ; iteration
      ;; (pretty-print (list 'qtrap: j))
      (let ((st (trapzd f a b ost j)))
        (if (or (> j JMAX)              ; prevent infinite recursion
                (and (<= (abs (- st ost)) (* EPS (abs ost)))
                     (> j 5)))          ; prevent early return
            st
            (loop st (+ 1 j)))))))


;;;
;;; This function, (qsimp f a b), integrates f from a to b:
;;;
(define (qsimp f a b)
  "(qsimp f a b) integrates f from a to b."
  (let ((JMAX 20)                       ; max iteration
        (EPS 1.0e-8))                   ; relative accuracy
    (let loop ((ost 0.0)                ; old trapzd result
               (os 0.0)                 ; old Simpson result
               (j 1))                   ; iteration
      ;; (pretty-print (list 'qsimp: j))
      (let* ((st (trapzd f a b ost j))
             (s (/ (- (* 4.0 st) ost) 3.0))) ; (4 * st - ost) / 3
        (if (or (> j JMAX)              ; prevent infinite recursion
                (and (<= (abs (- s os)) (* EPS (abs os)))
                     (> j 5)))          ; prevent early return
            s
            (loop st s (+ 1 j)))))))

;;;
;;; (pretty-print (- (sin 1.0) (qtrap cos 0.0 1.0)))
;;; =>
;;; 1.04490815999725e-9
;;;
;;; (pretty-print (- (sin 10.0) (qtrap cos 0.0 10.0)))
;;; =>
;;; -1.05554243123862e-9
;;;
;;; (pretty-print (- (sin 1.0) (qsimp cos 0.0 1.0)))
;;; =>
;;; -2.78650214013965e-10
;;;
;;; Here for the integration limit of 1.0 "qsimp" iterates to j = 7 to
;;; get 1.0e-8 realtive accuracy, whereas "qtrap" iterates to j = 14.
;;;

(define (dfridr func x h)
  "(dfridr f x h) -> (fprime . err)

Returns the  derivative of a function  `f' at a point  `x' by Ridders'
method  of polynomial  extrapolation.  The  value `h'  is input  as an
estimated initial stepsize; it need not be small, but rather should be
an increment in `x' over which `f' changes substantially.  An estimate
of the error in the derivative is returned as `err'.

Parameters: Stepsize  is decreased by  CON=1.4 at each  iteration. Max
size of tableau is set by NTAB=10. Return when error is SAFE=2.0 worse
than the best so far.

Translated from Numerical Recipes."
  ;;
  ;; Central difference approximation for derivative at x:
  ;;
  (define (fp h)
    (/ (- (func (+ x h))
          (func (- x h)))
       (* 2.0 h)))
  ;;
  (let* ((CON 1.4)
         (CON2 (* CON CON))
         (NTAB 10)
         (SAFE 2.0)
         (a (make-array 0.0 NTAB NTAB))             ; A[NTAB, NTAB]
         (a! (lambda (f i j) (array-set! a f i j))) ; A[i, j] := f
         (a (lambda (i j) (array-ref a i j)))       ; A[i, j]
         (break #f)
         (res #f)         ; result will be updated here
         (err (inf)))     ; very big number, initially
    ;;
    ;; Initial approximation for the derivative:
    ;;
    (a! (fp h) 0 0)
    (do ((i 1 (+ i 1))
         (h h (/ h CON)))
        ((or (<= NTAB i) break)         ; do while i < NTAB ...
         (cons res err))    ; ... return a pair of result and an error
      ;; Successive columns in the Neville tableau will go to smaller
      ;; stepsizes and higher orders of extrapolation.  Try new,
      ;; smaller stepsize:
      (a! (fp h) 0 i)                   ; A[0, i] := fp(h)
      (do ((j 1 (+ j 1))
           (fac CON2 (* fac CON2)))
          ((<= i j))                    ; do while j < i
        ;; Compute  extrapolations of  various orders, requiring  no new
        ;; function evaluations:
        (let* ((f1 (a (- j 1) i))       ; f1 := A[j-1, i]
               (f2 (a (- j 1) (- i 1))) ; f2 := A[j-1, i-1]
               (f3 (/ (- (* fac f1) f2) ; f3 := (fac * f1 - f2) / (fac - 1)
                      (- fac 1.0))))
          (a! f3 j i))
        ;; The error strategy is to compare each new extrapolation to
        ;; one order lower, both at the present stepsize and the
        ;; previous one:
        (let* ((norm0 (abs (- (a j i) (a (- j 1) i))))
               (norm1 (abs (- (a j i) (a (- j 1) (- i 1)))))
               (errt (max norm0 norm1)))
          (if (<= errt err)
              (begin
                (set! err errt)
                (set! res (a j i))))))
      ;;
      ;; If higher order is worse by a significant factor SAFE, then
      ;; signal an early quit:
      ;;
      (set! break (>= (abs (- (a i i)
                              (a (- i 1) (- i 1))))
                      (* SAFE err))))))

;; (dfridr sin 1.57 0.1)
;; =>
;; (7.96326710740725e-4 . 8.87148427636264e-15)
;;
;; (cos 1.57)
;; =>
;; 7.96326710733263e-4

;;;
;;; Traverse the nested structure x and apply f to every element.
;;;
(define (fmap f x)
  (cond
   ((null? x) x)
   ((pair? x) (cons (fmap f (car x)) (fmap f (cdr x))))
   (else (f x))))


;;;
;;; Given a procedure (d f x) to compute derivatives of a univariate f
;;; at x use it to  compute partial derivatives of the multivariate (f
;;; x)  at  x  where  x  is  represented by  a  nested  list  of  real
;;; numbers. Here  "dd" stays for "do derivatives".   FIXME: what kind
;;; of operation is this?
;;;
(define (dd d f x)
  (let go ((f f) (x x))
    (cond
     ((null? x) x)
     ((pair? x) (let ((h (car x))
                      (t (cdr x)))
                  (cons (go (lambda (x) (f (cons x t))) h)
                        (go (lambda (x) (f (cons h x))) t))))
     (else (d f x)))))

;;;
;;;
;;;
(define (ddd f x)
  (define (d f x)
    (car (dfridr f x 0.1)))        ; numerical derivative of f(x) at x
  ;;
  (dd d f x))

;; (let ((x 1.0)
;;       (f sin))
;;   (pretty-print (ddd f x)))

;; (let ((x 2.0)
;;       (f cos))
;;   (pretty-print (ddd f x)))

;; (let ((x '(1.0))
;;       (f (lambda (x) (sin (car x)))))
;;   (pretty-print (ddd f x)))

;; (let ((x '((1.0)))
;;       (f (lambda (x) (sin (caar x)))))
;;   (pretty-print (ddd f x)))


;; (let ((x '(1.0 2.0))
;;       (f (lambda (x) (+ (sin (car x)) (cos (cadr x))))))
;;   (pretty-print (ddd f x)))

;; (let ((x '(1.0 (2.0)))
;;       (f (lambda (x) (+ (sin (car x)) (cos (caadr x))))))
;;   (pretty-print (ddd f x)))


(define (derivatives f args h)
  "(derivatives f x h) ->
                       list of df / dx
                                      i
For scalar valued function `f' differentiates numerically with respect
to all components of `x'."
  (let loop ((f f)
             (args args)
             (acc '()))                 ; accumulator for derivatives
    (if (null? args)
        (reverse acc)
        (let* ((x (car args))           ; first arg ...
               (args (cdr args))        ; ... and the rest of them
               (fy (lambda (y) (apply f (+ x y) args))) ; function of the first arg only
               (f1 (car (dfridr fy 0.0 h)))) ; df/dx, the error is discarded
          ;;
          ;; Continue iterations with a function of less arguments:
          ;;
          (loop (lambda args (apply f x args))
                args
                (cons f1 acc))))))

;; (define (f x y) (+ x (* 10 y)))
;; (derivatives f '(1.57 3.0) 0.1)
;; =>
;; (0.999999999999961 9.99999999999999)

(define (list-derivatives f args h)
  "(list-derivatives f x h) ->
                            list of lists with df  / dx
                                                 i     j
For  list valued function  `f' differentiates  numerically all  of its
components at point `x'."
  (let loop ((f f)
             (y (apply f args)) ; to get the shape of the result, apply f here
             (acc '()))
    (if (null? y)
        (reverse acc)
        (let ((g (lambda args (car (apply f args)))) ; first component of f
              (f (lambda args (cdr (apply f args))))) ; other components
          (loop f
                (cdr y)
                (cons (derivatives g args h) acc)))))) ; FIXME: list-derivatives?

;; (define (sincos x) (list (sin x) (cos x)))
;;
;; (list-derivatives sincos '(1.0) 0.1)
;; =>
;; ((0.540302305868135) (-0.841470984807894))

;;
;; Newton  method for  integer  square root.  Shamelessly copied  from
;; http://programmingpraxis.com/contents/standard-prelude/
;;
(define (isqrt n)
  "(isqrt n) => integer square root"
  (if (not (and (positive? n) (integer? n)))
      (error "must be positive integer" n)
      (let loop ((x n))
        (let ((y (quotient (+ x (quotient n x)) 2)))
          (if (< y x) (loop y) x)))))

;;;
;;; This form  evaluates all sub-expressions in order  but returns the
;;; value(s) of the first one:
;;;
;;;   (begin0
;;;     (sin 1)
;;;     (clean-up))
;;;
;;;   => sin(1):
;;;
(define-syntax-rule (begin0 exp exp* ...)
  (call-with-values
      (lambda () exp)
    (lambda vals
      exp* ...
      (apply values vals))))

;;;
;;; The first attempt on macros:
;;;
;;; (getenv "asdf")
;;; => "999"
;;;
;;; (let-env (("asdf" "123"))
;;;   (getenv "asdf"))
;;; => "123"
;;;
;;; (getenv "asdf")
;;; => "999"
;;;
(define-syntax-rule (let-env ((var val) ...)
                      exp exp* ...)
  (let ((old-env (environ)))
    (setenv var val) ...
    (begin0
      (begin exp exp* ...)
      (environ old-env))))


