#!@bindir@/guile-qm
;;; -*- mode: scheme; -*- vim: set syntax=scheme:
!#
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
;;;
;;; Make sure  the shebang line  and %load-path are set  correctly. In
;;; particular "bindir"  and "pkgdatadir" between  @-markers should be
;;; replaced by, say, a path to PG repository.
;;;
;;; Load  user-specific config, this  may also  set the  %load-path to
;;; search for files/modules.  This should work with both 1.8 and 2.0.
;;; Use  absolute paths  when  extending %load-path  ---  there is  no
;;; globbing for ~/ at compile time.
;;;
(cond-expand
 ((not guile-2) (use-modules (ice-9 syncase))) ; eval-when for 1.8
 (else))                                ; nothing

(eval-when
 (eval load compile)      ; extend %load-path at compile time with 2.0
 (set! %load-path (cons "@pkgdatadir@" %load-path))
 (let ((conf (string-append (getenv "HOME") "/.qmrc")))
   (if (file-exists? conf)
       (load conf))))

;;;
;;; Home-grown modules need %load-path to be extended e.g. as in
;;;
;;; (set! %load-path (cons "~/darcs/ttfs-mac" %load-path))
;;;
;;; It may be  convenient to do this in ~/.qmrc,  see the load command
;;; above.
;;;
;;; One of  the procedures, qm-trace-hook, if defined  or imported, is
;;; called from scheme_trace_hook() in PG.  So import it if you want a
;;; trace.   Note that PG  code may  and probably  should look  up the
;;; names directly  in the guile modules. An  example is qm-mpts->npts
;;; from (guile scheduling) which is imported in schedeig-code.
;;;
(use-modules ((guile paragauss)
              #:select
              (qm-main)))

;;;
;;; This expression is a list of energies, one per input:
;;;
(qm-main (command-line)) ; first argument, the program name, is
                         ; handled by getopt
