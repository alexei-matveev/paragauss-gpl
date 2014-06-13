/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  ParaGauss, a program package for high-performance computations
  of molecular systems
  Copyright (C) 2014
  T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
  M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
  A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
  T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
  M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
  M. Roderus, N. Rösch

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License version 2 as published
  by the Free Software Foundation [1].

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  [1] http://www.gnu.org/licenses/gpl-2.0.html

  Please see the accompanying LICENSE file for further information.

  Copyright (c) 2011-2013 Alexei Matveev
*/

#include <stdlib.h>
#include <libguile.h>
#include <mpi.h>                /* used in libguile-comm.h */
#include "libguile-comm.h"

static void
inner_main (void *data, int argc, char **argv)
{
  (void) data;
  /*
    Module setup alone does not lead to MPI initialization. Instead it
    defines comm-init procedure for the user to issue

        (comm-init (command-line))
        ...
        (comm-finalize)

    in Scheme code.

    Note  that  the  names  defined  here  not as  part  of  a  module
    definitions (currently  none) were put into  the private namespace
    of (guile-user) module.  If you want to call any  of these you may
    need to "steal"  it from there by dereferencing  them as e.g.  (@@
    (guile-user) gsubr).

    Calling this  will define comm-*  gsubrs in (guile  comm internal)
    module:
  */
  scm_c_define_module ("guile comm internal", guile_comm_module_init, NULL);

  scm_shell (argc, argv);       /* never returns */
}

int
main (int argc, char **argv)
{
  scm_boot_guile (argc, argv, inner_main, NULL);
  return 0; /* never reached */
}
