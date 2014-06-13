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
#include <mpi.h>
#include <assert.h>

/* Here we assume that  communicators are implemented in fortran style
   as   SCM   integers,   libguile-comm   should  be   consistent   to
   interoperate: */
#include "libguile-comm.h"

#ifdef WITH_BGY3D
#  include "bgy3d-guile.h"      /* bgy3d_guile_init */
#endif

/*
 * FIXME: when we find a way to make Intel compiler return type(scm_t)
 * aka  struct  with  an  intptr_t  the  wrappers  for  these  Fortran
 * functions  may  be  removed.  Run  ParaGauss by  calling  these  in
 * sequence:
 */
int qm_init (void);
void qm_run (int world);
void qm_finalize (int world);

static SCM
guile_qm_init (void)
{
    int world = qm_init ();
    return scm_from_int (world);
}

static SCM
guile_qm_run (SCM world)
{
    int fworld = scm_to_int (world);
    qm_run (fworld);
    return SCM_UNSPECIFIED;
}

static SCM
guile_qm_finalize (const SCM world)
{
    int fworld = scm_to_int (world);
    qm_finalize (fworld);
    return SCM_UNSPECIFIED;
}

void qm_init_scheme (void);

static SCM
guile_paragauss_module_init (void)
{
    scm_c_define_gsubr ("qm-init", 0, 0, 0, guile_qm_init);
    scm_c_define_gsubr ("qm-run", 1, 0, 0, guile_qm_run);
    scm_c_define_gsubr ("qm-finalize", 1, 0, 0, guile_qm_finalize);

    /* See ../modules/paragauss.f90: */
    qm_init_scheme ();

    return SCM_UNSPECIFIED;
}

static void
guile_main (void *data, int argc, char **argv)
{
  /*
    Note  that  the  names defined  here  (not  as  a part  of  module
    definiitons) are  put into  the private namespace  of (guile-user)
    module. If you  want to call any of these you  may need to "steal"
    it  from there  by dereferencing  them as  e.g.   (@@ (guile-user)
    guile-paragauss-module-init).

    Calling this  will define comm-*  gsubrs in (guile  comm internal)
    module:
   */
    scm_c_define_module ("guile comm internal", guile_comm_module_init, NULL);

  /*
   * Calling this  will define  a few qm-*  gsubrs defined  in Fortran
   * sources:
   */
    scm_c_define_gsubr ("guile-paragauss-module-init", 0, 0, 0, guile_paragauss_module_init);

#ifdef WITH_BGY3D
    /* The function bgy3d_guile_init() initialzes Petsc, defines a few
       gsubrs with  bgy3d-prefix and returns. The tricky  part is that
       it registers  an atexit() function  that calls PetscFinalize().
       This function both, uses MPI and invokes MPI_Finalize(), if MPI
       was initialized from  PetscInitialize(). If you call MPI_Init()
       earlier than this  point, you seem to have  to register another
       atexit() handler that does MPI_Finalize(). */
    bgy3d_guile_init (argc, argv);
#endif

    scm_shell (argc, argv); // never returns
}

int
main (int argc, char **argv)
{
  /* void scm_boot_guile (int argc, char **argv, void (*main_func)
     (void *data, int argc, char **argv), void *data) */
  scm_boot_guile (argc, argv, guile_main, 0);
  return 0; /* never reached */
}
