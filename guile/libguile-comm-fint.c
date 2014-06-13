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
*/
#include <libguile.h>
#include <mpi.h>

#include "libguile-comm.h"

//
// Here we use the one-to-one relation between MPI_Comm instances
// and their integer Fortran representations of type MPI_Fint
//

//
// This converts an MPI_Comm to a SCM int:
//
SCM scm_from_comm (const MPI_Comm comm) // comm-init wants to return MPI_COMM_WORLD
{
    return scm_from_int (MPI_Comm_c2f (comm));
}

//
// This converts a SCM int to an MPI_Comm:
//
MPI_Comm scm_to_comm (const SCM comm)
{
    return MPI_Comm_f2c (scm_to_int (comm));
}

void guile_comm_smob_init (void)
{
    // nothing
}
