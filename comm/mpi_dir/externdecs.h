/*
  ParaGauss,  a program package  for high-performance  computations of
  molecular systems

  Copyright (C) 2014     T. Belling,     T. Grauschopf,     S. Krüger,
  F. Nörtemann, M. Staufer,  M. Mayer, V. A. Nasluzov, U. Birkenheuer,
  A. Hu, A. V. Matveev, A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman,
  D. I. Ganyushin,   T. Kerdcharoen,   A. Woiterski,  A. B. Gordienko,
  S. Majumder,     M. H. i Rotllant,     R. Ramakrishnan,    G. Dixit,
  A. Nikodem, T. Soini, M. Roderus, N. Rösch

  This program is free software; you can redistribute it and/or modify
  it under  the terms of the  GNU General Public License  version 2 as
  published by the Free Software Foundation [1].

  This program is distributed in the  hope that it will be useful, but
  WITHOUT  ANY   WARRANTY;  without  even  the   implied  warranty  of
  MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE. See  the GNU
  General Public License for more details.

  [1] http://www.gnu.org/licenses/gpl-2.0.html

  Please see the accompanying LICENSE file for further information.
*/
#include "mpi.h"

/*
 * These are only the extern declarations of the variables defined in
 * comm_variables.c. More documentation is available there.
 */

extern int my_index_buf;

extern int length_buf;
extern int count_buf;
extern char *pointer_buf;
extern char *recv_buf;

extern MPI_Comm comm_world;

/*
 * A few prototypes to make compiler happy.  These functions are also
 * called from mpipack.c, so export their prototypes:
 */
void errout (char *message);
void comm_send_buf ();
void comm_save_recv_buf ();
