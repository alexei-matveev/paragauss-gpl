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
#include <mpi.h>
#include <assert.h>
// #include <stdio.h>

#include "pi.h"

double pi (MPI_Comm world, int n)
{
    int rank, size, rc;

    rc = MPI_Comm_size (world, &size);
    assert(MPI_SUCCESS==rc);

    rc = MPI_Comm_rank (world, &rank);
    assert(MPI_SUCCESS==rc);

    double h = 1.0 / n;
    double s = 0.0;
    int i;
    for (i = rank; i < n; i += size) {
        double x = h * (i + 0.5);
        s += 4.0 / (1.0 + x * x);
    };

    double partial = s * h;
    double sum;
    rc = MPI_Allreduce (&partial, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
    assert(MPI_SUCCESS==rc);

    // printf("%d of %d pi = %f computed from %d terms\n", rank, size, sum, n);

    return sum;
}

