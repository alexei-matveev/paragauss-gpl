#include "mpi.h"
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

/* All variables in this file are needed */
/* in more than one file. Therefore they are made */
/* global. To distinguish them from other variables, */
/* these and the variables in the file comm.c */
/* have an appended string "buf".  */

/* All variables in this "module" are declared */
/* extern in header file externdecs.h. This is */
/* necessary to get external linkage for them, i.e */
/* to make them globally available. */

/* Contains index of processor. */
int my_index_buf=999;

/* Length of send and receive buffers. */
int length_buf=40960;
/* Counter which gives position in actually used buffer. */
int count_buf=0;
/* Pointer to the buffer used by isend/pack und irecv. */
char *pointer_buf;
/* Buffer used by unpack */
char *recv_buf;

/*
 * This is a replacment for literal MPI_COMM_WORLD used originally
 * everywhere. As a preparation for the general case where
 *
 *   comm_world != MPI_COMM_WORLD
 *
 * introduce this variable. Initialze it by an invalid communicator,
 * will be initialized by comm_init_buffer_data():
 */
MPI_Comm comm_world = MPI_COMM_NULL;
