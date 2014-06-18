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
/*
 * It is version of MPI C-subroutines using only one buffer for
 * receiving.  This version run without errors.  Using now.
 *
 * If you (like me) wonder why there is not a single call to MPI_Wait
 * in this file, this is what MPI 2.2 says in Section 3.7.4:
 *
 *   "If an MPI_TEST that completes a send is repeatedly called with
 *    the same arguments, and a matching receive has been started,
 *    then the call will eventually return flag = true, unless the
 *    receive is satisfied by another send."
 */

#include <stdio.h>
#include "mpi.h"
#include <malloc.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

/* Global variables and prototypes with scope bigger than this file */
#include "externdecs.h"

/* Contains the actual message tag. */
int msgtag_buf = 999;

/* Contains the sender of the actual message. */
int source_buf = 999;

/* Contains the target of the actual message. */
int target_buf = 999;


int comm_sending_host_buf () {
  /*
   * There is long story related to this gem, you can imagine the joy
   * of searching related bugs. Since those days the two values for
   * MPI tags were discouraged. Enforce that and set them free some
   * day (FXIME) ...
   *
   *   if (msgtag_buf==29) msgtag_buf=30;
   */
  assert (msgtag_buf != 29);
  assert (msgtag_buf != 30);
  return source_buf;
}


int comm_msgtag_buf () {
  return msgtag_buf;
}


int comm_target_buf () {
  return target_buf;
}

/*
 * The following variables are global for this file only.  As it is
 * the case for the variables in comm_variables.c, the appendix "buf"
 * separates the global from other name spaces.
 */

/* Number of processors. */
static int n_procs_buf;

/* MPI needs this. */
static MPI_Status *statuses_buf;

/* This is the array of communication buffers. */
static char **arr_buf;

/*
 * This are handles, n_procs_buf for each buffer, which allow to query
 * the status of a buffer (f.e. whether it is readily received or not)
 */
static MPI_Request **req_buf;

/*
 * Maximum number of used communication buffers.  The length is
 * defined in comm_variables.c.
 */
static int max_buf = 3000;

/*
 * Number of communication buffers which are kept allocated all the
 * time.
 */
static int static_buf=3;

/*
 * This pointer references the request handles associated with the
 * buffer to which pointer_buf points.
 */
static MPI_Request *pointer_req_buf;


void comm_init_buffer_data (MPI_Fint *fworld) {
  /*
   * Initialization routine for the C communication layer, called once
   * from comm_enroll(). The first thing it does is converting the
   * *fworld fortran comunicator to MPI_Comm and sets the global
   * "comm_world".  Stores the number of processors and 1-based index
   * as used by fortran in file-global variables.  Furthermore it
   * allocates above defined arrays.
   */
  int ierr;
  int i, j;
  int n_procs, rank;

  // convert the fortran world communicator to MPI_Comm,
  // and place it into global variable used everywhere
  // in the C-layer:
  comm_world = MPI_Comm_f2c (*fworld);

  assert (comm_world != MPI_COMM_NULL);
  // this does not hold if an agent is running:
  // assert(comm_world == MPI_COMM_WORLD);

  // get the process count:
  ierr = MPI_Comm_size (comm_world, &n_procs);
  assert (ierr == MPI_SUCCESS);

  // get the process rank:
  ierr = MPI_Comm_rank (comm_world, &rank);
  assert (ierr == MPI_SUCCESS);

  // set file globals, caching, one might think of it:
  n_procs_buf = n_procs;

  // this is also assumend by comm_enroll():
  my_index_buf = rank + 1;

  statuses_buf = (MPI_Status*) malloc (n_procs_buf * sizeof (MPI_Status));
  if (!statuses_buf)
    errout ("statuses_buf couldn't be allocated!\n");

  arr_buf = (char**) malloc (max_buf * sizeof (char*));
  if (!arr_buf)
    errout ("arr_buf couldn't be allocated!\n");

  recv_buf = (char*) malloc (length_buf * sizeof (char));
  if (!recv_buf)
    errout ("recv_buf couldn't be allocated!\n");

  req_buf = (MPI_Request**) malloc (max_buf * sizeof (MPI_Request*));
  if (!req_buf)
    errout("req_buf couldn't be allocated!\n");

  /* Allocate only the static part of the communication buffers. */
  for (i = 0; i < static_buf; i++) {
    arr_buf[i] = (char*) malloc (length_buf * sizeof (char));
    if (!arr_buf[i])
      errout ("arr_buf[%d] couldn't be allocated!\n");//,i);

    req_buf[i] = (MPI_Request*) malloc (n_procs_buf * sizeof (MPI_Request));
    if (!req_buf[i])
      errout ("req_buf[%d] couldn't be allocated!\n");//,i);

    /*
     * By setting the handles to the value of MPI_REQUEST_NULL it is
     * indicated that they are currently not in use.
     */
    for (j = 0; j < n_procs_buf; j++)
      req_buf[i][j] = MPI_REQUEST_NULL;
  }
  /*
   * Explicitly setting the pointers to 0 to initialize but not
   * allocate them.
   */
  for (i = static_buf; i < max_buf; i++) {
    arr_buf[i] = 0;
    req_buf[i] = 0;
  }
}


static void comm_dealloc_unused_buffers () {
  /*
   * Frees buffers currently not in use.  Is called from init_send and
   * recv.
   */

  int i, ierr, ready = 0;

  for (i = static_buf; i < max_buf; i++) {
    /* Is buffer allocated? */
    if (arr_buf[i] != 0) {
      /* Is an isend pending on this buffer? */
      ierr = MPI_Testall (n_procs_buf, req_buf[i], &ready, statuses_buf);
      assert (ierr == MPI_SUCCESS);

      if (ready) {
	/* Being currently not necessary, this buffer can be freed. */
	free (arr_buf[i]);
	arr_buf[i] = 0;

	free ((char*) req_buf[i]);
	req_buf[i] = 0;
      }
    }
  }
}


static void comm_next_free_buffer () {
  /*
   * Sets pointer_buf to next free buffer, pointer_req_buf to array of
   * request handles associated with this buffer.  Calles as well from
   * send as from receive routines.  NOT ANYMORE: Returns index of the
   * found free buffer. This return value is needed for the receive
   * queue.
   */

  int i, j, ierr, ready=0;

  /*
   * Search a free buffer in an infinite loop.  A timeout should be
   * implemented here one day.
   */
  for (;;) {
    /* Loop over all buffers. */
    for (i = 0; i < max_buf; i++) {
      /*
       * Same procedure as above: First check whether buffer is
       * allocated, then check pending handles
       */
      if (arr_buf[i] != 0) {
	ierr = MPI_Testall (n_procs_buf, req_buf[i], &ready, statuses_buf);
	assert (ierr == MPI_SUCCESS);

	if (ready) {
	  /* Buffer is free for reuse. Setting pointers accordingly. */
	  pointer_buf = arr_buf[i];
	  pointer_req_buf = req_buf[i];
	  /* NOT ANYMORE: Return index of found buffer. */
	  return;
	}
      } else {
	/* Buffer is not allocated, so it is clear that it is free. */
	arr_buf[i] = (char*) malloc (length_buf * sizeof (char));
	if (!arr_buf[i])
	  errout ("arr_buf[%d] couldn't be allocated!\n");//,i);

	req_buf[i] = (MPI_Request*) malloc (n_procs_buf * sizeof (MPI_Request));
	if (!req_buf[i])
	  errout ("req_buf[%d] couldn't be allocated!\n");//,i);

	/* Initialize request handles. */
	for (j = 0; j < n_procs_buf; j++)
	  req_buf[i][j] = MPI_REQUEST_NULL;

	/* Buffer is ready for use. Setting pointers accordingly. */
	pointer_buf = arr_buf[i];
	pointer_req_buf = req_buf[i];
	/* NOT ANYMORE: Return index of found buffer. */
	return;
      }
    }
  }
}

void comm_init_send_buf (int *target, int *msgtag) {
  /*
   * Sets msgtag and target for internal sends and receives.
   * Furthermore sets counter to zero and pointer_buf to next free
   * buffer.
   */
  msgtag_buf = *msgtag;
  target_buf = *target;
  count_buf = 0;
  comm_next_free_buffer ();
}

static void comm_send_host_buf (int *index, int *msgtag) {
  /*
   * Sends actual buffer and sets pointer_buf to next free buffer This
   * version is used of opaque and transparent sends as well.
   */
  int ierr;
  ierr = MPI_Isend (pointer_buf, count_buf, MPI_PACKED, *index - 1, *msgtag,
		    comm_world, &pointer_req_buf[*index - 1]);
  assert (ierr == MPI_SUCCESS);
  count_buf = 0;
  comm_next_free_buffer ();
}


static void comm_mcast_others_buf (int *msgtag) {
  /*
   * Trivial multicast. Sends to all other nodes one after
   * another. Opaque and transparent.
   */
  int ierr, index;
  for (index = 1; index <= n_procs_buf; index++) {
    if (index != my_index_buf) {
      ierr = MPI_Isend (pointer_buf, count_buf, MPI_PACKED, index - 1, *msgtag,
			comm_world, &pointer_req_buf[index - 1]);
      assert (ierr == MPI_SUCCESS);
    }
  }
  count_buf = 0;
  comm_next_free_buffer ();
}


void comm_send_buf () {
  /* Actually sends the buffer */
  if (target_buf < 0) {
    comm_mcast_others_buf (&msgtag_buf);
  }
  else {
    comm_send_host_buf (&target_buf, &msgtag_buf);
  }
}

void comm_send () {
  comm_send_buf ();
  comm_dealloc_unused_buffers ();
}

void comm_save_recv_buf () {
  /*
   * This is a misnomer, I guess it should have been
   * comm_safe_recv_buf. "Safe" is there because this function tried
   * to abort the execution when detecting an error tag. By now this
   * functionality is largely unused: use MPI_Abort instead.
   */

  int ierr, time, ready=0;
  int timeout = 9999999, wait_time = 100;
  int timelast;
  MPI_Status status;
  int rank; // for debug prints only

  ierr = MPI_Comm_rank (comm_world, &rank); // for debug prints only
  assert (ierr == MPI_SUCCESS);

  timelast = 0;
  for (time = 0; time < timeout; time++) {
    ierr = MPI_Iprobe (source_buf - 1, msgtag_buf, comm_world, &ready, &status);
    assert (ierr == MPI_SUCCESS);

    if (ready) {
      if ((source_buf - 1) == MPI_ANY_SOURCE)
	source_buf = status.MPI_SOURCE + 1;

      if (msgtag_buf == MPI_ANY_TAG)
	msgtag_buf = status.MPI_TAG;

      if (msgtag_buf == 255)
	errout ("Error went from slaves. Stop 1\n");

      /* Before this was implemented as MPI_Irecv/MPI_Wait */
      ierr = MPI_Recv (recv_buf, length_buf, MPI_PACKED,
		       source_buf - 1, msgtag_buf, comm_world,
		       &status);
      assert (ierr == MPI_SUCCESS);

      count_buf = 0;
      return;
    }

    if ( (time - timelast) > 50000 ){
      timelast = time;
      // let us use base-1 worker indices in this debug print:
      printf ("WARNING: Worker %d waiting for tag %d from %d\n already %d sec!\n",
	      rank + 1, msgtag_buf, source_buf, time / 10000);
    }

    ierr = MPI_Iprobe (MPI_ANY_SOURCE, 255, comm_world, &ready, &status);
    assert (ierr == MPI_SUCCESS);

    if (ready) errout ("Error went from slaves. Stop 2\n");
    /* SGI:    ierr = usleep(wait_time); */
    usleep (wait_time);
  }
  printf ("Could not successfully receive message %d from host %d\n", msgtag_buf, source_buf);
  errout ("Could not successfully receive! Timeout?\n");
}

void comm_save_recv_c (int *index, int *msgtag, int *info) {
  /*
   * This is a part of the Fortran interface. Misnomer "save" ==
   * "safe".
   */

  if (*index < 0) {
    source_buf = MPI_ANY_SOURCE + 1;
  }
  else {
    source_buf = *index;
  }

  if (*msgtag < 0) {
    msgtag_buf = MPI_ANY_TAG;
  }
  else {
    msgtag_buf = *msgtag;
  }

  comm_dealloc_unused_buffers ();
  comm_save_recv_buf ();

  // FIXME: real error handling, maybe?
  *info = 0;
}

int comm_save_recv_nonblocking_buf (int *index, int *msgtag) {
  int ierr, ready = 0;
  MPI_Status status;

  if (*index < 0) {
    source_buf = MPI_ANY_SOURCE + 1;
  } else {
    source_buf = *index;
  }

  if (*msgtag < 0) {
    msgtag_buf = MPI_ANY_TAG;
  } else {
    msgtag_buf = *msgtag;
  }
  comm_dealloc_unused_buffers ();

  ierr = MPI_Iprobe (source_buf - 1, msgtag_buf, comm_world, &ready, &status);
  assert (ierr == MPI_SUCCESS);

  if (ready) {
    if ((source_buf - 1) == MPI_ANY_SOURCE)
      source_buf = status.MPI_SOURCE + 1;

    if (msgtag_buf == MPI_ANY_TAG)
      msgtag_buf = status.MPI_TAG;

    /* Before this was implemented as MPI_Irecv/MPI_Wait */
    ierr = MPI_Recv (recv_buf, length_buf, MPI_PACKED,
		     source_buf - 1, msgtag_buf, comm_world,
		     &status);
    assert (ierr == MPI_SUCCESS);

    count_buf = 0;
    return 1;
  }
  return 0;
}
