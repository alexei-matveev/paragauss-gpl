/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
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

  Copyright (c) 2011-2013 Alexei Matveev
*/

#include <libguile.h>
#include <mpi.h>
#include <assert.h>

/* Also declares wrappers for communicator smob: */
#include "libguile-comm.h"

/* Example parallel code: */
#include "pi.h"

#define MAX_BUF_LENGTH 512


/* size_t varies between 32/64 platforms, will be redefined later: */
static MPI_Datatype MPI_SIZE_T = MPI_UNSIGNED;

static SCM object_to_string (SCM obj);
static SCM string_to_object (SCM obj);

static char *to_byte_string (SCM obj, size_t *lenp);
static SCM from_byte_string (const char *buf, size_t len);


/* Set a name on a communicator: */
static SCM
guile_comm_set_name (SCM world, SCM name)
{
  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  /* Some communicators have names associated with them: */
  char cname[MPI_MAX_OBJECT_NAME];

  /* Does not null-terninate, returns the number of bytes
     necessary: */
  int len = scm_to_locale_stringbuf (name, cname, MPI_MAX_OBJECT_NAME);

  if (len + 1 > MPI_MAX_OBJECT_NAME) /* bytes ++ \0 dont fit */
    len = MPI_MAX_OBJECT_NAME - 1;

  cname[len] = '\0';

  int ierr = MPI_Comm_set_name (comm, cname);
  assert (MPI_SUCCESS==ierr);

  return scm_from_int (len);
}

static SCM
guile_comm_init (SCM args)      /* MPI_Init */
{
  int argc, i;
  char **argv;

  /* Count number of arguments: */
  argc = scm_to_int (scm_length (args));

  argv = malloc ((argc + 1) * sizeof (char *));

  argv[argc] = NULL;

  for (i = 0; i < argc; i++)
    {
      argv[i] = scm_to_locale_string (scm_car (args));
      args = scm_cdr (args);
    }

  int ierr = MPI_Init (&argc, &argv);
  assert (MPI_SUCCESS==ierr);

  /*
    FIXME:  In  fact  we  dont  know if  MPI_Init  replaced  the  argv
    completely and who is  responsible for freeing these resources. So
    we do not attempt to free them.
  */

  return scm_from_comm (MPI_COMM_WORLD);
}

static SCM
guile_comm_finalize (void)      /* MPI_Finalize */
{
  int ierr = MPI_Finalize ();
  assert (MPI_SUCCESS==ierr);

  return scm_from_int (ierr);
}

static SCM
guile_comm_rank (SCM world)     /* MPI_Comm_rank (world, ...) */
{
  int rank;

  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int ierr = MPI_Comm_rank (comm, &rank);
  assert (MPI_SUCCESS==ierr);

  return scm_from_int (rank);
}

static SCM
guile_comm_size (SCM world)     /* MPI_Comm_size (world, ...) */
{
  int size;

  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int ierr = MPI_Comm_size (comm, &size);
  assert (MPI_SUCCESS==ierr);

  return scm_from_int (size);
}

static SCM
guile_comm_barrier (SCM world)  /* MPI_Barrier (world, ...) */
{
  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int ierr = MPI_Barrier (comm);
  assert (MPI_SUCCESS==ierr);

  return scm_from_int (ierr);
}

/* MPI_Bcast, note argument order */
static SCM
guile_comm_bcast (const SCM world, const SCM root, const SCM obj)
{
  size_t len;
  char *sendbuf;
  char recvbuf[MAX_BUF_LENGTH];

  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int iroot = scm_to_int (root);

  int rank;
  int ierr = MPI_Comm_rank (comm, &rank);
  assert (MPI_SUCCESS==ierr);

  if (rank == iroot)
    {
      /* Searialize the object, dont forget to free() later: */
      sendbuf = to_byte_string (obj, &len);
      /* FIXME: recv buffer has finite length: */
      assert (len <= MAX_BUF_LENGTH);
    }

  /* Broadcast the size, or should we always send MAX_BUF_LENGTH? */
  ierr = MPI_Bcast (&len, 1, MPI_SIZE_T, iroot, comm);
  assert (MPI_SUCCESS==ierr);

  /* FIXME: recv buffer has finite length: */
  assert (len <= MAX_BUF_LENGTH);

  if (rank == iroot)
    {
      ierr = MPI_Bcast (sendbuf, len, MPI_CHAR, iroot, comm);
      assert (MPI_SUCCESS==ierr);
      free (sendbuf);
    }
  else
    {
      ierr = MPI_Bcast (recvbuf, len, MPI_CHAR, iroot, comm);
      assert (MPI_SUCCESS==ierr);
    }

  if (rank == iroot)
    return obj;          /* FIXME: should we return a copy instead? */
  else
    return from_byte_string (recvbuf, len);
}

/* MPI_Sendrecv, note argument order */
static SCM
guile_comm_send_recv (SCM world, SCM dst, SCM src, SCM tag, SCM obj)
{
  size_t len;
  char *sendbuf;
  char recvbuf[MAX_BUF_LENGTH];

  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int idst = scm_to_int (dst);
  int isrc = scm_to_int (src);

  /* FIXME: the same tag for send and recv: */
  int itag = scm_to_int (tag);

  /* Searialize the object, dont forget to free() later: */
  sendbuf = to_byte_string (obj, &len);
  assert (len <= MAX_BUF_LENGTH); /* here: <= */

  MPI_Status stat;

  /* Send just enough elements: */
  int ierr = MPI_Sendrecv (sendbuf,            len, MPI_CHAR, idst, itag, \
                           recvbuf, MAX_BUF_LENGTH, MPI_CHAR, isrc, itag, \
                           comm, &stat);
  assert (MPI_SUCCESS==ierr);

  free (sendbuf);

  /* Get the size of the received data: */
  int ilen;
  ierr = MPI_Get_count (&stat, MPI_CHAR, &ilen);
  assert (MPI_SUCCESS==ierr);
  assert (ilen <= MAX_BUF_LENGTH); /* redundant, as MPI would
                                      fail */
  return from_byte_string (recvbuf, ilen);
}


/*
  Send as  a procedure, receive  as a function that  returns arbitrary
  types unrelated to input is an ugly abstraction:
*/
static SCM
guile_comm_send (SCM world, SCM dst, SCM tag, SCM obj)
{
  size_t len;
  char *buf;

  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int idst = scm_to_int (dst);
  int itag = scm_to_int (tag);

  /* Searialize the object, dont forget to free() later: */
  buf = to_byte_string (obj, &len);
  assert (len < MAX_BUF_LENGTH);

  /* Send just enough elements: */
  int ierr = MPI_Send (buf, len, MPI_CHAR, idst, itag, comm);
  assert (MPI_SUCCESS==ierr);

  free (buf);

  return scm_from_int (ierr);
}

static SCM
guile_comm_recv (SCM world, SCM src, SCM tag)
{
  char buf[MAX_BUF_LENGTH];

  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int isrc = scm_to_int (src);
  int itag = scm_to_int (tag);

  MPI_Status stat;

  int ierr = MPI_Recv (buf, MAX_BUF_LENGTH, MPI_CHAR, isrc, itag, comm, &stat);
  assert (MPI_SUCCESS==ierr);

  /* Get the size of the received data: */
  int ilen;
  ierr = MPI_Get_count (&stat, MPI_CHAR, &ilen);
  assert (MPI_SUCCESS==ierr);

  return from_byte_string (buf, ilen);
}

/* MPI_Comm_split (world, color, ...) */
static SCM
guile_comm_split (SCM world, SCM color)
{
  int ierr;

  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  /* Color will define the coutries: */
  int icolor = scm_to_int (color);

  /* Key defines  the rank assignment within the  country, use world
     ranks for that: */
  int key;
  ierr = MPI_Comm_rank (comm, &key);
  assert (MPI_SUCCESS==ierr);

  MPI_Comm country;        /* Part of the world of the same color */

  ierr = MPI_Comm_split (comm, icolor, key, &country);
  assert (MPI_SUCCESS==ierr);

  return scm_from_comm (country);
}


/* I  am afraid we  cannot delegate  freeing communicators  to garbage
   collector. Do it explicitly: */
static SCM
guile_comm_free (SCM world)     /* MPI_Comm_free */
{
  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int ierr = MPI_Comm_free (&comm);
  assert (MPI_SUCCESS==ierr);

  return scm_from_int (ierr);
}


/* double guile_comm_pi (MPI_Comm world, int n); */
static SCM
guile_comm_pi (SCM world, SCM n)
{
  /* Extract MPI_Comm, verifies the type: */
  MPI_Comm comm = scm_to_comm (world);

  int N = scm_to_int (n);

  /* Compute PI in parallel: */
  double dbl = pi (comm, N);

  return scm_from_double (dbl);
}


/* These write/read scheme objects to scheme strings: */
static SCM
string_to_object (SCM str)
{
  SCM port = scm_open_input_string (str);

  SCM obj = scm_read (port);

  scm_close_port (port);

  return obj;
}

static SCM
object_to_string (SCM obj)      /* variant 2 */
{
  SCM str = scm_object_to_string (obj, SCM_UNDEFINED);

  return str;
}

/*
  These serialize/deserialize objects to char buffers.  Dont forget to
  free() the result of to_byte_string() when finished:
*/
static char *
to_byte_string (SCM obj, size_t *lenp)
{
  SCM str = object_to_string (obj);

  return scm_to_locale_stringn (str, lenp);
}

static SCM
from_byte_string (const char *buf, size_t len)
{
  SCM str = scm_from_locale_stringn (buf, len);

  return string_to_object (str);
}

#define EXPORT(name, req, opt, rst, func)               \
  (scm_c_define_gsubr (name, req, opt, rst, func),      \
   scm_c_export (name, NULL))

void
guile_comm_module_init (void *unused)
{
  (void) unused;

  /* size_t varies between 32/64 platforms, set MPI_SIZE_T here: */
  if (sizeof (size_t) == sizeof (unsigned long))
    MPI_SIZE_T = MPI_UNSIGNED_LONG;
  else if (sizeof (size_t) == sizeof (unsigned int))
    MPI_SIZE_T = MPI_UNSIGNED;
  else
    assert (0);
  guile_comm_smob_init();

  EXPORT ("comm-init", 1, 0, 0, guile_comm_init);
  EXPORT ("comm-finalize", 0, 0, 0, guile_comm_finalize);
  EXPORT ("comm-rank", 1, 0, 0, guile_comm_rank);
  EXPORT ("comm-size", 1, 0, 0, guile_comm_size);
  EXPORT ("comm-barrier", 1, 0, 0, guile_comm_barrier);

  EXPORT ("comm-bcast", 3, 0, 0, guile_comm_bcast);

  EXPORT ("comm-send-recv", 5, 0, 0, guile_comm_send_recv);
  EXPORT ("comm-send", 4, 0, 0, guile_comm_send);
  EXPORT ("comm-recv", 3, 0, 0, guile_comm_recv);

  EXPORT ("comm-split", 2, 0, 0, guile_comm_split);
  EXPORT ("comm-free", 1, 0, 0, guile_comm_free);
  EXPORT ("comm-set-name", 2, 0, 0, guile_comm_set_name);
  EXPORT ("comm-pi", 2, 0, 0, guile_comm_pi);
}
