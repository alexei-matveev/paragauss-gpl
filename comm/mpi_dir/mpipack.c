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
#include <stdio.h>
#include <assert.h>
#include "mpi.h"
#include "externdecs.h"

// static int c1_size_buf = 1;
static int i4_size_buf = 4;
static int r8_size_buf = 8;

void mpix_pkdouble_vec (double *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  istr = *stride;
  act_item = 0;
  suiting_items = (length_buf - count_buf) / r8_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr = MPI_Pack (&p[act_item], suiting_items, MPI_DOUBLE,
			 pointer_buf, length_buf, &count_buf,
			 comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_send_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    ierr = MPI_Pack (&p[act_item], rem_items, MPI_DOUBLE, pointer_buf,
		     length_buf, &count_buf, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	  ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_DOUBLE,
			   pointer_buf, length_buf, &count_buf,
			   comm_world);
	  assert (ierr == MPI_SUCCESS);
	}
      }
      comm_send_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
      ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_DOUBLE, pointer_buf,
		       length_buf, &count_buf, comm_world);
      assert (ierr == MPI_SUCCESS);
    }
  }
  *ierror = 0;
}

void mpix_upkdouble_vec (double *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  act_item = 0;
  istr = *stride;
  suiting_items = (length_buf - count_buf) / r8_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item], suiting_items, MPI_DOUBLE,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item], rem_items, MPI_DOUBLE, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	}
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item + i * istr], 1, MPI_DOUBLE,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item + i * istr], 1, MPI_DOUBLE,
		       comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  *ierror = 0;
}

void mpix_pkint_vec (int *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  istr = *stride;
  act_item = 0;
  suiting_items = (length_buf - count_buf) / i4_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr=MPI_Pack (&p[act_item], suiting_items, MPI_INT,
		       pointer_buf, length_buf, &count_buf,
		       comm_world);
      }
      comm_send_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    ierr = MPI_Pack (&p[act_item], rem_items, MPI_INT, pointer_buf,
		     length_buf, &count_buf, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	  ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_INT,
			   pointer_buf, length_buf, &count_buf,
			   comm_world);
	  assert (ierr == MPI_SUCCESS);
	}
      }
      comm_send_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
      ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_INT, pointer_buf,
		       length_buf, &count_buf, comm_world);
      assert (ierr == MPI_SUCCESS);
    }
  }
  *ierror = 0;
}

void mpix_upkint_vec (int *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  act_item = 0;
  istr = *stride;
  suiting_items = (length_buf - count_buf) / i4_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item], suiting_items, MPI_INT,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item], rem_items, MPI_INT, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	}
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item + i * istr], 1, MPI_INT,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item + i * istr], 1, MPI_INT,
		       comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  *ierror = 0;
}

void mpix_pkdouble_scalar (double *p, int *ierror) {
  int ierr;
  if (length_buf - count_buf < r8_size_buf)
    comm_send_buf ();
  ierr = MPI_Pack (p, 1, MPI_DOUBLE, pointer_buf,
		   length_buf, &count_buf, comm_world);
  assert (ierr == MPI_SUCCESS);
  *ierror = 0;
}

void mpix_upkdouble_scalar (double *p, int *ierror) {
  int ierr;
  if (length_buf - count_buf < r8_size_buf)
    comm_save_recv_buf ();
  ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		     p, 1, MPI_DOUBLE, comm_world);
  assert (ierr == MPI_SUCCESS);
  *ierror = 0;
}

void mpix_pkint_scalar (int *p, int *ierror) {
  int ierr;
  if (length_buf - count_buf < i4_size_buf)
    comm_send_buf ();
  ierr = MPI_Pack (p, 1, MPI_INT, pointer_buf,
		   length_buf, &count_buf, comm_world);
  assert (ierr == MPI_SUCCESS);
  *ierror = 0;
}

void mpix_upkint_scalar (int *p, int *ierror) {
  int ierr;
  if (length_buf - count_buf < i4_size_buf)
    comm_save_recv_buf ();
  ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		     p, 1, MPI_INT, comm_world);
  assert (ierr == MPI_SUCCESS);
  *ierror = 0;
}

void mpix_pkdouble_vecsc (double *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  istr = *stride;
  act_item = 0;
  suiting_items = (length_buf - count_buf) / r8_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr = MPI_Pack (&p[act_item], suiting_items, MPI_DOUBLE,
			 pointer_buf, length_buf, &count_buf,
			 comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_send_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    ierr = MPI_Pack (&p[act_item], rem_items, MPI_DOUBLE, pointer_buf,
		     length_buf, &count_buf, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	  ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_DOUBLE,
			   pointer_buf, length_buf, &count_buf,
			   comm_world);
	  assert (ierr == MPI_SUCCESS);
	}
      }
      comm_send_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
      ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_DOUBLE, pointer_buf,
		       length_buf, &count_buf, comm_world);
      assert (ierr == MPI_SUCCESS);
    }
  }
  *ierror = 0;
}

void mpix_upkdouble_vecsc (double *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  act_item = 0;
  istr = *stride;
  suiting_items = (length_buf - count_buf) / r8_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item], suiting_items, MPI_DOUBLE,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item], rem_items, MPI_DOUBLE, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	}
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item + i * istr], 1, MPI_DOUBLE,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / r8_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item + i * istr], 1, MPI_DOUBLE,
		       comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  *ierror = 0;
}

void mpix_pkint_vecsc (int *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  istr = *stride;
  act_item = 0;
  suiting_items = (length_buf - count_buf) / i4_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr = MPI_Pack (&p[act_item], suiting_items, MPI_INT,
			 pointer_buf, length_buf, &count_buf,
			 comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_send_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    ierr = MPI_Pack (&p[act_item], rem_items, MPI_INT, pointer_buf,
		     length_buf, &count_buf, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	  ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_INT,
			   pointer_buf, length_buf, &count_buf,
			   comm_world);
	  assert (ierr == MPI_SUCCESS);
	}
      }
      comm_send_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
      ierr = MPI_Pack (&p[act_item + i * istr], 1, MPI_INT, pointer_buf,
		       length_buf, &count_buf, comm_world);
      assert (ierr == MPI_SUCCESS);
    }
  }
  *ierror = 0;
}

void mpix_upkint_vecsc (int *p, int *nitem, int *stride, int *ierror) {
  int ierr, i, rem_items, act_item, suiting_items, istr;
  rem_items = *nitem;
  act_item = 0;
  istr = *stride;
  suiting_items = (length_buf - count_buf) / i4_size_buf;
  if (istr == 1) {
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item], suiting_items, MPI_INT,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item], rem_items, MPI_INT, comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  else {
    rem_items = *nitem / istr;
    if (rem_items*istr < *nitem) rem_items++ ;
    while (suiting_items < rem_items) {
      if (suiting_items > 0) {
	for (i = 0; i < suiting_items; i++) {
	}
	ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
			   &p[act_item + i * istr], 1, MPI_INT,
			   comm_world);
	assert (ierr == MPI_SUCCESS);
      }
      comm_save_recv_buf ();
      act_item += suiting_items*istr;
      rem_items -= suiting_items;
      suiting_items = length_buf / i4_size_buf;
    }
    for (i = 0; i < rem_items; i++) {
    }
    ierr = MPI_Unpack (recv_buf, length_buf, &count_buf,
		       &p[act_item + i * istr], 1, MPI_INT,
		       comm_world);
    assert (ierr == MPI_SUCCESS);
  }
  *ierror = 0;
}
