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

/* There is no corresponding header  file with prototypes as these are
   only used  from Fortran. In C  use directly the  macros provided by
   libguile.h.  The prototypes are  encoded in the interface blocks in
   scm.f90 though and have to be consistent with these definitions: */
int guile_macro_scm_to_int (SCM obj)
{
  return scm_to_int(obj);
}

int (scm_is_true) (SCM obj)
{
  return scm_is_true(obj);
}

int (scm_is_symbol) (SCM obj)
{
  return scm_is_symbol(obj);
}

int (scm_is_null) (SCM obj)
{
  return scm_is_null(obj);
}

/* FIXME: these return constants: */
SCM scm_eol ()
{
  return SCM_EOL;
}

SCM scm_undefined ()
{
  return SCM_UNDEFINED;
}

/*
  Guile  API function  scm_c_export() works  with a  variable argument
  list which is not callable from Fortran by standard means:
*/
void
scm_c_export_1 (const char *name)
{
  scm_c_export (name, NULL);
}
