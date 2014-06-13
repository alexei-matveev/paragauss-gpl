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
/******************************************************************
*                                                                 *
*  Purpose: The standard PVM packing and unpacking C functions    *
*           are called by envelope functions of the same name.    *
*                                                                 *
*  This code is platform and compiler dependent !!!               *
*                                                                 *
*  Called by: module commpack_module.f90                           *
*                                                                 *
*  Author: TB                                                     *
*  Date 06/28/95                                                  *
*                                                                 *
******************************************************************/

#include "pvm3.h"
#include "f77names.h"

void F77_pkbyte_scalar(p,info)
char *p;
int *info;
{ *info =  pvm_pkbyte(p,1,1); }

void F77_pkstring(p,length,info)
char *p;
int *length, *info;
{ *info =  pvm_pkbyte(p,*length,1); }

void F77_pkcplx_scalar(p,info)
float *p;
int *info;
{ *info =  pvm_pkcplx(p,1,1); }

void F77_pkdcplx_scalar(p,info)
double *p;
int *info;
{ *info =  pvm_pkdcplx(p,1,1); }

void F77_pkdouble_scalar(p,info)
double *p;
int *info;
{ *info =  pvm_pkdouble(p,1,1); }

void F77_pkfloat_scalar(p,info)
float *p;
int *info;
{ *info =  pvm_pkfloat(p,1,1); }

void F77_pkint_scalar(p,info)
int *p;
int *info;
{ *info =  pvm_pkint(p,1,1); }

void F77_pkshort_scalar(p,info)
short *p;
int *info;
{ *info =  pvm_pkshort(p,1,1); }

void F77_upkbyte_scalar(p,info)
char *p;
int *info;
{ *info =  pvm_upkbyte(p,1,1); }

void F77_upkstring(p,length,info)
char *p;
int *length, *info;
{ *info =  pvm_upkbyte(p,*length,1); }

void F77_upkcplx_scalar(p,info)
float *p;
int *info;
{ *info =  pvm_upkcplx(p,1,1); }

void F77_upkdcplx_scalar(p,info)
double *p;
int *info;
{ *info =  pvm_upkdcplx(p,1,1); }

void F77_upkdouble_scalar(p,info)
double *p;
int *info;
{ *info =  pvm_upkdouble(p,1,1); }

void F77_upkfloat_scalar(p,info)
float *p;
int *info;
{ *info =  pvm_upkfloat(p,1,1); }

void F77_upkint_scalar(p,info)
int *p;
int *info;
{ *info =  pvm_upkint(p,1,1); }

void F77_upkshort_scalar(p,info)
short *p;
int *info;
{ *info =  pvm_upkshort(p,1,1); }

void F77_pkbyte_vec(p,nitem,stride,info)
char *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkbyte(p,*nitem,*stride); }

void F77_pkcplx_vec(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkcplx(p,*nitem,*stride); }

void F77_pkdcplx_vec(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkdcplx(p,*nitem,*stride); }

void F77_pkdouble_vec(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkdouble(p,*nitem,*stride); }

void F77_pkfloat_vec(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkfloat(p,*nitem,*stride); }

void F77_pkint_vec(p,nitem,stride,info)
int *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkint(p,*nitem,*stride); }

void F77_pkshort_vec(p,nitem,stride,info)
short *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkshort(p,*nitem,*stride); }

void F77_upkbyte_vec(p,nitem,stride,info)
char *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkbyte(p,*nitem,*stride); }

void F77_upkcplx_vec(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkcplx(p,*nitem,*stride); }

void F77_upkdcplx_vec(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkdcplx(p,*nitem,*stride); }

void F77_upkdouble_vec(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkdouble(p,*nitem,*stride); }

void F77_upkfloat_vec(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkfloat(p,*nitem,*stride); }

void F77_upkint_vec(p,nitem,stride,info)
int *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkint(p,*nitem,*stride); }

void F77_upkshort_vec(p,nitem,stride,info)
short *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkshort(p,*nitem,*stride); }

void F77_pkbyte_vecsc(p,nitem,stride,info)
char *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkbyte(p,*nitem,*stride); }

void F77_pkcplx_vecsc(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkcplx(p,*nitem,*stride); }

void F77_pkdcplx_vecsc(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkdcplx(p,*nitem,*stride); }

void F77_pkdouble_vecsc(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkdouble(p,*nitem,*stride); }

void F77_pkfloat_vecsc(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkfloat(p,*nitem,*stride); }

void F77_pkint_vecsc(p,nitem,stride,info)
int *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkint(p,*nitem,*stride); }

void F77_pkshort_vecsc(p,nitem,stride,info)
short *p;
int *nitem, *stride, *info;
{ *info =  pvm_pkshort(p,*nitem,*stride); }

void F77_upkbyte_vecsc(p,nitem,stride,info)
char *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkbyte(p,*nitem,*stride); }

void F77_upkcplx_vecsc(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkcplx(p,*nitem,*stride); }

void F77_upkdcplx_vecsc(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkdcplx(p,*nitem,*stride); }

void F77_upkdouble_vecsc(p,nitem,stride,info)
double *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkdouble(p,*nitem,*stride); }

void F77_upkfloat_vecsc(p,nitem,stride,info)
float *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkfloat(p,*nitem,*stride); }

void F77_upkint_vecsc(p,nitem,stride,info)
int *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkint(p,*nitem,*stride); }

void F77_upkshort_vecsc(p,nitem,stride,info)
short *p;
int *nitem, *stride, *info;
{ *info =  pvm_upkshort(p,*nitem,*stride); }
