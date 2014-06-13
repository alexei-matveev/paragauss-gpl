#
# ParaGauss, a program package for high-performance computations
# of molecular systems
# Copyright (C) 2014
# T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
# M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
# A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
# T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
# M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
# M. Roderus, N. Rösch
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 as published
# by the Free Software Foundation [1].
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# [1] http://www.gnu.org/licenses/gpl-2.0.html
#
# Please see the accompanying LICENSE file for further information.
#
#
# Make has no scopes, not to pollute the namespace
# use $(GUILE)/..., or GUILE- prefix for targets and $(GUILE-...)
# prefix for variables. A notable exception so far are the
# top-level tragets
#
#       $(libguile-comm.a)
#       $(guile-qm.o)
#

#
# We expect the before "include"ing this file the
# variable $(GUILE) is set to a suitable prefix.
#
# CURDIR is set on every make or $(MAKE) -C dir:
#
GUILE ?= $(CURDIR)

#
# These may be used to refer to the targets from outside:
#
libguile-comm.a = $(GUILE)/libguile-comm.a
guile-qm.o = $(GUILE)/guile-qm.o

#
# One implementation uses Fortran integers for communicators,
# another (more complicated) uses the C MPI_Comm wrapped into
# SMOB:
#
GUILE-libguile-comm-impl.o = $(GUILE)/libguile-comm-fint.o
#UILE-libguile-comm-impl.o = $(GUILE)/libguile-comm-smob.o

#
# libguile.so does not export many  macros which are part of the Guile
# API.   The file  guile-api.c implements  some of  them  as functions
# with guile_macro_ prefix.
#
GUILE-libguile-comm-objs = \
        $(GUILE)/libguile-comm.o \
        $(GUILE-libguile-comm-impl.o) \
        $(GUILE)/guile-api.o \
        $(GUILE)/pi.o \

GUILE-objs = $(GUILE-libguile-comm-objs) \
        $(guile-qm.o)

GUILE-fobjs = $(GUILE)/scm.o

$(libguile-comm.a): $(GUILE-libguile-comm-objs) $(GUILE-fobjs)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

GUILE-clean:
	rm -f $(libguile-comm.a)

.PHONY: GUILE-clean

#
# Below we modify "global" variables and prerequisites of
# top-level targets ...
#

#
# This (global) variable (cobjs) is used in Make.rules to build and
# include dependencies:
#
cobjs += $(GUILE-objs)
f90objs += $(GUILE-fobjs)

#
# This is also a top-level (global) target, extend the list
# of dependencies:
#
clean: GUILE-clean
