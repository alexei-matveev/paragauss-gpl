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

all: example

BASEDIR=$(CURDIR)/..
include $(BASEDIR)/machine.inc

#
# Set it before include $(BASEDIR)/Make.rules:
#
f90objs = \
        $(BASEDIR)/modules/type_module.o \
        $(BASEDIR)/assert_failed.o \
        $(CURDIR)/error_handler.o \
        $(CURDIR)/example.o

AD=.
include module.mk

include $(BASEDIR)/Make.rules

$(AD)/example: $(example.o) $(assert_failed.o) $(error_handler.o) $(libad.a)
	$(FC) -o $(@) $(example.o) -L$(AD) -lad $(assert_failed.o) $(error_handler.o)

clean:
	rm -f $(f90objs) $(f90objs:.o=.mod)
	rm -f $(f90objs:.o=.m.dep) $(f90objs:.o=.o.dep)
	rm -f $(f90objs:.o=.F90)
	rm -f $(AD)/example
