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
include ../../../machine.inc

COMPFLAGS =
LINKFLAGS =  $(SCLALIB) $(BLCSLIB) $(LAPACKLIBS)


all : rrec

clean : 
	(rm *.o *.mod)

rrec : options_module.o time_module.o runrec_module.o runtimerecord.o
	${FC} -o rrec options_module.o time_module.o runrec_module.o runtimerecord.o ${LINKFLAGS}
	
runtimerecord.o : runtimerecord.f90
	${FC} ${COMPFLAGS} -c -o runtimerecord.o runtimerecord.f90
	
runrec_module.o : runrec_module.f90
	${FC} ${COMPFLAGS} -c -o runrec_module.o runrec_module.f90
	
options_module.o : options_module.f90
	${FC} ${COMPFLAGS} -c -o options_module.o options_module.f90

time_module.o : time_module.f90
	${FC} ${COMPFLAGS} -c -o time_module.o time_module.f90
