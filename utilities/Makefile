# -*- makefile -*-
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
#  ### Makefile for utilites

#### Home Directory ####
ifndef HOME
HOME = /home/$(USER)
endif


#### Directory where executables are installed #####
INSTALLDIR = $(HOME)/bin

#### NAME OF EXECUTABLE ####
EXE_GX2XYZ       = gx2xyz
EXE_GX2XYZ_A     = gx2xyz_a
EXE_EPE2XYZ      = epe2xyz
EXE_XYZ2GX       = xyz2gx
EXE_EWALD        = ewald
EXE_EWALD_NEW    = ewald_new
EXE_GET_EPER     = get_epe_r
EXE_GENLAT       = genlat
EXE_MKZEO        = make_zeo
EXE_MKZEO_MM     = make_zeo_mm
EXE_TIP_FROM_XYZ = def_tip_from_xyz
EXE_TIP_WATERS   = def_tip_waters
EXE_GET_WATER_FF = get_water_forcefield
EXE_EPE2PC       = epe_env2ext_pc
EXE_CHECK_EPE    = check_epe
EXE_SOL_IN_WATER = solute_in_water


ALL_EXE = $(EXE_GX2XYZ) $(EXE_EPE2XYZ) $(EXE_PCR2XYZ) \
	  $(EXE_PCS2XYZ) $(EXE_XYZ2GX) $(EXE_EWALD) $(EXE_GET_EPER) \
	  $(EXE_GENLAT) $(EXE_MKZEO) $(EXE_MKZEO_MM) $(EXE_EWALD_NEW) \
	  $(EXE_GX2XYZ_A) $(EXE_TIP_FROM_XYZ) \
	  $(EXE_TIP_WATERS) $(EXE_GET_WATER_FF) \
	  $(EXE_EPE2PC) $(EXE_CHECK_EPE) $(EXE_SOL_IN_WATER)

#### COMPILER ####
FC = f90

#### COMPILER FLAGS ####
# These are not critical:
F90FLAGS = -g -O0
LINKFLAGS = 

#### LDFLAGS, LIBRARY-PATH ####
# Now only for Absoft
LIBS = -lU77 -lV77


all:	$(ALL_EXE)

%: %.f90
	$(FC) $(LINKFLAGS) $(LIBS) -o $(@) $(^)

##### SPECIAL COMMANDS #####
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -f $(ALL_EXE)

install:
	for f in $(ALL_EXE); do \
		if [ -f $f ] ; then \
			mv -i $f $(INSTALLDIR); \
		fi \
	done

uninstall:
	for f in $(ALL_EXE); do \
		rm -f $(INSTALLDIR)/$f; \
	done
