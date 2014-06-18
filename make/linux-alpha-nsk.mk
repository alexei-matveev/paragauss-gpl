#-*- makefile -*-
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
WITH_OLD_INPUT = 1
### Machine Keyword for FPP
MACH =	LINUX

#### Home Directory ####
HOME = /home/icct/shor

#### Directory for MPI Executable ####
MPIBIN = $(HOME)/exe

#### Directory where Sources are stored ####
ifndef BASEDIR
  BASEDIR = $(PWD)
endif

#### Directory where Scripts are stored #####
BINDIR =	$(BASEDIR)/bin

#### Executables are stored in:
EXEDIR = $(MPIBIN)

#### NAME OF EXECUTABLE ####
exe = mainscf_$(VERS)
EXE = $(EXEDIR)/$(exe)

#### COMPILER ####
## Fortran 90/77

MPIDIR = /common/mpich
FC  = $(MPIDIR)/bin/mpif90
F77 = $(FC)

#### Extensions compiler understands:
src_ext = F90

#### Does compiler UPPERCASE module files:
UPPERCASE.mod = 0

#### COMPILER FLAGS ####
## F90BASEFLAGS free  form Fortran 90 flags (always used)
## FBASEFLAGS   fixed form Fortran 90 flags (always used)
## F90FLAGS    additonal Fortran 90 flags
## F90ALTFLAGS additonal Fortran 90 alternate flags (see makro ALTFLAGS below to see for which files)
## CCFLAGS  C compiler flags

F90BASEFLAGS = -c -assume underscore -assume no2underscores -assume buffered_io
FBASEFLAGS   = -c -assume underscore -assume no2underscores -assume buffered_io
optflags     = -arch host -fast -tune host -O4
optaltflags  = -arch host  -tune host -O0
dbgflags     = $(optflags) -check bounds -check format -check overflow -g
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
staticflags  = -non_shared
LINKFLAGS    = $(if $(STATIC), $(staticflags), )
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS)
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS))

# this is a 32-bit version
ext = $(if $(STATIC),32s,32)
VERS := $(VERS)-$(ext)

## CPP/FPP preprocessor:

CPP = /usr/bin/cpp
CPPUNDEF = -undef
CPPOPTIONS = $(CPPUNDEF) -traditional -P -C -I$(BASEDIR)/include

FPP = $(CPP)
FPPOPTIONS = $(CPPOPTIONS) -D_$(MACH) -D_COMPAC_FORTRAN
#
# type "make DEBUG=1" to compile in debug features
# or add it to $(debug_targets) in "misc.inc"
#
FPPOPTIONS += $(if $(DEBUG), -DFPP_DEBUG=$(DEBUG), )
FPPOPTIONS += -DFPP_PARAGAUSS_VERS="\"$(paragauss_vers)\""
FPPOPTIONS += -DMAX_PATH=256


SETMTIME = touch -r


## C compiler:
CC = $(MPIDIR)/bin/mpicc
CCFLAGS = 


#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = ranlib

MPIINCLUDE =	/nowhere/mpi/include # mpif90 knows it!
#MPILIBS =	-lfmpich -lmpich 
MPILIBS =
LAPACKLIBS =	-lcxml
BLASLIBS =	
# NAG, Absoft provide them for system()/etime() calls:
SYSTEMLIBS =	$(MPILIBS) -lm
#PVMINCLUDE =	/usr/lib/pvm3/include
#PVMLIBS =	-lfpvm3 -lpvm3

# dynamic load balancing (DLB) module comes with several variants adapted to
# different needs, see dlb module README
# The default variant 1 uses RMA objects for storage purposes
DLB_VARIANT = 0

ifeq ($(WITH_SCALAPACK),1)
  BLACS_LIB = -lblacs
  SCALAPACK_LIB = -lscalapack
endif

COMMINCLUDE =	$(MPIINCLUDE)
COMMLIBS =	$(MPILIBS)
COMMDIR =	comm/mpi_dir

#
# These are, eventually, absolute paths:
#
DIRS = $(patsubst %,$(BASEDIR)/%,$(dirs))

# (note -p/-I syntax on different platforms/compilers)
# -I$(COMMINCLUDE) only needed for mpif.h/fpvm3.h,
# not needed with mpif90 script (as on linux):
INCLUDE = $(patsubst %,-I%,$(DIRS)) $(COMMINCLUDE)

# include path for C files:
CINCLUDE =	-I$(BASEDIR)/$(COMMDIR) $(COMMINCLUDE)

# The one with underscores:
COMMCDIR =	comm/mpi_dir

# libraries, possibly classified by type:
LIBS       = \
	$(SCALAPACK_LIB) $(BLACS_LIB) \
	$(LAPACKLIBS) $(BLASLIBS) \
	$(COMMLIBS) $(SYSTEMLIBS)

ifeq ($(WITH_GUILE),1)
    INCLUDE += -I$(BASEDIR)/guile
    LIBS += -lguile
endif

# FIXME: objects that require F90ALTFLAGS:
f90objs_altflags = modules/gradient_2c_fit_ch_module.o \
		modules/int_send_2cob3c_module.o \
		calc_3center.o \
		modules/symm_positions.o \
		modules/occupied_levels_module.o \
		modules/density_calc_module.o \
		modules/fitcontract_module.o \
		modules/gradient_data_module.o \
		optimizer/step_module.o
