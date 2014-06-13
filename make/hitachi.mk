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
### Machine Keyword for FPP
MACH =	SR8000

CROSS_COMPILER = 1
ifdef CROSS_COMPILER
	usr = /usr/SR8000/USR
else
	usr = /usr
endif


#### Home Directory ####

ifndef HOME
 HOME =  /home/h/h0351ag
endif

#### Directory for MPI Executable ####
MPIBIN = $(HOME)/exe

EXEDIR = $(MPIBIN)

#### Directory where Sources are stored ####
ifndef BASEDIR
  BASEDIR = $(PWD)
endif

#### Directory where Scripts are stored #####
BINDIR = $(BASEDIR)/bin

#### NAME OF EXECUTABLE ####
exe = mainscf_$(VERS)
EXE = $(EXEDIR)/$(exe)

#### COMPILER ####
## Fortran 90/77

#FC  = $(usr)/ccs/bin/xf90
FC  = mpif90 
F77 = $(FC)

#### Extensions compiler understands:
src_ext = F90

#### COMPILER FLAGS ####
## F90BASEFLAGS free  form Fortran 90 flags (always used)
## FBASEFLAGS   fixed form Fortran 90 flags (always used)
## F90FLAGS    additonal Fortran 90 flags
## F90ALTFLAGS additonal Fortran 90 alternate flags (see makro ALTFLAGS below to see for which files)
## CCFLAGS  C compiler flags
#
# type "make NOOPT=1 modules/anything.compile"
# or add it to $(noopt_targets) in "misc.inc"
#
F90BASEFLAGS = -64 -model=F1 -Xpcomp -e -i,P -w -loglist -nohugeary
FBASEFLAGS   = -64 -model=F1 -Xpcomp -e -i,P -w -loglist -nohugeary
optflags     = -Os -pvdiag -parallel=0
optaltflags  = -opt=0 -parallel=0
dbgflags     = $(optflags) -g
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
staticflags  = 
LINKFLAGS    = $(if $(STATIC), $(staticflags), )
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS)
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS))

CCFLAGS = -64 -parallel=0
# Brehm: LINKFLAGS = -64 +SBTLB -lc
LINKFLAGS = -64 +SBTLB

## C compiler:

ifdef CROSS_COMPILER
	CPP = /lib/cpp
	CPPUNDEF = -undef
	CPPOPTIONS = $(CPPUNDEF) -traditional-cpp -P -C -I$(BASEDIR)/include
else
	CPP  = $(usr)/ccs/bin/cpp
	CPPUNDEF =  
	CPPOPTIONS = $(CPPUNDEF) -P -C -I$(BASEDIR)/include
endif

FPP = $(CPP)
FPPOPTIONS += $(CPPOPTIONS) -D_$(MACH)
#
# type "make DEBUG=1" to compile in debug features
# or add it to $(debug_targets) in "misc.inc"
#
FPPOPTIONS += $(if $(DEBUG), -DFPP_DEBUG=$(DEBUG), )
FPPOPTIONS += -DFPP_PARAGAUSS_VERS="\"$(paragauss_vers)\""
FPPOPTIONS += -DMAX_PATH=256

SETMTIME = touch -r

# CC = $(usr)/ccs/bin/cc
CC = mpicc

#### LDFLAGS, LIBRARY-PATH ####

AR = xar
RANLIB = xranlib

MPIINCLUDE = 
MPILIBS = 
# MPIINCLUDE = $(usr)/include -I$(usr)/mpi/include  
# MPILIBS = -L$(usr)/mpi/lib/lib64 -lfmpi -lmpi -lf90 -lf90c 

LAPACKLIBS = -L$(usr)/lib/LAPACK/lib64s -llapack
BLASLIBS  = -L$(usr)/lib/BLAS/lib64s -lblas
# SYSTEMLIBS = -L$(usr)/local/lib -llrz64p -llrz64s 
SYSTEMLIBS = 	-L/usr/SR8000/USR/mpi/lib/lib64 \
		-lfmpi \
		-lmpi \
		-lf90c \
		-lhf90pvmath \
		-lhf90math \
		-lf90 \
		-lcpvmath \
		-lhxb \
		-lc \
		-lmach

# dynamic load balancing (DLB) module comes with several variants adapted to
# different needs, see dlb module README
# The default variant 1 uses RMA objects for storage purposes
DLB_VARIANT = 0

ifeq ($(WITH_SCALAPACK),1)
  BLACS_LIB = -lblacs
  SCALAPACK_LIB = -lscalapack
endif

ifeq ($(serial),1)
COMMINCLUDE =
COMMLIBS =
COMMDIR =	comm/serial_dir
else
COMMINCLUDE =	-I./
COMMLIBS =	
COMMDIR =	comm/mpi_dir
endif

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
LIBS       = $(LAPACKLIBS) $(BLASLIBS) $(COMMLIBS) $(SYSTEMLIBS)

ifeq ($(WITH_GUILE),1)
    INCLUDE += -I$(BASEDIR)/guile
    LIBS += -lguile
endif

# FIXME: objects that require F90ALTFLAGS:
f90objs_altflags = ss_calculate_grads.o 
