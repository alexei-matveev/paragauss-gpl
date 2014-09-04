#-*- makefile -*-
#
# ParaGauss,  a program package  for high-performance  computations of
# molecular systems
#
# Copyright (C) 2014     T. Belling,     T. Grauschopf,     S. Krüger,
# F. Nörtemann, M. Staufer,  M. Mayer, V. A. Nasluzov, U. Birkenheuer,
# A. Hu, A. V. Matveev, A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman,
# D. I. Ganyushin,   T. Kerdcharoen,   A. Woiterski,  A. B. Gordienko,
# S. Majumder,     M. H. i Rotllant,     R. Ramakrishnan,    G. Dixit,
# A. Nikodem, T. Soini, M. Roderus, N. Rösch
#
# This program is free software; you can redistribute it and/or modify
# it under  the terms of the  GNU General Public License  version 2 as
# published by the Free Software Foundation [1].
#
# This program is distributed in the  hope that it will be useful, but
# WITHOUT  ANY   WARRANTY;  without  even  the   implied  warranty  of
# MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE. See  the GNU
# General Public License for more details.
#
# [1] http://www.gnu.org/licenses/gpl-2.0.html
#
# Please see the accompanying LICENSE file for further information.
#
### Machine Keyword for FPP
MACH =	LINUX

#### Home Directory ####
HOME = /home

#### Directory for MPI Executable ####
MPIBIN = $(HOME)/exe
BASEBIN = $(HOME)/bin

#### Directory for different libraries (force fields, basis sets(?))
LIBDIR = $(HOME)/lib

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

MPIDIR = /nowhere/
GCCDIR = /usr/bin
ifeq ($(serial),1)
  FC = $(GCCDIR)/gfortran
else
  FC = /usr/local/openmpi-1.4.2/bin/mpif90
endif
F77 = $(FC)

ifeq ($(WITH_OLD_INPUT),0)
$(error NEW INPUT doesnt yet work with Gfortran 4.3, set WITH_OLD_INPUT=1 in Makefile)
endif


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

#
# type "make NOOPT=1 modules/anything.compile"
# or add it to $(noopt_targets) in "misc.inc"
#
# -x f95 tells to not apply preprocessor to *.F90 files
F90BASEFLAGS = -x f95 -Wall #-fbounds-check #-pedantic-errors -Werror
FBASEFLAGS   = -x f95 -Wall #-fbounds-check #-pedantic-errors -Werror
optflags     = -O1 -g -std=gnu
optaltflags  = -O1 -g -std=gnu
dbgflags     = -O0 -g -std=gnu
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
extraflags   =
#staticflags  = -static
LINKFLAGS    = -static-libgfortran
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS)
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS)) $(extraflags)

## CPP/FPP preprocessor:
CPP = /usr/bin/cpp
CPPUNDEF = -undef
CPPOPTIONS = $(CPPUNDEF) -traditional -P -I$(BASEDIR)/include

FPP = $(CPP)
FPPOPTIONS += $(CPPOPTIONS) -D_$(MACH) 
#
# type "make DEBUG=1" to compile in debug features
# or add it to $(debug_targets) in "misc.inc"
#
FPPOPTIONS += $(if $(DEBUG), -DFPP_DEBUG=$(DEBUG), )
FPPOPTIONS += -DFPP_PARAGAUSS_VERS="\"$(paragauss_vers)\""
FPPOPTIONS += -DFPP_GFORTRAN_BUGS
# OpenMPI provides f90 interface, one can directly "use mpi":
FPPOPTIONS += -DUSE_MPI_MODULE
FPPOPTIONS += -DINTRINSIC_ISNAN
FPPOPTIONS += -DMAX_PATH=256
FPPOPTIONS += -DINTRINSIC_SECNDS -DINTRINSIC_ETIME


SETMTIME = touch -r


## C compiler:
ifeq ($(serial),1)
  CC = $(GCCDIR)/gcc
else
  CC = /usr/local/openmpi-1.4.2/bin/mpicc
endif
CCFLAGS = -Wall


#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = ranlib

LAPACKLIBS =	-llapack
BLASLIBS =	-lf77blas -lblas -latlas 

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
COMMDIR =       comm/serial_dir
else
COMMINCLUDE =   $(MPIINCLUDE)
COMMLIBS =      $(MPILIBS)
COMMDIR =       comm/mpi_dir
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
LIBS       = \
	$(SCALAPACK_LIB) $(BLACS_LIB) \
	$(LAPACKLIBS) $(BLASLIBS) \
	$(COMMLIBS) $(SYSTEMLIBS)

ifeq ($(WITH_GUILE),1)
    INCLUDE += -I$(BASEDIR)/guile
    LIBS += -lguile
endif

# FIXME: objects that require F90ALTFLAGS:
f90objs_altflags = $(libttfs_opt.a)
