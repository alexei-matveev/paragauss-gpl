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
WITH_NANINFCHK = 0
### Machine Keyword for FPP
MACH =	LINUX

#### Home Directory ####
ifndef HOME
HOME = /home/icct/shor
endif

#### Directory for different libraries (force fields, basis sets(?))
LIBDIR = $(HOME)/lib

# dynamic load balancing (DLB) module comes with several variants adapted to
# different needs, see dlb module README
# The default variant 1 uses RMA objects for storage purposes
DLB_VARIANT = 0

ifeq ($(WITH_SCALAPACK),1)
  BLACS_LIB = -lblacs
  SCALAPACK_LIB = -lscalapack
endif

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

MPIDIR =
FC  = mpiifort
F77 = $(FC)

#### Extensions compiler understands:
src_ext = F90

#### Does compiler UPPERCASE module files:
UPPERCASE.mod = 0

# Command to compare MOD-files,
# Intel seems to put "random" number
# byte positions ~49:52 (decimal)
# so skip 52 bytes from comparison:
CMPMOD = cmp -s -i 52

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
# -B80 : print name of the subroutine when entering
# -et  : backtrace stack
#  -Rb -Rc : chack array boundaries
F90BASEFLAGS = 
FBASEFLAGS   = 
#optflags     = -O -mcmodel=large -no-global-hoist -auto -heap-arrays 1000
optflags     = -O1 -ftz -align all -traceback -heap-arrays 1000
optaltflags  = -O1 -ftz -align all -traceback -heap-arrays 1000
dbgflags     = $(optflags) -g
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
LINKFLAGS    = 
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS)
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS))

## CPP/FPP preprocessor:
CPP = /usr/bin/cpp
CPPUNDEF = -undef
CPPOPTIONS = $(CPPUNDEF) -traditional-cpp -P -C -I$(BASEDIR)/include

FPP = $(CPP)
FPPOPTIONS += $(CPPOPTIONS) -D_$(MACH) -D_NKSG6_NSK
#
# type "make DEBUG=1" to compile in debug features
# or add it to $(debug_targets) in "misc.inc"
#
FPPOPTIONS += $(if $(DEBUG), -DFPP_DEBUG=$(DEBUG), )
FPPOPTIONS += -DFPP_PARAGAUSS_VERS="\"$(paragauss_vers)\""
FPPOPTIONS += -DMAX_PATH=256
FPPOPTIONS += -DUSE_MPI_MODULE
FPPOPTIONS += -DFPP_NO_BIG_AUTOMATIC_ARRAYS


SETMTIME = touch -r


## C compiler:
CC = mpiicc
CCFLAGS = -DF77_EXT_NAMES=lowercase_


#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = ranlib

MPILIBS =
MKLPATH =	$(MKLROOT)/lib/intel64
MKLINCLUDE =	$(MKLROOT)/include

LAPACKLIBS =    -L$(MKLPATH) -I$(MKLINCLUDE) -I$(MKLINCLUDE)/intel64/lp64 \
		-lmkl_lapack95_lp64 -lmkl_blas95_lp64 \
		-Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a \
		$(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm
BLASLIBS =	
# NAG, Absoft provide them for system()/etime() calls:
SYSTEMLIBS =	$(MPILIBS)

ifeq ($(serial),1)
COMMINCLUDE =
COMMLIBS =
COMMDIR =	comm/serial_dir
else
COMMINCLUDE =	$(MPIINCLUDE)
COMMLIBS =	$(MPILIBS)
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
LIBS       = \
	$(SCALAPACK_LIB) $(BLACS_LIB) \
	$(LAPACKLIBS) $(BLASLIBS) \
	$(COMMLIBS) $(SYSTEMLIBS)

ifeq ($(WITH_GUILE),1)
    INCLUDE += -I$(BASEDIR)/guile
    LIBS += -lguile
endif

# FIXME: objects that require F90ALTFLAGS:
f90objs_altflags = modules/matrix_types.o \
	modules/matrix_methods.o \
	modules/reltrafo.o \
	modules/relgrads.o \
	modules/relgrads_store.o
