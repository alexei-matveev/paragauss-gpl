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
HOME = /users/ttfs

#### Directory for MPI Executable ####
MPIBIN = $(HOME)/exe/mpich
MPIBIN = /XEON/exe

#### Directory where Sources are stored ####
ifndef BASEDIR
  BASEDIR = $(PWD)
endif

#### Directory where Scripts are stored #####
BINDIR =	$(BASEDIR)/bin

#### Executables are stored in:
EXEDIR = $(MPIBIN)

#### NAME OF EXECUTABLE ####
#EXE = $(EXEDIR)/mainscf_$(VERS)
exe = mainscf_$(VERS)
EXE = $(EXEDIR)/$(exe)


#### COMPILER ####
## Fortran 90/77

MPIDIR = /opt/mpich/ch-p4
MPIDIR = /usr/local/mpich-1.2.2.3
FC  = $(MPIDIR)/bin/mpif90
F77 = $(FC)

#### Extensions compiler understands:
src_ext = F90

#### Does compiler UPPERCASE module files:
UPPERCASE.mod = 1

#### COMPILER FLAGS ####
## F90BASEFLAGS free  form Fortran 90 flags (always used)
## FBASEFLAGS   fixed form Fortran 90 flags (always used)
## F90FLAGS    additonal Fortran 90 flags
## F90ALTFLAGS additonal Fortran 90 alternate flags (see makro ALTFLAGS below to see for which files)
## CCFLAGS  C compiler flags

#
# type "make NOOPT=1 modules/anything.compile"
#
# -B80 : print name of the subroutine when entering
# -et  : backtrace stack
#  -Rb -Rc : chack array boundaries
ifdef NOOPT
  F90BASEFLAGS =  -YEXT_NAMES=LCS -YEXT_SFX=_
  FBASEFLAGS   =  -YEXT_NAMES=LCS -YEXT_SFX=_
  F90FLAGS     = -g -et #-Rb -Rc #-B80
  F90ALTFLAGS  = -g -et #-Rb -Rc #-B80
  F77FLAGS     = -YDEALLOC=ALL -YEXT_NAMES=LCS -YEXT_SFX=_ -g
  LINKFLAGS    =  -m4 -YEXT_NAMES=LCS -YEXT_SFX=_ -g
else
  F90BASEFLAGS =  -YEXT_NAMES=LCS -YEXT_SFX=_
  FBASEFLAGS   =  -YEXT_NAMES=LCS -YEXT_SFX=_
  F90FLAGS     =  -O -B100
  F90ALTFLAGS  =  -O -B100
  F77FLAGS     =  -O -B100  -YDEALLOC=ALL -YEXT_NAMES=LCS -YEXT_SFX=_
  LINKFLAGS    =  -YEXT_NAMES=LCS -YEXT_SFX=_
endif

## CPP/FPP preprocessor:
CPP = /lib/cpp
CPPUNDEF = -undef
CPPOPTIONS = $(CPPUNDEF) -traditional-cpp -P -C -I$(BASEDIR)/include

FPP = $(CPP)
#
# type "make DEBUG=1" to compile in debug features
#
FPPOPTIONS = $(CPPOPTIONS) -D_$(MACH)
ifdef DEBUG
	FPPOPTIONS += -DFPP_DEBUG=$(DEBUG)
endif
#FPPOPTIONS += -DFPP_FAST_COMPILE
FPPOPTIONS += -DFPP_PARAGAUSS_VERS=\"$(VERS)\"
FPPOPTIONS += -DFPP_O3_ABSOFT


SETMTIME = touch -r


## C compiler:
CC = $(MPIDIR)/bin/mpicc
CCFLAGS = 


#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = ranlib

MPIINCLUDE =	/nowhere/mpi/include # mpif90 knows it!
MPILIBS =	-lmpich
LAPACKLIBS =	-L$(HOME)/atlas/LAPACK -llapack1
BLASLIBS =	-L$(HOME)/atlas/ATLAS/lib/Linux_PIII -lf77blas -latlas
# NAG, Absoft provide them for system()/etime() calls:
SYSTEMLIBS =	$(MPILIBS) -lU77 -lV77 
PVMINCLUDE =	/usr/lib/pvm3/include
PVMLIBS =	-lfpvm3 -lpvm3

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
INCLUDE = $(patsubst %,-p%,$(DIRS)) $(COMMINCLUDE)

# include path for C files:
CINCLUDE =	-I$(BASEDIR)/$(COMMDIR) -I$(COMMINCLUDE)

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
f90objs_altflags = modules/disp_rep_module.o

