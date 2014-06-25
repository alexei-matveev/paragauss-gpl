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
MACH =	POWER4

#### Home Directory ####

HOME =	~/matveev

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

MPIDIR = /gpfs/apps/OPENMPI/1.4.3/64/smp
FC = $(MPIDIR)/bin/mpif90

F77 = $(FC)

ifeq ($(WITH_OLD_INPUT),0)
#(error NEW INPUT doesnt work on Regatta)
endif

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
#
# extraflags = -qinitauto=7FF7FFFF -qflttrap=inv:nanq:zero:enable -qsigtrap
F90BASEFLAGS = -qsuffix=f=F90  -q64 $(extraflags)
FBASEFLAGS   = -qsuffix=f=f -qfixed -q64 $(extraflags)
optflags     = -O3 -qstrict
optaltflags  = -O3 -qstrict
dbgflags     = $(optflags) -g
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
LINKFLAGS    = -q64 $(extraflags)
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS)
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS))

## C compiler
CPP = /lib/cpp
CPPOPTIONS = $(CPPUNDEF) -traditional-cpp -P -C -I$(BASEDIR)/include \
	-DFPP_AIX_XLF
SETMTIME = touch -r

FPP = $(CPP)
FPPOPTIONS += $(CPPOPTIONS) -D_$(MACH)
#
# type "make DEBUG=1" to compile in debug features
# or add it to $(debug_targets) in "misc.inc"
#
FPPOPTIONS += $(if $(DEBUG), -DFPP_DEBUG=$(DEBUG), )
FPPOPTIONS += -DFPP_PARAGAUSS_VERS="\"$(paragauss_vers)\""

CC = $(MPIDIR)/bin/mpicc

CCFLAGS = -q64
CCFLAGS += 

#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = ranlib

MPIINCLUDE =	
MPILIBS =	
LAPACKLIBS =	-L/gpfs/apps/LAPACK/lib64 -llapack -lessl # FIXME: $$BSC_LDFLAGS
BLASLIBS  =	
SYSTEMLIBS =	

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

# not comm/mpi_dir/vpp!
COMMCDIR =      comm/mpi_dir

# libraries, possibly classified by type:
LIBS       = $(LAPACKLIBS) $(BLASLIBS) $(COMMLIBS) $(SYSTEMLIBS)

ifeq ($(WITH_GUILE),1)
    INCLUDE += -I$(BASEDIR)/guile
    LIBS += -lguile
endif

