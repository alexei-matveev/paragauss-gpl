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
MACH =	LINUX

#### Directory for MPI Executable ####
MPIBIN = $(HOME)/exe/openmpi

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

#
# If this  is included  too late in  the scope it  will spoil/redefine
# some of the make variables. In fact  all we need is to link to PETSC
# in  this case. $(PETSC_LIB)  is the  only interesting  variable.  PG
# sources do not  use that beast of a  library. To define $(PETSC_DIR)
# on Debian derivatives do
#
#   export PETSC_DIR=/usr/lib/petsc
#
# On  Debian Lenny  the  location  of the  include  file setting  MAKE
# variables is $(PETSC_DIR)/bmake/common/variables
#
ifeq ($(WITH_BGY3D),1)
	include $(PETSC_DIR)/conf/variables
endif

#### COMPILER ####
## Fortran 90/77

MPIDIR = /nowhere/
GCCDIR = /usr/lib/gcc-snapshot/bin
GCCDIR = /usr/bin
ifeq ($(serial),1)
  FC = $(GCCDIR)/gfortran
else
  FC = mpif90
endif
F77 = $(FC)

ifeq ($(WITH_OLD_INPUT),0)
$(error NEW INPUT doesnt yet work with Gfortran 4.3, set WITH_OLD_INPUT=1 in Makefile)
endif

#### Extensions compiler understands:
src_ext = F90

#### Does compiler UPPERCASE module files:
UPPERCASE.mod = 0

# Command to compare MOD-files,
# GNU FORTRAN puts a date in the first line (approx 70 chars)
CMPMOD = cmp -s -i 100

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
F90BASEFLAGS = -x f95 -Wall -freg-struct-return #-fbounds-check -finit-logical=false -finit-real=nan -ffpe-trap=invalid,zero #-pedantic-errors -Werror
FBASEFLAGS   = -x f95 -Wall -freg-struct-return #-fbounds-check -finit-logical=false -finit-real=nan -ffpe-trap=invalid,zero #-pedantic-errors -Werror
# Polyhedron flags: -march=native -ffast-math -funroll-loops -O3
optflags     = -O3 -g
optaltflags  = -O3 -g
dbgflags     = -O3 -g
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
extraflags   =
staticflags  = -static
LINKFLAGS    = $(if $(STATIC), $(staticflags), )
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS)
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS)) $(extraflags)

#modules/filename_module.o: extraflags=-std=gnu

## CPP/FPP preprocessor:
#
# -undef: Do not predefine any system-specific or GCC-specific macros.
#   The standard predefined macros remain defined. If you pass that to
#   GFortran then __FILE__  and __LINE__ will not be  defined.  Not so
#   with CPP. See if omitting it breaks anything.
#
# -P: some compilers do not  understand #line directives. But then the
#   location of the error is not precise. Omit it when possible.
#
# -traditional-cpp:  I think  it is  here  for a  relaxed matching  of
#   quotes in Fortran comments.
#
# -C: must  be there, otherwise everything  after string concatenation
#   operator // gets discarded.
#
# -ffreestanding:  Disable  stdc-predef.h  preinclude.  Otherwise  the
#   preprocessed files will contain C-sources such as /**/-comments.
#
CPP = /lib/cpp
CPPOPTIONS = -traditional-cpp -C -I$(BASEDIR)/include -ffreestanding

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
FPPOPTIONS += -DFPP_HAVE_ISNAN
FPPOPTIONS += -DINTRINSIC_ISNAN
FPPOPTIONS += -DMAX_PATH=256
FPPOPTIONS += -DINTRINSIC_SECNDS -DINTRINSIC_ETIME
FPPOPTIONS += -DFPP_HAVE_FLUSH
FPPOPTIONS += -DFPP_HAVE_F2003_FLUSH


SETMTIME = touch -r


## C compiler:
ifeq ($(serial),1)
  CC = $(GCCDIR)/gcc
else
  CC = mpicc
endif
CCFLAGS = -DF77_EXT_NAMES=lowercase_ -Wall


#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = ranlib

MPIINCLUDE =	# mpif90 knows it!
MPILIBS =	#-lmpich
#LAPACKLIBS =	-L/home/matveev/cvs/lib/LAPACK -llapack_g95 -lblas_g95
LAPACKLIBS =	-llapack -lblas
#LAPACKLIBS =	-llapack -lg2c
#BLASLIBS =	-L$(HOME)/lib/ATLAS -lf77blas -latlas
#BLASLIBS =	-L$(ABSOFT)/lib -L$(ABSOFT)/extras/ATLAS/Linux_INTEL32SSE2 -lf77blas -latlas
BLASLIBS =	#-lblas
# NAG, Absoft provide them for system()/etime() calls:
#YSTEMLIBS =	$(MPILIBS) -L$(ABSOFT)/lib -lU77 -lV77 -lm
SYSTEMLIBS =
PVMINCLUDE =	/usr/lib/pvm3/include
PVMLIBS =	-lfpvm3 -lpvm3
# dynamic load balancing (DLB) module comes with several variants adapted to
# different needs, see dlb module README
# The default variant 1 uses RMA objects for storage purposes
DLB_VARIANT = 0

ifeq ($(WITH_SCALAPACK),1)
  BLACS_LIB = -lblacs-openmpi
  SCALAPACK_LIB = -lscalapack-openmpi
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

# The one with underscores:
COMMCDIR =	comm/mpi_dir

# libraries, possibly classified by type:
LIBS       = \
	$(SCALAPACK_LIB) $(BLACS_LIB) \
	$(LAPACKLIBS) $(BLASLIBS) \
	$(COMMLIBS) $(SYSTEMLIBS)

#
# BGY3D  is linked  in  just as  any  other library.   It  is not  the
# responsibility of (this) Make to  assemble the library. So beware of
# version incompatibilities:
#
ifeq ($(WITH_BGY3D),1)
	BGY3D ?= $(BASEDIR)/bgy3d
	LIBS += -L$(BGY3D) -lbgy3d -lfftw3_mpi -lfftw3 -lm $(PETSC_LIB)
	CCFLAGS += -DWITH_BGY3D -I$(BGY3D)
endif

ifeq ($(WITH_GUILE),1)
    LIBS += $(shell guile-config link)
    CINCLUDE += $(shell guile-config compile)
endif

# FIXME: objects that require F90ALTFLAGS:
f90objs_altflags = $(libttfs_opt.a)
