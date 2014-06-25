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
MACH =	VPP

#### Home Directory ####
ifndef HOME
 HOME = /t/t3831ak
endif

#### Directory for MPI Executable ####
MPIBIN = $(HOME)/mpi/bin/VPP300

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

FC  = frt
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

#
# type "make NOOPT=1 modules/anything.compile"
#
# -Of full optimization
# -Wv,-m3 extended Vectorisation messages
# -Eeilmpu -Post extent of compiler messages
# -NI no external procedure inlining
# -Z <file> output file for compiler listings
# -Fixed -w fixed form sources with extended lines
# -X9 f90 language level
# -Am generate .mod file
# Emergency strategy: change optimization flags as follows:
# Normal : -Of -Wv,-m3 -Wv,-m3     -NI
# 1. Step: -Oe,-P,-E   -Wv,-qs,-m3
# 2. Step: -Ob
# 3, Step: -On
F90BASEFLAGS = -Eeilmp -Posta -NI -Am
FBASEFLAGS   = -Eeilmp -Posta -NI -Fixed -w -X9 -Am
optflags     = -Of -Wv,-m3,-noalias
optaltflags  = -Oe,-P,-E -Wv,-qs,-m3
dbgflags     = -Wv,-m3,-noalias
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
staticflags  = 
LINKFLAGS    = $(if $(STATIC), $(staticflags), )
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS)
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS))

## CPP/FPP preprocessor:
CPP = /lib/cpp
CPPOPTIONS = $(CPPUNDEF) -P -C -I$(BASEDIR)/include -I$(BASEDIR)/include

FPP = $(CPP)
#
# type "make DEBUG=1" to compile in debug features
#
FPPOPTIONS += $(CPPOPTIONS) -D_$(MACH)
FPPOPTIONS += $(if $(DEBUG), -DFPP_DEBUG=$(DEBUG), )
FPPOPTIONS += -DFPP_PARAGAUSS_VERS="\"$(paragauss_vers)\""
FPPOPTIONS += -DMAX_PATH=256

SETMTIME = touch -r


## C compiler:
CC = cc
CCFLAGS = 


#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = echo 'WARNING: VPP seems not to have ranlib, skipping:'

MPIINCLUDE =	-I/usr/lang/mpi/include
MPILIBS	=	-L/usr/lang/mpi/lib -lmpi
LAPACKLIBS =	-L/opt/lib -llapack
BLASLIBS =	-lblas
# SYSTEMLIBS =	-Wl,-P,-J -dn -lmp -lpx -lm -lelf -lc
SYSTEMLIBS =	$(MPILIBS) -Wl,-P,-J -dn -lmp -lpx -lm -lelf -lc
# -Wl,-P,-J -dn mark as parllel Vector executable

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
INCLUDE = $(patsubst %,-p%,$(DIRS)) $(COMMINCLUDE)

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

# memstat used only on VPP:
libttfs_callcomm.a += memstat.o

# FIXME: objects that require F90ALTFLAGS:
f90objs_altflags = modules/symm_module.o \
                modules/grid_module.o \
                modules/pbe_ggcxc_module.o \
                ll_calculate_grads.o \
                ls_calculate_grads.o \
                integral_calc_quad_2cob3c.o

# FIXME: objects that have been compiled with
#        different flags historically:
modules/bessel_module.o:		FFLAGS := $(F90BASEFLAGS) -Of -Wv,-m3 -Wv,-m3 -NI
modules/efm_module.o:			FFLAGS := $(F90BASEFLAGS) -On
modules/gradient_data_module.o:		FFLAGS := $(F90BASEFLAGS) -On
modules/group_module.o:			FFLAGS := $(F90BASEFLAGS) -On
modules/input_module.o:			FFLAGS := $(F90BASEFLAGS) -On
modules/solvation_module.o:		FFLAGS := $(F90BASEFLAGS) -Oe,-P,-E -Wv,-qs,-m3
modules/symm_module.o:			FFLAGS := $(F90BASEFLAGS) -On
optimizer/coordinates_data_module.o:	FFLAGS := $(F90BASEFLAGS) -On
# AM added :
ll_calculate_hfc.o:			FFLAGS := $(F90BASEFLAGS) -On
