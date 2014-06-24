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
#
# NOTE: To set all environment variables properly call
#
# /work/System/ifort/composer_xe_2013.5.192/composer_xe_2013.5.192/mkl/bin/mklvars.sh intel64
# /work/System/ifort/composer_xe_2013.5.192/composer_xe_2013.5.192/mkl/bin/intel64/mklvars_intel64.sh
#
# before compiling
#
### Machine Keyword for FPP
MACH =	LINUX

#### Home Directory ####
HOME = /home/ttfs

#### Directory for MPI Executable ####
MPIBIN = $(HOME)/exe/openmpi

#### Directory where Sources are stored ####
BASEDIR := $(PWD)

#### Directory where Scripts are stored #####
BINDIR =	$(BASEDIR)/bin

#### Executables are stored in:
EXEDIR = $(MPIBIN)

#### NAME OF EXECUTABLE ####
exe = mainscf_$(VERS)
EXE = $(EXEDIR)/$(exe)

#### COMPILER ####
## Fortran 90/77

MPIDIR = /usr/include/mpi
FC  = ifort
F77 = $(FC)

ifeq ($(WITH_OLD_INPUT),0)
$(error NEW INPUT doesnt work with Intel compiler, set WITH_OLD_INPUT=1 in Makefile)
endif

#### Extensions compiler understands:
src_ext = F90

#### Does compiler UPPERCASE module files:
UPPERCASE.mod = 0

# Command to compare MOD-files,
# Intel seems to put "random" number
# byte positions 45:46 (decimal)
# version 9.1 at 49:50 (decimal)
# so skip 50 bytes from comparison:
CMPMOD = cmp -s -i 50

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
F90BASEFLAGS = -stand f03 -heap-arrays #-check -ftrapuv #-fpe0 -ftz
FBASEFLAGS   = -stand f03 -heap-arrays #-check -ftrapuv #-fpe0 -ftz
optflags     = -O2 #-g -traceback # -warn all # -check bounds,format,output_conversion,uninit,noarg_temp_created
optaltflags  = -O2 #-g -traceback # -warn all # -check bounds,format,output_conversion
dbgflags     = -O2 #-g -traceback -check bounds,format,output_conversion,uninit
F90FLAGS     = $(if $(NOOPT), $(dbgflags), $(optflags))
F90ALTFLAGS  = $(if $(NOOPT), $(dbgflags), $(optaltflags))
F77FLAGS     = $(optflags) $(FBASEFLAGS)
LINKFLAGS    =
# these are used in Make.rules:
FFLAGS       = $(F90BASEFLAGS) -nofor-main
FFLAGS      += $(if $(altflags), $(F90ALTFLAGS), $(F90FLAGS))

## CPP/FPP preprocessor:
CPP = /lib/cpp
CPPUNDEF = -undef
CPPOPTIONS = $(CPPUNDEF) -traditional-cpp -P -C -I$(BASEDIR)/include

FPP = $(CPP)
#
# type "make DEBUG=1" to compile in debug features
# or add it to $(debug_targets) in "misc.inc"
#
FPPOPTIONS += $(CPPOPTIONS) -D_$(MACH)
FPPOPTIONS += $(if $(DEBUG), -DFPP_DEBUG=$(DEBUG), )
#FPPOPTIONS += -DFPP_FAST_COMPILE
FPPOPTIONS += -DFPP_PARAGAUSS_VERS="\"$(paragauss_vers)\""
#
# Dont use automatic arrays, better yet use "-heap-arrays" or just
# "unlimit -s unlimited":
#
FPPOPTIONS += -DFPP_NO_BIG_AUTOMATIC_ARRAYS


SETMTIME = touch -r


## C compiler:
CC = mpicc
CCFLAGS = 


#### LDFLAGS, LIBRARY-PATH ####

AR = ar
RANLIB = ranlib

MPIINCLUDE =	-I/usr/include/mpi
MPILIBS =	-lmpi_f77 -lmpi
mkl_lapacklib =        -lmkl_lapack
threading_lib =        -lmkl_sequential
iface_lib =    -lmkl_intel_lp64
core_lib =     -lmkl_core
#LAPACKLIBS =    $(MKLROOT)/lib/intel64/libmkl_blas95_lp64 $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64 -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm
LAPACKLIBS = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group $(MKLROOT)/lib/intel64/libmkl_blacs_openmpi_lp64.a -lpthread -lm
BLASLIBS =	
SYSTEMLIBS =	
#PVMINCLUDE =	/usr/lib/pvm3/include
#PVMLIBS =	-lfpvm3 -lpvm3

# dynamic load balancing (DLB) module comes with several variants adapted to
# different needs, see dlb module README
# The default variant 1 uses RMA objects for storage purposes
DLB_VARIANT = 0

ifeq ($(WITH_SCALAPACK),1)
  BLACSLIB = $$BLACS_LIB
  SCALAPACKLIB = $$SCALAPACK_LIB
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
INCLUDE = $(patsubst %,-I%,$(DIRS)) $(COMMINCLUDE)  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include

# include path for C files:
CINCLUDE =	-I$(BASEDIR)/$(COMMDIR) $(COMMINCLUDE) $(LAPACKLIBS)

# The one with underscores:
COMMCDIR =	comm/mpi_dir

# libraries, possibly classified by type:
LIBS       = \
	$(SCALAPACK_LIB) $(BLACS_LIB) \
	$(LAPACKLIBS) $(BLASLIBS) \
	$(COMMLIBS) $(SYSTEMLIBS)

##ifeq ($(WITH_GUILE),1)
##    INCLUDE += -I$(BASEDIR)/guile
##    LIBS += -lguile
##    LINKFLAGS += -nofor-main
##endif

ifeq ($(WITH_GUILE),1)
    LIBS += $(shell guile-config link)
    CINCLUDE += $(shell guile-config compile)
    LINKFLAGS += -nofor-main
endif

# FIXME: objects that require F90ALTFLAGS:
f90objs_altflags =
