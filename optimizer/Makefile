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
# Makefile for standalone optimizer
#

all: optimizer

libttfs_data.a = \
		atom_data_module.o \
		../modules/constants.o \
		../gaussq/gaussq_data.o


libttfs_opt.a = opt_data_module.o \
		hesse_module.o \
		coordinates_module.o \
		coortype_module.o \
		frequency_module.o \
		geo_operations_module.o \
		ts_module.o \
		valence_coord_module.o \
		step_module.o \
		line_search_module.o \
		gradient_module.o \
		slspar_module.o \
		allocopt_module.o \
		optimizer.o \
                vff_hessian.o \
                thermodyn_prop_module.o \

libttfs_eis.a =	../eis/ch.o \
		../eis/htribk.o \
		../eis/htridi.o \
		../eis/rebak.o \
		../eis/reduc.o \
		../eis/rs.o \
		../eis/rsg.o \
		../eis/scopy.o \
		../eis/sscal.o \
		../eis/tql2.o \
		../eis/tqlrat.o \
		../eis/tred1.o \
		../eis/tred2.o

optimizer_objs = \
	optimizer_main.o \
	../modules/type_module.o \
	../modules/filename_module.o \
	../modules/iounitadmin_module.o \
	../modules/readwriteblocked_module.o \
	../modules/math_module.o \
	../modules/matrix_types.o \
	../modules/matrix_check.o \
	../modules/matrix_eigenval.o \
	../modules/io.o \
	gxfile.o \
	eigensolver.o \
	../assert_failed.o \
	../modules/strings.o \
	../lapack/f77_lapack.o \
	error_handler.o

# $(f90objs) needs to be set before including Make.rules
# as the latter is using this var to auto-generate dependencies:
f90objs = $(libttfs_opt.a) $(libttfs_data.a) $(libttfs_eis.a) $(optimizer_objs)

# machine.inc must not reset BASEDIR, if we do it here.
# Fix machine.inc if the build process doesnt find scripts
# from ../bin/ subdirectory:
BASEDIR=$(PWD)/..

# local settings here:
include ../machine.inc

# build rules here:
include ../Make.rules

# possibly avoid MPI-specific compiler options
serial = 1

# -DNO_COMM is not to compile and depend on MPI-staff
# -DFPP_OPTIMIZER is to disentangle PG-specific staff
#  from standalone optimizer:
FPPOPTIONS += -DNO_COMM -DFPP_OPTIMIZER

# this will create *copies* of the *.a libs in the
# optimizer dierectories. They are not guaranteed to
# be reusable as the compile options for PG and standalone
# optimizer differ:
libttfs_opt.a: $(libttfs_opt.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_data.a: $(libttfs_data.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_eis.a: $(libttfs_eis.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

# name of the executable (relative to PWD):
optimizer.exe = ./optimizer.exe

optimizer: $(optimizer.exe)

optimizer_main.o: optimizer.o opt_data_module.o

$(optimizer.exe): $(optimizer_objs) libttfs_opt.a libttfs_data.a libttfs_eis.a
	$(FC) $(LINKFLAGS) \
	-o $(optimizer.exe) \
	$(optimizer_objs) \
	-L. -lttfs_opt -lttfs_data -lttfs_eis \
	$(LAPACKLIBS) $(BLASLIBS) $(SYSTEMLIBS)
	@echo DONT FORGET TO INSTALL optimizer/optimizer.exe

clean:
	rm -f $(optimizer_objs) \
	      $(optimizer_objs:.o=.mod) \
	      $(optimizer_objs:.o=.$(src_ext)) \
	      libttfs_opt.a \
	      $(libttfs_opt.a) \
	      $(libttfs_opt.a:.o=.mod) \
	      $(libttfs_opt.a:.o=.$(src_ext)) \
	      libttfs_data.a \
	      $(libttfs_data.a) \
	      $(libttfs_data.a:.o=.mod) \
	      $(libttfs_data.a:.o=.$(src_ext)) \
	      libttfs_eis.a \
	      $(libttfs_eis.a) \
	      $(libttfs_eis.a:.o=.mod) \
	      $(libttfs_eis.a:.o=.$(src_ext))
