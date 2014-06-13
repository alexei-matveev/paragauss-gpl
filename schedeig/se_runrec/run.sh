#!/bin/sh
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
MPIPROCESSES=36
MINPROCCONF=1
MAXPROCCONF=36
BLOCKSIZE=64
MINMATSIZE=200
MAXMATSIZE=2500
STEPSIZE="25"


cd bin
mpirun -np $MPIPROCESSES ./rrec minprocconf=$MINPROCCONF maxprocconf=$MAXPROCCONF blocksize=$BLOCKSIZE minmatsize=$MINMATSIZE maxmatsize=$MAXMATSIZE stepsize=$STEPSIZE | tee run.log
./log2stat.pl run.log run.stat
./curvefit.pl run.stat run.coeff
cp run.coeff ../../cf_coefficients
mv run.coeff  run.log  run.stat ../
cd ..
