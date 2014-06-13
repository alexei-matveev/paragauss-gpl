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

include ../machine.inc

src =			\
	f77_x_b3.f	\
	f77_x_b88.f	\
	f77_x_ft97b.f	\
	f77_x_lda.f	\
	f77_x_pbe.f	\
	f77_x_pw91.f	\
	f77_c_ft97.f	\
	f77_c_lyp.f	\
	f77_c_p86.f	\
	f77_c_pbe.f	\
	f77_c_pw91.f	\
	f77_c_pw92.f	\
	f77_c_pz81.f	\
	f77_c_vwn5.f	\
	f77_c_vwn5rpa.f	\
	f77_xc_b3lyp.f	\
	f77_xc_b97_1.f	\
	f77_xc_b97_2.f	\
	f77_xc_b97.f	\
	f77_xc_edf1.f	\
	f77_xc_ft97.f	\
	f77_xc_hcth120.f\
	f77_xc_hcth147.f\
	f77_xc_hcth407.f\
	f77_xc_hcth.f	\
	f77_xc_pw91.f

libdftauto.a    = $(src:.f=.o)

libdftauto.a: $(libdftauto.a) libdftauto.o
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libdftauto.o: libdftauto.f90

clean:
	rm -f $(libdftauto.a) libdftauto.o libdftauto.a

#rl =  http://www.cse.clrc.ac.uk/qcg/dft
url =  http://www.cse.scitech.ac.uk/ccg/dft

.PRECIOUS : $(src)

%.f :
	wget $(url)/$*.html
	perl -ne 'if(   ( m|<PRE>| .. m|</PRE>| ) ){print;}'   $*.html | \
	perl -ne 'if( ! ( m|<PRE>| || m|</PRE>| ) ){print;}' > $*.f
	rm -f $*.html

%.o : %.f
	$(FC) $(F77FLAGS) -c $<

%.o : %.f90
	$(FC) $(F90BASEFLAGS) $(F90FLAGS) -c $<
