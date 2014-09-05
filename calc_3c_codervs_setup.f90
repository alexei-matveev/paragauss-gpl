!
! ParaGauss,  a program package  for high-performance  computations of
! molecular systems
!
! Copyright (C) 2014     T. Belling,     T. Grauschopf,     S. Krüger,
! F. Nörtemann, M. Staufer,  M. Mayer, V. A. Nasluzov, U. Birkenheuer,
! A. Hu, A. V. Matveev, A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman,
! D. I. Ganyushin,   T. Kerdcharoen,   A. Woiterski,  A. B. Gordienko,
! S. Majumder,     M. H. i Rotllant,     R. Ramakrishnan,    G. Dixit,
! A. Nikodem, T. Soini, M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify
! it under  the terms of the  GNU General Public License  version 2 as
! published by the Free Software Foundation [1].
!
! This program is distributed in the  hope that it will be useful, but
! WITHOUT  ANY   WARRANTY;  without  even  the   implied  warranty  of
! MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE. See  the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
         subroutine calc_3c_codervs_setup(i,la,lb)
#        include "def.h"
         use calc_3center_module,i_renamed=>i
         use  cpksdervs_matrices, only: cpks_coul_grad_size
         USE_MEMLOG
         implicit none
         integer(kind=i4_kind),intent(in)::i,la,lb
         ! *** end of interface ***

!         local vars:
          integer(kind=i4_kind):: k_gr,k2dr,i_grad

          ab_coul_gr: do k_gr=1,6

           allocate (coul_int_grad(k_gr)%l(-1:lmax_ch),stat=alloc_stat(57))
             ASSERT(alloc_stat(57).eq.0)
             alloc_stat(57)=1

             ncexps = unique_atoms(i)%r2_ch%n_exponents
             allocate(coul_int_grad(k_gr)%l(-1)%m(num,ncexps,&
                1,2*lb+1,2*la+1),stat=alloc_stat(56)) ! 1
             ASSERT(alloc_stat(56).eq.0)
             MEMLOG(+size(coul_int_grad(k_gr)%l(-1)%m))
#if 1
             cpks_coul_grad_size=cpks_coul_grad_size+&
               size(coul_int_grad(k_gr)%l(-1)%m)
#endif
             ncexps = unique_atoms(i)%l_ch(0)%n_exponents
             allocate(coul_int_grad(k_gr)%l(0)%m(num,ncexps,&
                1,2*lb+1,2*la+1), stat=alloc_stat(56)) !2
             ASSERT(alloc_stat(56).eq.0)
             MEMLOG(+size(coul_int_grad(k_gr)%l(0)%m))
#if 1
             cpks_coul_grad_size=cpks_coul_grad_size+&
                     size(coul_int_grad(k_gr)%l(0)%m)
#endif
                  alloc_stat(56)=1
             coul_int_grad(k_gr)%l(0)%m=0.0_r8_kind
             coul_int_grad(k_gr)%l(-1)%m=0.0_r8_kind

            ab_dervs: if(integralpar_dervs) then

             do k2dr=1,6
             allocate (coul_int_dervs(k_gr,k2dr)%l(-1:lmax_ch),stat=alloc_stat(67))
             ASSERT(alloc_stat(67).eq.0)
             alloc_stat(67)=1

             ncexps = unique_atoms(i)%r2_ch%n_exponents
             allocate(coul_int_dervs(k_gr,k2dr)%l(-1)%m(num,ncexps,&
                1,2*lb+1,2*la+1),stat=alloc_stat(68)) ! 1
             MEMLOG(+size(coul_int_dervs(k_gr,k2dr)%l(-1)%m))
             ASSERT(alloc_stat(68).eq.0)
             ncexps = unique_atoms(i)%l_ch(0)%n_exponents
             allocate(coul_int_dervs(k_gr,k2dr)%l(0)%m(num,ncexps,&
              1,2*lb+1,2*la+1), stat=alloc_stat(68)) !2
             MEMLOG(+size(coul_int_dervs(k_gr,k2dr)%l(0)%m))
             ASSERT(alloc_stat(68).eq.0)
                  alloc_stat(68)=1
             coul_int_dervs(k_gr,k2dr)%l(0)%m=0.0_r8_kind
             coul_int_dervs(k_gr,k2dr)%l(-1)%m=0.0_r8_kind
              enddo

            endif ab_dervs
           

          enddo ab_coul_gr



          allocate(coul_int_grad_totsym(grad_dim),stat=alloc_stat(59)) 
           ASSERT(alloc_stat(59).eq.0)
           alloc_stat(59)=1

     if(integralpar_dervs) then
       allocate(ca_dervs_mat(nbexps,naexps,2*lb+1,2*la+1,grad_dim,6), &
                coul_int_ca_dervs(grad_dim,6), &
                coul_int_dervs_totsym(grad_dim,grad_dim), &
                stat=alloc_stat(73))
          ASSERT(alloc_stat(73).eq.0)
          MEMLOG(size(ca_dervs_mat))
                 alloc_stat(73)=1
                 alloc_stat(69)=1
       ca_dervs_mat=0.0_r8_kind
     endif

            gr_totsym: do i_grad=1,grad_dim ! only if moving_c
             allocate(coul_int_grad_totsym(i_grad)%l(-1:lmax_ch),stat=alloc_stat(58))
               ASSERT(alloc_stat(58).eq.0)
               alloc_stat(58)=1 

             ncexps = unique_atoms(i)%r2_ch%n_exponents
             allocate(coul_int_grad_totsym(i_grad)%l(-1)%m&
                          (num,ncexps,1,2*lb+1,2*la+1),stat=alloc_stat(60))
             ASSERT(alloc_stat(60).eq.0)
             MEMLOG(+size(coul_int_grad_totsym(i_grad)%l(-1)%m))
#if 1
             cpks_coul_grad_size=cpks_coul_grad_size+size(coul_int_grad_totsym(i_grad)%l(-1)%m)
#endif

             ncexps = unique_atoms(i)%l_ch(0)%n_exponents
             allocate(coul_int_grad_totsym(i_grad)%l(0)%m&
                  (num,ncexps,1,2*lb+1,2*la+1),stat=alloc_stat(60))
             ASSERT(alloc_stat(60).eq.0)
             MEMLOG(+size(coul_int_grad_totsym(i_grad)%l(0)%m))
#if 1
             cpks_coul_grad_size=cpks_coul_grad_size+size(coul_int_grad_totsym(i_grad)%l(0)%m)
#endif
             coul_int_grad_totsym(i_grad)%l(0)%m=0.0_r8_kind
             coul_int_grad_totsym(i_grad)%l(-1)%m=0.0_r8_kind
             alloc_stat(60)=1
 
             if(integralpar_dervs) then

              dervs_tsalloc: do k2dr=1,grad_dim ! only if moving_c
               allocate(coul_int_dervs_totsym(i_grad,k2dr)%l(-1:lmax_ch), &
                                                       stat=alloc_stat(70))
               ASSERT(alloc_stat(70).eq.0)
               alloc_stat(70)=1 

             ncexps = unique_atoms(i)%r2_ch%n_exponents
             allocate(coul_int_dervs_totsym(i_grad,k2dr)%l(-1)%m(num,ncexps,1,2*lb+1,2*la+1), &
                                                                    stat=alloc_stat(71))
             MEMLOG(+size(coul_int_dervs_totsym(i_grad,k2dr)%l(-1)%m))
             ASSERT(alloc_stat(71).eq.0)
             ncexps = unique_atoms(i)%l_ch(0)%n_exponents
             allocate(coul_int_dervs_totsym(i_grad,k2dr)%l(0)%m(num,ncexps,1,2*lb+1,2*la+1), &
                                                                    stat=alloc_stat(71))
             MEMLOG(+size(coul_int_dervs_totsym(i_grad,k2dr)%l(0)%m))
             ASSERT(alloc_stat(71).eq.0)
             coul_int_dervs_totsym(i_grad,k2dr)%l(0)%m=0.0_r8_kind
             coul_int_dervs_totsym(i_grad,k2dr)%l(-1)%m=0.0_r8_kind
             alloc_stat(71)=1

             enddo dervs_tsalloc

             ca_alc: do k2dr=1,6
               allocate(coul_int_ca_dervs(i_grad,k2dr)%l(-1:lmax_ch), &
                                                    stat=alloc_stat(74))
               ASSERT(alloc_stat(74).eq.0)
               alloc_stat(74)=1 

             ncexps = unique_atoms(i)%r2_ch%n_exponents
             allocate(coul_int_ca_dervs(i_grad,k2dr)%l(-1)%m(num,ncexps,1,2*lb+1,2*la+1), &
                      stat=alloc_stat(75))
             MEMLOG(+size(coul_int_ca_dervs(i_grad,k2dr)%l(-1)%m))
             ASSERT(alloc_stat(75).eq.0)

             ncexps = unique_atoms(i)%l_ch(0)%n_exponents
             allocate(coul_int_ca_dervs(i_grad,k2dr)%l(0)%m(num,ncexps,1,2*lb+1,2*la+1), &
                      stat=alloc_stat(75))
             MEMLOG(+size(coul_int_ca_dervs(i_grad,k2dr)%l(0)%m))
             ASSERT(alloc_stat(75).eq.0)
             coul_int_ca_dervs(i_grad,k2dr)%l(0)%m=0.0_r8_kind
             coul_int_ca_dervs(i_grad,k2dr)%l(-1)%m=0.0_r8_kind
             alloc_stat(75)=1
             
             enddo ca_alc

            endif

            enddo gr_totsym
      end subroutine calc_3c_codervs_setup
