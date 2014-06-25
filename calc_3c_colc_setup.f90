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
        subroutine calc_3c_colc_setup(la,lb)        
#include "def.h"
         use calc_3center_module
         use  cpksdervs_matrices, only: cpks_coul_grad_size
         USE_MEMLOG
         implicit none

         integer(kind=i4_kind),intent(in):: la,lb
         ! *** end of interface ***

         ! local variables
         integer(kind=i4_kind):: k_gr,k2dr,i_grad

         allocate(coul_int(lc)%m(num, ncexps, n_independent_fcts, 2*lb+1, 2*la+1), &
                  stat=alloc_stat(14)) ! 1
           MEMLOG(+size(coul_int(lc)%m))
           ASSERT(alloc_stat(14).eq.0)
           alloc_stat(14)=1
           pointer_coul => coul_int(lc)%m
           pointer_coul=0.0_r8_kind

         if(integralpar_3cob_grad) then

          do k_gr=1,6

           allocate(coul_int_grad(k_gr)%l(lc)%m &
                   (num,ncexps,n_independent_fcts,2*lb+1,2*la+1),stat=alloc_stat(56))
           ASSERT(alloc_stat(56).eq.0)
           alloc_stat(56)=1 !2
#if 1
         cpks_coul_grad_size=cpks_coul_grad_size+size(coul_int_grad(k_gr)%l(lc)%m)
#endif
           MEMLOG(+size(coul_int_grad(k_gr)%l(lc)%m))
           coul_int_grad(k_gr)%l(lc)%m=0.0_r8_kind

          if(integralpar_dervs) then
           do k2dr=1,6
            allocate(coul_int_dervs(k_gr,k2dr)%l(lc)%m &
                     (num,ncexps,n_independent_fcts,2*lb+1,2*la+1),stat=alloc_stat(68))
!           print*,'coul_int_dervs alloc',k_gr,k2dr,lc
            MEMLOG(+size(coul_int_dervs(k_gr,k2dr)%l(lc)%m))
            ASSERT(alloc_stat(68).eq.0)
            alloc_stat(68)=1 !2
            coul_int_dervs(k_gr,k2dr)%l(lc)%m=0.0_r8_kind
           enddo
          endif
          enddo

             lc_gr_tots: do i_grad=1,grad_dim ! only if moving_c

                allocate(coul_int_grad_totsym(i_grad)%l(lc)%m(num,&
                         ncexps,n_independent_fcts,2*lb+1,2*la+1),&
                     stat=alloc_stat(60))
                     ASSERT(alloc_stat(60).eq.0)
                            alloc_stat(60)=1
                MEMLOG(+size(coul_int_grad_totsym(i_grad)%l(lc)%m))
#if 1
         cpks_coul_grad_size=cpks_coul_grad_size+&
          size(coul_int_grad_totsym(i_grad)%l(lc)%m)
#endif
                coul_int_grad_totsym(i_grad)%l(lc)%m=0.0_r8_kind

              if(integralpar_dervs) then

               do k2dr=1,grad_dim
                allocate(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m(num&
                        ,ncexps,n_independent_fcts,2*lb+1,2*la+1),&
                         stat=alloc_stat(71))
               MEMLOG(+size(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m))
                     ASSERT(alloc_stat(71).eq.0)
                            alloc_stat(71)=1
                coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m=0.0_r8_kind
               enddo

               do k2dr=1,6
                allocate(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m(num,&
                           ncexps,n_independent_fcts,2*lb+1,2*la+1),&
                     stat=alloc_stat(75))
                MEMLOG(+size(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m))
                     ASSERT(alloc_stat(75).eq.0)
                            alloc_stat(75)=1
                coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m=0.0_r8_kind
               enddo

              endif
             enddo lc_gr_tots

        endif

        end subroutine calc_3c_colc_setup     
