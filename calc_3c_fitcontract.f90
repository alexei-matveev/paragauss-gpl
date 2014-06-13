!
! ParaGauss, a program package for high-performance computations
! of molecular systems
! Copyright (C) 2014
! T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
! M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
! A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
! T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
! M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
! M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License version 2 as published
! by the Free Software Foundation [1].
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
      subroutine calc_3c_fitcontract(unique_c,equalb)
#include "def.h"
     use calc_3center_module,i_renamed=>i,index_renamed=>index!!!, & unique_c_renamed=>unique_c
     use calc3c_switches
     use fitcontract_module
     use fit_coeff_module, only: fit_coeff_n_ch
     use gradient_data_module, only: cpks_gradient_totalsym
     use cpksdervs_matrices,only: max_cpks_coul_grad_size, cpks_coul_grad_size
      use ll_calculate_grads_module, only: model_density &
                                         , spin_index &
                                         , split_gradients &
                                         , spin_polarized &
                                         , grad_mat_spin &
                                         , add_ca_dervs
      USE_MEMLOG

      implicit none
      integer(kind=i4_kind), intent(in):: unique_c,equalb
      ! *** end of interface ***

      integer(i4_kind):: k_gr, k2dr, index, i, i_grad

      i=unique_c


      fitcont_close: if(integralpar_3c_co.or.integralpar_3cob_grad) then

       if(new_3c_co.and.integralpar_3c_co) &
                       call fitcontract('ch', num, i, cutoff, coul_int)

          close_coul_m: if(new_3c_co) then
              do lc = -1, 0
              if (allocated(coul_int(lc)%m)) then
                 MEMLOG(-size(coul_int(lc)%m))
                 deallocate(coul_int(lc)%m, stat=alloc_stat(14))
                 ASSERT(alloc_stat(14).eq.0)
                endif
              enddo

              do lc =  1, lmax_ch
                n_independent_fcts  = &
                      ua_pointer%symadapt_partner(1,Lc)%n_independent_fcts
              if(n_independent_fcts.ne.0) then
               MEMLOG(-size(coul_int(lc)%m))
               deallocate(coul_int(lc)%m,stat=alloc_stat(14))
                 ASSERT(alloc_stat(14).eq.0)
                endif
               enddo
           endif close_coul_m

               deallocate (coul_int, STAT=alloc_stat(13))
               ASSERT(alloc_stat(13).eq.0)
      endif fitcont_close


      gr_fitcont_close: if(integralpar_3cob_grad) then


        ab_dervs: do k_gr=k_gr_0,k_gr_1

         if(integralpar_cpksdervs) &
#ifdef no_cpks_coul_grads
          cpks_grad_mat_p=>cpks_grad_mat(:,k_gr)      ! 1st contrib(a,b)
#else
          cpks_grad_mat_p=>cpks_grad_mat(:,:,:,:,:,k_gr)      ! 1st contrib(a,b)
#endif

         sp_mda: if (model_density .and. spin_polarized) then
          if(new_3c_fc.and.new_3c_co_grad)  &
          call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr)%l,&
                           grad_mat(:,:,:,:,k_gr), grad_mat_spin(:,:,:,:,k_gr) )

         else sp_mda ! i.e regular nonmda run

          if(new_3c_fc.and.new_3c_co_grad) then

          if(integralpar_cpksdervs) then
           call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr)%l, &   ! cpks_gradients, 2nd contrib(a,b)
                            grad_mat(:,:,:,:,k_gr),cpks_gradients=cpks_grad_mat_p)

          else ! i.e. no integralpar_cpksdervs
           call fitcontract('grad',num,i,cutoff,coul_int_grad(k_gr)%l,grad_mat(:,:,:,:,k_gr))
          endif

!                                    **  coul ab dervs on
      if(integralpar_dervs) then
      do k2dr=1,6
      call fitcontract &
             ('grad',num,i,cutoff,coul_int_dervs(k_gr,k2dr)%l, dervs_mat(:,:,:,:,k_gr,k2dr))
      enddo
      endif

         endif
        endif sp_mda
      enddo ab_dervs

!     print*,'sum grad_mat', sum(grad_mat)

!         print*,'coul_int_dervs_totsym unique_c ',unique_c

       fitcont_c: do i_grad=1,grad_dim

        index = gradient_index(imc) + i_grad - 1
         mda:  if (model_density) then
          if (spin_polarized) then
          spin_index = index + gradient_data_n_gradients
          if (split_gradients) then
              if(new_3c_fc.and.new_3c_co_grad) &
               call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad)%l, &
                              prim_int_3cob_coul_grad(index)%m, & ! V_H
                              prim_int_3cob_grad(spin_index)%m, & ! V_X,spin
                              prim_int_3cob_grad(index)%m )       ! V_X,tot

           else ! total_gradient only
                        if(new_3c_fc.and.new_3c_co_grad) &
                         call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad)%l, &
                              prim_int_3cob_grad(index)%m, &      ! V_H+V_X,tot
                              prim_int_3cob_grad(spin_index)%m )  ! V_X,spin
           endif

          else ! spin_restricted
                      if (split_gradients) then
                        if(new_3c_fc.and.new_3c_co_grad) &
                         call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad)%l, &
                              prim_int_3cob_coul_grad(index)%m, & ! V_H
                              mda_xcpot_gradients = &
                              prim_int_3cob_grad(index)%m )       ! V_X,tot
                      else ! total_gradient only
                       if(new_3c_fc.and.new_3c_co_grad) &
                         call fitcontract('grad',num,i,cutoff, &
                              coul_int_grad_totsym(i_grad)%l, &
                              prim_int_3cob_grad(index)%m )       ! V_H+V_X,tot
                      endif
          endif

        else mda ! standard SCF

         if(new_3c_fc.and.new_3c_co_grad) then

          cpks_dervs: if(integralpar_cpksdervs) then

          !---------------------------------------
          !** here (c) coul grads are summed up in prim_int_3cob_grad

!      prim_int_3cob_grad(index)%m=0.0_r8_kind !!!!!!

#ifdef no_cpks_coul_grads
      call fitcontract('grad',num,i,cutoff,coul_int_grad_totsym(i_grad)%l, &
                       prim_int_3cob_grad(index)%m, &
                       cpks_gradients=cpks_gradient_totalsym(:,index) ) ! 1st contrib to prim_int_cpks_coul_grad, (c) contrib
#else
!                                          ________input______________
      call fitcontract('grad',num,i,cutoff,coul_int_grad_totsym(i_grad)%l, &
                       prim_int_3cob_grad(index)%m, &
                       cpks_gradients=prim_int_cpks_coul_grad(index)%m) ! 1st contrib to prim_int_cpks_coul_grad, (c) contrib
                      !--------------------------------------------
#endif

          else cpks_dervs
          call fitcontract('grad',num,i,cutoff,coul_int_grad_totsym(i_grad)%l, &
                        prim_int_3cob_grad(index)%m )
          endif cpks_dervs

         if(integralpar_dervs)  then

          do k2dr=1,grad_dim

           ! contribs coul_int_dervs_totsym  collected for all equal atoms  of one unique (c)
           ! are now contracted and  remaped to prim_int_coul_dervs

          call fitcontract &
          ('grad',num,i,cutoff,coul_int_dervs_totsym(i_grad,k2dr)%l, &
           prim_int_coul_dervs(index,gradient_index(imc)-1+k2dr)%m )    !(1) - first of 3
          enddo

     do k2dr=1,6
      call fitcontract('grad',num,i,cutoff,coul_int_ca_dervs(i_grad,k2dr)%l, &
                                     ca_dervs_mat(:,:,:,:,i_grad,k2dr) )
     enddo

         endif

         endif ! new integrals

        endif mda
        enddo fitcont_c


       if(integralpar_dervs) then
           call add_ca_dervs(equalb,ca_dervs_mat,ima,imb,imc,prim_int_coul_dervs) !(2) - second of 3 contribs
        endif

         lcloop: do lc=-1,0

          gr_k: do k_gr=1,6
           if(allocated(coul_int_grad(k_gr)%l(lc)%m)) then
            MEMLOG(-size(coul_int_grad(k_gr)%l(lc)%m))
#if 1 /* (1) */
            if(max_cpks_coul_grad_size.lt.cpks_coul_grad_size) then
                  max_cpks_coul_grad_size=cpks_coul_grad_size
!            print*, 'calc_3c_fitcontract 1:  max_cpks_coul_grad_size', max_cpks_coul_grad_size
            endif
            cpks_coul_grad_size=cpks_coul_grad_size &
                      -size(coul_int_grad(k_gr)%l(lc)%m )
#endif
            deallocate(coul_int_grad(k_gr)%l(lc)%m,stat=alloc_stat(56)) ! 1 2
            ASSERT(alloc_stat(56).eq.0)
           endif

           if(integralpar_dervs) then

            do k2dr=1,6
            if(allocated(coul_int_dervs(k_gr,k2dr)%l(lc)%m)) then
              MEMLOG(-size(coul_int_dervs(k_gr,k2dr)%l(lc)%m))
              deallocate( coul_int_dervs(k_gr,k2dr)%l(lc)%m, &
                          stat=alloc_stat(68)) ! 1 2
              ASSERT(alloc_stat(68).eq.0)
           endif
            enddo
           endif
          enddo gr_k


          !--------------------------------
          ! for each uniqe (abc) combination
          ! just calculated block is symmetrized
          ! and remaped to prim_int_coul_dervs

          grtots: do i_grad=1,grad_dim
           if(allocated(coul_int_grad_totsym(i_grad)%l(lc)%m)) then
           MEMLOG(-size(coul_int_grad_totsym(i_grad)%l(lc)%m))
#if 1 /* lc=-1,0 */
            if(max_cpks_coul_grad_size.lt.cpks_coul_grad_size) then
                  max_cpks_coul_grad_size=cpks_coul_grad_size
             print*, 'calc_3c_fitcontract 2: max_cpks_coul_grad_size', max_cpks_coul_grad_size
            endif
            cpks_coul_grad_size=cpks_coul_grad_size- &
                     size(coul_int_grad_totsym(i_grad)%l(lc)%m)
#endif
           deallocate(coul_int_grad_totsym(i_grad)%l(lc)%m,STAT=alloc_stat(60))
           ASSERT(alloc_stat(60).eq.0)
          endif

         if(integralpar_dervs) then

          do k2dr=1,grad_dim
           if(allocated(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m)) then
           MEMLOG(-size(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m))
           deallocate(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m,STAT=alloc_stat(71))
           ASSERT(alloc_stat(71).eq.0)
            endif
          enddo

          do k2dr=1,6
          if(allocated(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m)) then
           MEMLOG(-size(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m))
           deallocate(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m,STAT=alloc_stat(75))
           ASSERT(alloc_stat(75).eq.0)
          endif
          enddo

         endif

          enddo grtots
         enddo lcloop

        if(new_3c_co) then
        looplc: do lc=1,lmax_ch
                n_independent_fcts  = &
                      ua_pointer%symadapt_partner(1,Lc)%n_independent_fcts

       indx: if(n_independent_fcts.ne.0) then
          do k_gr=1,6
          if(allocated(coul_int_grad(k_gr)%l(lc)%m)) then
            MEMLOG(-size(coul_int_grad(k_gr)%l(lc)%m))
#if 1 /* (4) */
            if(max_cpks_coul_grad_size.lt.cpks_coul_grad_size) then
                  max_cpks_coul_grad_size=cpks_coul_grad_size
             print*, 'calc_3c_fitcontract 3: max_cpks_coul_grad_size', max_cpks_coul_grad_size
            endif
            cpks_coul_grad_size=cpks_coul_grad_size- &
                     size(coul_int_grad(k_gr)%l(lc)%m)
#endif
            deallocate(coul_int_grad(k_gr)%l(lc)%m,stat=alloc_stat(56))
            ASSERT(alloc_stat(56).eq.0)
           endif

           if(integralpar_dervs) then
            do k2dr=1,6
           if(allocated(coul_int_dervs(k_gr,k2dr)%l(lc)%m)) then
              MEMLOG(-size(coul_int_dervs(k_gr,k2dr)%l(lc)%m))
       deallocate(coul_int_dervs(k_gr,k2dr)%l(lc)%m,stat=alloc_stat(68))
       ASSERT(alloc_stat(68).eq.0)
             endif
            enddo
           endif

          enddo

          do i_grad=1,grad_dim

          if(allocated(coul_int_grad_totsym(i_grad)%l(lc)%m)) then
          MEMLOG(-size(coul_int_grad_totsym(i_grad)%l(lc)%m))
#if 1 /* (3)l=1,lc_max */
            if(max_cpks_coul_grad_size.lt.cpks_coul_grad_size) then
                  max_cpks_coul_grad_size=cpks_coul_grad_size
             print*, 'calc_3c_fitcontract 4: max_cpks_coul_grad_size', max_cpks_coul_grad_size
            endif
            cpks_coul_grad_size=cpks_coul_grad_size- &
                     size(coul_int_grad_totsym(i_grad)%l(lc)%m)
#endif
           deallocate(coul_int_grad_totsym(i_grad)%l(lc)%m,STAT=alloc_stat(60))
           ASSERT(alloc_stat(60).eq.0)
          endif

          if(integralpar_dervs) then

          do k2dr=1,grad_dim
          if(allocated(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m)) then
           MEMLOG(-size(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m))
           deallocate(coul_int_dervs_totsym(i_grad,k2dr)%l(lc)%m,STAT=alloc_stat(71))
           ASSERT(alloc_stat(71).eq.0)
          endif
          enddo

          do k2dr=1,6
          if(allocated(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m)) then
          MEMLOG(-size(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m))
          deallocate(coul_int_ca_dervs(i_grad,k2dr)%l(lc)%m,STAT=alloc_stat(75))
           ASSERT(alloc_stat(75).eq.0)
          endif
          enddo

         endif
        enddo

        endif indx
       enddo looplc
       endif

            do i_grad=1,grad_dim ! only if moving_c
             deallocate(coul_int_grad_totsym(i_grad)%l,stat=alloc_stat(58))
             ASSERT(alloc_stat(58).eq.0)

             if(integralpar_dervs) then
              do k2dr=1,grad_dim ! only if moving_c
               deallocate(coul_int_dervs_totsym(i_grad,k2dr)%l,stat=alloc_stat(70))
               ASSERT(alloc_stat(70).eq.0)
              enddo

              do k2dr=1,6
               deallocate(coul_int_ca_dervs(i_grad,k2dr)%l,stat=alloc_stat(74))
               ASSERT(alloc_stat(74).eq.0)
              enddo

             endif
            end do

            deallocate(coul_int_grad_totsym,stat=alloc_stat(59))
            ASSERT(alloc_stat(59).eq.0)

            if(integralpar_dervs) then
             MEMLOG(-size(ca_dervs_mat))
             deallocate(coul_int_dervs_totsym,ca_dervs_mat, &
                        coul_int_ca_dervs,stat=alloc_stat(69))
             ASSERT(alloc_stat(69).eq.0)
             alloc_stat(73)=0
            endif

          do k_gr=1,6
           deallocate (coul_int_grad(k_gr)%l,stat=alloc_stat(57))
           ASSERT(alloc_stat(57).eq.0)
           if(integralpar_dervs) then
           do k2dr=1,6   ! A x B dervs
            deallocate(coul_int_dervs(k_gr,k2dr)%l,stat=alloc_stat(67))
           enddo
           endif
          enddo

      endif gr_fitcont_close


     end subroutine calc_3c_fitcontract
