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
    subroutine calc_radial_r2(L_a,L_b,L_c)

     !--------------------------------------------------
     ! result: radial3cmat_g     -  r2 radial factors
     !         radial3cmatG      -  gradients of  r2 radial factors
     !         radial3cmatG      -  2nd dervs of  r2 radial factors
     !--------------------------------

     use calc_3center_module,j1m=>j1,j2m=>j2,j3m=>j3,i1m=>i1,i2m=>i2,m_upM=>m_up

      implicit none

      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************


      integer(kind=i4_kind) :: Lam,i3,m,jj,LLam
      integer(kind=i4_kind) :: j1,j2,j3,i1,i2,m_up
      integer(kind=i4_kind) :: l_alph,l_beta

      !--------------------------------------------------------
      ! local temps:
      ! gamma_help_r2   -  (:,1:L_a+L_b+2,:)    - energy factors
      !                    (:,1:L_a+L_b+3,:)    - first derivs
      !                    (:,1:L_a+L_b+4,:)    - second derivs
      !
      ! jj_fac_sa                               - energy factors 
      ! jj_facG                                 - first derivs
      ! jj_facGG                                - second rerivs
      !
      ! gamma_help_gG  - required when first derivs calculated
      !
      !------------------------------------------------------------

      ghelp_shift_fac(:)=-fact0(:)/(cexps(k)*(fact0(:)+cexps(k)))

      exp_arg(:,0)=1.0_r8_kind
      exp_arg(:,1)=cexps(k)/(fact0(:)+cexps(k))

      do j1=2,L_a+L_b+L_c
        exp_arg(:,j1)=exp_arg(:,j1-1)*exp_arg(:,1)
      enddo

   gr1: if(integralpar_dervs) then

     do j=1,n_equals

      gamma_arg2(:)= gamma_arg2_vec(:,j) * exp_arg(:,1)

      gamma_help_r2(:,1:L_a+L_b+4,j)=gamma(L_a+L_b+4,gamma_arg2(:))
                   !            ^ - shape on alloc ^ - increment for 2nd dervs 

      r2_jjj_fac(:)=gamma_arg2(:)*fact0(:)/(cexps(k)*(fact0(:)+cexps(k)))

      do j1=1,1+L_a+L_b 

       gamma_help_g(:,j1,j)=gamma_help_r2(:,j1+1,j)*r2_jjj_fac(:)
                                         !  ^   -  L_a+L_b+2

       gamma_help_gG(:,j1,j)=gamma_help_r2(:,j1+2,j)*r2_jjj_fac(:) &
                    !  ^                  !  ^  -  L_a+L_b+3 for first dervs
          -gamma_help_r2(:,j1+1,j)*fact0(:)/(cexps(k)*(fact0(:)+cexps(k)))
                        !     ^ here for grads increment of energy factor, L_a+L_b+2


       gamma_help_gGG(:,j1,j)=gamma_help_r2(:,j1+3,j)*r2_jjj_fac(:) &
                     !  ^                  !  ^  -  L_a+L_b+4 for 2nd dervs
          -gamma_help_r2(:,j1+2,j)*fact0(:)/(cexps(k)*(fact0(:)+cexps(k))) &
                        !     ^ here for 2nd dervs L_a+L_b+3
          -gamma_help_r2(:,j1+2,j)*fact0(:)/(cexps(k)*(fact0(:)+cexps(k)))
      
      enddo 
    enddo 

   elseif(integralpar_gradients) then
        ! calculate both gamma_help_g and gamma_help_gG used
        ! to calculate  r2 radial factors and thei derivs

     do j=1,n_equals

      gamma_arg2(:)= gamma_arg2_vec(:,j) * exp_arg(:,1)

      gamma_help_r2(:,1:L_a+L_b+3,j)=gamma(L_a+L_b+3,gamma_arg2(:))
                   !            ^                  ^ - increment for first derivs 

      r2_jjj_fac(:)=gamma_arg2(:)*fact0(:)/(cexps(k)*(fact0(:)+cexps(k)))

      do j1=1,1+L_a+L_b 
       gamma_help_g(:,j1,j)=gamma_help_r2(:,j1+1,j)*r2_jjj_fac(:)
                                         !  ^   -  L_a+L_b+2

       gamma_help_gG(:,j1,j)=gamma_help_r2(:,j1+2,j)*r2_jjj_fac(:) &
                    !  ^                  !  ^  -  L_a+L_b+3 for first dervs
          -gamma_help_r2(:,j1+1,j)*fact0(:)/(cexps(k)*(fact0(:)+cexps(k)))

      enddo 
    enddo 

   else gr1 !i.e. no grads and dervs
     do j=1,n_equals
      gamma_arg2(:)= gamma_arg2_vec(:,j) * exp_arg(:,1)

      gamma_help_r2(:,1:L_a+L_b+2,j)=gamma(L_a+L_b+2,gamma_arg2(:))
                   !            ^                  ^    - increment for energy contribs

      r2_jjj_fac(:)=gamma_arg2(:)*fact0(:)/(cexps(k)*(fact0(:)+cexps(k)))

      do j1=1,1+L_a+L_b 
       gamma_help_g(:,j1,j)=gamma_help_r2(:,j1+1,j)*r2_jjj_fac(:)
                                         !  ^   - L_a+L_b+2
      enddo 
    enddo 
   endif gr1


        jj_fac_sa=0.0_r8_kind

        dervs_jjf: if(integralpar_dervs) then

          !   calculate jj_fac_sa,  jj_facG and jj_facGG
          ! - temps to calculate radial factors and their derivs

        jj_facG=0.0_r8_kind
        jj_facGG=0.0_r8_kind

         do LLam=llam_m,l_a+l_b+l_c
           do m_up=0,mup_llam(LLam)
           do m=0,min(m_up,LLam-llam_m)

          gamma_help_fac(:)=bin_fac(m_up, m) * exp_arg(:,LLam-m)
          gamma_help_fac_g(:)=gamma_help_fac(:)*(LLam-m)*ghelp_shift_fac(:)

          do j=1,n_equals

          jj_fac_sa(:,m_up,LLam,j)=jj_fac_sa(:,m_up,LLam,j)+gamma_help_fac(:) &
              !  ^                        ^
              *(gamma_help_r2(:,1+LLam-m,j)*g_shift_fac(:)+gamma_help_g(:,1+LLam-m,j)) & 
                            !   ^     la+lb+1 for energy contrib                  ^
              +gamma_help_fac_g(:)*gamma_help_r2(:,1+LLam-m,j)
                                                !  ^  - la+lb+1 for energy contribs

          jj_facG(:,m_up,LLam,j)=jj_facG(:,m_up,LLam,j)+gamma_help_fac(:) &
              ! ^                      ^
              *(gamma_help_r2(:,2+LLam-m,j)*g_shift_fac(:)+gamma_help_gG(:,1+LLam-m,j)) & 
                            !   ^      la+lb+2 for first dervs                 ^ - la+lb+1 for first dervs
              +gamma_help_fac_g(:)*gamma_help_r2(:,2+LLam-m,j)
                                                !  ^  - la+lb+2 for first dervs

          jj_facGG(:,m_up,LLam,j)=jj_facGG(:,m_up,LLam,j)+gamma_help_fac(:) &
              ! ^                      ^    ! temp for second dervs
              *(gamma_help_r2(:,3+LLam-m,j)*g_shift_fac(:)+gamma_help_gGG(:,1+LLam-m,j)) & 
                            !   ^      la+lb+3 for 2nd dervs           ^ - la+lb+1 for 2nd dervs
              +gamma_help_fac_g(:)*gamma_help_r2(:,3+LLam-m,j)
                                                !  ^  - la+lb+3 for 2nd dervs

          enddo 
          enddo
        enddo
        enddo

        elseif(integralpar_gradients) then

          !   calculate both jj_fac_sa and jj_facG 
          ! - temps to calculate radial factors and their derivs

        jj_facG=0.0_r8_kind

         do LLam=llam_m,l_a+l_b+l_c
           do m_up=0,mup_llam(LLam)
           do m=0,min(m_up,LLam-llam_m)

          gamma_help_fac(:)=bin_fac(m_up, m) * exp_arg(:,LLam-m)
          gamma_help_fac_g(:)=gamma_help_fac(:)*(LLam-m)*ghelp_shift_fac(:)

          do j=1,n_equals

          jj_fac_sa(:,m_up,LLam,j)=jj_fac_sa(:,m_up,LLam,j)+gamma_help_fac(:) &
              !  ^                        ^
              *(gamma_help_r2(:,1+LLam-m,j)*g_shift_fac(:)+gamma_help_g(:,1+LLam-m,j)) & 
                            !   ^     la+lb+1 for energy contrib                  ^
              +gamma_help_fac_g(:)*gamma_help_r2(:,1+LLam-m,j)
                                                !  ^  - la+lb+1 for energy contribs

          jj_facG(:,m_up,LLam,j)=jj_facG(:,m_up,LLam,j)+gamma_help_fac(:) &
              ! ^                      ^
              *(gamma_help_r2(:,2+LLam-m,j)*g_shift_fac(:)+gamma_help_gG(:,1+LLam-m,j)) & 
                            !   ^      la+lb+2 for first dervs                 ^ - la+lb+1 for first dervs
              +gamma_help_fac_g(:)*gamma_help_r2(:,2+LLam-m,j)
                                                !  ^  - la+lb+2 for first dervs
          enddo 
          enddo
        enddo
        enddo

        else dervs_jjf

        do LLam=llam_m,l_a+l_b+l_c
        do m_up=0,mup_llam(LLam)
           do m=0,min(m_up,LLam-llam_m)

          gamma_help_fac(:)=bin_fac(m_up, m) * exp_arg(:,LLam-m)
          gamma_help_fac_g(:)=gamma_help_fac(:)*(LLam-m)*ghelp_shift_fac(:)

          do j=1,n_equals
          jj_fac_sa(:,m_up,LLam,j)=jj_fac_sa(:,m_up,LLam,j)+gamma_help_fac(:) &
              *(gamma_help_r2(:,1+LLam-m,j)*g_shift_fac(:)+gamma_help_g(:,1+LLam-m,j)) &
                             !  ^                                      !  ^
              +gamma_help_fac_g(:)*gamma_help_r2(:,1+LLam-m,j)
                                                !  ^
          enddo 
          enddo
        enddo
        enddo
        endif dervs_jjf

      jj=1
      gr3: if(integralpar_dervs) then
      do j1=0,L_a
         do j2=0,L_b
           do j3=0,L_c
            if(.NOT. even_triangle(j1,j2,j3) ) cycle
            Lam = (j1 + j2 + j3)/2
            LLam=L_a+L_b+L_c-Lam
            do i1=0,L_a-j1
               do i2=0,L_b-j2
                m_up=L_a-i1+i2+j2-Lam
                l_alph=L_c+j1+i1-Lam
                l_beta=L_b-i2+j3-Lam
                   do i3=0,L_c-j3
                   do j=1,n_equals
                   radial3cmat_g(:,jj+i3,j)= &
                        i3_fac(:,m_up,l_alph-i3,l_beta+i3)*jj_fac_sa(:,m_up,LLam,j)
                   radial3cmatG(:,jj+i3,j)= &
                        i3_fac(:,m_up,l_alph-i3,l_beta+i3)*jj_facG(:,m_up,LLam,j)

                   radial3cmatGG(:,jj+i3,j)= &
                        i3_fac(:,m_up,l_alph-i3,l_beta+i3)*jj_facGG(:,m_up,LLam,j)
                   enddo 
                  enddo
                jj=jj+L_c-j3+1
               end do
            end do

           enddo
         end do
      end do


      elseif(integralpar_gradients) then
      do j1=0,L_a
         do j2=0,L_b
           do j3=0,L_c
            if(.NOT. even_triangle(j1,j2,j3) ) cycle
            Lam = (j1 + j2 + j3)/2
            LLam=L_a+L_b+L_c-Lam
            do i1=0,L_a-j1
               do i2=0,L_b-j2
                m_up=L_a-i1+i2+j2-Lam
                l_alph=L_c+j1+i1-Lam
                l_beta=L_b-i2+j3-Lam
                   do i3=0,L_c-j3
                   do j=1,n_equals
                   radial3cmat_g(:,jj+i3,j)= &
                        i3_fac(:,m_up,l_alph-i3,l_beta+i3)*jj_fac_sa(:,m_up,LLam,j)
                   radial3cmatG(:,jj+i3,j)= &
                        i3_fac(:,m_up,l_alph-i3,l_beta+i3)*jj_facG(:,m_up,LLam,j)
                   enddo 
                  enddo
                jj=jj+L_c-j3+1
               end do
            end do

           enddo
         end do
      end do

     else  gr3
      do j1=0,L_a
         do j2=0,L_b
           do j3=0,L_c
            if(.NOT. even_triangle(j1,j2,j3) ) cycle
            Lam = (j1 + j2 + j3)/2
            LLam=L_a+L_b+L_c-Lam
            do i1=0,L_a-j1
               do i2=0,L_b-j2
                m_up=L_a-i1+i2+j2-Lam
                l_alph=L_c+j1+i1-Lam
                l_beta=L_b-i2+j3-Lam
                   do i3=0,L_c-j3
                   do j=1,n_equals
                   radial3cmat_g(:,jj+i3,j)= &
                        i3_fac(:,m_up,l_alph-i3,l_beta+i3)*jj_fac_sa(:,m_up,LLam,j)
                   enddo 
                  enddo
                jj=jj+L_c-j3+1
               end do
            end do

           enddo
         end do
      end do
      endif gr3

    end subroutine calc_radial_r2
