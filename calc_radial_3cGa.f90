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
!   first and second derivatives of the l-dependent radial function

    subroutine calc_radial_3cGa(L_a,L_b,L_c)

      ! calculates first gamma derivative radial3cmatG
      ! using jj_fac_sa as local temp

      ! input:   exp_arg        - calculated with  calc_radial_3cSA
      !          gamma_help_co  - calculated with  calc_radial_3cSA
      ! result:  radial3cmatG

      use calc_3center_module, &
          j1m=>j1,j2m=>j2,j3m=>j3,i1m=>i1,i2m=>i2,m_upM=>m_up,jm=>j
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      integer(kind=i4_kind) :: Lam,j1,j2,j3,i1,i2,i3,m,jj,m_up,LLam,i3m,i3p,j

      jj_fac_sa=0.0_r8_kind
      do LLam=llam_m,l_a+l_b+l_c
        do m_up=0,mup_llam(LLam)
           do m=0,min(m_up,LLam-llam_m)
          gamma_help_fac(:)=bin_fac(m_up, m) * exp_arg(:,LLam-m)
          do j=1,n_equals
          jj_fac_sa(:,m_up,LLam,j)=jj_fac_sa(:,m_up,LLam,j)  &
                               +gamma_help_fac(:)*gamma_help_co(:,2+LLam-m,j)
                      !                                           ^ plus 1 as compared to energy integral
          enddo
          enddo
        enddo
        enddo


      jj=1
      do j1=0,L_a
         do j2=0,L_b
           do j3=0,L_c
            if(.NOT. even_triangle(j1,j2,j3) ) cycle
            Lam = (j1 + j2 + j3)/2
            LLam=L_a+L_b+L_c-Lam
            do i1=0,L_a-j1
               do i2=0,L_b-j2
                m_up=L_a-i1+i2+j2-Lam
                   i3m=L_c+j1+i1-Lam
                   i3p=L_b-i2+j3-Lam
                   do i3=0,L_c-j3
                   do j=1,n_equals
                    radial3cmatG(:,jj+i3,j)= &
                      jj_fac_sa(:,m_up,LLam,j)*i3_fac(:,m_up,i3m-i3,i3p+i3)
                   enddo
                   enddo
                jj=jj+L_c-j3+1
               end do
            end do

           enddo
         end do
      end do

    end subroutine calc_radial_3cGa

    subroutine calc_radial_3cGGa(L_a,L_b,L_c)
      ! calculates radial3cmatGG  using jj_fac_sa as auxilary
      use calc_3center_module, &
          j1m=>j1,j2m=>j2,j3m=>j3,i1m=>i1,i2m=>i2,m_upM=>m_up,jm=>j
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c
      !** End of interface *****************************************
      integer(kind=i4_kind) :: Lam,j1,j2,j3,i1,i2,i3,m,jj,m_up,LLam,i3m,i3p,j

      jj_fac_sa=0.0_r8_kind
      do LLam=llam_m,l_a+l_b+l_c
        do m_up=0,mup_llam(LLam)
        do m=0,min(m_up,LLam-llam_m)
          gamma_help_fac(:)=bin_fac(m_up, m) * exp_arg(:,LLam-m)
          do j=1,n_equals
          jj_fac_sa(:,m_up,LLam,j)=jj_fac_sa(:,m_up,LLam,j)  &
                               +gamma_help_fac(:)*gamma_help_co(:,3+LLam-m,j)
                                                               !  ^ for second derivs
                                                               !  2 for first derivs
                                                               !  1 for energy
          enddo
          enddo
        enddo
       enddo


      jj=1
      do j1=0,L_a
         do j2=0,L_b
           do j3=0,L_c
            if(.NOT. even_triangle(j1,j2,j3) ) cycle
            Lam = (j1 + j2 + j3)/2
            LLam=L_a+L_b+L_c-Lam
            do i1=0,L_a-j1
               do i2=0,L_b-j2
                m_up=L_a-i1+i2+j2-Lam
                   i3m=L_c+j1+i1-Lam
                   i3p=L_b-i2+j3-Lam
                   do i3=0,L_c-j3
                   do j=1,n_equals
                    radial3cmatGG(:,jj+i3,j)= &
                      jj_fac_sa(:,m_up,LLam,j)*i3_fac(:,m_up,i3m-i3,i3p+i3)
                   enddo
                   enddo
                jj=jj+L_c-j3+1
               end do
            end do

           enddo
         end do
      end do

    end subroutine calc_radial_3cGGa
