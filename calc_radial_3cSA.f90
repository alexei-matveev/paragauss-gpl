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
    subroutine calc_radial_3cSA(L_a,L_b,L_c,N_max)

     !  result:    radial3cmat  
     !             gamma_help_co   -  to be used in calc_radial_3cGa
     !             exp_arg         -  to be used in calc_radial_3cGa

     use calc3c_switches
     use calc_3center_module,j1m=>j1,j2m=>j2,j3m=>j3,i1m=>i1,i2m=>i2,m_upM=>m_up, &
                             N_maxM=>N_max
      implicit none

      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b,L_c,N_max
      !** End of interface *****************************************

      integer(kind=i4_kind) :: Lam,i3,m,jj,LLam,i3m,i3p
      integer(kind=i4_kind) :: j1,j2,j3,i1,i2,m_up

      exp_arg(:,0)=1.0_r8_kind
      exp_arg(:,1)=cexps(k)/(fact0(:)+cexps(k))

     do j=1,n_equals
      gamma_arg2(:)= gamma_arg2_vec(:,j) * exp_arg(:,1)
      gamma_help_co(:,1:N_max+1,j)=gamma(N_max+1,gamma_arg2(:))
      !                 la+lb+lc+1 in this sub
      !                 la+lb+lc+2 in calc_radial_3cGa
      !                 la+lb+lc+3 in calc_radial_3cGGa
    enddo 

!    if(gamma_var_fixed) gamma_help_co(:,1:N_max+1,j)=1.0_r8_kind
   ! it can be fixed with gamma_arg2_vec fixed

        do j1=2,L_a+L_b+L_c
         exp_arg(:,j1)=exp_arg(:,j1-1)*exp_arg(:,1)
        enddo

      jj_fac_SA=0.0_r8_kind   ! locally used temp
      do LLam=llam_m,l_a+l_b+l_c
        do m_up=0,mup_llam(LLam)
           do m=0,min(m_up,LLam-llam_m)
          gamma_help_fac(:)=bin_fac(m_up, m) * exp_arg(:,LLam-m)
          do j=1,n_equals
          jj_fac_SA(:,m_up,LLam,j)=jj_fac_SA(:,m_up,LLam,j)  &
                               +gamma_help_fac(:)*gamma_help_co(:,1+LLam-m,j)
                                                              !   ^ energy contrib
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
                    radial3cmat(:,jj+i3,j)=jj_fac_SA(:,m_up,LLam,j)*i3_fac(:,m_up,i3m-i3,i3p+i3)
                   enddo 
                  enddo
                jj=jj+L_c-j3+1
               end do
            end do

           enddo
         end do
      end do

    end subroutine calc_radial_3cSA
