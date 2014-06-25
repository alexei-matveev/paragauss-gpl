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
    subroutine calc_radial_nucP(L_a,L_b)

    ! used to calc Nuc ewpc and solvation contribs 
        use calc_3center_module, &
          j1m=>j1,j2m=>j2,i1m=>i1,i2m=>i2,m_upM=>m_up 
      implicit none
      !------------ Declaration of formal parameters ---------------
      integer(kind=i4_kind), intent(in) :: L_a,L_b
      !** End of interface *****************************************
      integer(kind=i4_kind)  :: Lam,j1,j2,i1,i2,m,jj,llam,m_up


      gamma_arg2(:)=( gamma_argXC(:,1)**2+gamma_argXC(:,2)**2 &
                     +gamma_argXC(:,3)**2 ) * fact0

      gamma_help(:,1:N_max+1)=gamma(N_max+1,gamma_arg2(:))*z

      jj_fac=0.0_r8_kind
      do LLam=llam_m,l_a+l_b
        do m_up=0,mup_llam(LLam)
        do m=0,min(m_up,LLam-llam_m)
        jj_fac(:,m_up,LLam)=jj_fac(:,m_up,LLam)  &
                               +bin_fac(m_up, m)*gamma_help(:,1+LLam-m)
        enddo
        enddo
      enddo
     
      jj=0
      do j1=0,L_a
         do j2=0,L_b
            if(.NOT. even_triangle(j1,j2,0) ) cycle
            Lam = (j1 + j2)/2
            LLam=L_a+L_b-Lam
            do i1=0,L_a-j1
               do i2=0,L_b-j2
                  jj=jj+1
                  m_up=L_a-i1+i2+j2-Lam
                  radialNuc_mat(:,jj)=jj_fac(:,m_up,LLam) &
                                     *nucR_fac(:,m_up,j1+i1-Lam,L_b-i2-Lam)
               end do
            end do
         end do
      end do

    end subroutine calc_radial_nucP
