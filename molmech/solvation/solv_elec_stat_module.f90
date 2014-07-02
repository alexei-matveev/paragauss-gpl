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
! Public interface of module
!=====================================================================
module solv_elec_stat_module
!
!  This module are used to calculate electrostatic
!  contribution to the solvent effect of molecule
!
!  The module was prepared by extracting corresponding routines from
!  cpcm_module
!== Interrupt end of public interface of module =====              

!----modules used ------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use common_data_module
  use tasks_main_options_module
  use species_module
  use slab_module, only: slab_calc
  use energy_and_forces_module
  use cavity_module
  use ewald2d_module, only: calc_2d_poten
  use ewald_module, only: calc_3d_poten
  use ewald_solv_module

  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------
  public matrix_generation, el_solv_energy, solv_el_grad
  public matrix_grad,dealloc_cavity_mm,dealloc_a_inv

  real(kind=r8_kind), public, allocatable ::A_matrix_inv(:,:)
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of private constants and variables ----
  real(r8_kind), allocatable :: A_matrix(:,:)
  real(r8_kind) :: zeta7(15)=(/ &
                               1.2020462359240087, & ! 1 -  3
                               1.0369272488039460, & ! 2 -  5
                               1.0083492062659691, & ! 3 -  7  
                               1.0020083561897211, & ! 4 -  9
                               1.0004941851778486, & ! 5 - 11
                               1.0001227124391356, & ! 6 - 13
                               1.0000305872698443, & ! 7 - 15
                               1.0000076371380555, & ! 8 - 17
                               1.0000019082090243, & ! 9 - 19
                               1.0000004773717979, & !10 - 21
                               1.0000001163199652, & !11 - 23
                               1.0000000298023224, & !12 - 25
                               1.0000000074505806, & !13 - 27
                               1.0000000018626451, & !14 - 29
                               0.99999999953433871 & !15 - 31
                              /)

contains
  !********************************************************
    subroutine matrix_generation
      !generate the "interaction" matrix of charged surface areas
      !inverse = matrix for calculating the induced surface charge from the
      !molecule potential on the surface
      !** End of interface *****************************************
      use math_module, only : invert_matrix

      integer(kind=i4_kind) :: i,k,status
      real(kind=r8_kind) :: distance,A_matrix_fs
      real(kind=r8_kind) :: vect(3)
!!$      real(r8_kind) :: tt,tt1,AA

      ! generating direct matrix A for the main COSMO equation
      allocate(A_matrix_inv(n_size,n_size),stat=status)
      if ( status /= 0) call error_handler( &
           "MolMech:points_on_cavity_surface: allocation of A_MATRIX is failed")


      A_matrix_inv=zero

!!$call cpu_time(tt)
      do i=1,n_size
         do k=i,n_size
            f_s: if (i==k) then
               A_matrix_fs= &
                    1.07_r8_kind*sqrt(four*pi/tessarea(i)%area)
            else
               vect=tessarea(i)%xyz(1,:)-tessarea(k)%xyz(1,:)
               if(lattice_calc .and. use_ewald) then
                  A_matrix_fs=ewald_3d_contrib(vect)
!!$if(i==20.and.k==21) then
!!$AA=ewald_3d_contrib(vect)
!!$print*,AA
!!$print*,'#######################'
!!$end if
               else
                  distance=sqrt(dot_product(vect,vect))
                  A_matrix_fs= one/distance
                  if(slab_calc) A_matrix_fs=A_matrix_fs+periodic_2d_contrib(vect,i,k)
                  if(lattice_calc) A_matrix_fs=A_matrix_fs+periodic_3d_contrib(vect,i,k)+ &
                       distance*distance*pi*two/(three*volume)
!!$if(i==24.and.k==25) then
!!$AA=periodic_2d_contrib(vect,i,k)
!!$AA=periodic_3d_contrib(vect,i,k)+distance*distance*pi*two/(three*volume)
!!$print*,AA
!!$print*,'#######################'
!!$end if
               end if
            endif f_s
            A_matrix_inv(i,k)=A_matrix_inv(i,k)+ &
                 A_matrix_fs
            A_matrix_inv(k,i)=A_matrix_inv(i,k)
         enddo
      enddo

!!$call cpu_time(tt1)
!!$print*,'matrix_generation: ',tt,tt1,tt1-tt

      call invert_matrix(n_size,A_matrix_inv)

    end subroutine matrix_generation
  !******************************************************

  !******************************************************
  function periodic_2d_contrib(rij,i,j) result(Aij)
    !------------ Modules used -----------------------------------
    use slab_module, only: vect_s,cel_s,image_slab
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: Aij
    real(r8_kind) :: rij(3)
    integer(i4_kind) :: i,j
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: n1,n2,n3,m
    real(r8_kind) :: drij,drij_2,a12(3),da12,da12_2,sum,cos_theta,zeta,rij_a12
    real(r8_kind) :: P0,P1,P2,sum1,Aij1,ab_min,r2a2,buf
    !------------ Executable code ------------------------------------

    call image_slab(rij)

    a12=zero
    Aij=zero; Aij1=1.0_r8_kind

    drij_2=dot_product(rij,rij); drij=sqrt(drij_2)
    if(drij >=cel_s%a .or.drij >=cel_s%b ) then
       ab_min=min(cel_s%a,cel_s%b)
       print*,'The distance ',drij,' between ',i,' and ',j,'point centers'
       print*,'larger the minimal possible 2D translation: ',ab_min
       call error_handler('Solvent matrix cannot be calculated :-(')
    end if

    n1=-1
    l_n1: do
       n1=n1+1
       l_n2: do n2=-n1,n1
          if(n1==0 .and. n2==0) then
             cycle l_n2
          else
             if(gcd_rec(n1,abs(n2)) /= 1) cycle l_n2
          end if

          a12(1:2)=n1*vect_s%v1+n2*vect_s%v2
          rij_a12=dot_product(rij,a12)
          da12_2=dot_product(a12,a12); da12=sqrt(da12_2)
          cos_theta=rij_a12/(drij*da12)
          if     (cos_theta < -one) then
             cos_theta=-one
          else if(cos_theta >  one) then
             cos_theta= one
          end if
          r2a2=drij_2/da12_2

          P0=one; P1=cos_theta
          m=0; sum=zero; sum1=10_r8_kind; buf=one
          l_ma: do
             m=m+1

             if(m > 15) then
                zeta=one
             else
                zeta=zeta7(m)
             end if
 
             buf=buf*r2a2 !r2a2**m

             P2=Pol_legendre(2*m,P1,P0,cos_theta)
             sum=sum+buf*zeta*P2
             if(abs(sum-sum1) <= small01) exit l_ma
             sum1=sum
             P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
          end do l_ma

          Aij=Aij+sum*two/da12
       end do l_n2

       l_n3: do n3=0,n1-1
          if(gcd_rec(n3,n1) /= 1) cycle l_n3

          a12(1:2)=n3*vect_s%v1+n1*vect_s%v2
          rij_a12=dot_product(rij,a12)
          da12_2=dot_product(a12,a12); da12=sqrt(da12_2)
          cos_theta=rij_a12/(drij*da12)
          if     (cos_theta < -one) then
             cos_theta=-one
          else if(cos_theta >  one) then
             cos_theta= one
          end if
          r2a2=drij_2/da12_2

          P0=one; P1=cos_theta
          m=0; sum=zero; sum1=10_r8_kind; buf=one
          l_mb: do
             m=m+1
             if(m > 15) then
                zeta=one
             else
                zeta=zeta7(m)
             end if

             buf=buf*r2a2 !r2a2**m

             P2=Pol_legendre(2*m,P1,P0,cos_theta)
             sum=sum+buf*zeta*P2
             if(abs(sum-sum1) <= small01) exit l_mb
             sum1=sum
             P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
          end do l_mb

          Aij=Aij+sum*two/da12

          if(n3==0) cycle l_n3 !??

          a12(1:2)=n3*vect_s%v1-n1*vect_s%v2
          rij_a12=dot_product(rij,a12)
          da12_2=dot_product(a12,a12); da12=sqrt(da12_2)
          cos_theta=rij_a12/(drij*da12)
          if     (cos_theta < -one) then
             cos_theta=-one
          else if(cos_theta >  one) then
             cos_theta= one
          end if
          r2a2=drij_2/da12_2

          P0=one; P1=cos_theta
          m=0; sum=zero; sum1=10_r8_kind; buf=one
          l_mc: do
             m=m+1

             if(m > 15) then
                zeta=one
             else
                zeta=zeta7(m) !2m+1
             end if

             buf=buf*r2a2 !r2a2**m

             P2=Pol_legendre(2*m,P1,P0,cos_theta)

             sum=sum+buf*zeta*P2
             if(abs(sum-sum1) <= small01) exit l_mc
             sum1=sum
             P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
          end do l_mc

          Aij=Aij+sum*two/da12
       end do l_n3

       if(abs(Aij-Aij1) <= small0*ten) exit l_n1
       Aij1=Aij
    end do l_n1

  end function periodic_2d_contrib
  !****************************************************

  !******************************************************
  function periodic_3d_contrib(rij,i,j) result(Aij)
    !------------ Modules used -----------------------------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: Aij
    real(r8_kind) :: rij(3)
    integer(i4_kind) :: i,j
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: n1,n2,n3,n4,n5,n6,m
    real(r8_kind) :: drij,drij_2,a12(3),da12,da12_2,sum,cos_theta,zeta,rij_a12
    real(r8_kind) :: P0,P1,P2,sum1,Aij1,ab_min,r2a2,buf
    !------------ Executable code ------------------------------------

    call image(rij)

    drij_2=dot_product(rij,rij); drij=sqrt(drij_2)
    if(drij >=cel%a .or.drij >=cel%b .or. drij >= cel%c ) then
       ab_min=min(cel%a,cel%b,cel%c)
       print*,'The distance ',drij,' between ',i,' and ',j,'point centers'
       print*,'larger the minimal possible 3D translation: ',ab_min
       call error_handler('Solvent matrix cannot be calculated :-(')
    end if

    Aij=zero; Aij1=1.0_r8_kind

    n3=-1
    l_n3: do
       n3=n3+1
       l_n1: do n1=-n3,n3
          l_n2: do n2=-n3,n3
             if(n1==0 .and. n2==0 .and.n3==0) then
                cycle l_n2
             else
                if(gcd_rec(gcd_rec(n3,abs(n1)),abs(n2)) /= 1) cycle l_n2
             end if

             a12=n1*vect%v1+n2*vect%v2+n3*vect%v3
             rij_a12=dot_product(rij,a12)
             da12_2=dot_product(a12,a12); da12=sqrt(da12_2)
             cos_theta=rij_a12/(drij*da12)
             if     (cos_theta < -one) then
                cos_theta=-one
             else if(cos_theta >  one) then
                cos_theta= one
             end if
             r2a2=drij_2/da12_2

             P0=one; P1=cos_theta
             m=0; sum=zero; sum1=10_r8_kind; buf=one
             l_ma: do
                m=m+1

                if(m > 15) then
                   zeta=one
                else
                   zeta=zeta7(m)
                end if
 
                buf=buf*r2a2 !r2a2**m

                P2=Pol_legendre(2*m,P1,P0,cos_theta)
                sum=sum+buf*zeta*P2
                if(abs(sum-sum1) <= small01) exit l_ma
                sum1=sum
                P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
             end do l_ma

             Aij=Aij+sum*two/da12
          end do l_n2
       end do l_n1

       l_n6: do n6=0,n3-1
          l_n4: do n4=-n3,n3,2*n3
             if(n6==0.and.n4<0) cycle l_n4
             l_n5: do n5=-n3,n3
                if(gcd_rec(gcd_rec(n6,abs(n4)),abs(n5)) /= 1) cycle l_n5

                a12=n4*vect%v1+n5*vect%v2+n6*vect%v3
                rij_a12=dot_product(rij,a12)
                da12_2=dot_product(a12,a12); da12=sqrt(da12_2)
                cos_theta=rij_a12/(drij*da12)
                if     (cos_theta < -one) then
                   cos_theta=-one
                else if(cos_theta >  one) then
                   cos_theta= one
                end if
                r2a2=drij_2/da12_2

                P0=one; P1=cos_theta
                m=0; sum=zero; sum1=10_r8_kind; buf=one
                l_mb: do
                   m=m+1
                   if(m > 15) then
                      zeta=one
                   else
                      zeta=zeta7(m)
                   end if

                   buf=buf*r2a2 !r2a2**m

                   P2=Pol_legendre(2*m,P1,P0,cos_theta)
                   sum=sum+buf*zeta*P2
                   if(abs(sum-sum1) <= small01) exit l_mb
                   sum1=sum
                   P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
                end do l_mb

                Aij=Aij+sum*two/da12
             end do l_n5
          end do l_n4

          l_n5a: do n5=-n3,n3,2*n3
             l_n4a: do n4=-(n3-1),n3-1
                if(n6==0.and.n4<0) cycle l_n4a
                if(n6==0.and.n4==0.and.n5<0) cycle l_n4a
                if(gcd_rec(gcd_rec(n6,abs(n4)),abs(n5)) /= 1) cycle l_n4a

                a12=n4*vect%v1+n5*vect%v2+n6*vect%v3
                rij_a12=dot_product(rij,a12)
                da12_2=dot_product(a12,a12); da12=sqrt(da12_2)
                cos_theta=rij_a12/(drij*da12)
                if     (cos_theta < -one) then
                   cos_theta=-one
                else if(cos_theta >  one) then
                   cos_theta= one
                end if
                r2a2=drij_2/da12_2

                P0=one; P1=cos_theta
                m=0; sum=zero; sum1=10_r8_kind; buf=one
                l_mc: do
                   m=m+1

                   if(m > 15) then
                      zeta=one
                   else
                      zeta=zeta7(m) !2m+1
                   end if

                   buf=buf*r2a2 !r2a2**m

                   P2=Pol_legendre(2*m,P1,P0,cos_theta)

                   sum=sum+buf*zeta*P2
                   if(abs(sum-sum1) <= small01) exit l_mc
                   sum1=sum
                   P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
                end do l_mc

                Aij=Aij+sum*two/da12
             end do l_n4a
          end do l_n5a
       end do l_n6

       if(abs(Aij-Aij1) <= small0*ten) exit l_n3
       Aij1=Aij
    end do l_n3

  end function periodic_3d_contrib
  !****************************************************

  !****************************************************
  function Pol_legendre(n1,P_n,P_n_1,x) result(P_n1)
    !calulate P(x)_n+1 Legendre polynomial using P(x)_n and P(x)_n-1
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: P_n1
    real(r8_kind), intent(in) :: x, P_n,P_n_1
    integer(i4_kind), intent(in) :: n1
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: n
    !------------ Executable code ------------------------------------

    n=float(n1-1)

    P_n1=((two*n+one)*P_n*x-n*P_n_1)/(n+one)

  end function Pol_legendre
  !****************************************************

  !****************************************************
  function dPol_legendre(n,P_n,P_n_1,x) result(dP_n)
    !calulate P(x)_n+1 Legendre polynomial using P(x)_n and P(x)_n-1
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: dP_n
    real(r8_kind), intent(in) :: x, P_n,P_n_1
    integer(i4_kind), intent(in) :: n
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: n1
    integer(i4_kind) :: n2
    !------------ Executable code ------------------------------------

    n1=float(n); n2=n+1

    if(x == one) then
       dP_n=n1*(n1+one)/two
    elseif(x == -one) then
       dP_n=((-one)**n2)*n1*(n1+one)/two
    else
       dP_n=(P_n*x-P_n_1)*n1/(x*x-one)
    end if

  end function dPol_legendre
  !****************************************************

  !****************************************************
  recursive function gcd_rec(u1, v1) result(gcd)
    !Greatest_common_divisor
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind)             :: gcd
    integer(i4_kind), intent(in) :: u1, v1
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind)             :: u, v
    !------------ Executable code ------------------------------------
    if(u1==0.and.v1==0) then
       gcd=1
       return
    end if
    if(v1 > u1) then
       u=v1; v=u1
    else
       u=u1; v=v1
    end if
    if(u==0) then
       u=v
    elseif(v==0) then
       v=u
    end if

    if (mod(u, v) /= 0) then
       gcd = gcd_rec(v, mod(u, v))
    else
       gcd = v
    end if

   end function gcd_rec
  !****************************************************

  !****************************************************
  subroutine dealloc_A_inv
   ! deallocates the inverse of the cavity matrix A
   ! (A_inv is needed to calculate the induced surface charges)
   !** End of interface *****************************************
    integer(kind=i4_kind) :: status

    deallocate(A_matrix_inv,stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:dealloc_A_inv: deallocation of A_MATRIX_INV is failed")

  end subroutine dealloc_A_inv
  !***************************************************

  !***************************************************
  subroutine matrix_grad()
    !gradient of the cavity matrix A

    integer(kind=i4_kind) :: i,k,ng,lg
    real(kind=r8_kind) :: distance,help,eps_help,help_d,sum
    real(kind=r8_kind) :: vect(3),gradient(3),Qi,Qk,dvect(3,3)
!!$    real(r8_kind) :: tt,tt1,AA(3)

    eps_help=dielectric_constant/(two*(dielectric_constant-1.0_r8_kind))

!!$Grad=zero

!!$call cpu_time(tt)
    ng_N: do ng=1,n_species
       gradient(:)=zero
       i_n_size: do i=1,n_size
          Qi=Q_at(i)
          help_d=-1.07_r8_kind*Qi*Qi*sqrt(pi/tessarea(i)%area)/&
               tessarea(i)%area

          k_n_size: do k=i,n_size
             Qk=Q_at(k)

             f_s: if (i==k) then
                do lg=1,3
                   gradient(lg) = gradient(lg) + help_d* &
                        cagr%darea(lg,ng,1)%m(i)
                end do
             else
                vect=tessarea(i)%xyz(1,:)-tessarea(k)%xyz(1,:)
                distance=sqrt(dot_product(vect,vect))
                sum=zero
                do lg=1,3
                   dvect(lg,:)=cagr%dcenter(lg,ng,1)%m(:,i)-cagr%dcenter(lg,ng,1)%m(:,k)
                   sum=sum+dvect(lg,lg)
                end do
                if(sum == zero) cycle k_n_size
                if(lattice_calc .and. use_ewald) then
                   gradient=gradient+ewald_3d_contrib_grad(vect,dvect)*Qi*Qk*two
!!$if(i==20.and.k==21) then
!!$AA=ewald_3d_contrib_grad(vect,dvect)
!!$print*,ng,AA
!!$print*,'#######################'
!!$end if
                else
                   help=-one/distance**3
                   do lg=1,3
                      gradient(lg) = gradient(lg) + help*two*dot_product(vect,dvect(lg,:))*Qi*Qk
                   end do
                   if(slab_calc) then
                      gradient=gradient+periodic_2d_contrib_grad(vect,dvect)*Qi*Qk*two
!!$if(i==24.and.k==25) then
!!$AA=periodic_2d_contrib_grad(vect,dvect)
!!$print*,ng,AA
!!$print*,'#######################'
!!$end if
                   else if(lattice_calc) then
                      gradient=gradient+periodic_3d_contrib_grad(vect,dvect)*Qi*Qk*two+ &
                           (four*pi/(three*volume))*matmul(vect,dvect)*Qi*Qk*two
!!$if(i==24.and.k==25) then
!!$AA=periodic_3d_contrib_grad(vect,dvect)+(four*pi/(three*volume))*matmul(vect,dvect)
!!$print*,ng,AA
!!$print*,'#######################'
!!$end if
                   end if
                end if
             endif f_s
          enddo k_n_size
       enddo i_n_size
       gradient(:)=gradient(:)*eps_help*coulomb_factor
       Grad(:,ng)=Grad(:,ng)+gradient(:)
    enddo ng_N
!!$call cpu_time(tt1)
!!$print*,'matrix_grad: ',tt,tt1,tt1-tt

  end subroutine matrix_grad
  !***************************************************

  !************************************************************
  function periodic_2d_contrib_grad(rij,d_rij) result(dAij)
    !------------ Modules used -----------------------------------
    use slab_module, only: vect_s, image_slab
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: dAij(3)
    real(r8_kind) :: rij(3),d_rij(3,3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: Aij,Aij1,drij,drij_2,rij_a12,da12_2,da12
    real(r8_kind) :: cos_theta,rdr(3),a12(3),adr(3)
    real(r8_kind) :: sum,sum1,sum_d(3),aaaa,dP2,drda,buf,dr3da
    real(r8_kind) :: P0,P1,P2,r2a2,zeta,b2(3),db2(3)!,dA1(3),dA2(3),dA3(3),A1,A2,A3
    integer(i4_kind) :: n1,n2,n3,m,k
    !------------ Executable code ------------------------------------

    call image_slab(rij)

    a12=zero
    dAij=zero
!!$dA1=zero;dA2=zero;dA3=zero;A1=zero;A2=zero;A3=zero
    drij_2=dot_product(rij,rij); drij=sqrt(drij_2)

    do k=1,3
       rdr(k)=dot_product(rij,d_rij(k,:))
    end do

    Aij=zero; Aij1=1.0_r8_kind
    n1=-1
    l_n1: do
       n1=n1+1
       l_n2: do n2=-n1,n1
          if(n1==0 .and. n2==0) then
             cycle l_n2
          else
             if(gcd_rec(n1,abs(n2)) /= 1) cycle l_n2
          end if

          a12(1:2)=n1*vect_s%v1+n2*vect_s%v2
          rij_a12=dot_product(rij,a12)
          da12_2=dot_product(a12,a12); da12=sqrt(da12_2); drda=drij*da12
          dr3da=da12*drij*drij*drij
          cos_theta=rij_a12/drda
          if     (cos_theta < -one) then
             cos_theta=-one
          else if(cos_theta >  one) then
             cos_theta= one
          end if
          r2a2=drij_2/da12_2

          do k=1,3
             adr(k)=dot_product(a12,d_rij(k,:))
          end do
          b2=two*rdr/drij_2; db2=adr/drda-rij_a12*rdr/dr3da
          P0=one; P1=cos_theta
          m=0; sum=zero; sum1=10_r8_kind; sum_d=zero 
          buf=one
          l_ma: do
             m=m+1
             if(m > 15) then
                zeta=one
             else
                zeta=zeta7(m)
             end if

             buf=buf*r2a2 !r2a2**m

             P2=Pol_legendre(2*m,P1,P0,cos_theta); dP2=dPol_legendre(2*m,P2,P1,cos_theta)

             aaaa=buf*zeta
             sum=sum+aaaa*P2
             sum_d=sum_d+aaaa*(m*P2*b2+dP2*db2)
             if(abs(sum-sum1) <= small01) exit l_ma
             sum1=sum

             P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
          end do l_ma

          Aij=Aij+sum*two/da12
!!$A1=A1+sum*two/da12
!!$dA1=dA1+sum_d*two/da12
          dAij=dAij+sum_d*two/da12
       end do l_n2

       l_n3: do n3=0,n1-1
          if(gcd_rec(n3,n1) /= 1) cycle l_n3

          a12(1:2)=n3*vect_s%v1+n1*vect_s%v2
          rij_a12=dot_product(rij,a12)
          da12_2=dot_product(a12,a12); da12=sqrt(da12_2); drda=drij*da12
          dr3da=da12*drij*drij*drij
          cos_theta=rij_a12/drda
          if     (cos_theta < -one) then
             cos_theta=-one
          else if(cos_theta >  one) then
             cos_theta= one
          end if
          r2a2=drij_2/da12_2

          do k=1,3
             adr(k)=dot_product(a12,d_rij(k,:))
          end do
          b2=two*rdr/drij_2; db2=adr/drda-rij_a12*rdr/dr3da

          P0=one; P1=cos_theta
          m=0; sum=zero; sum1=10_r8_kind; sum_d=zero 
          buf=one
          l_mb: do
             m=m+1
             if(m > 15) then
                zeta=one
             else
                zeta=zeta7(m)
             end if

             buf=buf*r2a2 !r2a2**m

             P2=Pol_legendre(2*m,P1,P0,cos_theta); dP2=dPol_legendre(2*m,P2,P1,cos_theta)

             aaaa=buf*zeta
             sum=sum+aaaa*P2
             sum_d=sum_d+aaaa*(m*P2*b2+dP2*db2)
             if(abs(sum-sum1) <= small01) exit l_mb
             sum1=sum

             P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
          end do l_mb

          Aij=Aij+sum*two/da12
!!$A2=A2+sum*two/da12
!!$dA2=dA2+sum_d*two/da12
          dAij=dAij+sum_d*two/da12

          if(n3==0) cycle l_n3 !??

          a12(1:2)=n3*vect_s%v1-n1*vect_s%v2
          rij_a12=dot_product(rij,a12)
          da12_2=dot_product(a12,a12); da12=sqrt(da12_2); drda=drij*da12
          dr3da=da12*drij*drij*drij
          cos_theta=rij_a12/drda
          if     (cos_theta < -one) then
             cos_theta=-one
          else if(cos_theta >  one) then
             cos_theta= one
          end if
          r2a2=drij_2/da12_2

          do k=1,3
             adr(k)=dot_product(a12,d_rij(k,:))
          end do
          b2=two*rdr/drij_2; db2=adr/drda-rij_a12*rdr/dr3da

          P0=one; P1=cos_theta
          m=0; sum=zero; sum1=10_r8_kind; sum_d=zero 
          buf=one
          l_mc: do
             m=m+1
             if(m > 15) then
                zeta=one
             else
                zeta=zeta7(m)
             end if

             buf=buf*r2a2 !r2a2**m

             P2=Pol_legendre(2*m,P1,P0,cos_theta); dP2=dPol_legendre(2*m,P2,P1,cos_theta)

             aaaa=buf*zeta
             sum=sum+aaaa*P2
             sum_d=sum_d+aaaa*(m*P2*b2+dP2*db2)
             if(abs(sum-sum1) <= small01) exit l_mc
             sum1=sum

             P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
          end do l_mc

          Aij=Aij+sum*two/da12
!!$A3=A3+sum*two/da12
!!$dA3=dA3+sum_d*two/da12
          dAij=dAij+sum_d*two/da12
       end do l_n3

       if(abs(Aij-Aij1) <= small0*ten) exit l_n1
       Aij1=Aij
    end do l_n1

!!$print*,Aij,A1+A2+A3
!!$print*,A1,A2,A3
!!$print*,'======================='
!!$print*,dAij
!!$print*,dA1+dA2+dA3
!!$print*,'-----------------------'
!!$print*,dA1
!!$print*,dA2
!!$print*,dA3
!!$print*,'***********************'
  end function periodic_2d_contrib_grad
  !***************************************************

  !************************************************************
  function periodic_3d_contrib_grad(rij,d_rij) result(dAij)
    !------------ Modules used -----------------------------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: dAij(3)
    real(r8_kind) :: rij(3),d_rij(3,3)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    real(r8_kind) :: Aij,Aij1,drij,drij_2,rij_a12,da12_2,da12
    real(r8_kind) :: cos_theta,rdr(3),a12(3),adr(3)
    real(r8_kind) :: sum,sum1,sum_d(3),aaaa,dP2,drda,buf,dr3da
    real(r8_kind) :: P0,P1,P2,r2a2,zeta,b2(3),db2(3)!,dA1(3),dA2(3),dA3(3),A1,A2,A3
    integer(i4_kind) :: n1,n2,n3,n4,n5,n6,m,k
    !------------ Executable code ------------------------------------

    call image(rij)

    dAij=zero
!!$dA1=zero;dA2=zero;dA3=zero;A1=zero;A2=zero;A3=zero
    drij_2=dot_product(rij,rij); drij=sqrt(drij_2)

    do k=1,3
       rdr(k)=dot_product(rij,d_rij(k,:))
    end do

    Aij=zero; Aij1=1.0_r8_kind
    n3=-1
    l_n3: do
       n3=n3+1
       l_n1: do n1=-n3,n3
          l_n2: do n2=-n3,n3
             if(n1==0 .and. n2==0 .and. n3==0) then
                cycle l_n2
             else
                if(gcd_rec(gcd_rec(n3,abs(n1)),abs(n2)) /= 1) cycle l_n2
             end if

             a12=n1*vect%v1+n2*vect%v2+n3*vect%v3
             rij_a12=dot_product(rij,a12)
             da12_2=dot_product(a12,a12); da12=sqrt(da12_2); drda=drij*da12
             dr3da=da12*drij*drij*drij
             cos_theta=rij_a12/drda
             if     (cos_theta < -one) then
                cos_theta=-one
             else if(cos_theta >  one) then
                cos_theta= one
             end if
             r2a2=drij_2/da12_2

             do k=1,3
                adr(k)=dot_product(a12,d_rij(k,:))
             end do
             b2=two*rdr/drij_2; db2=adr/drda-rij_a12*rdr/dr3da
             P0=one; P1=cos_theta
             m=0; sum=zero; sum1=10_r8_kind; sum_d=zero 
             buf=one
             l_ma: do
                m=m+1
                if(m > 15) then
                   zeta=one
                else
                   zeta=zeta7(m)
                end if

                buf=buf*r2a2 !r2a2**m

                P2=Pol_legendre(2*m,P1,P0,cos_theta); dP2=dPol_legendre(2*m,P2,P1,cos_theta)

                aaaa=buf*zeta
                sum=sum+aaaa*P2
                sum_d=sum_d+aaaa*(m*P2*b2+dP2*db2)
                if(abs(sum-sum1) <= small01) exit l_ma
                sum1=sum

                P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
             end do l_ma

             Aij=Aij+sum*two/da12
!!$A1=A1+sum*two/da12
!!$dA1=dA1+sum_d*two/da12
             dAij=dAij+sum_d*two/da12
          end do l_n2
       end do l_n1

       l_n6: do n6=0,n3-1
          l_n4: do n4=-n3,n3,2*n3
             if(n6==0.and.n4<0) cycle l_n4
             l_n5: do n5=-n3,n3
                if(gcd_rec(gcd_rec(n6,abs(n4)),abs(n5)) /= 1) cycle l_n5

                a12=n4*vect%v1+n5*vect%v2+n6*vect%v3
                rij_a12=dot_product(rij,a12)
                da12_2=dot_product(a12,a12); da12=sqrt(da12_2); drda=drij*da12
                dr3da=da12*drij*drij*drij
                cos_theta=rij_a12/drda
                if     (cos_theta < -one) then
                   cos_theta=-one
                else if(cos_theta >  one) then
                   cos_theta= one
                end if
                r2a2=drij_2/da12_2

                do k=1,3
                   adr(k)=dot_product(a12,d_rij(k,:))
                end do
                b2=two*rdr/drij_2; db2=adr/drda-rij_a12*rdr/dr3da

                P0=one; P1=cos_theta
                m=0; sum=zero; sum1=10_r8_kind; sum_d=zero 
                buf=one
                l_mb: do
                   m=m+1
                   if(m > 15) then
                      zeta=one
                   else
                      zeta=zeta7(m)
                   end if

                   buf=buf*r2a2 !r2a2**m

                   P2=Pol_legendre(2*m,P1,P0,cos_theta); dP2=dPol_legendre(2*m,P2,P1,cos_theta)

                   aaaa=buf*zeta
                   sum=sum+aaaa*P2
                   sum_d=sum_d+aaaa*(m*P2*b2+dP2*db2)
                   if(abs(sum-sum1) <= small01) exit l_mb
                   sum1=sum

                   P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
                end do l_mb

                Aij=Aij+sum*two/da12
!!$A2=A2+sum*two/da12
!!$dA2=dA2+sum_d*two/da12
                dAij=dAij+sum_d*two/da12
             end do l_n5
          end do l_n4

          l_n5a: do n5=-n3,n3,2*n3
             l_n4a: do n4=-(n3-1),n3-1
                if(n6==0.and.n4<0) cycle l_n4a
                if(n6==0.and.n4==0.and.n5<0) cycle l_n4a
                if(gcd_rec(gcd_rec(n6,abs(n4)),abs(n5)) /= 1) cycle l_n4a

                a12=n4*vect%v1+n5*vect%v2+n6*vect%v3
                rij_a12=dot_product(rij,a12)
                da12_2=dot_product(a12,a12); da12=sqrt(da12_2); drda=drij*da12
                dr3da=da12*drij*drij*drij
                cos_theta=rij_a12/drda
                if     (cos_theta < -one) then
                   cos_theta=-one
                else if(cos_theta >  one) then
                   cos_theta= one
                end if
                r2a2=drij_2/da12_2

                do k=1,3
                   adr(k)=dot_product(a12,d_rij(k,:))
                end do
                b2=two*rdr/drij_2; db2=adr/drda-rij_a12*rdr/dr3da

                P0=one; P1=cos_theta
                m=0; sum=zero; sum1=10_r8_kind; sum_d=zero 
                buf=one
                l_mc: do
                   m=m+1
                   if(m > 15) then
                      zeta=one
                   else
                      zeta=zeta7(m)
                   end if

                   buf=buf*r2a2 !r2a2**m

                   P2=Pol_legendre(2*m,P1,P0,cos_theta); dP2=dPol_legendre(2*m,P2,P1,cos_theta)

                   aaaa=buf*zeta
                   sum=sum+aaaa*P2
                   sum_d=sum_d+aaaa*(m*P2*b2+dP2*db2)
                   if(abs(sum-sum1) <= small01) exit l_mc
                   sum1=sum

                   P0=P2; P1=Pol_legendre(2*m+1,P2,P1,cos_theta)
                end do l_mc

                Aij=Aij+sum*two/da12
!!$A3=A3+sum*two/da12
!!$dA3=dA3+sum_d*two/da12
                dAij=dAij+sum_d*two/da12
             end do l_n4a
          end do l_n5a
       end do l_n6

       if(abs(Aij-Aij1) <= small0*ten) exit l_n3
       Aij1=Aij
    end do l_n3

!!$print*,Aij,A1+A2+A3
!!$print*,A1,A2,A3
!!$print*,'======================='
!!$print*,dAij
!!$print*,dA1+dA2+dA3
!!$print*,'-----------------------'
!!$print*,dA1
!!$print*,dA2
!!$print*,dA3
!!$print*,'***********************'

  end function periodic_3d_contrib_grad
  !***************************************************

  !**********************************************************
  subroutine el_solv_energy()

    use inp_out_module, only: output_device

    real(kind=r8_kind), allocatable :: V(:),r(:,:),dV_dxyz(:,:,:)
    real(r8_kind) :: xa(3), xb(3), za, dist
    real(kind=r8_kind) :: Q_at_sum, Z_sum, correct_factor_at
    integer(kind=i4_kind) :: status,typ
    integer(kind=i4_kind) :: i,j

    allocate(V(n_size),stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:el_solv_energy: allocation of V is failed")
    V=zero
    if(slab_calc .or. lattice_calc) then
       allocate(r(n_size,3),stat=status)
       ASSERT(status == 0)
       do i=1,n_size
          r(i,:)=tessarea(i)%xyz(1,:)
       end do
    end if

    if(slab_calc) then
       call calc_2d_poten(n_size,r,V,dV_dxyz,do_grads=.false.)
    else if(lattice_calc) then
       call calc_3d_poten(n_size,r,V,dV_dxyz,do_grads=.false.)
    else
       do i=1,n_size
          xb=tessarea(i)%xyz(1,:)
          do j=1,n_species
             typ=atoms_cart(j)%type
             za=atoms(typ)%charge
             xa=atoms_cart(j)%r
             dist = sqrt(sum((xa(:)-xb(:))**2))
             if(dist < 0.001_r8_kind) dist=0.001_r8_kind
             V(i) = V(i)+za/dist
          end do
       end do
    end if
!!$print*,V(20)
!!$print*,sum(V)

    allocate(Q_at(n_size),stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:el_solv_energy: allocation of Q_n is failed")

    ! calculation point charges due to the nuclear part of the solute molecule
    Q_at= -((dielectric_constant-1.0_r8_kind)/(dielectric_constant))* &
         MATMUL(A_matrix_inv,V)

    if(output_level >= 3) then
       write (output_device,*) '*******************************************************************'
       write (output_device,*) '      Calculated point charges at electrostatic cavity surface     '
       write (output_device,*) '*******************************************************************'
       write(output_device,*)  'number                   coordinates(angs)              charge     '
       write (output_device,*) '-------------------------------------------------------------------'
       do i=1,n_size
          write(output_device,'(1x,i4,3x,3(1x,f13.9),1x,f14.9)') i,tessarea(i)%xyz(1,:),Q_at(i)
       end do
    end if
    !correction of point charges respect to Gauss equation
    Z_sum=zero
    do j=1,n_species
       typ=atoms_cart(j)%type
       Z_sum=Z_sum+atoms(typ)%charge
    enddo
    Z_sum=-Z_sum*((dielectric_constant-1.0_r8_kind)/dielectric_constant)

    Q_at_sum=0.0_r8_kind
    do j=1,n_size
       Q_at_sum=Q_at_sum+Q_at(j)
    enddo

    correct_factor_at=Z_sum/Q_at_sum
    correct_factor_at=one

    Q_at=correct_factor_at*Q_at

    ! E=(1/2)*Unn
    E_solv_el=zero
    do i=1,n_size
       E_solv_el=E_solv_el+Q_at(i)*V(i)*half*coulomb_factor
    enddo
    E_solv_tot=E_solv_tot+E_solv_el
    E_total=E_total+E_solv_el

!print*,'E_solv_tot=',E_solv_el

    deallocate(V,stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech:el_solv_energy: deallocation of V is failed")
    if(slab_calc) then
       deallocate(r,stat=status)
       ASSERT(status == 0) 
    end if

    call dealloc_A_inv()

  end subroutine el_solv_energy
  !**********************************************************

  !**************************************************
  subroutine solv_el_grad()
   !calculates and adds the part of gradient
   !(induced charges) * (gradient of electrostatic potential due to nuclei)
   !** End of interface *****************************************

    real(kind=r8_kind), allocatable :: V(:),r(:,:),dV_dxyz(:,:,:)
    real(kind=r8_kind)          :: za,dist
    integer(kind=i4_kind)       :: ma,j,l,n,i
    integer(kind=i4_kind)       :: typ,status
    real(kind=r8_kind)  :: xa(3),xb(3)
    real(kind=r8_kind) :: gradient(3),dr_n(3,3)

    if(slab_calc.or.lattice_calc) then
       allocate(dV_dxyz(n_size,n_species,3),stat=status)
       ASSERT(status == 0)
       allocate(r(n_size,3),stat=status)
       ASSERT(status == 0)
       do i=1,n_size
          r(i,:)=tessarea(i)%xyz(1,:)
       end do
    end if

    if(slab_calc.or.lattice_calc) then
       if(slab_calc)    call calc_2d_poten(n_size,r,V,dV_dxyz,do_grads=.true.)
       if(lattice_calc) call calc_3d_poten(n_size,r,V,dV_dxyz,do_grads=.true.)
!!$print*,dV_dxyz(20,6,:)
       do ma=1,n_species
          gradient=zero
          do j=1,n_size
             do l=1,n_species
                dr_n = zero
                if(l == ma) then
                   do n=1,3
                      dr_n(n,n) = one
                   enddo
                endif
                do n=1,3
                   gradient(n)=gradient(n)+coulomb_factor*Q_at(j)*dot_product(dV_dxyz(j,l,:), &
                        (dr_n(n,:) - cagr%dcenter(n,ma,1)%m(:,j)))
                end do
             end do
          end do
          Grad(:,ma)=Grad(:,ma)+gradient(:)
!print*,ma,Grad(:,ma)
       end do
    else
       unique1: do ma=1,n_species
          gradient=zero

          unique3: do j=1,n_size
             xb=tessarea(j)%xyz(1,:)

             unique5: do l=1,n_species
                typ=atoms_cart(l)%type
                za=atoms(typ)%charge
                xa=atoms_cart(l)%r

                dr_n = zero
                if(l == ma) then
                   do n=1,3
                      dr_n(n,n) = one
                   enddo
                endif
                dist = sqrt(sum((xa(:)-xb(:))**2))
                do n=1,3
                   gradient(n) = gradient(n) - coulomb_factor*za*Q_at(j)/dist**3* &
                        dot_product((xa(:)-xb(:)), &
                        (dr_n(n,:) - cagr%dcenter(n,ma,1)%m(:,j)))
                enddo
             enddo unique5
          enddo unique3
!!$print*,'el',ma,gradient
          Grad(:,ma)=Grad(:,ma)+gradient(:)
       enddo unique1
    end if

    if(slab_calc.or.lattice_calc) then
       deallocate(dV_dxyz,stat=status)
       ASSERT(status == 0)
       deallocate(r,stat=status)
       ASSERT(status == 0)
    end if

    deallocate(Q_at)

  end subroutine solv_el_grad
  !****************************************************************
end module solv_elec_stat_module
