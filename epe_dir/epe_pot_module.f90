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
!===============================================================
! Public interface of module
!===============================================================
module epe_pot_module
  !
  !  Purpose : Calculate a electrostatic potential of epe
  !            system on choosen plane grid
  !
  !  Author: AS
  !  Date: 02/2001
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
  use type_module
  use iounitadmin_module
  use epecom_module

  implicit none
  private
  save

  !== Interrupt end of public interface of module =================
  !------------ public functions and subroutines ------------------
  public get_plane_grid,calc_pot
!==================================================================
! End of public interface of module
!==================================================================
  !-- Declaration of privat constant and variable  ----------------
  type poten
     real(kind=r8_kind) :: r(3)
     real(kind=r8_kind) :: pot
  end type poten

  type (poten), allocatable :: pot_on_plane(:,:)
  character(len=256) :: input_dir
!------------ Subroutines -----------------------------------------

contains

!******************************************************************
  subroutine get_plane_grid
    !------------ Declaration of formal parameters ----------------
    !== End of interface ==========================================
    !------------ Declaration of local variables ------------------
    real(kind=r8_kind) :: center_coor(3),xaxis(3),yaxis(3),zaxis(3)
    real(r8_kind) :: length01, r01(3), rbuf(3), mat_rot(3,3)
    integer(kind=i4_kind) :: istat,i1,j1
    !------------ Executable code ---------------------------------

    allocate(pot_on_plane(nxgrid,nygrid),stat=istat)
    if(istat /= 0) call error_handler( &
         "epe_pot_module: get_plane_grid: allocate pot_on_plane failed")

    !saving X and Y points into ParaGauss INPUT directory
    call getenv("TTFSINPUTDIR",input_dir)
    if(index(output_format,"scilab")/=0 .or. &
         index(output_format,"SCILAB")/=0) call save_X_Y_points()

    !defining coordinate system on the grid plane
    length01=sqrt(dot_product(point01-point00,point01-point00))
    r01=point01-point00

    center_coor=point00
    call vect_product(point01-point00,point02-point00,zaxis)
    zaxis=zaxis/sqrt(dot_product(zaxis,zaxis))

    xaxis=r01/sqrt(dot_product(r01,r01))

    call vect_product(zaxis,xaxis,yaxis)

!!$    print*,xaxis
!!$    print*,yaxis
!!$    print*,zaxis

    !build grid on choosen plane
    do i1=1,3
       mat_rot(i1,1)=xaxis(i1)/sqrt(dot_product(xaxis,xaxis))
       mat_rot(i1,2)=yaxis(i1)/sqrt(dot_product(yaxis,yaxis))
       mat_rot(i1,3)=zaxis(i1)/sqrt(dot_product(zaxis,zaxis))
    end do

    rbuf(3)=zero
    do i1=1,nxgrid
       rbuf(1)=dxaxis(1)+(i1-1)*(dxaxis(2)-dxaxis(1))/(nxgrid-1)
       do j1=1,nygrid
          rbuf(2)=dyaxis(1)+(j1-1)*(dyaxis(2)-dyaxis(1))/(nygrid-1)
          pot_on_plane(i1,j1)%r=matmul(mat_rot,rbuf)+center_coor
!!$          print*,pot_on_plane(i1,j1)%r
       end do
    end do

  contains
    subroutine vect_product(av,bv,cv)
      !------------ Declaration of formal parameters ------------
      real(kind=r8_kind),intent(in) :: av(3),bv(3)
      real(kind=r8_kind),intent(out) :: cv(3)
      !== End of interface ======================================
      !------------ Declaration of local variables --------------
      !------------ Executable code -----------------------------
      cv(1)=av(2)*bv(3)-av(3)*bv(2)
      cv(2)=av(3)*bv(1)-av(1)*bv(3)
      cv(3)=av(1)*bv(2)-av(2)*bv(1)
    end subroutine vect_product
!------------------------------------------------------------------

!------------------------------------------------------------------
    subroutine save_X_Y_points
      use constants, only: angstrom
      implicit none
      !== End of interface ======================================

      !------------ Declaration of local variables --------------
      real(kind=r8_kind) :: x,y
      integer(kind=i4_kind) :: grid_unit,xunit,yunit,units,i
      !------------ Executable code -----------------------------

      units=1
      if(trim(output_units) == "a.u.") units=0

      grid_unit=openget_iounit(trim(input_dir)//'/xygrid',  &
           form='formatted', status='unknown')
      write(grid_unit,'(i5)') nxgrid
      write(grid_unit,'(i5)') nygrid
      write(grid_unit,'(i5)') units
      call returnclose_iounit(grid_unit)

      xunit=openget_iounit(trim(input_dir)//'/Xpoints',  &
           form='formatted', status='unknown')
      do i=1,nxgrid
         x=dxaxis(1)+(i-1)*(dxaxis(2)-dxaxis(1))/(nxgrid-1)
         if(index(output_units,"a.u.")/=0) x = x * angstrom
         write(xunit,'(f20.15)') x
      end do
      call returnclose_iounit(xunit)

      yunit=openget_iounit(trim(input_dir)//'/Ypoints',  &
           form='formatted', status='unknown')
      do i=1,nxgrid
         y=dyaxis(1)+(i-1)*(dyaxis(2)-dyaxis(1))/(nygrid-1)
         if(index(output_units,"a.u.")/=0) y = y * angstrom
         write(yunit,'(f20.15)') y
      end do
      call returnclose_iounit(yunit)

    end subroutine save_X_Y_points

  end subroutine get_plane_grid
!******************************************************************

!******************************************************************
  subroutine calc_pot
    use culon_module, only : erf, erfc, gauss_potential
    use constants, only: angstrom
    implicit none
    !== End of interface ==========================================

    !------------ Declaration of local variables ------------------
    real(kind=r8_kind) :: V_pot,x,y
    real(kind=r8_kind) :: vrr,vrs,vrc,rr(3)
    integer(kind=i4_kind) :: vunit
    integer(kind=i4_kind) :: i,j,k,istat
    !------------ Executable code ---------------------------------

    epe(1:reg_2a_n_ions)%s(1)=R_SH_ION(1:reg_2a_n_ions,1)
    epe(1:reg_2a_n_ions)%s(2)=R_SH_ION(1:reg_2a_n_ions,2)
    epe(1:reg_2a_n_ions)%s(3)=R_SH_ION(1:reg_2a_n_ions,3)
    epe(1:reg_2a_n_ions)%c(1)=R_nuc_ION(1:reg_2a_n_ions,1)
    epe(1:reg_2a_n_ions)%c(2)=R_nuc_ION(1:reg_2a_n_ions,2)
    epe(1:reg_2a_n_ions)%c(3)=R_nuc_ION(1:reg_2a_n_ions,3)

    do j=1,reg_2a_n_ions
       epe(j)%qs=q_shell(epe(j)%k)
       epe(j)%qc=q_nuclear(epe(j)%k)
       epe(j)%q=q_ion(epe(j)%k)
    enddo

    if(index(output_format,"scilab")/=0 .or. &
         index(output_format,"SCILAB")/=0) then
       vunit=openget_iounit(trim(input_dir)//'/Vpoints',  &
            form='formatted', status='unknown')
    else if(index(output_format,"worksheet")/=0 .or. &
         index(output_format,"WORKSHEET")/=0) then
       vunit=openget_iounit(trim(input_dir)//'/V.txt',  &
            form='formatted', status='unknown')
    else if(index(output_format,"gnuplot")/=0 .or. &
         index(output_format,"GNUPLOT")/=0) then
       vunit=openget_iounit(trim(input_dir)//'/V.gpl',  &
            form='formatted', status='unknown')
    endif

    nxgrd: do i=1,nxgrid
       if(index(output_format,"worksheet")/=0 .or. &
            index(output_format,"WORKSHEET")/=0 .or. &
            index(output_format,"gnuplot")/=0 .or. &
            index(output_format,"GNUPLOT")/=0) then
          x=dxaxis(1)+(i-1)*(dxaxis(2)-dxaxis(1))/(nxgrid-1)
          if(trim(output_units) == "a.u.") x = x * angstrom
       endif
       nygrd: do j=1,nygrid
          if(index(output_format,"worksheet")/=0 .or. &
               index(output_format,"WORKSHEET")/=0 .or. &
               index(output_format,"gnuplot")/=0 .or. &
               index(output_format,"GNUPLOT")/=0) then
             y=dyaxis(1)+(j-1)*(dyaxis(2)-dyaxis(1))/(nygrid-1)
             if(trim(output_units) == "a.u.") y = y * angstrom
          endif
          rr=pot_on_plane(i,j)%r
          V_pot=gauss_potential(rr,1,1)
          reg_2a_nions: do k=1,reg_2a_n_ions
             vrr=sqrt(dot_product(pot_on_plane(i,j)%r-epe(k)%r, &
                  pot_on_plane(i,j)%r-epe(k)%r))
             if(vrr <= 0.001_r8_kind) vrr=0.001_r8_kind
             vrs=sqrt(dot_product(pot_on_plane(i,j)%r-epe(k)%s, &
                  pot_on_plane(i,j)%r-epe(k)%s))
             if(vrs <= 0.001_r8_kind) vrs=0.001_r8_kind
             vrc=sqrt(dot_product(pot_on_plane(i,j)%r-epe(k)%c, &
                  pot_on_plane(i,j)%r-epe(k)%c))
             if(vrc <= 0.001_r8_kind) vrc=0.001_r8_kind

             V_pot=V_pot+(epe(k)%qs/vrs+epe(k)%qc/vrc- &
                  epe(k)%q*erf(ERROR_FUNCTION_PARAMETER*vrr)/vrr)
          end do reg_2a_nions
          V_pot=V_pot*qau_qepe
          if(V_abs_limit > 0.0_r8_kind .and. abs(V_pot) > V_abs_limit) then
             if(V_pot < 0.0_r8_kind) V_pot = -V_abs_limit
             if(V_pot > 0.0_r8_kind) V_pot = V_abs_limit
          end if
          if(index(output_units,"a.u.")/=0) V_pot=V_pot/27.211652_r8_kind
          if(index(output_format,"scilab")/=0 .or. &
               index(output_format,"SCILAB")/=0) then
             write(vunit,'(f25.15)') V_pot
          else if(index(output_format,"worksheet")/=0 .or. &
               index(output_format,"WORKSHEET")/=0 .or. &
               index(output_format,"gnuplot")/=0 .or. &
               index(output_format,"GNUPLOT")/=0) then
             write(vunit,'(3f25.15)') x,y,V_pot
             if((index(output_format,"gnuplot")/=0 .or. &
                  index(output_format,"GNUPLOT")/=0) .and. j==nygrid) write(vunit,'(a)') ""
          endif
       end do nygrd
    end do nxgrd

    call returnclose_iounit(vunit)

    deallocate(pot_on_plane,stat=istat)
    if(istat /= 0) call error_handler( &
         "epe_pot_module: calc_pot: deallocate pot_on_plane failed")

  end subroutine calc_pot

end module epe_pot_module

