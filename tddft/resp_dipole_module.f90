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
module resp_dipole_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  ! 
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

#include "def.h"
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  public resp_dipole_rewrite

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine resp_dipole_rewrite()
    use symmetry_data_module, only: symmetry_data_dip_irrep_mult, &
         symmetry_data_dip_components, &
         symmetry_data_dimension
    use clebsch_gordan, only: cg=>cg_eliminated, prod_bas
    use symmetry_data_module, only: symmetry_data_n_irreps, & 
         symmetry_data_n_partners, &
         symmetry_data_dimension, &
         symmetry_data_n_spin
    use constants, only: zero
    use resp_util_module, only: resp_util_borders, resp_util_calc_ou,&
         resp_util_fname, min_diff
    use iounitadmin_module,   only: openget_iounit,returnclose_iounit,&
         write_to_output_units
    USE filename_module, ONLY: tmpfile, resp_dir
    USE datatype, only: arrmat1, arrmat2
    USE debug
    use io
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind)      :: tm_xyz(3), dip, tma(3),tmb(3)
    real(kind=r8_kind)      :: coeff_abc, coeff_xyz
    integer(kind=i4_kind)   :: i_ir_c, nirr, i_pa_c, mult, isym, i_xyz
    integer(kind=i4_kind)   :: i_ir_a, i_pa_a, i_ir_b, i_pa_b
    integer(kind=i4_kind)   :: tmp_io_unit,io_stat,alloc_stat
    integer(kind=i4_kind)   :: na, nb, npaa, npab, npac
    integer(kind=i4_kind)   :: ispin, ia, ib, n_spin, i_dip_mlt
    integer(kind=i4_kind)   :: occs, occe, unoccs, unocce
    integer(kind=i4_kind)   :: i_mlt, dim2, is_count
    type(prod_bas),pointer  :: pcg    

    real(kind=r8_kind),allocatable :: dip_int(:,:)
    integer(i4_kind),  allocatable :: dm(:,:)

    type exponents3D
       real(kind = r8_kind),pointer :: expn3D(:,:,:,:)
!       real(kind = r8_kind),allocatable :: expn3D(:,:,:,:)
    end type exponents3D

    type partner2
       type(exponents3D),    pointer :: prtn(:,:)
    end type partner2

    type(partner2),          allocatable :: dip_matrix(:,:)

    integer(i4_kind), parameter:: X = 1, Y = 2, Z = 3

    !------------ Executable code --------------------------------

    call write_to_output_units('response_rewrite_dipoletape: preparations')

    nirr   = symmetry_data_n_irreps()
    n_spin = symmetry_data_n_spin()

    ALLOCATE(dm(nirr,n_spin),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    do ispin = 1, n_spin
       do i_ir_c = 1, nirr
          npac = symmetry_data_n_partners(i_ir_c)
          call resp_util_calc_ou(i_ir_c,ispin,dim2)
          dm(i_ir_c,ispin) = dim2
       end do
    end do

    !! dipole matrix dip_matrix(ira,irb)%partn(paa,pab)%expn3D(ea,eb,spin,3)
    allocate(dip_matrix(nirr,nirr),STAT = alloc_stat)
    ASSERT(alloc_stat==0)

    do i_ir_a = 1, nirr
       npaa = symmetry_data_n_partners(i_ir_a)
       do i_ir_b = 1, nirr
          npab = symmetry_data_n_partners(i_ir_b)

          allocate(dip_matrix(i_ir_a,i_ir_b)%prtn(npaa,npab),STAT=alloc_stat)
          ASSERT(alloc_stat==0)

          do i_pa_a = 1, npaa 
             do i_pa_b = 1, npab 
                na = symmetry_data_dimension(i_ir_a)
                nb = symmetry_data_dimension(i_ir_b)                
                allocate(dip_matrix(i_ir_a,i_ir_b)%&
                     prtn(i_pa_a,i_pa_b)%expn3D(na,nb,n_spin,3),STAT = alloc_stat)
                ASSERT(alloc_stat == 0)
                
             end do
          end do
       end do
    end do

    ! open the temporary file written in dipole_module
    tmp_io_unit=openget_iounit(status='unknown',form='unformatted',&
         & file=trim(tmpfile("resp_dipoles_tmp.dat")))

    do !while (.not.eof(tmp_io_unit))
       read(tmp_io_unit,iostat=io_stat) &
            ispin, & 
            i_ir_a, i_pa_a, ia,&
            i_ir_b, i_pa_b, ib,&
            (tm_xyz(i_xyz), i_xyz=1,3)
       if( io_stat == -1 ) EXIT ! the loop at EOF!
       ASSERT(io_stat==0)

       dip_matrix(i_ir_a,i_ir_b)%prtn(i_pa_a,i_pa_b)%&
            expn3D(ia,ib,ispin,1:3) = tm_xyz(1:3)

       dip_matrix(i_ir_b,i_ir_a)%prtn(i_pa_b,i_pa_a)%&
            expn3D(ib,ia,ispin,1:3) = tm_xyz(1:3)

    end do

    do ispin =1,n_spin
       do i_ir_a = 1, nirr
          do i_ir_b = 1, nirr
             do i_pa_a = 1, symmetry_data_n_partners(i_ir_a)
                do i_pa_b = 1, symmetry_data_n_partners(i_ir_b)
                   do ia = 1,  symmetry_data_dimension(i_ir_a)
                      do ib = 1, symmetry_data_dimension(i_ir_b)
                         tma = dip_matrix(i_ir_a,i_ir_b)%prtn(i_pa_a,i_pa_b)%&
                              expn3D(ia,ib,ispin,:)
                         tmb = dip_matrix(i_ir_b,i_ir_a)%prtn(i_pa_b,i_pa_a)%&
                              expn3D(ib,ia,ispin,:)
                         do i_xyz = X ,Z
         ASSERT(tma(i_xyz) == tmb(i_xyz))
                         end do

                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    close(tmp_io_unit,iostat=io_stat)
    ASSERT(io_stat==0)    

    ispin_: do ispin = 1, n_spin
       isym = 0
       i_ir_c_: do i_ir_c = 1, nirr
          npac = symmetry_data_n_partners(i_ir_c)
          allocate(dip_int(dm(i_ir_c,ispin),3),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          dip_int = zero
!         io_unit = openget_iounit(trim(resp_dir)//'/'//&
!              resp_util_fname('dipol',i_ir_c,ispin),form='unformatted',status='unknown')          

          i_dip_mlt_: do i_dip_mlt = 1,symmetry_data_dip_irrep_mult(i_ir_c)
!!$             print *," i_dip_mlt = ",i_dip_mlt
             is_count = 0

             i_ir_a_: do i_ir_a = 1, nirr
                npaa = symmetry_data_n_partners(i_ir_a)
                i_ir_b_: do i_ir_b = 1, nirr
!!$                   print *,"i_ir_c = ",i_ir_c," i_ir_a = ",i_ir_a," i_ir_b = ",i_ir_b
                   npab = symmetry_data_n_partners(i_ir_b)

                   occs = 0
                   occe = 0
                   unoccs = 0
                   unocce = 0
                   call resp_util_borders(i_ir_a,i_ir_b,ispin,occs,occe,unoccs,unocce)

                   if (occs == 0) then 
                      ASSERT(occe==0)
                      cycle
                   end if

                   if (unoccs == 0) then
                      ASSERT(unocce==0)
                      cycle
                   end if

                   mult = cg(i_ir_c,i_ir_a,i_ir_b)%mult
                   i_mlt_: do i_mlt = 1, mult
!!$                      print *,"i_mlt = ",i_mlt
                      pcg =>   cg(i_ir_c,i_ir_a,i_ir_b)%sub(i_mlt)

                      ia_: do ia = occs, occe
                         ib_: do ib = unoccs,unocce
                            is_count = is_count + 1

!!$                            print *,"ia =",ia," ib ",ib," is_count = ",is_count
                            !! HERE: IN: ( A -> B ) C, i_mlt, ia->ib
                            !!      OUT: < (A->B)C, ia,ib || rC > -> dip_int(is_count)

                            dip = zero
                             i_pa_c_: do i_pa_c = 1, npac
                               isym = isym_calc(i_ir_c,i_pa_c,i_dip_mlt)!!isym + 1
!!$                               print *,"isym = ",isym
                               i_pa_a_: do i_pa_a = 1, npaa
                                  i_pa_b_: do i_pa_b = 1, npab
!!$                                     print *,"i_pa_c = ",i_pa_c," i_pa_a = ",i_pa_a," i_pa_b = ",i_pa_b
                                     coeff_abc = pcg%c(i_pa_c,i_pa_a,i_pa_b)
                                     tm_xyz(X:Z) = dip_matrix(i_ir_a,i_ir_b)%prtn(i_pa_a,i_pa_b)%expn3D(ia,ib,ispin,X:Z)
!!$                                     print *," tm_xyz = ",tm_xyz
                                     do i_xyz = X,Z
                                        coeff_xyz = symmetry_data_dip_components(isym,i_xyz)
                                        dip_int(is_count,isym) = dip_int(is_count,isym) &
                                             &                 + coeff_abc * coeff_xyz * tm_xyz(i_xyz)
                                     end do
!!$                                     print *,"dip_int(",is_count,isym,") = ",dip_int(is_count,isym)
                                  end do i_pa_b_
                               end do i_pa_a_
                            end do i_pa_c_
                         end do ib_
                      end do ia_
                   end do i_mlt_
                end do i_ir_b_
             end do i_ir_a_
          end do i_dip_mlt_
!         io_unit = openget_iounit(trim(resp_dir)//'/'//&
!              resp_util_fname('dipol',i_ir_c,ispin),form='unformatted',status='unknown')          
          call write_buffer( trim(resp_dir)//'/'//resp_util_fname('dipol',i_ir_c,ispin) &
                           , dip_int &
                           )
!         call returnclose_iounit(io_unit)
          deallocate(dip_int,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do i_ir_c_
    end do ispin_

    do i_ir_a = 1, nirr
       npaa = symmetry_data_n_partners(i_ir_a)
       do i_ir_b = 1, nirr
          npab = symmetry_data_n_partners(i_ir_b)
          do i_pa_a = 1, npaa 
             do i_pa_b = 1, npab 

                deallocate(dip_matrix(i_ir_a,i_ir_b)%prtn(i_pa_a,i_pa_b)%&
                     expn3D,STAT = alloc_stat)
                ASSERT(alloc_stat==0)

             end do
          end do
          deallocate(dip_matrix(i_ir_a,i_ir_b)%prtn,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
    end do

    deallocate(dip_matrix,dm,STAT = alloc_stat)
    ASSERT(alloc_stat==0)

  end subroutine resp_dipole_rewrite

  integer(i4_kind) function isym_calc(ir,pa,mlt)
    use symmetry_data_module, only: symmetry_data_dip_irrep_mult,&
         symmetry_data_n_partners
    integer(i4_kind), intent(in) :: ir,pa,mlt
    !===================================================
    integer(i4_kind) :: i, p, m
    
    isym_calc  = 0
    ir_: do i = 1, ir-1
       pa_: do p = 1,symmetry_data_n_partners(i) 
          mlt_: do m = 1,symmetry_data_dip_irrep_mult(i)  
             isym_calc = isym_calc + 1
          end do mlt_
       end do pa_
    end do ir_

    do p = 1,pa 
       do m = 1,mlt  
          isym_calc = isym_calc + 1
       end do
    end do

    ASSERT(isym_calc<4)

  end function isym_calc
  !*************************************************************
  

  !--------------- End of module ----------------------------------
end module resp_dipole_module
