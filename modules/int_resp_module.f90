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
# include "def.h"
!===============================================================
! Public interface of module
!===============================================================
module int_resp_module
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

  use type_module ! type specification parameters
  use comm_module  
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public int_resp_Clb_3c

  !================================================================
  ! End of public interface of module
  !================================================================

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine int_resp_Clb_3c
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use ch_response_module,   only: dimension_of_fit_ch
    use int_send_3c_resp,     only: int_send_3c_resp_read
    use symmetry_data_module, only: symmetry_data_n_irreps
    use clebsch_gordan,       only: cg=>cg_eliminated
!   use iounitadmin_module,   only: openget_iounit,returnclose_iounit
    use eigen_data_module,    only: eigvec
    use symmetry_data_module
    use resp_util_module
    use filename_module
    use msgtag_module
    use debug
    use xpack
    use io
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)              :: n_ir, i_ir_c, i_ir_a, i_ir_b, na, nb, nc
    integer(kind=i4_kind)              :: i_spin
    integer(kind=i4_kind)              :: is_counter
    integer(kind=i4_kind)              :: n_procs!, io_unit 
    integer(kind=i4_kind)              :: k_index
    integer(kind=i4_kind)              :: occ_times_unocc_dim, status
    real(kind=r8_kind),    allocatable :: coulomb_3c_matrix(:,:,:)
    real(kind=r8_kind),    allocatable :: final_coulomb(:,:)
    integer(kind=i4_kind), allocatable :: nff(:)
    integer(kind=i4_kind)              :: imult, mult_irr
    integer(kind=i4_kind)              :: i_proc
    integer(kind=i4_kind)              :: dim_ma, dim_sl

    integer(kind=i4_kind)              :: DM1, iac, NS, NF
    real(kind=r8_kind),    allocatable :: RM(:,:)

    integer(kind=i4_kind)              :: n_spin, dim_factor
    
    !------------ Executable code --------------------------------
    if(comm_i_am_master()) then
       if(comm_parallel()) then 
          call comm_init_send(comm_all_other_hosts,msgtag_response_3Clb_start)
          call comm_send()
       end if
    end if

    n_ir       = ssym%n_irrep
    n_spin     = ssym%n_spin
    dim_factor = 2

    n_procs = comm_get_n_processors()

    allocate(nff(n_procs),STAT = status)
    ASSERT(status==0)

    i_spin_: do i_spin = 1, n_spin

       i_ir_c_: do i_ir_c = 1, n_ir

          NS = 0
          NF = 0

          call resp_util_calc_ou(i_ir_c,i_spin,DM1)          

          dim_sl = INT(dimension_of_fit_ch(i_ir_c)/n_procs)
          dim_ma = dimension_of_fit_ch(i_ir_c) &
               & - (n_procs-1) * dim_sl

!!$          nff(0) = 0
          nff(1) = dim_ma
          do i_proc = 2, n_procs
             nff(i_proc) = dim_sl
          end do

          nc = nff(comm_myindex())
!!          if (nc .eq. 0) exit

          if ( comm_i_am_master()) then
!             io_unit = openget_iounit(trim(resp_dir)//'/'//&
!                resp_util_fname('co_3c',i_ir_c,i_spin),form='unformatted',status='unknown')
          ALLOCATE(RM(DM1,dimension_of_fit_ch(i_ir_c)),STAT = status)
          ASSERT(status==0)
          end if
          i_ir_a_: do i_ir_a = 1, n_ir
             na = ssym%dim(i_ir_a)
             i_ir_b_: do i_ir_b = 1, n_ir
                nb = ssym%dim(i_ir_b)

                ! get the number of possible transitions
                call resp_util_calc_transitions_v2(i_ir_a, i_ir_b, i_ir_c, &
                     i_spin, occ_times_unocc_dim)

                if (occ_times_unocc_dim == 0) cycle 

                allocate(coulomb_3c_matrix(nc,nb,na), STAT = status)
                ASSERT(status==0)
                coulomb_3c_matrix = 0.0_r8_kind

!!$                allocate(final_coulomb(occ_times_unocc_dim,nc),STAT = status)
                allocate(final_coulomb(occ_times_unocc_dim,MAXVAL(nff)),STAT = status)
                ASSERT(status==0)

                mult_irr = cg(i_ir_c,i_ir_a,i_ir_b)%mult
                
                imult_: do imult = 1, mult_irr

                   call int_send_3c_resp_read(i_ir_c, i_ir_a, i_ir_b, imult, coulomb_3c_matrix, nff)

                   fit_charge_index_: do k_index = 1, nc

                      call coulomb_mult_eigvec(i_ir_a,i_ir_b, &
                           i_spin, na, nb, coulomb_3c_matrix(k_index,:,:), final_coulomb(:,k_index) )

                   end do fit_charge_index_

                   ! ******
                   ! collect results on master and save them to tape
                   ! ******
                   if(comm_i_am_master()) then
                      NS = NF + 1
                      NF = NF + size(final_coulomb,1)
                      RM(NS:NF,1:nc) = final_coulomb(:,1:nc)
                      iac = 1
                      do i_proc=2,n_procs
                         call comm_save_recv(i_proc,msgtag_response_3Clb_send)
!!$                         call upck(final_coulomb(:,1:nff(i_proc)))
                         call upck(final_coulomb)
                         iac = iac + nff(i_proc-1)
                         RM(NS:NF,iac:iac+nff(i_proc)-1) = final_coulomb(:,1:nff(i_proc))
                      end do

#if 0
                      where(abs(RM)<1.0D-16_r8_kind) RM = 0.0_r8_kind
                      print *,"=== irc, ira, irb = ",i_ir_c, i_ir_a, i_ir_b, " RM === "
                      call octave("RM",RM(NS:NF,:))
#endif

                   else                         
                      call comm_init_send(comm_master_host,msgtag_response_3Clb_send)
                      call pck(final_coulomb)
                      call comm_send()
                   end if

                end do imult_

                ! deallocate memory for final result for this spin and irrep
                deallocate(final_coulomb,STAT = status)
                ASSERT(status==0)

                deallocate(coulomb_3c_matrix,STAT=status)
                ASSERT(status==0)

             end do i_ir_b_
          end do i_ir_a_
          if (comm_i_am_master()) then

#if 0
             where(abs(RM)<1.0D-16) RM = 0.0_r8_kind
             print *,"=== i_ir_c = ",i_ir_c, " RM === "
             call octave("RM",RM)
#endif
!             io_unit = openget_iounit(trim(resp_dir)//'/'//&
!                resp_util_fname('co_3c',i_ir_c,i_spin),form='unformatted',status='unknown')
             call write_buffer( trim(resp_dir)//'/'//resp_util_fname('co_3c',i_ir_c,i_spin) &
                              , RM                                                           &
                              )
             DEALLOCATE(RM,STAT=status)
             ASSERT(status==0)
!            call returnclose_iounit(io_unit)
          end if
       end do i_ir_c_

    end do i_spin_

    deallocate(nff,STAT = status)
    ASSERT(status==0)

  contains

    subroutine coulomb_mult_eigvec(&
         i_ir_a,i_ir_b,    &
         i_spin, na, nb, coulomb_3c, coulomb_matrix)
      use constants,         only: zero
      use resp_util_module,  only: min_diff
      
      use linalg_module,     only: matmatmul
      implicit none

      integer(kind=i4_kind), intent(in)    :: i_ir_a, i_ir_b, i_spin 
      integer(kind=i4_kind), intent(in)    :: na, nb
      real(kind=r8_kind),    intent(inout) :: coulomb_3c(:,:)
      real(kind=r8_kind),    intent(inout) :: coulomb_matrix(:)
      ! --- local variables
      real(kind=r8_kind),allocatable :: aux_matrix(:)
      integer(kind=i4_kind)      :: i_occ
      integer(kind=i4_kind)      :: i_unocc, status
      integer(kind=i4_kind)      :: occs, occe, unoccs, unocce
      real(kind=r8_kind),pointer :: trafo(:)

      call resp_util_borders(i_ir_a,i_ir_b,i_spin,occs,occe,unoccs,unocce)
      ! allocate auxiliary matrix
      allocate(aux_matrix(nb), stat=status)
      ASSERT(status==0)

      aux_matrix = 0.0_r8_kind
      
      is_counter = 0   ! combined metaindex a->b
      
      i_occ_: do i_occ = occs, occe
         trafo=>eigvec(i_ir_a)%m(:,i_occ,i_spin)
         aux_matrix = matmul(coulomb_3c,trafo)
         i_unocc_: do i_unocc = unoccs, unocce
            trafo=>eigvec(i_ir_b)%m(:,i_unocc,i_spin)
            is_counter=is_counter + 1
            coulomb_matrix(is_counter) = dot_product(trafo, aux_matrix)
         end do i_unocc_
      end do i_occ_
      
      ! deallocate matrices for intermediate results
      deallocate(aux_matrix,  stat=status)
      ASSERT(status==0)

    end subroutine coulomb_mult_eigvec

  end subroutine int_resp_Clb_3c
  !*************************************************************


  !--------------- End of module ----------------------------------
end module int_resp_module
