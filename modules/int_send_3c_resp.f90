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
!===============================================================
! Public interface of module
!===============================================================
module  int_send_3c_resp
  !---------------------------------------------------------------
  ! 
  !  Purpose: Sending, receiving and storing of 
  !           3 center coulomb integrals.
  !
  !  Module called by: main_slave, integral_main_2cob3c, 
  !    integral_calc_quad_2cob3c, integral_interrupt_2cob3c,
  !    integral_setup_2cob3c, integral_shutdown_2cob3c
  !
  !  Author: SB
  !  Date: 3/05
  !
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
# include "def.h"
  use type_module
  use filename_module, only: tmpfile
  use symmetry_data_module, only: symmetry_data_n_irreps
  use readwriteblocked_module, only : readwriteblocked_startwrite, &
       readwriteblocked_write, &
       readwriteblocked_stopwrite, &
       readwriteblocked_tapehandle, &
       readwriteblocked_startread, &
       readwriteblocked_stopread, &
       readwriteblocked_read
  use int_data_2cob3c_module, only : symadapt_totsym_3c_int_type, multiplicity
  use ch_response_module, only: dimension_of_fit_ch, orb_position
  use unique_atom_module, only: unique_atoms, N_unique_atoms
  use clebsch_gordan, only: cg=>cg_eliminated
  use debug, only: show 

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private

  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public int_send_3c_resp_save, int_send_3c_resp_read

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine int_send_3c_resp_save(ua1, ua2, l1, l2, sa_3co_resp)
    !------------ Modules used ---------------------------------
    implicit none

    integer(kind=i4_kind), intent(in) :: ua1, ua2, l1, l2
    type(multiplicity), intent(in) :: sa_3co_resp(:,:,:) ! (n_irr, n_irr, n_irr)

    !------------ Declaration of local variables -----------------
    type(readwriteblocked_tapehandle)          :: th_3c_resp
    real(kind=r8_kind), pointer                :: work_int(:,:,:,:,:)
    integer(kind=i4_kind)  :: i_ir_c, i_ir_a, i_ir_b, ia, ib, ic
    integer(kind=i4_kind)  :: n_ifa, n_ifb, count, ifa, ifb, n_ir
    integer(kind=i4_kind)  :: nia, nib, dim, mult_irrep, status, imult
    real(kind=r8_kind), allocatable            :: tmp_matrix(:)
    !------------ Executable code --------------------------------

    n_ir = symmetry_data_n_irreps()
    i_ir_c_: do i_ir_c = 1, n_ir
       if (dimension_of_fit_ch(i_ir_c)==0) cycle
       i_ir_a_: do i_ir_a = 1, n_ir
          n_ifa = unique_atoms(ua1)%symadapt_partner(i_ir_a,l1)%N_independent_fcts !!$size(work_int,5)
          if (n_ifa .eq.0) cycle
          i_ir_b_: do i_ir_b = 1, n_ir
             n_ifb = unique_atoms(ua2)%symadapt_partner(i_ir_b,l2)%N_independent_fcts !!$size(work_int,4)
             if (n_ifb .eq. 0) cycle
             mult_irrep = cg(i_ir_c,i_ir_a,i_ir_b)%mult

             if (mult_irrep .eq. 0) cycle

             i_mlt_: do imult = 1, mult_irrep                

                work_int => sa_3co_resp(i_ir_c, i_ir_a, i_ir_b)%mult(imult)%int

                nib   = size(work_int,1)
                nia   = size(work_int,2)
                dim = (n_ifa * nia) * (n_ifb * nib)

                allocate(tmp_matrix(dim), STAT = status)
                ASSERT(status==0)
                tmp_matrix = 0.0_r8_kind

                !DBG:                print *,"file: ",fname(ua1,ua2,l1,l2,i_ir_c,i_ir_a,i_ir_b,imult),&
                !DBG:                     "will be saved"

                call readwriteblocked_startwrite(trim(tmpfile(fname(ua1,ua2,l1,l2,i_ir_c,i_ir_a,i_ir_b,imult))), th_3c_resp)

                i_fit_ch_: do ic = 1, dimension_of_fit_ch(i_ir_c)

                   count = 0

                   ifa_: do ifa = 1, n_ifa
                      ia_: do ia = 1, nia

                         ifb_: do ifb = 1, n_ifb
                            ib_: do ib =1, nib

                               count = count + 1
                               tmp_matrix (count) = work_int(ib, ia, ic, ifb, ifa)
                               !! if (tmp_matrix(count) > 1.0E+20) print *,"WARNING!!! tmp_matrix > 1.0E+20"

                            end do ib_
                         end do ifb_
                      end do ia_
                   end do ifa_

                   call readwriteblocked_write(tmp_matrix,th_3c_resp)

                end do i_fit_ch_

                call readwriteblocked_stopwrite(th_3c_resp)

                tmp_matrix = 0.0_r8_kind

                call readwriteblocked_startwrite(trim(tmpfile(fname(ua2,ua1,l2,l1,i_ir_c,i_ir_b,i_ir_a,imult))), th_3c_resp) 

                do ic = 1, dimension_of_fit_ch(i_ir_c)

                   count = 0

                   ifbT_: do ifb = 1, n_ifb
                      ibT_: do ib =1, nib

                         ifaT_: do ifa = 1, n_ifa
                            iaT_: do ia = 1, nia

                               count = count + 1
                               tmp_matrix(count) = work_int(ib, ia, ic, ifb, ifa)

                            end do iaT_
                         end do ifaT_
                      end do ibT_
                   end do ifbT_

                   call readwriteblocked_write(tmp_matrix,th_3c_resp)

                end do

                call readwriteblocked_stopwrite(th_3c_resp)

                deallocate(tmp_matrix, STAT=status)
                ASSERT(status==0)


             end do i_mlt_

          end do i_ir_b_
       end do i_ir_a_
    end do i_ir_c_



  end subroutine int_send_3c_resp_save
  !*************************************************************


  !*************************************************************
  subroutine int_send_3c_resp_read(i_ir_c, i_ir_a, i_ir_b, imult, &
       coulomb_3c_matrix, nff)
    !------------ Modules used ---------------------------------
    use symmetry_data_module, only: symmetry_data_n_partners, symmetry_data_n_irreps
    use comm_module
    use msgtag_module, only: msgtag_tmp_3co_send,msgtag_packed_message
    use debug
    use xpack
    use comm, only: comm_reduce, comm_bcast
    implicit none
    integer(kind=i4_kind), intent(in)    :: i_ir_c, i_ir_a, i_ir_b, imult, nff(:)
    real(kind=r8_kind),    intent(inout) :: coulomb_3c_matrix(:,:,:)
    !------------ Declaration of local types ---------------------
    integer(kind=i4_kind) :: ua1, ua2
    integer(kind=i4_kind) :: l1, l2, n_proc
    integer(kind=i4_kind) :: ic 
    integer(kind=i4_kind) :: status
    integer(kind=i4_kind) :: n_if_a, n_contr_a, n_uncontr_a, dim_ch_a
    integer(kind=i4_kind) :: n_if_b, n_contr_b, n_uncontr_b, dim_ch_b, i_proc
    integer(kind=i4_kind) :: ff, myindex
    type(readwriteblocked_tapehandle)          :: th_3c_resp
    real(kind=r8_kind), allocatable            :: tmp_matrix(:,:)
    logical                                    :: check_exist
    integer(i4_kind)                           :: ce ! (int) check_exists

    integer(kind=i4_kind) :: ira,irb,ua,ub,la,lb,dima,dimb
    !------------ Executable code --------------------------------

    n_proc = comm_get_n_processors()
    myindex = comm_myindex()
    ff = nff(myindex)

    ua2_: do ua2 = 1, N_unique_atoms         
       ua1_: do ua1 = 1, N_unique_atoms

          l2_: do l2 = 0, unique_atoms(ua2)%lmax_ob
             n_if_b      = unique_atoms(ua2)%symadapt_partner(i_ir_b,l2)%N_independent_fcts
             n_contr_b   = unique_atoms(ua2)%l_ob(l2)%N_contracted_fcts
             n_uncontr_b = unique_atoms(ua2)%l_ob(l2)%N_uncontracted_fcts
             dim_ch_b    = (n_contr_b + n_uncontr_b) * n_if_b
             if (dim_ch_b == 0) cycle

             l1_: do l1 = 0, unique_atoms(ua1)%lmax_ob
                n_if_a      = unique_atoms(ua1)%symadapt_partner(i_ir_a,l1)%N_independent_fcts 
                n_contr_a   = unique_atoms(ua1)%l_ob(l1)%N_contracted_fcts
                n_uncontr_a = unique_atoms(ua1)%l_ob(l1)%N_uncontracted_fcts
                dim_ch_a    = (n_contr_a + n_uncontr_a) * n_if_a
                if (dim_ch_a == 0) cycle

                inquire (file = trim(tmpfile(fname(ua1,ua2,l1,l2,i_ir_c,i_ir_a,i_ir_b,imult))), exist=check_exist)

                !! TEST IF THE FILE EXIST SOMEWHERE
                ce = 0
                if( check_exist ) ce = 1
                call comm_reduce(ce)
                call comm_bcast(ce)
                if ( ce == 0 ) cycle

                file_found_: if (check_exist) then

                   allocate(tmp_matrix(MAXVAL(nff),(dim_ch_a*dim_ch_b)), STAT = status)
                   ASSERT(status==0)

                   call readwriteblocked_startread(trim(tmpfile(fname(ua1,ua2,l1,l2,i_ir_c,i_ir_a,i_ir_b,imult))), th_3c_resp)

                   proc_cycle_: do i_proc = 1, n_proc
                      tmp_matrix = 0.0_r8_kind
                      do ic = 1, nff(i_proc)
                         call readwriteblocked_read(tmp_matrix(ic,:),th_3c_resp)
                      end do

                      if (myindex==i_proc) then
                         call make_3c(i_ir_a, ua1, l1, &
                              &       i_ir_b, ua2, l2, &
                              &       dim_ch_a, dim_ch_b, nff(i_proc), &
                              tmp_matrix,coulomb_3c_matrix)
                         cycle
                      end if

                      !! FIXME: I WILL THINK ABOUT IT
                      if (.not. comm_parallel()) cycle
                      call comm_init_send(i_proc,msgtag_tmp_3co_send)
                      call pck(i_ir_a)
                      call pck(i_ir_b)
                      call pck(ua1)
                      call pck(ua2)
                      call pck(l1)
                      call pck(l2)
                      call pck(dim_ch_a)
                      call pck(dim_ch_b)
                      call pck(tmp_matrix)                      
                      call comm_send()

                   end do proc_cycle_
                   call readwriteblocked_stopread(th_3c_resp)

                   deallocate(tmp_matrix,STAT = status)
                   ASSERT(status==0)

                else !! file_not_found_
                   if (.not. comm_parallel()) cycle
                   call comm_save_recv(comm_all_other_hosts,msgtag_tmp_3co_send)
                   call upck(ira)
                   call upck(irb)
                   call upck(ua)
                   call upck(ub)
                   call upck(la)
                   call upck(lb)
                   call upck(dima)
                   call upck(dimb)
                   allocate(tmp_matrix(MAXVAL(nff),(dima*dimb)),STAT=status)
                   ASSERT(status==0)
                   call upck(tmp_matrix)

                   call make_3c(ira, ua, la, &
                        &       irb, ub, lb, &
                        &       dima, dimb, ff, &
                        tmp_matrix,coulomb_3c_matrix)

                   deallocate(tmp_matrix,STAT=status)
                   ASSERT(status==0)

                end if file_found_

             end do l1_
          end do l2_

       end do ua1_
    end do ua2_

    ff = 1979
    call comm_bcast(ff)  !! needs to be here, because of synchronization

  contains

    subroutine make_3c(ira, ua, la, &
         &             irb, ub, lb, &
         &             da , db, ff, &
         &   from,to)
      integer(kind = i4_kind), intent(in   ) :: ira,irb,ua,ub,la,lb
      integer(kind = i4_kind), intent(in   ) :: da, db, ff
      real   (kind = r8_kind), intent(in   ) :: from(:,:)
      real   (kind = r8_kind), intent(inout) ::   to(:,:,:)
      !------------ Declaration of local types ---------------------
      integer(kind = i4_kind)                :: counter, ca, cb, ia,ib

      counter   = 0_i4_kind
      ca = orb_position(ira, ua, la-1)
      expa: do ia = 1, da
         ca = ca + 1
         cb = orb_position(irb, ub, lb-1)
         expb: do ib = 1, db
            cb = cb + 1
            counter = counter + 1
            to(:ff,cb,ca) = from(:ff,counter)
         end do expb
      end do expa
    end subroutine make_3c

  end subroutine int_send_3c_resp_read
  !*************************************************************

  character(len=45) function fname(ua1,ua2,l1,l2,i_ir_c,i_ir_a,i_ir_b,imult)
    implicit none
    integer(kind = i4_kind), intent(in) :: ua1,ua2,l1,l2
    integer(kind = i4_kind), intent(in) :: i_ir_c,i_ir_a,i_ir_b,imult
    !------------ Declaration of local types ---------------------
    character(len=4) :: ua1_char,ua2_char
    character(len=4) :: l1_char,l2_char
    character(len=4) :: irc_char,ira_char,irb_char
    character(len=4) :: imlt_char
    !------------ Executable code --------------------------------

    write (ua1_char, '(i4)') ua1
    write (ua2_char, '(i4)') ua2
    write (l1_char,  '(i4)') l1
    write (l2_char,  '(i4)') l2
    write (irc_char, '(i4)') i_ir_c
    write (ira_char, '(i4)') i_ir_a
    write (irb_char, '(i4)') i_ir_b
    write (imlt_char,'(i4)') imult
    ua1_char  = adjustl(ua1_char)
    ua2_char  = adjustl(ua2_char)
    l1_char   = adjustl(l1_char)
    l2_char   = adjustl(l2_char)
    irc_char  = adjustl(irc_char)
    ira_char  = adjustl(ira_char)
    irb_char  = adjustl(irb_char)
    imlt_char = adjustl(imlt_char)
    fname = 'resp_3c_'//trim(ua1_char)//"."//trim(ua2_char)//"."//trim(l1_char)//"."//trim(l2_char)//&
         '_'//trim(irc_char)//"."//trim(ira_char)//"."//trim(irb_char)//'_'//trim(imlt_char)

  end function fname

  !--------------- End of module ----------------------------------
end module int_send_3c_resp
