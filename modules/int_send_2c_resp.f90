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
!=====================================================================
! Public interface of module
!=====================================================================
module  int_send_2c_resp
!---------------------------------------------------------------
!
!  Purpose: Sending, receiving and storing of
!           2 center coulomb integrals.
!
!  Module called by: integral_main_2cff,
!
!  Based on: int_send_3c_resp
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
#include <def.h>
use type_module
use filename_module, only : tmpfile, resp_dir
use symmetry_data_module, only: symmetry_data_n_irreps
use readwriteblocked_module, only : readwriteblocked_startwrite, &
     readwriteblocked_write, &
     readwriteblocked_stopwrite, &
     readwriteblocked_tapehandle, &
     readwriteblocked_startread, &
     readwriteblocked_stopread, &
     readwriteblocked_read
use ch_response_module, only: dimension_of_fit_ch
use unique_atom_module, only: unique_atoms, N_unique_atoms

implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private

!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
public int_send_2c_resp_save, int_send_2c_resp_rewrite

  !===================================================================
  ! End of public interface of module
  !===================================================================

!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine int_send_2c_resp_save(ua1,ua2,l1,l2,sa_2co_resp,i_ir)
    use debug
    !------------ Modules used ---------------------------------
    implicit none

    integer(kind=i4_kind), intent(in) :: ua1, ua2, l1, l2
    real(kind=r8_kind),    intent(in) :: sa_2co_resp(:,:,:,:)

    !------------ Declaration of local variables ---------------------
    type(readwriteblocked_tapehandle)          :: th_2c_resp
    integer(kind=i4_kind)  :: i_ir, ia, ib
    integer(kind=i4_kind)  :: n_ifa, n_ifb, count, ifa, ifb, n_ir
    integer(kind=i4_kind)  :: nia, nib, dim, status
    real(kind=r8_kind), allocatable            :: tmp_matrix(:)
    !------------ Executable code ------------------------------------

    !! FIXME: if slave send to master, if not send to tape

    n_ir = symmetry_data_n_irreps()
    n_ifb = size(sa_2co_resp,3) !!$ua1%symadapt_partner(i_ir_a,l1)%N_independent_fcts
    n_ifa = size(sa_2co_resp,4) !!$ua2%symadapt_partner(i_ir_b,l2)%N_independent_fcts

    nib   = size(sa_2co_resp,1)
    nia   = size(sa_2co_resp,2)
    dim = (n_ifa * nia) * (n_ifb * nib)

    if (dim==0) return

    call readwriteblocked_startwrite(trim(tmpfile(fname(ua1,ua2,l1,l2,i_ir))), th_2c_resp)

    count = 0
    allocate(tmp_matrix(dim),STAT = status)
    ASSERT(status==0)

    ib_: do ib =1, nib
       ifb_: do ifb = 1, n_ifb

          ia_: do ia = 1, nia
             ifa_: do ifa = 1, n_ifa
                count = count + 1
                tmp_matrix(count) = sa_2co_resp(ib,ia,ifb,ifa)
             end do ifa_
          end do ia_

       end do ifb_
    end do ib_
    call readwriteblocked_write(tmp_matrix,th_2c_resp)
    deallocate(tmp_matrix,STAT = status)
    ASSERT(status==0)

    call readwriteblocked_stopwrite(th_2c_resp)

   end subroutine int_send_2c_resp_save
  !*************************************************************


  !*************************************************************
   subroutine int_send_2c_resp_rewrite()
     !
     ! Main purpose: When scf finished,
     !   int_send_2c_resp_rewrite called from response_module
     !   to rewrite 2c coulomb in the correct way for RESTDD
     !
     !
     !------------ Modules used ---------------------------------
     use unique_atom_module, only: n_unique_atoms
     use symmetry_data_module, only: symmetry_data_n_partners
     use debug, only: show
     use msgtag_module, only: msgtag_send_2c_Colmb
!    use iounitadmin_module,   only: openget_iounit,returnclose_iounit
     use ch_response_module,   only: fit_position
     use comm_module
     use xpack
     use io
     implicit none
     !------------ Declaration of local types ---------------------
     integer(kind=i4_kind) :: counter, i_ir
     integer(kind=i4_kind) :: ua1, ua2
     integer(kind=i4_kind) :: l1, l2
     integer(kind=i4_kind) :: ia, ib
     integer(kind=i4_kind) :: dim_ch
     integer(kind=i4_kind) :: n_if_a, n_contr_a, n_uncontr_a, dim_ch_a
     integer(kind=i4_kind) :: n_if_b, n_contr_b, n_uncontr_b, dim_ch_b
     type(readwriteblocked_tapehandle)          :: th_2c_resp
     real(kind=r8_kind), allocatable            :: tmp_matrix(:),coulomb_2c_matrix(:,:), receive_coulomb(:)
     logical                                    :: check_exist
     character(len=4)                           :: ir_char
     integer(kind=i4_kind) :: status, n_ir, i_proc, n_procs
     integer(kind=i4_kind) :: imeta_a, imeta_b, n_ua, dim_diag_ch
!    integer(kind=i4_kind) :: io_unit
     !------------ Executable code --------------------------------

     n_ir = symmetry_data_n_irreps()
     n_ua = n_unique_atoms

     i_ir_: do i_ir = 1, n_ir

        dim_ch      =  dimension_of_fit_ch(i_ir)

        allocate(coulomb_2c_matrix(dim_ch,dim_ch), STAT=status)
        ASSERT(status==0)
        coulomb_2c_matrix = 0.0D0

        ua1_: do ua1 = 1, N_unique_atoms
           ua2_: do ua2 = 1, N_unique_atoms

              l1_: do l1 = -1, unique_atoms(ua1)%lmax_ch

                 !! It's stupid, but it's so... sorry, not my idea!
                 if (l1 == -1) then
                    n_contr_a   = unique_atoms(ua1)%l_ch(0)%N_contracted_fcts
                    n_uncontr_a = unique_atoms(ua1)%l_ch(0)%N_uncontracted_fcts
                 elseif (l1==0) then
                    n_contr_a   = unique_atoms(ua1)%r2_ch%N_contracted_fcts
                    n_uncontr_a = unique_atoms(ua1)%r2_ch%N_uncontracted_fcts
                 else
                    n_contr_a   = unique_atoms(ua1)%l_ch(l1)%N_contracted_fcts
                    n_uncontr_a = unique_atoms(ua1)%l_ch(l1)%N_uncontracted_fcts
                 end if

                 dim_ch_a = n_contr_a+n_uncontr_a

                 if (dim_ch_a==0) cycle

                 l2_: do l2 = -1, unique_atoms(ua2)%lmax_ch


                    inquire (file = trim(tmpfile(fname(ua1,ua2,l1,l2,i_ir))), exist=check_exist)

!                    if (check_exist == .false.) then
                    if (.not.check_exist) then
                       cycle
                    end if

                    if (l1==-1) then
                       n_if_a = unique_atoms(ua1)%symadapt_partner(i_ir,0)%N_independent_fcts
                    else
                       n_if_a = unique_atoms(ua1)%symadapt_partner(i_ir,l1)%N_independent_fcts
                    end if

                    if (l2==-1) then
                       n_if_b = unique_atoms(ua2)%symadapt_partner(i_ir,0)%N_independent_fcts
                    else
                       n_if_b = unique_atoms(ua2)%symadapt_partner(i_ir,l2)%N_independent_fcts
                    end if

                    call readwriteblocked_startread(trim(tmpfile(fname(ua1,ua2,l1,l2,i_ir))), th_2c_resp)

                    !! The same headache !!
                    if (l2 == -1) then
                       n_contr_b   = unique_atoms(ua2)%l_ch(0)%N_contracted_fcts
                       n_uncontr_b = unique_atoms(ua2)%l_ch(0)%N_uncontracted_fcts
                    elseif (l2 == 0) then
                       n_contr_b   = unique_atoms(ua2)%r2_ch%N_contracted_fcts
                       n_uncontr_b = unique_atoms(ua2)%r2_ch%N_uncontracted_fcts
                    else
                       n_contr_b   = unique_atoms(ua2)%l_ch(l2)%N_contracted_fcts
                       n_uncontr_b = unique_atoms(ua2)%l_ch(l2)%N_uncontracted_fcts
                    end if

                    dim_ch_b    = n_contr_b + n_uncontr_b

                    if (dim_ch_b == 0) cycle

                    allocate(tmp_matrix((dim_ch_a*n_if_a)*(dim_ch_b*n_if_b)))

                    counter   = 0
                    call readwriteblocked_read(tmp_matrix,th_2c_resp)

                    imeta_a = fit_position(i_ir,ua1,l1-1)
                    imeta_b = fit_position(i_ir,ua2,l2-1)

                    expa_: do ia = 1, dim_ch_a*n_if_a
                       expb_: do ib = 1, dim_ch_b*n_if_b
                          counter = counter + 1
                          coulomb_2c_matrix(imeta_a+ia,imeta_b+ib) = tmp_matrix(counter)
                          coulomb_2c_matrix(imeta_b+ib,imeta_a+ia) = tmp_matrix(counter)
                       end do expb_
                    end do expa_

                    deallocate(tmp_matrix)
                    call readwriteblocked_stopread(th_2c_resp)

                    end do l2_
                 end do l1_

              end do ua2_
           end do ua1_

           dim_diag_ch = (dim_ch*(dim_ch+1))/2
           allocate(tmp_matrix(dim_diag_ch),STAT=status)
           ASSERT(status==0)
           counter = 0
           do ia = 1, dim_ch
              do ib = 1, ia
                 counter = counter + 1
                 tmp_matrix(counter) = coulomb_2c_matrix(ia,ib)
              end do
           end do

           deallocate (coulomb_2c_matrix,STAT = status)
           ASSERT(status==0)

           !! begin: send-receive block
           !! all 2center Coulombs will be:
           !! 1) sended from slave
           !! 2) received on the master
           !! 3) added to the main one on the master
           !! NOTE: only master coulomb will be complete, other not!!!
           n_procs = comm_get_n_processors()
           if (comm_i_am_master()) then
              ALLOCATE(receive_coulomb(dim_diag_ch),STAT=status)
              ASSERT(status==0)
              do i_proc=2,n_procs
                 call comm_save_recv (i_proc, msgtag_send_2c_Colmb)
                 call upck(receive_coulomb)
                 tmp_matrix = tmp_matrix + receive_coulomb
              end do
              DEALLOCATE(receive_coulomb,STAT=status)
              ASSERT(status==0)
           else
              call comm_init_send (comm_master_host, msgtag_send_2c_Colmb)
              call pck(tmp_matrix)
              call comm_send()
           end if
           !! end: send-receive block

           write (ir_char, '(i4)') i_ir
           ir_char = adjustl(ir_char)

           if (comm_i_am_master()) then
!             io_unit = openget_iounit(trim(resp_dir)//'/'&
!                  //'co_2c_'//trim(ir_char),form='unformatted',status='unknown')
              call write_buffer( trim(resp_dir)//'/'//'co_2c_'//trim(ir_char) &
                               , tmp_matrix                                   &
                               )
!             call returnclose_iounit(io_unit)
           end if

           deallocate(tmp_matrix,STAT=status)
           ASSERT(status==0)

        end do i_ir_

      end subroutine int_send_2c_resp_rewrite
  !*************************************************************

   character(len=32) function fname(ua1,ua2,l1,l2,i_ir)
    implicit none
    integer(kind = i4_kind), intent(in) :: ua1,ua2,l1,l2
    integer(kind = i4_kind), intent(in) :: i_ir
    !------------ Declaration of local types ---------------------
    character(len=4) :: ua1_char,ua2_char
    character(len=4) :: l1_char,l2_char
    character(len=4) :: ir_char
    !------------ Executable code ------------------------------------

    write (ua1_char,'(i4)') ua1
    write (ua2_char,'(i4)') ua2

    if (l1 == -1) then
       write (l1_char, '(i4)') 0_i4_kind
    elseif (l1 == 0) then
       l1_char = 'r'
    else
       write (l1_char, '(i4)') l1
    end if

    if (l2 == -1) then
       write (l2_char, '(i4)') 0_i4_kind
    elseif (l2 == 0 ) then
       l2_char = 'r'
    else
       write (l2_char, '(i4)') l2
    end if

    write (ir_char, '(i4)') i_ir
    ua1_char = adjustl(ua1_char)
    ua2_char = adjustl(ua2_char)
    l1_char  = adjustl(l1_char)
    l2_char  = adjustl(l2_char)
    ir_char = adjustl(ir_char)
    fname = 'resp_2c_'//trim(ua1_char)//"."//trim(ua2_char)//"."//trim(l1_char)//"."//trim(l2_char)//&
         '_'//trim(ir_char)

  end function fname

  !--------------- End of module -------------------------------------
end module int_send_2c_resp
