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
#include <def.h>
!===============================================================
! Public interface of module
!===============================================================
MODULE  read_module
  !---------------------------------------------------------------
  !
  !  Purpose:
  !
  !  Detailization:
  !  --------------
  !
  !  Module called by: closed_shell_module, open_shell_module
  !
  !
  !  References: ...
  !  Remarks   : As usual NOTHING is stored within this module.
  !              Any important information of needed by other parts
  !              parts of the program will be stored in the
  !              "global_data_module.f90"
  !
  !  Author: SB
  !  Date:   04/05
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  USE debug, only: show
  USE type_module          ! contains standard data types
! USE iounitadmin_module   ! to open output units
  USE output_module        ! defines amount of output
  USE datatype, only: ARRMAT2
  USE filename_module, ONLY: data_dir
! USE io

  IMPLICIT NONE
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  PUBLIC read_C, read_R, read_Q, read_noRI_XC

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
CONTAINS

  !*************************************************************
  SUBROUTINE read_C(i_ir,i_sp,eta,eps,N,NK,NAS,RM)
    !  Purpose:
    !  Distribute matrix C to all processors
    ! called by: calc_diag_send_L()
    !------------ Modules used ----------------------------------
    USE filename_module,   ONLY: data_dir
!   USE iounitadmin_module,ONLY: openget_iounit,returnclose_iounit
    USE global_module,     ONLY: gl_N_spin, gl_N_as
    USE comm_module
    USE debug
    USE io, only: read_buffer

    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    INTEGER(KIND=i4_kind), INTENT(in)    :: i_ir, i_sp, N, NK
    REAL   (KIND=r8_kind), INTENT(in)    :: eta(:),eps(:)
    INTEGER(KIND=i4_kind), INTENT(in)    :: NAS(:)
    REAL(KIND=r8_kind),    INTENT(INOUT) :: RM(:,:)
    !*************************************************************
    INTEGER(KIND=i4_kind)             :: status!!, k_loop
    REAL(KIND=r8_kind), ALLOCATABLE   :: TV(:,:), SV(:)
    INTEGER(KIND=i4_kind)             :: ias, CC !, io_unit
    LOGICAL                           :: file_exist

    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code --------------------------------

!!#if 0
    inquire(file = trim(data_dir)//"/"//trim(fname('co_3c',i_ir,i_sp)), exist=file_exist )
    if (.not.file_exist) return
!!#endif

    !** get memory for temp. column vector of matrix L

    ASSERT(N==gl_N_as(i_ir,i_sp))
    RM = 0.0D0

    ALLOCATE(TV(N,NK),&
         &   SV(size(eta)),&
         &   STAT=status)
    ASSERT(status==0)

    if (gl_N_spin == 1) then
       CC = 2
    else
       CC = 1
    end if
    SV = SQRT(REAL(CC * eps * eta,r8_kind))

!   io_unit = openget_iounit(trim(data_dir)//'/'&
!        //fname('co_3c',i_ir,i_sp),form='unformatted',status='unknown')
    CALL read_buffer( trim(data_dir)//'/'//fname('co_3c',i_ir,i_sp) &
                    , TV &
                    )
!   call returnclose_iounit(io_unit)

    do ias = 1, N
       if(NAS(ias) == -1) cycle
       RM(NAS(ias),:) = RM(NAS(ias),:) &
            &         + SV(NAS(ias)  ) * TV(    ias ,:)
    end do

    DEALLOCATE(TV,SV,stat=status)
    ASSERT(status==0)

  END SUBROUTINE read_C
  !*************************************************************

  !*************************************************************
  SUBROUTINE read_R(i_ir, XC)
    !  Purpose:
    !  Read 2-index-integrals <g_l|fxc_aa|g_k>, <g_l|fxc_ab|g_k>
    !  Read 2-index-integrals <g_l|fxc_ba|g_k>, <g_l|fxc_bb|g_k>
    !  from file.
    !  the output matrix is constructed as
    !
    !  open_shell case:
    !  tmp_ptr(:,1) : temp_vec = <g_l|fxc_aa |g_k>
    !  tmp_ptr(:,2) : temp_vec = <g_l|fxc_ab |g_k>
    !  tmp_ptr(:,3) : temp_vec = <g_l|fxc_ba |g_k>
    !  tmp_ptr(:,4) : temp_vec = <g_l|fxc_bb |g_k>
    !
    !  closed_shell case:
    !  tmp_ptr(:,1) : temp_vec = <g_l|fxc_aa |g_k>
    !  tmp_ptr(:,2) : temp_vec = <g_l|fxc_ab |g_k>
    !
    ! ** NOTE ** tmp_ptr(:,:) is assumed to be symmetrical
    !
    ! called by: open_shell_
    !------------ Modules used ----------------------------------
    USE global_module, ONLY: gl_N_k,gl_N_spin
    use io, only: read_buffer
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND = i4_kind), INTENT(in   ) :: i_ir
    TYPE(ARRMAT2),           INTENT(inout) :: XC(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    INTEGER(KIND=i4_kind)             :: alloc_stat,i_fit_k,i_fit_l,&
         & nmat,counter, n_k, i
    INTEGER(KIND=i4_kind)             :: n_spin
    REAL(KIND=r8_kind), ALLOCATABLE   :: help_matrix(:)
    LOGICAL                           :: file_exist
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code --------------------------------

    ! allocate memory for auxiliary array
    n_k = gl_N_k(i_ir)
    nmat = (n_k*(n_k+1))/2

    ALLOCATE( help_matrix(nmat),STAT=alloc_stat)
    IF(alloc_stat.ne.0) CALL error_handler&
         ("    df_data_read_R: allocation of help_matrix failed")

    ! read the matrix f_xc from input tape by first reading the help_matrix
    ! from tape and then redistributing the data on f_xc
    ! Note: we assume C1 symmetry => n_irrep =1
    help_matrix = 0.0_r8_kind

    ! read matrix fxc_aa, fxc_ab, fxc_ba, fxc_bb
    ! read the auxiliary matrix = lower half of the true matrix
    n_spin = gl_N_spin

    DO i = 1,n_spin * 2
!      io_unit = openget_iounit(trim(data_dir)//'/'//fname('xc_2c',i_ir,i),&
!           form='unformatted',status='unknown')

       inquire(file = trim(data_dir)//"/"//trim(fname('xc_2c',i_ir,i)), exist=file_exist )
       if (.not.file_exist) cycle

       CALL read_buffer( trim(data_dir)//'/'//fname('xc_2c',i_ir,i) &
                       , help_matrix &
                       )
       ! and recreate true matrix from aux. matrix
       counter=1
       DO i_fit_k=1,n_k
          DO i_fit_l=1,i_fit_k
             XC(i)%m(i_fit_k,i_fit_l) = help_matrix(counter)
             IF(i_fit_l/=i_fit_k) &
                  &  XC(i)%m(i_fit_l,i_fit_k) = help_matrix(counter)
             counter=counter+1
          END DO
       END DO
!      call returnclose_iounit(io_unit)
    END DO

    ! deallocate auxiliary array
    DEALLOCATE(help_matrix,stat=alloc_stat)
    IF(alloc_stat/=0) CALL error_handler(&
         'df_data_read_R: deallocating help_matrix failed')

  END SUBROUTINE read_R
  !*************************************************************

  !*************************************************************
  SUBROUTINE read_Q(i_ir, tmp_ptr)
    !  Purpose:
    !  Read 2-index-integrals <g_l|1/(r-r')|g_k> from file.
    !
    ! ** NOTE ** tmp_ptr(:,:) is assumed to be symmetrical
    !
    ! called by: open_shell_
    !------------ Modules used ----------------------------------
    use io, only: read_buffer
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind),intent(in)  :: i_ir
    REAL(KIND=r8_kind),intent(inout)  :: tmp_ptr(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    INTEGER(KIND=i4_kind)             :: alloc_stat, ia, ib, counter, dim
    INTEGER(KIND=i4_kind)             :: n_k
    character(len=4)                  :: ir_char
    REAL(KIND=r8_kind), ALLOCATABLE   :: coulomb_2c_matrix(:)
    LOGICAL                           :: file_exist
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code --------------------------------

    n_k = size(tmp_ptr,1)
    write (ir_char, '(i4)') i_ir
    ir_char = adjustl(ir_char)
    tmp_ptr = 0.0

    inquire(file = trim(data_dir)//"/co_2c_"//trim(ir_char), exist=file_exist )
    if (.not.file_exist) return

    dim = INT((n_k*(n_k+1))/2)
!!    print *,"dim = ",dim
    ALLOCATE(coulomb_2c_matrix(dim),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
!!    print *,"PATH = ",trim(data_dir)//'/co_2c_'//trim(ir_char)
!   io_unit = openget_iounit(trim(data_dir)//'/co_2c_'//trim(ir_char),&
!        form='unformatted',status='unknown')

    CALL read_buffer( trim(data_dir)//'/co_2c_'//trim(ir_char) &
                    , coulomb_2c_matrix &
                    )
    counter = 0
    ia_: do ia = 1, n_k
       ib_: do ib = 1, ia
          counter = counter + 1
          tmp_ptr(ia,ib) = coulomb_2c_matrix(counter)
          tmp_ptr(ib,ia) = coulomb_2c_matrix(counter)
       end do ib_
    end do ia_

!   call returnclose_iounit(io_unit)

    ! deallocate auxiliary array
    DEALLOCATE(coulomb_2c_matrix,stat=alloc_stat)
    ASSERT(alloc_stat==0)

  END SUBROUTINE read_Q
  !*************************************************************

  !*************************************************************
  SUBROUTINE read_noRI_XC(i_ir_c,eta_up,eps_up,XC,&
       &                         eta_dn,eps_dn)
    !  Purpose:
    !------------ Modules used ----------------------------------
    USE global_module,      ONLY: gl_N_spin
    USE io,                 ONLY: read_buffer
!   USE iounitadmin_module, ONLY: openget_iounit,returnclose_iounit
    USE filename_module,    ONLY: data_dir
    USE linalg_module,      ONLY: matvecmul
    USE datatype,           ONLY: arrmat1, arrmat2
    USE constants,          ONLY: ZERO,ONE,TWO

    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    INTEGER(KIND=i4_kind), INTENT(IN )          :: i_ir_c
    REAL   (KIND=r8_kind), INTENT(IN )          :: eta_up(:),eps_up(:)
    TYPE   (     ARRMAT2), INTENT(OUT)          :: XC(:)
    REAL   (KIND=r8_kind), INTENT(IN ),OPTIONAL :: eta_dn(:),eps_dn(:)
    !*************************************************************
    INTEGER(KIND=i4_kind)             :: status, NSP, CC, idx, i, j
    TYPE(ARRMAT1), ALLOCATABLE, target :: epsxeta(:)
    INTEGER(KIND=i4_kind)             :: CF(2)
    REAL(kind = r8_kind), POINTER     :: R(:), L(:)
    REAL(kind = r8_kind), ALLOCATABLE :: XCR(:)
    LOGICAL                           :: file_exist

    INTEGER(KIND = i4_kind), PARAMETER :: &
         &    UP = 1, &
         &    DN = 2

    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code --------------------------------

!!    print *,"i_ir_c = ",i_ir_c

    NSP = gl_N_spin

    ALLOCATE(epsxeta(NSP),STAT=status)
    ASSERT(status==0)

    ALLOCATE(epsxeta(UP)%m(size(eta_up)),STAT=status)
    ASSERT(status==0)
    CC = 1

    if (NSP==2) then
       ALLOCATE(epsxeta(DN)%m(size(eta_dn)),STAT=status)
       ASSERT(status==0)
       CC = 2
    end if

    epsxeta(UP)%m             = CC * eps_up * eta_up
    if (nsp==2) epsxeta(DN)%m = CC * eps_dn * eta_dn

    do idx = 1, size(XC,1)

       inquire(file = trim(data_dir)//"/"&
            //trim(co_4c_fname(i_ir_c,idx)), exist=file_exist )
       if (.not.file_exist) return

       call read_buffer( trim(data_dir)//'/'//trim(co_4c_fname(i_ir_c,idx)) &
                       , XC(idx)%m )

       select case (NSP)
       case (1) !! closed shell
          CF = (/UP,UP/)
       case (2) !! open   shell
          select case (idx)
          case (1)
             CF = (/UP,UP/)
          case (2)
             CF = (/UP,DN/)
          case (3)
             CF = (/DN,DN/)
          case (4)
             CF = (/DN,UP/)
          end select
       end select

       R => epsxeta(CF(1))%m
       L => epsxeta(CF(2))%m
       ALLOCATE(XCR(size(R)),STAT = status)
       ASSERT(status==0)
       call matvecmul(XC(idx)%m, R, XCR)
       !!FIXME: BLAS HERE!!!
       do i = 1, size(L)
          do j = 1, size(R)
             XC(idx)%m(i,j) = L(i) * XCR(j)
          end do
       end do
       DEALLOCATE(XCR,STAT=status)
       ASSERT(status==0)

    end do

!   call returnclose_iounit(unit)

  contains

    character(len=32) function co_4c_fname(i_ir_c,idx)
      implicit none
      integer(kind = i4_kind), intent(in) :: i_ir_c,idx
      !------------ Declaration of local types ---------------------
      character(len=4) :: irc_char
      character(len=4) :: idx_char
      !------------ Executable code --------------------------------

      write (irc_char, '(i4)') i_ir_c
      write (idx_char, '(i1)') idx
      irc_char = adjustl(irc_char)
      idx_char = adjustl(idx_char)
      co_4c_fname = 'XC_4i_'//trim(irc_char)//'_'//trim(idx_char)

    end function co_4c_fname

  END SUBROUTINE read_noRI_XC
  !*************************************************************

  !*************************************************************
  character(len=32) function fname(name,i_ir,i_sp)
    implicit none
    integer  (kind = i4_kind), intent(in) :: i_ir,i_sp
    character(len  = 5),       intent(in) :: name
    !------------ Declaration of local types ---------------------
    character(len=4) :: irc_char,isp_char
    character(len=5) :: fnm_char
    !------------ Executable code --------------------------------

    write (irc_char, '(i4)') i_ir
    write (isp_char, '(i1)') i_sp
    irc_char = adjustl(irc_char)
    isp_char = adjustl(isp_char)
    fnm_char = adjustl(name)
    fname    = trim(fnm_char)//"_"//trim(irc_char)//"_"//trim(isp_char)
  end function fname
  !*************************************************************


!!$  !*************************************************************
!!$  subroutine read_
!!$  !  Purpose: ..
!!$  !------------ Modules used ----------------------------------
!!$    use
!!$    implicit none
!!$    !------------ Declaration of formal parameters ---------------
!!$    integer(kind=i4_kind), intent(     ) ::
!!$    real(kind=r8_kind),    intent(     ) ::
!!$    logical,               intent(     ) ::
!!$    character,             intent(     ) ::
!!$    !** End of interface *****************************************
!!$    !------------ Declaration of local variables -----------------
!!$    integer(kind=i4_kind)                ::
!!$    real(kind=r8_kind)                   ::
!!$    logical                              ::
!!$    character                            ::
!!$    !------------ Executable code --------------------------------
!!$
!!$
!!$  end subroutine read_
!!$  !*************************************************************


  !--------------- End of module ----------------------------------
END MODULE read_module
