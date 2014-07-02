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
MODULE  init_tddft_module
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
#include <def.h>
  USE type_module          ! contains standard data types
  USE iounitadmin_module   ! to open output units
  USE debug
  USE msgtag_module, ONLY: msgtag_tddft_eps_eta
  USE xpack

  IMPLICIT NONE
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ public functions and subroutines ---------------------
  PUBLIC init_tddft_start

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of constants and variables ---------------


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
CONTAINS


  !*************************************************************
  SUBROUTINE init_tddft_start()
    !  Purpose:
    !    1. Read the df data header file (N_as,N_K etc.)
    !    2. Read epsilon, eta from input tape
    !------------ Modules used ----------------------------------
    USE comm_module
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    if(comm_i_am_master()) then
       if(comm_parallel()) then
          call comm_init_send(comm_all_other_hosts,msgtag_tddft_eps_eta)
          call comm_send()
       end if
    end if

    if (comm_i_am_master()) &
         CALL write_to_output_units(&
         & '  init_start: call df_data_read_header()')
    CALL init_send_eigenparam()

    if (comm_i_am_master()) then
       CALL write_to_output_units('  init_start: call df_data_read_header()')
    end if

    CALL init_read_header()

    if (comm_i_am_master()) then
       CALL write_to_output_units('  init_start: call df_data_read_eps_eta()')
    end if

    CALL init_read_eps_eta()

  END SUBROUTINE init_tddft_start
  !*************************************************************

  !*************************************************************
  SUBROUTINE init_send_eigenparam()
    !  Purpose:
    !  send to the slave eigensolveparameters.
    !  !!do not forget send it all!!
    !------------ Modules used ----------------------------------
    USE comm_module
!!    USE comm_bcast, only: bcast
    USE global_module, ONLY: gl_eig_crite
    USE msgtag_module, ONLY: msgtag_tddft_sendeig
    USE xpack

    INTEGER(KIND = i4_kind) :: i_proc
!!    INTEGER(KIND = i4_kind) :: num

!!#if 0
    if (comm_i_am_master()) then
       do i_proc = 2, comm_get_n_processors()
          call comm_init_send(i_proc,msgtag_tddft_sendeig)
          call pck(gl_eig_crite)
          call comm_send()
       end do
    else
       call comm_save_recv(comm_master_host,msgtag_tddft_sendeig)
       call upck(gl_eig_crite)
    end if
!!#endif

!!    num = gl_eig_crite
!!    call bcast(num)
!!    gl_eig_crite = num

  END SUBROUTINE init_send_eigenparam



  !*************************************************************
  SUBROUTINE init_read_header()
    !  Purpose:
    !  Read df data header file
    !------------ Modules used ----------------------------------
    USE filename_module, ONLY: data_dir
    USE comm_module
    USE global_module
    USE exchange
    IMPLICIT NONE
    !------------ Declaration of local variables ---------------------

    INTEGER(KIND=i4_kind) :: io_unit, i_ir_c, in_chfit
    INTEGER(KIND=i4_kind) :: io_stat
    INTEGER(KIND=i4_kind) :: i_spin, i_skip, in_spin,&
         & in_irrep_c, in_trans, n_procs

    ! 25 +  230 ==  255. FIXME:  the code here  does not  use namelist
    ! input, instead it tries to  parse the inputs itself (badly). The
    ! formats here  must skip the "NAME  =" part of  the namelist line
    ! and have a wide enough field for the actual value:
    character (len=255) :: input_line
    character (len=*), parameter :: ifmt = "(25X, I230)"
    character (len=*), parameter :: dfmt = "(25X, F230.10)"
    character (len=*), parameter :: lfmt = "(25X, L230)"

    ! This is for fancy output only:
    character (len=*), parameter :: prefix = "      * "

    REAL (r8_kind) :: LeV

    integer(i4_kind), parameter :: & !! VERY IMPORTANT
       & X_NONE  = 0, &
       & C_NONE  = 0, &
       & C_VWN   = 1, &
       & C_PWlda = 2, &
       & C_PBE   = 3, &
       & C_PW91  = 4, &
       & C_PERDEW = 5

    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    n_procs = comm_get_n_processors()

    if (comm_i_am_master()) then
       CALL write_to_output_units( "    df_data_read_header: read header of interface file...")
    end if

    ! first open the interfacefile from the DF program program
    io_unit = openget_iounit(file=TRIM(data_dir)//"/resp_header.out",&
         status="old",form="formatted",action="read")

    if (comm_i_am_master()) then
       CALL write_to_output_units("    TDDFT: HEADER  ")
    end if

    in_chfit = 0

    do i_skip = 1, 500

       READ (io_unit, FMT='(A)', IOSTAT=io_stat) input_line
       ASSERT(io_stat==0)

       call clean_default(input_line)

       if (0 /= index (input_line, 'NUM_SPINS')) then
          READ (input_line, fmt=ifmt, IOSTAT=io_stat) gl_N_spin
          ASSERT(io_stat==0)

          if (comm_i_am_master()) then
             if (gl_N_spin == 1) then
                CALL write_to_output_units(prefix // "    proceeding closed shell...")
             else
                CALL write_to_output_units(prefix // "    proceeding open   shell...")
             end if
          end if

       elseif (0 /= index (input_line, 'NUM_IRREPS')) then
          READ (input_line, ifmt, IOSTAT=io_stat) gl_N_irr
          ASSERT(io_stat==0)
          if (comm_i_am_master()) &
               CALL write_to_output_units(prefix // "    IRREPS             = ", gl_N_irr)
          call global_alloc('chr',gl_N_irr)
          call global_alloc('nas',gl_N_irr)
          call global_alloc('mlt',gl_N_irr)
          call global_alloc('trs',gl_N_irr)
       elseif (0 /= index (input_line, 'NUM_CHARGEFITFCTS')) then
          in_chfit = in_chfit+1
          READ (input_line, ifmt, IOSTAT=io_stat) gl_N_k(in_chfit)
          ASSERT(io_stat==0)

          ! RESPONSE CONTROL BLOCK
       elseif (0 /= index (input_line, "TARGET")) then
          gl_SS = .false.
          gl_ST = .false.
          if (0 /= index (input_line, "SS")) then
             gl_SS = .true.
          end if
          if (0 /= index (input_line, "ST")) then
             gl_ST = .true.
          end if
          if (comm_i_am_master()) then
             if (gl_N_spin == 1 ) then
                if (gl_SS) then
                   CALL write_to_output_units(prefix // "    singlet -> singlet")
                end if
                if (gl_ST) then
                   CALL write_to_output_units(prefix // "    singlet -> triplet")
                end if
             end if
          end if

       elseif (0/=index(input_line,'EXCHANGE')) then
          if     (0/=index(input_line,"NONE")) then
             gl_X  = X_NONE
          elseif (0/=index(input_line,"XA"  )) then
             gl_X  = X_XALPHA
          elseif (0/=index(input_line,"BECKE"  )) then
             gl_X  = X_BECKE88
          elseif (0/=index(input_line,"PBE" )) then
             gl_X  = X_PBE
          elseif (0/=index(input_line,"REVPBE" )) then
             gl_X  = X_REVPBE
          elseif (0/=index(input_line,"PBEN" )) then
             gl_X  = X_PBEN
          elseif (0/=index(input_line,"PW91")) then
             gl_X  = X_PW91
          end if
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0/=index(input_line,'CORRELATION')) then
          if     (0/=index(input_line,"NONE" )) then
             gl_C  = C_NONE
          elseif (0/=index(input_line,"VWN"  )) then
             gl_C  = C_VWN
          elseif (0/=index(input_line,"PERDEW")) then
             gl_C  = C_PERDEW
          elseif (0/=index(input_line,"PBE"  )) then
             gl_C  = C_PBE
          elseif (0/=index(input_line,"PW91" )) then
             gl_C  = C_PW91
          elseif (0/=index(input_line,"PWlda")) then
             gl_C  = C_PWlda
          end if
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if


       elseif (0 /= index (input_line, "CALC_ALL")) then
          READ (input_line, lfmt, IOSTAT=io_stat) gl_calcall
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0 /= index (input_line, "LANCZOS")) then
          READ (input_line, lfmt, IOSTAT=io_stat) gl_lanczos
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0 /= index (input_line, "S_APP")) then
          READ (input_line, lfmt, IOSTAT=io_stat) gl_S_App
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0 /= index (input_line, "noRI")) then
          READ (input_line, lfmt, IOSTAT=io_stat) gl_noRI
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0 /= index (input_line, "CALC_N_LOW")) then
          READ (input_line, ifmt,IOSTAT=io_stat) gl_Nlow
          ASSERT(io_stat==0)
          if (.not. gl_calcall) then
             if (comm_i_am_master()) then
                CALL write_to_output_units(trim(prefix // input_line))
             end if
          else
             gl_Nlow = 65535
          end if

       elseif (0 /= index (input_line, "LOWESTeV")) then
          READ (input_line, dfmt, IOSTAT=io_stat) LeV
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if
          gl_LOWESTeV = (LeV * 0.0367493090_r8_kind) ** 2

       elseif (0 /= index (input_line, "MAX_SP_TRANS")) then
          READ (input_line, ifmt, IOSTAT=io_stat) gl_max_sp_trans
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0 /= index (input_line, "MAX_ITER")) then
          READ (input_line, ifmt, IOSTAT=io_stat) gl_MaxIt
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0 /= index (input_line, "CALC_OSC_STR")) then
          READ (input_line, lfmt, IOSTAT=io_stat) gl_oscstr
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

       elseif (0 /= index (input_line, "NTO")) then ! MH: NTO Calculation
          READ (input_line, lfmt, IOSTAT=io_stat) gl_NTO
          ASSERT(io_stat==0)
          if (comm_i_am_master()) then
             CALL write_to_output_units(trim(prefix // input_line))
          end if

          ! END OF RESPONSE CONTROL BLOCK

       elseif (0 /= index (input_line, 'NUMBER OF TRANSITIONS')) then
          do i_spin = 1, gl_N_spin
             do i_ir_c = 1, gl_N_irr
!               if (eof(io_unit)) EXIT
                READ (io_unit, '(4X,I2,5X,I2,5X,I5)', IOSTAT=io_stat) in_irrep_c, in_spin, in_trans
                if (io_stat == -1) EXIT ! the loop at EOF!
                ASSERT(io_stat==0)

                gl_N_as(in_irrep_c, in_spin) = in_trans

                ! FIXME: What is it?
                if (comm_parallel()) then
                   gl_N_as_slave(in_irrep_c, in_spin) = INT(gl_N_as(in_irrep_c,in_spin)/n_procs)
                   gl_N_as_mastr(in_irrep_c, in_spin) = gl_N_as(in_irrep_c, in_spin) &
                        &                             - (n_procs-1) * gl_N_as_slave(in_irrep_c, in_spin)
                else
                   gl_N_as_slave(in_irrep_c, in_spin) = 0_i4_kind
                   gl_N_as_mastr(in_irrep_c, in_spin) = gl_N_as(in_irrep_c, in_spin)
                end if

             end do
          end do
          EXIT                  ! ... without checking io_stat
       end if
       ASSERT(io_stat==0)
    end do

    if ((gl_X == X_NONE) .and. (gl_C  == C_NONE)) then
       gl_XC = .false.
    else
       gl_XC = .true.
    end if

    CALL returnclose_iounit(io_unit)

  contains

    subroutine clean_default(input_line)
      implicit none
      CHARACTER(LEN=255),intent(inout)    :: input_line
      !------------ Declaration of local variables -----------------
      INTEGER(KIND=i4_kind) :: j
      !------------ Executable code --------------------------------
      j = index(input_line,"# (the default)")
      if (j.ne.0) then
         input_line = input_line(:j-1)
      end if

    end subroutine clean_default

  END SUBROUTINE init_read_header

  !*************************************************************
  SUBROUTINE init_read_eps_eta()
    !  Purpose:
    !  Read the differences of MO energy eigenvalues
    !       (e_a - e_s),    a=occupied, s=unoccupied
    !  and the differences of MO occupation numbers
    !       (n_a - n_s)
    !  and the MO level indices of the occupied and unoccupied
    !       MO's in the ground state spectrum
    !  from the binary data file "resp_eig_occ_ind.dat"
    !  Note:
    !  1) The order is the same as for the usual metaindex 'as'
    !  2) In the notation of the paper we need the MO energy
    !     eigenvalue differences just the other way round:
    !       (e_s - e_a)
    !     Therefore multiply the MO energy differences from tape
    !     by -1.
    !------------ Modules used ----------------------------------
    USE filename_module, ONLY: data_dir
    USE comm_module

    USE global_module, ONLY: &
         global_alloc, &
         gl_eps, &
         gl_eta, &
         gl_MO, &
         gl_N_irr, &
         gl_N_spin, &
         gl_N_as, &
         gl_N_as_slave, &
         gl_N_as_mastr, &
         gl_what_N_as, &
         gl_IRR

    USE phys_param_module, ONLY: hartree2eV

    IMPLICIT NONE
    !------------ Declaration of local variables ---------------------
    INTEGER(KIND=i4_kind) :: io_unit,io_stat, i_ir_c, i_spin
    INTEGER(KIND=i4_kind) :: as_c, dummy, ir_a, ir_b
    INTEGER(KIND=i4_kind) :: n_spin, n_irrep

    INTEGER(KIND=i4_kind) :: a_index, s_index

    REAL(KIND=r8_kind)    :: eps_read, eta_read

    INTEGER(KIND=i4_kind),POINTER :: nas(:)

!!    INTEGER(KIND=i4_kind) :: n_procs

    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    n_spin  = gl_N_spin
    n_irrep = gl_N_irr

    CALL global_alloc('prc',n_irrep,gl_N_as)

    if (comm_i_am_master()) then
       CALL global_alloc('ALL',n_irrep,gl_N_as_mastr)
    else
       CALL global_alloc('ALL',n_irrep,gl_N_as_slave)
    end if

    i_spin_: dO i_spin = 1, n_spin
       i_ir_c_: do i_ir_c = 1, n_irrep
          io_unit = openget_iounit(file=TRIM(data_dir)//'/'//fname('eg_oc',i_ir_c,i_spin),&
               form="unformatted",action="READ")

          as_c_: DO as_c = 1, gl_N_as(i_ir_c,i_spin)
             READ(io_unit, IOSTAT=io_stat) &
                  & dummy,   &
                  & ir_a,    &
                  & ir_b,    &
                  & a_index, &
                  & s_index, &
                  & eps_read,&
                  & eta_read
             ASSERT(io_stat==0)
             eps_read = -eps_read

             !! SPLIITING TO PROCESSORS
             nas => gl_what_N_as(i_ir_c,i_spin)%m
             if ((nas(as_c)) == -1) cycle
             gl_MO (i_ir_c,i_spin)%m(nas(as_c),1) = a_index
             gl_MO (i_ir_c,i_spin)%m(nas(as_c),2) = s_index
             gl_eps(i_ir_c,i_spin)%m(nas(as_c))   = eps_read
             gl_eta(i_ir_c,i_spin)%m(nas(as_c))   = eta_read
             gl_IRR(i_ir_c,i_spin)%m(nas(as_c),1) = ir_a
             gl_IRR(i_ir_c,i_spin)%m(nas(as_c),2) = ir_b

          end do as_c_

       end do i_ir_c_
    end do i_spin_

    ! close input tape
    CALL returnclose_iounit(io_unit)

#if 0
    !---begin debug
    PRINT*,""
    PRINT*,"**********************************"
    PRINT*,"*** Information about Spectrum ***"
    do i_ir_c = 1, n_irrep
       DO i_spin = 1,n_spin
          print *,"Spin = ",i_spin
          n_procs = comm_myindex()
          print *,"N_proc = ", n_procs
          PRINT*,"**********************************"
          DO as_c = 1, size(gl_eps(i_ir_c,i_spin)%m,1)
             WRITE(*,fmt='(3I5,2X,2(E30.15,2X))') as_c, &
                  & gl_MO (i_ir_c,i_spin)%m(as_c,1), &
                  & gl_MO (i_ir_c,i_spin)%m(as_c,2), &
                  & gl_eps(i_ir_c,i_spin)%m(as_c), &
                  & gl_eta(i_ir_c,i_spin)%m(as_c)

          ENDDO
       ENDDO
    end do
    PRINT*,"********"
    !---end debug
#endif

  END SUBROUTINE init_read_eps_eta
  !*************************************************************

  !*************************************************************
  character(len=32) function fname(name,i_ir,i_sp)
    implicit none
    integer  (kind = i4_kind), intent(in) :: i_ir,i_sp
    character(len  = 5),       intent(in) :: name
    !------------ Declaration of local types ---------------------
    character(len=4) :: irc_char,isp_char
    character(len=5) :: fnm_char
    !------------ Executable code ------------------------------------

    write (irc_char, '(i4)') i_ir
    write (isp_char, '(i1)') i_sp
    irc_char = adjustl(irc_char)
    isp_char = adjustl(isp_char)
    fnm_char = adjustl(name)
    fname    = trim(fnm_char)//"_"//trim(irc_char)//"_"//trim(isp_char)
  end function fname
  !*************************************************************



!!$  !*************************************************************
!!$  subroutine df_data_
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
!!$  end subroutine df_data_
!!$  !*************************************************************


  !--------------- End of module -------------------------------------
END MODULE init_tddft_module
