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
module hamiltonian_module
  !-------------------------------------------------------------------
  !-------------- Module specification ---------------------------
  !-------------------------------------------------------------------
  !
  !  Purpose: Contains the 'three-center' parts of the
  !           hamiltonian : HAM_XC,EN_XC
  !           and the total hamiltonian HAM_TOT
  !           provide the kinetic and nuclear parts
  !           of the hamiltonian and a an appropriate
  !           read routine to get it from the
  !           intermediate file written by 'convert'
  !
  !  Contents:
  !   Routines:
  !      - public variables : ham_xc, ham_tot (arrmat3)
  !
  !     (i) PUBLIC     reset_ham   -> allocates and initializes
  !                                   hams
  !                    free_ham    -> free one or more parts of
  !                                   ham
  !
  !            print_hamiltonian   -> print out the hamiltonian
  !
  !     (ii) PRIVATE   select_free -> does the deallocation
  !
  !
  !  Module called by: build_hamiltonian,ham_calc
  !                    prescf_module,scf
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification
  ! Author: TG
  ! Date:   10/96
  ! Description: removed unnecessary spin dependence of
  !              ham_kin and ham_nuc
  !
  ! Modification
  ! Author: FN,TG
  ! Date:   1/97
  ! Description: included ham_kin_nuc_module
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   14/7/97
  ! Description: Control variable OPTIONS_XCMODE() introduced
  !              Control variable OPTIONS_ETOTMODE() introduced
  !
  !
  ! Modification
  ! Author: MM
  ! Date:   10/97
  ! Description: extension to spin orbit
  !              hamiltonian_store and hamiltonian_recover are not yet
  !              adapted to spin orbit
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  !------------ Modules used -----------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype, only: arrmat2, arrmat3 ! user defined types
  use init_module, only: init ! initialization stuff
  use symmetry_data_module, only : sym, &
       symmetry_data_n_irreps,symmetry_data_n_proj_irreps, &
       symmetry_data_dimension, symmetry_data_dimension_proj, &
       symmetry_data_n_spin
  use options_module, only: options_spin_orbit, lvshift_mixing
  implicit none
  private
  save
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of public constants and variables -----
  public :: arrmat3,arrmat2
  type(arrmat3),allocatable,public :: ham_xc(:),ham_tot(:)
  type(arrmat3),allocatable,target,public :: ham_lsft(:) ! level shift hamiltonian

  ! complex hamiltonians for spin orbit
  type(arrmat2),allocatable,public :: ham_tot_real(:)
  type(arrmat2),allocatable,public :: ham_tot_imag(:)
  type(arrmat2),allocatable,public :: ham_xc_real(:)
  type(arrmat2),allocatable,public :: ham_xc_imag(:)

  integer(kind=i4_kind),public :: lsft_allocstat(3)=0,diis_allocstat(3)=0


  !------------ public functions and subroutines ---------------------
  public :: hamiltonian_setup!(ssym), allocate hamtiltonian
  public :: hamiltonian_shutdown!(), clean up the module
  public :: reset_ham, print_hamiltonian, &
            hamiltonian_store,hamiltonian_recover
  ! function read_ham_kin_nuc had to be moved to ham_calc_module


  !===================================================================
  ! End of public interface of module
  !===================================================================

  public :: sym

  !------------ Subroutines ------------------------------------------
contains

  subroutine hamiltonian_setup(ssym)
    !
    ! Purpose: allocate hamiltonian. Called from a parallel context in
    ! prescf_init()
    !
    implicit none
    type(sym), intent(in) :: ssym
    ! *** end of interface ***

    !
    ! Dont forget to use appropriate dimensions of irreps,
    ! (that is use projective irreps in case of spin orbit)
    !
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       if(.not.allocated(ham_tot_real))then
         call alloc_ham(ssym%dim_proj, 0, 'ham_tot_real')
         call alloc_ham(ssym%dim_proj, 0, 'ham_tot_imag')
       endif
    else
       !
       ! STANDARD SCF:
       !
       ! call free_ham('tot') has been moved out of scf-loop:
       if(.not.allocated(ham_tot))then
          call alloc_ham(ssym%dim, ssym%n_spin, 'ham_tot')
       endif

       if(lvshift_mixing.and. .not.allocated(ham_lsft)) then
          call alloc_ham(ssym%dim, ssym%n_spin, 'ham_lsft' )
          lsft_allocstat(2)=1
       endif
    endif ! options_spin_orbit
  end subroutine hamiltonian_setup

  subroutine reset_ham ()
    !
    ! Purpose:  Reset hamiltonian  to zeros.   Called from  a parallel
    ! context in ham_calc_main every SCF cycle.
    !
    implicit none
    ! *** end of interface ***

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       call init (ham_tot_real)
       call init (ham_tot_imag)
    else
       !
       ! STANDARD SCF:
       !
       call init (ham_tot)

       ! FIXME: why NOT allocated?
       if (lvshift_mixing .and. .not. allocated (ham_lsft)) then
          call init (ham_lsft)
       endif
    endif ! options_spin_orbit
  end subroutine reset_ham

  subroutine alloc_ham( dims, n_spin, name )
    !
    ! Allocate one of global module vars that hold hamiltonian (parts)
    ! as specified by "name" and dimensions "dims(1:n_irrep)" and "n_spin"
    !
    ! Called by "reset_ham()"
    !
    implicit none
    integer(i4_kind), intent(in) :: dims(:) ! (n_irrep)
    integer(i4_kind), intent(in) :: n_spin
    character(len=*), intent(in) :: name
    ! *** end of interface ***

    integer(i4_kind) :: status
    integer(i4_kind) :: n_irrep, i

    n_irrep = size(dims)

    select case ( name )

    case ( 'ham_tot' )

      allocate(ham_tot(n_irrep),STAT=status)
      ASSERT(status.eq.0)

      do i=1,n_irrep
         allocate(ham_tot(i)%m(dims(i), dims(i), n_spin), STAT=status)
         ASSERT(status.eq.0)
      enddo

    case ( 'ham_tot_real' )

      allocate(ham_tot_real(n_irrep),STAT=status)
      ASSERT(status.eq.0)

      do i=1,n_irrep
         allocate(ham_tot_real(i)%m(dims(i), dims(i)), STAT=status)
         ASSERT(status.eq.0)
      enddo

    case ( 'ham_tot_imag' )

      allocate(ham_tot_imag(n_irrep),STAT=status)
      ASSERT(status.eq.0)

      do i=1,n_irrep
         allocate(ham_tot_imag(i)%m(dims(i), dims(i)), STAT=status)
         ASSERT(status.eq.0)
      enddo

    case ( 'ham_xc' )

      allocate(ham_xc(n_irrep),STAT=status)
      ASSERT(status.eq.0)

      do i=1,n_irrep
         allocate(ham_xc(i)%m(dims(i), dims(i), n_spin), STAT=status)
         ASSERT(status.eq.0)
      enddo

    case ( 'ham_xc_real' )

      allocate(ham_xc_real(n_irrep),STAT=status)
      ASSERT(status.eq.0)

      do i=1,n_irrep
         allocate(ham_xc_real(i)%m(dims(i), dims(i)), STAT=status)
         ASSERT(status.eq.0)
      enddo

    case ( 'ham_xc_imag' )

      allocate(ham_xc_imag(n_irrep),STAT=status)
      ASSERT(status.eq.0)

      do i=1,n_irrep
         allocate(ham_xc_imag(i)%m(dims(i), dims(i)), STAT=status)
         ASSERT(status.eq.0)
      enddo

    case ( 'ham_lsft' )

      allocate(ham_lsft(n_irrep),STAT=status)
      ASSERT(status.eq.0)

      do i=1,n_irrep
         allocate(ham_lsft(i)%m(dims(i), dims(i), n_spin), STAT=status)
         ASSERT(status.eq.0)
      enddo

    case default
       print *,'ERROR: no such case >'//name//'<'
       ABORT('not implemented')
    end select
  end subroutine alloc_ham

  subroutine free_ham(char1, char2)
    !  Purpose: deallocate the requested part of the
    !           hamiltonian
    !           allowed values for char<i> are:
    !           "tot" and "lsft" (see select_free())
    !
    !           FN,10/95
    !------------ Declaration of formal parameters ---------------
    character(*),intent(in ) ::char1, char2
    optional :: char2
    !** End of interface *****************************************

    call select_free(char1)

    if(present(char2)) then
       call select_free(char2)
    endif
  end subroutine free_ham

  subroutine select_free(char)
    !
    !  Purpose: free hamiltonian
    !
    implicit none
    character(*),intent(in ) :: char
    integer :: i,alloc_stat
    !** End of interface *****************************************
    external error_handler
    !------------ Executable code ------------------------------------

    select case (char)
    case("tot")
       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          do i=1,size(ham_tot_real,1)
             deallocate(ham_tot_real(i)%m,ham_tot_imag(i)%m,STAT=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("FREE_HAM: deallocation (6) failed")
          enddo
          deallocate(ham_tot_real,ham_tot_imag,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("FREE_HAM: deallocation (6/6) failed")
       else ! options_spin_orbit
          do i=1,size(ham_tot,1)
             deallocate(ham_tot(i)%m,STAT=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("FREE_HAM: deallocation (6) failed")
          enddo
          deallocate(ham_tot,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("FREE_HAM: deallocation (6/6) failed")
       endif
    case("lsft")
     if(.not. options_spin_orbit) then
      do i=1,size(ham_lsft,1)
             deallocate(ham_lsft(i)%m,STAT=lsft_allocstat(2))
             ASSERT(lsft_allocstat(2).eq.0)
       enddo
             deallocate(ham_lsft,STAT=lsft_allocstat(1))
             ASSERT(lsft_allocstat(1).eq.0)
     endif

    case default
       call error_handler("FREE_HAM: unkown selection: "//char)
    end select

  end subroutine select_free


  subroutine hamiltonian_store(th)
    ! Purpose: stores the Kohn-Sham hamiltonian on a readwriteblocked
    !          file in case the 'save_ksmatrix' option is set.
    !
    ! routine called by: main_scf
    !** End of interface *****************************************
    !------------ Modules ----------------------------------------
    use readwriteblocked_module
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module     , only: output_main_scf, output_data_saved
    !------------ Declaration of formal parameters ---------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !------------ Declaration of local variables  ----------------
    integer(kind=i4_kind) :: i,s,j, n_spin, n_dim

    if (output_main_scf) call write_to_output_units &
         ("HAMILTONIAN_STORE: saving Kohn-Sham hamiltonian")

    if ( options_spin_orbit ) then
       WARN('SAVE_KSMATRIX not implemented for SO')
       RETURN
    endif

    n_spin = symmetry_data_n_spin()
    call readwriteblocked_write((/real(n_spin,r8_kind)/),th)
    if (output_data_saved) then
       write(output_unit,'(/ a     )')'Stored Kohn-Sham hamiltonian :'
       write(output_unit,'(  a     )')'n_spin'
       write(output_unit,'(4es20.13)')(/real(n_spin,r8_kind)/)
    endif
    do i=1,symmetry_data_n_irreps()
       n_dim = symmetry_data_dimension(i)
       do s=1,n_spin
          do j=1,n_dim
             call readwriteblocked_write(ham_tot(i)%m(:,j,s),th)
             if (output_data_saved) then
                1000 format('ham_tot(',i2,')%m(:,',i4,',',i1,')')
                write(output_unit, fmt = 1000 )i,j,s
                write(output_unit,'(4es20.13)')ham_tot(i)%m(:,j,s)
             endif
          enddo
       enddo
    enddo
  end subroutine hamiltonian_store

!*************************************************************
! record  1: n_spin
! FOREACH irrep i DO
! record  2: ham_tot(i)%m(:,:,1)
! record  3: ham_tot(i)%m(:,:,2)  [if n_spin > 1]
! DONE
!*************************************************************

  subroutine hamiltonian_recover(th)
    !
    ! Purpose:    recovers   the    Kohn-Sham    hamiltonian   on    a
    ! readwriteblocked file in case the 'read_ksmatrix' option is set.
    ! So far called from a  master only context in main_scf(). Assumes
    ! the storage for the  hamiltonian has already been allocated. See
    ! prescf_init(). FIXME: does not yet work for SO.
    !
    use readwriteblocked_module
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module, only: output_main_scf, output_data_read
    implicit none
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !** End of interface *****************************************

    allocatable           :: buffer
    integer(kind=i4_kind) :: i,s,j, n_spin, n_irrep, n_dim, spin_stored
    real(kind=r8_kind)    :: spin(1), buffer(:)
    real(r8_kind), parameter :: half = 0.5_r8_kind
    integer               :: status

    if (output_main_scf) call write_to_output_units &
         ("HAMILTONIAN_RECOVER: reading Kohn-Sham hamiltonian")

    if ( options_spin_orbit ) then
       WARN('READ_KSMATRIX not implemented for SO')
       RETURN
    endif

    ASSERT(allocated(ham_tot))

    n_spin = symmetry_data_n_spin()
    n_irrep = symmetry_data_n_irreps()

    call readwriteblocked_read(spin,th)
    if (output_data_read) then
       write(output_unit,'(/ a     )')'Recovered Kohn-Sham hamiltonian :'
       write(output_unit,'(  a     )')'n_spin'
       write(output_unit,'(4es20.13)')spin(1)
    endif
    spin_stored = int(spin(1),i4_kind)
    if (spin_stored > n_spin .and. output_main_scf) then
       call write_to_output_units &
            ("HAMILTONIAN_RECOVER: trying to convert spin-polarized data from")
       call write_to_output_units &
            ("                     tape into the spin-restricted data required")
       ! ham_tot(i)%m(:,:) := 1/2 * sum(s) ham_tot(i)%m(:,:,s)
    endif
    if (spin_stored < n_spin .and. output_main_scf) then
       call write_to_output_units &
            ("HAMILTONIAN_RECOVER: trying to convert spin-restricted data from")
       call write_to_output_units &
            ("                     tape into the spin-polarized data required")
       ! ham_tot(i)%m(:,:,s) := ham_tot(i)%m(:,:)
    endif

    do i=1,n_irrep
       n_dim = symmetry_data_dimension(i)

       do s=1,spin_stored
          if (s <= n_spin) then
             do j=1,n_dim
                call readwriteblocked_read(ham_tot(i)%m(:,j,s),th)
                if (output_data_read) then
                   1000 format('ham_tot(',i2,')%m(:,',i4,',',i1,')')
                   write(output_unit, fmt = 1000 )i,j,s
                   write(output_unit,'(4es20.13)')ham_tot(i)%m(:,j,s)
                endif
             enddo
          else ! s = spin_stored > n_spin = 1
             allocate(buffer(n_dim),STAT=status)
             if (status /= 0) call error_handler &
                  ("HAMILTONIAN_RECOVER: allocation (3) failed ")
             do j=1,n_dim
                call readwriteblocked_read(buffer,th)
                if (output_data_read) then
                   write(output_unit, fmt = 1000 )i,j,s
                   write(output_unit,'(4es20.13)')buffer
                endif
                ham_tot(i)%m(:,j,1) = ( ham_tot(i)%m(:,j,1) + buffer ) * half
             enddo
             deallocate(buffer,STAT=status)
             if (status /= 0) call error_handler &
                  ("HAMILTONIAN_RECOVER: deallocation (3) failed ")
          endif
       enddo
       if (n_spin > 1 .and. spin_stored == 1) then
          ham_tot(i)%m(:,:,2) = ham_tot(i)%m(:,:,1)
       endif
    enddo

  end subroutine hamiltonian_recover

  subroutine print_hamiltonian(char1,char2,char3,char4,char5,char6, &
       ham_xc_arr,ham_xc_arr_imag)
    !  Purpose: print out the selected part of the hamiltonian
    !           to the file $TTFSOUT/hamiltonian.
    !           allowed values for char<i> are:
    !           "kin", "nuc", "coul", "xc", "tot", and "full"
    !
    !           The keyword "full" indicates that all suceeding matrices
    !           in the parameter list are printed as quadratic matrices.
    !           This is the default if "direct_energy_calc" is turned off.
    !
    !           FN,1/96
    !------------ Modules ----------------------------------------
    !------------ Declaration of formal parameters ---------------
    character(*),intent(in ) ::char1,char2,char3,char4,char5,char6
    optional :: char2,char3,char4,char5,char6
    real(kind=r8_kind), optional :: ham_xc_arr(:)
    real(kind=r8_kind), optional :: ham_xc_arr_imag(:)
    !** End of interface *****************************************
    integer(kind=i4_kind)    :: count_loop
    logical                  :: quadratic
    data count_loop / 0 /

    !------------ Executable code ------------------------------------

    quadratic = .true. ! was /= etotmode_direct
    count_loop = count_loop + 1
    if(present(ham_xc_arr)) then
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             call select_print(char1,count_loop,quadratic,ham_xc_arr,ham_xc_arr_imag)
          else
             !
             ! STANDARD RUN
             !
             call select_print(char1,count_loop,quadratic,ham_xc_arr)
          endif
    else
       call select_print(char1,count_loop,quadratic)
    endif

    if(present(char2)) then
       if(present(ham_xc_arr)) then
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             call select_print(char2,count_loop,quadratic,ham_xc_arr,ham_xc_arr_imag)
          else
             !
             ! STANDARD RUN
             !
             call select_print(char2,count_loop,quadratic,ham_xc_arr)
          endif
       else
          call select_print(char2,count_loop,quadratic)
       endif
    else
       return
    endif

    if(present(char3)) then
       if(present(ham_xc_arr)) then
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             call select_print(char3,count_loop,quadratic,ham_xc_arr,ham_xc_arr_imag)
          else
             call select_print(char3,count_loop,quadratic,ham_xc_arr)
          endif
       else
          call select_print(char3,count_loop,quadratic)
       endif
    else
       return
    endif

    if(present(char4)) then
       if(present(ham_xc_arr)) then
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             call select_print(char4,count_loop,quadratic,ham_xc_arr,ham_xc_arr_imag)
          else
             call select_print(char4,count_loop,quadratic,ham_xc_arr)
          endif
       else
          call select_print(char4,count_loop,quadratic)
       endif
    else
       return
    endif

    if(present(char5)) then
       if(present(ham_xc_arr)) then
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             call select_print(char5,count_loop,quadratic,ham_xc_arr,ham_xc_arr_imag)
          else
             call select_print(char5,count_loop,quadratic,ham_xc_arr)
          endif
       else
          call select_print(char5,count_loop,quadratic)
       endif
    else
       return
    endif

    if(present(char6)) then
       if(present(ham_xc_arr)) then
          if (options_spin_orbit) then
             !
             ! SPIN ORBIT
             !
             call select_print(char6,count_loop,quadratic,ham_xc_arr,ham_xc_arr_imag)
          else
             call select_print(char6,count_loop,quadratic,ham_xc_arr)
          endif
       else
          call select_print(char6,count_loop,quadratic)
       endif
    else
       return
    endif

  end subroutine print_hamiltonian

  subroutine select_print(char,count,quadratic,ham_xc_arr,ham_xc_arr_imag)
    !  Purpose: print hamiltonian to file
    !           $TTFSOUT/hamiltonian. Works only
    !           for a system suited to the FORMAT specifiers
    !           In the most recent case this is H2.
    !
    !  subroutine called by :'print_hamiltonian'
    !
    !  FN,1/96
    !------------ Modules ----------------------------------------
    use iounitadmin_module
    use filename_module, only: outfile
    !------------ Declaration of formal parameters ---------------
    character(*),intent(in )          :: char
    integer(kind=i4_kind),intent(in)  :: count
    logical,intent(inout)             :: quadratic
    real(kind=r8_kind), optional      :: ham_xc_arr(:)
    ! in case of spin orbit
    real(kind=r8_kind), optional      :: ham_xc_arr_imag(:)
     !** End of interface *****************************************
    !------------ Declaration of local variables   ---------------
    integer(kind=i4_kind)   :: io_u,i_gamma,is,m,mm,&
         i_dim, j_dim, counter, counter0, i_spin
    real(kind=r8_kind), allocatable :: ham_num_xc(:,:),ham_num_xc_real(:,:),ham_num_xc_imag(:,:)
    logical :: use_ham_num_xc
    character(len=64) :: frmt
    !------------ Executable code ------------------------------------

    if (char == 'full') then
       quadratic = .true.
       return
    endif
    ! debug>>>
    quadratic = .true.

    i_spin = symmetry_data_n_spin()
    use_ham_num_xc = .false. ! was dependent on split_ham path...

    io_u=get_iounit()
    open(io_u,form='formatted',status='unknown', &
         position='append',&
         file=trim(outfile('hamiltonian')))

    write(io_u,*)' '
    write(io_u,*)' '
    write(io_u,*)' '
    write(io_u,*)'+++++++++++++ Loop ',count,' ++++++++++++++++'
    write(io_u,*)' Part of hamiltonian: ',char

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       counter0 = 0
       counter = 0
       do i_gamma=1,symmetry_data_n_proj_irreps()
          i_dim=symmetry_data_dimension_proj(i_gamma)
          write(io_u,*)' '
          write(io_u,*)' '
          write(io_u,*)'---------- Irrep ',i_gamma,' -------------'

          ! XC-ham first if present:
          if (use_ham_num_xc) then
             allocate(ham_num_xc_real(i_dim,i_dim),ham_num_xc_imag(i_dim,i_dim),STAT=is)
             if (is /= 0) call error_handler &
                  ("PRINT_HAMILTONIAN: allocation failed")

             do m=1,i_dim
                do mm=1,m
                   counter = counter + 1
                   ham_num_xc_real(m,mm) = ham_xc_arr(counter)
                   ham_num_xc_imag(m,mm) = ham_xc_arr_imag(counter)
                   ham_num_xc_real(mm,m) = ham_xc_arr(counter)
                   ham_num_xc_imag(mm,m) = ham_xc_arr_imag(counter)
                enddo
             enddo
          endif


          frmt = '(4(A,F16.13))'
             do m=1,i_dim
                if (quadratic) then
                   j_dim = i_dim
                else
                   j_dim = m
                endif
                select case (char)
                case ("tot")
                   write(io_u,'(10(A,G12.5E2))')&
                        ( ' ',ham_tot_real(i_gamma)%m(mm,m),mm=1,j_dim)
                   write(io_u,'(10(A,G12.5E2))')&
                        ( ' ',ham_tot_imag(i_gamma)%m(mm,m),mm=1,j_dim)
                case default
                   write(io_u,*)"PRINT_HAM: unknown selection: "//char
                end select
             enddo

!!$          enddo
          if (use_ham_num_xc) then
             deallocate(ham_num_xc_real,ham_num_xc_imag)
          endif
          counter0 = counter0 + i_spin*(i_dim*(i_dim+1))/2
       enddo
    else ! options_spin_orbit

       counter0 = 0
       do i_gamma=1,symmetry_data_n_irreps()
          i_dim=symmetry_data_dimension(i_gamma)
          write(io_u,*)' '
          write(io_u,*)' '
          write(io_u,*)'---------- Irrep ',i_gamma,' -------------'
          if (use_ham_num_xc) then
             allocate(ham_num_xc(i_dim,i_dim),STAT=is)
             if (is /= 0) call error_handler &
                  ("PRINT_HAMILTONIAN: allocation failed")
          endif
          do is=1,i_spin
             counter = counter0 + is
             write(io_u,*)' Spin : ',is
             write(io_u,*)' '
             if (use_ham_num_xc) then
                do m=1,i_dim
                   do mm=1,m
                      ham_num_xc(m,mm) = ham_xc_arr(counter)
                      ham_num_xc(mm,m) = ham_xc_arr(counter)
                      counter=counter+i_spin
                   enddo
                enddo
             endif
             do m=1,i_dim
                if (quadratic) then
                   j_dim = i_dim
                else
                   j_dim = m
                endif
                select case (char)

                case ("tot")
!!$                   write(io_u,'(4(A,F12.8))')&
!!$                        ( ' ',ham_tot(i_gamma)%m(mm,m,is),mm=1,j_dim)
                   write(io_u,'(10(A,G12.5E2))')&
                        ( ' ',ham_tot(i_gamma)%m(mm,m,is),mm=1,j_dim)
                case default
                   write(io_u,*)"PRINT_HAM: unknown selection: "//char
                end select
             enddo

          enddo
          if (use_ham_num_xc) then
             deallocate(ham_num_xc)
          endif
          counter0 = counter0 + i_spin*(i_dim*(i_dim+1))/2
       enddo
    endif! options_spin_orbit
    write(io_u,*)' '
    write(io_u,*)'+++++++++++ End of loop ',count,' +++++++++++++'
    close(io_u)
    call return_iounit(io_u)

  end subroutine select_print

  subroutine hamiltonian_shutdown()
    !
    ! Purpose:  deallocate all  parts of  hamiltonian (that  are still
    ! allocated) at the end of the SCF-cycles.  Called from a parallel
    ! context in prescf_finalize()
    !
    use comm, only: comm_rank
    implicit none
    !** End of interface *****************************************

    ! moved from inside of SCF loop:
    call free_ham('tot')

    ! moved from after the SCF loop:
    if(allocated(ham_lsft)) call free_ham('lsft')

    call print_lsft_diis_allocstat()
  end subroutine hamiltonian_shutdown

  subroutine print_lsft_diis_allocstat()
    implicit none
    ! *** end of interface ***

    integer(kind=i4_kind) i_al

    do i_al=1,size(lsft_allocstat)
       if(lsft_allocstat(i_al).ne.0) print*,i_al, 'error lsft_allocstat ne.0'
    enddo
  end subroutine print_lsft_diis_allocstat

!--------------- End of module ----------------------------------
  end module hamiltonian_module
