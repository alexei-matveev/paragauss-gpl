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
!================================================================
! Public interface of module
!================================================================
subroutine  integral_calc_quad_dipole()
  !----------------------------------------------------------------
  !
  !  Purpose: This routine is the main routine for the
  !           integral calculation of one quadrupel
  !           ( unique atom 1, l 1, unique atom 2, l 2 )
  !           for dipole integral calculation.
  !           The corresponding data (including
  !           the quadrupel to be calculated) are stored in
  !           int_data_dipole_module.
  !           Contraction and symmetry-adaption are directly
  !           done within this subroutine.
  !
  !  Subroutine called by: main_slave integral_main_dipole
  !  It contains: symadapt_add_nottotalsym_gten, 
  !                symadapt_add_nottotalsym_sporgt
  !                 renorm_spor_gt
  !                  and other procedures
  !
  !  References: Publisher Document: Concepts of Integral Part
  ! 
  !  Author: TB
  !  Date: 7/97
  !
  !================================================================
  ! End of public interface of module
  !================================================================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification 
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification 
  ! Author: DG
  ! Date:   20/03/2001
  ! Description: I must not do contractions now in case if kinematic factors is on in the gtensor stuff
  ! I will do it late in gtensor module (dipole_module.f90?)
  !
  ! Modification use DLB instead of master/slave for dipole integrals
  ! Author: AN
  ! Date:   4/11
  ! Description: use DLB for scheduling batches of dipole integrals
  !              remove reporting back of slaves, they get their new
  !              tasks via DLB
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------

  !------------ Modules used ------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use timer_module
  use time_module
  use integralpar_module
  use int_data_dipole_module
  use output_module
  use int_send_dipole_module
  use iounitadmin_module
  use operations_module, only: operations_gtensor, operations_dipole, operations_hfc
  use options_module, only:options_kinematic_factors
  implicit none
  !------------ Declaration of local variables ------------------
  integer(kind=i4_kind)  :: i_ea1, i_ea2 ! loop indices for equal atoms

  if ( output_int_loops .or. output_int_taskdistribution) then
     write(output_unit,*) "integral_calc_quad_dipole: start with quadrupel ", &
          quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
     write(stdout_unit,*) "integral_calc_quad_dipole: start with quadrupel ", &
          quadrupel%ua1, quadrupel%l1, quadrupel%ua2, quadrupel%l2
  endif

  if ( output_int_loops ) call write_to_output_units( &
       "integral_calc_quad_dipole: setup")
  call setup()

  eqal_atom_1: do i_ea1 = 1, ua1%N_equal_atoms
     eqal_atom_2: do i_ea2 = 1, ua2%N_equal_atoms

        if ( output_int_loops ) then
           write(output_unit,*) &
                "integral_calc_quad_dipole: processing pair of equal atoms ", &
                i_ea1, i_ea2
           write(stdout_unit,*) &
                "integral_calc_quad_dipole: processing pair of equal atoms ", &
                i_ea1, i_ea2
        endif

        center1 = ua1%position(:,i_ea1)
        center2 = ua2%position(:,i_ea2)

        ! calculate primitive integrals and contract them 
        if ( output_int_loops ) call write_to_output_units( &
             "integral_calc_quad_dipole: calc_primitives_and_contract")
        call calc_primitives_and_contract(quadrupel%ua1,i_ea1,quadrupel%l1 &
                                         ,quadrupel%ua2,i_ea2,quadrupel%l2)


        ! add contribution of one pair
        ! of equal atom to symmetry adapted integrals
        if ( output_int_loops ) call write_to_output_units( &
             "integral_calc_quad_dipole:&
             & adding results to symmetry adapted arrays")

        
        if( operations_dipole )then
           if( options_spin_orbit )then
              call symadapt_add_nottotalsym_spor()
           else
              call symadapt_add_nottotalsym()
           endif
        endif

#ifdef WITH_GTENSOR
        if( operations_gtensor )then
           ASSERT(options_spin_orbit)
           call symadapt_add_nottotalsym_sporgt()
        endif

        if( operations_hfc )then
           if (options_spin_orbit) then
              call symadapt_add_nottotalsym_sohfc()
           else
              call symadapt_add_nottotalsym_hfc()
           endif
        endif
#endif

        ! deallocate contracted integrals
        if ( output_int_loops ) call write_to_output_units( &
             "integral_calc_quad_dipole: deallocate_contracted")
        call deallocate_contracted() ! L291 dipol and gten

        if ( output_int_loops ) call write_to_output_units( &
             "integral_calc_quad_dipole: 1 pair of equal atoms done")

     end do eqal_atom_2
  end do eqal_atom_1

#if 0
  ! write debug output
  if ( output_int_data ) call symadapt_write()
#endif

  if ( output_int_loops ) call write_to_output_units( &
       "integral_calc_quad_dipole: calculation done ")

  if ( output_int_loops ) call write_to_output_units( &
       "integral_calc_quad_dipole: int_send_dipole_send")
  call int_send_dipole_send()


  if ( output_int_loops ) call write_to_output_units( &
       "integral_calc_quad_dipole: shutdown")
  call shutdown()


  if ( output_int_loops ) call write_to_output_units( &
       "integral_calc_quad_dipole: end")

  !--------------------------------------------------------------
  !------------ Private Subroutines -----------------------------
contains

  !*************************************************************
  subroutine setup
    !  Purpose: contains setup routines of various modules
    !------------ Executable code ------------------------------

    call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
    call init_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_quadrupel_2cob3c(integralpar_i_int_part))

    call int_data_dipole_setup()

  end subroutine setup
  !**************************************************************
  subroutine shutdown
    !  Purpose: contains shutdown routines of various modules
    !------------ Executable code ------------------------------

    call int_data_dipole_shutdown()

    call stop_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
    call stop_timer(timer_int_quadrupel_2cob3c(integralpar_i_int_part))
    call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
    call timer_small_to_large( &
         timer_int_calc_2cob3c(integralpar_i_int_part), &
         timer_int_calcsum_2cob3c(integralpar_i_int_part) )

  end subroutine shutdown
  !**************************************************************


  !**************************************************************
  subroutine allocate_primitives()
    ! Purpose: allocates storage for primitive integrals
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: status
    !------------ Executable code ------------------------------
    integralpar: if ( integralpar_2cob_dipole ) then
       dip:  if(operations_dipole) then
          allocate( prim_int_2cob_dipole(n_exp2,n_exp1,n_m2,n_m1,3), &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: allocate_primitives: 2cob_dipole")
          prim_int_2cob_dipole = 0.0_r8_kind
       end if dip

#ifdef WITH_GTENSOR
       gt:if (operations_gtensor) then
          allocate( prim_int_2cob_dipoleg(n_exp2,n_exp1,n_m2,n_m1,13), &    ! 4 + 9 DG
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: allocate_primitives: 2cob_dipoleg")
          prim_int_2cob_dipoleg = 0.0_r8_kind
          kf:   if(options_kinematic_factors) then
             allocate( cont_int_2cob_dipoleg(n_exp2,n_exp1,n_m2,n_m1,13), &  
                  ! 4 + 9 DG in case if kinematic factors I mustn`t run 
                  stat=status ) 
             if ( status .ne. 0 ) call error_handler( &
                  "integral_calc_quad_dipole: allocate_primitives: 2cob_dipoleg")
             !contract procedure, I`ll do it later in gtensor module
             cont_int_2cob_dipoleg = 0.0_r8_kind
          end if kf
       end if gt

       hfc:if (operations_hfc) then
          allocate( prim_int_2cob_hfc(n_exp2,n_exp1,n_m2,n_m1,10,n_unique_atoms), &    ! DG
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: allocate_primitives: 2cob_hfc")
          prim_int_2cob_hfc = 0.0_r8_kind
          kf_hfc:   if(options_kinematic_factors) then
             allocate( cont_int_2cob_hfc(n_exp2,n_exp1,n_m2,n_m1,10,n_unique_atoms), &     ! 
                  stat=status ) 
             if ( status .ne. 0 ) call error_handler( &
                  "integral_calc_quad_dipole: allocate_primitives: 2cob_hfc") !
             cont_int_2cob_hfc = 0.0_r8_kind
          end if kf_hfc
       end if hfc
#endif


    endif integralpar
  end subroutine allocate_primitives
  !**************************************************************


  !**************************************************************
  subroutine deallocate_contracted()
    ! Purpose: deallocates storage of contracted integrals
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: status
    !------------ Executable code -------------------------------
    integralpar: if (integralpar_2cob_dipole ) then
       if (operations_dipole) then
          deallocate( cont_int_2cob_dipole, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: deallocate_contracted: 2cob_dipole")
       endif

#ifdef WITH_GTENSOR
       if (operations_gtensor ) then
          deallocate( cont_int_2cob_dipoleg, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: deallocate_contracted: 2cob_dipoleg")
       endif

       if (operations_hfc ) then
          deallocate( cont_int_2cob_hfc, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: deallocate_contracted: 2cob_hfc")
       endif
#endif
    end if integralpar

  end subroutine deallocate_contracted
  !**************************************************************

 !**************************************************************
  subroutine deallocate_primitivies()
    ! Purpose: deallocates storage of contracted integrals
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: status
    !------------ Executable code -------------------------------
    integralpar: if (integralpar_2cob_dipole ) then
       if (operations_dipole) then
          deallocate( prim_int_2cob_dipole, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: deallocate_primitivies: 2cob_dipole")
       endif

#ifdef WITH_GTENSOR
       if (operations_gtensor ) then
          deallocate( prim_int_2cob_dipoleg, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: deallocate_primitivies: 2cob_dipoleg")
       endif

        if (operations_hfc ) then
          deallocate( prim_int_2cob_hfc, &
               stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integral_calc_quad_dipole: deallocate_primitivies: 2cob_hfc")
       endif
#endif
    end if integralpar

  end subroutine deallocate_primitivies
  !**************************************************************


  !**************************************************************
  subroutine contract_2center(prim,cont,n_xyz)
    !  Purpose: does contraction on two center integrals
    !  the result is in cont that is allocated by this subroutine
    !  prim is deallocated
    !------------ Declaration of formal parameters -------------
    real(kind=r8_kind), pointer, dimension(:,:,:,:,:) :: prim,cont
    integer(kind=i4_kind), intent(in) :: n_xyz
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind), allocatable :: int(:)
    real(kind=r8_kind)    :: coef
    integer(kind=i4_kind) :: i_m1, i_m2, i_exp1, i_exp2, i_c1, &
         i_c2, i_co1, i_co2, status, i_xyz
    !------------ Executable code ------------------------------
    allocate( int(n_c2), stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_dipole: contract_2center: allocate of int failed")
    allocate( cont(n_c2,n_c1,n_m2,n_m1,n_xyz), stat=status )                  !Allocate contracted integrals
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_dipole: contract_2center: allocate of cont failed")
    cont(:,n_uncontr1+1:,:,:,:) = 0.0_r8_kind
    xyz: do i_xyz = 1, n_xyz
       m1: do i_m1 = 1, n_m1
          m2: do i_m2 = 1, n_m2
             exp1_uc: do i_exp1 = 1, n_uncontr1
                do i_c2 =1, n_uncontr2
                   int(i_c2) = prim(i_c2,i_exp1,i_m2,i_m1,i_xyz)
                enddo
                i_co2 = 1
                do i_c2 = n_uncontr2 + 1, n_c2
                   int(i_c2) = 0.0_r8_kind
                   do i_exp2 = 1, n_exp2
                      int(i_c2) = int(i_c2) + &
                           contractions2(i_exp2,i_co2) * &
                           prim(i_exp2,i_exp1,i_m2,i_m1,i_xyz)
                   enddo
                   i_co2 = i_co2 + 1
                enddo
                i_co1 = 1
                do i_c1= n_uncontr1 + 1, n_c1
                   coef = contractions1(i_exp1,i_co1)
                   do i_c2 = 1, n_c2
                      cont(i_c2,i_c1,i_m2,i_m1,i_xyz) = &
                           cont(i_c2,i_c1,i_m2,i_m1,i_xyz) + &
                           coef * int(i_c2)
                   enddo
                   i_co1 = i_co1 + 1
                enddo
                do i_c2 = 1, n_c2
                   cont(i_c2,i_exp1,i_m2,i_m1,i_xyz) = int(i_c2)
                enddo
             enddo exp1_uc
             exp1_c: do i_exp1 = n_uncontr1 + 1, n_exp1
                do i_c2 =1, n_uncontr2
                   int(i_c2) = prim(i_c2,i_exp1,i_m2,i_m1,i_xyz)
                enddo
                i_co2 = 1
                do i_c2 = n_uncontr2 + 1, n_c2
                   int(i_c2) = 0.0_r8_kind
                   do i_exp2 = 1, n_exp2
                      int(i_c2) = int(i_c2) + &
                           contractions2(i_exp2,i_co2) * &
                           prim(i_exp2,i_exp1,i_m2,i_m1,i_xyz)
                   enddo
                   i_co2 = i_co2 + 1
                enddo
                i_co1 = 1
                do i_c1= n_uncontr1 + 1, n_c1
                   coef = contractions1(i_exp1,i_co1)
                   do i_c2 = 1, n_c2
                      cont(i_c2,i_c1,i_m2,i_m1,i_xyz) = &
                           cont(i_c2,i_c1,i_m2,i_m1,i_xyz) + &
                           coef * int(i_c2)
                   enddo
                   i_co1 = i_co1 + 1
                enddo
             enddo exp1_c
          enddo m2
       enddo m1
    enddo xyz
    !deallocation of primitivies removed to deallocate_primitivies() DG

    deallocate( int, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_dipole: contract_2center: deallocate of int failed")
  end subroutine contract_2center
  !**************************************************************
  !**************************************************************
  subroutine contract_2center_nua(prim,cont,n_xyz)
    !  Purpose: does contraction on two center integrals
    !  the result is in cont that is allocated by this subroutine
    !  prim is deallocated for n_unique_atoms
    !------------ Declaration of formal parameters -------------
    real(kind=r8_kind), pointer, dimension(:,:,:,:,:,:) :: prim,cont
    integer(kind=i4_kind), intent(in) :: n_xyz
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind), allocatable :: int(:)
    real(kind=r8_kind)    :: coef
    integer(kind=i4_kind) :: i_m1, i_m2, i_exp1, i_exp2, i_c1, &
         i_c2, i_co1, i_co2, status, i_xyz, i_ua
    !------------ Executable code ------------------------------
    allocate( int(n_c2), stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_dipole: contract_2center: allocate of int failed")
    allocate( cont(n_c2,n_c1,n_m2,n_m1,n_xyz,n_unique_atoms), stat=status )
    !Allocate contracted integrals
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_dipole: contract_2center: allocate of cont failed")
    
    cont(:,n_uncontr1+1:,:,:,:,:) = 0.0_r8_kind
    loopua: do i_ua = 1,n_unique_atoms
    xyz: do i_xyz = 1, n_xyz
       m1: do i_m1 = 1, n_m1
          m2: do i_m2 = 1, n_m2
             exp1_uc: do i_exp1 = 1, n_uncontr1
                do i_c2 =1, n_uncontr2
                   int(i_c2) = prim(i_c2,i_exp1,i_m2,i_m1,i_xyz,i_ua)
                enddo
                i_co2 = 1
                do i_c2 = n_uncontr2 + 1, n_c2
                   int(i_c2) = 0.0_r8_kind
                   do i_exp2 = 1, n_exp2
                      int(i_c2) = int(i_c2) + &
                           contractions2(i_exp2,i_co2) * &
                           prim(i_exp2,i_exp1,i_m2,i_m1,i_xyz,i_ua)
                   enddo
                   i_co2 = i_co2 + 1
                enddo
                i_co1 = 1
                do i_c1= n_uncontr1 + 1, n_c1
                   coef = contractions1(i_exp1,i_co1)
                   do i_c2 = 1, n_c2
                      cont(i_c2,i_c1,i_m2,i_m1,i_xyz,i_ua) = &
                           cont(i_c2,i_c1,i_m2,i_m1,i_xyz,i_ua) + &
                           coef * int(i_c2)
                   enddo
                   i_co1 = i_co1 + 1
                enddo
                do i_c2 = 1, n_c2
                   cont(i_c2,i_exp1,i_m2,i_m1,i_xyz,i_ua) = int(i_c2)
                enddo
             enddo exp1_uc
             exp1_c: do i_exp1 = n_uncontr1 + 1, n_exp1
                do i_c2 =1, n_uncontr2
                   int(i_c2) = prim(i_c2,i_exp1,i_m2,i_m1,i_xyz,i_ua)
                enddo
                i_co2 = 1
                do i_c2 = n_uncontr2 + 1, n_c2
                   int(i_c2) = 0.0_r8_kind
                   do i_exp2 = 1, n_exp2
                      int(i_c2) = int(i_c2) + &
                           contractions2(i_exp2,i_co2) * &
                           prim(i_exp2,i_exp1,i_m2,i_m1,i_xyz,i_ua)
                   enddo
                   i_co2 = i_co2 + 1
                enddo
                i_co1 = 1
                do i_c1= n_uncontr1 + 1, n_c1
                   coef = contractions1(i_exp1,i_co1)
                   do i_c2 = 1, n_c2
                      cont(i_c2,i_c1,i_m2,i_m1,i_xyz,i_ua) = &
                           cont(i_c2,i_c1,i_m2,i_m1,i_xyz,i_ua) + &
                           coef * int(i_c2)
                   enddo
                   i_co1 = i_co1 + 1
                enddo
             enddo exp1_c
          enddo m2
       enddo m1
    enddo xyz
 end do loopua
    !deallocation of primitivies removed to deallocate_primitivies() DG

    deallocate( int, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integral_calc_quad_dipole: contract_2center_nua: deallocate of int failed")
  end subroutine contract_2center_nua
  !**************************************************************

  !**************************************************************
  subroutine calc_primitives_and_contract(U1,E1,L1,U2,E2,L2)
    ! calculate primitive integrals and do contraction
    use shgi_dip, only: shgi_dip_drv
    implicit none
    integer(i4_kind), intent(in) :: U1,E1,L1,U2,E2,L2
    ! *** end of interface ***

    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind)  :: i_ua
!!$    integer(kind=i4_kind)  :: i_c1, i_c2, i_exp1, i_exp2
    ! help pointers
    !------------ Executable code ------------------------------

    ! allocate primitive integrals
    call allocate_primitives()

    ! calculate primitive integrals
    call start_timer(timer_int_prim_2cob3c(integralpar_i_int_part))

    if (operations_dipole) then
       ! calculate primitive integrals
       call shgi_dip_drv(U1,E1,L1,U2,E2,L2,unique_atoms &
                    , prim_int_2cob_dipole)
    endif

#ifdef WITH_GTENSOR
    gt:if(operations_gtensor) then
       call ll_calculate_dipoleg(U1,U2,L1,L2)
    end if gt

    hfc:if(operations_hfc) then
       call ll_calculate_hfc(U1,U2,L1,L2)
    end if hfc
#endif

    if ( output_int_loops ) call write_to_output_units( &
         "integral_calc_quad_dipole: calc_primitives_and_contract:&
         &primitive integrals done")

    call stop_timer(timer_int_prim_2cob3c(integralpar_i_int_part))

#if 0
    if ( output_int_data  ) then
       ! write primitive 2-center integrals
       write(debug_unit,*)
       write(debug_unit,*)
       write(debug_unit,*) "######### primitive 2-center integrals ##########"
       if (integralpar_2cob_dipole) then
          write(debug_unit,*) "i_exp1, i_exp2, dipole"
          do i_exp1 = 1,n_exp1
             do i_exp2 = 1,n_exp2
                write(debug_unit,'(2I4,100F20.8)') i_exp1, i_exp2, &
                     prim_int_2cob_dipole(i_exp2,i_exp1,:,:,:)   
             enddo
          enddo
       endif
       write(debug_unit,*) "######### primitive 2-center integrals end ##########"
       write(debug_unit,*)
       write(debug_unit,*)
       if ( output_int_loops ) call write_to_output_units( &
            "integral_calc_quad_dipole: write 2 center primitives done")
    endif
#endif
    
    ! do contraction
    call start_timer(timer_int_cont_2cob3c(integralpar_i_int_part))

    if ( output_int_loops ) call write_to_output_units( &
         "integral_calc_quad_dipole: calc_primitives_and_contract:&
         & contract 2cob_dipole")
    dip: if (operations_dipole) then

       call contract_2center(prim_int_2cob_dipole,cont_int_2cob_dipole,3)

    end if dip

#ifdef WITH_GTENSOR
    gten:if (operations_gtensor) then

       kf:if(options_kinematic_factors) then
          DPRINT  'don"t cotract now'
          cont_int_2cob_dipoleg(:,:,:,:,:) = prim_int_2cob_dipoleg(:,:,:,:,:)
          !DG I`m passing the array directly, w/o contraction
          !Here   I must deallocate primitives integrals because
          !I don`t like to change the name of array later in all
          !normaly they are deallocated in contact_2centre
          !places, so I remapping it one to one directly
       else
          DPRINT ' cotract now'
          call contract_2center(prim_int_2cob_dipoleg,cont_int_2cob_dipoleg,13) !DG
       end if kf
    end if gten

    hfc_c:if (operations_hfc) then

       hfckf:if(options_kinematic_factors) then
          DPRINT  'don"t cotract now'
          cont_int_2cob_hfc(:,:,:,:,:,:) = prim_int_2cob_hfc(:,:,:,:,:,:)
          !DG I`m passing the array directly, w/o contraction
          !Here   I must deallocate primitives integrals because
          !I don`t like to change the name of array later in all
          !normaly they are deallocated in contact_2centre
          !places, so I remapping it one to one directly
          ! this is not a  good idea, I showld  have just to change pointer
          ! to avoid data transfer
       else
          DPRINT " cotract now"
          Do i_ua = 1, n_unique_atoms

             call contract_2center_nua(prim_int_2cob_hfc,cont_int_2cob_hfc, 10 ) !DG
          end Do
       end if hfckf
    end if hfc_c
#endif

 !deallocation of primitivies integrals
 call deallocate_primitivies()

    call stop_timer(timer_int_cont_2cob3c(integralpar_i_int_part))

    if ( output_int_loops ) call write_to_output_units( &
         "integral_calc_quad_2cob3c: calc_primitives_and_contract:&
         & contraction done")

#if 0
    if ( output_int_data ) then
       ! write contracted integrals
       write(debug_unit,*)
       write(debug_unit,*)
       write(debug_unit,*) "######### comtracted 2-center integrals ##########"
       if (integralpar_2cob_dipole) then
          write(debug_unit,*) "i_c1, i_c2, dipole"
          do i_c1 = 1,n_c1
             do i_c2 = 1,n_c2
                write(debug_unit,'(2I4,100F20.8)') i_c1, i_c2, &
                     cont_int_2cob_dipole(i_c2,i_c1,:,:,:)   
             enddo
          enddo
       endif
       write(debug_unit,*) "######### contracted 2-center integrals end ##########"
       write(debug_unit,*)
       write(debug_unit,*)
       if ( output_int_loops ) call write_to_output_units( &
            "integral_calc_quad_dipole: write done")
    endif
#endif

  end subroutine calc_primitives_and_contract
  !**************************************************************


  !**************************************************************
  subroutine symadapt_add_nottotalsym()
    ! add contribution of one pair
    ! of equal atom to symmetry adapted integrals
    ! the case of not total symmetric integrals
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
         saint_2cob_dipole
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
         i_ip1, i_if1, i_if2, i_cf1, i_cf2, m1, m2
    integer(kind=i4_kind) :: i_ip2    
    real(kind=r8_kind) :: coef1, coef
    !------------ Executable code ------------------------------

    if(.not.integralpar_2cob_dipole) RETURN ! FIXME: dont enter in the first line

    call start_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
!   dipole: if ( integralpar_2cob_dipole ) then
       xyz: do i_xyz = 1, 3
          i_ip1 = 1
          irrep1: do i_ir1 = 1, symmetry_data_n_irreps()
             sap1 => ua1%symadapt_partner(i_ir1,quadrupel%l1)
             partner1: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                i_ip2 = 1
                irrep2: do i_ir2 = 1, i_ir1
                   sap2 => ua2%symadapt_partner(i_ir2,quadrupel%l2)
                   partner2: do i_pa2 = 1, symmetry_data_n_partners(i_ir2)
                      needed: if ( associated( &
                           symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int ) ) then
                         saint_2cob_dipole => symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int
                         ind_fct_1: do i_if1 = 1, sap1%N_independent_fcts
                            sat1 => sap1%sa_int(i_ea1,i_if1,i_pa1)
                            ind_fct_2: do i_if2 = 1, sap2%N_independent_fcts
                               sat2 => sap2%sa_int(i_ea2,i_if2,i_pa2)
                               m_sum_1: do i_cf1 = 1, sat1%N_fcts
                                  m1 = sat1%m(i_cf1)
                                  coef1 = sat1%c(i_cf1)
                                  m_sum_2: do i_cf2 = 1, sat2%N_fcts
                                     m2 = sat2%m(i_cf2)
                                     coef = sat2%c(i_cf2) * coef1
                                     saint_2cob_dipole(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipole(:,:,i_if2,i_if1) + &
                                          coef * cont_int_2cob_dipole(:,:,m2,m1,i_xyz)
                                  enddo m_sum_2
                               enddo m_sum_1
                            enddo ind_fct_2
                         enddo ind_fct_1
                      endif needed
                      i_ip2 = i_ip2 + 1
                   enddo partner2
                enddo irrep2
                i_ip1 = i_ip1 + 1
             enddo partner1
          enddo irrep1
       enddo xyz
!   endif dipole
    call stop_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
  end subroutine symadapt_add_nottotalsym
  !**************************************************************

#ifdef WITH_GTENSOR
  !**************************************************************
  subroutine symadapt_add_nottotalsym_hfc()
    ! add contribution of one pair
    ! of equal atom to symmetry adapted integrals
    ! the case of not total symmetric integrals
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
         saint_2cob_hfc
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
         i_ip1, i_if1, i_if2, i_cf1, i_cf2, m1, m2, i_ua
    real(kind=r8_kind) :: coef1, coef
    !------------ Executable code ------------------------------

    if(.not.integralpar_2cob_dipole) RETURN ! FIXME: dont enter at all

    call start_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
!   neednrhfc: if ( integralpar_2cob_dipole ) then
      nua : do i_ua = 1, n_unique_atoms
        xyznr:   do i_xyz = 1, 7
              i_ip1 = 1
              do i_ir1 = 1, symmetry_data_n_irreps()
                 sap1 => ua1%symadapt_partner(i_ir1,quadrupel%l1)
                 do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                    i_ir2 = i_ir1
                    i_pa2 = i_pa1
                 
                       sap2 => ua2%symadapt_partner(i_ir2,quadrupel%l2)
                     
                          if ( associated( &
                               symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz, i_ua)%int ) ) then
                             
                             saint_2cob_hfc => symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz, i_ua)%int

                             select case(i_xyz)
                                
                                
                             case(7)! Aiso  10->7
                             do i_if1 = 1, sap1%N_independent_fcts
                                sat1 => sap1%sa_int(i_ea1,i_if1,i_pa1)
                                do i_if2 = 1, sap2%N_independent_fcts
                                   sat2 => sap2%sa_int(i_ea2,i_if2,i_pa2)
                                   do i_cf1 = 1, sat1%N_fcts
                                      m1 = sat1%m(i_cf1)
                                      coef1 = sat1%c(i_cf1)
                                      do i_cf2 = 1, sat2%N_fcts
                                         m2 = sat2%m(i_cf2)
                                         coef = sat2%c(i_cf2) * coef1
                                         saint_2cob_hfc(:,:,i_if2,i_if1) = &
                                              saint_2cob_hfc(:,:,i_if2,i_if1) + &
                                              coef * cont_int_2cob_hfc(:,:,m2,m1,10, i_ua)
                                      enddo
                                   enddo
                                enddo
                             enddo

                             case(1:6)!  1->1, 2->2 etc. See description in the end of ll_calculate_hfc
                             do i_if1 = 1, sap1%N_independent_fcts
                                sat1 => sap1%sa_int(i_ea1,i_if1,i_pa1)
                                do i_if2 = 1, sap2%N_independent_fcts
                                   sat2 => sap2%sa_int(i_ea2,i_if2,i_pa2)
                                   do i_cf1 = 1, sat1%N_fcts
                                      m1 = sat1%m(i_cf1)
                                      coef1 = sat1%c(i_cf1)
                                      do i_cf2 = 1, sat2%N_fcts
                                         m2 = sat2%m(i_cf2)
                                         coef = sat2%c(i_cf2) * coef1
                                         saint_2cob_hfc(:,:,i_if2,i_if1) = &
                                              saint_2cob_hfc(:,:,i_if2,i_if1) + &
                                              coef * cont_int_2cob_hfc(:,:,m2,m1,i_xyz, i_ua)
                                      enddo
                                   enddo
                                enddo
                             enddo

                          end select
                       end if
                       i_ip1 = i_ip1 + 1
                    enddo
              
               
             enddo

          enddo xyznr
          
       end do nua
!   endif neednrhfc

    call stop_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
  end subroutine symadapt_add_nottotalsym_hfc
#endif

  subroutine symadapt_add_nottotalsym_spor()
    ! add contribution of one pair
    ! of equal atom to symmetry adapted integrals
    ! the case of not total symmetric integrals
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
         saint_2cob_dipole_real,saint_2cob_dipole_imag
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
         i_ip1, i_ip2, i_if1, i_if2, i_cf1, i_cf2, m1, m2,alpha
    real(kind=r8_kind) :: coef1_real, coef2_real, coef_real, coef1_imag, coef2_imag, coef_imag
    real(kind=r8_kind) :: coef_real_arr(9,9),coef_imag_arr(9,9)
    !------------ Executable code ------------------------------

    if(.not.integralpar_2cob_dipole) RETURN ! FIXME: dont enter at all

    call start_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
!   dipole: if ( integralpar_2cob_dipole ) then
       xyz: do i_xyz = 1, 3
          i_ip1 = 1
          irrep1: do i_ir1 = 1, symmetry_data_n_proj_irreps()
             sap1 => ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)
             partner1: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                i_ip2 = 1
                irrep2: do i_ir2 = 1, i_ir1
                   sap2 => ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)
                   partner2: do i_pa2 = 1, symmetry_data_n_partners_proj(i_ir2)
                      needed: if ( associated( &
                           symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int ) ) then
                         saint_2cob_dipole_real => symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int
                         saint_2cob_dipole_imag => symadapt_int_2cob_dipole_p(i_ip2,i_ip1,i_xyz)%int_imag
                         ind_fct_1: do i_if1 = 1, sap1%N_independent_fcts
                            ind_fct_2: do i_if2 = 1, sap2%N_independent_fcts
                               coef_real_arr = 0.0_r8_kind
                               coef_imag_arr = 0.0_r8_kind
                               spin_up_and_down: do alpha =1,2
                                  sat1 => sap1%sa_spor_int(i_ea1,alpha,i_if1,i_pa1)
                                  sat2 => sap2%sa_spor_int(i_ea2,alpha,i_if2,i_pa2)
                                  m_sum_1: do i_cf1 = 1, sat1%N_fcts
                                     m1 = sat1%m(i_cf1)
                                     coef1_real = sat1%re(i_cf1)
                                     coef1_imag = sat1%im(i_cf1)
                                     m_sum_2: do i_cf2 = 1, sat2%N_fcts
                                        m2 = sat2%m(i_cf2)
                                        coef2_real = sat2%re(i_cf2)
                                        coef2_imag = sat2%im(i_cf2)
                                        coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag
                                        coef_imag = - coef1_real * coef2_imag + coef1_imag * coef2_real
                                        coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                        coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                                     enddo m_sum_2
                                  enddo m_sum_1
                               enddo spin_up_and_down
                               do m2 = 1,quadrupel%l2*2+1
                                  do m1 = 1,quadrupel%l1*2+1
                                     if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then
                                        saint_2cob_dipole_real(:,:,i_if2,i_if1) = &
                                             saint_2cob_dipole_real(:,:,i_if2,i_if1) + &
                                             coef_real_arr(m2,m1) * cont_int_2cob_dipole(:,:,m2,m1,i_xyz)
                                     endif
                                     if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then
                                        saint_2cob_dipole_imag(:,:,i_if2,i_if1) = &
                                             saint_2cob_dipole_imag(:,:,i_if2,i_if1) + &
                                             coef_imag_arr(m2,m1) * cont_int_2cob_dipole(:,:,m2,m1,i_xyz)
                                     endif
                                  enddo
                               enddo
                            enddo ind_fct_2
                         enddo ind_fct_1
                      endif needed
                    
                      i_ip2 = i_ip2 + 1
                   enddo partner2
                enddo irrep2
                i_ip1 = i_ip1 + 1
             enddo partner1
          enddo irrep1
       enddo xyz
!   endif dipole
    
    call stop_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
  end subroutine symadapt_add_nottotalsym_spor
  !**************************************************************

#ifdef WITH_GTENSOR
  !**************************************************************
  subroutine symadapt_add_nottotalsym_sporgt()
    ! add contribution of one pair
    ! of equal atom to symmetry adapted integrals
    ! the case of not total symmetric integrals
    ! calculates all associated symmetry adapted integrals
    !------------ Declaration of local variables ---------------

    real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
         saint_2cob_dipoleg_real,saint_2cob_dipoleg_imag, &
         saint_2cob_dipoleg_offdiag_real,&
         saint_2cob_dipoleg_offdiag_imag
   
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
         i_ip1, i_if1, i_if2, i_cf1, i_cf2, m1, m2,alpha
    real(kind=r8_kind) :: coef1_real, coef2_real, coef_real, coef1_imag, coef2_imag, coef_imag
    real(kind=r8_kind) :: coef_real_arr(9,9),coef_imag_arr(9,9)

    integer(i4_kind)          :: i
    integer(i4_kind), pointer :: cc_coupling(:) ! Coupling of irreps by complex conjugation

    !------------ Executable code ------------------------------

    if(.not.integralpar_2cob_dipole) RETURN ! FIXME: dont enter at all

#if _DPRINT
print *,'symadapt_add_nottotalsym_sporgt runs...'
print *,'symadapt_add_nottotalsym_sporgt: symadapt_int_2cob_dipoleg_p:'
print *,'symadapt_add_nottotalsym_sporgt  +shape=',shape(symadapt_int_2cob_dipoleg_p)
print *,'symadapt_add_nottotalsym_sporgt  +shape(%int)',&
     & shape(symadapt_int_2cob_dipoleg_p(1,1,1)%int)
#endif


    call start_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
!   gten: if ( integralpar_2cob_dipole ) then
    !!!!???? Diagonal matrix elements are calculated
       xyzd: do i_xyz = 1, 7 ! Lx,Ly,Lz,S,sigma_x,sigma_y,sygma_z + 3 correction terms in the future
          DPRINT 'symadapt_add_nottotalsym_sporgt: xyz=',i_xyz
          i_ip1 = 1
          irrep1_d: do i_ir1 = 1, symmetry_data_n_proj_irreps()
             DPRINT 'symadapt_add_nottotalsym_sporgt: ir1=',i_ir1
             sap1 => ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)
             partner1_d: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                i_ir2 = i_ir1
                i_pa2 = i_pa1
                sap2 => ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)

                needed_d: if (  associated( &
                     symadapt_int_2cob_dipoleg_p(1,i_ip1,i_xyz)%int )) then

                   saint_2cob_dipoleg_real => &
                        symadapt_int_2cob_dipoleg_p(1,i_ip1,i_xyz)%int
                   saint_2cob_dipoleg_imag => &
                        symadapt_int_2cob_dipoleg_p(1,i_ip1,i_xyz)%int_imag
#if _DPRINT
print *,'symadapt_add_nottotalsym_sporgt  +shape(%int)', shape(saint_2cob_dipoleg_real)
print *,'symadapt_add_nottotalsym_sporgt  +shape(%int)', shape(saint_2cob_dipoleg_imag)
#endif
                   select case(i_xyz)
                   case(1:4)! Lx,Ly,Lz,S
                      offind_fct_1_L1_d: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_L1_d: do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind
                            do alpha =1,2
                               sat1 => sap1%sa_spor_int(i_ea1,alpha,i_if1,i_pa1)
                               sat2 => sap2%sa_spor_int(i_ea2,alpha,i_if2,i_pa2)
                               do i_cf1 = 1, sat1%N_fcts
                                  m1 = sat1%m(i_cf1)
                                  coef1_real = sat1%re(i_cf1)
                                  coef1_imag = sat1%im(i_cf1)
                                  do i_cf2 = 1, sat2%N_fcts
                                     m2 = sat2%m(i_cf2)
                                     coef2_real = sat2%re(i_cf2)
                                     coef2_imag = sat2%im(i_cf2)
                                     coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag
                                     coef_imag = - coef1_real * coef2_imag + coef1_imag * coef2_real
                                     coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                     coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                                  enddo
                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_dipoleg_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_real(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)

                                  endif
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) &
                                          - coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)
                                  endif
                               enddo
                            enddo
                         enddo offind_fct_2_L1_d
                      enddo offind_fct_1_L1_d
                   case(7) !Sigmaz
                      offind_fct_1_Sz1_d: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_Sz1_d: do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            spin_up_and_down_sz1_d: do alpha =1,2
                               sat1 => sap1%sa_spor_int(i_ea1,alpha,i_if1,i_pa1)
                               sat2 => sap2%sa_spor_int(i_ea2,alpha,i_if2,i_pa2)
                               do i_cf1 = 1, sat1%N_fcts
                                  m1 = sat1%m(i_cf1)
                                  coef1_real = sat1%re(i_cf1)
                                  coef1_imag = sat1%im(i_cf1)
                                  do i_cf2 = 1, sat2%N_fcts
                                     m2 = sat2%m(i_cf2)
                                     coef2_real = sat2%re(i_cf2)
                                     coef2_imag = sat2%im(i_cf2)
                                 
                                     if(alpha == 2) then
                                        coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag
                                        coef_imag =  -coef1_real * coef2_imag + coef1_imag * coef2_real
                                     else
                                        coef_real = -(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                        coef_imag =  coef1_real * coef2_imag - coef1_imag * coef2_real
                                     end if

                                     coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                     coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                                  enddo
                               enddo
                            enddo spin_up_and_down_sz1_d

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
                                  endif
                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)         
                                  endif

                               enddo
                            enddo
                         end do offind_fct_2_Sz1_d
                      end do offind_fct_1_Sz1_d
                   case (5) !Sigmax
                      offind_fct_1_Sx_d: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_Sx_d: do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind
                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
                                  coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag !ok
                                  coef_imag = (- coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_real = (coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  coef_imag = (-coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
                                  endif

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)

                                  endif
                               enddo
                            enddo
                         end do offind_fct_2_Sx_d
                      end do offind_fct_1_Sx_d
                   case(6) !Sigmay
                      offind_fct_1_Sy_d: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_Sy_d: do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind
                            ! spin_up_and_down_sx: do alpha =1,2
                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =-(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  coef_real =-coef1_real * coef2_imag + coef1_imag * coef2_real  !ok
                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =coef1_real * coef2_real + coef1_imag * coef2_imag
                                  coef_real =-(-coef1_real * coef2_imag + coef1_imag * coef2_real) !DG spicious sign
                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_dipoleg_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
                                  endif
                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
                                  endif

                               enddo
                            enddo
                         end do offind_fct_2_Sy_d
                      end do offind_fct_1_Sy_d
                   end select
                end if needed_d
                i_ip1 = i_ip1 + 1
             end do partner1_d
          end do irrep1_d

       end do xyzd

       ! Overlap between different irreps:
       cc_coupling => symmetry_data_get_cccoupling()
       if( any(cc_coupling.NE.(/(i,i=1,size(cc_coupling))/)) )then
!!$          WARN("FIXME: why to do it at all?")
       i_xyz = 4
       i_ip1 = 1
!!$          irrep4: do i_ir1 = 1, symmetry_data_n_proj_irreps(), 2
       irrep4: do i_ir1 = 1, symmetry_data_n_proj_irreps()
             sap1 => ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)
             partner4: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)

!!$                i_ir2 = i_ir1 + 1 
                i_ir2 = cc_coupling(i_ir1)
                if( i_ir2<i_ir1 )then
                   WARN("FIXME: what about i_ip1?")
                   cycle irrep4 ! already processed
                endif

                i_pa2 = i_pa1 
                sap2 => ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)

                 mainif4: if ( associated( symadapt_int_2cob_dipoleg_p&
                     (2,i_ip1,i_xyz)%int)) then

                   saint_2cob_dipoleg_offdiag_real => &
                        symadapt_int_2cob_dipoleg_p(2,i_ip1,i_xyz)%int
                   saint_2cob_dipoleg_offdiag_imag => &
                        symadapt_int_2cob_dipoleg_p(2,i_ip1,i_xyz)%int_imag

                      offind_fct_1_S1: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_S1: do i_if2 = 1, sap2%N_independent_fcts! 
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
!!$                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)  !OK
!!$                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)  !OK
                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)  !OK
!!$                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)  !OK
                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo
                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) &
                                          +coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)
                                  endif

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)

                                  endif

                               enddo
                            enddo

                         end do offind_fct_2_S1
                      end do offind_fct_1_S1
                   end if mainif4
                   ! i_ip1 = i_ip1 + 1
                    i_ip1 = i_ip1 + 2
                end do partner4
             end do irrep4 
          endif! cc_coupling
!!$!----------------------------------------------------------------------------------------------------------------
       xyz: do i_xyz = 1, 7 ! Lx,Ly,Lz,S,sigma_x,sigma_y,sygma_z  Here offdiagonal matrix elements are calculated
          if(i_xyz == 4) cycle
          i_ip1 = 1
          irrep1: do i_ir1 = 1, symmetry_data_n_proj_irreps()
             sap1 => ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)
             partner1: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                i_ir2 = i_ir1
                i_pa2 = i_pa1
                sap2 => ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)

                 mainif: if ( associated( symadapt_int_2cob_dipoleg_p&
                     (2,i_ip1,i_xyz)%int)) then

                   saint_2cob_dipoleg_offdiag_real => &
                        symadapt_int_2cob_dipoleg_p(2,i_ip1,i_xyz)%int
                   saint_2cob_dipoleg_offdiag_imag => &
                        symadapt_int_2cob_dipoleg_p(2,i_ip1,i_xyz)%int_imag

                   select case(i_xyz)
                       case(4)

                      offind_fct_1_S: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_S: do i_if2 = 1, sap2%N_independent_fcts! 
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
                                  coef_real = (coef1_real * coef2_real + coef1_imag * coef2_imag)  !OK
                                  coef_imag = (coef1_real * coef2_imag - coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag) !OK
                                  coef_imag = (coef1_real * coef2_imag - coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo
                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) &
                                          +coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)
                                  endif

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)

                                  endif

                               enddo
                            enddo

                         end do offind_fct_2_S
                      end do offind_fct_1_S

                   case(1:3)

                      offind_fct_1_L: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_L: do i_if2 = 1, sap2%N_independent_fcts! 
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
                                  coef_imag = (coef1_real * coef2_real - coef1_imag * coef2_imag)  
                                  coef_real = -(coef1_real * coef2_imag + coef1_imag * coef2_real)
!!$
!!$                                  coef_imag = -(coef1_real * coef2_real - coef1_imag * coef2_imag)  
!!$                                  coef_real = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag = -(coef1_real * coef2_real - coef1_imag * coef2_imag) !OK
                                  coef_real = (coef1_real * coef2_imag + coef1_imag * coef2_real)
!!$
!!$                                  coef_imag = +(coef1_real * coef2_real - coef1_imag * coef2_imag) 
!!$                                  coef_real = -(coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo
                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) &
                                          +coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)
                                  endif

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,i_xyz)

                                  endif

                               enddo
                            enddo

                         end do offind_fct_2_L
                      end do offind_fct_1_L

                   case(5) !sigma_x
                      !--------------------------------------------------
                      offind_fct_1_Sxd1: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_Sxd1: do i_if2 = 1, sap2%N_independent_fcts! 
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
!!$                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)

!                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)
!                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)
!Old version
                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)  !OK
                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)
                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)

                       !           coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
                       !           coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)
!Old version
                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) &
                                          +coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4) !!sign!!!
                                  endif

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)

                                  endif

                               enddo
                            enddo

                         end do offind_fct_2_Sxd1
                      end do offind_fct_1_Sxd1

                   case(6) !sigma_y

                      offind_fct_1_Sy: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_Sy: do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_imag =-(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_real =(coef1_real * coef2_imag + coef1_imag * coef2_real)  !ok
                           !       coef_real =-(coef1_real * coef2_imag + coef1_imag * coef2_real) 
                           !       coef_imag =-(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  !Old version
                                  coef_imag =-(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_real =(coef1_real * coef2_imag + coef1_imag * coef2_real)  !ok
                             !     coef_imag =(coef1_real * coef2_real - coef1_imag * coef2_imag)
                             !     coef_real =-(coef1_real * coef2_imag + coef1_imag * coef2_real)  !OK
                                   
                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_imag =-(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_real =(coef1_real * coef2_imag + coef1_imag * coef2_real) 

!                                  coef_real =-(coef1_real * coef2_imag + coef1_imag * coef2_real) 
!                                  coef_imag =-(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  !Old version
                                  coef_imag =-(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_real =(coef1_real * coef2_imag + coef1_imag * coef2_real)

                   !               coef_imag =(coef1_real * coef2_real - coef1_imag * coef2_imag) !OK
                   !               coef_real =-(coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
                                  endif
                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
                                  endif
                               enddo
                            enddo
                         end do offind_fct_2_Sy
                      end do offind_fct_1_Sy

                   case(7) !sigma_z

                      offind_fct_1_Sz: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_Sz: do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2) 
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_real = -2*(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = -2*(coef1_real * coef2_imag + coef1_imag * coef2_real)
!!$                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag) !OK
                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real) 

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo
                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2) 
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_real = -2*(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = -2*(coef1_real * coef2_imag + coef1_imag * coef2_real)
                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real) !OK

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4) !mixing real and
                                  endif!imag part
                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_dipoleg_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)         
                                  endif

                               enddo
                            enddo

                         end do offind_fct_2_Sz
                      end do offind_fct_1_Sz

                   end select
                endif mainif
                i_ip1 = i_ip1 + 1
             enddo partner1
          enddo irrep1
       enddo xyz
!   endif gten

    DPRINT 'symadapt_add_nottotalsym_sporgt...done  ' 
    call stop_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))

  end subroutine symadapt_add_nottotalsym_sporgt!!!!!!! to be changed
  !**************************************************************
  subroutine symadapt_add_nottotalsym_sohfc()
    ! add contribution of one pair
    ! of equal atom to symmetry adapted integrals
    ! the case of not total symmetric integrals
    ! calculates all associated symmetry adapted integrals
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind), pointer, dimension(:,:,:,:) :: saint_2cob_hfc_real,saint_2cob_hfc_imag,&
    & saint_2cob_hfc_offdiag_real,saint_2cob_hfc_offdiag_imag     
    type(unique_atom_partner_type), pointer :: sap1, sap2
    type(unique_atom_sa_int_type), pointer :: sat1, sat2
    integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
         i_ip1, i_if1, i_if2, i_cf1, i_cf2, m1, m2,alpha, i_ua
    real(kind=r8_kind) :: coef1_real, coef2_real, coef_real, coef1_imag, coef2_imag, coef_imag
    real(kind=r8_kind) :: coef_real_arr(9,9),coef_imag_arr(9,9),coef_real_arr_x(9,9),coef_imag_arr_x(9,9),&
         coef_real_arr_y(9,9),coef_imag_arr_y(9,9),coef_real_arr_z(9,9),coef_imag_arr_z(9,9) 

    !------------ Executable code ------------------------------

    if(.not.integralpar_2cob_dipole) RETURN ! FIXME: dont enter at all

    call start_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))
    DPRINT "symadapt_add_nottotalsym_sohfc runs..."

    !--------------------------------------
!   hfc1: if ( integralpar_2cob_dipole ) then
    !!!!???? Diagonal matrix elements are calculated
        ua_loop: Do i_ua = 1 , n_unique_atoms
       xyz_hfc: do i_xyz = 1, 9 ! Lx/r^3,Ly/r^3,Lz/r^3,dipol_sigma_x,dipol_sigma_y,dipol_sygma_z
                                !   isotropic_term: sigmax * mux, sigmay * muy, sugmaz * muz 
          i_ip1 = 1
          irrep1_hfc: do i_ir1 = 1, symmetry_data_n_proj_irreps()
             sap1 => ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)
             partner1_hfc: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                i_ir2 = i_ir1
                i_pa2 = i_pa1
                sap2 => ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)
               
                needed_hfc: if (  associated( &
                     symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz,i_ua)%int )) then

                   saint_2cob_hfc_real => &
                        symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz,i_ua)%int
                   saint_2cob_hfc_imag => &
                        symadapt_int_2cob_hfc_p(1,i_ip1,i_xyz,i_ua)%int_imag
                   select case(i_xyz)
!!$                       case(7)! iso a "unit matrix" variant
!!$                      !10 -> 7 
!!$                          offind_fct_1_diso: do i_if1 = 1, sap1%N_independent_fcts
!!$                             offind_fct_2_diso: do i_if2 = 1, sap2%N_independent_fcts
!!$                            coef_real_arr = 0.0_r8_kind
!!$                            coef_imag_arr = 0.0_r8_kind
!!$                            do alpha =1,2
!!$                               sat1 => sap1%sa_spor_int(i_ea1,alpha,i_if1,i_pa1)
!!$                               sat2 => sap2%sa_spor_int(i_ea2,alpha,i_if2,i_pa2)
!!$                               do i_cf1 = 1, sat1%N_fcts
!!$                                  m1 = sat1%m(i_cf1)
!!$                                  coef1_real = sat1%c(i_cf1)
!!$                                  coef1_imag = sat1%c_imag(i_cf1)
!!$                                  do i_cf2 = 1, sat2%N_fcts
!!$                                     m2 = sat2%m(i_cf2)
!!$                                     coef2_real = sat2%c(i_cf2)
!!$                                     coef2_imag = sat2%c_imag(i_cf2)
!!$                                     coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag
!!$                                     coef_imag = - coef1_real * coef2_imag + coef1_imag * coef2_real
!!$                                     coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
!!$                                     coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
!!$                                  enddo
!!$                               enddo
!!$                            enddo
!!$
!!$                            do m2 = 1,quadrupel%l2*2+1
!!$                               do m1 = 1,quadrupel%l1*2+1
!!$
!!$                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
!!$                                          + coef_real_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)
!!$
!!$                                  endif
!!$                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
!!$                                          + coef_imag_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)
!!$                                  endif
!!$                               enddo
!!$                            enddo
!!$                         enddo offind_fct_2_diso
!!$                      enddo offind_fct_1_diso

                   case(7)!Sigmax * mux
                      do i_if1 = 1, sap1%N_independent_fcts
                          do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind
                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
                                  coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag !ok
                                  coef_imag = (- coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_real = (coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  coef_imag = (-coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
!!$                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$
!!$                                     saint_2cob_dipoleg_real(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_dipoleg_real(:,:,i_if2,i_if1) &
!!$                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
!!$                                  endif
!!$
!!$                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$
!!$                                     saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) &
!!$                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
!!$
!!$                                  endif
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)

                                  endif
                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)
                                  endif
                               enddo
                            enddo
                         end do
                      end do
                   case (8) ! Sigmay * muy
                       do i_if1 = 1, sap1%N_independent_fcts
                          do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind
                            ! spin_up_and_down_sx: do alpha =1,2
                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =-(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  coef_real =-coef1_real * coef2_imag + coef1_imag * coef2_real  !ok
                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =coef1_real * coef2_real + coef1_imag * coef2_imag
                                  coef_real =-(-coef1_real * coef2_imag + coef1_imag * coef2_real) !DG spicious sign
                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
!!$                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$                                     saint_2cob_dipoleg_real(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_dipoleg_real(:,:,i_if2,i_if1) &
!!$                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
!!$                                  endif
!!$                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$                                     saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) &
!!$                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
!!$                                  endif
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)

                                  endif
                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)
                                  endif

                               enddo
                            enddo
                         end do
                      end do

                   case (9) !9 Sigmaz * muz
                     do i_if1 = 1, sap1%N_independent_fcts
                          do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                             do alpha =1,2
                               sat1 => sap1%sa_spor_int(i_ea1,alpha,i_if1,i_pa1)
                               sat2 => sap2%sa_spor_int(i_ea2,alpha,i_if2,i_pa2)
                               do i_cf1 = 1, sat1%N_fcts
                                  m1 = sat1%m(i_cf1)
                                  coef1_real = sat1%re(i_cf1)
                                  coef1_imag = sat1%im(i_cf1)
                                  do i_cf2 = 1, sat2%N_fcts
                                     m2 = sat2%m(i_cf2)
                                     coef2_real = sat2%re(i_cf2)
                                     coef2_imag = sat2%im(i_cf2)
                                 
                                     if(alpha == 2) then
                                        coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag
                                        coef_imag =  -coef1_real * coef2_imag + coef1_imag * coef2_real
                                     else
                                        coef_real = -(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                        coef_imag =  coef1_real * coef2_imag - coef1_imag * coef2_real
                                     end if

                                     coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                     coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                                  enddo
                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
!!$                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$
!!$                                     saint_2cob_dipoleg_real(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_dipoleg_real(:,:,i_if2,i_if1) &
!!$                                          + coef_real_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)
!!$                                  endif
!!$                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$
!!$                                     saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_dipoleg_imag(:,:,i_if2,i_if1) &
!!$                                          + coef_imag_arr(m2,m1) * cont_int_2cob_dipoleg(:,:,m2,m1,4)         
!!$                                  endif
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)

                                  endif
                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,10,i_ua)
                                  endif

                               enddo
                            enddo
                         end do
                      end do
                      
                   case(1:3)! Lx,Ly,Lz
                      ! Mapping to symadapt integrals
                      ! 1 =  Lx/r^3
                      ! 2 =  Ly/r^3
                      ! 3 =  Lz/r^3
                      ! Mapping to contracted primitive integrals
                      ! 7 =  Lx/r^3
                      ! 8 =  Ly/r^3
                      ! 9 =  Lz/r^3
                       do i_if1 = 1, sap1%N_independent_fcts
                          do i_if2 = 1, sap2%N_independent_fcts
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind
                            do alpha =1,2
                               sat1 => sap1%sa_spor_int(i_ea1,alpha,i_if1,i_pa1)
                               sat2 => sap2%sa_spor_int(i_ea2,alpha,i_if2,i_pa2)
                               do i_cf1 = 1, sat1%N_fcts
                                  m1 = sat1%m(i_cf1)
                                  coef1_real = sat1%re(i_cf1)
                                  coef1_imag = sat1%im(i_cf1)
                                  do i_cf2 = 1, sat2%N_fcts
                                     m2 = sat2%m(i_cf2)
                                     coef2_real = sat2%re(i_cf2)
                                     coef2_imag = sat2%im(i_cf2)
                                     coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag
                                     coef_imag = - coef1_real * coef2_imag + coef1_imag * coef2_real
                                     coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                     coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
                                  enddo
                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,i_xyz+6,i_ua)
                                     !1,2,3 ->7,8,9

                                  endif
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          - coef_real_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,i_xyz+6,i_ua)
                                  endif
                               enddo
                            enddo
                         enddo
                      enddo
                   case(4:6) !Sigmaz,Sigmay,Sigmax
                      ! Mapping to symadapt integrals
                      ! 4 = 3*sigma_x * (x^2-r^2)/r^5 + 3 * sigma_y * xy/r^5 + 3*sigma_z * xz/r^5
                      ! 5 = 3*sigma_y * (y^2-r^5)/r^5 + 3 * sigma_x * xy/r^5 + 3*sigma_z * zy/r^5
                      ! 6 = 3*sigma_z * (z^2-r^2)/r^5 + 3 * sigma_x * xz/r^5 + 3*sigma_y * yz/r^5
                      ! Mapping to contracted primitive integrals
                      ! 1 = d2/dx2   4 = d2/dxdy
                      ! 2 = d2/dy2   5 = d2/dxdz
                      ! 3 = d2/dz2   6 = d2/dydz
                      offind_fct_1_Sz1_d: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_Sz1_d: do i_if2 = 1, sap2%N_independent_fcts
                            ! contribution from sigma_z
                            coef_real_arr_z = 0.0_r8_kind
                            coef_imag_arr_z = 0.0_r8_kind
!!$                            coef_real_arr_x = 0.0_r8_kind
!!$                            coef_imag_arr_x = 0.0_r8_kind
!!$                            coef_real_arr_y = 0.0_r8_kind
!!$                            coef_imag_arr_y = 0.0_r8_kind
                            spin_up_and_down_sz1_d: do alpha =1,2
                               sat1 => sap1%sa_spor_int(i_ea1,alpha,i_if1,i_pa1)
                               sat2 => sap2%sa_spor_int(i_ea2,alpha,i_if2,i_pa2)
                               do i_cf1 = 1, sat1%N_fcts
                                  m1 = sat1%m(i_cf1)
                                  coef1_real = sat1%re(i_cf1)
                                  coef1_imag = sat1%im(i_cf1)
                                  do i_cf2 = 1, sat2%N_fcts
                                     m2 = sat2%m(i_cf2)
                                     coef2_real = sat2%re(i_cf2)
                                     coef2_imag = sat2%im(i_cf2)
                                     
!                                     if(alpha == 1) then
                                      if(alpha == 2) then
                                        coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag
                                        coef_imag =  -coef1_real * coef2_imag + coef1_imag * coef2_real
                                     else
                                        coef_real = -(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                        coef_imag =  coef1_real * coef2_imag - coef1_imag * coef2_real
                                     end if

                                     coef_real_arr_z(m2,m1) = coef_real_arr_z(m2,m1) + coef_real
                                     coef_imag_arr_z(m2,m1) = coef_imag_arr_z(m2,m1) + coef_imag
!!$                                     !contribution from sigma_x
!!$                                     coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag !ok
!!$                                     coef_imag = (- coef1_real * coef2_imag + coef1_imag * coef2_real)
!!$
!!$                                     coef_real_arr_x(m2,m1) = coef_real_arr_x(m2,m1) + coef_real
!!$                                     coef_imag_arr_x(m2,m1) = coef_imag_arr_x(m2,m1) + coef_imag
!!$                                     !contribution from sigma_y
!!$                                     if(alpha == 1) then
!!$                                     coef_imag =-(coef1_real * coef2_real + coef1_imag * coef2_imag)
!!$                                     coef_real =-coef1_real * coef2_imag + coef1_imag * coef2_real  !ok
!!$                                     else
!!$                                        coef_imag =(coef1_real * coef2_real + coef1_imag * coef2_imag)
!!$                                        coef_real =-(-coef1_real * coef2_imag + coef1_imag * coef2_real)  
!!$                                     end if
!!$                                     coef_real_arr_y(m2,m1) = coef_real_arr_y(m2,m1) + coef_real
!!$                                     coef_imag_arr_y(m2,m1) = coef_imag_arr_y(m2,m1) + coef_imag
                                  enddo
                               enddo
                            enddo spin_up_and_down_sz1_d

                            !Sigmax coeffitients
                      
                            coef_real_arr_x = 0.0_r8_kind
                            coef_imag_arr_x = 0.0_r8_kind
                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
                                  coef_real = coef1_real * coef2_real + coef1_imag * coef2_imag !ok
                                  coef_imag = (- coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr_x(m2,m1) = coef_real_arr_x(m2,m1) + coef_real
                                  coef_imag_arr_x(m2,m1) = coef_imag_arr_x(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_real = (coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  coef_imag = (-coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr_x(m2,m1) = coef_real_arr_x(m2,m1) + coef_real
                                  coef_imag_arr_x(m2,m1) = coef_imag_arr_x(m2,m1) + coef_imag
                               enddo
                            enddo

                            !Sigmay coeffitients
                            coef_real_arr_y = 0.0_r8_kind
                            coef_imag_arr_y = 0.0_r8_kind
                           
                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =-(coef1_real * coef2_real + coef1_imag * coef2_imag)
                                  coef_real =-coef1_real * coef2_imag + coef1_imag * coef2_real  !ok
                                  coef_real_arr_y(m2,m1) = coef_real_arr_y(m2,m1) + coef_real
                                  coef_imag_arr_y(m2,m1) = coef_imag_arr_y(m2,m1) + coef_imag
                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =coef1_real * coef2_real + coef1_imag * coef2_imag
                                  coef_real =-(-coef1_real * coef2_imag + coef1_imag * coef2_real) !DG spicious sign
                                  coef_real_arr_y(m2,m1) = coef_real_arr_y(m2,m1) + coef_real
                                  coef_imag_arr_y(m2,m1) = coef_imag_arr_y(m2,m1) + coef_imag
                               enddo
                            enddo


                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  ! Mapping to symadapt integrals
                                  ! 4 = 3*sigma_x * (x^2-r^2)/r^5 + 3 * sigma_y * xy/r^5 + 3*sigma_z * xz/r^5
                                  ! 5 = 3*sigma_y * (y^2-r^5)/r^5 + 3 * sigma_x * xy/r^5 + 3*sigma_z * zy/r^5
                                  ! 6 = 3*sigma_z * (z^2-r^2)/r^5 + 3 * sigma_x * xz/r^5 + 3*sigma_y * yz/r^5
                                  ! Mapping to contracted primitive integrals
                                  ! 1 = d2/dx2   4 = d2/dxdy
                                  ! 2 = d2/dy2   5 = d2/dxdz
                                  ! 3 = d2/dz2   6 = d2/dydz
                                xyzeq6:  if (i_xyz == 6 ) then ! Term  6 
                                  !contributions from sigma z
                                  if (abs(coef_real_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,3,i_ua) !(z^2-r^2)/r^5
                                  endif
                                  if (abs(coef_imag_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,3,i_ua) !  
                                  endif
                                  !contributions from sigma y
                                  if (abs(coef_real_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !5 yz/r^5
                                  endif
                                  if (abs(coef_imag_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !5 
                                  endif
                                  !contributions from sigma x
                                  if (abs(coef_real_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) !xz/r^5 6
                                  endif
                                  if (abs(coef_imag_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) ! 6
                                  endif
                               end if xyzeq6

                               xyzeq5:  if (i_xyz == 5 ) then !term 5
                                  !contributions from sigma z
                                  if (abs(coef_real_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !zy/r^5
                                  endif
                                  if (abs(coef_imag_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !       
                                  endif
                                  !contributions from sigma y
                                  if (abs(coef_real_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,2,i_ua) !(y^2-r^2)/r^5
                                  endif
                                  if (abs(coef_imag_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,2,i_ua) !
                                  endif
                                  !contributions from sigma x
                                  if (abs(coef_real_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !xy/r^5
                                  endif
                                  if (abs(coef_imag_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !
                                  endif
                               end if xyzeq5
                               
                                xyzeq4:  if (i_xyz == 4 ) then !term 4
                                  !contributions from sigma z
                                  if (abs(coef_real_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) !xz/r^5
                                  endif
                                  if (abs(coef_imag_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) !       
                                  endif
                                  !contributions from sigma y
                                  if (abs(coef_real_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !xy/r^5
                                  endif
                                  if (abs(coef_imag_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !
                                  endif
                                  !contributions from sigma x
                                  if (abs(coef_real_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,1,i_ua) !(x^2-r^2)/r^5
                                  endif
                                  if (abs(coef_imag_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,1,i_ua) !
                                  endif
                               end if xyzeq4
                               
                               enddo
                            enddo
                         end do offind_fct_2_Sz1_d
                      end do offind_fct_1_Sz1_d   
                   end select
                end if needed_hfc
            
                i_ip1 = i_ip1 + 1
             end do partner1_hfc
          end do irrep1_hfc

       end do xyz_hfc
    end Do ua_loop

    ualoop: do i_ua = 1,n_unique_atoms
       xyz_off: do i_xyz = 1, 7 ! Lx,Ly,Lz,sigma_x,sigma_y,sygma_z  Here offdiagonal matrix elements are calculated

          i_ip1 = 1
          irrep1_off: do i_ir1 = 1, symmetry_data_n_proj_irreps()
             sap1 => ua1%symadapt_spor_partner(i_ir1,quadrupel%l1)
             partner1_off: do i_pa1 = 1, symmetry_data_n_partners_proj(i_ir1)
                i_ir2 = i_ir1
                i_pa2 = i_pa1
                sap2 => ua2%symadapt_spor_partner(i_ir2,quadrupel%l2)
              
                 mainif_off: if ( associated( symadapt_int_2cob_hfc_p&
                     (2,i_ip1,i_xyz,i_ua)%int)) then

                   saint_2cob_hfc_offdiag_real => &
                        symadapt_int_2cob_hfc_p(2,i_ip1,i_xyz,i_ua)%int
                   saint_2cob_hfc_offdiag_imag => &
                        symadapt_int_2cob_hfc_p(2,i_ip1,i_xyz,i_ua)%int_imag

                   select case(i_xyz)
!!$                    !We do not need it  case (7)! iso
!!$                           offind_fct_1_iso: do i_if1 = 1, sap1%N_independent_fcts
!!$                             offind_fct_2_iso: do i_if2 = 1, sap2%N_independent_fcts! 
!!$                            coef_real_arr = 0.0_r8_kind
!!$                            coef_imag_arr = 0.0_r8_kind
!!$
!!$                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
!!$                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
!!$
!!$                            do i_cf1 = 1, sat1%N_fcts
!!$                               m1 = sat1%m(i_cf1)
!!$                               coef1_real = sat1%c(i_cf1)
!!$                               coef1_imag = sat1%c_imag(i_cf1)
!!$                               do i_cf2 = 1, sat2%N_fcts
!!$                                  m2 = sat2%m(i_cf2)
!!$                                  coef2_real = sat2%c(i_cf2)
!!$                                  coef2_imag = sat2%c_imag(i_cf2)
!!$                             !     coef_imag = (coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                             !     coef_real = -(coef1_real * coef2_imag + coef1_imag * coef2_real)
!!$                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)
!!$
!!$                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
!!$                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
!!$
!!$                               enddo
!!$                            enddo
!!$
!!$                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
!!$                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
!!$
!!$                            do i_cf1 = 1, sat1%N_fcts
!!$                               m1 = sat1%m(i_cf1)
!!$                               coef1_real = sat1%c(i_cf1)
!!$                               coef1_imag = sat1%c_imag(i_cf1)
!!$                               do i_cf2 = 1, sat2%N_fcts
!!$                                  m2 = sat2%m(i_cf2)
!!$                                  coef2_real = sat2%c(i_cf2)
!!$                                  coef2_imag = sat2%c_imag(i_cf2)
!!$
!!$                                  coef_imag = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_real = (coef1_real * coef2_imag + coef1_imag * coef2_real)
!!$
!!$                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
!!$                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag
!!$
!!$                               enddo
!!$                            enddo
!!$                            do m2 = 1,quadrupel%l2*2+1
!!$                               do m1 = 1,quadrupel%l1*2+1
!!$                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$
!!$                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
!!$                                          +coef_real_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,i_xyz,i_ua)
!!$                                  endif
!!$
!!$                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 
!!$
!!$                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
!!$                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
!!$                                          + coef_imag_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,i_xyz,i_ua)
!!$
!!$                                  endif
!!$
!!$                               enddo
!!$                            enddo
!!$
!!$                         end do offind_fct_2_iso
!!$                      end do offind_fct_1_iso

                   case(1:3)
                      ! Mapping to symadapt integrals
                      ! 1 =  Lx/r^3
                      ! 2 =  Ly/r^3
                      ! 3 =  Lz/r^3
                      ! Mapping to contracted primitive integrals
                      ! 7 =  Lx/r^3
                      ! 8 =  Ly/r^3
                      ! 9 =  Lz/r^3
                      offind_fct_1_L: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_L: do i_if2 = 1, sap2%N_independent_fcts! 
                            coef_real_arr = 0.0_r8_kind
                            coef_imag_arr = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
                                  coef_imag = (coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_real = -(coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_real = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr(m2,m1) = coef_real_arr(m2,m1) + coef_real
                                  coef_imag_arr(m2,m1) = coef_imag_arr(m2,m1) + coef_imag

                               enddo
                            enddo
                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                  if (abs(coef_real_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          +coef_real_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,i_xyz+6,i_ua)
                                  endif

                                  if (abs(coef_imag_arr(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,i_xyz+6,i_ua)

                                  endif

                               enddo
                            enddo

                         end do offind_fct_2_L
                      end do offind_fct_1_L
                      
                        case(4:6) !Sigmaz,Sigmay,Sigmax
                      ! Mapping to symadapt integrals
                      ! 4 = 3*sigma_x * d2/dx2 + 3 * sigma_y * d2/dydz + 3*sigma_z * d2/dzdx
                      ! 5 = 3*sigma_x * d2/dy2 + 3 * sigma_x * d2/dydx + 3*sigma_z * d2/dzdy
                      ! 6 = 3*sigma_x * d2/dz2 + 3 * sigma_y * d2/dydz + 3*sigma_x * d2/dzdx
                      ! Mapping to contracted primitive integrals
                      ! 1 = d2/dx2   4 = d2/dxdy
                      ! 2 = d2/dy2   5 = d2/dxdz
                      ! 3 = d2/dz2   6 = d2/dydz
                      offind_fct_1_S: do i_if1 = 1, sap1%N_independent_fcts
                         offind_fct_2_S: do i_if2 = 1, sap2%N_independent_fcts
                          !  Sigmaz
                            coef_real_arr_z = 0.0_r8_kind
                            coef_imag_arr_z = 0.0_r8_kind
!!$                            coef_real_arr_x = 0.0_r8_kind
!!$                            coef_imag_arr_x = 0.0_r8_kind
!!$                            coef_real_arr_y = 0.0_r8_kind
!!$                            coef_imag_arr_y = 0.0_r8_kind
                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2) 
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)
                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr_z(m2,m1) = coef_real_arr_z(m2,m1) + coef_real
                                  coef_imag_arr_z(m2,m1) = coef_imag_arr_z(m2,m1) + coef_imag

                               enddo
                            enddo
                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2) 
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

!!$                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
!!$                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)
                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr_z(m2,m1) = coef_real_arr_z(m2,m1) + coef_real
                                  coef_imag_arr_z(m2,m1) = coef_imag_arr_z(m2,m1) + coef_imag

                               enddo
                            enddo

                               !sigmax 
                               coef_real_arr_x = 0.0_r8_kind
                               coef_imag_arr_x = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)
                                  coef_real = (coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_imag = (coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr_x(m2,m1) = coef_real_arr_x(m2,m1) + coef_real
                                  coef_imag_arr_x(m2,m1) = coef_imag_arr_x(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)

                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_real = -(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_imag = -(coef1_real * coef2_imag + coef1_imag * coef2_real)

                                  coef_real_arr_x(m2,m1) = coef_real_arr_x(m2,m1) + coef_real
                                  coef_imag_arr_x(m2,m1) = coef_imag_arr_x(m2,m1) + coef_imag

                               enddo
                            enddo
                            !sigmay
                            coef_real_arr_y = 0.0_r8_kind
                            coef_imag_arr_y = 0.0_r8_kind

                            sat1 => sap1%sa_spor_int(i_ea1,1,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,1,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =-(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_real =(coef1_real * coef2_imag + coef1_imag * coef2_real)  !ok

                                  coef_real_arr_y(m2,m1) = coef_real_arr_y(m2,m1) + coef_real
                                  coef_imag_arr_y(m2,m1) = coef_imag_arr_y(m2,m1) + coef_imag

                               enddo
                            enddo

                            sat1 => sap1%sa_spor_int(i_ea1,2,i_if1,i_pa1)
                            sat2 => sap2%sa_spor_int(i_ea2,2,i_if2,i_pa2)
                            do i_cf1 = 1, sat1%N_fcts
                               m1 = sat1%m(i_cf1)
                               coef1_real = sat1%re(i_cf1)
                               coef1_imag = sat1%im(i_cf1)
                               do i_cf2 = 1, sat2%N_fcts
                                  m2 = sat2%m(i_cf2)
                                  coef2_real = sat2%re(i_cf2)
                                  coef2_imag = sat2%im(i_cf2)

                                  coef_imag =-(coef1_real * coef2_real - coef1_imag * coef2_imag)
                                  coef_real =(coef1_real * coef2_imag + coef1_imag * coef2_real) 
                                  coef_real_arr_y(m2,m1) = coef_real_arr_y(m2,m1) + coef_real
                                  coef_imag_arr_y(m2,m1) = coef_imag_arr_y(m2,m1) + coef_imag

                               enddo
                            enddo

                            do m2 = 1,quadrupel%l2*2+1
                               do m1 = 1,quadrupel%l1*2+1
                                xyzdeq6:  if (i_xyz == 6 ) then !Term 6
                                  !contributions from sigma z 
                                  if (abs(coef_real_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,3,i_ua) !!(z^2-r^2)/r^5
                                  endif
                                  if (abs(coef_imag_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,3,i_ua) ! 
                                  endif
                                  !contributions from sigma y
                                  if (abs(coef_real_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !yz/r^5
                                  endif
                                  if (abs(coef_imag_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !
                                  endif
                                  !contributions from sigma x
                                  if (abs(coef_real_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) !xz/r^5
                                  endif
                                  if (abs(coef_imag_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) !
                                  endif
                               end if xyzdeq6

                               xyzdq5:  if (i_xyz == 5 ) then
                                  !contributions from sigma z
                                  if (abs(coef_real_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !d2/dzdy
                                  endif
                                  if (abs(coef_imag_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,6,i_ua) !d2/dzdy       
                                  endif
                                  !contributions from sigma y
                                  if (abs(coef_real_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,2,i_ua) !d2/dy2
                                  endif
                                  if (abs(coef_imag_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,2,i_ua) !d2/dy2
                                  endif
                                  !contributions from sigma x
                                  if (abs(coef_real_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !d2/dxdy
                                  endif
                                  if (abs(coef_imag_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !d2/dxdy
                                  endif
                               end if xyzdq5
                               
                                xyzdeq4:  if (i_xyz == 4 ) then
                                  !contributions from sigma z
                                  if (abs(coef_real_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) !d2/dzdx
                                  endif
                                  if (abs(coef_imag_arr_z(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_z(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,5,i_ua) !d2/dzdx       
                                  endif
                                  !contributions from sigma y
                                  if (abs(coef_real_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !d2/dydx
                                  endif
                                  if (abs(coef_imag_arr_y(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_y(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,4,i_ua) !d2/dydx
                                  endif
                                  !contributions from sigma x

                                  if (abs(coef_real_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_real(:,:,i_if2,i_if1) &
                                          + coef_real_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,1,i_ua) !d2/dx2
                                  endif
                                  if (abs(coef_imag_arr_x(m2,m1)).gt.0.0_r8_kind) then 

                                     saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) = &
                                          saint_2cob_hfc_offdiag_imag(:,:,i_if2,i_if1) &
                                          + coef_imag_arr_x(m2,m1) * cont_int_2cob_hfc(:,:,m2,m1,1,i_ua) !d2/dx2
                                  endif
                               end if xyzdeq4
                               
                            enddo
                         enddo
                      end do offind_fct_2_S
                   end do offind_fct_1_S
                end select
             endif mainif_off
            
                i_ip1 = i_ip1 + 1
             enddo partner1_off
          enddo irrep1_off
       enddo xyz_off
    end do ualoop
!   endif hfc1
 
    DPRINT "symadapt_add_nottotalsym_sohfc...done  "
    call stop_timer(timer_int_symadapt_2cob3c(integralpar_i_int_part))

  end subroutine symadapt_add_nottotalsym_sohfc!!!!!!! to be changed
#endif

#if 0
  !**************************************************************
  subroutine symadapt_write()
    ! output to debug unit
    !------------ Declaration of local variables ---------------
    real(kind=r8_kind), pointer, dimension(:,:,:,:) :: &
         saint_2cob_dipole
    type(unique_atom_partner_type), pointer :: sap1, sap2
    integer(kind=i4_kind) :: i_xyz, i_ir1, i_ir2, i_pa1, i_pa2, &
         i_ip1, i_ip2, i_if1, i_if2, i_c1, i_c2
    !------------ Executable code ------------------------------
    dipole: if ( integralpar_2cob_dipole ) then
       write(debug_unit,*)
       write(debug_unit,*)
       write(debug_unit,*) "######### symmetry-adapted dipole integrals ##########"
       write(debug_unit,*) "i_c1,  i_c2,  dipole"
       xyz: do i_xyz = 1, 3
          write(debug_unit,*) "######## i_xyz ", i_xyz
          i_ip1 = 1
          irrep1: do i_ir1 = 1, symmetry_data_n_irreps()
             sap1 => ua1%symadapt_partner(i_ir1,quadrupel%l1)
             partner1: do i_pa1 = 1, symmetry_data_n_partners(i_ir1)
                i_ip2 = 1
                irrep2: do i_ir2 = 1, i_ir1
                   sap2 => ua2%symadapt_partner(i_ir2,quadrupel%l2)
                   partner2: do i_pa2 = 1, symmetry_data_n_partners(i_ir2)
                      needed: if ( associated( &
                           symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int ) ) then
                         saint_2cob_dipole => symadapt_int_2cob_dipole(i_ip2,i_ip1,i_xyz)%int
                         ind_fct_1: do i_if1 = 1, sap1%N_independent_fcts
                            ind_fct_2: do i_if2 = 1, sap2%N_independent_fcts
                               write(debug_unit,*) "##### Irrep1 ", i_ir1, " Irrep2 ", i_ir2, &
                                    " Partner1 ", i_pa1, " Partner2 ", i_pa2, &
                                    " ind_fct_1 ", i_if1, " ind_fct_2 ", i_if2
                               do i_c1 = 1, ubound(saint_2cob_dipole,2)
                                  do i_c2 = 1, ubound(saint_2cob_dipole,1)
                                     write(debug_unit,'(2I4,3F20.8)') i_c1, i_c2, &
                                          saint_2cob_dipole(i_c2,i_c1,i_if2,i_if1)
                                  enddo
                               enddo
                            enddo ind_fct_2
                         enddo ind_fct_1
                      endif needed
                      i_ip2 = i_ip2 + 1
                   enddo partner2
                enddo irrep2
                i_ip1 = i_ip1 + 1
             enddo partner1
          enddo irrep1
       enddo xyz
    endif dipole
    write(debug_unit,*) "######### symmetry-adapted integrals end ##########"
    write(debug_unit,*)
    write(debug_unit,*)
  end subroutine symadapt_write
  !**************************************************************
#endif

end subroutine integral_calc_quad_dipole
