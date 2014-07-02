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
!=====================================================================
! Public interface of module
!=====================================================================
module resp_util_module
  !-------------------------------------------------------------------
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
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
  use iounitadmin_module
  use symmetry_data_module
  use eigen_data_module, only : eigvec, eigval
  use occupation_module, only: occ_num,n_occo
  use phys_param_module
  use comm_module
  use xpack
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of constants and variables ---------------
  ! integer flag array indicating how the MO "i_mo" shall be used in the
  ! calculation of 3-index-integrals, where "i_mo=1,...,dim(i_irrep)".
  ! "MO_status(i_mo, i_irrep, i_spin)":
  !  = -1   :  do not use this MO for the 3-index-integrals 
  !            (useful for neglecting unphysically high MO`s)
  !  =  0   :  use as an fully occupied MO 
  !  =  1   :  use as an partially occupied MO
  !  =  2   :  use as an MO which contains absolutely NO CHARGE (=:empty)
  ! Note: this variable is calculated on the MASTER and then distributed
  !       in subroutine "send_level_index"
  ! further variables:
  ! begin_index(i,irrep,spin), i=1,2,3
  !    i=1: first index in MO_status() with value 0
  !    i=2: first index in MO_status() with value 1
  !    i=3: first index in MO_status() with value 2
  ! end_index(i,irrep,spin), i=1,2,3: dito for last index
  ! n_full    - full orbitals
  ! n_partial - partially filled orbitals
  ! n_empty   - empty orbitals
  ! max_level_energy_au : input parameter converted from eV to hartrees
  integer(kind=i4_kind),allocatable,public :: &
       & MO_status  (:,:,:),&
       & begin_index(:,:,:),&
       & end_index  (:,:,:),&
       & n_full     (:,:),  &
       & n_partial  (:,:),  &
       & n_empty    (:,:)

  ! Difference between levels: if gap between two levels smaller than min_diff,
  ! such transitions will be ignored.
  real(kind=r8_kind), parameter, public :: min_diff = 0.00000001_r8_kind 

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public resp_util_set_level_index, &
       resp_util_upck_level_index, &
       resp_util_set_trans_bounds, &
       resp_util_calc_ou, &
       resp_util_calc_transitions_v2, &
       resp_util_fname, &
       resp_util_borders, &
       resp_util_buildi, &
       resp_util_buildr

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************
  subroutine resp_util_set_level_index(&
       & unoccupied_level_criterion,&
       & limit_unoccupied_levels,&
       & max_level_index,&
       & min_level_energy_au,&
       & max_level_energy_au, &
       & min_unocc_level_energy_au,&
       & max_unocc_level_energy_au, &
       & num_spectrum_levels)
    !  Purpose: 
    !  Determine the value of the integer flag array "MO_status"
    !  using the user input parameters
    !    "max_empty_level_index" 
    !    "max_empty_level_energy"
    !    "unoccupied_level_criterion"
    !  and distribute it to the slaves.
    !  The input parameters "max_empty_level_index" and
    !  "max_empty_level_energy" are mutually exclusive: if the
    !  parameter "max_empty_level_index" is specified,
    !  then the parameter "max_empty_level_energy" is ignored.
    !  The following rules are used to determine the
    !  possible values of "MO_status(i_mo, i_irrep, i_spin)",
    !  where "i_mo=1,...,dim(i_irrep)":
    !  = -1   :  do not use this MO for the 3-index-integrals 
    !            This value is assigned to all MO`s with an index
    !            larger than "max_empty_level_index", or, if the 
    !            former is unspecified, with an energy eigenvalue
    !            larger than "max_empty_level_energy".
    !  =  0   :  use as an fully occupied MO 
    !            All MO`s with a charge larger than
    !               N_spin - "unoccupied_level_criterion"
    !  =  1   :  use as an partially occupied MO
    !            All MO`s which do not fall into any of the other 
    !            categories.
    !  =  2   :  use as an MO which contains absolutely NO 
    !            CHARGE (=:empty).
    !            All MO`s with a charge less or equal than
    !               "unoccupied_level_criterion"
    !
    !  NOTE: This subroutine is called by MASTER        
    !------------ Modules used ------------------- ---------------
    use msgtag_module, only: msgtag_response_level_send
    implicit none
    real(kind=r8_kind),intent(in) :: &
         & unoccupied_level_criterion, &
         & min_level_energy_au, &
         & max_level_energy_au, &
         & min_unocc_level_energy_au, &
         & max_unocc_level_energy_au
    logical,intent(in)            :: limit_unoccupied_levels
    integer(kind=i4_kind),intent(in) :: & 
         & max_level_index, &
         & num_spectrum_levels
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    logical               :: was_found(3)
    integer(kind=i4_kind) :: alloc_stat, i_irrep, i_spin, i_dim,&
         & print_dim
    integer(kind=i4_kind) :: max_dim, new_max_index
    integer(kind=i4_kind) :: i_empty, i_rest, i_full, i_partial
    real(kind=r8_kind)    :: full_charge, empty_charge, new_max_energy, new_min_energy
    real(kind=r8_kind)    :: new_unocc_max_energy, new_unocc_min_energy
    ! symmetry information

    integer(kind=i4_kind) :: n_spin, n_irrep 
    !------------ Executable code ------------------------------------

    !## then allocate memory for the array MO_status

    ! find largest dimension of all irreps
    max_dim = maxval(ssym%dim(:))
    n_spin  = ssym%n_spin      ! number of spins
    n_irrep = ssym%n_irrep     ! number of irreps

    allocate(MO_status(max_dim,n_irrep,n_spin),&
         & begin_index(3,n_irrep,n_spin),&
         & end_index(3,n_irrep,n_spin),&
         & n_full(n_irrep,n_spin),n_partial(n_irrep,n_spin),&
         & n_empty(n_irrep,n_spin),stat=alloc_stat)
    ASSERT(alloc_stat==0)

    !## find status of all MO`s according to rules documented above

    ! this will be the minimum charge for a fully occupied level
    !  it should be a little less then 2.0 for spin=1
    !           and a little less then 1.0 for spin=2 
    full_charge = 2.0_r8_kind/real(n_spin,kind=r8_kind)&
         & - unoccupied_level_criterion
    ! this will be the maximum charge allowed in an totally empty level
    empty_charge = unoccupied_level_criterion

    ! initialize MO_status with -1 
    ! note: size(MO_status(:,i_irrep, i_spin)) may be larger than
    !       ssym%dim(i_irrep) !!!
    MO_status = -1_i4_kind

    ! Note: "occ_num(i_irrep)%m(i_dim,i_spin)" is ... ???
    do i_irrep=1, n_irrep
       do i_spin=1, n_spin

          ! consistency check of limit levels criteria
          if(limit_unoccupied_levels) then  ! such a limit requested ?
             if(max_level_index/=-1_i4_kind) then
                if(max_level_index> ssym%dim(i_irrep) ) then
                   call write_to_output_units("response_send_level_index: *** warning ***: &
                        & max_level_index > dimension of irrep ",inte=i_irrep)
                   new_max_index = ssym%dim(i_irrep)
                else
                   new_max_index = max_level_index
                end if
             else
                !! cutoff for occupied
                new_min_energy = max(min_level_energy_au,minval(eigval(i_irrep)%m(:, i_spin)))
                new_max_energy = min(max_level_energy_au,maxval(eigval(i_irrep)%m(:, i_spin)))

                !! cutoff for unoccupied
                new_unocc_min_energy = max(min_unocc_level_energy_au,&
                     minval(eigval(i_irrep)%m(:, i_spin)))
                new_unocc_max_energy = min(max_unocc_level_energy_au,&
                     maxval(eigval(i_irrep)%m(:, i_spin)))
             end if
          end if

          do i_dim=1, ssym%dim(i_irrep)

             MO_status(i_dim, i_irrep,i_spin) = -1

             if(limit_unoccupied_levels) then
                if(max_level_index /= -1) then
                   ! check if the index of this level is to high
                   if( i_dim < new_max_index) &
                        MO_status(i_dim,i_irrep,i_spin) = stOfMO(i_irrep,i_spin,i_dim)
                else
                   ! occupied
                   if(  (eigval(i_irrep)%m(i_dim, i_spin) >= new_min_energy) .and. &
                        (eigval(i_irrep)%m(i_dim, i_spin) <= new_max_energy)) then
                      MO_status(i_dim,i_irrep,i_spin) = stOfMO(i_irrep,i_spin,i_dim)
                   end if
                   ! unoccupied
                   if(  (eigval(i_irrep)%m(i_dim, i_spin) >= new_unocc_min_energy) .and. &
                        (eigval(i_irrep)%m(i_dim, i_spin) <= new_unocc_max_energy)) then
                      MO_status(i_dim,i_irrep,i_spin) = stOfMO(i_irrep,i_spin,i_dim)
                   end if
                end if
             else
                MO_status(i_dim, i_irrep,i_spin) = stOfMO(i_irrep,i_spin,i_dim)
             end if

          end do
       end do
    end do

    ! ## determine first and last index of array elements with values 0,1,2
    do i_irrep=1, n_irrep
       do i_spin=1, n_spin
          begin_index(:,i_irrep,i_spin)=-1_i4_kind
          was_found=.false.
          do i_dim = 1, ssym%dim(i_irrep)
             if((MO_status(i_dim, i_irrep,i_spin) == 0) .and..not.was_found(1)) then
                was_found(1)=.true.
                begin_index(1,i_irrep,i_spin)=i_dim
             else if ((MO_status(i_dim, i_irrep,i_spin) == 1) .and..not.was_found(2)) then
                was_found(2)=.true.
                begin_index(2,i_irrep,i_spin)=i_dim
             else if ((MO_status(i_dim, i_irrep,i_spin) == 2) .and..not.was_found(3)) then
                was_found(3)=.true.
                begin_index(3,i_irrep,i_spin)=i_dim
             end if
          end do
          end_index(:,i_irrep,i_spin)=-1_i4_kind
          was_found=.false.
          do i_dim = ssym%dim(i_irrep), 1, -1
             if((MO_status(i_dim, i_irrep,i_spin) == 0) .and..not.was_found(1)) then
                was_found(1)=.true.
                end_index(1,i_irrep,i_spin)=i_dim
             else if ((MO_status(i_dim, i_irrep,i_spin) == 1) .and..not.was_found(2)) then
                was_found(2)=.true.
                end_index(2,i_irrep,i_spin)=i_dim
             else if ((MO_status(i_dim, i_irrep,i_spin) == 2) .and..not.was_found(3)) then
                was_found(3)=.true.
                end_index(3,i_irrep,i_spin)=i_dim
             end if
          end do
       end do
    end do

    ! calculate total number of
    !  - full orbitals               n_full, 
    !  - partially filled orbitals   n_partial
    !  - empty orbitals              n_empty
    n_full    = 0_i4_kind
    n_partial = 0_i4_kind
    n_empty   = 0_i4_kind
    do i_irrep=1, n_irrep
       do i_spin=1, n_spin
          i_full    = 0_i4_kind
          i_partial = 0_i4_kind
          i_empty   = 0_i4_kind
          do i_dim=1, ssym%dim(i_irrep)
             if(MO_status(i_dim, i_irrep, i_spin)==0) i_full    = i_full    + 1
             if(MO_status(i_dim, i_irrep, i_spin)==1) i_partial = i_partial + 1
             if(MO_status(i_dim, i_irrep, i_spin)==2) i_empty   = i_empty   + 1
          end do

          i_rest = ssym%dim(i_irrep) - i_full - i_partial - i_empty

          n_full   (i_irrep, i_spin) = i_full
          n_partial(i_irrep, i_spin) = i_partial
          n_empty  (i_irrep, i_spin) = i_empty

       end do
    end do

    if (comm_i_am_master()) then
       write(output_unit,*) ""
       write(output_unit,*) " First and last indices [a,b] of "
       write(output_unit,*) "  1) fully     occupied levels"
       write(output_unit,*) "  2) partially occupied levels"
       write(output_unit,*) "  3) empty              levels"
       write(output_unit,*) ""
       write(output_unit,fmt=1000)
       do i_spin=1, n_spin
          do i_irrep=1, n_irrep
             write(output_unit,1010) i_spin,i_irrep,&
                  &begin_index(1,i_irrep,i_spin),&
                  &end_index  (1,i_irrep,i_spin),&
                  &begin_index(2,i_irrep,i_spin),&
                  &end_index  (2,i_irrep,i_spin),&
                  &begin_index(3,i_irrep,i_spin),&
                  &end_index  (3,i_irrep,i_spin)
          end do
       end do
1000   format(5X,"spin  irrep           (1)            (2)            (3)")
1010   format(7X,I1,5X,I1,5X,3(4X,"[",I4,",",I4,"]"))
       write(output_unit,*) ""
       write(output_unit,*) ""
       write(output_unit,*) " MO occupation as calculated in this module:"
       write(output_unit,*) "   status = -1 : level will be ignored"
       write(output_unit,*) "   status =  0 : fully occupied level"
       write(output_unit,*) "   status =  1 : partially occupied level"
       write(output_unit,*) "   status =  2 : empty level"
       write(output_unit,*) ""
       write(output_unit,*) "    spin  irrep    index          energy [eV]&
            &      occupation        status"
       do i_spin=1, n_spin
          do i_irrep=1, n_irrep
             print_dim=min(ssym%dim(i_irrep),num_spectrum_levels)
             do i_dim=print_dim,1,-1
                write(output_unit,fmt=1020) i_spin,i_irrep,i_dim,&
                     &eigval(i_irrep)%m(i_dim, i_spin)*hartree2ev,&
                     &occ_num(i_irrep)%m(i_dim,i_spin),&
                     &MO_status(i_dim,i_irrep,i_spin)
             end do
          end do
       end do
       write(output_unit,*) ""
1020   format(7X,I1,5X,I1,5X,I4,4X,F20.10,5X,ES10.4,11X,I2)
    end if

    if (comm_parallel()) then
       call comm_init_send(comm_all_other_hosts,msgtag_response_level_send)
       call pck(MO_status)
       call pck(begin_index)
       call pck(end_index)
       call pck(n_full)
       call pck(n_partial)
       call pck(n_empty)
       call comm_send()
    end if

  contains

    integer(i4_kind) function stOfMO(ir,is,id) result(r)
      integer(i4_kind), intent(in)  :: ir,is,id
      if     (occ_num(ir)%m(id,is) >= full_charge ) then
         r = 0 ! occupied
      elseif (occ_num(ir)%m(id,is) <= empty_charge ) then
         r = 2 ! empty
      else
         r = 1 ! partially
      end if
    end function stOfMO

  end subroutine resp_util_set_level_index

  !*************************************************************
  subroutine resp_util_upck_level_index()
    !  Purpose: 
    !  Receive the indices of the highest FULLY occupied level
    !  and of the lowest EMPTY level to all slaves.
    !  NOTE: This subroutine is called by "main_slave()" only !         
    !------------ Modules used ------------------- ---------------
    use constants, only: zero
    use msgtag_module, only: msgtag_response_level_send
    implicit none
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(i4_kind) :: max_dim, n_spin, n_irrep, alloc_stat
    !------------ Executable code ------------------------------------

    !## then allocate memory for the array MO_status

    ! find largest dimension of all irreps
    max_dim = maxval(ssym%dim(:))
    n_spin  = ssym%n_spin      ! number of spins
    n_irrep = ssym%n_irrep     ! number of irreps
    allocate(MO_status(max_dim,n_irrep,n_spin),&
         & begin_index(3,n_irrep,n_spin),&
         & end_index(3,n_irrep,n_spin),&
         & n_full(n_irrep,n_spin),n_partial(n_irrep,n_spin),&
         & n_empty(n_irrep,n_spin),stat=alloc_stat)
    ASSERT(alloc_stat==0)

    MO_status   = 0
    begin_index = 0
    end_index   = 0
    n_full      = 0
    n_partial   = 0
    n_empty     = 0

    call comm_save_recv(comm_master_host,msgtag_response_level_send)
    call upck(MO_status)
    call upck(begin_index)
    call upck(end_index)
    call upck(n_full)
    call upck(n_partial)
    call upck(n_empty)

#if 0
    print *,MO_status
    print *,begin_index
    print *,end_index
    print *,n_full
    print *,n_partial
    print *,n_empty
#endif

  end subroutine resp_util_upck_level_index
  !*************************************************************


  !*************************************************************
  subroutine resp_util_set_trans_bounds(i_case, i_irrep, i_spin,&
       & n_dim, occ_start, occ_end, unocc_start, unocc_end)
    ! Purpose: 
    ! For fixed irrep and spin direction calculate the correct 
    ! start and stop indices for the 3 cases "I_CASE":
    !       "occupied"                  "unoccupied"
    ! 1   full           MOs    ->   part. + empty  MOs
    ! 2   partially occ. MOs    ->   empty          MOs
    ! 3   partially occ. MOs    ->   partially occ. MOs 

    implicit none
    ! --- declaration of formal parameters ---------------------
    integer(kind=i4_kind), intent(in) :: i_case,i_irrep, i_spin,&
         & n_dim
    integer(kind=i4_kind), intent(out)           :: occ_start, occ_end
    integer(kind=i4_kind), intent(out)           :: unocc_start, unocc_end
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) ::  i_full, i_partial, i_empty, i_mo
    !------------ Executable code ------------------------------------

    ! define some internal abbreviations
    i_full    = n_full   (i_irrep, i_spin)
    i_partial = n_partial(i_irrep, i_spin)
    i_empty   = n_empty  (i_irrep, i_spin)

    select case (i_case)
       ! CASE 1: i=fully occupied levels
       !         s=partially occupied levels + empty levels
    case(1)
       occ_start   = begin_index(1,i_irrep,i_spin)
       occ_end     = end_index(1,i_irrep,i_spin)
       if( (i_partial/=0) .AND. (i_empty/=0)) then
          unocc_start = minval(begin_index(2:3,i_irrep,i_spin))
          unocc_end   = maxval(end_index(2:3,i_irrep,i_spin)) 
       else if( (i_partial/=0) .and. (i_empty==0)) then
          unocc_start = begin_index(2,i_irrep,i_spin)
          unocc_end   = end_index(2,i_irrep,i_spin)
       else if( (i_partial==0) .and. (i_empty/=0))  then
          unocc_start = begin_index(3,i_irrep,i_spin)
          unocc_end   = end_index(3,i_irrep,i_spin)
       else
          call error_handler("response_set_trans_bounds:case(1): there are no partially filled or empty orbitals !")
       end if

       ! CASE 2: i=partially occupied levels
       !         s=empty levels levels
    case(2)
       if ( (i_partial/=0) .AND. (i_empty/=0)) then
          occ_start   = begin_index(2,i_irrep,i_spin)
          occ_end     = end_index  (2,i_irrep,i_spin)
          unocc_start = begin_index(3,i_irrep,i_spin)
          unocc_end   = end_index  (3,i_irrep,i_spin)
       else
          call error_handler("response_set_trans_bounds:case(2): there are no partially filled or empty orbitals !")
       end if

       ! CASE 3: i=partially occupied levels
       !         s=partially occupied levels
       ! -> calculate upper block diagonal without main diagonal
       !    (there is no transition from one state in itself)
    case(3)
       if ( i_partial/=0 ) then
          occ_start   = begin_index(2,i_irrep,i_spin)
          ! usually we have occ_end = end_index(2,i_irrep,i_spin) - 1
          ! but just in case that there are holes between partially occupied states
          ! we search for the correct by scanning all indices.
          ! after the loop  occ_end holds that MO index which is the 
          ! last but one (=Vorletzter Index) with MO_status=1 i.e. partially filled
          do i_mo=begin_index(2,i_irrep,i_spin),end_index(2,i_irrep,i_spin)-1
             if(MO_status(i_mo,i_irrep, i_spin)/=0) occ_end = i_mo
          end do
          unocc_start = -1  ! not needed here
          unocc_end   = end_index(2,i_irrep,i_spin)
       else
          call error_handler("response_set_trans_bounds:case(3) there are no partially filled orbitals !")
       end if

       ! Otherwise we have an error
    case default
       call error_handler("response_set_trans_bounds: unknown value of i_case")

    end select

  end subroutine resp_util_set_trans_bounds

  !*************************************************************
  subroutine resp_util_calc_ou(i_ir_c,i_spin,dim_ou)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use symmetry_data_module, only: symmetry_data_n_irreps
    use clebsch_gordan,       only: cg=>cg_eliminated
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: i_ir_c, i_spin
    integer(kind=i4_kind), intent(OUT  ) :: dim_ou
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: occs, occe,unoccs, unocce,i_occ, i_unocc
    integer(kind=i4_kind) :: i_ir_a,i_ir_b, n_irr
    integer(kind=i4_kind) :: imlt, nmult
    !------------ Executable code ------------------------------------

    n_irr  = symmetry_data_n_irreps()
    dim_ou = 0_i4_kind
    do i_ir_a = 1, n_irr
       do i_ir_b = 1, n_irr

          call resp_util_borders(i_ir_a,i_ir_b,i_spin,occs,occe,unoccs,unocce)

          if ((occs == 0) .or. (unoccs==0)) cycle
          if ((occe == 0) .or. (unocce==0)) cycle

          nmult =  cg(i_ir_c,i_ir_a,i_ir_b)%mult

          do imlt = 1, nmult
             do i_occ = occs, occe
                do i_unocc = unoccs, unocce
                   if (abs(eigval(i_ir_a)%m(i_occ,i_spin)&
                        -eigval(i_ir_b)%m(i_unocc,i_spin))<min_diff) cycle
                   dim_ou = dim_ou + 1
                end do
             end do
          end do

       end do
    end do

  end subroutine resp_util_calc_ou
  !*************************************************************

  !*************************************************************
  subroutine resp_util_borders(i_ira,i_irb,i_spin,& !! in
       &                   occs,occe,unoccs,unocce) !! out
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN ) :: i_ira, i_irb, i_spin
    integer(kind=i4_kind), intent(OUT) :: occs,occe,unoccs,unocce
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)              :: shift, i_shift
    !------------ Executable code ------------------------------------

    occs = 0
    occe = 0
    unoccs = 0
    unocce = 0

    !! IMPORTANT!! LOGIC SHOULD BE IN THIS WAY

    if(begin_index(1,i_ira,i_spin) /=-1) then
       occs = begin_index(1,i_ira,i_spin)  !! start with full
    elseif(begin_index(2,i_ira,i_spin) /=-1) then
       occs = begin_index(2,i_ira,i_spin)  !! start with partial
    end if

    if (begin_index(2,i_irb,i_spin) /= -1) then 
       unoccs = begin_index(2,i_irb,i_spin) !! start with partial
    elseif(begin_index(3,i_irb,i_spin) /= -1) then
       unoccs = begin_index(3,i_irb,i_spin) !! start with empty
    endif

    if (end_index(2,i_ira,i_spin) /=-1) then
       occe = end_index(2,i_ira,i_spin) !! finished with partial
    elseif(end_index(1,i_ira,i_spin) /=-1) then
       occe = end_index(1,i_ira,i_spin) !! finished with full
    end if

    if (end_index(3,i_irb,i_spin) /=-1) then
       unocce = end_index(3,i_irb,i_spin) !! finished with empty
    elseif(end_index(2,i_irb,i_spin) /=-1) then
       unocce = end_index(2,i_irb,i_spin) !! finished with partial
    end if

    shift = unoccs
    do i_shift = unoccs, unocce
       if (abs(eigval(i_ira)%m(occe,i_spin)-eigval(i_irb)%m(i_shift,i_spin))<min_diff) then 
          shift = shift + 1
       end if
    end do
    unoccs = shift

    if (occs<0) occs=0
    if (occe<0) occe=0
    if (unoccs<0) unoccs=0
    if (unocce<0) unocce=0

  end subroutine resp_util_borders
  !*************************************************************


  !*************************************************************
  character(len=32) function resp_util_fname(name,i_ir,i_sp)
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
    resp_util_fname    = trim(fnm_char)//"_"//trim(irc_char)//"_"//trim(isp_char)
  end function resp_util_fname
  !*************************************************************

  subroutine resp_util_calc_transitions_v2(i_ira, i_irb, i_irc, i_spin, &
       & occ_times_unocc_dim, print_info)

    use clebsch_gordan, only: cg=>cg_eliminated, prod_bas

    !  Purpose: 
    ! For fixed irrep and spin direction calculate total number of 
    ! occupied times unoccupied orbitals by taking fractional 
    ! occupation numbers into account.
    !       "occupied"                  "unoccupied"
    ! 1   full           MOs    ->   part. + empty  MOs
    ! 2   partially occ. MOs    ->   empty          MOs
    ! 3   partially occ. MOs    ->   partially occ. MOs 
    !
    ! note: 
    ! for type (1) we have (N_full)      going into (N_partial+N_empty)
    ! for type (2) we have (N_partial)   going into (N_empty)
    ! for type (3) we have (N_partial-1) going into (N_partial-1)
    implicit none
    ! --- declaration of formal parameters ---------------------
    integer(kind=i4_kind), intent(in) :: i_ira, i_irb, i_irc, i_spin !!$, na, nb
    integer(kind=i4_kind), intent(out):: occ_times_unocc_dim
    logical, optional                 :: print_info
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: i_full, i_parta, i_partb, i_empty !!, i_resta, i_restb
    integer(kind=i4_kind) :: occs, occe, unocce, unoccs, i_occ, i_unocc
    !------------ Executable code ------------------------------------

    if (cg(i_irc,i_ira,i_irb)%mult==0) then
       occ_times_unocc_dim = 0
       return
    end if

    ! define some internal abbreviations
    
    i_full    = n_full   (i_ira, i_spin)
    i_parta   = n_partial(i_ira, i_spin)

    i_partb   = n_partial(i_irb, i_spin)
    i_empty   = n_empty  (i_irb, i_spin)

!!$    i_resta = na - i_full  - i_parta 
!!$    i_restb = nb - i_partb - i_empty 

    occ_times_unocc_dim = 0

    call resp_util_borders(i_ira,i_irb,i_spin,occs,occe,unoccs,unocce)
    
    if ((occs == 0) .or. (unoccs==0)) return
    if ((occe == 0) .or. (unocce==0)) return

    do i_occ = occs, occe
       do i_unocc = unoccs, unocce
          if (abs(eigval(i_ira)%m(i_occ,i_spin)&
               -eigval(i_irb)%m(i_unocc,i_spin))<min_diff) cycle
          occ_times_unocc_dim = occ_times_unocc_dim + 1
       end do
    end do

    if (comm_i_am_master()) then
       if(present(print_info)) then
          if(print_info) then
             call write_to_output_units("")
             call write_to_output_units("response_calc_transitions_v2: transitions for ")
             call write_to_output_units("   from ira ",inte=i_ira)
             call write_to_output_units("     to irb ",inte=i_irb)
             call write_to_output_units("     as irc ",inte=i_irc)
             call write_to_output_units("   spin  ",inte=i_spin)
#if 0
             call write_to_output_units("")
             call write_to_output_units("   n_full    =",inte=i_full)
             call write_to_output_units("   n_part a  =",inte=i_parta)
             call write_to_output_units("   n_part b  =",inte=i_partb)
             call write_to_output_units("   n_empty   =",inte=i_empty)
#endif
             call write_to_output_units("   Number of transitions =",inte=occ_times_unocc_dim)
             call write_to_output_units("")
          end if
       end if
    end if



  end subroutine resp_util_calc_transitions_v2


  !*************************************************************
  subroutine resp_util_buildr(NM,NS,WM,RM)
    !  Purpose: ..
    !------------ Modules used ----------------------------------
    use comm_module
    use msgtag_module, ONLY: msgtag_tddft_dipolem
    use xpack
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: NM, NS
    real(kind=r8_kind),    intent(IN   ) :: WM(:)
    real(kind=r8_kind),    intent(OUT  ) :: RM(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)           :: NST,NFN, i_proc, status
    real(kind=r8_kind), ALLOCATABLE :: AM(:)
    !------------ Executable code ------------------------------------
    if (comm_i_am_master()) then
       RM(1:NM) = WM

       ALLOCATE(AM(NS), STAT = status)
       ASSERT(status==0)

       do i_proc = 2, comm_get_n_processors()
          NST = NM + (i_proc-2) * NS + 1
          NFN = NM + (i_proc-1) * NS
          call comm_save_recv(i_proc,msgtag_tddft_dipolem)
          call upck(AM)
          RM(NST:NFN) = AM
       end do
       DEALLOCATE(AM, STAT=status)
       ASSERT(status==0)
    else
       call comm_init_send(comm_master_host,msgtag_tddft_dipolem)
       call pck(WM)
       call comm_send()
    end if
  end subroutine resp_util_buildr
  !***********************************************************

  !*************************************************************
  subroutine resp_util_buildi(NM,NS,WM,RM)
    !  Purpose: ..
    !------------ Modules used ----------------------------------
    use comm_module
    use msgtag_module, ONLY: msgtag_tddft_dipolem
    use xpack
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: NM, NS
    integer(kind=i4_kind), intent(IN   ) :: WM(:)
    integer(kind=i4_kind), intent(OUT  ) :: RM(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)               :: NST,NFN, i_proc, status
    integer(kind=i4_kind) , ALLOCATABLE :: AM(:)
    !------------ Executable code ------------------------------------

    if (comm_i_am_master()) then
       RM(1:NM) = WM

       ALLOCATE(AM(NS), STAT = status)
       ASSERT(status==0)

       do i_proc = 2, comm_get_n_processors()
          NST = NM + (i_proc-2) * NS + 1
          NFN = NM + (i_proc-1) * NS
          call comm_save_recv(i_proc,msgtag_tddft_dipolem)
          call upck(AM)
          RM(NST:NFN) = AM
       end do
       DEALLOCATE(AM, STAT=status)
       ASSERT(status==0)
    else
       call comm_init_send(comm_master_host,msgtag_tddft_dipolem)
       call pck(WM)
       call comm_send()
    end if
  end subroutine resp_util_buildi
  !***********************************************************



  !--------------- End of module -------------------------------------
end module resp_util_module
