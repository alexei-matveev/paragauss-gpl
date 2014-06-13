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
module  symmetry_data_module
!----------------------------------------------------------------
!
!  Purpose: contains information about symmetry as used in SCF-Part:
!             Number of Irreps, Point Group,
!             for each Irrep:
!                Name, Number of Partners,
!                Number of Basis Functions ("dimension")
!
!  Routines to read this information from input file and
!  to write it in formatted form to output are included.
!
!  The modul is filled by the symmetry part.
!
!  Symmetry adaption coefficients are contained in
!  unique_atom_module.
!
!  Data are module private and made accesible via Subroutines.
!
!  The informaion is available both on master and slaves.
!  Subroutines for packing and unpacking data are included.
!
!  A Metaindex i_ip that combines irrep and partner indices i_ir
!  and i_pa is also included. Mapping of i_ir and i_pa to i_ip
!  and vive versa is done by function symmetry_data_i_ip(i_ir,i_pa)
!  aand subroutine symmetry_data_unmap_i_ip(i_ip,i_ir,i_pa).
!  Value range of if is 1 to symmetry_data_n_ip().
!  Subroutine symmetry_data_setup() must be called before.
!
!  A logical array symmetry_data_dipoles_exist(i_ir,i_ir,i_xyz)
!  is also included that describes which blocks (i_ir,i_ir)
!  between the x, y and z component of dipole integrals are not
!  equal zero.
!
!  Author: Folke Noertemann
!  Date: 10/95
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: TB
! Date:   10/95
! Description: packing and unpacking routines added, restructured
!
! Modification (Please copy before editing)
! Author: TB
! Date:   1/96
! Description: routine for elliminating empty ireps added
!
! Modification (Please copy before editing)
! Author: TB
! Date:   9/96
! Description: removed n_spin and dimension from input namelists
!              Dimensions are calculated automatically.
!              The new equivalent to n_spin is now
!              spin_restricted in namelist main_options
!              However, for sake of compatibability, n_spin
!              remains in sym as well as in options module.
!
! Modification (Please copy before editing)
! Author: TB
! Date:   7/97
! Description: Added Metaindex i_ip that combines irrep and partner
!              indices and related Routines. (See above)
!              Added logical array symmetry_data_dipoles_exist.
!              Haeder of module totally rewritten.
!
! Modification (Please copy before editing)
! Author: MS
! Date:   8/97
! Description: The logical pseudo was added to ssym, to indicate if a
!              irrep is pseudo 2D
!
! Modification (Please copy before editing)
! Author: MM
! Date:   6/98
! Description: Projective Irreps were added
!
! Modification (Please copy before editing)
! Author:
! Date:
! Description:
!
!
!----------------------------------------------------------------

!------------ Modules used --------------------------------------
#include "def.h"
use type_module, only:&
     & IK=>i4_kind,&
     & RK=>r8_kind ! type specification parameters
use operations_module

implicit none
save
private         ! by default, everything is private
!== Interrupt end of public interface of module =================

!------------ Declaration of types ------------------------------
type, public :: sym
   ! vector irreps
   integer(IK)          :: n_irrep    ! number of IRREPS
   integer(IK),pointer  :: dim(:)     ! dim of IRREP i
   integer(IK),pointer  :: partner(:) ! number of partners
   logical, pointer               :: pseudo(:)  ! .true. if we have pseudo 2D irrep
   character(len=12), pointer :: name(:)    ! name of IRREP
   ! projective irreps
   integer(IK)          :: n_proj_irrep    ! number of projective IRREPS
   integer(IK),pointer  :: dim_proj(:)     ! dim of proj. IRREP i
   integer(IK),pointer  :: partner_proj(:) ! number of partners
   character(len=12), pointer :: name_proj(:)    ! name of IRREP
   real(RK),pointer  :: jz(:)           ! for axial groups: jz
   ! general information about pointgroup
   character(len=4) :: point_group

   ! FIXME: see the strange logic with n_spin_set
   integer(IK)          :: n_spin = 1 ! spin-restricted or not (1 or 2 )

   integer(IK)          :: totalsymmetric_irrep  ! index
   logical              :: data_allocated
end type sym


!------------ Declaration of public variables -------------------

!
! One way to  get symmetry info, such as the  number and dimensions of
! irreps is to import this struct:
!
!       use symmetry_data_module, only: ssym
!
type(sym), public, protected :: ssym

integer(IK)             , public :: number_of_irreps = 0
integer(IK), allocatable, public :: irrep_dimensions(:) ! (>=number_of_irreps), includes empty ones!

integer(IK)         , public :: symmetry_data_dip_irrep_mult(20)=0 ! (:n_irr) is used
real(RK)            , public :: symmetry_data_dip_components(3,3)
logical, allocatable, public :: symmetry_data_dipoles_exist(:,:,:)
logical, allocatable, public :: symmetry_data_pdipoles_exist(:,:,:)
!(N_irrep,Nirrep,3) , third index: x,y,z
! describes which blocks (i_ir,i_ir) between the x, y and z component
! of dipole integrals are not equal zero.

!interface symmetry_data_set
!   module procedure symmetry_data_set_irreps
!   module procedure symmetry_data_set_spin_orbit
!end interface

!------------ public functions and subroutines ------------------

public :: symmetry_data_setup
public :: symmetry_data_close

public &
     symmetry_data_n_partners, &
     symmetry_data_irrepname, &
     symmetry_data_irrepname_proj, &
     symmetry_data_dimension, &
     get_totalsymmetric_irrep, &
     symmetry_data_set,&
     symmetry_data_set_pcoupling,&
     symmetry_data_set_cccoupling,&
     find_totalsymmetric_irrep, &
     symmetry_data_n_spin, &
     symmetry_data_elim_irreps, &
     symmetry_data_group_read, &
     symmetry_data_group_write, &
     symmetry_data_point_group, &
     symmetry_data_write_formatted,&
     symmetry_data_i_ip, &
     symmetry_data_i_ip_proj, &
     symmetry_data_n_ip, &
     symmetry_data_n_ip_proj, &
     symmetry_data_unmap_i_ip, &
     symmetry_data_pseudo,&
     symmetry_data_dimension_proj, &
     symmetry_data_n_partners_proj, &
     symmetry_data_jz, &
     symmetry_data_elim_pirreps,&
     & symmetry_data_n_irreps,&          ! NO spin-orbit
     & symmetry_data_n_proj_irreps,&     !    spin-orbit
     & symmetry_data_n_irr,&             ! always, if op_spin_orbit is set
     & symmetry_data_get_pcoupling,&
     & symmetry_data_get_cccoupling

!================================================================
! End of public interface of module
!================================================================

!------------ Declaration of private constants and variables ----

integer(IK), pointer, private :: pcoupling_all(:); !...(number_of_irreps)
integer(IK), pointer, private :: pcoupling(:);     !...(number_of_not_empty_irreps)
integer(IK), pointer, private :: cccoupling_all(:); !...(number_of_irreps)
integer(IK), pointer, private :: cccoupling(:);     !...(number_of_not_empty_irreps)

logical, private :: op_spin_orbit = .false.

character(len=12), private :: point_group ! length 12 kept for compatibility
character(len=4), parameter, private :: df_point_group = "C1  "

namelist /symmetry_group/ point_group

logical, private :: n_spin_set = .false. ! FIXME: must die!
logical, private :: projective = .false.

integer, allocatable, private :: &
     first_ip_of_ir(:), & ! (ssym%n_irrep)
     ir_and_pa_of_ip(:,:), & ! (2,n_ip)
     first_ip_of_ir_proj(:), & ! (ssym%n_irrep)
     ir_and_pa_of_ip_proj(:,:) ! (2,n_ip)
integer, private :: n_ip


!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine symmetry_data_setup(options_spin_orbit)
    ! purpose : setup work:
    !           Mapping to Metaindex i_ip that combines irrep and
    !           partner indices i_ir and i_pa is prepared.
    ! Must be called after symmetry_data_set.
    ! Called by: initialize master, initialize slave
    implicit none
    logical, intent(in) :: options_spin_orbit
    !** End of interface ***************************************

    integer :: i_ir, i_pa, i_ip, alloc_stat

    call symmetry_data_bcast()

    allocate(first_ip_of_ir(ssym%n_irrep),stat=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("symmetry_data_setup: allocataion (1) failed")
    i_ip = 1
    do i_ir = 1, ssym%n_irrep
       first_ip_of_ir(i_ir) = i_ip
       i_ip = i_ip + ssym%partner(i_ir)
    enddo
    n_ip = i_ip - 1
    allocate(ir_and_pa_of_ip(2,n_ip), stat=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("symmetry_data_setup: allocataion (2) failed")
    i_ip = 1
    do i_ir = 1, ssym%n_irrep
       do i_pa = 1, ssym%partner(i_ir)
          ir_and_pa_of_ip(1,i_ip) = i_ir
          ir_and_pa_of_ip(2,i_ip) = i_pa
          i_ip = i_ip + 1
       enddo
    enddo

    if ( options_spin_orbit ) then
        call symmetry_data_setup_proj()
    endif
  end subroutine symmetry_data_setup
  !*************************************************************

  !*************************************************************
  subroutine symmetry_data_setup_proj()
    ! purpose : setup work:
    !           Mapping to Metaindex i_ip that combines irrep and
    !           partner indices i_ir and i_pa is prepared.
    ! Must be called after symmetry_data_set.
    ! Called by: initialize master, initialize slave
    !** End of interface ***************************************
    integer :: i_ir, i_pa, i_ip, alloc_stat
    allocate(first_ip_of_ir_proj(ssym%n_proj_irrep),stat=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("symmetry_data_setup: allocataion (1) failed")
    i_ip = 1
    do i_ir = 1, ssym%n_proj_irrep
       first_ip_of_ir_proj(i_ir) = i_ip
       i_ip = i_ip + ssym%partner_proj(i_ir)
    enddo
    n_ip = i_ip - 1
    allocate(ir_and_pa_of_ip_proj(2,n_ip), stat=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("symmetry_data_setup: allocataion (2) failed")
    i_ip = 1
    do i_ir = 1, ssym%n_proj_irrep
       do i_pa = 1, ssym%partner_proj(i_ir)
          ir_and_pa_of_ip_proj(1,i_ip) = i_ir
          ir_and_pa_of_ip_proj(2,i_ip) = i_pa
          i_ip = i_ip + 1
       enddo
    enddo
  end subroutine symmetry_data_setup_proj
  !*************************************************************

  !*************************************************************
  integer function symmetry_data_i_ip(i_ir,i_pa)
    ! purpose : returns Metaindex i_ip that combines irrep and
    !           partner indices i_ir and i_pa
    integer(IK),intent(in)  :: i_ir,i_pa
    !** End of interface ***************************************
    symmetry_data_i_ip = first_ip_of_ir(i_ir) + i_pa - 1
  end function symmetry_data_i_ip
  !*************************************************************

  !*************************************************************
  integer function symmetry_data_i_ip_proj(i_ir,i_pa)
    ! purpose : returns Metaindex i_ip that combines irrep and
    !           partner indices i_ir and i_pa
    integer(IK),intent(in)  :: i_ir,i_pa
    !** End of interface ***************************************
    symmetry_data_i_ip_proj = first_ip_of_ir_proj(i_ir) + i_pa - 1
  end function symmetry_data_i_ip_proj
  !*************************************************************

  !*************************************************************
  subroutine symmetry_data_unmap_i_ip(i_ip,i_ir,i_pa)
    ! purpose : returns irrep and partner indices i_ir and i_pa
    !           for a given Metaindex i_ip
    integer(IK),intent(in)  :: i_ip
    integer(IK),intent(out) :: i_ir,i_pa
    !** End of interface ***************************************
    i_ir = ir_and_pa_of_ip(1,i_ip)
    i_pa = ir_and_pa_of_ip(2,i_ip)
  end subroutine symmetry_data_unmap_i_ip
  !*************************************************************

  !*************************************************************
  subroutine symmetry_data_unmap_i_ip_proj(i_ip,i_ir,i_pa)
    ! purpose : returns irrep and partner indices i_ir and i_pa
    !           for a given Metaindex i_ip
    integer(IK),intent(in)  :: i_ip
    integer(IK),intent(out) :: i_ir,i_pa
    !** End of interface ***************************************
    i_ir = ir_and_pa_of_ip_proj(1,i_ip)
    i_pa = ir_and_pa_of_ip_proj(2,i_ip)
  end subroutine symmetry_data_unmap_i_ip_proj
  !*************************************************************

  !*************************************************************
  integer function symmetry_data_n_ip()
    ! purpose : returns maximum of Metaindex i_ip that combines
    !           irrep and partner indices i_ir and i_pa
    !** End of interface ***************************************
    symmetry_data_n_ip = n_ip
  end function symmetry_data_n_ip
  !*************************************************************

  !*************************************************************
  integer function symmetry_data_n_ip_proj()
    ! purpose : returns maximum of Metaindex i_ip that combines
    !           irrep and partner indices i_ir and i_pa
    !** End of interface ***************************************
    symmetry_data_n_ip_proj = n_ip
  end function symmetry_data_n_ip_proj
  !*************************************************************

  !*************************************************************
  subroutine free_sym(ssym_dummy,proj_only,vec_only)
    ! purpose: free the allocated space in ssym
    type(sym)         :: ssym_dummy
    logical,optional  :: proj_only
    logical,optional  :: vec_only
    !** End of interface ***************************************
    integer(IK)   :: dealloc_stat

    if (present(proj_only)) then
       ! free only projective irreps
       deallocate(ssym_dummy%dim_proj,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation of ssym_dummy%dim_proj failed")

       deallocate(ssym_dummy%partner_proj,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation of ssym_dummy%partner_projfailed")

       deallocate(ssym_dummy%name_proj,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation of ssym_dummy%name_proj failed")

       deallocate(ssym_dummy%jz,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation of ssym_dummy%jz failed")
    elseif (present(vec_only)) then
       ! free only vector irreps
       deallocate(ssym_dummy%dim,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation failed (1)")

       deallocate(ssym_dummy%partner,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation failed (2)")

       deallocate(ssym_dummy%name,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation failed (3)")

       deallocate(ssym_dummy%pseudo,STAT=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_SYM: deallocation pseudo failed")
    else
       ! free either irreps
       if(associated(ssym_dummy%dim)) then
          deallocate(ssym_dummy%dim,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation failed (1)")
       end if

       if(associated(ssym_dummy%partner)) then
          deallocate(ssym_dummy%partner,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation failed (2)")
       end if

       if(associated(ssym_dummy%name)) then
          deallocate(ssym_dummy%name,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation failed (3)")
       end if


       if (projective) then
          deallocate(ssym_dummy%dim_proj,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation of ssym_dummy%dim_proj failed")

          deallocate(ssym_dummy%partner_proj,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation of ssym_dummy%partner_projfailed")

          deallocate(ssym_dummy%name_proj,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation of ssym_dummy%name_proj failed")

          deallocate(ssym_dummy%jz,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation of ssym_dummy%jz failed")
       endif


       if(associated(ssym_dummy%pseudo)) then
          deallocate(ssym_dummy%pseudo,STAT=dealloc_stat)
          if(dealloc_stat.ne.0) call error_handler &
               ("FREE_SYM: deallocation pseudo failed")
       end if
    endif

    ssym_dummy%data_allocated = .false.
  end subroutine free_sym
  !*************************************************************

  !*************************************************************
  integer(IK) function symmetry_data_n_spin()
    ! purpose: returns n_spin ( 1: closed shell, 2: open shell)
    !** End of interface ***************************************
    if ( .not. n_spin_set ) then
       ssym%n_spin = 1
       n_spin_set = .true.
    endif
    symmetry_data_n_spin = ssym%n_spin
   end function symmetry_data_n_spin
   !*************************************************************

  !*************************************************************
   character(len=4) function symmetry_data_point_group()
   ! purpose: returns pointgroup
     !** End of interface ***************************************
     symmetry_data_point_group = ssym%point_group
   end function symmetry_data_point_group
   !*************************************************************

  !*************************************************************
   integer(IK) function symmetry_data_n_irreps()
     ! purpose: returns number of irreps
     !** End of interface ***************************************
     symmetry_data_n_irreps = ssym%n_irrep
   end function symmetry_data_n_irreps
   !*************************************************************

  !*************************************************************
   integer(IK) function symmetry_data_n_proj_irreps()
     ! purpose: returns number of irreps
     !** End of interface ***************************************
   symmetry_data_n_proj_irreps = ssym%n_proj_irrep
   end function symmetry_data_n_proj_irreps

   function symmetry_data_n_irr() result(n_irr)
     integer(IK) :: n_irr !<<< result
     ! *** end of interface ***

     if(op_spin_orbit)then
        n_irr = symmetry_data_n_proj_irreps()
     else
        n_irr = symmetry_data_n_irreps()
     endif
   end function symmetry_data_n_irr
   !*************************************************************

   !*************************************************************
   integer(IK) function symmetry_data_n_partners(index)
     ! purpose: returns number of partners of irreps
     !------------ Declaration of formal parameters -------------
     integer(IK),intent(in) :: index ! of irrep
     !** End of interface ***************************************
     character(len=10) number
     if (index .lt. 1 .or. index .gt. ssym%n_irrep) then
        write(number,*) index
        call error_handler( &
             "symmetry_data_n_partners: wrong index " &
             // trim(number) )
     endif
     symmetry_data_n_partners = ssym%partner(index)
   end function symmetry_data_n_partners
   !*************************************************************

   !*************************************************************
   integer(IK) function symmetry_data_n_partners_proj(index)
     ! purpose: returns number of partners of irreps
     !------------ Declaration of formal parameters -------------
     integer(IK),intent(in) :: index ! of irrep
     !** End of interface ***************************************
     character(len=10) number
     if (index .lt. 1 .or. index .gt. ssym%n_proj_irrep) then
        write(number,*) index
        call error_handler( &
             "symmetry_data_n_partners_proj: wrong index " &
             // trim(number) )
     endif
     symmetry_data_n_partners_proj = ssym%partner_proj(index)
   end function symmetry_data_n_partners_proj
   !*************************************************************

   !*************************************************************
   real(RK) function symmetry_data_jz(index)
   ! purpose: returns number of partners of irreps
   !------------ Declaration of formal parameters -------------
   integer(IK),intent(in) :: index ! of irrep
   !** End of interface ***************************************
   character(len=10) number
   if (index .lt. 1 .or. index .gt. ssym%n_proj_irrep) then
      write(number,*) index
      call error_handler( &
           "symmetry_data_n_partners_proj: wrong index " &
           // trim(number) )
   endif
   symmetry_data_jz = ssym%jz(index)
   end function symmetry_data_jz
  !*************************************************************

   !*************************************************************
   character(len=12) function symmetry_data_irrepname(index)
     ! purpose: returns name of irrep
     !------------ Declaration of formal parameters -------------
     integer(IK),intent(in) :: index ! of irrep
     !** End of interface ***************************************
     character(len=10) number
     if (index .lt. 1 .or. index .gt. ssym%n_irrep) then
        write(number,*) index
        call error_handler( &
             "symmetry_data_irrepname: wrong index " &
             // trim(number) )
     endif
     symmetry_data_irrepname = ssym%name(index)
   end function symmetry_data_irrepname
   !*************************************************************

   !*************************************************************
   character(len=12) function symmetry_data_irrepname_proj(index)
     ! purpose: returns name of irrep
     !------------ Declaration of formal parameters -------------
     integer(IK),intent(in) :: index ! of irrep
     !** End of interface ***************************************
     character(len=10) number
     if (index .lt. 1 .or. index .gt. ssym%n_proj_irrep) then
        write(number,*) index
        call error_handler( &
             "symmetry_data_irrepname: wrong index " &
             // trim(number) )
     endif
     symmetry_data_irrepname_proj = ssym%name_proj(index)
   end function symmetry_data_irrepname_proj
   !*************************************************************

   !*************************************************************
   integer function get_totalsymmetric_irrep()
     ! purpose : returns index of totalsymmetric irrep
     !** End of interface ***************************************
     get_totalsymmetric_irrep = ssym%totalsymmetric_irrep
   end function get_totalsymmetric_irrep
   !*************************************************************

   !*************************************************************
   integer(IK) function symmetry_data_dimension(index)
     ! purpose: returns dimension of irreps
     !------------ Declaration of formal parameters -------------
     integer(IK),intent(in) :: index ! of irrep
     !** End of interface ***************************************
     character(len=10) number
     if (index .lt. 1 .or. index .gt. ssym%n_irrep) then
        write(number,'(i5)') index
        call error_handler( &
             "symmetry_data_dimension: wrong index " &
             // trim(number) )
     endif
     symmetry_data_dimension = ssym%dim(index)
   end function symmetry_data_dimension
   !*************************************************************

   !*************************************************************
   integer(IK) function symmetry_data_dimension_proj(index)
     ! purpose: returns dimension of irreps
     !------------ Declaration of formal parameters -------------
     integer(IK),intent(in) :: index ! of irrep
     !** End of interface ***************************************
     character(len=10) number
     if (index .lt. 1 .or. index .gt. ssym%n_proj_irrep) then
        write(number,*) index
        call error_handler( &
             "symmetry_data_dimension: wrong index " &
             // trim(number) )
     endif
     symmetry_data_dimension_proj = ssym%dim_proj(index)
   end function symmetry_data_dimension_proj
   !*************************************************************

   !*************************************************************
   logical function symmetry_data_pseudo(index)
     ! purpose: returns true if irrep is pseudo 2d and false else
     !------------ Declaration of formal parameters -------------
     integer(IK),intent(in) :: index ! of irrep
     !** End of interface ***************************************
     character(len=10) number
     if (index .lt. 1 .or. index .gt. ssym%n_irrep) then
        write(number,*) index
        call error_handler( &
             "symmetry_data_dimension: wrong index " &
             // trim(number) )
     endif
     symmetry_data_pseudo = ssym%pseudo(index)
   end function symmetry_data_pseudo
   !************************************************************

  subroutine symmetry_data_bcast()
    use comm, only: comm_rank
    use comm_module, only: comm_init_send, comm_send, &
        comm_all_other_hosts, comm_save_recv, &
        comm_master_host
    use msgtag_module, only: msgtag_packed_message
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        call comm_init_send(comm_all_other_hosts, msgtag_packed_message)
        call symmetry_data_pack()
        call comm_send()
    else
        call comm_save_recv(comm_master_host, msgtag_packed_message)
        call symmetry_data_unpack()
    endif
  end subroutine symmetry_data_bcast

   !*************************************************************
   subroutine symmetry_data_pack
     ! purpose: packs spvm into pvm buffer
     !** End of interface ***************************************
     !------------ Modules used -----------------------------------
     use commpack_module
     use xpack
     !------------ Declaration of local variables -----------------
     integer(IK)  :: info, i, j
     logical :: alloc_sym_dip
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
     if ( .not. n_spin_set ) then
        ssym%n_spin = 1
        n_spin_set = .true.
     endif

     ! op_spin_orbit:
     call pck(op_spin_orbit)

     if(op_spin_orbit)then
        ! pcoupling_all:
        call pck(size(pcoupling_all))
        ! P coupling:
        call pck(pcoupling_all)
        ! CC coupling:
        call pck(cccoupling_all)

        ! pcoupling, cccoupling size:
        call pck(size(pcoupling))
     endif

     call pck(size(irrep_dimensions))
     call pck(irrep_dimensions)
     call pck(number_of_irreps) ! of non-empty irreps!

     ! dipole symmetry adaption coeffs and multiplicities:
     call pck(symmetry_data_dip_irrep_mult)
     call pck(symmetry_data_dip_components)

     call commpack(ssym%n_irrep,info)
     if (info .ne. 0) call error_handler("symmetry_data_pack: n_irrep")
     call commpack(ssym%totalsymmetric_irrep,info)
     if (info .ne. 0) call error_handler("symmetry_data_pack: totalsymmetric_irrep")
     call commpack(ssym%dim,ssym%n_irrep,1,info)
     if (info .ne. 0) call error_handler("symmetry_data_pack: dimension")
     call commpack(ssym%partner,ssym%n_irrep,1,info)
     if (info .ne. 0) call error_handler("symmetry_data_pack: patrner")
     do i=1,ssym%n_irrep
        call commpack(ssym%name(i),info)
        if (info .ne. 0) call error_handler("symmetry_data_pack: name")
     enddo
     do i=1,ssym%n_irrep
        call commpack(ssym%pseudo(i),info)
        if (info .ne. 0) call error_handler("symmetry_data_pack: pseudo")
     end do
     if (info .ne. 0) call error_handler("symmetry_data_pack: pseudo")
     call commpack(ssym%point_group,info)
     if (info .ne. 0) call error_handler("symmetry_data_pack: point_group")
     call commpack(ssym%n_spin,info)
     if (info .ne. 0) call error_handler("symmetry_data_pack: n_spin")
     ! now pack projective irreps
     call commpack(projective,info)
     if (info .ne. 0) call error_handler("symmetry_data_pack: projective")
     if (projective) then
        call commpack(ssym%n_proj_irrep,info)
        if (info .ne. 0) call error_handler("symmetry_data_pack: n_irrep")
        call commpack(ssym%dim_proj,ssym%n_proj_irrep,1,info)
        if (info .ne. 0) call error_handler("symmetry_data_pack: dimension")
        call commpack(ssym%partner_proj,ssym%n_proj_irrep,1,info)
        if (info .ne. 0) call error_handler("symmetry_data_pack: patrner")
        do i=1,ssym%n_proj_irrep
           call commpack(ssym%jz(i),info)
           if (info .ne. 0) call error_handler("symmetry_data_pack: jz")
        enddo
        do i=1,ssym%n_proj_irrep
           call commpack(ssym%name_proj(i),info)
           if (info .ne. 0) call error_handler("symmetry_data_pack: name")
        enddo
     endif


     alloc_sym_dip= .false.
     if ( allocated(symmetry_data_dipoles_exist) ) then
        alloc_sym_dip= .true.
        call commpack(alloc_sym_dip,info)
        if (info .ne. 0) call error_handler("symmetry_data_pack: n_spin")
        do j = 1, 3
           do i=1,ssym%n_irrep
              call commpack(symmetry_data_dipoles_exist(:,i,j),ssym%n_irrep,1,info)
              if (info .ne. 0) call error_handler("symmetry_data_pack: dipoles_exist")
           enddo
           if (projective) then
              do i=1,ssym%n_proj_irrep
                 call commpack(symmetry_data_pdipoles_exist(:,i,j),ssym%n_proj_irrep,1,info)
                 if (info .ne. 0) call error_handler("symmetry_data_pack: dipoles_exist")
              enddo
           endif
        enddo
     else
        call commpack(alloc_sym_dip,info)
     endif
   end subroutine symmetry_data_pack
   !*************************************************************

   !*************************************************************
   subroutine symmetry_data_unpack
     ! purpose: packs ssym into pvm buffer
     !** End of interface ***************************************
     !------------ Modules used -----------------------------------
     use error_module
     use commpack_module
     use xpack
     !------------ Declaration of local variables -----------------
     integer(IK)  :: info, i, stat, j, n
     logical :: alloc_sym_dip
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
!     call buggy("symmetry_data_unpack: ENTERED")
     n_spin_set = .true.

     ! op_spin_orbit:
     call upck(op_spin_orbit)

     if(op_spin_orbit)then
        ! pcoupling_all:
        call upck(n)

        allocate(&
             & pcoupling_all(n),&
             & cccoupling_all(n),&
             & STAT=stat)
        call error(stat,"symmetry_data_unpack: alloc pcoupling_all failed")

        call upck(pcoupling_all)
        call upck(cccoupling_all)

        ! pcoupling:
        call upck(n)

        pcoupling => pcoupling_all(1:n)
        cccoupling => cccoupling_all(1:n)
     endif

     call upck(number_of_irreps)
     allocate(irrep_dimensions(number_of_irreps),stat=stat)
     ASSERT(stat==0)
     call upck(irrep_dimensions)
     call upck(number_of_irreps) ! of non-empty irreps!

     ! dipole symmetry adaption coeffs and multiplicities:
     call upck(symmetry_data_dip_irrep_mult)
     call upck(symmetry_data_dip_components)

     call communpack(ssym%n_irrep,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: n_irrep")
     allocate ( ssym%dim(ssym%n_irrep), &
                ssym%partner(ssym%n_irrep), &
                ssym%name(ssym%n_irrep), &
                ssym%pseudo(ssym%n_irrep), &
                stat = stat )
     if ( stat .ne. 0 ) call error_handler( &
          "symmetry_data_unpack: allocate failed" )
     call communpack(ssym%totalsymmetric_irrep,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: totalsymmetric_irrep")
     call communpack(ssym%dim,ssym%n_irrep,1,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: dimension")
     call communpack(ssym%partner,ssym%n_irrep,1,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: patrner")
     do i=1,ssym%n_irrep
        call communpack(ssym%name(i),info)
        if (info .ne. 0) call error_handler("symmetry_data_unpack: name")
     enddo
     do i=1,ssym%n_irrep
        call communpack(ssym%pseudo(i),info)
        if (info .ne. 0) call error_handler("symmetry_data_unpack: pseudo")
     enddo
     call communpack(ssym%point_group,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: point_group")
     call communpack(ssym%n_spin,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: n_spin")

     ! now read projective irreps
!     call buggy("symmetry_data_unpack: unpack projective")
     call communpack(projective,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: projective")
     if (projective) then
!        call buggy("symmetry_data_unpack: now unpack projective irreps")
        call communpack(ssym%n_proj_irrep,info)
!        call buggy("symmetry_data_unpack: dimension of proj. irrep: ",int=ssym%n_proj_irrep)
        if (info .ne. 0) call error_handler("symmetry_data_unpack: n_irrep")
!        call buggy("symmetry_data_unpack: allocate stuff")
        allocate ( ssym%dim_proj(ssym%n_proj_irrep), &
             ssym%partner_proj(ssym%n_proj_irrep), &
             ssym%name_proj(ssym%n_proj_irrep), &
             ssym%jz(ssym%n_proj_irrep), &
             stat = stat )
        if ( stat .ne. 0 ) call error_handler( &
             "symmetry_data_unpack: allocate failed" )
!        call buggy("symmetry_data_unpack: unpack sizes")
        call communpack(ssym%dim_proj,ssym%n_proj_irrep,1,info)
        if (info .ne. 0) call error_handler("symmetry_data_unpack: dimension_proj")
        call communpack(ssym%partner_proj,ssym%n_proj_irrep,1,info)
        if (info .ne. 0) call error_handler("symmetry_data_unpack: patrner_proj")

        do i=1,ssym%n_proj_irrep
           call communpack(ssym%jz(i),info)
           if (info .ne. 0) call error_handler("symmetry_data_unpack: jz")
        enddo
        do i=1,ssym%n_proj_irrep
           call communpack(ssym%name_proj(i),info)
           if (info .ne. 0) call error_handler("symmetry_data_unpack: name_proj")
        enddo
     endif
!     call buggy("symmetry_data_unpack: completed unpacking projective")
     call communpack(alloc_sym_dip,info)
     if (info .ne. 0) call error_handler("symmetry_data_unpack: n_spin")
     if (alloc_sym_dip) then
!        call buggy("symmetry_data_unpack: unpack dipole information")
        allocate(symmetry_data_dipoles_exist(ssym%n_irrep,ssym%n_irrep,3),stat=stat)
        if (projective) then
           allocate(symmetry_data_pdipoles_exist(ssym%n_proj_irrep,ssym%n_proj_irrep,3),stat=stat)
        endif

        if (stat .ne. 0) call error_handler("symmetry_data_unpack: allocate dipoles_exist failed")
        do j = 1, 3
           do i=1,ssym%n_irrep
              call communpack(symmetry_data_dipoles_exist(:,i,j),ssym%n_irrep,1,info)
              if (info .ne. 0) call error_handler("symmetry_data_unpack: dipoles_exist")
           enddo
           if (projective) then
              do i=1,ssym%n_proj_irrep
                 call communpack(symmetry_data_pdipoles_exist(:,i,j),ssym%n_proj_irrep,1,info)
                 if (info .ne. 0) call error_handler("symmetry_data_unpack: dipoles_exist")
              enddo
           endif
        enddo
     endif
!     call buggy("symmetry_data_unpack: LEAVING")
   end subroutine symmetry_data_unpack

   !*************************************************************

   subroutine symmetry_data_close()
     ! purpose: deallocates arrays in ssym
     !** End of interface ***************************************

     call free_sym(ssym)

     if ( allocated(symmetry_data_dipoles_exist) ) &
          deallocate(symmetry_data_dipoles_exist)
     if ( allocated(symmetry_data_pdipoles_exist) ) &
          deallocate(symmetry_data_pdipoles_exist)
     if ( allocated(first_ip_of_ir) ) &
          deallocate(ir_and_pa_of_ip,first_ip_of_ir)
     if ( allocated(first_ip_of_ir_proj) ) &
          deallocate(ir_and_pa_of_ip_proj,first_ip_of_ir_proj)
     if ( allocated(irrep_dimensions) )then
        deallocate(irrep_dimensions)
     endif

     ! restore default values, did we forget something?:
     op_spin_orbit = .false.
     projective = .false.
     n_spin_set = .false.
     number_of_irreps = 0
   end subroutine symmetry_data_close

   !*************************************************************

   subroutine symmetry_data_group_read
     ! purpose: reads pointgroup from input file
     !** End of interface ***************************************
     !------------ Modules used -----------------------------------
     use input_module
#ifdef WITH_EFP
     use qmmm_interface_module, only: efp
#endif
     !------------ Declaration of local variables -----------------
     integer(IK)                :: status, unit
     !------------ Declaration of subroutines used ----------------
     external error_handler
     !------------ Executable code --------------------------------
     DPRINT 'sdm::symmetry_data_group_read: entered'
     point_group = df_point_group ! missing blanks are added by =
     DPRINT 'sdm::symmetry_data_group_read: df_point_group=',point_group
     if ( input_line_is_namelist("symmetry_group") ) then
        call input_read_to_intermediate
        unit = input_intermediate_unit()
        read(unit, nml=symmetry_group, iostat=status)
        if (status .gt. 0) call input_error( &
             "symmetry_data_read: namelist symmetry_group")
     endif

#ifdef WITH_EFP
     if (efp .and. trim(point_group) /= "C1") call error_handler &
        ("symmetry_data_read: EFP method currently runs without symmetry")
#endif

     ssym%point_group = adjustl(point_group) ! trailing blanks ignored by =
     DPRINT 'sdm::symmetry_data_group_read: read>'//point_group//'<'

   end subroutine symmetry_data_group_read

   !*************************************************************

   subroutine symmetry_data_group_write(iounit)
     ! purpose: writes point group to iounit in input format
     !------------ Modules used -----------------------------------
     use echo_input_module
     !use operations_module, only: operations_echo_input_level
     !------------ Declaration of formal parameters -------------
     integer, intent(in) :: iounit
     !** End of interface ***************************************
     !------------ Executable code ------------------------------

     word_format = '("    ",a," = ",a6:" # ",a)' ! including quotes

     call start("SYMMETRY_GROUP","SYMMETRY_DATA_GROUP_WRITE", &
          iounit,operations_echo_input_level)
     call word("POINT_GROUP",ssym%point_group,df_point_group)
     call stop()

   end subroutine symmetry_data_group_write
   !*************************************************************

   !*************************************************************
   subroutine symmetry_data_write_formatted(output_unit)
     ! purpose: writes data to output_unit in a format suited
     !          for the output file
     !------------ Declaration of formal parameters -------------
     integer, intent(in) :: output_unit
     !** End of interface ***************************************
     !------------ Declaration of local variables ---------------
     integer(IK) :: status, i, n_ob, i_xyz, i_ir1, i_ir2
     !------------ Declaration of subroutines used --------------
     external error_handler
     !------------ Executable code ------------------------------
     write(output_unit, *, iostat=status)
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     write(output_unit, *, iostat=status) "The point group ", &
          ssym%point_group, " results in ", ssym%n_irrep, &
          "non_empty IRREPs"
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     write(output_unit, *, iostat=status) &
          "Index, Name, Number of Partners and &
          &Number of basis functions for each IRREP"
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     n_ob = 0
     do i=1,ssym%n_irrep

        write(output_unit, "(I5,2X,A,2I6)", iostat=status) &
             i, ssym%name(i), ssym%partner(i), ssym%dim(i)
        if (status .gt. 0) call error_handler( &
             "symmetry_data_write_formatted")

        n_ob = n_ob + ssym%partner(i) * ssym%dim(i)
     enddo

     write(output_unit, *, iostat=status) &
          "Number of totally symmetric IRREP: ", ssym%totalsymmetric_irrep
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     write(output_unit, *, iostat=status)
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     write(output_unit, *, iostat=status) &
          "Total number of orbitals: ", n_ob
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     write(output_unit, *, iostat=status)
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")


     if (ssym%n_spin .eq. 1) then
        write(output_unit, *, iostat=status) &
             "This is a closed shell calculation"
     else
        write(output_unit, *, iostat=status) &
             "This is a open shell calculation"
     endif
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     write(output_unit, *, iostat=status)
     if (status .gt. 0) call error_handler( &
          "symmetry_data_write_formatted")

     if (allocated(symmetry_data_dipoles_exist)) then
        write(output_unit,'("Dipol integrals exist between the following Irreps:")',iostat=status)
        if (status .gt. 0) call error_handler( &
             "symmetry_data_write_formatted")
        write(output_unit, '("Irrep1/2  X Y Z",100(5X,I3))', iostat=status) &
             (i_ir2, i_ir2=1,ssym%n_irrep)
        do i_ir1 = 1, ssym%n_irrep
           write(output_unit,'(I3,13X,100(2X,3L2))',iostat=status) &
                i_ir1, &
                ((symmetry_data_dipoles_exist(i_ir2,i_ir1,i_xyz), i_xyz=1,3), i_ir2=1,ssym%n_irrep)
           if (status .gt. 0) call error_handler( &
                "symmetry_data_write_formatted")
        enddo
        write(output_unit, *, iostat=status)
        if (status .gt. 0) call error_handler( &
             "symmetry_data_write_formatted")
        write(output_unit, *, iostat=status)
        if (status .gt. 0) call error_handler( &
             "symmetry_data_write_formatted")
     endif
     if (allocated(symmetry_data_pdipoles_exist)) then
        write(output_unit,'("Dipol integrals exist between the following projective Irreps:")',iostat=status)
        if (status .gt. 0) call error_handler( &
             "symmetry_data_write_formatted")
        write(output_unit, '("Irrep1/2  X Y Z",100(5X,I3))', iostat=status) &
             (i_ir2, i_ir2=1,ssym%n_proj_irrep)
        do i_ir1 = 1, ssym%n_proj_irrep
           write(output_unit,'(I3,13X,100(2X,3L2))',iostat=status) &
                i_ir1, &
                ((symmetry_data_pdipoles_exist(i_ir2,i_ir1,i_xyz), i_xyz=1,3), i_ir2=1,ssym%n_proj_irrep)
           if (status .gt. 0) call error_handler( &
                "symmetry_data_write_formatted")
        enddo
        write(output_unit, *, iostat=status)
        if (status .gt. 0) call error_handler( &
             "symmetry_data_write_formatted")
        write(output_unit, *, iostat=status)
        if (status .gt. 0) call error_handler( &
             "symmetry_data_write_formatted")
     endif

   end subroutine symmetry_data_write_formatted
   !*************************************************************

   !*************************************************************
   subroutine find_totalsymmetric_irrep()
     ! purpose: determines index of totalalsymmetric irrep
     !          when other data are read in
     !** End of interface ***************************************
     integer :: i
     if ( ssym%totalsymmetric_irrep .ne. 0 ) return
     do i = 1,ssym%N_irrep
        if ( trim(adjustl(ssym%name(i))) .eq. "A1"  .or. &
             trim(adjustl(ssym%name(i))) .eq. "a1"  .or. &
             trim(adjustl(ssym%name(i))) .eq. "A1G" .or. &
             trim(adjustl(ssym%name(i))) .eq. "A1g" .or. &
             trim(adjustl(ssym%name(i))) .eq. "a1g"      ) then
           ssym%totalsymmetric_irrep = i
           return
        endif
     enddo
     ssym%totalsymmetric_irrep = 1
   end subroutine find_totalsymmetric_irrep
   !*************************************************************

   !*************************************************************
   subroutine symmetry_data_set( n_irrep, n_proj_irrep, point_group, &
        n_spin, totalsymmetric_irrep, index, dim, dim_proj, partner, partner_proj,&
        name, name_proj, pseudo, proj, jz  &
        , spin_orbit &
        )
     ! purpose: sets the private data acordind to arguments passed.
     ! if n_irrep is given, allocation is done
     ! the arguments dim, partner and name refer to irrep index
     !------------ Declaration of formal parameters -------------
     integer, intent(in), optional :: n_irrep
     integer, intent(in), optional :: n_proj_irrep
     character(len=*),    optional :: point_group
     integer, intent(in), optional :: n_spin
     integer, intent(in), optional :: totalsymmetric_irrep
     integer, intent(in), optional :: index
     logical, intent(in), optional :: proj
     ! the following arguments refer to irrep index
     integer, intent(in), optional :: dim
     integer, intent(in), optional :: dim_proj
     integer, intent(in), optional :: partner
     integer, intent(in), optional :: partner_proj
     real(RK), intent(in), optional :: jz
     character(len=*),    optional :: name
     character(len=*),    optional :: name_proj
     logical, intent(in), optional :: pseudo
     logical, intent(in), optional :: spin_orbit
     !** End of interface ***************************************
     character(len=10) :: number
     integer :: stat
     if ( present(n_irrep) ) then
        ssym%n_irrep = n_irrep
        allocate ( ssym%dim(ssym%n_irrep), &
             ssym%partner(ssym%n_irrep), &
             ssym%name(ssym%n_irrep), &
             ssym%pseudo(ssym%n_irrep), &
             stat = stat )
        if ( stat .ne. 0 ) call error_handler( &
             "symmetry_data_set: allocate failed" )

        ! SET GLOBALS:
        number_of_irreps = n_irrep
        allocate(irrep_dimensions(number_of_irreps),stat=stat)
        ASSERT(stat==0)
     endif
     if ( present(n_proj_irrep) ) then
        ssym%n_proj_irrep = n_proj_irrep
        allocate ( ssym%dim_proj(ssym%n_proj_irrep), &
             ssym%partner_proj(ssym%n_proj_irrep), &
             ssym%name_proj(ssym%n_proj_irrep), &
             ssym%jz(ssym%n_proj_irrep), &
             stat = stat )
        if ( stat .ne. 0 ) call error_handler( &
             "symmetry_data_set: allocate failed" )
     endif
     if ( present(proj) ) then
        projective = .true.
     endif
     if ( present(point_group) ) then
        ssym%point_group = adjustl(point_group) ! trailing blanks ignored by =
     endif
     if ( present(n_spin) ) then
        ssym%n_spin  = n_spin
        n_spin_set = .true.
     endif
     if ( present(totalsymmetric_irrep) ) then
        ssym%totalsymmetric_irrep  = totalsymmetric_irrep
     endif
     if ( present(index) ) then
        if (present(dim).or.present(partner).or.present(name).or.present(pseudo)) then
           if (index .lt. 1 .or. index .gt. ssym%n_irrep) then
              write(number,*) index
              call error_handler( &
                   "symmetry_data_set: wrong index "// trim(number) )
           endif
        else
           if (index .lt. 1 .or. index .gt. ssym%n_proj_irrep) then
              write(number,*) index
              call error_handler( &
                   "symmetry_data_set: wrong index "// trim(number) )
           endif
        endif
        if ( present(dim) ) then
           ssym%dim(index) = dim

           ! SET GLOBALS:
ASSERT(allocated(irrep_dimensions))
           irrep_dimensions(index) = dim
        endif
        if ( present(dim_proj) ) then
           ssym%dim_proj(index) = dim_proj
        endif
        if ( present(partner) ) then
           ssym%partner(index) = partner
        endif
        if ( present(partner_proj) ) then
           ssym%partner_proj(index) = partner_proj
        endif
        if ( present(jz) ) then
           ssym%jz(index) = jz
        endif
        if ( present(name) ) then
           ssym%name(index) = name
        endif
        if ( present(name_proj) ) then
           ssym%name_proj(index) = name_proj
        endif
        if ( present(pseudo) ) ssym%pseudo(index) = pseudo
     elseif ( present(dim) .or. present(partner) .or. present(name)&
          .or. present(pseudo) .or. present(dim_proj) .or. present(partner_proj)&
          .or. present(name_proj)) then
        call error_handler( &
             "symmetry_data_set: arguments given require index" )
     endif
     if ( present(spin_orbit) ) then
        op_spin_orbit = spin_orbit
     endif
   end subroutine symmetry_data_set

   subroutine symmetry_data_set_pcoupling(pseudo_coupling)
     integer(IK),intent(in) :: pseudo_coupling(:)
     ! *** end of interface ***

     integer :: memstat

     allocate(pcoupling_all(size(pseudo_coupling)),STAT=memstat)
     if(memstat/=0) call error_handler("sdm/symmetry_data_set_pcoupling: alloc failed")

     pcoupling_all = pseudo_coupling
   end subroutine symmetry_data_set_pcoupling

   function symmetry_data_get_pcoupling () result (pc)
     implicit none
     integer (IK) :: pc(size (pcoupling)) ! FIXME: n_irreps, but which one?
     ! *** end of interface ***

     pc = pcoupling
   end function symmetry_data_get_pcoupling

   subroutine symmetry_data_set_cccoupling(cc_coupling)
     integer(IK),intent(in) :: cc_coupling(:)
     ! *** end of interface ***

     integer :: memstat

     allocate(cccoupling_all(size(cc_coupling)),STAT=memstat)
     if(memstat/=0) call error_handler("sdm/symmetry_data_set_cccoupling: alloc failed")

     cccoupling_all = cc_coupling
   end subroutine symmetry_data_set_cccoupling

   function symmetry_data_get_cccoupling() result(ccc)
     integer(IK),pointer :: ccc(:); !<<< result
     ! *** end of interface ***

     ccc => cccoupling
   end function symmetry_data_get_cccoupling

   !*************************************************************

   !*************************************************************
   subroutine symmetry_data_elim_irreps()
     ! purpose : elliminates all irreps with dimension 0
     !** End of interface **************************************
     !------------ Declaration of local variables ---------------
     integer(IK) :: stat, n_ir_new, i_ir, ts
     type(sym) :: sym_dummy
     logical, allocatable :: dipoles_exist(:,:,:)
     logical  :: vec_only
     integer(IK) :: n_ir
     integer(IK) :: ord(20) ! (:n_irr) is used
     !------------ Executable code ------------------------------
     vec_only = .true.

     n_ir = symmetry_data_n_irreps()
     ASSERT(n_ir<=20)

     ! new number of irreps
     n_ir_new = 0
     ord = -1
     do i_ir = 1, n_ir
        if ( ssym%dim(i_ir) .gt. 0 ) then
           n_ir_new = n_ir_new + 1
           ord(n_ir_new) = i_ir
        endif
     enddo

     number_of_irreps = n_ir_new
     irrep_dimensions(:n_ir_new) = irrep_dimensions(ord(:n_ir_new))
     irrep_dimensions(n_ir_new+1:) = 0

     allocate ( sym_dummy%dim(n_ir_new), &
                sym_dummy%partner(n_ir_new), &
                sym_dummy%name(n_ir_new), &
                sym_dummy%pseudo(n_ir_new), &
                stat = stat )
     if ( stat .ne. 0 ) call error_handler( &
          "symmetry_data_elim_irreps: allocate failed" )

     ! copy to new arrays
     sym_dummy%dim(:)     = ssym%dim    (ord(:n_ir_new))
     sym_dummy%partner(:) = ssym%partner(ord(:n_ir_new))
     sym_dummy%name(:)    = ssym%name   (ord(:n_ir_new))
     sym_dummy%pseudo(:)  = ssym%pseudo (ord(:n_ir_new))

     ts = 0
     do i_ir = 1, n_ir_new
        if (ssym%totalsymmetric_irrep .eq. ord(i_ir) ) ts = i_ir
     enddo


     if (allocated(symmetry_data_dipoles_exist)) then
        allocate(dipoles_exist(n_ir_new,n_ir_new,3))

        dipoles_exist(:,:,:) = &
             symmetry_data_dipoles_exist(ord(:n_ir_new),ord(:n_ir_new),:)

        deallocate(symmetry_data_dipoles_exist)
        allocate(symmetry_data_dipoles_exist(n_ir_new,n_ir_new,3))
        symmetry_data_dipoles_exist = dipoles_exist
        deallocate(dipoles_exist)
     endif


     call free_sym(ssym,vec_only=vec_only)

     ssym%dim => sym_dummy%dim
     ssym%partner => sym_dummy%partner
     ssym%name => sym_dummy%name
     ssym%pseudo => sym_dummy%pseudo
     if ( ts .gt. 0 ) ssym%totalsymmetric_irrep = ts
     ssym%n_irrep = n_ir_new

     symmetry_data_dip_irrep_mult(:n_ir_new) = &
          symmetry_data_dip_irrep_mult(ord(:n_ir_new))

   end subroutine symmetry_data_elim_irreps
   !*************************************************************

   !*************************************************************
   subroutine symmetry_data_elim_pirreps()
     use error_module
     implicit none
     ! purpose : elliminates all projective irreps with dimension 0
     !** End of interface **************************************
     !------------ Declaration of local variables ---------------
     integer(IK) :: stat, n_ir_old, n_ir_new, n_ir_down, i_ir, i_ir_new,i_ir2, i_ir2_new
     logical, allocatable :: pdipoles_exist(:,:,:)
     type(sym) :: sym_dummy
     logical  :: proj_only
     integer             :: memstat
     integer(IK),pointer :: pforw(:),pback(:)
     ! forward & back permutations
     !------------ Executable code ------------------------------
     proj_only = .true.

     !mdf>>>
     n_ir_old = symmetry_data_n_proj_irreps()


     if(n_ir_old/=size(pcoupling_all))&
          & call error("symmetry_data_elim_pirreps: size(pcoupling_all) ???")

     allocate(pforw(n_ir_old),pback(n_ir_old),STAT=memstat)
     if(memstat/=0)call error("symmetry_data_elim_pirreps: alloc failed")
     !<<<mdf

     ! new number of irreps
     n_ir_new  = 0
     n_ir_down = n_ir_old !<<mdf
     do i_ir = 1, symmetry_data_n_proj_irreps()
        !mdf>>>
        if ( ssym%dim_proj(i_ir) .gt. 0 )then
           n_ir_new        = n_ir_new + 1
           pforw(n_ir_new) = i_ir
           pback(i_ir)     = n_ir_new
        else
           pforw(n_ir_down) = i_ir
           pback(i_ir)      = n_ir_down
           n_ir_down        = n_ir_down - 1
        endif
        !<<<mdf
     enddo

     !mdf>>>
     DPRINT '...'
     DPRINT 'symmetry_data_elim_pirreps: P/CC before >>>'
     DWRITE(*,'(10I3)') (i_ir,i_ir=1,n_ir_old)
     DWRITE(*,'(10I3)') pcoupling_all(1:n_ir_old)
     DWRITE(*,'(10I3)') cccoupling_all(1:n_ir_old)
     DPRINT '...'
     DPRINT 'symmetry_data_elim_pirreps: found',n_ir_new,'not empty irreps:'
     DWRITE(*,'(10I3)') pforw(1:n_ir_new)
     DPRINT '...'
     DPRINT 'symmetry_data_elim_pirreps: new order of irreps is:'
     DWRITE(*,'(10I3)') pforw(1:n_ir_old)

     if(all(pcoupling_all.gt.0))then
        pcoupling_all(1:n_ir_old) = pback(pcoupling_all(pforw(1:n_ir_old)))
     else
        print *,'ERROR!: not all(pcoupling_all.gt.0) (mailto:matveev@ch.tum.de)'
     endif
     pcoupling => pcoupling_all(1:n_ir_new)

     if(all(cccoupling_all.gt.0))then
        cccoupling_all(1:n_ir_old) = pback(cccoupling_all(pforw(1:n_ir_old)))
     else
        print *,'ERROR!: not all(cccoupling_all.gt.0) (mailto:matveev@ch.tum.de)'
     endif
     cccoupling => cccoupling_all(1:n_ir_new)

     DPRINT '...'
     DPRINT 'symmetry_data_elim_pirreps: P/CC after  >>>:'
     DWRITE(*,'(10I3)') (i_ir,i_ir=1,n_ir_old)
     DWRITE(*,'(10I3)') pcoupling_all(1:n_ir_old)
     DWRITE(*,'(10I3)') cccoupling_all(1:n_ir_old)
     DPRINT '...'
     DPRINT 'symmetry_data_elim_pirreps: the following section will be used:'
     DWRITE(*,'(10I3)') (i_ir,i_ir=1,n_ir_new)
     DWRITE(*,'(10I3)') pcoupling
     DWRITE(*,'(10I3)') cccoupling
     !<<<mdf

     allocate ( sym_dummy%dim_proj(n_ir_new), &
                sym_dummy%partner_proj(n_ir_new), &
                sym_dummy%name_proj(n_ir_new), &
                sym_dummy%jz(n_ir_new), &
                stat = stat )
     if ( stat .ne. 0 ) call error( &
          "symmetry_data_elim_pirreps: allocate failed" )

     ! copy to new arrays
     i_ir_new = 0
     do i_ir = 1, symmetry_data_n_proj_irreps()
        if ( ssym%dim_proj(i_ir) .gt. 0 ) then
           i_ir_new = i_ir_new + 1
           sym_dummy%dim_proj(i_ir_new) = ssym%dim_proj(i_ir)
           sym_dummy%partner_proj(i_ir_new) = ssym%partner_proj(i_ir)
           sym_dummy%name_proj(i_ir_new) = ssym%name_proj(i_ir)
           sym_dummy%jz(i_ir_new) = ssym%jz(i_ir)
        endif
     enddo

     if (allocated(symmetry_data_pdipoles_exist)) then
        allocate(pdipoles_exist(n_ir_new,n_ir_new,3))
        i_ir_new = 0
        do i_ir = 1, symmetry_data_n_proj_irreps()
           if ( ssym%dim_proj(i_ir) .gt. 0 ) then
              i_ir_new = i_ir_new + 1
              i_ir2_new = 0
              do i_ir2 = 1, symmetry_data_n_proj_irreps()
                 if ( ssym%dim_proj(i_ir2) .gt. 0 ) then
                    i_ir2_new = i_ir2_new + 1
                    pdipoles_exist(i_ir2_new,i_ir_new,:) = &
                         symmetry_data_pdipoles_exist(i_ir2,i_ir,:)
                 endif
              enddo
           endif
        enddo
        deallocate(symmetry_data_pdipoles_exist)
        allocate(symmetry_data_pdipoles_exist(n_ir_new,n_ir_new,3))
        symmetry_data_dipoles_exist = pdipoles_exist
        deallocate(pdipoles_exist)
     endif

     if (allocated(symmetry_data_pdipoles_exist)) then
        allocate(pdipoles_exist(n_ir_new,n_ir_new,3))
        i_ir_new = 0
        do i_ir = 1, symmetry_data_n_proj_irreps()
           if ( ssym%dim_proj(i_ir) .gt. 0 ) then
              i_ir_new = i_ir_new + 1
              i_ir2_new = 0
              do i_ir2 = 1, symmetry_data_n_proj_irreps()
                 if ( ssym%dim_proj(i_ir2) .gt. 0 ) then
                    i_ir2_new = i_ir2_new + 1
                    pdipoles_exist(i_ir2_new,i_ir_new,:) = &
                         symmetry_data_pdipoles_exist(i_ir2,i_ir,:)
                 endif
              enddo
           endif
        enddo
        deallocate(symmetry_data_pdipoles_exist)
        allocate(symmetry_data_pdipoles_exist(n_ir_new,n_ir_new,3))
        symmetry_data_dipoles_exist = pdipoles_exist
        deallocate(pdipoles_exist)
     endif

     call free_sym(ssym,proj_only=proj_only)

     ssym%dim_proj => sym_dummy%dim_proj
     ssym%partner_proj => sym_dummy%partner_proj
     ssym%name_proj => sym_dummy%name_proj
     ssym%jz => sym_dummy%jz
     ssym%n_proj_irrep = n_ir_new

     !mdf>>>
     deallocate(pforw,pback,STAT=memstat)
     if(memstat/=0)call error("symmetry_data_elim_pirreps: dealloc failed")
     !<<<mdf
   end subroutine symmetry_data_elim_pirreps
   !*************************************************************

!--------------- End of module ----------------------------------
end module symmetry_data_module
