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
!==============================================================a
module coordinates_module
  !---------------------------------------------------------------
  !  Purpose: routines to convert coordinates from
  !          cartesian to internal. The coordinates themselves
  !          are also kept in this module.
  !
  !  Module called by: main_opt
  !
  !  References: Wilson-Buch
  !
  !----------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use math_module
  use opt_data_module
  use coortype_module
  use atom_data_module
  use geo_operations_module
  use iounitadmin_module
  use allocopt_module
#ifdef WITH_EFP
  use qmmm_interface_module, only: efp
  use efp_module, only: n_efp, efp_fixed, qm_fixed
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private

  !------------ Declaration of constants and variables ------------
  integer(kind=i4_kind),pointer,public :: sym_type(:)
  ! Wilsons B-Matrix
  real(kind=r8_kind),allocatable,public    :: bmat(:,:),bmat_prim(:,:)
  real(kind=r8_kind),allocatable,public    :: gmat_sqr(:,:),gmat_sqr_inv(:,:)
  ! general right inverse of Wilsons B-Matrix
  real(kind=r8_kind),allocatable,public    :: bmat_inv(:,:)

#if 1 /* mass_weighted */
  real(kind=r8_kind),allocatable,public    :: gmat(:,:)
#endif

  ! transformation matrix from primitive internals to delocalized internals
  real(kind=r8_kind),allocatable,public    :: umat(:,:),umat_trans(:,:)
  ! matrix used for constrained optimization in Z-Matrix coordinates
  ! maybe this is a waste of memory, but it simplifies the steps
  ! and together with 'ParaGauss' memory is not really a concern
  ! for the geometry optimizer

  real(kind=r8_kind),pointer,public    :: q_old(:) !(n_internal)
  ! holds values of internal coordinates

  real(kind=r8_kind),allocatable,public    :: tmat(:,:)
  ! reduction and expand-matrices to go from symmetry-reduced coordinates
  ! to the full set of internal coordinates and back. Just for the
  ! case zmat_coordinates.and.zmat_format

  real(kind=r8_kind),allocatable,public    :: reduc_mat(:,:),expand_mat(:,:)
  ! constraint-vectors used for constrained optimzation in delocalized
  ! coordinates
  real(kind=r8_kind),allocatable,public    :: constraint(:,:)
  integer(kind=i4_kind),public             :: n_constraint


  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public coordinates_setup,cart_to_internal,generate_bmat, &
       invert_bmat, fake_internals,internal_to_cart,sym_check,cart2cart, &
       print_internal,set_internal,val_types,alloc_bmat,free_bmat
  public dealloc_intcoor,dealloc_reduce,dealloc_constraint,dealloc_cartcoor
  !------------ Subroutines ---------------------------------------
contains
  !*************************************************************

  subroutine coordinates_setup(geo_loop)
    ! Purpose: main routine for the projection of cartesian
    !          coordinates to internal coordinates. At the
    !          same time the Matrix B and its inverse are
    !          generated.
    !          The variables n_internal and n_primitive have been
    !           set as follows:
    !           a) zmat_coordinates
    !              -zmat_format:
    !               n_internal = maxval(numx)
    !               n_primitive = 3*N-6
    !              -free_format:
    !               n_internal = 3*N-6
    !               n_primitive user specified
    !           b) delocalized_coordinates
    !              -zmat_format:
    !               n_internal = 3N-6 = n_primitive
    !              -free_format:
    !               n_internal = 3*N-6
    !               n_primitive user specified
    !              -valence_format (= 'blanket approach')
    !               n_internal = 3N-6 - N_CONSTRAINT
    !               n_primitive determined in 'val_types'
    ! Subroutine called by: main_opt
    !----------------------------------------------------
!  use allocopt_module
    ! --- Declaration of formal parameters --------------
    integer(kind=i4_kind), optional :: geo_loop
    ! --- Declaration of local variables ----------------
    integer(kind=i4_kind) :: alloc_stat,i
!!$! These variables can be used to test the transformation properties
!!$! of the B-Matrix
!!$    real(kind=r8_kind)             :: summ
!!$    real(kind=r8_kind)             :: summ_vec(3)
    ! --- Executable code -------------------------------

    ! find out the number of possible primitive coordinates

    call set_cartesian(atom)
    if(tsscan_sphere) call set_cartesian(atom_reactant,xyz_reactant)
    if(tsscan_sphere.and.exist_product) call set_cartesian(atom_product,xyz_product)

    if (delocalized_coordinates.and.valence_format) then
       call val_types(n_coor=n_primitive,deloc=.true., &
            calc_dimension=.true.,keep_coor=keep_valence)
    endif

    if(cart_format) then
       call alloc_cartcoor()
    else
       call alloc_intcoor(s_prim,q_prim,s,sym_type)
       !    and its inverse for n_internal coordinates
       if(tsscan_sphere) call alloc_intcoor(s_prim_reactant,q_prim_reactant,s_reactant)
       if(tsscan_sphere) call alloc_intcoor(s_prim_pointonmep,q_prim_pointonmep,s_pointonmep)
       if(qst_step.or.(tsscan_sphere.and.exist_product)) call alloc_intcoor(s=s_product)
       if(tsscan_sphere.and.exist_product) call alloc_intcoor(s_prim_product,q_prim_product)
    end if

    deloc: if (delocalized_coordinates) then

       if (zmat_format) then
          call zmat_types(s_prim,n_primitive)
       elseif (free_format) then
          call freeformat_types(s_prim)
       elseif (valence_format) then
          call val_types(coor=s_prim,n_coor=n_primitive,deloc=.true.,&
               calc_dimension=.false.,keep_coor=keep_valence)
       endif

       call  set_internal(s_prim(:)%value,n_primitive,s_prim,atom) !(1)

       if (print_debug) then
          write(OPT_STDOUT,*)" Underlying set of primitive internals is :"
          call print_internal(s_prim,n_primitive,full=.true.)
       endif

       call constraint_setup(s_prim,n_primitive) !(1) delocalized

       call alloc_bmat(n_internal) !this allocates the B-Matrix, its Transpose (1)
       allocate(s(n_internal),STAT=alloc_stat)
       if (alloc_stat/=0 ) call error_handler("coordinates_setup: allocation s failed")
       s(:)%value = zero
       s(:)%typ = 0_i4_kind

    else   deloc
       if (zmat_format) then
          call zmat_types(s_prim,n_primitive)
          if(tsscan_sphere) call zmat_types(s_prim_reactant,n_primitive)
          if(tsscan_sphere.and.exist_product) call zmat_types(s_prim_product,n_primitive)
          call alloc_bmat(n_primitive)  !(2)
       elseif (free_format) then
          call freeformat_types(s)
          call alloc_bmat(n_internal)  !(3)
       elseif (cart_format) then
       endif
    endif deloc

    allocate(q_old(n_internal),q(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    ! Do extra allocations for the case of delocalized coordinates
    if (delocalized_coordinates) then
       allocate(umat(n_primitive,n_internal),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            (" coordinates_setup : allocation of umat failed")
       allocate(umat_trans(n_internal,n_primitive),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            (" coordinates_setup : allocation of umat failed")
       umat=zero
       umat_trans=zero
    endif

    ! comment this out only if you doubt that the B-Matrix (in primitive internals)
    ! is correct, FN 2/98
    !    call check_bmat(bmat)

    if (zmat_coordinates) then
       if (zmat_format) then
          call generate_bmat(s_prim,n_primitive,bmat)
          call set_internal(q_prim,n_primitive,s_prim,atom) !(2)
          if(tsscan_sphere) then
           if(exist_product) call set_internal(q_prim_product,n_primitive,s_prim_product,atom_product)
           call set_internal(q_prim_reactant,n_primitive,s_prim_reactant,atom_reactant)
          endif

          if(present(geo_loop)) then
             call reduce_setup(geo_loop)
          else
             call reduce_setup()
          end if

          if(tsscan_sphere) then
           call reduce(s_reactant,s_prim_reactant,q,q_prim_reactant)
           if(exist_product) call reduce(s_product,s_prim_product,q,q_prim_product)
          endif
          call reduce(s,s_prim,q,q_prim)
          if(tsscan_sphere.and.exist_product) then
           sphere_dependent_var= &
               dependent_var(n_internal,s,s_reactant,distance_to_reactant,distance_to_product)
          elseif(tsscan_sphere) then
           sphere_dependent_var=   dependent_var(n_internal,s,s_reactant,distance_to_reactant)
          endif

       else
          call generate_bmat(s,n_internal,bmat)
          call set_internal(q,n_internal,s,atom) !(3)
       endif

       call constraint_setup(s,n_internal) !(2)

    else if( cart_coordinates) then
       call set_cart_for_opt(q_prim,n_primitive)
       call cart_constraint_setup(n_internal) !?????????????
       q=q_prim(1:n_internal)
#ifdef WITH_EFP
       if(qm_fixed) q=q_prim(3*n_atoms+1:3*n_atoms+n_internal)
#endif
    elseif ( delocalized_coordinates) then
       DPRINT 'generate_bmat'
       call generate_bmat(s_prim(1:n_primitive),n_primitive,bmat_prim)
       DPRINT 'generate_bmat done'
       ! first test if equation (1) is valid for both B_prim and B

!      call generate_bmat_reduced(umat_flag=.true.)

       ! generate the reduction matrix:
       call generate_umat(bmat_prim,umat)
       umat_trans = transpose(umat)

       ! reduce the bmat_prim into bmat:
       bmat = matmul(umat_trans,bmat_prim)

!!$! Uncomment these lines to test the transformation properties of
!!$! the B-Matrix, i.e. conservation of total momentum and total angular
!!$! momentum for the update procedure.
!!$       do i=1,n_primitive
!!$          summ=zero
!!$          do j=1,3*(n_atoms+n_dummy)
!!$             summ = summ + bmat_prim(i,j)
!!$          enddo
!!$          print*,summ
!!$          summ_vec=zero
!!$          start=1
!!$          do j=1,n_atoms+n_dummy
!!$             summ_vec =summ_vec + cross_product(atom(j)%x,bmat_prim(i,start:start+2))
!!$             start=start+3
!!$          enddo
!!$          print*,summ_vec
!!$       enddo
!!$       do i=1,n_internal
!!$          summ=zero
!!$          do j=1,3*(n_atoms+n_dummy)
!!$             summ = summ + bmat(i,j)
!!$          enddo
!!$          print*,summ
!!$          summ_vec=zero
!!$          start=1
!!$          do j=1,n_atoms+n_dummy
!!$             summ_vec =summ_vec + cross_product(atom(j)%x,bmat(i,start:start+2))
!!$             start=start+3
!!$          enddo
!!$          print*,summ_vec
!!$       enddo

       ! generate a set of delocalized coordinates
       q_prim = s_prim(1:n_primitive)%value
       q = matmul(umat_trans,q_prim)

       if (print_debug) then
          write(OPT_STDOUT,*)" Delocalized Internals at Start :"
          do i=1,n_internal
             write(OPT_STDOUT,'(i4,3x,f9.5)')i,q(i)
          enddo
          write(OPT_STDOUT,*)" "
       endif
    endif

#if 0 /* mass_weighted */
    allocate(gmat_sqr(ubound(bmat,1),ubound(bmat,1)), &
             gmat_sqr_inv(ubound(bmat,1),ubound(bmat,1)),stat=allocopt_stat(11))
    ASSERT(allocopt_stat(11).eq.0)
#endif

    if(.not. cart_format) then
       call invert_bmat(bmat,bmat_inv,.false.) ! (2)

       if (print_bmat) then
          write(OPT_STDOUT,*)" B-Matrix at Start:"
          call print_matrix(bmat,ubound(bmat,1),3*(n_atoms+n_dummy),10_i4_kind)
          write(OPT_STDOUT,*)" Generalized Inverse of B-Matrix at Start "
          call print_matrix(bmat_inv,3*(n_atoms+n_dummy),ubound(bmat,1),10_i4_kind)
          write(OPT_STDOUT,*)" "
       endif
    end if

#if 0 /* mass_weighted */
    allocate(gmat(n_primitive,n_primitive),stat=allocopt_stat(8))
    ASSERT(allocopt_stat(8).eq.0)
    call generate_gmat(gmat)
    deallocate(gmat,stat=allocopt_stat(8))
    ASSERT(allocopt_stat(8).eq.0)
    allocopt_stat(8)=1
#endif

    q_old=q

  end subroutine coordinates_setup

  subroutine set_cartesian(atom,xyz)
    ! Purpose: allocate and set the variable 'atom' by reading
    !          it from x,y,z. These in turn may be read in
    !          in different formats and thus the routine
    !          is decoupled from the actual reading in of
    !          cartesian coordinates.
    !          At the same time the variable containing the internal
    !          Z-matrix-type coordinates is allocated and
    !          initialized.
    !
    ! subroutine called by: 'coordinates_setup'
    ! ---------------------------------------------------------
    integer(kind=i4_kind)   :: alloc_stat,i,n_frag_atoms,j
    type(atom_type), pointer   :: atom(:)
    real(kind=r8_kind),intent(in),optional :: xyz(:,:)
#ifdef WITH_EFP
    integer(kind=i4_kind)   :: k
#endif

    n_frag_atoms=0
#ifdef WITH_EFP
    if(efp .and. n_efp > 0) n_frag_atoms=3*n_efp
#endif
    allocate(atom(n_atoms+n_dummy+n_frag_atoms),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

   if(present(xyz)) then
    do i=1,n_atoms+n_dummy
       atom(i)%x(1) = real(xyz(1,i),kind=r8_kind)
       atom(i)%x_old(1) = real(xyz(1,i),kind=r8_kind)
       atom(i)%x(2) = real(xyz(2,i),kind=r8_kind)
       atom(i)%x_old(2) = real(xyz(2,i),kind=r8_kind)
       atom(i)%x(3) = real(xyz(3,i),kind=r8_kind)
       atom(i)%x_old(3) = real(xyz(3,i),kind=r8_kind)
       atom(i)%charge = charge(i)
       atom(i)%dummy=.false.
       atom(i)%heavy=.false.
       if(charge(i)==dummy_charge) atom(i)%dummy=.true.
    enddo
   else
    do i=1,n_atoms+n_dummy+n_frag_atoms
       atom(i)%x(1) = real(x(i),kind=r8_kind)
       atom(i)%x_old(1) = real(x(i),kind=r8_kind)
       atom(i)%x(2) = real(y(i),kind=r8_kind)
       atom(i)%x_old(2) = real(y(i),kind=r8_kind)
       atom(i)%x(3) = real(z(i),kind=r8_kind)
       atom(i)%x_old(3) = real(z(i),kind=r8_kind)
       atom(i)%charge = charge(i)
       atom(i)%dummy=.false.
       atom(i)%heavy=.false.
       if(charge(i)==dummy_charge) atom(i)%dummy=.true.
    enddo
    endif

    do i=1,n_tot_atoms
       atom(i)%mass=nuc_mass(int(charge(i),i4_kind))
    end do
    do i=1,n_infinite_masses
       if(atom(infinite_masses_index(i))%dummy) &
          call error_handler('set_cartesian: you try to set a infinfinite mass to a dummy atom')
       atom(infinite_masses_index(i))%heavy=.true.
    end do
    do i=1,n_defined_masses
       if(atom(defined_masses_index(i))%dummy) &
            call error_handler('set_cartesian: you try to set a mass to a dummy atom')
       if(atom(defined_masses_index(i))%heavy) &
            call error_handler('set_cartesian: you try to set a mass to a infinite heavy atom')
       write(OPT_STDOUT,*) 'setting mass of atom number',i,'to',defined_masses(i),'au'
       atom(defined_masses_index(i))%mass = defined_masses(i)
    end do

    ! output
    if (.not.atomic_units) then
       write(OPT_STDOUT,*)
       write(OPT_STDOUT,*)" --- Cartesian coordinates (Angstroms) ---"
       j=0
       if(n_atoms+n_dummy > 0) then
          write(OPT_STDOUT,'("Atom",3x,"Element",5x,"x      ",3x,"y      ",3x,"z       ")')
          do i=1,n_atoms+n_dummy
             write(OPT_STDOUT,1103)i,adjustl(trim(symbol(nint(atom(i)%charge)))),&
                  atom(i)%x(1),&
                  atom(i)%x(2),&
                  atom(i)%x(3)
          enddo
       end if
#ifdef WITH_EFP
       if(efp .and. n_efp > 0) then
          write(OPT_STDOUT,'("Frag",3x,"Element",5x,"x           ",3x,"y           ",3x,"z            ")')
          do i=1,n_efp
             do k=1,3
                j=j+1
                write(OPT_STDOUT,1103)i+n_atoms+n_dummy,adjustl(trim(symbol(nint(atom(j)%charge)))),&
                     atom(j)%x(1),&
                     atom(j)%x(2),&
                     atom(j)%x(3)
             end do
          end do
       end if
#endif
    else
       write(OPT_STDOUT,*)
       write(OPT_STDOUT,*)" --- Cartesian coordinates (au) -----------"
       j=0
       if(n_atoms+n_dummy > 0) then
          write(OPT_STDOUT,'("Atom",3x,"Element",5x,"x           ",3x,"y           ",3x,"z            ")')
          do i=1,n_atoms+n_dummy
             j=j+1
             write(OPT_STDOUT,1103)i,adjustl(trim(symbol(nint(atom(i)%charge)))),&
                  atom(i)%x(1),&
                  atom(i)%x(2),&
                  atom(i)%x(3)
          enddo
       end if
#ifdef WITH_EFP
       if(efp .and. n_efp > 0) then
          write(OPT_STDOUT,'("Frag",3x,"Element",5x,"x           ",3x,"y           ",3x,"z            ")')
          do i=1,n_efp
             do k=1,3
                j=j+1
                write(OPT_STDOUT,1103)i+n_atoms+n_dummy,adjustl(trim(symbol(nint(atom(j)%charge)))),&
                     atom(j)%x(1),&
                     atom(j)%x(2),&
                     atom(j)%x(3)
             end do
          end do
       end if
#endif
       write(OPT_STDOUT,*)
    endif
1103 format(i3,6x,a3,4x,f12.8,3x,f12.8,3x,f12.8)


!1101 format(3(2x,f13.7))
  end subroutine set_cartesian


  subroutine zmat_types(coor,n_coor)
    ! Purpose : set the type of internal coordinates to
    !           one of:
    !           bond_length
    !           bond_angle
    !           dihedral_angle
    !           This routine is based on the 'gxfile'-Option
    !           where a gxfile according to Vladimirs
    !           format from 'simol' is read in.
    ! ----------------------------------------------
    type(int_coor),intent(inout)     :: coor(:)
    integer(kind=i4_kind),intent(in) :: n_coor
    ! --- Declaration of local variables -----------
    integer(kind=i4_kind)  :: num_coord,n,i,counter,j,&
         l_eq, sign_flag(n_coor)
    sign_flag=1_i4_kind ! sign flag is used if you want to build
    ! antisymmetric modes
    ! if in the second deck two symmetrieaequivalent coordinates are marked
    ! with different signs, the resulting mode will be built according
    ! to the signs
    num_coord=0_i4_kind
    counter=0_i4_kind
    !  zmat = 1. deck of gxfile-input
    !  numx = 2. deck of gxfile-input
    if(tsscan_sphere) select_sphere_vars=minval(zmat).lt.0
    comp_1: do i=1,3
       atoms: do n=1,n_atoms+n_dummy
          if (zmat(i,n).eq.0) then
             ! if coordinate is omitted , skip this
             cycle atoms

          elseif (numx(i,n).eq.0) then
             !numx /= 0 indicates that this coordinate is to
             !be optimized. The value of numx refers to the type of
             !symmetry equivalent internal coordinate.
             num_coord=num_coord+1
             coor(num_coord)%var = .false.
             sym_type(num_coord) = 0
          elseif(numx(i,n) /= 0 ) then
             if(numx(i,n)<0) then
                sign_flag(num_coord+1)=-1
             end if
             num_coord=num_coord+1
             sym_type(num_coord) = abs(numx(i,n))
             coor(num_coord)%var = .true.
          endif
          ! the sym_type array contains the symmetry-type of the
          ! internal coordinate with index num_coord

          if (zmat_coordinates) then
          coor(num_coord)%sphere=zmat(i,n).lt.0
             if (sym_type(num_coord)==0 ) then
                coor(num_coord)%unique = .false.
                coor(num_coord)%var = .false.
             else
                coor(num_coord)%unique=.true.
                do j=1,num_coord-1
                   if (sym_type(num_coord)==sym_type(j) ) then
                      ! this means that a coordinate of this
                      ! symmetry type already exists
                      coor(num_coord)%unique=.false.
                      exit
                   endif
                enddo
             endif
          endif

          if (i==1) then
             coor(num_coord)%typ = b_length
             coor(num_coord)%length%partner1 = n
             coor(num_coord)%length%partner2 = abs(zmat(i,n))
          elseif( i==2) then
             coor(num_coord)%typ = b_angle
             coor(num_coord)%angle%partner1 = n
             coor(num_coord)%angle%partner2 = abs(zmat(i,n))
             coor(num_coord)%angle%apex = abs(zmat(1,n))
          elseif (i==3) then
             coor(num_coord)%typ = d_angle
             coor(num_coord)%dihedral%partner1 = n !8
             coor(num_coord)%dihedral%partner2 = abs(zmat(i,n)) !2
             coor(num_coord)%dihedral%base1 = abs(zmat(1,n)) !3
             coor(num_coord)%dihedral%base2 = abs(zmat(2,n))
          else
             call error_handler("zmat_types: something went totally wrong.Check gxfile")
          endif

       enddo atoms
    enddo comp_1

    do i=1,n_coor
       if (coor(i)%unique) then
          l_eq=0
          do j=1,n_coor
             if(sym_type(j)==sym_type(i)) then
                l_eq=l_eq+1
                coor(i)%equal(l_eq) = j
                coor(i)%equal_sign(j-i+1) = sign_flag(j)
                if (i==j) cycle
                coor(i)%n_equal = coor(i)%n_equal + 1
             endif
          enddo
       endif
    enddo
  end subroutine zmat_types


  subroutine freeformat_types(coor)
    ! Purpose: set the types for the internal coordinate
    !          variable 's'. As a preliminary maesure the
    !          component 'unique' is set to true for
    !          all internal variables at the moment.
    !
    ! subroutine called by: coordinates_setup
    ! ---------------------------------------------------------
    type(int_coor),intent(inout)   :: coor(:)
    integer(kind=i4_kind) :: num_coord,i


    do num_coord=1,n_primitive
       coor(num_coord)%var = .true. ! no constraints
       coor(num_coord)%unique = .true. ! no symmetry conservation
    enddo

    num_coord=1
    do i=1,n_stretch
       coor(num_coord)%typ = b_length
       coor(num_coord)%length%partner1 = stretch(i,1)
       coor(num_coord)%length%partner2 = stretch(i,2)
       if (.not.var_stretch(i)) then
          coor(num_coord)%var = .false.
       endif
       num_coord=num_coord+1
    enddo

    do i=1,n_bend
       coor(num_coord)%typ = b_angle
       coor(num_coord)%angle%partner1 = bend(i,1)
       coor(num_coord)%angle%partner2 = bend(i,2)
       coor(num_coord)%angle%apex     = bend(i,3)
       if (.not.var_bend(i)) then
          coor(num_coord)%var = .false.
       endif
       num_coord=num_coord+1
    enddo

    do i=1,n_torsion
       coor(num_coord)%typ = d_angle
       coor(num_coord)%dihedral%partner1 = torsion(i,1)
       coor(num_coord)%dihedral%partner2 = torsion(i,2)
       coor(num_coord)%dihedral%base1 = torsion(i,3)
       coor(num_coord)%dihedral%base2 = torsion(i,4)
       if (.not.var_tors(i)) then
          coor(num_coord)%var = .false.
       endif
       num_coord=num_coord+1
    enddo

  end subroutine freeformat_types
   !*************************************************************
  subroutine val_types(coor,n_coor,deloc,calc_dimension,keep_coor)
    ! Purpose: compute valence coordinates of the current molecule.
    !          Attention: in the current implementation, dummy
    !          atoms are simply SKIPPED. This is because the valence
    !          coordinates are currently only used for establishing
    !          an initial Hessian using force constants and
    !          empirical rules (see also valence_coord_module).
    !
    ! 26.5.98: implement a further extension that will save the
    !          valence coordinates to file after the first run.
    !          On subsequent calls (or subsequent geometry loops)
    !          the routine will check if the file is present. If so,
    !          the coordinates defined in this file will be used as
    !          the underlying primitives to the delocalized
    !          coordinates.
    !-----------------------------------------------------------
    use constants, only: angstrom
    implicit none
    type(int_coor),optional                     :: coor(:)
    integer(kind=i4_kind),optional              :: n_coor
    logical,optional,intent(in)                 :: deloc
    logical,optional,intent(in)                 :: calc_dimension
    logical,optional,intent(in)                 :: keep_coor
    ! --- Declaration of local variables -----------------------
    logical                :: neighbour(n_atoms+n_dummy,n_atoms+n_dummy)
    integer(kind=i4_kind)  :: n_tot, counter, i, j, k, l, &
         p1,b1,b2,p2,i_perm,n_cons,summi,summo,i_cons,io_val
    integer(kind=i4_kind)  :: array(4),arr(4)
    real(kind=r8_kind),dimension(3)  :: e12, e23, e34, e31, e32
    real(kind=r8_kind)               :: c_phi, fac
    real(kind=r8_kind)               :: rij, max_length
    logical                          :: delocalized,calc,keep,exist
    real(kind=r8_kind),parameter     :: close_to_pi = 0.01_r8_kind
    ! --- executable code -------------------------------------
    io_val=60
    if (present(deloc)) then
       if (deloc) then
          delocalized=.true.
       else
          delocalized=.false.
       endif
    else
       delocalized=.false.
    endif
    if (present(calc_dimension)) then
       if (calc_dimension) then
          calc=.true.
       else
          calc=.false.
       endif
       if(calc.and.present(coor)) &
          call error_handler &
          ("val_types: EITHER use this routine to calculate&
          & n_primitive OR setup the coordinates")
    else
       calc=.false.
    endif
    if (present(keep_coor)) then
       if (keep_coor) then
          keep=.true.
       else
          keep=.false.
       endif
    else
       keep=.false.
    endif

    if (keep) then
       exist=.false.
       inquire(EXIST=exist,FILE=trim(opt_data_dir)//'/val_coor.dat')
       if (exist) then
          write(OPT_STDOUT,*)" "
          write(OPT_STDOUT,*)"val_types: reading valence coordinates from file"
          write(OPT_STDOUT,*)" "
          call read_val()
          return
       else
          write(OPT_STDOUT,*)" "
          write(OPT_STDOUT,*)"val_types: setting up valence coordinates "
          write(OPT_STDOUT,*)" "
       endif
    endif


    counter = 0
    if (delocalized) then
       ! determine the number and type of constrained
       ! variables
       n_cons = ubound(const_stretch,1) + ubound(const_bend,1) + &
            ubound(const_torsion,1)
       if (calc) then
          counter=counter+n_cons
       else
          do i=1,ubound(const_stretch,1)
             counter=counter+1
             coor(counter)%typ = b_length
             coor(counter)%length%partner1 = const_stretch(i,1)
             coor(counter)%length%partner2 = const_stretch(i,2)
             coor(counter)%unique = .true.
             coor(counter)%var = .false.
             coor(counter)%n_equal = 1
          enddo
          do i=1,ubound(const_bend,1)
             counter=counter+1
             coor(counter)%typ = b_angle
             coor(counter)%angle%partner1 = const_bend(i,1)
             coor(counter)%angle%partner2 = const_bend(i,2)
             coor(counter)%angle%apex = const_bend(i,3)
             coor(counter)%unique = .true.
             coor(counter)%n_equal = 1
             coor(counter)%var = .false.
          enddo
          do i=1,ubound(const_torsion,1)
             counter=counter+1
             coor(counter)%typ = d_angle
             coor(counter)%dihedral%partner1 = const_torsion(i,1)
             coor(counter)%dihedral%partner2 = const_torsion(i,2)
             coor(counter)%dihedral%base1 = const_torsion(i,3)
             coor(counter)%dihedral%base2 = const_torsion(i,4)
             coor(counter)%unique = .true.
             coor(counter)%n_equal = 1
             coor(counter)%var = .false.
          enddo
       endif
    endif


    ! first establish the neighbour matrix thereby setting up the
    ! variables of type 'bond_length'
    neighbour = .false.
    n_tot=n_atoms+n_dummy

    do i=1,n_tot
       if (.not.delocalized) then
          if (atom(i)%charge == dummy_charge) cycle
       endif
       stretch_j: do j=i+1,n_tot

          if (.not.delocalized) then
             if (atom(j)%charge == dummy_charge ) cycle
          endif
          rij = abs_value(atom(i)%x - atom(j)%x)
          if (atomic_units) then
             fac = alpha_valence * angstrom
          else
             fac = alpha_valence
          endif
          max_length = (covrad(nint(atom(i)%charge)) + &
               covrad(nint(atom(j)%charge))) * fac
          ! the above values are taken from T.A.Halgrens program 'NABOR' -> check
          if (rij <= max_length) then
             neighbour(i,j) = .true.
             neighbour(j,i) = .true.
             if (delocalized) then
                summi = i+j
                c_stretch: do i_cons=1,ubound(const_stretch,1)
                   summo = sum(const_stretch(i_cons,:))
                   if (summo/=summi) cycle c_stretch
                   if (i==const_stretch(i_cons,1).or.&
                        i==const_stretch(i_cons,2)) then
                      cycle stretch_j
                   endif
                enddo c_stretch
             endif
             counter= counter+1
             if(.not.calc) then
                coor(counter)%typ = b_length
                coor(counter)%length%partner1 = i
                coor(counter)%length%partner2 = j
                coor(counter)%unique=.true.
                coor(counter)%var = .true.
                coor(counter)%n_equal=1
             endif
          endif

       enddo stretch_j
    enddo

    if (print_debug) then
       do i=1,n_tot
          print*,(neighbour(i,j),j=1,n_tot)
       enddo
    endif
    ! now 'bond_angles'
    do i=1,n_tot
       angle_j: do j=1,n_tot

          if (i==j) cycle angle_j
          if (.not.neighbour(i,j)) cycle angle_j

          angle_k: do k=j+1,n_tot
             if (k==i.or.k==j) cycle angle_k
             if (.not.neighbour(i,k)) cycle angle_k
             if (delocalized) then
                ! if this coordinate is among the constraint variables
                ! skip it
                summi=i+j+k
                c_bend: do i_cons=1,ubound(const_bend,1)
                   summo=sum(const_bend(i_cons,:))
                   if (summi/=summo) cycle c_bend
                   if (i/=const_bend(i_cons,3)) cycle c_bend
                   if (j==const_bend(i_cons,1).or.&
                        j==const_bend(i_cons,2)) then
                      cycle angle_k
                   endif
                enddo c_bend
             endif

             ! first see if this bend coordinate is /= 180 degrees.
             ! If so, throw it away. In the Paper by Baker, it is
             ! replaced by sth. else. This is only required for
             ! really pathological cases.
             e31 = atom(j)%x - atom(i)%x
             e31 = e31/abs_value(e31)
             e32 = atom(k)%x - atom(i)%x
             e32 = e32/abs_value(e32)
             c_phi = dot_product(e31,e32)
             if ( abs(abs(c_phi)-one)<=close_to_pi) then ! this angle is near pi
                if (print_debug) then
                   write(OPT_STDOUT,*)"val_types: the angle formed by atoms :",j,i,k
                   write(OPT_STDOUT,*)"           is close to 180/0 degrees (pi/0) :",acos(c_phi)*convert1," degrees"
                   write(OPT_STDOUT,*)"           This coordinate will be skipped "
                   write(OPT_STDOUT,*)" If you think this angle should be included anyway "
                   write(OPT_STDOUT,*)" change the variable CLOSE_TO_PI in routine val_types. "
                   write(OPT_STDOUT,*)" "
                endif
                cycle angle_k
             endif
             counter = counter+1
             if (.not.calc) then
                coor(counter)%typ=b_angle
                coor(counter)%angle%apex = i
                coor(counter)%angle%partner1 = j
                coor(counter)%angle%partner2 = k
                coor(counter)%unique=.true.
                coor(counter)%var=.true.
                coor(counter)%n_equal=1
             endif
          enddo angle_k
       enddo angle_j
    enddo
    ! now dihedral angles
    do i=1,n_tot-3
       dihedral_j: do j=i+1,n_tot-2
          if (i==j) cycle dihedral_j

          dihedral_k: do k=j+1,n_tot-1

             if (k==j.or.k==i) cycle dihedral_k
             dihedral_l: do l=k+1,n_tot
                ! decide if the three vectors defining the dihedral
                ! angle are pairwise colinear. If so, skip this coordinate
                ! because the dihedral angle is not defined for this
                ! situation.
                arr = (/i,j,k,l/)
                array=0_i4_kind
                permutation: do i_perm=1,12
                   call perm(i_perm,arr,array)
                   p1=array(1)
                   b1=array(2)
                   b2=array(3)
                   p2=array(4)
                   if (.not.neighbour(p1,b1)) then
                      cycle permutation
                   endif
                   if (.not.neighbour(b1,b2)) then
                      cycle permutation
                   endif
                   if (.not.neighbour(b2,p2)) then
                      cycle permutation
                   endif

                   if (delocalized) then
                      ! now find out if this dihedral angle is already among
                      ! the constraint variables. If so, skip it.
                      summi=i+j+k+l
                      c_torsion: do i_cons=1,ubound(const_torsion,1)
                         summo=sum(const_torsion(i_cons,:))
                         if (summi/=summo) cycle c_torsion
                         summi=b1+b2
                         summo=const_torsion(i_cons,3)+const_torsion(i_cons,4)
                         if (summi/=summo) cycle c_torsion
                         summi=p1+p2
                         summo=const_torsion(i_cons,1)+const_torsion(i_cons,2)
                         if (summi/=summo) cycle c_torsion
                         if ( (p1==const_torsion(i_cons,1).or.p1==const_torsion(i_cons,2)) .and. &
                              (b1==const_torsion(i_cons,3).or.b1==const_torsion(i_cons,4)) ) then
                            cycle dihedral_l
                         endif
                      enddo c_torsion
                   endif

                   e12 = atom(b1)%x - atom(p1)%x
                   e12 = e12/abs_value(e12)

                   e23 = atom(b2)%x - atom(b1)%x
                   e23 = e23/abs_value(e23)

                   e34 = atom(p2)%x - atom(b2)%x
                   e34 = e34/abs_value(e34)

                   c_phi = dot_product(e12,e23)
                   ! if atoms 1,2 and 3 make up a line, skip this combination
                   if (abs(abs(c_phi)-one)<=small) then
                      cycle dihedral_l
                   endif
                   c_phi = dot_product(e23,e34)
                   ! dire for atoms 2,3 and 4
                   if (abs(abs(c_phi)-one)<=small) then
                      cycle dihedral_l
                   endif

                   counter = counter+1
                   if (.not.calc) then
                      coor(counter)%typ = d_angle
                      coor(counter)%dihedral%partner1 = p1
                      coor(counter)%dihedral%base1 = b1
                      coor(counter)%dihedral%base2 = b2
                      coor(counter)%dihedral%partner2 = p2
                      coor(counter)%unique=.true.
                      coor(counter)%var = .true.
                      coor(counter)%n_equal=1
                   endif
                enddo permutation
             enddo dihedral_l
          enddo dihedral_k
       enddo dihedral_j

    enddo
    ! now the value of the counter variable is = n_valence
    n_coor = counter

    if (.not.calc .and. (keep.and. .not.exist)) then
       write(OPT_STDOUT,*)" "
       write(OPT_STDOUT,*)"val_types: writing valence coordinates to file"
       write(OPT_STDOUT,*)" "
       call write_val()
    endif

  contains
    subroutine write_val()
      ! writes the valence coordinate information to the file
      ! 'val_coor.dat'
      ! The following informattion is on the file:
      ! 1. The number of valence coordinates: n_coor
      ! 2. The valence coordinates themselves:
      !    stretches: typ,p1,p2,unique,var,n_equal
      !    angle bends: typ,p1,ap,p1,unique,var,n_equal
      !    torsions: typ,p1,b1,b2,p2,unique,var,n_equal
      ! ----------------------------------------------------
      integer(kind=i4_kind)  :: i,typ

      io_val=openget_iounit(status='replace',form='unformatted',&
           file=trim(opt_data_dir)//'/val_coor.dat')
      write(io_val)n_coor,n_cons
      do i=1,n_coor
         typ = coor(i)%typ
         write(io_val)typ
         select case (typ)
         case ( b_length )
            write(io_val)coor(i)%length%partner1,&
                         coor(i)%length%partner2, &
                         coor(i)%unique,&
                         coor(i)%var,&
                         coor(i)%n_equal
         case ( b_angle )
            write(io_val)coor(i)%angle%partner1,&
                         coor(i)%angle%apex,&
                         coor(i)%angle%partner2,&
                         coor(i)%unique,&
                         coor(i)%var,&
                         coor(i)%n_equal
         case ( d_angle )
            write(io_val)coor(i)%dihedral%partner1,&
                         coor(i)%dihedral%base1,&
                         coor(i)%dihedral%base2,&
                         coor(i)%dihedral%partner2,&
                         coor(i)%unique,&
                         coor(i)%var,&
                         coor(i)%n_equal
         case default
         end select
      enddo
      call returnclose_iounit(io_val)
    end subroutine write_val
    subroutine read_val()
      ! reads the valence coordinate information from the file
      ! 'val_coor.dat'
      ! The following information is on the file:
      ! 1. The number of valence coordinates: n_coor
      ! 2. The valence coordinates themselves:
      !    stretches: typ,p1,p2,unique,var,n_equal
      !    angle bends: typ,p1,ap,p1,unique,var,n_equal
      !    torsions: typ,p1,b1,b2,p2,unique,var,n_equal
      ! ----------------------------------------------------
      integer(kind=i4_kind)  :: i,typ

      io_val=openget_iounit(status='old',form='unformatted',&
           file=trim(opt_data_dir)//'/val_coor.dat')
      read(io_val)n_coor,n_cons
      if (calc) then
         write(OPT_STDOUT,*)" val_types/read_val: reading the dimensions only from file"
         rewind io_val
         call returnclose_iounit(io_val)
         return
      endif
      do i=1,n_coor
         read(io_val)typ
         coor(i)%typ = typ
         select case (typ)
         case ( b_length )
            read(io_val) coor(i)%length%partner1,&
                         coor(i)%length%partner2, &
                         coor(i)%unique,&
                         coor(i)%var,&
                         coor(i)%n_equal
         case ( b_angle )
            read(io_val) coor(i)%angle%partner1,&
                         coor(i)%angle%apex,&
                         coor(i)%angle%partner2,&
                         coor(i)%unique,&
                         coor(i)%var,&
                         coor(i)%n_equal
         case ( d_angle )
            read(io_val) coor(i)%dihedral%partner1,&
                         coor(i)%dihedral%base1,&
                         coor(i)%dihedral%base2,&
                         coor(i)%dihedral%partner2,&
                         coor(i)%unique,&
                         coor(i)%var,&
                         coor(i)%n_equal
         case default
         end select
      enddo
      call returnclose_iounit(io_val)
    end subroutine read_val

    subroutine perm(p,arr_in,arr_out)
      ! outputs the p=th permutation of the
      ! for indices i,j,k,l
      ! ----------------------------------
      integer(kind=i4_kind),intent(in) :: p
      integer(kind=i4_kind),intent(in):: arr_in(4)
      integer(kind=i4_kind),intent(out):: arr_out(4)
      !----------------------------------------
      arr_out = 0_i4_kind
      if (p<=6) then
         arr_out(2)=arr_in(1)
         if (p<=2) then
            arr_out(3)=arr_in(2)
            if (p==1) then
               arr_out(1)=arr_in(3)
               arr_out(4)=arr_in(4)
            else
               arr_out(1)=arr_in(4)
               arr_out(4)=arr_in(3)
            endif
         elseif(p<=4) then
            arr_out(3) = arr_in(3)
            if (p==3) then
               arr_out(1) = arr_in(2)
               arr_out(4) = arr_in(4)
            else
               arr_out(1) = arr_in(4)
               arr_out(4) = arr_in(2)
            endif
         else
            arr_out(3) = arr_in(4)
            if (p==5) then
               arr_out(1) = arr_in(2)
               arr_out(4) = arr_in(3)
            else
               arr_out(1) = arr_in(3)
               arr_out(4) = arr_in(2)
            endif
         endif
      elseif(p<=10) then
         arr_out(2) = arr_in(2)
         if (p<=8) then
            arr_out(3) = arr_in(3)
            if (p==7) then
               arr_out(1) = arr_in(1)
               arr_out(4) = arr_in(4)
            else
               arr_out(1) = arr_in(4)
               arr_out(4) = arr_in(1)
            endif
         else
            arr_out(3) = arr_in(4)
            if (p==9) then
               arr_out(1) = arr_in(1)
               arr_out(4) = arr_in(3)
            else
               arr_out(1) = arr_in(3)
               arr_out(4) = arr_in(1)
            endif
         endif
      else
         arr_out(2) = arr_in(3)
         arr_out(3) = arr_in(4)
         if (p==11) then
            arr_out(1) = arr_in(1)
            arr_out(4) = arr_in(2)
         else
            arr_out(1) = arr_in(2)
            arr_out(4) = arr_in(1)
         endif
      endif
    end subroutine perm
  end subroutine val_types

  !*************************************************************

  subroutine set_cart_for_opt(q_dummy,n_coor)
    !------------ Modules used -----------------------------------
#ifdef WITH_EFP
    use pointcharge_module, only: rcm
#endif
    implicit none
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind),intent(out)   :: q_dummy(:)
    integer(kind=i4_kind),intent(in) :: n_coor
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i
#ifdef WITH_EFP
    integer(i4_kind) :: j
#endif
    !------------ Executable code --------------------------------

    do i=1,n_atoms
       q_dummy(3*i-2)=atom(i)%x(1)
       q_dummy(3*i-1)=atom(i)%x(2)
       q_dummy(3*i)  =atom(i)%x(3)
    end do
#ifdef WITH_EFP
100 if(efp .and. n_efp > 0 .and. .not. efp_fixed) then
       j=1+n_atoms
       do i=1,n_efp
          q_dummy(3*j-2)=rcm(1,i)
          q_dummy(3*j-1)=rcm(2,i)
          q_dummy(3*j)  =rcm(3,i)
          j=j+1
          q_dummy(3*j-2)=xyz_torque(i,1)
          q_dummy(3*j-1)=xyz_torque(i,2)
          q_dummy(3*j)  =xyz_torque(i,3)
          j=j+1
       end do
    end if
#endif

  end subroutine set_cart_for_opt

  !*************************************************************
  subroutine set_internal(q_dummy,n_coor,coor,atom)
    !  Purpose: transform cartesian coordinates atom(i)%x from
    !           the 'main_opt' program to internal coordinates
    !
    !------------ Modules used ----------------------------------
    use strings, only: itoa
    implicit none
    !------------ Declaration of local variables -----------------
    type(atom_type), intent(in) :: atom(:)
    real(kind=r8_kind),intent(out)   :: q_dummy(:)
    integer(kind=i4_kind),intent(in) :: n_coor
    type(int_coor),intent(inout)        :: coor(:)
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)        :: num_coord, p1, p2, b1, b2, ap, n_sym
    real(kind=r8_kind)           :: r13,r23,r34,r12,fac3
    real(kind=r8_kind)           :: c_phi,s_phi,c_phi2,s_phi2,c_phi3,s_phi3,&
         c_tau
    real(kind=r8_kind)           :: e31(3), e32(3), e34(3), e23(3), e12(3)
    real(kind=r8_kind)           :: help1(3),help2(3)
    !------------ Executable code --------------------------------

    integer(i4_kind) :: errors

    errors = 0

    n_sym=1
    do num_coord=1,n_coor
       select case(coor(num_coord)%typ)

       case (b_length)

          p1 = coor(num_coord)%length%partner1
          p2 = coor(num_coord)%length%partner2
          r12 = abs_value(atom(p1)%x - atom(p2)%x)
          if(r12<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",p1,"and",p2,"have identical positions"
             write(OPT_STDOUT,*) "in definition of the bond ",p1,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             ! dont CYCLE ! to next coord
          end if
          q_dummy(num_coord) = r12
          if (coor(num_coord)%unique) then
             n_sym=n_sym+1
          endif

       case (b_angle)

          p1 =  coor(num_coord)%angle%partner1
          p2 =  coor(num_coord)%angle%partner2
          ap =  coor(num_coord)%angle%apex
          e31 = atom(p1)%x - atom(ap)%x
          e32 = atom(p2)%x - atom(ap)%x
          r13 = abs_value(e31)
          r23 = abs_value(e32)
          if(r13<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",p1,"and",ap,"have identical positions"
             write(OPT_STDOUT,*) "in definition of the bond angle",p1,"-",ap,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             CYCLE ! to next coord
          end if
          if(r23<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",p2,"and",ap,"have identical positions"
             write(OPT_STDOUT,*) "in definition of the bond angle",p1,"-",ap,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             CYCLE ! to next coord
          end if
          e31 = e31/r13
          e32 = e32/r23
          fac3=sqrt(r13*r23)
          c_phi = dot_product(e31,e32)
          if ( c_phi > one .and. c_phi < one+small ) c_phi = one
          if ( c_phi < -one .and. c_phi > -one-small ) c_phi = -one
          s_phi = sqrt(one - c_phi**2)
          q_dummy(num_coord) = acos(c_phi)
          if (coor(num_coord)%unique) then
             n_sym=n_sym+1
          endif

       case (d_angle)
          p1 = coor(num_coord)%dihedral%partner1
          p2 = coor(num_coord)%dihedral%partner2
          b1 = coor(num_coord)%dihedral%base1
          b2 = coor(num_coord)%dihedral%base2
          ! r12 = distance partner1 - base1
          ! r34 = distance partner2 - base2
          ! r23 = distance base1 - base2
          ! definition of the vectors e12,e34,e23 analogous

          ! c_phi2 = e12*e23
          ! c_phi3 = e34*e23

          e34 = atom(p2)%x - atom(b2)%x
          r34 = abs_value(e34)
          if(r34<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",p2,"and",b2,"have identical positions"
             write(OPT_STDOUT,*) "in definition of the dihedral angle",p1,"-",b1,"-",b2,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             CYCLE ! to next coord
          end if
          e34 = e34/r34

          e12 = atom(b1)%x-atom(p1)%x
          r12 = abs_value(e12)
          if(r12<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",b1,"and",p1,"have identical positions"
             write(OPT_STDOUT,*) "in definition of the dihedral angle",p1,"-",b1,"-",b2,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             CYCLE ! to next coord
          end if

          e12 = e12/r12

          e23 = atom(b2)%x - atom(b1)%x
          r23 = abs_value(e23)
          if(r23<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",b2,"and",b1,"have identical positions"
             write(OPT_STDOUT,*) "in definition of the dihedral angle",p1,"-",b1,"-",b2,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             CYCLE ! to next coord
          end if

          e23 = e23/r23
          fac3 = sqrt(r12*r34)
          c_phi2 = dot_product(-e12,e23)
          if ( abs(one-c_phi2)<=small ) then
             c_phi2=one
          elseif( abs(one+c_phi)<=small) then
             c_phi2=-one
          endif

          c_phi3 = dot_product(e34,-e23 )
          if (abs(one-c_phi3)<=small) then
             c_phi3=one
          elseif(abs(one+c_phi3)<=small) then
             c_phi3=-one
          endif

          help1 = cross_product(e12,e23)
          s_phi2 = abs_value(help1)
          help2 = cross_product(e23,e34)
          s_phi3 = abs_value(help2)
          if(s_phi2<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",p1,b2,"and",b1,"form a straight line"
             write(OPT_STDOUT,*) "in definition of the dihedral angle",p1,"-",b1,"-",b2,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             CYCLE ! to next coord
          end if
          if(s_phi3<small) then
             errors = errors + 1
             write(OPT_STDOUT,*) "Error",errors,":"
             write(OPT_STDOUT,*) "atoms",p2,b2,"and",b1,"form a straight line"
             write(OPT_STDOUT,*) "in definition of the dihedral angle",p1,"-",b1,"-",b2,"-",p2
             write(OPT_STDOUT,*) "proposed as internal coordinate number",num_coord
             CYCLE ! to next coord
          end if
          help1 = help1/s_phi2
          help2 = help2/s_phi3
          c_tau = dot_product(help1,help2)
          if ( abs(one-c_tau)<=small ) c_tau = one
          if ( abs(one+c_tau)<=small ) c_tau = -one
          ! now check if cross_product(-e12,e34) is parallel to e23
          help1 = cross_product(-e12,e34)
          if (dot_product(help1,e23)>=zero) then
             q_dummy(num_coord) = acos(c_tau)
          else
             q_dummy(num_coord) = -acos(c_tau)
          endif
          if (delocalized_coordinates) then
             ! This is new to Version 19: dihedral angles of values -pi are
             ! changed to +pi, which is physically identical.
             if ((q_dummy(num_coord)+pi)<small)then
                q_dummy(num_coord) = -q_dummy(num_coord)
             endif
          endif
          if(q_dummy(num_coord)<-pi.or.q_dummy(num_coord)>pi) call error_handler&
               ("sth. fishy with dihedral angle")
          if (coor(num_coord)%unique) then
             n_sym=n_sym+1
          endif

       case default

          write(OPT_STDOUT,*)" cart_to_internal : rubbish here "
          ABORT('no such coordinate type')

       end select
    enddo

    coor(:)%value = q_dummy

    if( errors > 0 )then
      call error_handler("Error: found "//trim(itoa(errors))//" bad coordinate(s), see optimizer output")
    endif

  end subroutine set_internal


  function dependent_var(n_coor,s,s_reactant,distance_to_reactant,distance_to_product)
   integer(kind=i4_kind) :: dependent_var
    real(kind=r8_kind), intent(out):: distance_to_reactant
    real(kind=r8_kind), intent(out), optional:: distance_to_product
   integer(kind=i4_kind), intent(in):: n_coor
    type(int_coor), intent(in) :: s(:),s_reactant(:)
    real(kind=r8_kind):: max_contib,contib,sqt,sqmt

   integer(kind=i4_kind):: i

   distance_to_reactant=0.0_r8_kind
   if(present(distance_to_product)) distance_to_product=0.0_r8_kind
   max_contib=0.0
   dependent_var=0
   do i=1,n_coor
    if(select_sphere_vars.and..not.s(i)%sphere) cycle
    contib=(s(i)%value-s_reactant(i)%value)**2
    if(contib.gt.max_contib) then
     max_contib=contib
     dependent_var=i
    endif
    distance_to_reactant=distance_to_reactant+contib
    if(present(distance_to_product)) then
     distance_to_product=distance_to_product+(s(i)%value-s_product(i)%value)**2
     write(io_flepo,*) i,contib,s(i)%value,s_reactant(i)%value,s_product(i)%value,&
                         'contib s s_r s_p',(s(i)%value-s_product(i)%value)**2
    else
     write(io_flepo,*) i,contib,s(i)%value,s_reactant(i)%value,'contib s s_r'
    endif
   enddo
ASSERT(distance_to_reactant.gt.small)
   distance_to_reactant=sqrt(distance_to_reactant)
   if(present(distance_to_product)) then
     distance_to_product=sqrt(distance_to_product)
     tsscan_rp_var=distance_to_product/(distance_to_product+distance_to_reactant)
     write(io_flepo,*) 'tsscan_rp_var ', tsscan_rp_var,distance_to_reactant,distance_to_product
   sqt=tsscan_rp_var**2
   sqmt=(1.0_r8_kind-tsscan_rp_var)**2
   max_contib=0.0_r8_kind
   do i=1,n_coor
    if(select_sphere_vars.and..not.s(i)%sphere) cycle
    contib=abs((s(i)%value-s_reactant(i)%value)*sqt-sqmt*(s(i)%value-s_product(i)%value))
    if(contib.gt.max_contib) then
     max_contib=contib
     dependent_var=i
    endif
   enddo
   endif
  end function dependent_var

  subroutine print_internal (coor, n_coor, n_val, full)
    ! Purpose: provide a nice output of the internal coordinates
    !          lengths in au/angstrom
    !          angles in radians/degrees
    !
    !--------------------------------------------------------
    use constants, only: angstrom
    implicit none
    type(int_coor), intent(in) :: coor(:)
    integer(i4_kind), intent(in) :: n_coor
    integer(i4_kind), optional :: n_val
    logical, optional, intent(in) :: full
    ! *** end of interface ***

    real(kind=r8_kind)          :: fac_length1,fac_length2,&
                                   fac_angle1,fac_angle2
    integer(kind=i4_kind),allocatable :: sym_loc(:)
    integer(kind=i4_kind)             :: alloc_stat
    ! --- Declaration of local variables --------------------
    integer(kind=i4_kind) :: i
    logical :: full_loc

    full_loc=.false.
    if (present(full) ) then
       if (full) then
          full_loc=.true.
       else
          full_loc=.false.
       endif
    else
       full_loc=.false.
    endif

    if (present(n_val)) then
       allocate(sym_loc(n_val),STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler&
            ("print_internals: allocation (1) failed")
       sym_loc=0_i4_kind
    else
       allocate(sym_loc(n_primitive),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            ("print_internals: allocatin (2) failed")
       sym_loc = sym_type
    endif

    if (.not.full_loc) then
       write(OPT_STDOUT,*)"   No.      au/rad                  angstr./deg.&
            &                        type  mult"
       write(OPT_STDOUT,*)"------------------------------------------------&
            &--------------------------------------"
    else
       write(OPT_STDOUT,*)"   No.      au/rad                  angstr./deg.&
            &                                    type"
       write(OPT_STDOUT,*)"------------------------------------------------&
            &------------------------------------------"
    endif


    if (atomic_units) then
       fac_length1 = one
       fac_length2 = one / angstrom
       fac_angle1=one
       fac_angle2=convert1
    else
       fac_length1 = angstrom
       fac_length2 = one
       fac_angle1 = one
       fac_angle2 = convert1
    endif
    do i=1,n_coor
       if (.not.full_loc.and. .not.delocalized_coordinates) then
          if ( coor(i)%unique  ) then
             if (coor(i)%typ == b_length) then
                if (coor(i)%var ) then
                   write(OPT_STDOUT,1000)i,coor(i)%value*fac_length1,coor(i)%value*fac_length2,&
                        adjustr(trim(bond_type(coor(i)%typ))),sym_loc(i),coor(i)%n_equal
                endif
             else
                if (coor(i)%var ) then
                   write(OPT_STDOUT,1000)i,coor(i)%value*fac_angle1,coor(i)%value*fac_angle2,&
                        adjustr(trim(bond_type(coor(i)%typ))),sym_loc(i),coor(i)%n_equal
                endif
             endif
          endif
       else
          if (coor(i)%typ == b_length) then
             if (coor(i)%var ) then
                write(OPT_STDOUT,1001)i,coor(i)%value*fac_length1,coor(i)%value*fac_length2,&
                     adjustr(trim(bond_type(coor(i)%typ))),sym_loc(i)
             else
                write(OPT_STDOUT,1011)i,coor(i)%value*fac_length1,coor(i)%value*fac_length2,&
                     adjustr(trim(bond_type(coor(i)%typ))),sym_loc(i)
             endif
          else
             if (coor(i)%var ) then
                write(OPT_STDOUT,1001)i,coor(i)%value*fac_angle1,coor(i)%value*fac_angle2,&
                     adjustr(trim(bond_type(coor(i)%typ))),sym_loc(i)
             else
                write(OPT_STDOUT,1011)i,coor(i)%value*fac_angle1,coor(i)%value*fac_angle2,&
                     adjustr(trim(bond_type(coor(i)%typ))),sym_loc(i)
             endif
          endif
       endif
    enddo
    if (.not.full_loc) then
       write(OPT_STDOUT,*)"------------------------------------------------&
            &--------------------------------------"
    else
       write(OPT_STDOUT,*)"------------------------------------------------&
            &------------------------------------------"
    endif
    write(OPT_STDOUT,*)" "

    if (full_loc) then
       do i=1,n_coor
          if (coor(i)%typ == b_length) then
             write(OPT_STDOUT,'("Bond Length ",i3,"      ",i3,"  ",i3)')i,&
                  coor(i)%length%partner1,&
                  coor(i)%length%partner2
          elseif(coor(i)%typ == b_angle ) then
             write(OPT_STDOUT,'(" Bond Angle ",i3,"      ",i3,"  ",i3,"  ",i3)')i,&
                  coor(i)%angle%partner1,&
                  coor(i)%angle%apex,&
                  coor(i)%angle%partner2
          elseif(coor(i)%typ == d_angle ) then
             write(OPT_STDOUT,'("   Dihedral ",i3,"      ",i3,"  ",i3,"  ",i3,"  ",i3)')i,&
                  coor(i)%dihedral%partner1,&
                  coor(i)%dihedral%base1,&
                  coor(i)%dihedral%base2,&
                  coor(i)%dihedral%partner2
          endif
       enddo
    endif

    deallocate(sym_loc,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" print_internals: deallocation (1) failed")

1000 format(3x,i3,'   ',f18.12,5x,f18.12,5x,A16,2x,i3,2x,i3)
!1010 format(3x,i3,'   ',f18.12,5x,f18.12,5x,A16,'   constant',2x,i3,2x,i3)
1001 format(3x,i3,'   ',f18.12,5x,f18.12,5x,A16,'           ',2x,i3)
1011 format(3x,i3,'   ',f18.12,5x,f18.12,5x,A16,'   constant',2x,i3)
  end subroutine print_internal

!*************************************************************
  subroutine generate_bmat(coor,n_coor,bmat_in)
    !  Purpose: generate Wilsons B-Matrix
    !           This routine assumes the the types of internal
    !           coordinates are already set ( -> 'zmat_types').
    !           The 'bmat' is generated using the set of coordinates
    !           'q'.
    !           All angles are multiplied with sqrt(r12*r23) to
    !           give the same unit as bond lengths.
    !
    ! References: E.B.Wilson Jr, J.C. Decius and P.C.Cross,
    !             'Molecular Vibrations',(McGraw-Hill,New York, 1955)
    !            pages: 54 ff.
    !------------ Modules used ----------------------------------
    ! --- Declaration of formal parameters ----------------------
    type(int_coor),intent(in)          :: coor(:)
    integer(kind=i4_kind),intent(in)   :: n_coor
    real(kind=r8_kind),intent(out)   :: bmat_in(:,:)
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)        :: num_coord, p1, p2, b1, b2, ap
    real(kind=r8_kind)           :: r13,r23,r34,r12
    real(kind=r8_kind)           :: c_phi, s_phi, c_phi2, s_phi2, c_phi3, s_phi3
    real(kind=r8_kind)           :: e31(3), e32(3), e34(3), e23(3), e12(3)
    real(kind=r8_kind)           :: help1(3), help2(3)
    real(kind=r8_kind)           :: fac1,fac2,fac3,test
    integer(kind=i4_kind)        :: start1, start2, start3, start4
    integer(kind=i4_kind),parameter  :: n_col = 10
    ! number of columns for b-matrix output
    !------------ Executable code --------------------------------

    ! ATTENTION: for this subroutine it might be the case that
    ! ubound(coor) /= n_coor, namely when th B-matrix for the
    ! valence coordinates is done. There we dont know from the
    ! beginning on how many valence coordinates we will have
    ! and thus the coordinate vector is allocated with a fixed maximum
    ! number.
    bmat_in=zero

    do num_coord=1,n_coor
       if (coor(num_coord)%typ.eq.b_length ) then

          p1 = coor(num_coord)%length%partner1
          p2 = coor(num_coord)%length%partner2
          start1=3*(p1-1)+1
          start2=3*(p2-1)+1
          r12 = abs_value(atom(p1)%x-atom(p2)%x)
          bmat_in(num_coord,start1:start1+2) = (atom(p1)%x-atom(p2)%x)/r12
          bmat_in(num_coord,start2:start2+2) = &
               -bmat_in(num_coord,start1:start1+2)

       elseif (coor(num_coord)%typ == b_angle) then

          p1 =  coor(num_coord)%angle%partner1
          p2 =  coor(num_coord)%angle%partner2
          ap =  coor(num_coord)%angle%apex
          e31 = atom(p1)%x - atom(ap)%x
          e32 = atom(p2)%x - atom(ap)%x
          r13 = abs_value(e31)
          r23 = abs_value(e32)
          fac3 = sqrt(r13*r23)
          e31 = e31/r13
          e32 = e32/r23
          c_phi = dot_prod(e31,e32)
          if ( c_phi > one .and. c_phi < one+small ) c_phi = one
          if ( c_phi < -one .and. c_phi > -one-small ) c_phi = -one
          s_phi = sqrt(one - c_phi**2)
          start1=3*(p1-1)+1
          start2=3*(p2-1)+1
          start3=3*(ap-1)+1


          bmat_in(num_coord,start1:start1+2) = (c_phi*e31 - e32)/(r13*s_phi)
          bmat_in(num_coord,start2:start2+2) = (c_phi*e32 - e31)/(r23*s_phi)
          bmat_in(num_coord,start3:start3+2) = ( (r13-r23*c_phi)*e31 + &
               (r23-r13*c_phi)*e32 ) / ( r13*r23*s_phi )

       elseif (coor(num_coord)%typ == d_angle) then

          p1 = coor(num_coord)%dihedral%partner1
          p2 = coor(num_coord)%dihedral%partner2
          b1 = coor(num_coord)%dihedral%base1
          b2 = coor(num_coord)%dihedral%base2
          start1 = 3*(p1-1)+1
          start2 = 3*(b1-1)+1
          start3 = 3*(b2-1)+1
          start4 = 3*(p2-1)+1

          ! r12 = distance partner1 - base1
          ! r34 = distance partner2 - base2
          ! r23 = distance base1 - base2
          ! definition of the vectors e12,e34,e23 analogous

          e34 = atom(p2)%x - atom(b2)%x
          r34 = abs_value(e34)
          e34 = e34/r34

          e12 = atom(b1)%x-atom(p1)%x
          r12 = abs_value(e12)

          e12 = e12/r12

          e23 = atom(b2)%x - atom(b1)%x
          r23 = abs_value(e23)

          e23 = e23/r23
          fac3 = sqrt(r12*r34)

          c_phi2 = dot_product(-e12,e23)
          if ( abs(one-c_phi2)<=small) then
             c_phi2 = one
          elseif (abs(one+c_phi2)<=small ) then
             c_phi2 = -one
          endif

          c_phi3 = dot_product(-e23,e34 )
          if ( abs(one-c_phi3)<=small ) then
             c_phi3 = one
          elseif( abs(one+c_phi3)<=small ) then
             c_phi3 = -one
          endif

          help1 = cross_product(e12,e23)
          help2 = cross_product(e23,e34)
          s_phi2 = abs_value(help1)
          s_phi3 = abs_value(help2)

          ! stelle fest ob beide Winkel < oder > pi sind ODER
          ! ob einer < pi und einer > pi ist.
          test = dot_product(help1,-help2)
          if ( test < zero ) then
             s_phi3=-s_phi3
          endif
          help1=help1/s_phi2
          help2=help2/s_phi3

          bmat_in(num_coord,start1:start1+2) = -cross_product(e12,e23)/(r12*s_phi2**2)
          bmat_in(num_coord,start4:start4+2) = -cross_product(e34,e23)/(r34*s_phi3**2)

          fac1 = (r23-r12*c_phi2)/(r12*r23*s_phi2)
          fac2 = c_phi3/(r23*s_phi3)

          bmat_in(num_coord,start2:start2+2) = (fac1*help1 - fac2*help2)

          fac1 = (r23 - r34*c_phi3) / (r23*r34*s_phi3)
          fac2 = c_phi2 / (r23*s_phi2)

          bmat_in(num_coord,start3:start3+2) = -(fac1*help2 - &
               fac2*help1)
       else
          write(OPT_STDOUT,*) " generate_bmat: unknown type",coor(num_coord)%typ &
                            , " of coordinate number",num_coord                  &
                            , " having the value",coor(num_coord)%value
          call error_handler("generate_bmat: rubbish here, see optimizer output")
       endif
    enddo
    if(print_bmat) then
       write(OPT_STDOUT,*)" Primitive B-Matrix  "
       call print_matrix(bmat_in,n_coor,3*(n_atoms+n_dummy),n_col)
       write(OPT_STDOUT,*)" "
    endif
!1201 format(i4,20(2x,f9.5))
!1202 format(4x,20(8x,i3))

  end subroutine generate_bmat

     !*************************************************************
  subroutine invert_bmat(bmat,bmat_inv,spectroscopic)
    ! Purpose: construct the general right-inverse of bmat.
    !          if 'spectroscopic' is set to TRUE, the inverse is
    !          built using the atomic masses, such that the
    !          by-product matrix G = B*M*Bt is the spectroscopic
    !          matrix used in a frequency calculation.
    !          Otherwise simply the unit matrix is taken.
    !
    ! Extension for fixed orientation (8.6.98) FN. Here 6 additional
    ! constraints that fix the orientation of the molecule
    ! are added to the bmatrix. Then, being quadratic, it is inverted
    ! using the standard invert_matrix() procedure.
    !
    ! ----- Modules used ---------------------------------------
    use atom_data_module
    ! ----- Declaration of formal parameters -------------------
    implicit none
    logical,intent(in) :: spectroscopic
    real(r8_kind), intent(in)  :: bmat(:,:)
    real(r8_kind), intent(out) :: bmat_inv(:,:)
    ! *** end of interface ***

    ! ----- Declaration of local variables ---------------------
    real(kind=r8_kind),allocatable :: m(:,:),gmat(:,:), gmat_inv(:,:),&
         help_mat(:,:), help_vec(:,:), help_vec_dummy(:,:), help_vec_heavy(:,:), &
         x_heavy(:), y_heavy(:), z_heavy(:), m_heavy(:)
    integer(kind=i4_kind)          :: alloc_stat, i, j, n_local, i_atom, counter
    integer(kind=i4_kind)          :: i_help
    real(kind=r8_kind)             :: coord(3), &
         center(3), diff_vec(3),rotaxes(3,3)
    real(kind=r8_kind),allocatable :: constraint_mat(:,:)
    ! all things for pc
    real(kind=r8_kind),dimension(3):: con_vec,vec_1,vec_2,normal_vec
    real(kind=r8_kind),parameter,dimension(3) :: e1=(/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind/),&
                                                 e2=(/0.0_r8_kind,1.0_r8_kind,0.0_r8_kind/),&
                                                 e3=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
    real(kind=r8_kind)             :: abs_e
    logical :: crash ! for debug only
    ! ----- executable code ------------------------------------

    n_local = size(bmat,1) ! this will be either n_primitive or n_internal

    allocate(m(3*(n_atoms+n_dummy),3*(n_atoms+n_dummy)),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    allocate(gmat(n_local,n_local),gmat_inv(n_local,n_local),STAT=allocopt_stat(9))
    ASSERT(allocopt_stat(9).eq.0)

    m = zero

    if (.not.fixed_orientation) then

       if (.not.spectroscopic) then
          allocate(help_mat(3*(n_atoms+n_dummy),n_local),stat=alloc_stat)
          ASSERT(alloc_stat.eq.0)

          do i=1,3*(n_atoms+n_dummy)
             m(i,i) = one
          enddo

          help_mat = matmul(m,transpose(bmat)) ! m==1??
          gmat = matmul(bmat,help_mat)
          gmat_inv = gmat ! just as a helparray for the LAPACK routines

          call invert_matrix(n_local,gmat_inv)

          ! now test this by eigensolving ...
          ! -> dgeev.f LAPACK

          bmat_inv = matmul(help_mat,gmat_inv)

#if 0 /* mass_weigted */
!     gmat_sqr=gmat_inv
!      gmat_sqr=matmul(bmat_inv,transpose(bmat_inv))
!      gmat_sqr=matmul(bmat,bmat_inv)

!       gmat_sqr=matmul(gmat,gmat_inv)

       gmat_sqr=matmul(transpose(bmat),bmat)
       gmat_sqr_inv=gmat_sqr
!       call invert_matrix(3*(n_atoms+n_dummy),gmat_sqr_inv)

!      write(OPT_STDOUT,*) sum(matmul(gmat_sqr_inv,transpose(gmat_sqr))),' gmat_sqr'

     call invert_matrix(n_local,gmat_sqr)
    call  rs(n_local,n_local,gmat_inv,eig_real,1_i4_kind,eig_vec,fv1,fv2,ierr)
#endif

!!$!   This is a test for conservation of geometrical (or mass) center
!!$    ! test, if sum(i)[bmat_inv(i,j_internal)] = 0
!!$    do j=1,n_internal
!!$       summ_1 = zero
!!$       summ_2 = zero
!!$       summ_3 = zero
!!$       jj=1
!!$       do i=1,3*(n_atoms+n_dummy),3
!!$          summ_1 = summ_1 + bmat_inv(i,j)*mass(nint(atom(jj)%charge))
!!$          summ_2 = summ_2 + bmat_inv(i+1,j)*mass(nint(atom(jj)%charge))
!!$          summ_3 = summ_3 + bmat_inv(i+2,j)*mass(nint(atom(jj)%charge))
!!$          jj=jj+1
!!$       enddo
!!$       write(OPT_STDOUT,'(" Internal ",i2,"    Sum(i)[bmat_inv(i,",i2,")    ",3(2x,f9.5))')j,j,summ_1,summ_2,summ_3
!!$    enddo
          ! if no frequency calculation is desired, gmat can be thrown away
          if (print_gmat) then
             write(OPT_STDOUT,*)" G-Matrix (B * M * Bt ) "
             call print_matrix(gmat,n_local,n_local,10_i4_kind)
             write(OPT_STDOUT,*)" "
          endif
          ! dito for m, help_mat
          deallocate(help_mat,STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

!1000      format(40(f11.7))
!2000      format(20(2x,f8.4))
       else ! i.e. spectroscopic
          ! calculate the inverse by extending Bmat by 6 additional
          ! lines, which represent additional conditions for translational
          ! and rotational invariance. Atoms with vanishing masses (dummy atoms)
          ! are ommitted in these conditions. In the limit of atoms with infinite
          ! masses, the masses of the infenitely heavy atoms are formaly set to 1,
          ! wheras the mass of all other atoms is set to zero.
          if(n_infinite_masses>0) then
             ! we have havy atoms
             allocate(x_heavy(n_infinite_masses), y_heavy(n_infinite_masses),&
                  z_heavy(n_infinite_masses), m_heavy((n_infinite_masses)), &
                  stat = alloc_stat)
             m_heavy=1.0_r8_kind
             counter=1
             do i_atom=1,n_tot_atoms
                if(atom(i_atom)%heavy) then
                   x_heavy(counter)=atom(i_atom)%x(1)
                   y_heavy(counter)=atom(i_atom)%x(2)
                   z_heavy(counter)=atom(i_atom)%x(3)
                   counter=counter+1
                end if
             end do
!             if(mass_center_print) &
             call mass_center(x_heavy,y_heavy,z_heavy,&      !(1) in invert_bmat
                  m_heavy,n_infinite_masses,center,rotaxes)

             deallocate(x_heavy, y_heavy, z_heavy, m_heavy, &
                  stat = alloc_stat)
             if (alloc_stat/=0) call error_handler&
                  (" invert_bmat : deallocation (2) failed")
          else
             ! calculate center of mass and main rotational axes
!             if(mass_center_print) &
             call mass_center(atom(1:n_tot_atoms)%x(1),atom(1:n_tot_atoms)%x(2),atom(1:n_tot_atoms)%x(3),&  !(2) in invert_bmat
                  atom(1:n_tot_atoms)%mass,n_tot_atoms,center,rotaxes)
          end if
          ! now copy b_mat to help_mat
          allocate(help_mat(3*n_tot_atoms,3*n_tot_atoms), &
               help_vec(6,3*n_tot_atoms), &
               help_vec_dummy(6,3*n_tot_atoms), &
               help_vec_heavy(6,3*n_tot_atoms), &
               stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               (" invert_bmat : allocation (3) failed")

          help_mat=0.0_r8_kind
          help_mat(1:n_local,:)=bmat(1:n_local,:)
          ! now we have to add the additional conditions for
          ! translational and rotational invariance after update
          ! first translations
          help_vec_dummy=0.0_r8_kind
          help_vec_heavy=0.0_r8_kind
          help_vec=0.0_r8_kind
          do j=1,3
             do i_atom=1,n_tot_atoms
                if(atom(i_atom)%heavy) then
                   help_vec_heavy(j,(i_atom-1)*3+j)=1.0_r8_kind
                elseif(atom(i_atom)%dummy) then
                   help_vec_dummy(j,(i_atom-1)*3+j)=1.0_r8_kind
                else
                   help_vec(j,(i_atom-1)*3+j)=atom(i_atom)%mass
                end if
             enddo
          enddo
          ! now rotations
          do i_atom=1,n_tot_atoms
             diff_vec=atom(i_atom)%x-center(:)
             ! now lets determine coordinates in coordinate system
             ! defined in main axis coordinate system
             coord(1)=sum(rotaxes(:,1)*diff_vec)
             coord(2)=sum(rotaxes(:,2)*diff_vec)
             coord(3)=sum(rotaxes(:,3)*diff_vec)
             if(atom(i_atom)%dummy) then ! dummy atoms
                help_vec_dummy(4,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     (coord(1)*rotaxes(:,2)-coord(2)*rotaxes(:,1))
                help_vec_dummy(5,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     (coord(2)*rotaxes(:,3)-coord(3)*rotaxes(:,2))
                help_vec_dummy(6,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     (coord(1)*rotaxes(:,3)-coord(3)*rotaxes(:,1))
             elseif(atom(i_atom)%heavy) then ! heavy atoms
                help_vec_heavy(4,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     (coord(1)*rotaxes(:,2)-coord(2)*rotaxes(:,1))
                help_vec_heavy(5,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     (coord(2)*rotaxes(:,3)-coord(3)*rotaxes(:,2))
                help_vec_heavy(6,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     (coord(1)*rotaxes(:,3)-coord(3)*rotaxes(:,1))
             else ! "normal" atoms
                help_vec(4,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     atom(i_atom)%mass*(coord(1)*rotaxes(:,2)-coord(2)*rotaxes(:,1))
                help_vec(5,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     atom(i_atom)%mass*(coord(2)*rotaxes(:,3)-coord(3)*rotaxes(:,2))
                help_vec(6,(i_atom-1)*3+1:(i_atom-1)*3+3)= &
                     atom(i_atom)%mass*(coord(1)*rotaxes(:,3)-coord(3)*rotaxes(:,1))
             end if
          enddo

          !
          ! FIXME: code below is broken, it behaves differently for diatomics
          ! oriented along Z-, X-, or Y-axis and randomly oriented molecule.
          !
          if ( n_tot_atoms == 2 ) then
            WARN("Code for linear molecules is broken!")
          endif

          i_help=1
          do i=1,6
             ! FIXME: these if's are evil:
             if(sum(help_vec_heavy(i,:)**2)<small) then
                ! these constraint can not be used, use normal atoms instead
                if(sum(help_vec(i,:)**2)<small) then
                   if(sum(help_vec_dummy(i,:)**2)<small) then
                      ! dont increment i_help!
                      if(n_tot_atoms/=2)  call error_handler('invert_bmat: something fishy')
                   else
                      crash = n_local+i_help > size(help_mat,1)
                      ASSERT(.not.crash)
                      help_mat(n_local+i_help,:)=help_vec_dummy(i,:)
                      i_help=i_help+1
                   end if
                else
                   crash = n_local+i_help > size(help_mat,1)
                   ASSERT(.not.crash)
                   help_mat(n_local+i_help,:)=help_vec(i,:)
                   i_help=i_help+1
                end if
             else
                crash = n_local+i_help > size(help_mat,1)
                ASSERT(.not.crash)
                help_mat(n_local+i_help,:)=help_vec_heavy(i,:)
                WARN("FIXME: for sure i and not 1?")
                i_help=i_help+i ! FIXME: for sure i and not 1?
             end if
          end do

          ! now we can already invert help_mat
          call invert_matrix(3*n_tot_atoms,help_mat)
          ! and copy staff into b_mat_inv
          bmat_inv(:,:)=help_mat(:,1:n_primitive)
          ! b_mat_inv contains the derivatives of the cartesian coordinates with
          ! respect to the internal coordinates
          deallocate(help_mat, help_vec, help_vec_dummy, help_vec_heavy,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               'b_mat_inv: deallocating (3) failed')
       endif ! .not.spectroscopic

    else  ! fixed_orientation
!   print*,'fixed_orientation fixed_orientation',fixed_orientation
       allocate(constraint_mat(3*(n_atoms+n_dummy),3*(n_atoms+n_dummy)),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            (" invert_bmat : allocation (4) failed")
       constraint_mat=zero
       constraint_mat(1:n_local,:) = bmat(1:n_local,:)

       ! fixation of atom 'fixed_atom_1' to its original position --------
       j=(fixed_atom_1-1)*3+1
       do i = n_local+1,n_local+3
          ! erster Spaltenindex fuer atom i_atom ist j=(i_atom-1)*3+1
          constraint_mat(i,j)=one
          j=j+1
       enddo
       ! -----------------------------------------------------------------

       con_vec = atom(fixed_atom_2)%x - atom(fixed_atom_1)%x
       con_vec=con_vec/abs_value(con_vec)
       ! Using vectors e1=(1,0,0),e2=(0,1,0) and e3=(0,0,1)
       ! find two vectors orthogonal to con_vec with a little Schmidt-procedure---
       vec_1 = e1 - dot_product(e1,con_vec)*con_vec
       abs_e = abs_value(vec_1)
       if (abs_e<small) then ! e1 was (anti)parallel to con_vec and thus cannont be used
          !                    which means that e2 and e3 should be ok
          ! test this hypothesis:
          abs_e = dot_product(e2,con_vec)
          if ( .not. (dot_product(e2,con_vec)<small .and. dot_product(e3,con_vec)<small) ) then
             call error_handler("invert_bmat: sth. fishy here")
          endif
          vec_1 = e2
          vec_2 = e3
       else ! this means e1 is not parallel (or antiparallel) to the connecvtion vector
          vec_1=vec_1/abs_e
          ! now try to construct vec_2 with e2
          vec_2 = e2 - dot_product(e2,con_vec)*con_vec - dot_product(e2,vec_1)*vec_1
          abs_e = abs_value(vec_2)
          if ( abs_e<small ) then ! this means that only e3 is left to make up vec_2
             vec_2 = e3 - dot_product(e3,con_vec)*con_vec - dot_product(e3,vec_1)*vec_1
             abs_e = abs_value(vec_2)
             if (abs_e < small ) then
                write(OPT_STDOUT,*)" invert_bmat: Vector VEC_2 orthogonal to the connection line"
                write(OPT_STDOUT,*)"              of atoms ",fixed_atom_1," and ",fixed_atom_2
                write(OPT_STDOUT,*)"              could not be constructed. Seek technical assistance"
                call error_handler(" ")
             endif
             vec_2 = vec_2/abs_e
          else !this means e2 was ok
             vec_2=vec_2/abs_e
          endif
       endif

       if (print_debug) then
          write(OPT_STDOUT,*)" invert_bmat: vectors for fixation of the ",fixed_atom_2," atoms movement"
          write(OPT_STDOUT,*)"              along the connection vector from atom ",fixed_atom_1," to "
          write(OPT_STDOUT,*)"              atom ",fixed_atom_2
          write(OPT_STDOUT,'("Connection Vector :",3(2x,f8.4))')(con_vec(i),i=1,3)
          write(OPT_STDOUT,'("Vector vec_1      :",3(2x,f8.4))')(vec_1(i),i=1,3)
          write(OPT_STDOUT,'("Vector vec_2      :",3(2x,f8.4))')(vec_2(i),i=1,3)
       endif

       ! now find out column j of the b-matrix belonging to i_atom
       j=(fixed_atom_2-1)*3+1
       do i=1,3
          constraint_mat(n_local+4,j) = vec_1(i)
          constraint_mat(n_local+5,j) = vec_2(i)
          j=j+1
       enddo

       normal_vec = cross_product(atom(fixed_atom_2)%x-atom(fixed_atom_1)%x,atom(fixed_atom_3)%x-atom(fixed_atom_1)%x)
       abs_e = abs_value(normal_vec)
       if (abs_e<small) then
          write(OPT_STDOUT,*)" invert_bmat: the normal vector of the plane spanned by atoms"
          write(OPT_STDOUT,*)"             ",fixed_atom_1,", ",fixed_atom_2," and ",fixed_atom_3
          write(OPT_STDOUT,*)"              is rubbish. Seek technical assistance"
          call error_handler(" ")
       endif

       normal_vec=normal_vec/abs_e
       ! find out column for atom i_atom
       j=(fixed_atom_3-1)*3+1
       do i=1,3
          constraint_mat(n_local+6,j) = normal_vec(i)
          j=j+1
       enddo

       ! finally constraint_mat can be inverted using standard methods
       call invert_matrix(3*(n_atoms+n_dummy),constraint_mat)
       ! ... and the important stuff is copied to bmat_inv
       bmat_inv(:,:) = constraint_mat(:,1:n_local)
       deallocate(constraint_mat,STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            ("invert_bmat: deallocation (4) failed")

    endif

    if(allocated(gmat)) then
     deallocate(gmat,gmat_inv,m,STAT=allocopt_stat(9))
     ASSERT(allocopt_stat(9).eq.0)
     allocopt_stat(9)=1
    endif

  if(allocated(m)) then
   deallocate(m,STAT=alloc_stat)
   ASSERT(alloc_stat.eq.0)
  endif

  end subroutine invert_bmat

     !*************************************************************
  subroutine fake_internals(new_internal, sym_int, q)
    ! Purpose: preliminary ... invent *some* updated internal
    !          coordinates to test the update procedure.
    real(r8_kind), intent(in) :: new_internal(:)  ! (num_int)
    integer(i4_kind), intent(in) :: sym_int(:)  ! (num_int)
    real(r8_kind), intent(inout) :: q(:)   ! (n_internal)
    ! *** end of interface ***

    integer(kind=i4_kind)  :: i, num_int

    DPRINT 'inside fake_internals'

    num_int=ubound(new_internal,1)
    print*,'This is the new list of (primitive) internals',new_internal

    if (zmat_coordinates) then
       do i=1,num_int
       if(ts_scan) then
        q(sym_int(i))=q(sym_int(i))+new_internal(i)
       else
        q(sym_int(i))=new_internal(i)
       endif
       enddo
    elseif (delocalized_coordinates) then
       print*,' Ubounds of umat_trans   :',ubound(umat_trans,1),ubound(umat_trans,2)

       q = matmul(umat_trans,new_internal)
    endif
    write(OPT_STDOUT,*)
    write(OPT_STDOUT,*)" fake_internals : q, q_old "
    do i=1,n_internal
       write(OPT_STDOUT,'(i4, 2(3x,F9.4),3x,"Step  ",f9.4)')i,q(i),q_old(i),q(i)-q_old(i)
    enddo
  end subroutine fake_internals

  subroutine internal_to_cart(q_old,q,bmat,bmat_inv)
    !
    ! FIXME: add cartesian geometry to the argument list
    !
    ! Purpose: transform the internal displacement coordinates
    !          Qold-Qnew to cartesian displacements Xold-Xnew
    !          using the general right inverse of  Wilsons
    !          B-Matrix:
    !          (1) Xnew = Xold + bmat_inv[Xold]*(Qnew-Qold)
    !          (2) construct a new bmat[Xnew] and a new
    !              bmat_inv[Xnew].
    !          (3) use dX = Xnew-Xold to determine an intermediate
    !              dQhelp = Qnew_help - Qold
    !          (4) compare this to the original dQ=Qnew-Qold.
    !              If not sufficient go to step (1).
    !
    ! References: P.Pulay and G.Fogarasi:
    !        'Geometry optimization in redundant internal coordinates'
    !        J.Chem.Phys. 96, (1992), 2856
    !        and references therein.
    !
    ! Extension to delocalized internal coordinates:
    !          Step (2): construct a new B-Matrix for each
    !                    SCF-step but leave the transformation
    !                    Matrix umat_trans CONSTANT. Otherwise
    !                    the backtransformation will not converge.
    !
    ! References: J.Baker, A.Kessi and B.Delley:
    !             ' The generation and use of delocalized internal'
    !             '  coordinates in geometry optimization'
    !             J.Chem.Phys. 105, (1) July 1996
    !
    ! Attention: n_coor is set to the upper bound of bmat in order
    !            to be sure that all upper llimits are coorect for the
    !            case:
    !            -zmat_format AND zmat_coordinates, where the n_primitive
    !             internal coordinates are symmetry-reduced to n_internal
    !             ones.
    !            -delocalized_coordinates or zmat_coordinates AND free_format
    !             where all internal coordinates (symmetry-equivalent or
    !             not) form a set of n_internal coordinates.
    !            -frequency_calculation OR ts_search in delocalized_coordinates
    !             where the Hessian has to be set up in symmetry-unique
    !             coordinates. FN 3/98
    !
    ! ---- Modules ----------------------------------------------
    implicit none
    real(kind=r8_kind), intent(in)    :: q_old(:) ! (n_internal?) values of internal coordinates
    real(kind=r8_kind), intent(in)    :: q(:)     ! (n_internal?) updated values
    real(kind=r8_kind), intent(inout) :: bmat(:,:)     ! IN: at point q_old, OUT: at point q
    real(kind=r8_kind), intent(inout) :: bmat_inv(:,:) ! IN: at point q_old, OUT: at point q
    ! *** end of interface ***

    ! ---- Declaration of local variables ------------------------
    logical  :: converged
    integer(kind=i4_kind)    :: i,start,loop,&
         alloc_stat,j,sumi,n_coor
    integer(kind=i4_kind),parameter :: max_iteration=100
    real(kind=r8_kind),allocatable  :: dx(:,:),q_ist(:),q_new(:),&
         q_prim_old(:),q_keep(:)
    real(kind=r8_kind)              :: delta,delta_max
    real(kind=r8_kind)              :: xgeo(3),xmass(3),rotaxes(3,3)
    type(int_coor),pointer          :: s_local(:)
    logical                         :: contract
    integer(kind=i4_kind)           :: contract_counter
    ! ---- executable code ---------------------------------------

    contract=.false.
    contract_counter=1

    DPRINT 'internal_to_cart: shape(bmat)=',shape(bmat)

    n_coor = size(bmat,1)
    allocate(dx(3*(n_atoms+n_dummy),3),&
         q_ist(n_coor),q_new(n_coor),q_keep(n_coor),&
         STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    if (zmat_format.and.zmat_coordinates) then
       s_local => s_prim
    else
       s_local => s
       allocate(q_prim_old(n_primitive),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       q_prim_old=q_prim
    endif

    ! calculate geometrical center:

    call geo_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),(n_atoms+n_dummy),xgeo)
    call mass_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),atom(:)%mass,(n_atoms+n_dummy),xmass,rotaxes) !(3) in internal_to_cart

    if(mass_center_print)then
     write(OPT_STDOUT,'("Geometrical center BEFORE update ",3(2x,f9.4))')(xgeo(i),i=1,3)
     write(OPT_STDOUT,'("Mass center BEFORE update        ",3(2x,f9.4))')(xmass(i),i=1,3)
     write(OPT_STDOUT,*)" "
    endif

    DPRINT 'internal_to_cart: internal coordinates q=',q
    converged = .false.
    loop = 0
    if (zmat_format.and.zmat_coordinates) then
       ! at this point:
       q_ist = q_prim
       call expand(q_prim,q,q_new)
       ! q_prim =  full (symmrtry-expanded ) vector of coords as before the optimization step
       ! q = symmetry-reduced set of coords already including the optimization step
       ! q_new = full vector of coords including the optimization step
    else
       q_ist = q_old
       q_new = q
       ! q_ist = vector of delocalized coords as before the optimization step
       ! q_new = vector of delocalized coords including the optimization step
    endif

    scf : do while (.not.converged)
       loop = loop + 1
       ! Step (1)
       ASSERT(.not.converged)
!      if (.not.converged) then
       dx=zero
!      endif
       if(contract) q_ist=(q_ist+q_keep)*0.5_r8_kind
       do i=1,n_atoms+n_dummy
          start = 3*(i-1)+1
          do j=1,n_coor
             delta =q_new(j) - q_ist(j)
             if (zmat_coordinates) then
                if (s_local(j)%typ  == d_angle ) then
                   if (delta > pi ) delta = -two*pi+delta
                   if (delta < -pi ) delta = two*pi+delta
                   if (abs(abs(delta)-two*pi)<=small) delta = zero
                endif
             endif
             dx(i,1) = dx(i,1) + bmat_inv(start,j)*delta
             dx(i,2) = dx(i,2) + bmat_inv(start+1,j)*delta
             dx(i,3) = dx(i,3) + bmat_inv(start+2,j)*delta
          enddo
          atom(i)%x = atom(i)%x + dx(i,:)
!          write(OPT_STDOUT,1012)i,dx(i,1),dx(i,2),dx(i,3)
       enddo
!1012   format ('Atom',i3,'dx(1..3)  :',3(f12.8,2x))
       ! --------- end if step (1) ------------

       ! Step (2)
       if(contract) then
          contract=.false.
          if (print_debug) print *,'restarting after ',max_iteration,'cycles'
       end if
       if (print_debug) print*,' ----------- LOOP ',loop,' --------------------'

       q_keep=q_ist
       if (delocalized_coordinates) then
          bmat_prim=zero
          bmat = zero
          call generate_bmat(s_prim,n_primitive,bmat_prim)
!         call generate_bmat_reduced(umat_flag=.false.)

          ! reduce bmat_prim into bmat, umat_trans must have been generated before:
          bmat = matmul(umat_trans,bmat_prim)

          call set_internal(q_prim,n_primitive,s_prim,atom) !(4)
          do i=1,n_primitive
             if(s_prim(i)%typ==d_angle) then
                if (q_prim(i)-q_prim_old(i) > pi ) then
                   delta = q_prim(i)-q_prim_old(i) - two*pi
                   q_prim(i)=q_prim_old(i)+delta
                endif
                if (q_prim(i)-q_prim_old(i) < -pi ) then
                   delta = two*pi + (q_prim(i)-q_prim_old(i))
                   q_prim(i)=q_prim_old(i)+delta
                endif
!!$                if (abs(q_prim(i)-q_prim_old(i)) == pi ) then
!!$                   delta=zero
!!$                   q_prim(i)=q_prim_old(i)
!!$                endif
             endif
          enddo
          s_prim(1:n_primitive)%value = q_prim
          q_ist = matmul(umat_trans,q_prim)
       elseif ( zmat_coordinates) then
          call generate_bmat(s_local,n_coor,bmat)
          call set_internal(q_ist,n_coor,s_local,atom) ! (5)
       endif

       call invert_bmat(bmat,bmat_inv,.false.) !!!(1)

       ! Step (3)+(4)

       ! test the difference
       sumi=0
       if (print_debug) &
            write(OPT_STDOUT,*)" Difference for Internals with respect to previous cycle"
       delta_max=zero
       do i=1,n_coor
          delta = q_new(i)-q_ist(i)
          if (print_debug ) &
               write(OPT_STDOUT,'(2(3x,f20.12),3x,es10.3,4x,i5)')q_new(i),q_ist(i),q_new(i)-q_ist(i),i
          if (zmat_coordinates) then
             if (s_local(i)%typ == d_angle ) then
                if (delta > pi ) delta = -two*pi+delta
                if (delta < -pi ) delta = two*pi+delta
             endif
          endif
          ! use small*0.01 as criterion to make sure that all geometrical
          ! constraints are fullfilled in the next gemetrie with a accuracy higher
          ! then small
          !mdf>>> (am)
          delta = abs(delta)
!!$          if ( q_new(i)==zero .or. abs(delta)<=small*0.01_r8_kind ) then
!!$             delta = abs(delta)
!!$          else
!!$             delta = abs(delta/q_new(i))
!!$          endif
          !<<<mdf
          if (delta>delta_max) &
               delta_max=delta

          if ( delta <= small*0.1_r8_kind ) then
             sumi=sumi+1
          endif
       enddo
       if (sumi.eq.n_coor) then
          converged=.true.
          ! set the internal variables a last time
          if (delocalized_coordinates) then
             call set_internal(q_prim,n_primitive,s_prim,atom) !(6)
          elseif (zmat_coordinates) then
             call set_internal(q_ist,n_coor,s_local,atom)   !(7)
          endif
       else
          converged=.false.
       endif

       if (loop == max_iteration ) then
          if(contract_counter==3) then
             write(OPT_STDOUT,*)"  internal_to_cart: no self-consistent transformation to cartesians"
             write(OPT_STDOUT,*)"                    could be found"
             write(OPT_STDOUT,'("  maximal deviation:     ",4x,2es10.3)')delta_max,small*0.1_r8_kind
             write(OPT_STDOUT,*)" "
             ABORT('no self-cons trafo found, see tty')
          else
             loop = 0 ! which means try once again
             contract_counter = contract_counter + 1
             contract = .true.
          end if
       endif


    enddo scf
    ! map the converged geometry to the x(:),y(:) and z(:)
    ! coordinates ... it is nicer for the output routine

    x(1:(n_atoms+n_dummy)) = atom(:)%x(1)
    y(1:(n_atoms+n_dummy)) = atom(:)%x(2)
    z(1:(n_atoms+n_dummy)) = atom(:)%x(3)

    ! calculate geometrical center again
    call geo_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),(n_atoms+n_dummy),xgeo)
    call mass_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),atom(:)%mass,(n_atoms+n_dummy),xmass,rotaxes) !(4) in intern to cart

    if(mass_center_print) then
      write(OPT_STDOUT,'("Geometrical center AFTER update  ",3(2x,f9.4))')(xgeo(i),i=1,3)
      write(OPT_STDOUT,'("Mass center AFTER update         ",3(2x,f9.4))')(xmass(i),i=1,3)
      write(OPT_STDOUT,*)" "
    endif

    s_local(:)%value = q_ist

    if(print_updated_geometry) call print_geometry()

    nullify(s_local)
    deallocate(q_ist,dx,q_new,q_keep,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" internal_to_cart: deallocation (0) failed")
    if(delocalized_coordinates) then

       deallocate(q_prim_old,STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            ("internal_to_cart: deallocation (1) failed")
    endif

!  if(present(bmat_inv_d)) then
!    bmat_inv_d=matmul(bmat_inv,transpose(reduc_mat))
!  endif

!1102 format(3(2x,f13.7))
!1020 format(' Interne Koord. No.',i4,'   Ist-Wert :',f12.8,3x,' Soll_wert :',f12.8)
  end subroutine internal_to_cart


  subroutine cart_to_internal(converged)
    ! Purpose: transform cartesian displacement coordinates
    !          atom(i)%x-atom(i)%x_old to  displacements
    !          of internal coordinates dq(j):
    !          dq = bmat*(atom(:)%x-atom(:)%x_old).
    !          These are compared
    !          to the actual displacements q(j)-q_old(j). If
    !          the difference is sufficiently small the value
    !          of 'converged' is set to TRUE and FALSE
    !          otherwise.
    !
    ! subroutine called by: internal_to_cart
    !
    !-----  Declaration of formal parameters --------------------
    logical,intent(inout) :: converged
    !----- Declaration of local variables -----------------------
    real(kind=r8_kind), allocatable :: q_test(:),q_test_p(:)
    integer(kind=i4_kind)    :: alloc_stat, i, sumi
    real(kind=r8_kind)       :: delta
    !---- executable code ---------------------------------------

    allocate(q_test(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop ' cart_to_internal : allocation (1a) failed'
    q_test=zero

    if (zmat_coordinates) then
       call set_internal(q_test,n_internal,s,atom)
    elseif (delocalized_coordinates) then
       allocate(q_test_p(n_primitive),STAT=alloc_stat)
       if (alloc_stat/=0) &
            stop ' cart_to_internal : allocation (1b) failed'
       q_test_p=zero
       call set_internal(q_test_p,n_primitive,s_prim,atom)
       q_test = matmul(umat_trans,q_test_p)
       deallocate(q_test_p,STAT=alloc_stat)
       if (alloc_stat/=0) &
            stop 'internal_to_cart: deallocation (1) failed'
    endif

    ! test the difference
    sumi=0
    do i=1,n_internal
       delta = abs((q_test(i)-q(i))/q(i))
       if ( delta <= 0.0005 ) then
          sumi=sumi+1
       endif
    enddo
    if (sumi.eq.n_internal) then
       converged=.true.
    else
       converged=.false.
    endif

    deallocate(q_test,STAT=alloc_stat)
    if (alloc_stat/=0) &
         stop ' cart_to_internal : deallocation (1) failed'
  end subroutine cart_to_internal

  subroutine cart2cart()
    ! ---- Modules ----------------------------------------------
#ifdef WITH_EFP
    use efp_module, only: update_efps
    use iounitadmin_module, only: returnclose_iounit
    use pointcharge_module, only: rcm
#endif
    !-----  Declaration of formal parameters --------------------
    !----- Declaration of local variables -----------------------
    integer(i4_kind) :: i
    real(r8_kind)    :: xgeo(3),xmass(3),rotaxes(3,3)
#ifdef WITH_EFP
    integer(i4_kind) :: alloc_stat, j, l
    real(r8_kind)    :: drot(3),dtrans(3)
#endif
    !---- executable code ---------------------------------------

    call  geo_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),(n_atoms+n_dummy),xgeo)
    call mass_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),atom(:)%mass,(n_atoms+n_dummy), &
         xmass,rotaxes)

    if(mass_center_print) then
       write(OPT_STDOUT,'("Geometrical center BEFORE update ",3(2x,f9.4))')(xgeo(i),i=1,3)
       write(OPT_STDOUT,'("Mass center BEFORE update        ",3(2x,f9.4))')(xmass(i),i=1,3)
       write(OPT_STDOUT,*)" "
    endif

#ifdef WITH_EFP
    if(qm_fixed) goto 100
#endif

    do i=1,n_atoms
       atom(i)%x(1)=q(3*i-2)
       atom(i)%x(2)=q(3*i-1)
       atom(i)%x(3)=q(3*i)
    end do

100 do i=1,n_atoms
       x(i) = atom(i)%x(1)
       y(i) = atom(i)%x(2)
       z(i) = atom(i)%x(3)
    end do

#ifdef WITH_EFP
    if(efp .and. n_efp > 0 .and. .not.efp_fixed) then
       l=n_atoms+1
       if(qm_fixed) l=1
       do i=1,n_efp
          dtrans(1)=q(3*l-2)-q_old(3*l-2)!rcm(1,i)
          dtrans(2)=q(3*l-1)-q_old(3*l-1)!rcm(2,i)
          dtrans(3)=q(3*l)  -q_old(3*l)!rcm(3,i)

          q(3*l-2)=dtrans(1)
          q(3*l-1)=dtrans(2)
          q(3*l)  =dtrans(3)
          l=l+1

          drot(1)=q(3*l-2)-q_old(3*l-2)!xyz_torque(i,1)
          drot(2)=q(3*l-1)-q_old(3*l-1)!xyz_torque(i,2)
          drot(3)=q(3*l)  -q_old(3*l)!xyz_torque(i,3)

          xyz_torque(i,1)=q(3*l-2)
          xyz_torque(i,2)=q(3*l-1)
          xyz_torque(i,3)=q(3*l)
          write(i_tq,*) xyz_torque(i,:)

          q(3*l-2)=drot(1)
          q(3*l-1)=drot(2)
          q(3*l)  =drot(3)
          l=l+1
       end do
       call returnclose_iounit(i_tq)

       deallocate(xyz_torque,stat=alloc_stat)
       ASSERT(alloc_stat==0)

       j=n_atoms
       if(qm_fixed) j=0
       call update_efps(x(n_atoms+1:n_atoms+3*n_efp), &
                        y(n_atoms+1:n_atoms+3*n_efp), &
                        z(n_atoms+1:n_atoms+3*n_efp), &
                        q(3*j+1:3*j+6*n_efp))
    end if
#endif

    ! calculate geometrical center again
    call  geo_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),(n_atoms+n_dummy),xgeo)
    call mass_center(atom(:)%x(1),atom(:)%x(2),atom(:)%x(3),atom(:)%mass,(n_atoms+n_dummy), &
         xmass,rotaxes)

    if(mass_center_print) then
       write(OPT_STDOUT,'("Geometrical center AFTER update  ",3(2x,f9.4))')(xgeo(i),i=1,3)
       write(OPT_STDOUT,'("Mass center AFTER update         ",3(2x,f9.4))')(xmass(i),i=1,3)
       write(OPT_STDOUT,*)" "
    endif
  end subroutine cart2cart

  !*************************************************************

  subroutine sym_check(step)
    ! Purpose: check and preserve the symmetry of the actual step.
    !          The data structures s(:)%unique and s(:)%equal are
    !          used to preserve the symmetry specified in the
    !          gxfile-format. Actually the update procedure should
    !          preserve it anyway since the gradients are
    !          symmetric.
    !
    ! Subroutine called by: newton_step,line_search_main
    !------------ Modules used ---------------------------------
    !------------ Declaration of formal parameters -------------
    real(kind=r8_kind),intent(inout)   :: step(:)
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind)   :: i,l_eq,counter
    ! test symmetry of step using the s-structure

    counter=0
    do i=1,n_internal
       if (s(i)%unique) then
          l_eq=1
          do while (s(i)%equal(l_eq)/=0)
             if ( (step(i)-step(s(i)%equal(l_eq))).gt.small ) then
                write(OPT_STDOUT,*)" Difference in symmetry for internal coord.",i
                write(OPT_STDOUT,1000)step(i),step(s(i)%equal(l_eq))
                counter=counter+1
             endif
             ! PRESERVE Symmetry at any cost
             step(s(i)%equal(l_eq)) = step(i)
             l_eq=l_eq+1
          enddo
       endif
    enddo
    if (counter/=0) then
       write(OPT_STDOUT,*)"newton_step : inconsistency in symmetry found"
    endif
    write(OPT_STDOUT,1100)(step(i),i=1,n_internal)
    write(OPT_STDOUT,*)" "
    write(OPT_STDOUT,*)" "
1000 format(3x,f14.11,5x,f14.11)
1100 format("  Step :",3x,7(f10.7,3x),//,8(f10.7,3x),//20(f10.7,3x))
  end subroutine sym_check
  !*************************************************************
  subroutine check_bmat(bmat_approx)
    ! generate Bmatrix using small variations in cartesian
    ! geometry.
    ! This routine is only intended as a numerical check for the
    ! *analytical* generation of the B-Matrix.
    ! ------------------------------------------------------
    use geo_operations_module
    real(kind=r8_kind)    :: bmat_approx(:,:)
    real(kind=r8_kind),allocatable    :: q_local(:),q_ref(:)
    integer(kind=i4_kind)             :: i,j,n,nn,ncol
    real(kind=r8_kind),parameter      :: delta=1.0e-5_r8_kind
    real(kind=r8_kind)                :: keep,signum

    allocate(q_local(n_internal),q_ref(n_internal))
    q_local=zero
    q_ref=zero
    call set_internal(q_ref,n_internal,s,atom)
    do j=1,n_internal
       i=0
       do n=1,n_atoms+n_dummy
          do nn=1,3
             i=i+1
             keep=atom(n)%x(nn)
             atom(n)%x(nn) = atom(n)%x(nn) + delta
             q_local = zero
             signum=one
             call set_internal(q_local,n_internal,s,atom)
             ! for the dihedral angle there is a discontinuity pi since
             ! it is defined between -pi and pi.
             if (abs(pi-q_ref(j))<=small) then !q_ref=pi within  accuracy
                if (q_local(j)<zero) then
                   q_local(j) = -q_local(j)
                   signum=-one
                endif
             elseif (abs(pi+q_ref(j))<=small) then ! q_ref=-pi within accuracy
                if (q_local(j)>zero) then
                   q_local(j)=-q_local(j)
                   signum=-one
                endif
             endif
             bmat_approx(j,i) = signum*(q_local(j)-q_ref(j))/delta
             atom(n)%x(nn) = keep
          enddo
       enddo! atoms+dummies
    enddo! internals

    ncol=10_i4_kind
    write(OPT_STDOUT,*)" --- approximate B-Matrix -----------------------"
    call print_matrix(bmat,n_internal,3*(n_atoms+n_dummy),ncol)
    write(OPT_STDOUT,*)" "
    deallocate(q_ref,q_local)
  end subroutine check_bmat

  ! ***************************************************************************
  subroutine generate_umat(bmat,umat)
    ! Purpose: generate the matrix for reducing the B-Matrix by
    !          (1) constructing B*Bt
    !          (2) diagonalizing B*Bt
    !          (3) sort out those eigenvectors belonging to
    !              non-zero eigenvalues: active set U
    !          (4) bmat_red = Ut*b_mat
    !
    !          The matrix containing the linear combination coefficients
    !          of the primitive internals (umat) is biult up.
    !          For instance for the backtransformation from internal
    !          to cartesian coordinates.
    !
    ! Subroutine called by: coordinates_setup,internal_to_cart
    ! --------------------------------------------------------------------
    use readwriteblocked_module
    implicit none
!   logical,intent(in)              :: umat_flag
    real(r8_kind), intent(in)  :: bmat(:,:) ! primitive B-matrix (bmat_prim?)
    real(r8_kind), intent(out) :: umat(:,:)
    ! *** end of interface ***

    ! ---- Declaration of local variables --------------------------------
    ! this logical switch determines, if previous eigenvectors are
    ! available to adjust the signs of the current eigenvectors
    logical,save                    :: eigvec_init=.false.
    ! ----------------------------------------------------
    real(kind=r8_kind),allocatable    :: eigval(:),eigvec(:,:),&
         help_mat(:,:),eigvec_prev(:,:),eigval_prev(:),umat_tmp(:,:)
    integer(kind=i4_kind)             :: alloc_stat,i,n_max,counter,j,dimen
    integer(kind=i4_kind),parameter   :: n_col = 10, io_deloc = 25
    logical                           :: exist = .false.,stop_flag=.false.
    real(kind=r8_kind)                :: check,dimen_real(1)
    type(readwriteblocked_tapehandle) :: th   ! tape_handle for 'delocalized.dat'

!   if ( .not.umat_flag ) then
!      bmat = matmul(umat_trans,bmat_prim)
!      return
!   else

       inquire(EXIST=exist,FILE=trim(opt_data_dir)//'/delocalized.dat')
       if (print_debug) then
          if (exist) write(OPT_STDOUT,*)" Taking previous eigenvectors from file"
          if (.not.exist) write(OPT_STDOUT,*)" No file containing previous eigenvectors found"
       endif
       if (update_delocalized) exist =.false.

       allocate(eigvec(n_primitive,n_primitive),eigval(n_primitive),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            (" generate_umat: allocation (1) failed")
       if (n_constraint>0_i4_kind) then
          allocate(umat_tmp(n_primitive,n_internal+n_constraint),STAT=alloc_stat)
          if (alloc_stat/=0 ) call error_handler&
               (" generate_umat: allocation (2) failed")
       endif

       if ( exist ) then
          allocate(eigvec_prev(n_primitive,n_primitive),&
                  eigval_prev(n_primitive),STAT=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               (" generate_umat: allocation (3) failed")
          call readwriteblocked_startread(trim(opt_data_dir)//'/delocalized.dat',th)
          call readwriteblocked_read(dimen_real,th)
          dimen=int(dimen_real(1),kind=i4_kind)

          if (dimen /= n_primitive ) call error_handler&
               (" generate_umat: eigvecs on delocalized.dat do not have right dimension")
          call readwriteblocked_read(eigval_prev(:),th)
          do i=1,dimen
             call readwriteblocked_read(eigvec_prev(:,i),th)
          enddo
          call readwriteblocked_stopread(th)
          eigval = eigval_prev
          eigvec = eigvec_prev
          eigvec_init = .true.
       endif

       allocate(help_mat(n_primitive,n_primitive),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            (" generate_umat: allocation (2) failed")
       help_mat = matmul(bmat,transpose(bmat))

       if (.not. eigvec_init ) then
          eigvec=zero
          eigval=zero
          call eigensolver(help_mat,n_primitive,eigval,eigvec)
       elseif ( eigvec_init .and. update_delocalized ) then
          eigvec=zero
          eigval=zero
          call eigensolver(help_mat,n_primitive,eigval,eigvec)
          call eigvec_sign_check(n_primitive,eigvec,eigvec_prev)
       endif

       if (print_debug) &
            write(OPT_STDOUT,*)" generate_umat: writing eigenvecs to file delocalized.dat"

       call readwriteblocked_startwrite(trim(opt_data_dir)//'/delocalized.dat',th)
       dimen_real(1) = real(n_primitive,kind=r8_kind)
       call readwriteblocked_write(dimen_real,th)
       call readwriteblocked_write(eigval(:),th)
       do i=1,n_primitive
          call readwriteblocked_write(eigvec(:,i),th)
       enddo
       call readwriteblocked_stopwrite(th)

       ! sort out the active subset
       n_max = n_internal
       counter=1_i4_kind
       umat=zero
       do i=1,n_primitive
          if (eigval(i)>=small) then
             ! this belongs to the active subset
             if (n_constraint>0_i4_kind) then
                umat_tmp(:,counter) = eigvec(:,i)
             else
                umat(:,counter) = eigvec(:,i)
             endif
             check=zero
             do j=1,n_primitive
                check=check+abs(eigvec(j,i))
             enddo
             if(check<=small) then
                write(OPT_STDOUT,*)"     DEFINITION OF INTERNAL COORDINATES  "
                write(OPT_STDOUT,*)"     DOES NOT SPAN AN 3*N-6 SUBSPACE     "
                write(OPT_STDOUT,*)"     REDEFINE YOUR INTERNALS             "
                stop_flag=.true.
             endif
             counter=counter+1_i4_kind
          else
!             if (print_debug) then
!                write(OPT_STDOUT,*)"generate_umat: The ",i,"th Eigenvalue is equal to zero"
!                write(OPT_STDOUT,*)"                       The corresponding eigenvector is "
!                do j=1,n_primitive
!                   write(OPT_STDOUT,'("component ",i4,"      ",f9.5)')j,eigvec(j,i)
!                enddo
!             endif
          endif
       enddo
       counter=counter-1_i4_kind
       write(OPT_STDOUT,*)" generate_umat: The active subset contains ",counter," eigenvectors"
       if (counter/=n_max+n_constraint) then
          write(OPT_STDOUT,*)" generate_umat: wrong number of  eigenvectors in the active subset"
          stop_flag=.true.
       endif
       if (n_constraint>0_i4_kind) then
          if (print_debug) then
             print*,' Old UMAT '
             call print_matrix(umat_tmp,n_primitive,n_internal,10_i4_kind)
             print*,' '
          endif
! This is an older version of the same procedure. Actually it is not
! necessary any more.
!          call schmidt(n_constraint,constraint,umat_tmp,umat)
          call ortho2(constraint,umat_tmp,umat)
!          if (.not.stop_flag) &
!               call ortho_constraint(n_constraint,constraint,umat_tmp,umat)
          deallocate(umat_tmp,STAT=alloc_stat)
          if (alloc_stat/=0) &
               stop ' generate_umat: deallocation (1) failed'
       endif

!      umat_trans = transpose(umat)

       if (print_debug) then
          write(OPT_STDOUT,*)" This is the matrix containing the eigenvectors of the active subset"
          call print_matrix(umat,n_primitive,n_max,n_col)
          write(OPT_STDOUT,*)" "
       endif

       if (stop_flag) then
          stop 1
       endif
!      bmat = matmul(umat_trans,bmat_prim)

       deallocate(help_mat,eigval,eigvec,STAT=alloc_stat)
       ASSERT(alloc_stat==0)

       if (exist ) then
          deallocate(eigvec_prev,eigval_prev,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
       endif
!   endif
  end subroutine generate_umat
  ! ***************************************************************************
  subroutine alloc_bmat(n_coor)
    ! Purpose: just a brief wrapper for the few lines necessary to allocate
    !          the B-Matrix 'bmat'
    ! ----------------------------------------------------------------------
    integer(kind=i4_kind),intent(in)   :: n_coor
    integer(kind=i4_kind) :: alloc_stat


    if(.not.allocated(bmat)) then
    allocate(bmat(n_coor,3*(n_atoms+n_dummy)),STAT=allocopt_stat(12))
    ASSERT(allocopt_stat(12).eq.0)
    endif
    allocate(bmat_inv(3*(n_atoms+n_dummy),n_coor),STAT=allocopt_stat(10))
    ASSERT(allocopt_stat(10).eq.0)
    bmat=zero
    bmat_inv = zero

    allocate(bmat_prim(n_primitive,3*(n_atoms+n_dummy)),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    bmat_prim=zero
  end subroutine alloc_bmat
  ! ***************************************************************************
  subroutine free_bmat()
    ! Purpose: do the opposite of 'alloc_bmat'
    ! ------------------------------------------------------------------------
    integer(kind=i4_kind)    :: alloc_stat
    deallocate(bmat,STAT=allocopt_stat(12))
    ASSERT(allocopt_stat(12).eq.0)
    allocopt_stat(12)=1
    deallocate(bmat_inv,STAT=allocopt_stat(10))
    ASSERT(allocopt_stat(10).eq.0)
    allocopt_stat(10)=1

    deallocate(bmat_prim,STAT=alloc_stat)
    if ( alloc_stat/=0) then
       stop ' alloc_bmat : deallocation of bmat_inv failed'
    endif
  end subroutine free_bmat

  ! ***************************************************************************

  subroutine cart_constraint_setup(n_coor)
    ! Purpose:  the matrix that eliminates constant coordinates from the
    !          optimization step.
    ! -----------------------------------------------------------------------
    implicit none
    integer(kind=i4_kind),intent(in)   :: n_coor
    ! *** end of interface ***

    integer(kind=i4_kind) :: alloc_stat, i

    n_constraint=0

    n_internal = n_internal-n_constraint
    allocate(tmat(n_internal,n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    tmat = zero

    do i=1,n_internal
       tmat(i,i) = one
    enddo
  end subroutine cart_constraint_setup

  subroutine constraint_setup(coor,n_coor)
    ! Purpose: set the constraint vectors (delocalized coordinates only) and
    !          the matrix that eliminates constant internals from the
    !          optimization step.
    ! -----------------------------------------------------------------------
    type(int_coor),intent(in)          :: coor(:)
    integer(kind=i4_kind),intent(in)   :: n_coor
    ! --- Declaration of local variables ------------------------------------
    integer(kind=i4_kind) :: alloc_stat,i,c_count

    ASSERT(n_coor==size(coor))
    n_constraint=0
    do i=1,n_coor
       if (.not.coor(i)%var) then
          n_constraint=n_constraint+1
       endif
    enddo

    n_internal = n_internal-n_constraint
    allocate(tmat(n_internal,n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    tmat = zero

    do i=1,n_internal
       tmat(i,i) = one
    enddo

    if (delocalized_coordinates) then
       allocate(constraint(n_coor,n_constraint),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       constraint = zero
       c_count=1
       do i=1,n_coor
          if (.not.coor(i)%var) then
             constraint(i,c_count) = one
             ! tmat(c_count,c_count) = zero
             c_count = c_count+1
          endif
       enddo
    elseif ( zmat_coordinates) then
       if (n_coor/=n_internal) then
          write(OPT_STDOUT,*)" constraint_setup: Number of internals for Z-Matrix Coordinates wrong"
          stop 1
       endif
       do i=1,n_internal
          if (.not.coor(i)%var ) tmat(i,i) = zero
       enddo
    endif
  end subroutine constraint_setup

  subroutine dealloc_constraint()
integer(kind=i4_kind) :: alloc_stat
    if(allocated(tmat)) then
       deallocate(tmat,STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    end if
    if (delocalized_coordinates) then
       deallocate(constraint,STAT=alloc_stat)
       if (alloc_stat/=0) &
            stop 'constraint_setup: deallocation (1) failedd'
    endif
  end subroutine dealloc_constraint

  subroutine alloc_intcoor(s_prim,q_prim,s,sym_type)

    ! This is a brief wrapper for the allocation of the (primitive)
    ! internal coordinates. We habe to distinguish the following cases:
    ! (1) Z-Matrix coordinates are used
    !     We do not have to allocate 's_prim', but can directly
    !     proceed to allocate s(3*N-6=n_internal).
    ! (2) Delocalized internal coordinates are used
    !     We need to allocate s_prim(n_primitive) for the primitive
    !     internals. n_primitive is either
    !     n_internal : if zmat-connectivities are taken
    !     n_primitive: if freeformat-connectivities are taken
    !     n_valence  : if valence-format connectivities are used.
    !     Then we need to allocate s(n_internal=3*N-6).
    ! ***************************************************************************
    type(int_coor),optional, pointer:: s_prim(:)
    type(int_coor), optional, pointer:: s(:)
    real(r8_kind),optional, pointer:: q_prim(:)
    integer(kind=i4_kind),optional, pointer   :: sym_type(:)

    integer(kind=i4_kind)   :: alloc_stat,num,i

    if (zmat_format.or.free_format) then
       num=n_primitive
    else
       num=int(binomial_coeff(n_atoms+n_dummy,2),i4_kind) + &
            3*int(binomial_coeff(n_atoms+n_dummy,3),i4_kind) + &
            6* int(binomial_coeff(n_atoms+n_dummy,4),i4_kind)
       ! this is the maximum number of bond stretches,
       ! bond angles and dihedral angles
       write(OPT_STDOUT,*)" alloc_intcoor: maximum number of variables = ",num
       num=n_primitive
       write(OPT_STDOUT,*)"                actual number of variables  = ",num
       write(OPT_STDOUT,*)" "
    endif

   if(present(s_prim)) then
    allocate(s_prim(num),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    ! retain this - it may become a fossil soon
    do i=1,num
       allocate(s_prim(i)%equal(max_equal),&
                s_prim(i)%equal_sign(max_equal),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       s_prim(i)%equal = 0_i4_kind
       s_prim(i)%n_equal = 1_i4_kind
       s_prim(i)%equal_sign = 0_i4_kind
       s_prim(i)%value = zero
       s_prim(i)%var = .false.
       s_prim(i)%typ = 0_i4_kind
       s_prim(i)%sphere= .false.
    enddo
   endif

   if(present(sym_type)) then
    allocate(sym_type(num),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    sym_type = 0_i4_kind
   endif

    if(.not.delocalized_coordinates.and.present(s)) then
       allocate(s(n_internal),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
       s(:)%value = zero
       s(:)%typ = 0_i4_kind
    end if

   if(present(q_prim)) then
    allocate(q_prim(n_primitive),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    q_prim = zero
   endif

  end subroutine alloc_intcoor

  subroutine dealloc_intcoor(s_prim,q_prim,s,sym_type,q,q_old)
    ! This is a brief wrapper for the allocation of the (primitive)
    ! internal coordinates. We habe to distinguish the following cases:
    ! (1) Z-Matrix coordinates are used
    !     We do not have to allocate 's_prim', but can directly
    !     proceed to allocate s(3*N-6=n_internal).
    ! (2) Delocalized internal coordinates are used
    !     We need to allocate s_prim(n_primitive) for the primitive
    !     internals. n_primitive is either
    !     n_internal : if zmat-connectivities are taken
    !     n_primitive: if freeformat-connectivities are taken
    !     n_valence  : if valence-format connectivities are used.
    !     Then we need to allocate s(n_internal=3*N-6).
    ! ***************************************************************************
    type(int_coor), optional, pointer:: s_prim(:)
    type(int_coor), optional, pointer:: s(:)
    real(r8_kind), optional, pointer:: q_prim(:)
    real(r8_kind), pointer, optional:: q(:),q_old(:)
    integer(kind=i4_kind),optional, pointer   :: sym_type(:)

    integer(kind=i4_kind)   :: alloc_stat,num,i

    if (zmat_format.or.free_format) then
       num=n_primitive
    else
       num=int(binomial_coeff(n_atoms+n_dummy,2),i4_kind) + &
            3*int(binomial_coeff(n_atoms+n_dummy,3),i4_kind) + &
            6* int(binomial_coeff(n_atoms+n_dummy,4),i4_kind)
       ! this is the maximum number of bond stretches,
       ! bond angles and dihedral angles
       write(OPT_STDOUT,*)" alloc_intcoor: maximum number of variables = ",num
       num=n_primitive
       write(OPT_STDOUT,*)"                actual number of variables  = ",num
       write(OPT_STDOUT,*)" "
    endif

   if(present(s_prim)) then
    if(associated(s_prim) ) then
    DPRINT 'deallocate s_prim'
    do i=1,num
       deallocate(s_prim(i)%equal,s_prim(i)%equal_sign,STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    enddo
    deallocate(s_prim,STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    DPRINT 'done deallocate s_prim'
   endif
   endif

    if(present(sym_type)) then
     if(associated(sym_type)) then
     deallocate(sym_type,STAT=alloc_stat)
     ASSERT(alloc_stat.eq.0)
    endif
    endif

    if(.not.delocalized_coordinates.and.present(s)) then
    if(associated(s)) then
       deallocate(s,STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    endif
    end if

   if(present(q_prim)) then
   if(associated(q_prim)) then
    deallocate(q_prim,STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
   endif
   endif

  if(present(q)) then
   if(associated(q)) then
    deallocate(q_old,q,STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
  endif
  endif

  end subroutine dealloc_intcoor

  ! ***************************************************************************

  subroutine alloc_cartcoor()
    implicit none
    ! *** end of interface ***

    integer(kind=i4_kind) :: alloc_stat

    allocate(q_prim(n_primitive),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    q_prim = zero
  end subroutine alloc_cartcoor

  subroutine dealloc_cartcoor()
#ifdef WITH_EFP
    use pointcharge_module, only: rcm
#endif
    implicit none
    integer(kind=i4_kind) :: alloc_stat
    ! *** end of interface ***

    if(associated(q_old)) then
       deallocate(q_prim,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(q_prim)

       deallocate(q_old,q,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
       nullify(q_old)
    end if

#ifdef WITH_EFP
    if(allocated(rcm) .and. .not. qm_fixed) then
       deallocate(rcm,stat=alloc_stat)
       ASSERT(alloc_stat==0)
    end if
#endif
  end subroutine dealloc_cartcoor

  ! *********************************************************************

  subroutine reduce_setup(geo_loop)
    ! Purpose: construct the matrices 'reduc' and 'expand'.
    !          q = reduc * q_prim
    !          q_prim = expand * q + q' where q' contains
    !          the original coordinate set for those
    !          components which are not varied.
    !          for the case of zmat_coordinates specified in zmat_format.
    ! --------------------------------------------------------------
    implicit none
    integer(kind=i4_kind), optional :: geo_loop ! during
    ! the calculation of the hessian it is possible, that two
    ! symmetryequivalent coordinates (number geo_loop) have different
    ! values (antisymmetric modes!)
    ! --------------------------------------------------------------
    integer(kind=i4_kind)   :: alloc_stat,i,j,int_count,i_equal
    real(kind=r8_kind) :: ivalue,jvalue

!   print *,'reduce_setup entered'
    int_count=0
    do i=1,n_primitive
       if ( .not.s_prim(i)%unique) cycle ! i primitive
          int_count=int_count+1
          equals: do i_equal=2,s_prim(i)%n_equal

             j = s_prim(i)%equal(i_equal)

             ivalue = s_prim(i)%value
             jvalue = s_prim(j)%value

             if (abs(ivalue + jvalue)<=small .and. abs(ivalue)>small ) then
                if (present(geo_loop)) then
                   if(geo_loop==sym_type(i)) then
                      if(abs(abs(ivalue - jvalue) - 2.0*step_size)>small) then
                         write(OPT_STDOUT,*)" reduce_setup: Unique coordinate ",int_count," and symmetry-equivalent"
                         write(OPT_STDOUT,*)"               coordinate ",j," have opposite sign "
                         write(OPT_STDOUT,*)"               They will be considered as symmetry-equivalent"
                         s_prim(i)%equal_sign(i_equal) = -1_i4_kind*s_prim(i)%equal_sign(i_equal)
!                         write(OPT_STDOUT,*) s_prim(i)%equal_sign(i_equal)
                      end if
                   else
                      write(OPT_STDOUT,*)" reduce_setup: Unique coordinate ",int_count," and symmetry-equivalent"
                      write(OPT_STDOUT,*)"               coordinate ",j," have opposite sign "
                      write(OPT_STDOUT,*)"               They will be considered as symmetry-equivalent"
                      s_prim(i)%equal_sign(i_equal) = -1_i4_kind*s_prim(i)%equal_sign(i_equal)
!                         write(OPT_STDOUT,*) s_prim(i)%equal_sign(i_equal)
                   endif
                else
                   write(OPT_STDOUT,*)" reduce_setup: Unique coordinate ",int_count," and symmetry-equivalent"
                   write(OPT_STDOUT,*)"               coordinate ",j," have opposite sign "
                   write(OPT_STDOUT,*)"               They will be considered as symmetry-equivalent"
                   s_prim(i)%equal_sign(i_equal) = -1_i4_kind*s_prim(i)%equal_sign(i_equal)
!                         write(OPT_STDOUT,*) s_prim(i)%equal_sign(i_equal)
                end if
             elseif (abs(ivalue - jvalue)>small ) then
                if (present(geo_loop)) then
                   if(geo_loop==sym_type(i)) then

                      if(abs(abs(ivalue + jvalue)-2.0*step_size)<=small &
                           .and.abs(ivalue - jvalue)>small) then
                         write(OPT_STDOUT,*)" reduce_setup: Unique coordinate ",int_count," and symmetry-equivalent"
                         write(OPT_STDOUT,*)"               coordinate ",j," have opposite sign "
                         write(OPT_STDOUT,*)"               They will be considered as symmetry-equivalent"
                         s_prim(i)%equal_sign(i_equal) = -1_i4_kind*s_prim(i)%equal_sign(i_equal)
                         write(OPT_STDOUT,*) s_prim(i)%equal_sign(i_equal)

                      elseif(abs(abs(ivalue - jvalue)-2.0*step_size)>small) then
                         write(OPT_STDOUT,*)" reduce_setup: Unique Coordinate ",int_count," Value ",ivalue
                         write(OPT_STDOUT,*)"               Value of wrong coordinate   ",jvalue
                         write(OPT_STDOUT,*)"               Coordinate declared as symmetry-equivalent does not"
                         write(OPT_STDOUT,*)"               have the same value"
                         call error_handler(" --- Please correct your gxfile (2.1), see flepo --- ")
                      end if
                   else
                         write(OPT_STDOUT,*)" reduce_setup: Unique Coordinate ",int_count," Value ",ivalue
                         write(OPT_STDOUT,*)"               Value of wrong coordinate   ",jvalue
                         write(OPT_STDOUT,*)"               Coordinate declared as symmetry-equivalent does not"
                         write(OPT_STDOUT,*)"               have the same value"
                         write(OPT_STDOUT,*)"               primitive  =", i
                         write(OPT_STDOUT,*)"               equivalent =", j,'(',i_equal,'-th of',s_prim(i)%n_equal,')'
                         call error_handler(" --- Please correct your gxfile (2.2), see flepo --- ")
                   endif
                else
                   write(OPT_STDOUT,*)" reduce_setup: Unique Coordinate ",int_count," Value ",ivalue
                   write(OPT_STDOUT,*)"               Value of wrong coordinate   ",jvalue
                   write(OPT_STDOUT,*)"               Coordinate declared as symmetry-equivalent does not"
                   write(OPT_STDOUT,*)"               have the same value"
                   write(OPT_STDOUT,*)"               Coordinate definition from Z-Matrix: "

                   if (s_prim(j)%typ == b_length ) then
                      write(OPT_STDOUT,*)"         bond-length, Partner1, Partner2 = ",&
                           s_prim(s_prim(i)%equal(i_equal))%length%partner1,&
                           s_prim(s_prim(i)%equal(i_equal))%length%partner2

                   elseif(s_prim(j)%typ == b_angle ) then
                      write(OPT_STDOUT,*)"         bond-angle, Partner1,Apex,Partner2:",&
                           s_prim(s_prim(i)%equal(i_equal))%angle%partner1,&
                           s_prim(s_prim(i)%equal(i_equal))%angle%apex,&
                           s_prim(s_prim(i)%equal(i_equal))%angle%partner2

                   elseif(s_prim(j)%typ == d_angle ) then
                      write(OPT_STDOUT,*)"         torsion, Partner1,Base1,Bas2,Partner2",&
                           s_prim(s_prim(i)%equal(i_equal))%dihedral%partner1,&
                           s_prim(s_prim(i)%equal(i_equal))%dihedral%base1,&
                           s_prim(s_prim(i)%equal(i_equal))%dihedral%base2,&
                           s_prim(s_prim(i)%equal(i_equal))%dihedral%partner2
                   else
                      call error_handler("reduce : Types are messed up")
                   endif
                   call error_handler(" --- Please correct your gxfile (2.3), see flepo --- ")
                end if
             end if
          enddo equals
    enddo

    allocate(reduc_mat(n_internal,n_primitive),&
         expand_mat(n_primitive,n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    reduc_mat = zero
    expand_mat = zero
    int_count = 0_i4_kind

    do i=1,n_primitive
       if ( .not. s_prim(i)%unique .or. sym_type(i)==0 ) cycle
       int_count = int_count + 1
       if (int_count > n_internal ) call error_handler &
            ("too many internal coordinates in reduce_setup")

       do j=1,s_prim(i)%n_equal
          reduc_mat(int_count,s_prim(i)%equal(j)) = float(s_prim(i)%equal_sign(j))
       enddo
    enddo
    expand_mat = transpose(reduc_mat)
  end subroutine reduce_setup

  subroutine dealloc_reduce()
        if(allocated(reduc_mat)) then
        deallocate(reduc_mat,expand_mat,stat=allocopt_stat(31))
        ASSERT(allocopt_stat(31).eq.0)
        endif
  end subroutine dealloc_reduce

  subroutine reduce(s,s_prim,q,q_prim)
    ! Purpose: reduce the full set  of internal coordinates to those
    !          coordinates which are unique (-> s_prim(i)%unique )
    ! ------------------------------------------------------------
    implicit none
    type(int_coor), intent(in):: s_prim(:)
    type(int_coor), intent(out):: s(:)
    real(kind=r8_kind), intent(in):: q_prim(:)
    real(kind=r8_kind), intent(out):: q(:)
    ! *** end of interface ***

    integer(kind=i4_kind)   :: i, int_count, n_sym

    q = matmul(reduc_mat,q_prim)
    do i=1,n_internal
       n_sym = sum(abs(reduc_mat(i,:))) !take the absolute value because of
       if (n_sym == 0) call error_handler&
            ("reduce: n_sym=0. Seek technical assistance")
       q(i) = q(i)/real(n_sym,kind=r8_kind) ! symmetry-equivalent internals
    enddo!                                    with opposite sign.

    s(:)%value = q
    s(:)%var = .true.
    int_count=0_i4_kind
    do i=1,n_primitive
       if (s_prim(i)%unique) then
          int_count=int_count+1
          s(int_count)%typ = s_prim(i)%typ
          s(int_count)%sphere=s_prim(i)%sphere
          if (s_prim(i)%typ == b_length ) then
             s(int_count)%length%partner1 = s_prim(i)%length%partner1
             s(int_count)%length%partner2 = s_prim(i)%length%partner2
          elseif (s_prim(i)%typ == b_angle ) then
             s(int_count)%angle%partner1 = s_prim(i)%angle%partner1
             s(int_count)%angle%apex = s_prim(i)%angle%apex
             s(int_count)%angle%partner2 = s_prim(i)%angle%partner2
          elseif (s_prim(i)%typ == d_angle ) then
             s(int_count)%dihedral%partner1 = s_prim(i)%dihedral%partner1
             s(int_count)%dihedral%base1 = s_prim(i)%dihedral%base1
             s(int_count)%dihedral%base2 = s_prim(i)%dihedral%base2
             s(int_count)%dihedral%partner2 = s_prim(i)%dihedral%partner2
          else
             call error_handler (" sth. fishy when setting  types in reduce")
          endif
       endif
    enddo

  end subroutine reduce

  ! ******************************************************************

  subroutine expand(q_primitive,q_in,q_expand)
    ! Purpose: expand the coordinate vector again
    ! --------------------------------------------------------------
    real(kind=r8_kind),intent(in)   :: q_in(:),q_primitive(:)
    real(kind=r8_kind),intent(out)  :: q_expand(:)
    integer(kind=i4_kind)           :: i
    !---------------------  local variables  -----------------------
    real(kind=r8_kind)              :: q_reduc(size(q_in)), n_sym
    q_reduc=matmul(reduc_mat,q_primitive)
    do i=1,n_internal
       ! now take into account the multiplicity of the
       ! coordinates
       n_sym = sum(abs(reduc_mat(i,:)))  ! take the absolute value because of
       q_reduc(i) = q_reduc(i)/n_sym     ! symmetry-equivalent internals
    enddo                                ! with opposite sign.
    q_expand=q_primitive+matmul(expand_mat,(q_in-q_reduc))
  end subroutine expand

  ! ******************************************************************

end module coordinates_module
