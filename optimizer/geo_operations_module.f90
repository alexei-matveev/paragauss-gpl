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
module  geo_operations_module
  !---------------------------------------------------------------
  !
  !  Purpose: contains geometric operations such as
  !           translation and rotation of the whole molecule
  !
  !  Module called by: coordinates_module, ...
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use opt_data_module
  use math_module
  use atom_data_module
  use iounitadmin_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================
  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ public functions and subroutines ------------------
  public translate_molecule,xmol_output,rotate_molecule,&
       print_geometry,geo_center,mass_center,minerttensor,centerofmass

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of constants and variables ----
  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains
  !*************************************************************

  subroutine translate_molecule(vec,index_atom)
    ! Purpose: translate the whole molecule by moving atom
    !          'index_atom' to 'vec'.
    ! ----------------------------------------------------------
    ! ----- Declaration of formal parameters ------------------
    real(kind=r8_kind),intent(in)    :: vec(3)
    integer(kind=i4_kind),intent(in) :: index_atom
    ! ---- Declaration of local variables ---------------------
    real(kind=r8_kind)    :: dx(3)
    integer(kind=i4_kind) :: i
    ! ---- executable code ------------------------------------
    dx = vec-atom(index_atom)%x
    do i=1,(n_atoms+n_dummy)
       atom(i)%x = atom(i)%x + dx
    enddo
  end subroutine translate_molecule

  !*************************************************************

  subroutine rotate_molecule(axis,angle)
    ! Purpose: rotate the whole molecule about the
    !          desired 'axis' by angle 'angle'.
    !----------------------------------------------------
    use math_module
    integer(kind=i4_kind),intent(in)  :: axis
    ! 1 = x-axis, 2 = y-axis, 3 = z-axis
    real(kind=r8_kind)  :: angle
    ! Declaration of local variables -------
    real(kind=r8_kind)     :: vec1(3),vec_zero(3),help(3)
    real(kind=r8_kind)     :: m(3,3)
    integer(kind=i4_kind)  :: i, n
    ! first store the position of atom 1 and moce the
    ! molecule to the origin of the cartesian coordinate
    ! system:
    vec1 = atom(1)%x
    vec_zero = zero
    call translate_molecule(vec_zero,1)

    m = rotmat(axis,angle)
    do n=1,n_atoms+n_dummy
       help = zero
       do i=1,3
          help(i) = dot_product(m(i,:),atom(n)%x)
       enddo
       atom(n)%x = help
    enddo

    call translate_molecule(-vec1,1)

  end subroutine rotate_molecule
  !*************************************************************

  subroutine print_geometry()
    ! Purpose: preliminary output of actual geometry
    integer(kind=i4_kind)  :: i,j
    character(len=3)       :: sym
    write(OPT_STDOUT,*)'The updated geometry (au) is'
    do i=1,n_atoms+n_dummy
       sym = adjustl(trim(symbol(nint(atom(i)%charge))))
       write(OPT_STDOUT,1000)i,sym,(atom(i)%x(j),j=1,3)
    enddo
    write(OPT_STDOUT,*)" "
1000 format(i4,6x,a3,4x,3(2x,f13.9))
  end subroutine print_geometry

  !*************************************************************

  subroutine geo_center(x,y,z,n,xgeo)
    ! Purpose: calculate the geometrical center of mass of
    !          the molecule using:
    ! x_g = 1/(N_atoms+N_dummy) sum(i) [ x(i) ]
    !
    ! x,y,z(N_atoms+N_dummy) are the x,y and z-coordinates
    ! of all atoms, including dummy-atoms
    ! --------------------------------------------------------
    integer(kind=i4_kind),intent(in)          :: n
    real(kind=r8_kind),target,intent(in)      :: x(n),y(n),z(n)
    real(kind=r8_kind),intent(out)            :: xgeo(3)
    ! --- Declaration of local variables ---------------------
    real(kind=r8_kind),pointer                :: xpoint(:)
    integer(kind=i4_kind)                     :: comp,i
    integer                                   :: ierr
    real(kind=r8_kind)                        :: xtensor(3,3),rot(3,3),eigrot(3)
    ! geometrical tensor of rank 2, main rotation axes and rotation moments
    real(kind=r8_kind)                        :: fv1(3),fv2(3) ! help_arrays

    xgeo=zero

    do comp=1,3
       select case (comp)
       case (1)
          xpoint => x
       case (2)
          xpoint => y
       case (3)
          xpoint => z
       case default
          stop 'rubbish in routine GEO_CENTER'
       end select
       do i=1,n
          xgeo(comp) = xgeo(comp) + &
                       xpoint(i)
       enddo
    enddo
    xgeo = one/real(n,kind=r8_kind) * xgeo
    ! now we additionaly calculate the three main rotation axis, in order to
    ! check if they have changed.
    xtensor(1,1)=sum((x(1:n)-xgeo(1))**2)
    xtensor(2,2)=sum((y(1:n)-xgeo(2))**2)
    xtensor(3,3)=sum((z(1:n)-xgeo(3))**2)
    xtensor(2,1)=sum((x(1:n)-xgeo(1))*(y(1:n)-xgeo(2)))
    xtensor(3,1)=sum((x(1:n)-xgeo(1))*(z(1:n)-xgeo(3)))
    xtensor(3,2)=sum((y(1:n)-xgeo(2))*(z(1:n)-xgeo(3)))
    xtensor(1,2)=xtensor(2,1)
    xtensor(1,3)=xtensor(3,1)
    xtensor(2,3)=xtensor(3,2)
    call rs(3,3,xtensor,eigrot,1,rot,fv1,fv2,ierr)
    if(ierr/=0)  stop 'geo_center: diagonalization failed'
!VVP: I think this information is not necessary for ParaGauss's users.
!    write(OPT_STDOUT,*) 'Rotational moments and corresponding eigenvectors:'
!    write(OPT_STDOUT,*)"          moment         eigvec 1       eigvec 2       eigvec 3"
!    write(OPT_STDOUT,*)" --------------------------------------------------------------"
!    do comp=1,3
!       write(OPT_STDOUT,'("     ",4F15.10)') eigrot(comp),rot(comp,:)
!    end do
!    write(OPT_STDOUT,*)" --------------------------------------------------------------"
!    write(OPT_STDOUT,*)" "
  end subroutine geo_center

  !*************************************************************
  subroutine mass_center(x,y,z,mass,n,xmass,rotaxes)
    ! Purpose: calculate the center of mass of the molecule
    !          using:
    !          xmass = 1/(Mges) sum(i) [ Mi Xi ]
    !          Input parameters as in 'geo_center', the only
    !          difference being that, the atomic masses are also
    !          supplied
    ! -------------------------------------------------------
    integer(kind=i4_kind),intent(in)        :: n
    real(kind=r8_kind),intent(in)           :: x(:),y(:),z(:) !(n)
    real(kind=r8_kind),intent(in)           :: mass(:)        !(n)
    real(kind=r8_kind),intent(out)          :: xmass(3)
    real(kind=r8_kind),intent(out),optional :: rotaxes(3,3)
    ! --- Declaration of local variables --------------------
    integer(kind=i4_kind)                     :: comp,ierr
    real(kind=r8_kind)                        :: fv1(3),fv2(3)
    real(kind=r8_kind)                        :: mtensor(3,3),eigrot(3)

    ASSERT(n==size(x))
    ASSERT(n==size(y))
    ASSERT(n==size(z))
    ASSERT(n==size(mass))

    xmass = centerofmass(x,y,z,mass)

    if(present(rotaxes)) then
      !*** Calculate tensor of moment of inertia with the function 'minerttensor' ********
      mtensor = minerttensor(x,y,z,mass,xmass)
      call rs(3,3,mtensor,eigrot,1,rotaxes,fv1,fv2,ierr)
      if(ierr/=0)  stop 'geo_center: diagonalization failed'

      if(mass_center_print) then
        write(OPT_STDOUT,*) 'Rotational moments and corresponding eigenvectors:'
        write(OPT_STDOUT,*)" --------------------------------------------------------------"
        write(OPT_STDOUT,*)"     principial        eigvec 1       eigvec 2       eigvec 3"
        write(OPT_STDOUT,*)"      moments"
        write(OPT_STDOUT,*)" --------------------------------------------------------------"
        do comp=1,3
        write(OPT_STDOUT,'("     ",4F15.7)') eigrot(comp),rotaxes(comp,:)
        end do
        write(OPT_STDOUT,*)" --------------------------------------------------------------"
      endif
    endif
  end subroutine mass_center
  !***************************************************************************************

  function minerttensor(x,y,z,mass,xmass) result(tensor)
    ! Purpose: calculation of the moment of inertia tensor
    ! --------------------------------------------------------
    real(kind=r8_kind),intent(in)           :: x(:),y(:),z(:)
    real(kind=r8_kind),intent(in)           :: mass(:)
    real(kind=r8_kind),intent(in)           :: xmass(3)
    real(kind=r8_kind),dimension(3,3)       :: tensor
    !*** End of interface ******************************************

    integer(kind=i4_kind)                   :: i_atom = 0

    !*** executable code *****************************************************************
    tensor=0.0_r8_kind
    do i_atom=1,size(x)
      !Please see vib.pdf from gaussian.com
      tensor(1,1)=tensor(1,1)+((z(i_atom)-xmass(3))**2&
        +(y(i_atom)-xmass(2))**2)*mass(i_atom)
      tensor(2,2)=tensor(2,2)+((x(i_atom)-xmass(1))**2&
        +(z(i_atom)-xmass(3))**2)*mass(i_atom)
      tensor(3,3)=tensor(3,3)+((x(i_atom)-xmass(1))**2&
        +(y(i_atom)-xmass(2))**2)*mass(i_atom)
      tensor(1,2)=tensor(1,2)-(x(i_atom)-xmass(1))*&
        (y(i_atom)-xmass(2))*mass(i_atom)
      tensor(1,3)=tensor(1,3)-(x(i_atom)-xmass(1))*&
        (z(i_atom)-xmass(3))*mass(i_atom)
      tensor(2,3)=tensor(2,3)-(y(i_atom)-xmass(2))*&
        (z(i_atom)-xmass(3))*mass(i_atom)
    end do
    tensor(2,1)=tensor(1,2)
    tensor(3,1)=tensor(1,3)
    tensor(3,2)=tensor(2,3)

  end function minerttensor
  !***************************************************************************************

  function centerofmass(x,y,z,mass) result(masscenter)
    ! Purpose: calculation of the molecules center of mass
    ! --------------------------------------------------------
    real(kind=r8_kind),intent(in)           :: x(:),y(:),z(:)
    real(kind=r8_kind),intent(in)           :: mass(:)
    real(kind=r8_kind),dimension(3)         :: masscenter
    !*** End of interface ******************************************

    integer(i4_kind) :: i_atom

    masscenter = ZERO
    DO i_atom=1,size(x)
      masscenter = masscenter + (/ x(i_atom), y(i_atom), z(i_atom) /) * mass(i_atom)
    END DO
    masscenter = masscenter / sum(mass)
  end function centerofmass
  !***************************************************************************************

  subroutine xmol_output(filename)
    ! Purpose: create file that can be read in by the 'xmol'
    !          program.
    ! --------------------------------------------------------
    ! use atom_data_module
    use constants, only: angstrom
    character(len=*), intent(in) :: filename
    integer(kind=i4_kind) :: i,io_xmol=11
    character(len=3)      :: sym
    ! --- executable code ------------------------------------
    io_xmol=openget_iounit(status='unknown',form='formatted',&
         file=trim(filename)//".xyz")
    write(io_xmol,*)n_atoms+n_dummy
    write(io_xmol,*)"interesting molecule"
    do i=1,n_atoms+n_dummy
       sym = symbol(nint(charge(i)))
       write (io_xmol, 1000) &
            sym, x(i) / angstrom, y(i) / angstrom, z(i) / angstrom
    enddo
1000 format(a3,3(1x,f10.6))
    call returnclose_iounit(io_xmol)
  end subroutine xmol_output

!--------------- End of module ----------------------------------
end module geo_operations_module
