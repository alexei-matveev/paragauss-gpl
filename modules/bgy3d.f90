!===============================================================
! Public interface of module
!===============================================================
module bgy3d
  !---------------------------------------------------------------
  !
  ! Copyright (c) 2014 Alexei Matveev
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
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
  use type_module, only: r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------


  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: rism_term ! (unique_atoms, energy, gradient_cartesian)

  !================================================================
  ! End of public interface of module
  !================================================================

  interface
     subroutine bgy3d_molmech (n, x, e, g) bind(c)
       use iso_c_binding, only: c_int, c_double
       implicit none
       integer (c_int), intent (in), value :: n
       real (c_double), intent (in) :: x(3, n)
       real (c_double), intent (out) :: e, g(3, n)
     end subroutine bgy3d_molmech
  end interface

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !------------ Subroutines ---------------------------------------
contains

  subroutine rism_term (unique_atoms, energy, gradient_cartesian)
    !
    ! Initiates computation  of "mechanic" gradients and  adds them to
    ! the accumulator.
    !
    use iso_c_binding, only: c_ptr
    use constants, only: angstrom, kcal
    use unique_atom_module, only: unique_atom_type
    use datatype, only: arrmat2
    implicit none
    type (unique_atom_type), intent(in) :: unique_atoms(:) ! (nua)
    real (r8_kind), intent (inout) :: energy
    type (arrmat2), intent (inout) :: gradient_cartesian(:) ! (nua) % m(1:3, nea)
    ! *** end of interface ***

    ! Nuclear coordinates in a (3 x natm) aray:
    real (r8_kind) :: x(3, sum (unique_atoms % n_equal_atoms))

    ! Gradients in a (3 x natm) aray:
    real (r8_kind) :: g(3, sum (unique_atoms % n_equal_atoms))
    real (r8_kind) :: e

    integer :: i, j, k

    ! First collect  atomic coordinates for nuclei in  a single array.
    ! Convert a.u.  to angstrom:
    k = 0
    do i = 1, size (unique_atoms)
       do j = 1, unique_atoms(i) % N_equal_atoms
          k = k + 1

          ! Coordinates in angstrom:
          x(:, k) = unique_atoms(i) % position(:, j) / angstrom
       enddo
    enddo
    ! At this point "k" is the total number of atoms:
    call bgy3d_molmech (k, x, e, g)
    print *, "rism_term: e=", e, "maxval(abs(g))=", maxval(abs(g))

    ! Increment  accumulators.  MM  energy  in  kcal/mol,  convert  to
    ! working units:
    energy = energy + e * kcal

    k = 0
    do i = 1, size (gradient_cartesian)
       do j = 1, size (gradient_cartesian(i) % m, 2)
          k = k + 1

          ! Gradients  in  units of  kcal/angstrom,  scale to  working
          ! units. In other if  the energy term is K * f  (x / A) then
          ! the derivatives wrt x is (K / A) * f'(x / A):
          associate (v => gradient_cartesian (i) % m (1:3, j))
            v = v + g(1:3, k) * (kcal / angstrom)
          end associate
       enddo
    enddo
  end subroutine rism_term

  !--------------- End of module ----------------------------------
end module bgy3d
