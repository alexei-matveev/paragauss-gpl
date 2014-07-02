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
module scheme
  !-------------------------------------------------------------------
  !
  !  Purpose:  Contains  bits  that  interoperate  with  the  embedded
  !  Scheme/Guile interpreter.
  !
  ! Copyright (c) Alexei Matveev
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
# include "def.h"
  ! Interoprability with C assumes RK  is the double precision, and IK
  ! is the default C-int:
  use type_module, only: IK => i4_kind, RK => r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  interface scheme_make_list
     module procedure scheme_make_list_int_1
     module procedure scheme_make_list_double_1
  end interface

  public :: scheme_make_list

  interface scheme_make_array
     module procedure scheme_make_array_double_1
     module procedure scheme_make_array_double_2
  end interface

  public :: scheme_make_array

  !------------ public functions and subroutines ---------------------

  public :: scheme_define       ! (key, val)
  public :: scheme_defined_p    ! (key) -> logical
  public :: scheme_lookup       ! (key) -> SCM val
  public :: scheme_find_basis   ! (atom, type)
  public :: scheme_trace_hook   ! (key, file, line, time)
  public :: scheme_atomic_sites ! (unique_atoms) -> SCM sites

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  function scheme_make_list_int_1 (array) result (list)
    use scm, only: scm_t, scm_list, scm_cons, scm_from
    implicit none
    integer(IK), intent(in) :: array(:)
    type(scm_t) :: list
    ! *** end of interface ***

    ! FIXME: Intel compiler has  troubles assembling the final list in
    ! the result directly. Introduce temp var:
    type(scm_t) :: temp
    integer :: i

    temp = scm_list ()
    do i = size (array), 1, -1
       temp = scm_cons (scm_from (array (i)), temp)
    enddo
    list = temp
  end function scheme_make_list_int_1

  function scheme_make_list_double_1 (array) result (list)
    use scm, only: scm_t, scm_list, scm_cons, scm_from
    implicit none
    real(RK), intent(in) :: array(:)
    type(scm_t) :: list
    ! *** end of interface ***

    ! FIXME: Intel compiler has  troubles assembling the final list in
    ! the result directly. Introduce temp var:
    type(scm_t) :: temp
    integer :: i

    temp = scm_list ()
    do i = size (array), 1, -1
       temp = scm_cons (scm_from (array (i)), temp)
    enddo
    list = temp
  end function scheme_make_list_double_1

  function scheme_make_array_double_1 (n, list) result (arr)
    use scm, only: scm_t, scm_car, scm_cdr, assignment(=)
    implicit none
    integer, intent(in) :: n
    type(scm_t), value :: list
    real(RK) :: arr(n)
    ! *** end of interface ***

    integer :: i

    do i = 1, n
       arr(i) = scm_car (list)  ! defined assignment

       ! passed by value, redefine:
       list = scm_cdr (list)
    enddo
  end function scheme_make_array_double_1

  function scheme_make_array_double_2 (n, m, list) result (arr)
    use scm, only: scm_t, scm_car, scm_cdr
    implicit none
    integer, intent(in) :: n, m
    type(scm_t), value :: list
    real(RK) :: arr(n, m)
    ! *** end of interface ***

    integer :: i

    do i = 1, m
       arr(:, i) = scheme_make_array (n, scm_car (list))

       ! passed by value, redefine:
       list = scm_cdr (list)
    enddo
  end function scheme_make_array_double_2


  subroutine scheme_define (key, val)
    !
    ! Purpose: (re)defines a Scheme variable. Here a double variable.
    !
    use scm, only: scm_t, scm_define, scm_from
    implicit none
    character(len=*), intent(in) :: key
    real(RK), intent(in) :: val
    ! *** end of interface ***

    type(scm_t) :: var

    var = scm_define (scheme_symbol (key), scm_from (val))
  end subroutine scheme_define


  function scheme_defined_p (key) result (yes)
    !
    ! String -> Logical
    !
    use scm, only: scm_is_true, scm_defined_p
    implicit none
    character(len=*), intent(in) :: key
    logical :: yes
    ! *** end of interface ***

    yes = scm_is_true (scm_defined_p (scheme_symbol (key)))
  end function scheme_defined_p


  subroutine scheme_find_basis(atom, type, uab)
    !
    ! Purpose: find a basis from a library.
    !
    use scm, only: scm_t, scm_call, scm_from, scm_length, scm_to_int, &
         scm_car, scm_cdr
    use unique_atom_module, only: unique_atom_basis_type, &
         unique_atom_make_empty_basis
    implicit none
    character(len=*), intent(in) :: atom, type
    type(unique_atom_basis_type), allocatable, intent(out) :: uab(:) ! array of shells
    ! *** end of interface ***

    type(scm_t) :: find_basis, basis
    integer :: lmax, L

    if (type == "fake") then
       !
       ! FIXME: PG is  broken with a zero-sized array  of shells. Make
       ! an empty s-shell:
       !
       allocate(uab(0:0))
       uab(0) = unique_atom_make_empty_basis()
    else
       !
       ! Lookup  Scheme  procedure qm-find-basis  that  does the  library
       ! search for us:
       !
       find_basis = scheme_lookup ("qm-find-basis")

       !
       ! Call the Scheme function passing the arguments (converted to
       ! Scheme values first):
       !
       basis = scm_call (find_basis, scm_from (atom), scm_from (type))

       !
       ! Number of shells, s-, p-, d-, etc. Zero or more:
       !
       lmax = scm_to_int (scm_length (basis)) - 1

       allocate (uab(0:lmax))
       do L = 0, lmax
          uab(L) = parse (L, scm_car (basis))
          basis = scm_cdr (basis)
       enddo
    endif
  contains

    function parse (L, shell) result (uab)
      use scm, only: scm_t, scm_call, scm_length, scm_car, &
           scm_cdr, assignment(=)
      use unique_atom_module, only: unique_atom_make_basis
      implicit none
      integer, intent(in) :: L
      type(scm_t), value :: shell
      type(unique_atom_basis_type) :: uab
      ! *** end of interface ***

      type(scm_t) :: symbol, exponents, contractions
      integer :: ne, nc

      ! symbol is unused, we assumed a proper spd sequence:
      symbol = scm_car (shell)

      ! passed by value, redefine:
      shell = scm_cdr (shell)

      exponents = scm_car (shell)
      contractions = scm_cdr (shell)

      ! number of exponents, contractions:
      ne = scm_length (exponents)    ! defined assignment
      nc = scm_length (contractions) ! defined assignment

      uab = unique_atom_make_basis (L, scheme_make_array (ne, exponents), &
           scheme_make_array (ne, nc, contractions))
    end function parse
  end subroutine scheme_find_basis

  function scheme_lookup (key, mod) result (val)
    !
    ! Mostly to remind one that Fortran  string is not the same as SCM
    ! string, a  symbol is not a  string, and a variable  is an object
    ! that needs to be "dereferenced" to get its value.
    !
    use scm, only: scm_t, scm_resolve_module, scm_lookup, &
         scm_variable_ref
    implicit none
    character (len=*), intent (in) :: key
    character (len=*), intent (in), optional :: mod
    type(scm_t) :: val
    ! *** end of interface **

    type(scm_t) :: var

    if (present (mod)) then
       var = scm_lookup (scm_resolve_module (mod), scheme_symbol (key))
    else
       var = scm_lookup (scheme_symbol (key))
    endif
    val = scm_variable_ref (var)
  end function scheme_lookup

  function scheme_symbol (str) result (sym)
    use scm, only: scm_t, scm_string_to_symbol, assignment(=)
    implicit none
    character(len=*), intent(in) :: str
    type(scm_t) :: sym
    ! *** end of interface **

    type(scm_t) :: string

    string = str                ! defined assignment
    sym = scm_string_to_symbol (string)
  end function scheme_symbol

  subroutine scheme_trace_hook (key, file, line, time)
    use scm, only: scm_t, scm_is_true, scm_defined_p, scm_call, &
         scm_variable_ref, scm_lookup, scm_from
    implicit none
    character(len=*), intent(in) :: key, file
    integer, intent(in) :: line
    double precision, intent(in) :: time
    ! *** end of interface ***

    type(scm_t) :: trace, ret

    trace = scheme_symbol ("qm-trace-hook")

    if (scm_is_true (scm_defined_p (trace))) then
       ret = scm_call ( &
            scm_variable_ref (scm_lookup (trace)), &
            scm_from (key), &
            scm_from (file), &
            scm_from (line), &
            scm_from (time))
    endif
  end subroutine scheme_trace_hook

  function scheme_atomic_sites (unique_atoms) result (sites)
    !
    ! Return molecular  description as a Scheme stucture  such as this
    ! one for Pd dimer:
    !
    ! (("Pd" (0.0 0.0 1.2) 0.0 0.0 46.0)
    !  ("Pd" (0.0 0.0 -1.2) 0.0 0.0 46.0))
    !
    ! The numbers are the cartesian  coordinates of a site in angstrom
    ! followed by a triple of (sigma, epsilon, charge) in angstrom, eV
    ! and "e" which are the force-field parameters.  The charge is the
    ! charge  of   the  nucleus,  eventually  screened   by  the  core
    ! electrons. This is  not the net charge according  to any kind of
    ! population  analysis. The Lennard-Jones  parameters are  fake so
    ! far.
    !
    use scm, only: scm_t, scm_from, scm_list, scm_cons
    use unique_atom_module, only: unique_atom_type
    use constants, only: angstrom
    implicit none
    type(unique_atom_type), intent(in) :: unique_atoms(:)
    type(scm_t) :: sites
    ! *** end of interface ***

    type(scm_t) :: site
    real(RK), parameter :: sigma = 0.0, epsilon = 0.0
    real(RK) :: position(3), charge
    integer :: ua, ea

    sites = scm_list ()               ! empty
    do ua = size(unique_atoms), 1, -1 ! cons them in reverse
       do ea = unique_atoms(ua)%n_equal_atoms, 1, -1
          position = unique_atoms(ua)%position(:, ea) / angstrom

          ! FIXME:  BGY3d code  uses two  Lennard-Jones  parameters in
          ! addition to  the site charge.  Where am I supposed  to get
          ! them?
          charge = unique_atoms(ua)%Z - unique_atoms(ua)%ZC

          ! We  deliver  coordinates in  atomic  units.  Make sure  to
          ! convert them to angstroms when calling BGY3d:
          site = scm_list (scm_from (trim (adjustl (unique_atoms(ua)%name))), &
               scheme_make_list (position), &
               scm_from (sigma), &
               scm_from (epsilon), &
               scm_from (charge))

          ! Cons the new site onto the head of the list:
          sites = scm_cons (site, sites)
       enddo
    enddo
  end function scheme_atomic_sites

  !--------------- End of module -------------------------------------
end module scheme
