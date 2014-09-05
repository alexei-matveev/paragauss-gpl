!=====================================================================
! Public interface of module
!=====================================================================
module bgy3d
  !-------------------------------------------------------------------
  !
  ! Copyright (c) 2014 Bo Li
  ! Copyright (c) 2014 Alexei Matveev
  !
  ! This  program is  free software;  you can  redistribute  it and/or
  ! modify  it under  the  terms  of the  GNU  General Public  License
  ! version 2 as published by the Free Software Foundation [1].
  !
  ! This program  is distributed in the  hope that it  will be useful,
  ! but  WITHOUT ANY WARRANTY;  without even  the implied  warranty of
  ! MERCHANTABILITY or  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  ! General Public License for more details.
  !
  ! [1] http://www.gnu.org/licenses/gpl-2.0.html
  !
  ! Please see the accompanying LICENSE file for further information.
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
  use type_module, only: i4_kind, r8_kind ! type specification parameters
#ifdef WITH_GUILE
  use scm, only: scm_t
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------


  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: bgy3d_init_scheme!()
  public :: bgy3d_finalize!()
  public :: bgy3d_term !(iter, unique_atoms, densmat, ham_tot, e)
  public :: rism_term ! (unique_atoms, energy, gradient_cartesian)

  !===================================================================
  ! End of public interface of module
  !===================================================================

  interface
     subroutine bgy3d_pot_destroy (ctx) bind(c)
       use iso_c_binding, only: c_ptr
       implicit none
       type(c_ptr), intent(in), value :: ctx
     end subroutine bgy3d_pot_destroy

     function bgy3d_pot_get_value (ctx, n, x, v, p) result (offset) bind (c)
       use iso_c_binding, only: c_ptr, c_int, c_double, c_bool
       implicit none
       type (c_ptr), intent (in), value :: ctx
       integer (c_int), intent (in), value :: n
       real (c_double), intent (out) :: x(3, n)
       real (c_double), intent (out) :: v(n)
       integer (c_int), intent (out) :: p
       logical (c_bool) :: offset
     end function bgy3d_pot_get_value

     subroutine bgy3d_pot_interp (ctx, n, x, v) bind (c)
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type (c_ptr), intent (in), value :: ctx
       integer (c_int), intent (in), value :: n
       real (c_double), intent (in) :: x(3, n)
       real (c_double), intent (out) :: v(n)
     end subroutine bgy3d_pot_interp

     subroutine bgy3d_molmech (n, x, e, g) bind(c)
       use iso_c_binding, only: c_int, c_double
       implicit none
       integer (c_int), intent (in), value :: n
       real (c_double), intent (in) :: x(3, n)
       real (c_double), intent (out) :: e, g(3, n)
     end subroutine bgy3d_molmech
  end interface

  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

#ifdef WITH_GUILE
  type(scm_t), private :: solvent_hook
  type(scm_t), private :: solute_hook
  type(scm_t), private :: restart_destroy_hook
  type(scm_t), private :: restart
  logical, private :: once = .false.
#endif
  logical, private :: bgy3d_in_use = .false.

  logical, parameter, private :: debug = .false.
  !------------ Subroutines ------------------------------------------
contains

  function density (x, unit) result(rho)
    !
    ! Simple  interface  to  compute  the density  of  electrons  with
    ! arbitrary length  units. The "unit" specifies the  length of the
    ! unit vectors corresponding the  coordinates "x" in atomic units:
    ! when coordinates are measured in  angstroms unit == 1 angstrom ~
    ! 1.9 au.  When  coordinates are in au, unit  = 1, obviously. This
    ! identity holds:
    !
    !   density (x, A) == A**3 * density (A * x, 1):
    !
    use orbitalstore_module, only: orbital_type
    use orbital_module, only: orbital_setup, orbital_allocate, &
        orbital_free, orbital_shutdown, orbital_calculate
    use density_calc_module, only: density_calc_setup, density_calc_close, &
        density_calc_nl
    use machineparameters_module, only: nmax => machineparameters_veclen
    use symmetry_data_module, only: ssym
    use options_module, only: options_spin_orbit
    implicit none
    real(r8_kind), intent(in) :: x(:, :) ! (3, n), coordinates
    real(r8_kind), intent(in) :: unit    ! unit length in au
    real(r8_kind) :: rho(size(x, 2))     ! (n), output densities
    ! *** end of interface ***

    real(r8_kind) :: grdpts(nmax, 3) ! grid points
    type(orbital_type), pointer :: orbs_ob(:) ! see orbital_allocate()
    integer(i4_kind) :: i, nstart, nportion, nend

    !
    ! We  only need the  total density.   Subroutine density_calc_nl()
    ! returns the total when ispin = 1 or two densities for alpha- and
    ! beta- electrons when ispin = 2:
    !
    real(r8_kind) :: s_rho(nmax, ssym%n_spin) ! (polarized) density
    integer(i4_kind) :: n

    n = size (rho)

    ASSERT(size(x,1)==3)

    ! This code does not work for SO yet:
    ASSERT(.not. options_spin_orbit)

    !
    ! Preparations  begin, first  orbital_module.  Nmax  is  the chunk
    ! size,  by using  a different  value here  one can  change memory
    ! requirements and performance, but not the result:
    !
    call orbital_setup (nmax)

    ! Allocating  orbital objects.  The size  is proportional  to nmax
    ! times the number of basis orbitals:
    call orbital_allocate (orbs_ob = orbs_ob)

    call density_calc_setup()

    !
    ! Done with preparations ...
    !

    ! Start point
    nstart = 1
    do while (nstart <= n)

       ! Get possible maximum portion, less or equal to nmax
       nportion = min (1 + n - nstart, nmax)
       nend = nstart + nportion - 1

       ! Input  coordinates to  grid change  with portion.  BGY3d code
       ! provides  coordinates in  Angstroms but  QM codes  use atomic
       ! units. Convert coordinates to au:
       do i = 1, nportion
          grdpts(i, :) = x(:, nstart + (i - 1)) * unit
       end do

       ! calculate primitive, contracted, symmetry-adapted orbitals:
       call orbital_calculate (grdpts(1:nportion, 1:3), nportion, orbs_ob)

       ! calculate density
       call density_calc_nl (nportion, s_rho, orbs_ob = orbs_ob)

       ! rho for output, only updated the portion used, SUM(ARRAY,
       ! DIM, MASK). Convert the density to unit^-3:
       rho(nstart:nend) = sum(s_rho(1:nportion, :), 2) * unit**3

       ! Move to next portion
       nstart = nend + 1
    end do
    ASSERT(nstart==n+1)

    ! Clean up and exit:
    call orbital_free(orbs_ob = orbs_ob)
    call orbital_shutdown()

    call density_calc_close()
  end function density

  subroutine qm_density (n, x, rho) bind(c)
    !
    ! Simple interface  for use in BGY3d  code.  FIXME: it  is an ugly
    ! inconsistency  with  the rest  of  PG  code,  but this  function
    ! assumes the coordinates in  angstroms A and delivers the density
    ! in A^-3 as expected by BGY3d code.
    !
    use iso_c_binding, only: c_int, c_double
    use constants, only: angstrom
    implicit none
    integer(c_int), intent(in), value :: n ! number of points
    real(c_double), intent(in) :: x(3, n)  ! their coordinates
    real(c_double), intent(out) :: rho(n)  ! output densities
    ! *** end of interface ***

    !
    ! density (x, A) == A**3 * density (A * x, 1):
    !
    rho = density (x, angstrom)
  end subroutine qm_density

  subroutine bgy3d_pot_test (ctx)
    !
    ! test for interface usage
    !
    use iso_c_binding, only: c_ptr
    use machineparameters_module, only: nmax => machineparameters_veclen
    use comm, only: comm_allreduce
    implicit none
    type(c_ptr), intent(in) :: ctx
    ! *** end of interface ***

    real(r8_kind) :: x(3, nmax)
    real(r8_kind) :: v(nmax)
    real(r8_kind) :: m0, m1(3)
    real(r8_kind), parameter :: vol = 8000.0d0 ! use as parameter
    integer(i4_kind) :: n, times

    print *, "moments from fortran: "
    do times = 1, 3
      m0 = 0.0d0
      m1 = 0.0d0
      !
      ! FIXME: the noop  with .true.  .and.  is there  to workaround a
      ! regression in GCC 4.7 which otherwise fails at -O0:
      !
      do while (.true. .and. bgy3d_pot_get_value (ctx, nmax, x, v, n))
        m0 = m0 + sum(v(:n))
        m1 = m1 + matmul(x(:, :n), v(:n))
      end do

      call comm_allreduce(m0)
      call comm_allreduce(m1)
      print '("Moments divided by cell volume V = ", F11.6)',  vol
      print '("<1 * v> =", F11.6)', m0 / vol
      print '("<x * v> =", F11.6)', m1(1) / vol
      print '("<y * v> =", F11.6)', m1(2) / vol
      print '("<z * v> =", F11.6)', m1(3) / vol
    end do
  end subroutine bgy3d_pot_test


  subroutine integrate (ctx, densmat, ham_tot, e_elt)
    !
    ! Integrate  matrix elements <i|V|j>  of solvent  field V  for all
    ! pairs of orbitals i and j over (cell) volume.
    !
    use orbitalstore_module, only: orbital_type
    use orbital_module, only: orbital_setup, orbital_allocate, &
        orbital_free, orbital_shutdown, orbital_calculate
    use iso_c_binding, only: c_ptr
    use machineparameters_module, only: nmax => machineparameters_veclen
    use constants, only: angstrom, kcal ! units used by BGY code
    use options_module, only: options_spin_orbit
    use symmetry_data_module, only: ssym
    use f77_blas, only: dgemm
    use datatype, only: arrmat2, arrmat3
    use eigen_data_module, only: eigen_data_solve ! debug only
    use comm, only: comm_rank, comm_allreduce
    USE_MPI, only: MPI_Wtime
    implicit none
    type(c_ptr), intent(in) :: ctx             ! Context*
    type(arrmat3), intent(in) :: densmat(:)    ! for trace2()
    type(arrmat3), intent(inout) :: ham_tot(:) ! incremented!
    real(r8_kind), intent(out) :: e_elt        ! electron energy
    ! *** end of interface ***

    real(r8_kind) :: x(3, nmax) ! coordinates from bgy

    real(r8_kind) :: grdpts(nmax, 3) ! grid points used for orbital
    real(r8_kind) :: v(nmax) ! value of solvent field
    real(r8_kind), allocatable :: orbs_help(:,:)
    real(r8_kind) :: alpha, beta
    integer(i4_kind) :: n, i, j, k
    integer(i4_kind) :: n_irrep ! irreps
    integer(i4_kind), allocatable :: dims(:), partners(:)
    integer(i4_kind) :: maxdim ! maxval(dims(:))
    integer(i4_kind) :: alloc_stat
    type(orbital_type), pointer :: orbs_ob(:) ! see orbital_allocate()
    type(arrmat2), allocatable :: ham(:)
    integer(i4_kind) :: n_spin

    integer :: s
    ! variables to store time estimates
    real(r8_kind) :: time_tot, time_orb, time_intgr, time_before

    real(r8_kind), parameter :: one = 1.0

    ! Does not work for SO yet
    ASSERT (.not. options_spin_orbit)

    n_irrep = ssym % n_irrep
    n_spin = ssym % n_spin

    !
    ! Preparation work
    !

    time_tot = MPI_Wtime()      ! record start time

    ! preparations for intermediate arrays
    allocate (dims(n_irrep), partners(n_irrep), stat = alloc_stat)
    ASSERT (alloc_stat == 0)

    dims = ssym % dim           ! dimension of irreps
    partners = ssym % partner   ! partners

    ! dimension for intermediate arrays, max of irrep dimensions:
    maxdim = maxval (dims)

    allocate (orbs_help(nmax, maxdim), stat = alloc_stat)
    ASSERT (alloc_stat == 0)

    allocate (ham(n_irrep), stat = alloc_stat)
    ASSERT (alloc_stat == 0)

    do i = 1, size (ham)
       allocate (ham(i) % m(dims(i), dims(i)), stat=alloc_stat)
       ASSERT(alloc_stat==0)

       ham(i) % m = 0.0
    enddo

    ! preparation work for orbital
    call orbital_setup (nmax)
    call orbital_allocate (orbs_ob = orbs_ob)

    !
    ! End of preparation
    !

    !
    ! FIXME: the  noop with  .true.  .and.  is  there to  workaround a
    ! regression in GCC 4.7 which otherwise fails at -O0:
    !
    ! Fill x and v with coordinates and (weighted) solvent field:
    !
    do while (.true. .and. bgy3d_pot_get_value (ctx, nmax, x, v, n))

      ! convert coordinates and value of solvent field to au
      do i = 1, n
        grdpts(i, :) = x(:, i) * angstrom

        ! Potential is  returned in  kcal/mol premultiplied by  a grid
        ! weight  dV  in  A^3.   QM  code operates  in  atomic  units.
        ! Convert kcal/mol  * A^3 to  atomic units. XXX:  here solvent
        ! field acts on 'negative' unit electron
        v(i) = -1.0 * v(i) * kcal * angstrom**3
      end do

      time_before = MPI_Wtime()
      ! calculate primitive, contracted, symmetry-adapted orbitals:
      call orbital_calculate (grdpts(1:n, 1:3), n, orbs_ob)
      time_orb = time_orb + (MPI_Wtime() - time_before)

      !
      ! integral below
      !

      time_before = MPI_Wtime()
      do i = 1, n_irrep ! loop over irreps

        ! Make  sure ham(:)%m  is initialized  to zero.   DGEMM() will
        ! increment ham(i)%m on every call (beta == 1):
        alpha = one / partners(i)
        beta = one

        do j = 1, partners(i) ! loop partners

          ! work array: |V_solv|j>
          do k = 1, dims(i)
            orbs_help(1:n, k) = orbs_ob(i) % o(1:n, k, j) * v(1:n)
          end do

          ! call bals dgemm to calculate <i|V_solv|j>
          call dgemm ('t', 'n', dims(i), dims(i), n, alpha, &
                orbs_help(:, :), size (orbs_help, 1), &
                orbs_ob(i) % o(:, :, j), size(orbs_ob(i) % o, 1), beta, &
                ham(i) % m(:, :), size (ham(i) % m, 1))
        end do ! end loop over partners
      end do ! end loop over irreps
      time_intgr = time_intgr + (MPI_Wtime() - time_before)

      !
      ! end of interval
      !

    end do ! end loop to get values

    time_tot = MPI_Wtime() - time_tot ! duration

    ! These  will  be   cumulative  timings,  ideally  independent  of
    ! processor count:
    call comm_allreduce (time_tot)
    call comm_allreduce (time_orb)
    call comm_allreduce (time_intgr)

    if (comm_rank() == 0) then
      print *, 'bgy3d: TIMING total      =', time_tot
      print *, 'bgy3d: TIMING |-orbitals =', time_orb
      print *, 'bgy3d: TIMING |-integral =', time_intgr
      print *, 'bgy3d: TIMING (cumulative time)'
    endif

    do i = 1, size (ham)
       call comm_allreduce (ham(i) % m)
    enddo

    ! Electrostatic   interaction  of   the  electrons   with  solvent
    ! charge. Since  ham(:) has been  just all-reduced this  number is
    ! the same on all workers:
    e_elt = 0.0
    do i = 1, n_irrep
       do s = 1, n_spin
          e_elt = e_elt + trace2 (densmat(i) % m(:, :, s), ham(i) % m)
       enddo
    enddo

    ! Show   ham_tot,  eigval   and  eigvec   before   adding  solvent
    ! contribution
    if (debug .and. comm_rank() == 0) then
       print *, "Before adding solvent contribution"
       call show_ham()
    endif

    ! Add solvent  contribution to ham_tot(:). It is  probably not the
    ! most  efficient  way  to  reduce  evey  contribution  on  master
    ! separately. It  happens that  at diagonalization stage  only the
    ! value of the hamiltonian on master matters.
    if (comm_rank() == 0) then
       do i = 1, n_irrep
          do s = 1, n_spin
             ! ham_tot(i) % m(:, :, s) = ham_tot(i) % m(:, :, s) + ham(i) % m
          enddo
       enddo
    endif

    if (debug) then
       ! Diagonize 'new' ham_tot(:). This is a parallel context here:
       WARN("overwriting eigenvectors!")
       call eigen_data_solve()

       if (comm_rank() == 0) then
          print *, "After adding solvent contribution"
          call show_ham()
       endif
    endif

    ! Clean up and exit:
    call orbital_free(orbs_ob = orbs_ob)
    call orbital_shutdown()

    !
    ! Deallocation work
    !
    deallocate (dims, partners, orbs_help, stat = alloc_stat)
    ASSERT (alloc_stat == 0)

    do i = 1, size(ham)
       deallocate (ham(i) % m, stat = alloc_stat)
       ASSERT (alloc_stat == 0)
    enddo

    deallocate (ham, stat = alloc_stat)
    ASSERT (alloc_stat == 0)

    contains

      function trace2 (A, B) result (tr)
        !
        ! return the trace of multiplation of two symmetric matrices
        !
        implicit none
        real(r8_kind), intent(in) :: A(:, :), B(:, :)
        real(r8_kind) :: tr
        ! *** end of interface ***

        integer(i4_kind) :: i, j

        ! make sure input matrices are symmetric
        ASSERT (size (A, 1) == size (A, 2))
        ASSERT (size (B, 1) == size (B, 2))

        ! and two matrices are in the same dimension
        ASSERT (size (A, 1) == size (B, 1))

        tr = 0.0
        do i = 1, size (A, 1)
           do j = 1, size (A, 2)
              tr = tr + A(i, j) * B(j, i)
           enddo
        enddo
      end function trace2

      subroutine show_ham ()
        !
        ! Debug prints only.
        !
        use eigen_data_module, only: eigvec, eigval
        use overlap_module, only: overlap
        use debug, only: show
        implicit none
        ! *** end of interface ***

        integer :: i, s

        do i = 1, n_irrep
           call show ("ham", ham(i) % m)
           do s = 1, n_spin
              print *, "Matrix rep in basis of spin=", s, "eigenfunctions"
              call show ("eig", matmul (transpose (eigvec(i) % m(:, :, s)), &
                   matmul (ham(i) % m, &
                   eigvec(i) % m(:, :, s))))
              call show ("eig val", eigval(i) % m(:, s))
              call show("ham_tot", ham_tot(i) % m(:, :, s))

              ! debug only, show number of electrons
              print *, "number of electrons:", trace2(densmat(i) % m(:, :, s), &
                   overlap(i) % m(:, :))
           enddo
        enddo
      end subroutine show_ham

  end subroutine integrate


  function nuclear_energy (potential, unique_atoms) result (e_nuc)
    !
    ! Get the interaction between nuclei and solvent in Hartree.
    !
    use iso_c_binding, only: c_ptr
    use constants, only: angstrom, kcal
    use unique_atom_module, only: unique_atom_type
    implicit none
    type(c_ptr), intent(in) :: potential ! Context*
    type(unique_atom_type), intent(in) :: unique_atoms(:)
    real(r8_kind) :: e_nuc               ! result
    ! *** end of interface ***

    ! Nuclear coordinates in a (3 x natm) aray:
    real(r8_kind) :: x_nuc(3, sum (unique_atoms % n_equal_atoms))

    ! Proton number, or rather charge of the screened nuclei:
    real(r8_kind) :: chg_nuc(sum (unique_atoms % n_equal_atoms))

    ! Solvent field acts on nuclei site:
    real(r8_kind) :: v_nuc(sum (unique_atoms % n_equal_atoms))

    integer :: i, j, k

    ! First collect  atomic coordinates for nuclei in  a single array.
    ! Convert a.u.  to angstrom:
    k = 0
    do i = 1, size (unique_atoms)
       do j = 1, unique_atoms(i) % N_equal_atoms
          k = k + 1

          ! Coordinates in angstrom:
          x_nuc(:, k) = unique_atoms(i) % position(:, j) / angstrom

          ! Screened nuclear charge:
          chg_nuc(k) = unique_atoms(i) % Z - unique_atoms(i) % ZC
       enddo
    enddo

    ! At this point  "k" is the total number  of atoms.  Interpolating
    ! solvent field on nuclei sites:
    call bgy3d_pot_interp (potential, k, x_nuc, v_nuc)

    ! calculate energy between sovlent and nuclei
    e_nuc = 0.0
    do i = 1, k
      e_nuc = e_nuc + v_nuc(i) * chg_nuc(i)
    enddo

    ! Interpolated potential is in kcal/mol, convert to hartree:
    e_nuc = e_nuc * kcal
  end function nuclear_energy


  subroutine bgy3d_term (iter, unique_atoms, densmat, ham_tot, e_solv)
    !
    ! Increments the hamiltonian instance  ham_tot(:) on master by the
    ! solvation term.   The instances  of ham_tot(:) on  other workers
    ! are not affected. The energy contribution is, on the other hand,
    ! valid on all workers and does not need further reduction. FIXME:
    ! maybe one should avoid reduce operations for each and every term
    ! and do that once for the final hamiltonian/energy.
    !
    use iso_c_binding, only: c_funloc, c_funptr, c_ptr
    use scm, only: scm_t, scm_call, scm_car, scm_cdr, assignment(=)
    use scheme, only: scheme_atomic_sites
    use unique_atom_module, only: unique_atom_type
    use datatype, only: arrmat3
    use constants, only: kcal   ! for debug printing only
    use comm, only: comm_rank
    implicit none
    integer, intent(in) :: iter                ! starts at zero
    type(unique_atom_type), intent(in) :: unique_atoms(:)
    type(arrmat3), intent(in) :: densmat(:)    ! (n_irr)
    type(arrmat3), intent(inout) :: ham_tot(:) ! (n_irr), incremented
    real(r8_kind), intent(out) :: e_solv       ! solvation energy
    ! *** end of interface ***

    type(c_funptr) :: temp
    type(scm_t) :: solute, funptr, pair, unspecified
    type(c_ptr) :: field
    real(r8_kind) :: e_elt, e_nuc

    if (iter >= 0 .and. comm_rank() == 0) print *, "bgy3d: iter=", iter

    ! Module initialization:
    if (.not. once) then
       once = .true.

       ! This sets module global bgy3d_in_use:
       call bgy3d_init ()

    endif

    !
    ! FIXME: We need a better  way to tell PG compiled WITH_BGY3D that
    ! this is a  plain old "gas phase" calculation  that does not need
    ! BGY functionality.  The input file may need to be augmented. See
    ! how bgy3d_term() is unconditionally called from ham_calc_main().
    !
    ! This flag was set in bgy3d_init():
    if (.not. bgy3d_in_use) return

    if (iter < 0) then
       !
       ! Do nothing.  The caller didnt want to waste time on solvation
       ! term at the initial stage of SCF. Otherwise he/she would have
       ! supplied a non-negative number.
       !
    else if (iter == 0) then
       !
       ! Density  may   not  yet  be   available,  as  in   first  SCF
       ! iteration.  Run  pure  solvent  calculation.  Call  the  hook
       ! without  arguments  and  hope  the  solver  picks  the  right
       ! solvent.  Return value is unspecified as of yet:
       !
       unspecified = scm_call (solvent_hook)
    else
       ! GFortran   4.3  nees   a   temp  var   for   the  result   of
       ! c_funloc(). Next is a defined assignmend of a c_funptr to SCM
       ! exact:
       temp = c_funloc (qm_density)
       funptr = temp            ! SCM ptr <- c_funptr

       ! Solute  description. The  force field  parameters,  except of
       ! core charges, are not meaningful. Geometry of the solute does
       ! not change during SCF but may well become variable later:
       solute = scheme_atomic_sites (unique_atoms)

       ! Call  the hook  with  the  description of  QM  sites and  the
       ! address  of  qm_density()  and  get  back  a  solvent  medium
       ! description.
       pair = scm_call (solute_hook, solute, funptr, restart)

       ! At the moment the CAR position  of the result is a pointer to
       ! a C-struct describing  the solvent electrostatic field.  Here
       ! a defined assignment of an SCM value to a c_ptr is involved:
       field = scm_car (pair)   ! c_ptr <- SCM ptr

       ! The CDR position is an SCM  ptr to be passed to the next call
       ! of the  hook.  FIXME:  Here we save  it in the  module global
       ! variable  in order  to be  able to  pass it  back to  the BGY
       ! solver at  the next  call. If  you are not  going to  pass it
       ! back,   because,   say,    the   SCF   is   converged,   call
       ! bgy3d-restart-destroy on it.
       restart = scm_cdr (pair) ! SCM ptr

       ! Compute  matrix elements of  the solvent  field and  add them
       ! (eventually) to the Fock matrix. Also compute the interaction
       ! energy of the electrons with solvent charge:
       call integrate (field, densmat, ham_tot, e_elt)

       ! Interaction of nuclei with the solvent:
       e_nuc = nuclear_energy (field, unique_atoms)

       ! Electrostatic  interaction  energy  of  the solute  with  the
       ! solvent charge density, the total:
       e_solv = e_elt + e_nuc

       ! Only print the values on root
       if (comm_rank() == 0) then
          print *, "Interaction of solvent and nuclei:"
          print *, "e_nuc =", e_nuc, "a.u. =", e_nuc / kcal, "kcal/mol"
          print *, "Expectation value Tr(P_uv * ham):"
          print *, "e_elt =", e_elt, "a.u. =", e_elt / kcal, "kcal/mol"
          print *, "Total energy of solute with solvent:"
          print *, "e_slv =", e_solv, "a.u. =", e_solv / kcal, "kcal/mol"
       endif

       ! Here  is why we  treat field-  and restart  objects not  as a
       ! single  entity.  After  the  matrix elements  of the  solvent
       ! field have been computed, we destroy the field object to free
       ! memory for QM code.  The restart object needs to be preserved
       ! until the  next call to the  solute_hook and is  saved in the
       ! module global variable (see above):
       call bgy3d_pot_destroy (field)
    endif
  end subroutine bgy3d_term


  subroutine bgy3d_init ()
    !
    ! Set module global variables.
    !
    use iso_c_binding, only: c_null_ptr
    use scm, only: scm_car, assignment(=), scm_fluid_ref, scm_is_true
    use scheme, only: scheme_defined_p, scheme_lookup
    implicit none
    ! *** end of interface ***

    integer :: version


    !
    ! This  is a  short-curcuit  for PG  compiled  WITH_BGY3D but  run
    ! through  the  normal PG  startup  sequence  that  does not  know
    ! anything  about BGY.  The variable  bgy3d_in_use remains  set to
    ! false:
    !
    if (.not. scheme_defined_p ("bgy3d-api-version")) return

    !
    ! libbgy3d.a is  develeoped out of  PG tree.  Check if  the API
    ! version is compatible.  The Scheme variable bgy3d-api-version
    ! is a 3-list of (major  minor extra) integer numbers. Here int
    ! <- SCM int:
    !
    version = scm_car (scheme_lookup ("bgy3d-api-version"))
    ASSERT(version==2)

    !
    ! We do not (yet) interpret the BGY3d input here. We only check if
    ! it was  supplied. Set  the module global,  if false  this module
    ! should have no effect.
    !
    bgy3d_in_use = scm_is_true (scm_fluid_ref (scheme_lookup ("*input*")))

    restart = c_null_ptr     ! SCM ptr <- c_ptr
    solvent_hook = scheme_lookup ("bgy3d-solvent")
    solute_hook = scheme_lookup ("bgy3d-solute")
    restart_destroy_hook = scheme_lookup ("bgy3d-restart-destroy")
  end subroutine bgy3d_init


  subroutine bgy3d_finalize ()
    !
    ! Reset module global state.
    !
    use scm, only: scm_call
    implicit none
    ! *** end of interface ***

    ! See comments in bgy3d_term():
    if (.not. bgy3d_in_use) return

    ! We are  not going to restart  anymore. Maybe it  will later make
    ! sense to preserve the restart info across geometry iterations.
    restart = scm_call (restart_destroy_hook, restart)
  end subroutine bgy3d_finalize


#ifdef WITH_GUILE
  subroutine bgy3d_init_scheme() bind(c)
    !
    ! Exports a  few auxiliary methods  to the Scheme  interpreter. At
    ! the moment  bgy3d_init_scheme() is run from  qm_init().  So that
    ! the  gsubrs it  exposes  only become  available after  executing
    ! (qm-init)  in  Scheme interpreter  OR  upon (use-modules  (guile
    ! paragauss)).
    !
    use iso_c_binding, only: c_funptr, c_funloc
    use scm, only: scm_t, scm_define_gsubr
    implicit none
    ! *** end of interface ***

    type(c_funptr) :: fun
    type(scm_t) :: proc

    !
    ! FIXME: GFortran 4.3 fails if  you inline call to c_funloc(), use
    ! a temp var as a workaround:
    !
    fun = c_funloc (guile_qm_density)
    proc = scm_define_gsubr ("qm-density", 1, 0, 0, fun)
  end subroutine bgy3d_init_scheme


  function guile_qm_density (x) result(rho) bind(c)
    !
    ! Guile interface for use in scripts.
    !
    use scm, only: scm_t, scm_from
    use scheme, only: scheme_make_array
    use iso_c_binding, only: c_double
    implicit none
    type(scm_t), intent(in), value :: x ! list of three numbers
    type(scm_t) :: rho
    ! *** end of interface ***

    real(c_double) :: xs(3, 1), rhos(1)

    xs(:, 1) = scheme_make_array (3, x)

    call qm_density (1, xs, rhos)

    rho = scm_from (rhos(1))
  end function guile_qm_density
#endif

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

  !--------------- End of module -------------------------------------
end module bgy3d
