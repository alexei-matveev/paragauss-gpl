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
subroutine chargefit(loop, coeff_dev, coulomb_dev)
  !
  !  Purpose: This subroutine performs the charge fit.
  !
  !  This subroutine runs in mater-only context so far.
  !
  !  Variables used:
  !
  !  Name:              Description/Range:
  !  ssym               symmetry flags  -> datatypes
  !  eigvec             Coefficients for the density matrix
  !  densmat            Density matrix, type(arrmat3)
  !  eig                Eigenvalues
  !
  !
  !  Varaiables used and modified:
  !
  !  mat_charge         Two center-integrals for density fit
  !  coeff_charge:      fit coeffs for total charge distribution
  !  coeff_spin  :      fit coeffs for spin  charge distribution
  !
  !
  !  Subroutine called by: mainscf
  !
  !  Author: T. Grauschopf
  !  Date: 8/95
  !
  !--------------------------------------------------------------
  ! Modifications:
  ! 1. Interfaces for subroutines inserted -> otherwise dimension
  !    of assumed-shape arrays are wrong in the subroutines
  !
  ! 2. Information that is sent to the slaves via comm is now
  !    packed in ONE buffer and sent with only ONE message tag,
  !    namely the one which  results in the call to the appropriate
  !    subroutine in 'decide'. No other message is sent with an
  !    extra message tag to the slaves
  !
  ! FN 29.8.95
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   2/7/97
  ! Description: 1. calculation of density matrix moved to new subroutine
  !                 CALC_DENSMAT being a member of the DENSITY_DATA_MODULE
  !              2. introduction of the OPTIONS_PERTURBATION_THEORY switch
  !              3. restructured to separate the perturbation theory
  !                 from the (constrained) charge fit.
  !              4. perturbation theory now based on occupied_levels_module
  !                 and virtual_levels_module to avoid multiple sending of
  !                 eigenstates
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   8/97
  ! Description: Calculation of bounds moved to prescf
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: Calculation of EXT-MDA projection <rho_spin|g_l> introduced
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   6/98
  ! Description: Extension to Spin Orbit
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   3/99
  ! Description: I deleted commpack(n_pairs,info) - n_pairs is the global PARAMETER
  !              which is known for any part of programm it is used
  !
  ! Modification (Please copy before editing)
  ! Author:      ...
  ! Date:        ...
  ! Description: ...
  !
  !== Interrupt of public interface of module ===================
  !--------------------------------------------------------------
  ! MODULES
  !--------------------------------------------------------------

# include "def.h"
  use type_module    ! contains standard data types
  use comm, only: comm_rank, comm_bcast
  use comm_module, only: comm_init_send, comm_send, comm_all_other_hosts
  use msgtag_module, only: msgtag_charge_fit
  use time_module, only: start_timer, stop_timer
  use timer_module, only: timer_scf_chfit
  use symmetry_data_module, only: symmetry_data_n_irreps, &
    symmetry_data_n_spin, symmetry_data_dimension, &
    symmetry_data_n_proj_irreps, symmetry_data_dimension_proj
  use iounitadmin_module, only: output_unit
  use occupation_module, ONLY : get_n_elec, get_spin_diff, get_n_core_elec, &
                                n_rot, eigvec_kept, eigen_kept, keep_eigen, &
                                occupation_spindiff
  use eigen_data_module, only: eigval, eigvec, eigvec_real, eigvec_imag
  use fit_coeff_module,       only : coeff_charge                              &
                                   , coeff_charge_old                          &
                                   , coeff_charge_veff                         &
                                   , coeff_deltarho                            &
                                   , coeff_spin                                &
                                   , coeff_core                                &
                                   , charge_norm                               &
                                   , ch_initialized                            &
                                   , dr_initialized                            &
                                   , copy_coeff_charge_old                     &
                                   , free_coeff_charge_old                     &
                                   , fit_coeff_charge_norm                     &
                                   , fit_coeff_n_ch                            &
                                   , fit_coeff_n_xc                            &
                                   , fit_coeff_n_cd                            &
                                   , fit_coeff_set_ch                          &
                                   , fitfct_map
  ! FIXME: how comes ptpair2 is unused?
  use pairs_module, only: n_pairs0, n_pairs, ptpair1, deallocate_pairs
  use pert_coeff_module, only: pcalc_pert_coeff, update_charge_overlap_matrix
  use pert_coeff_module, only: pert_vec, pert_vec_core, pert_vec_xc, pert_vec2, pert_mat
  use pert_coeff_module, only: pert_mat_real, pert_mat_imag
  use pert_coeff_module, only: pert_coeff_free
  use linsys_module, only: decompose_linsys, solve_linsys
  use mat_charge_module, only: mat_charge, matinv_charge, dual_charge_norm
  use mat_charge_module, only: mat_xc_ch
  ! the following use-statement is for debug purposes only
  use output_module, only: output_chargefit, output_n_coeff_dev
  use mixing_module, only: mixing_ch, mixing_ch_ignore_coeff_old
  use convergence_module, only: convergence_check_coeff_dev, &
                                convergence_check_coulomb_dev, &
                                convergence_load_coeff_dev, &
                                convergence_load_metric_dev
  use xcmda_hamiltonian, only: mda_constrain_rhospin
  use options_module, only: options_perturbation_theory, options_xcmode, &
                            xcmode_model_density, xcmode_extended_mda, &
                            options_spin_orbit
  use operations_module, only: operations_core_density
  use occupied_levels_module, only: update_eigvec_occ
  use virtual_levels_module, only: eigvec_vir_dealloc
  use integralpar_module, only: integralpar_cpksdervs
  use spin_orbit_module, only: whatis, op_BackTrafo
  use interfaces, only: IMAST
  implicit none
  !------------ Declaration of formal parameters ----------------
  integer(i4_kind), intent(in) :: loop
  real(r8_kind), intent(out) :: coeff_dev, coulomb_dev
  ! *** end of interface ***


  !------------ Declaration of local variables --------------------
  integer(kind=i4_kind) :: n, ni,&
       i_s, i, j, k, l, i_r, ii, alloc_stat, &
       n_spin, n_irrep, n_ch, n_xc, n_mat, n_proj, n_cd
  real(r8_kind) :: coeff_norm
  real(kind=r8_kind) :: c1,c2,f,fac,a,b,help_real,beta,deig,target_charge,&
       help_real_real,help_real_imag
  real(kind=r8_kind) :: c2_real, c2_imag, a_real, a_imag
  real(kind=r8_kind), dimension(:,:), allocatable :: al, bl
  real(kind=r8_kind), dimension(:), allocatable :: al_real, al_imag, bl_real
  real(kind=r8_kind), dimension(:), allocatable :: linsys_mat, linsys_dual
  ! the core fit functions are distributed precisely as the charge fit functions
  logical :: eigvec_rotated
  logical                             :: perturbation_theory,ignore_coeff_old,&
                                         extended_mda
  external error_handler
  ! ----------------------------------------------------------------
  ! first check if the norm is still ok
  ! This way the exact charge density of the cirrent cycle is compared
  ! with the fitted on of the  p r e v i o u s  cycle.
  ! call fit_coeff_charge_norm()

  DPRINT 'chargefit: entered'

  if ( comm_rank() == 0 ) then
      !
      ! Tell slaves to enter chargefit():
      !
      call comm_init_send(comm_all_other_hosts, msgtag_charge_fit)
      call comm_send()
  endif

  !
  ! From here on all workers run in a parallel context!
  !

  ! This timer  was previously in  main_scf() then executed  at master
  ! only.   Other  steps like  broadcasting  eigenvectors and  density
  ! matrix generation were also included back then:
  call start_timer (timer_scf_chfit)

  n_ch = fit_coeff_n_ch()
  n_xc = fit_coeff_n_xc()
  n_mat = (n_ch*(n_ch+1))/2
  if (operations_core_density) then
     n_cd = fit_coeff_n_cd()
  endif
  n_spin = symmetry_data_n_spin()

  if (options_spin_orbit) then
     !
     ! SPIN ORBIT
     !
     n_irrep = symmetry_data_n_proj_irreps()
  else
     n_irrep = symmetry_data_n_irreps()
  endif

  ! ATTENTION perturbation_theory flag may be changed during the SCF run
  ! only the master has the most recent value
  perturbation_theory = options_perturbation_theory()
  !
  ! FIXME: please make sure that options_perturbation_theory() is consistent
  !        across processes.
  !
  call comm_bcast(perturbation_theory)

  extended_mda = options_xcmode() == xcmode_extended_mda
  if (options_xcmode() == xcmode_model_density .or. extended_mda) then
     n_proj = n_spin
  else
     n_proj = 1
  endif

  !
  ! Calculate the charge density projection [rho|f_k], the
  ! extended MDA exchange fit function variant <rho|g_l> and
  ! all related quantities required for the perturbation theory
  !
  ! One of the first things the master does when entering the
  ! subroutine pcalc_pert_coeff() is sending a message to slaves
  ! telling them to also call pcalc_pert_coeff(). From that
  ! moment on all workers cooperate on computing the RHS for
  ! the linear equation system below (and more).
  !
  ! See pert_coeff_module.f90 for details ...
  !
  call pcalc_pert_coeff()

  !
  ! Historically, chargefit.f90 was executed by master only,
  ! slaves cooperated while constructing the RHS, master did
  ! the rest:
  !
  if ( comm_rank() /= 0 ) GOTO 999 ! finalize and exit

  !
  ! From here on these variables imported from
  ! pert_coeff_module.f90 should contains meaningfull data
  ! (ok on master only, maybe not so meaningfull on slaves)
  !
  ! pert_vec(1:n_ch, 1:n_proj)
  ! pert_vec2(1:n_ch, *)
  ! pert_mat(1:n_ch, *)
  ! pert_mat_real(1:n_ch, *)
  ! pert_mat_imag(1:n_ch, *)
  !

  if (n_proj > 1) then ! spin polarized model density calculation
     pert_vec(:,1) = pert_vec(:,1) + pert_vec(:,2)                 ! rho_tot
     pert_vec(:,2) = pert_vec(:,1) - pert_vec(:,2) - pert_vec(:,2) ! rho_spin
  endif

  if(output_chargefit) then
     pert_output: if (perturbation_theory) then
     write(output_unit,*) 'Anzahl Paare ',n_pairs
     write(output_unit,*) 'Groesse Stoerungstheoriegleichungssystem',n_ch
     do l=1,n_irrep
        do i=1,n_spin
           if (options_spin_orbit .and. i==1) then
              write(output_unit,*) 'pert_mat'
           else
              write(output_unit,*) 'pert_mat Spin ',i
           endif
           do j=1,n_pairs
              write(output_unit,*) 'Koeffizienten zu Paarnr. ',j
              if (options_spin_orbit .and. i==1) then
                 write(output_unit,'(  A   )') 'Real Part:'
                 write(output_unit,'(10F8.4)') pert_mat_real(:,l,j)
                 write(output_unit,'(  A   )') 'Imaginary Part:'
                 write(output_unit,'(10F8.4)') pert_mat_imag(:,l,j)
              else
                 write(output_unit,'(10F8.4)') pert_mat(:,l,i,j)
              endif
              write(output_unit,*)
           enddo
           write(output_unit,*)
        enddo
     enddo
     endif pert_output
     write(output_unit,*) 'PERT_VEC(tot)'
     write(output_unit,'(10F8.4)') pert_vec(:,1)
     write(output_unit,*)
     if (n_proj > 1) then ! spin polarized model density calculation
        write(output_unit,*) 'PERT_VEC(spin)'
        write(output_unit,'(10F8.4)') pert_vec(:,2)
        write(output_unit,*)
     endif
     write(output_unit,*) 'PERT_VEC_core(tot)'
     write(output_unit,'(10F8.4)') pert_vec_core(:)
     write(output_unit,*)
     if (extended_mda) then
        write(output_unit,*) 'PERT_VEC_XC(tot)'
        write(output_unit,'(10F8.4)') pert_vec_xc(:,1)
        write(output_unit,*)
        if (n_proj > 1) then ! spin polarized model density calculation
           write(output_unit,*) 'PERT_VEC_XC(spin)'
           write(output_unit,'(10F8.4)') pert_vec_xc(:,2)
           write(output_unit,*)
        endif
     endif
  endif

  !
  ! Save old charge fit coefficients before updating them
  !
  call copy_coeff_charge_old()

  !
  ! mixing_ch_ignore_coeff_old = loop <= IDYSTR - 1
  !     .and. chold_initialized .and. discard_init
  !
  ignore_coeff_old = mixing_ch_ignore_coeff_old(loop)

  !
  ! Ultimate target for the density fit norm:
  !
  call get_n_elec(target_charge)

  !
  ! Current (soon to become old) density fit norm:
  !
  coeff_norm = dot_product(coeff_charge, charge_norm)

  DPRINT   "chargefit: current norm", coeff_norm
  DPRINT   "chargefit:  target norm", target_charge
  DPRINT   "chargefit:         diff", coeff_norm - target_charge

  !
  ! Current norm may be deficient (e.g. zero inital coeff_charge or by numerics)
  ! use this difference as the target charge:
  !
  if ( ignore_coeff_old ) then
    !
    ! Solve for a in F * a = b  with  Q * a = N
    !
    ! Keep target charge N and rhs b
    !                               k
  else
    !
    ! Solve for da in F * da = db  with  Q * da = dN instead
    !
    target_charge = target_charge - coeff_norm

    DPRINT   "chargefit:   new target", target_charge

    !
    ! Prepare right hand side: b_k = b_k - Sum(k) [f_k|f_l] a_l(old)
    !
    k = 0
    do i = 1, n_ch
       do j = 1, i - 1
          k = k + 1
          pert_vec(i, 1) = pert_vec(i, 1) - mat_charge(k) * coeff_charge(j)
          pert_vec(j, 1) = pert_vec(j, 1) - mat_charge(k) * coeff_charge(i)
       enddo
       k = k + 1
       pert_vec(i, 1) = pert_vec(i, 1) - mat_charge(k) * coeff_charge(i)
    enddo
  endif

  if( output_chargefit) then
     write(output_unit,*)
     write(output_unit,*) 'Werte von LINSYS_TARGET'
     write(output_unit,'(F12.8)') target_charge
     write(output_unit,*)
     write(output_unit,*) 'Werte von LINSYS_RHS (including -F*a_old)'
     write(output_unit,'(F12.8)') pert_vec(:,1)
  endif

  pert_solve: if ( perturbation_theory ) then
    !
    ! Set up and solve
    !
    !   (F + dF_pert) * da = b - F * a_old
    !
    ! with
    !
    !   Q * da = 0
    !

    !
    ! Setup F_ij = F_ij(orig) + dF_ij(pert)
    !
    allocate(linsys_mat((n_ch*(n_ch+1))/2),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    linsys_mat = mat_charge

    if (.not.ignore_coeff_old) then
       !  assemble correction matrix to density fit
       call update_charge_overlap_matrix(linsys_mat)
    endif

    if(output_chargefit) then
       write(output_unit,*) 'D1'
       k = 0
       do i=1,n_ch
          do j=1,i
             k = k+1
             write(output_unit,*) 'linsys_mat(',i,',',j,')=',linsys_mat(k)
          enddo
       enddo
    endif

    !
    ! Linear solver
    ! note that size of linsys_dual is n_ch + 1 !
    !
    allocate(linsys_dual(n_ch + 1), stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    call decompose_linsys(n_ch, linsys_mat, charge_norm, linsys_dual) ! (1)

    call solve_linsys(n_ch, linsys_mat, pert_vec(:, 1), linsys_dual, target_charge)

    deallocate(linsys_dual, linsys_mat, stat=alloc_stat)
    ASSERT(alloc_stat.eq.0)

  else pert_solve
    !
    ! Set up and solve: F * a = b  with  Q * a = N
    !
    ! NOTE: matinv_charge is the FACTORIZED matrix for the CONSTRAINED (dual)
    ! linear equation system; dual_charge_norm is the corresponding transformed
    ! version of charge_norm.
    !
    ! In this branch one optimizes the call to decompose_linsys() away
    ! because the matrix defining the linear equations is always the same.
    !
    call solve_linsys(n_ch, matinv_charge, pert_vec(:, 1), dual_charge_norm, target_charge)
  endif pert_solve

  if ( ignore_coeff_old ) then
    call fit_coeff_set_ch( pert_vec(:, 1) )

    !
    ! FIXME: pert_vec(:, 1) will be used below for perturbation
    ! theory of orbitals. The change from 0 to total is too big to be treated
    ! as perturbation. Thus make PT of orbitals a NOOP (?) by zeroing:
    !
    pert_vec(:, 1) = 0.0_r8_kind
  else
    !
    ! Regular branch:
    !
    call fit_coeff_set_ch( coeff_charge + pert_vec(:, 1) )
  endif

  DPRINT   "chargefit:     new norm", dot_product(coeff_charge, charge_norm)
  DPRINT   "chargefit: coeff_charge", coeff_charge
! DPRINT   "chargefit:   coeff_norm", charge_norm

  if(output_chargefit) then
     write(output_unit,*) 'D2'
     do i=1,n_ch
        write(output_unit,*) 'coeff_charge(',i,')=',coeff_charge(i)
     enddo
  endif

  projected: if (n_proj > 1) then ! spin polarized model density calculation
!    Set up and solve : F*a_spin = b_spin  with or without  Q*a_spin = N_spin
     if (mda_constrain_rhospin().and.loop>1) then
!!! MF wrong?! Only default magnetic momnet, not actual
!       call get_spin_diff(target_charge)
        target_charge=occupation_spindiff(signum=.true.)
        if( output_chargefit) then
           write(output_unit,*) 'Werte von LINSYS_TARGET(rho_spin)'
           write(output_unit,'(F12.8)') target_charge
        endif
        call solve_linsys(n_ch,matinv_charge,pert_vec(:,2),dual_charge_norm, &
                          target_charge)
     else
        call solve_linsys(n_ch,matinv_charge,pert_vec(:,2))
     endif
     coeff_spin = pert_vec(:,2)

     if(output_chargefit) then
        write(output_unit,*) 'D2'
        do i=1,n_ch
           write(output_unit,*) 'coeff_spin(',i,')=',coeff_spin(i)
        enddo
     endif

  endif projected

  ch_initialized = .false.

  coredens:if (operations_core_density) then
     ! Set up and solve : F*a_core = b_core  with  Q*a_core = N_core
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! pack mat_charge into linsys_mat by using mapping array fitfct_map
     allocate(linsys_mat((n_cd+1)*n_cd/2), stat = alloc_stat )
     ASSERT(alloc_stat==0)
     k = 0
     l = 0
     do i=1,n_ch
        do j=1,i
           k = k + 1
           if ( fitfct_map(i) .and. fitfct_map(j) ) then
              l = l + 1
              linsys_mat(l) = mat_charge(k)
           endif
        end do
     end do
     ! pack charge_norm into coeff_core by using mapping array fitfct_map
     ! (here, coeff_core is just a working array)
     l = 0
     do i = 1, n_ch
        if ( fitfct_map(i) ) then
           l = l + 1
           coeff_core(l) = charge_norm(i)
        end if
     end do
     ! pack pert_vec_core into itself using mapping array fitfct_map
     l = 0
     do i = 1, n_ch
        if ( fitfct_map(i) ) then
           l = l + 1
           pert_vec_core(l) = pert_vec_core(i)
        end if
     end do
     ! get the number of core electrons
     call get_n_core_elec(target_charge)
     if( output_chargefit) then
        write(output_unit,*) 'Werte von LINSYS_TARGET(core)'
        write(output_unit,'(F12.8)') target_charge
     endif

     ! now solve the linear equation system
     ! note that size of linsys_dual is n_cd + 1 !
     allocate(linsys_dual(n_cd+1), stat = alloc_stat )
     ASSERT(alloc_stat.eq.0)
     call decompose_linsys(n_cd,linsys_mat,coeff_core(1:n_cd),linsys_dual) !(2)
     call solve_linsys(n_cd,linsys_mat,pert_vec_core(1:n_cd),&
                       linsys_dual,target_charge)
     deallocate(linsys_dual,linsys_mat,STAT=alloc_stat)
     ASSERT(alloc_stat.eq.0)

     ! finally unpack pert_vec_core into coeff_core using mapping array fitfct_map
     coeff_core = 0.0_r8_kind
     l = 0
     do i = 1, n_ch
        if ( fitfct_map(i) ) then
           l = l + 1
           coeff_core(i) = pert_vec_core(l)
        end if
     end do
     if(output_chargefit) then
        write(output_unit,*) 'core density fit coefficients'
        do i=1,n_cd
           write(output_unit,*) 'coeff_core(',i,')=',pert_vec_core(i)
        enddo
     endif
     deallocate( pert_vec_core, stat=alloc_stat )
     ASSERT(alloc_stat==0)
  endif coredens

  ! now check if the norm is still ok
  call fit_coeff_charge_norm()
  if(output_chargefit) then
     write(output_unit,*) 'chargefit coefficients after normalization'
     do i=1,n_ch
        write(output_unit,*) 'coeff_charge(',i,')=',coeff_charge(i)
     enddo
     if (n_proj > 1) then
        write(output_unit,*) 'spinfit coefficients after normalization'
        do i=1,n_ch
           write(output_unit,*) 'coeff_spin(',i,')=',coeff_spin(i)
        enddo
     endif
  endif

  if(whatis(op_BackTrafo).eq.2)then
     ASSERT(.not.perturbation_theory)
     DPRINT  '2) coeff_charge=',sum(coeff_charge),maxval(abs(coeff_charge))

     WARN('OFFSET coeff_charge by coeff_charge_veff')
     call fit_coeff_set_ch( coeff_charge - coeff_charge_veff )
  endif

  !
  ! Mixing of charge coefficients and getting beta, the mixing ratio beta is
  ! used in perturbation theory code below:
  !
  call mixing_ch(loop, beta, mat_charge)

  if(output_chargefit) then
     write(output_unit,*) 'chargefit coefficients after mixing'
     do i=1,n_ch
        write(output_unit,*) 'coeff_charge(',i,')=',coeff_charge(i)
     enddo
     if (n_proj > 1) then
        write(output_unit,*) 'spinfit coefficients after mixing'
        do i=1,n_ch
           write(output_unit,*) 'coeff_spin(',i,')=',coeff_spin(i)
        enddo
     endif
  endif
  ! find the five largest deviations of the charge fit coefficients
  if (convergence_check_coeff_dev()) then
     coeff_dev = convergence_load_coeff_dev(coeff_charge_old,coeff_charge, &
                 output_n_coeff_dev,"MAX. DIFF. IN CHARGEFIT COEFFICIENTS")
  else
     coeff_dev = huge(0.0_r8_kind)
  endif
  ! find the self-interaction of the fitted charge deviation
  if( convergence_check_coulomb_dev()) then
     coulomb_dev = convergence_load_metric_dev(coeff_charge_old,coeff_charge, &
                   mat_charge,"COULOMB SELF-INTERACTION NORM OF DELTA RHO_FIT")
  else
     coulomb_dev = huge(0.0_r8_kind)
  endif
  ! the old charge coefficients are not needed any more, thus:
  call free_coeff_charge_old

  pert_orbitals: if (perturbation_theory) then
     if(.not.integralpar_cpksdervs)  call eigvec_vir_dealloc(IMAST)

     if (.not.eigen_kept) then
        ! allocate n_rot as temporary working array
        allocate(n_rot(n_irrep,n_spin),STAT=alloc_stat)
        ASSERT(alloc_stat==0)
     endif

     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        ! allocate space for the matrixelements of the perturbation operator
        allocate(al_real(n_irrep), al_imag(n_irrep), bl_real(n_irrep), STAT=alloc_stat)
        ASSERT(alloc_stat==0)

        ! prepare_final_rotation
        ! here pert_vec holds the actual change of the charge fit coefficients
        do i_r=1,n_irrep
              bl_real(i_r) = sum(pert_vec2(:, i_r, 1) * pert_vec(:, 1))
              al_real(i_r)=sum(pert_mat_real(:,i_r,n_pairs)*pert_vec(:,1))
              al_imag(i_r)=sum(pert_mat_imag(:,i_r,n_pairs)*pert_vec(:,1))
        enddo

        ! solve the saekular equation for the highest occupied and lowest unoccupied
        ! orbital for every irrep and spin
        bl_real = bl_real / 2.0_r8_kind
        n_rot=-1 ! must be < 0 to ensure ni = (n=n_rot) + 1 to be less then 1
        do i_r=1,n_irrep
                 !
                 ! If list of pairs is empty, cycle to the next irrep:
                 !
                 if ( n_pairs0(1, i_r) > n_pairs ) cycle

                 n = ptpair1(i_r)%m(n_pairs, 1)
                 ASSERT(n>0)

                 ! FIXME:  is there  a guarantee  this index  is still  in the
                 ! range? Isnt that the purpose of ptpair2?
                 ni = n + 1

                 ASSERT(n<=size(eigval(i_r)%m,1))
                 ASSERT(ni<=size(eigval(i_r)%m,1))

                 deig=max((eigval(i_r)%m(ni,1)-eigval(i_r)%m(n,1))/2.0_r8_kind,&
                      0.01_r8_kind)
                 c1 = deig - beta * bl_real(i_r)
                 c2_real=beta*al_real(i_r)
                 c2_imag=beta*al_imag(i_r)
                 f = c1 - sqrt(c1 * c1 + c2_real * c2_real + c2_imag * c2_imag)
                 fac = sqrt(f * f + c2_real * c2_real + c2_imag * c2_imag)
                 if(fac<=1.0E-5_r8_kind) cycle
                 n_rot(i_r,1) = n
                 a_real=c2_real/fac
                 a_imag=c2_imag/fac
                 b = f / fac
                 ! build up the new orbitals
                 do ii=1,symmetry_data_dimension_proj(i_r)
                    help_real_real=eigvec_real(i_r)%m(ii,n)
                    help_real_imag=eigvec_imag(i_r)%m(ii,n)
                    eigvec_real(i_r)%m(ii,n)=a_real*eigvec_real(i_r)%m(ii,n)-a_imag*eigvec_imag(i_r)%m(ii,n)+&
                         b * eigvec_real(i_r)%m(ii, ni)
                    eigvec_imag(i_r)%m(ii,n)=a_real*eigvec_imag(i_r)%m(ii,n)+a_imag*eigvec_real(i_r)%m(ii,n)+&
                         b * eigvec_imag(i_r)%m(ii, ni)
                    eigvec_real(i_r)%m(ii, ni) = -b * help_real_real + &
                         a_real*eigvec_real(i_r)%m(ii,ni)-a_imag*eigvec_imag(i_r)%m(ii,ni)
                    eigvec_imag(i_r)%m(ii, ni) = -b * help_real_imag + &
                         a_real*eigvec_imag(i_r)%m(ii,ni)+a_imag*eigvec_real(i_r)%m(ii,ni)
                 end do
        enddo

        deallocate(al_real, al_imag, bl_real, STAT=alloc_stat)
        ASSERT(alloc_stat==0)
        !? print*," did pertubation theory of orbitals"

     else ! options_spin_orbit
        !
        ! STANDARD SCF (NO SPIN ORBIT)
        !
        ! allocate space for the matrixelements of the perturbation operator
        allocate(al(n_irrep,n_spin),bl(n_irrep,n_spin),STAT=alloc_stat)
        ASSERT(alloc_stat==0)

        ! prepare_final_rotation
        ! here pert_vec holds the actual change of the charge fit coefficients
        do i_r=1,n_irrep
           do i_s=1,n_spin
              bl(i_r,i_s)=sum(pert_vec2(:,i_r,i_s)*pert_vec(:,1))
              al(i_r,i_s)=sum(pert_mat(:,i_r,i_s,n_pairs)*pert_vec(:,1))
           enddo
        enddo

        ! solve the saekular equation for the highest occupied and lowest unoccupied
        ! orbital for every irrep and spin
        bl=bl/2.0_r8_kind
        n_rot=-1 ! must be < 0 to ensure ni = (n=n_rot) + 1 to be less then 1
        do i_r=1,n_irrep
           do i_s=1,n_spin
                 !
                 ! If list of pairs is  empty, cycle to the next (irrep, spin)
                 ! block:
                 !
                 if ( n_pairs0(i_s, i_r) > n_pairs ) cycle

                 n = ptpair1(i_r)%m(n_pairs, i_s)
                 ASSERT(n>0)

                 ! FIXME:  is there  a guarantee  this index  is still  in the
                 ! range? Isnt that the purpose of ptpair2?
                 ni = n + 1

                 ASSERT(n<=size(eigval(i_r)%m,1))
                 ASSERT(ni<=size(eigval(i_r)%m,1))

                 deig=max((eigval(i_r)%m(ni,i_s)-eigval(i_r)%m(n,i_s))/2.0_r8_kind,&
                      0.01_r8_kind)
                 c1=deig-beta*bl(i_r,i_s)
                 c2=beta*al(i_r,i_s)
                 f=c1-sqrt(c1*c1+c2*c2)
                 fac=sqrt(f*f+c2*c2)
                 if(fac<=1.0E-5_r8_kind) cycle
                 n_rot(i_r,i_s) = n
                 ! save the original orbitals (if required)
                 if (keep_eigen) then
                    eigvec_kept(i_r)%m(:,1,i_s) = eigvec(i_r)%m(:,n,i_s)
                    eigvec_kept(i_r)%m(:,2,i_s) = eigvec(i_r)%m(:,ni,i_s)
                 endif
                 a=c2/fac
                 b=f/fac
                 ! build up the new orbitals
                 do ii=1,symmetry_data_dimension(i_r)
                    help_real=eigvec(i_r)%m(ii,n,i_s)
                    eigvec(i_r)%m(ii,n,i_s)=a*eigvec(i_r)%m(ii,n,i_s)+&
                         b*eigvec(i_r)%m(ii,ni,i_s)
                    eigvec(i_r)%m(ii,ni,i_s)=-b*help_real+&
                         a*eigvec(i_r)%m(ii,ni,i_s)
                 end do
           enddo
        enddo
        keep_eigen = .false.

        deallocate(bl,al,STAT=alloc_stat)
        ASSERT(alloc_stat==0)
     endif ! options_spin_orbit
  endif pert_orbitals

  ! load < rho_s - rhofit_s | g_l > for s = up and down
  if (extended_mda) then
     if (n_spin == 1) then
        do i=1,n_ch
           pert_vec_xc(:,1) = pert_vec_xc(:,1) - &
                              mat_xc_ch(:,i)*coeff_charge(i)
        end do
        if (output_chargefit) then
           write(output_unit,*) 'PERT_VEC_XC(tot) including <rho_fit|g_l>'
           write(output_unit,'(10F8.4)') pert_vec_xc(:,1)
        endif
     else
        ! load charge fit coefficients for spin up and down
        pert_vec(:,1) = 0.5_r8_kind*( coeff_charge(:) + coeff_spin(:) )
!!!MF
!       pert_vec(:,2) = 0.5_r8_kind*( coeff_charge(:) - coeff_spin(:) )
        do i_s=1,n_proj
           do i=1,n_ch
              pert_vec_xc(:,i_s) = pert_vec_xc(:,i_s) - &
                                   mat_xc_ch(:,i)*pert_vec(:,i_s)
           end do
        end do
        if (output_chargefit) then
           write(output_unit,*) 'PERT_VEC_XC(up) including <rho_fit|g_l>'
           write(output_unit,'(10F8.4)') pert_vec_xc(:,1)
           write(output_unit,*) 'PERT_VEC_XC(down) including <rho_fit|g_l>'
           write(output_unit,'(10F8.4)') pert_vec_xc(:,2)
        endif
     end if

     coeff_deltarho = pert_vec_xc
     dr_initialized = .false.
  endif

999 CONTINUE
    !
    ! Here again, all workers run in a parallel context ...
    !

    if ( perturbation_theory ) then
        !
        ! Rotated eigenvectors need to be updated on all workers:
        !
        if ( comm_rank() == 0 ) then
            eigvec_rotated = any(n_rot > 0)
        endif

        ! This flag was set to a meaniful value only on master:
        call comm_bcast(eigvec_rotated)

        if ( eigvec_rotated ) then
            !
            ! Perturbation theory modified (?) eigenvectors, update
            ! them everywhere (another question is WHY we do this):
            !
            call update_eigvec_occ()
        endif

        !
        ! FIXME: we are in master only context here, doing manipulaitons
        !        in the global state of other modules:
        !
        if ( comm_rank() == 0 ) then
            if (.not.eigen_kept) then
                deallocate(n_rot, STAT=alloc_stat)
                ASSERT(alloc_stat==0)
            endif
        endif
    endif

    !
    ! Dealocate module variables in pert_coeff_module:
    !
    call pert_coeff_free()

    if ( perturbation_theory ) then
        !
        ! This cleans up pairs_module()
        !
        call deallocate_pairs()
    endif

  ! See the top for the corresponding start:
  call stop_timer (timer_scf_chfit)

  DPRINT 'chargefit: exit'
end subroutine chargefit
