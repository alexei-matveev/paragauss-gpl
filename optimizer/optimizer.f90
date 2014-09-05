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
module optimizer
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
# include "def.h"
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public main_opt

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine main_opt( task, converged,stop_after_eperelaxation,convert_internal, &
                       cross_boundary_3b)
    use type_module
    use constants, only: angstrom
    use math_module
    use coordinates_module
    use gradient_module
    use hesse_module
    use atom_data_module
    use opt_data_module
    use step_module, only: newton_step,step_max_comp,step_mean_square,number_nolines,noline
    use line_search_module, only: line_search_main,line_search_main2
    use valence_coord_module, only: valence_setup
    use geo_operations_module, only: xmol_output
    use ts_module, only: ts_main
    use frequency_module, only: frequency_main
    use slspar_module!, n_types_central_atoms_3b=>n_types_central_atoms_3body
    use iounitadmin_module
    use filename_module, only: inpfile
    use allocopt_module
#ifndef FPP_OPTIMIZER
    use options_module, only: update_hessian_iteration
    use operations_module, only: operations_qm_epe
#endif
#ifdef NEW_EPE
    use ewaldpc_module, only: opt_and_relax
    use qm_epe_interface_module, only: read_interface_param,free_ff_arrays
    use qm_epe_interface_module, only: at_type,interface_ff
#endif
#ifdef WITH_EFP
    use qmmm_interface_module, only: efp
    use efp_module, only: n_efp
#endif

    implicit none
    character(len=*), intent(in)            :: task
    logical         , intent(out)           :: converged
    logical         , intent(out), optional :: stop_after_eperelaxation
    logical         , intent(inout), optional :: convert_internal
    real(kind=r8_kind), intent(in), optional :: cross_boundary_3b
    ! *** end of interface ***


    real(r8_kind), parameter :: angsau = 1 / angstrom
  integer(kind=i4_kind)             :: io_gx,geo_loop,io_dip,io_gxstore,io_eq

  real(kind=r8_kind)                :: energy_prev=-1.0E15_r8_kind
  real(kind=r8_kind)                :: iwork, dipole_moment(3)
  real(kind=r8_kind) :: z_dummy
  real(kind=r8_kind)::small_dist=0.35_r8_kind,ant
  real(kind=r8_kind),dimension(3):: pos_temp

  logical                           :: do_newton,hesse_complete,do_step
!  logical :: epe_pc_shells
  logical :: gx_pos_regular=.false.  !regular pos of gx-file
  logical :: epe_pos_regular=.false.   !regular pos of epe lattice
  logical :: ewald_PCs=.false.         !ewald point charges
!!$  logical :: opened_stat
  ! to store OPT_STDOUT when it is redirected to io_flepo
  integer(kind=i4_kind) :: old_stdout = -1
  integer(kind=i4_kind):: io_epe_r=2
  integer(kind=i4_kind)  :: io_epe_pcs=14
  integer(kind=i4_kind)  :: io_epe_pcr=14
  integer(kind=i4_kind)  :: io_ewald_PCs=14
  integer(kind=i4_kind):: nucen_pcr,nucen_ewa
  integer(kind=i4_kind):: i,istat,n_counter,kl,kl1
  integer(kind=i4_kind):: k_counter,ieqt
  integer(kind=i4_kind), dimension(5)::ind_deck
  integer(kind=i4_kind):: ewa_low_limit
  real(kind=r8_kind) :: V_ewa,Coul_ewa
  logical :: have_hesse_cartesian


  type gx_reg
     real(kind=r8_kind),dimension(3)::r
     real(kind=r8_kind),dimension(3)::s
     real(kind=r8_kind),dimension(3)::eq
     real(kind=r8_kind)::q
     integer::  ieq
  end type gx_reg

  type(gx_reg), dimension(:),allocatable:: gx
  character(len=4):: char_geo_loop

  real(r8_kind),parameter:: gepe_const=evau*angsau**2

  ! List of I/O-Units
  ! io_hessecartes = 50
  ! io_val = 60
  ! io_conv = 35
  ! io_hesse = 15
  ! io_input = 30
  ! io_linmin = 20
  ! io_deloc = 25
  ! io_eigmod = 40
  ! io_epe_r =2
  ! io_epepar =14
  ! io_ffh_epepar = 12

 integer(kind=i4_kind):: internal_counter

 DPRINT  'main_opt(',trim(task),'): entered'

 ! what is this loop for?
 internal_loop: do internal_counter=1,999
 n_types_central_atoms_3body=0

  ! set the global from opt_data_module:
  ! will be re-set to io_flepo, after opening
  ! the latter:
#ifndef FPP_OPTIMIZER
  OPT_STDOUT = output_unit ! named constant (1) from iounitadmin_module
#else
  OPT_STDOUT = stdout_unit ! named constant (6) from iounitadmin_module
#endif

  write (OPT_STDOUT, *) "--- Optimizer started ---"

  !
  ! Optimizer seems to be using some SAVEed variables to store
  ! persistent data over geometry iterations. This state is lost
  ! after crash/restart cycle.
  !
  ! PLEASE DONT (AB)USE PERSISTENT STATE IF POSSIBLE!
  !
  ! Here the state of persistent optimizer variables should be restored
  ! if a valid (non empty) "optimizer.state" is found in the input directory:

  call persistent_state('restore')

  ! User have to create an empty (or valid) "optimizer.state" to initiate
  ! saving/restoring persistent variables on disk.
  ! By default, persistent vars are SAVEd in memory and get lost upon crash/termination.

  ! assume that not converged:
  converged = .false.
  if(present(stop_after_eperelaxation) ) stop_after_eperelaxation=.true.
  if(present(cross_boundary_3b)) crossboundary_3b=cross_boundary_3b


  ! read the input ----------------------------------------
  if(present(convert_internal)) then
    call opt_read_input(task,convert_internal=convert_internal)
    if(ts_scan.and.convert_internal) write(io_flepo,*) '  TS-SCAN CONVERT INTERNAL'
  else
    call opt_read_input(task)
  endif

#ifndef FPP_OPTIMIZER
  ! FIXME: what is the purpose of "analitic_hessian_calculated"?
  ! it seems to be always set to "false" independently of "task"
  ! in task_optimizer() ???
  ! it may be set to true with namelist for optimizer_only run
  ! when required
  if(analitic_hessian_calculated) update_hessian_iteration=1
#endif

  ! note different spelling of the two:
  if(present(stop_after_eperelaxation) )  &
              stop_after_eperelaxation=stop_after_epe_relaxation

#ifndef FPP_OPTIMIZER
  if(operations_qm_epe.and.optimization.and..not.cart_format) fixed_orientation=.true.
#endif
  if(fixed_orientation) then
     if(fixed_atom_1==0 .or. fixed_atom_2==0 .or. fixed_atom_3==0) &
          call error_handler("FIXED_ORIENTATION=TRUE, but not all fixed_atoms are defined")
  end if

#ifdef NEW_EPE
  if(operations_qm_epe) then
     epe_interfaced_mode=.true.
     epe_forces=.true.
     epeparameters=.false.
     if(opt_and_relax .and. present(stop_after_eperelaxation)) stop_after_eperelaxation=.false.
  end if
#endif

  if(epe_forces) then
     write(OPT_STDOUT,*)"start epe part"
     if(calc_epeff_hessian)  then
        call init_slsp_ffh
        call get_ffh_parameters
 !      print*,'get_ffh_parameters done'
        if(epeparameters.and.list_epepar) call list_ffh_par
     else
#ifndef NEW_EPE
        call init_slsp()       !
        !initialize known pair potentials
        if(epeparameters) call make_epe_namelist
        if(epeparameters.and.list_epepar) call list_epe_par
#endif
     endif
  endif

  DPRINT "main_opt: trying "//">"//adjustl(trim(opt_data_dir))//'/gxfile'//"<"

!  if(.not.linear_transit) then
   io_gx=openget_iounit(status='unknown',form='formatted', &
       file=trim(opt_data_dir)//'/gxfile')
   allocopt_stat(22)=0
   DPRINT "main_opt: opened on unit",io_gx
!  endif

  call coordinates_read(io_gx,iwork,xyz,step_counter) ! (1)
 !print*,'impu', impu(1:10)


  if(tsscan_mix.and..not.optimization) then
    io_gxstore=openget_iounit(file=trim(opt_data_dir)//'/gx.rp_mix', &
                              status='unknown',form='formatted')
   allocopt_stat(23)=0
    call coordinates_write(io_gxstore,-1.0_r8_kind,rp_mix=.true.)
    ! FIXME: shouldnt we also close io_gxstore before exit?
    goto 9999 ! finalize and then RETURN
  endif

  if(save_epe_r) then
     ! save regular positions of ions for use in epe interface mode
     ! to make proper relations between atoms

     io_epe_r=openget_iounit( file=adjustl(trim(inpfile('epe.r'))), &
          status='new',form='formatted')
   allocopt_stat(24)=0
     print*, 'file epe.r created'
     do i=1,n_atoms+n_dummy
        if(impu(i).ne.0) write(io_epe_r,*) x(i),y(i),z(i)
     enddo !i=1,kolat
     call returnclose_iounit(io_epe_r)
     allocopt_stat(24)=1
  end if

  if (print_xmol)  call xmol_output("xmol_start")

     geo_loop = int(-iwork,kind=i4_kind)
     write(char_geo_loop,'(i4)') geo_loop

    ! default value, see also the "if" below:
    analitic_hessian_calculated = .false.
#ifndef FPP_OPTIMIZER
    if( update_hessian_iteration > 0 )then
       analitic_hessian_calculated = mod(step_counter,update_hessian_iteration).eq.0
       DPRINT 'since update_hessian_iteration = ',update_hessian_iteration
       DPRINT 'and step_counter = ',step_counter
       DPRINT 'setting analitic_hessian_calculated to ',analitic_hessian_calculated
    endif
#endif

  io_flepo=openget_iounit(form='formatted',status='unknown', &
                          file=trim(opt_dir)//'/flepo.'//adjustl(char_geo_loop))
   allocopt_stat(25)=0

  ! re-set the global from opt_data_module:
  write(OPT_STDOUT,*) "main_opt: the rest of output goes to flepo ..."
  old_stdout = OPT_STDOUT ! save, to restore later
  OPT_STDOUT = io_flepo
  write (OPT_STDOUT, *) "--- Optimizer running ---"

  io_gxstore=openget_iounit(status='unknown',form='formatted', &
       file=trim(opt_dir)//'/gx.in.'//adjustl(char_geo_loop))
   allocopt_stat(23)=0
     call coordinates_write(io_gxstore,iwork)
  DPRINT 'main_opt: coordinates_write done'

  if(gx_test)then
    call make_gx_test()
    call returnclose_iounit(io_gxstore)
   allocopt_stat(23)=1
    goto 9999 ! finalize and RETURN
  endif

  DPRINT 'main_opt: main optimization block',optimization,line_search2

  if(.not.convert_internals) then
    if (optimization.or.calc_hessian.or.ts_search.or.calc_cart_hess.or.calc_cart_grad) then
      if (optimization .and. line_search2) then
         if(mod(geo_loop,2)==0) then
            read(io_gx,'(2F24.12)') energy,energy_ph
            write(io_gxstore,'(2F24.12)') energy,energy_ph
         else
            call gradient_read(io_gx,IO_GXSTORE=io_gxstore)                 !!!!!!!!!!! (2)
         end if
      else if(calc_cart_grad) then
         if(geo_loop==1) then
            call gradient_read(io_gx,IO_GXSTORE=io_gxstore)
         else
            read(io_gx,'(2F24.12)') energy,energy_ph
            write(io_gxstore,'(2F24.12)') energy,energy_ph
         end if
      else
         if( .not. tsscan_mix )then
            DPRINT 'main_opt: call gradient_read optimization calc_hessian ts_search calc_cart_hess'
            call gradient_read(io_gx,IO_GXSTORE=io_gxstore)  !!! (3)
            DPRINT 'main_opt: gradient_read done, dervs_cartes allocated',allocated(dervs_cartes)
         endif
      end if

    elseif(frequency_calculation) then
      DPRINT '(frequency_calculation) call gradient_read(...)'
      call gradient_read(io_gx,IO_GXSTORE=io_gxstore)
    elseif(epe_forces) then
!      read(io_gx,'(2F24.12)',iostat=istat) energy,energy_ph
!      ASSERT(istat.eq.0)
!      write(io_gxstore,'(2F24.12)') energy,energy_ph
       call gradient_read(io_gx,IO_GXSTORE=io_gxstore)
    else
       call gradient_read(io_gx,IO_GXSTORE=io_gxstore)
    endif

    write(opt_stdout,'(F24.12,a13)') energy,' TOTAL ENERGY'

    call returnclose_iounit(io_gxstore)
    allocopt_stat(23)=1
    DPRINT 'main_opt: epe_forces = ', epe_forces
    if(epe_forces)  call setup_calc_ff()
    if(allocated(grad_cartes)) call print_cartesian_gradients()  ! when they are complete
  endif ! of if(.not.convert_internals)

  call returnclose_iounit(io_gxstore)
   allocopt_stat(23)=1

  if (calculate_intensities) then
    if( .not. frequency_calculation .or. .not. calc_cart_hess )then
      WARN('calculating intensities without hessian?')
    endif
    call setup_intensities()
  endif

  call returnclose_iounit(io_gx)
  allocopt_stat(22)=1

  if(linear_transit) then
   write(io_flepo,*) energy, 'Transit energy, point', step_counter+1
   converged=.false.
  if(allocated(grad_cartes)) then
   deallocate(grad_cartes,stat=allocopt_stat(3))
   ASSERT(allocopt_stat(3).eq.0)
   allocopt_stat(3)=1
  endif
  if(allocated(dervs_cartes)) then
   deallocate(dervs_cartes,stat=allocopt_stat(4))
   ASSERT(allocopt_stat(4).eq.0)
   allocopt_stat(4)=1
  endif
  if(allocated(epe)) then
   deallocate(epe,stat=allocopt_stat(16))
   ASSERT(allocopt_stat(16).eq.0)
   allocopt_stat(16)=1
  endif
   return
  endif

  DPRINT 'main_opt: now do coordinates_setup'
  ! This routine sets the internal coordinates including B-Matrix
  ! and transformations of gradients and hessian
  if (print_internals .or. convert_internals .or.&
       optimization.or.ts_search.or.calc_hessian.or.frequency_calculation) then
     if(.not.calc_hessian) then
        call coordinates_setup()
     else
        ! pass the number of the actually displaced coordinate to
        ! coordinates_setup
        if(single_step) then
           call coordinates_setup(geo_loop-1)
        else
           call coordinates_setup(geo_loop/2)
        end if
     end if
  end if
  DPRINT 'main_opt: done coordinates_setup'

  if(tsscan_mix) energy=eg_rpmix(epe_kl,epe_nucen,epe)

  if (print_internals ) then
     write(OPT_STDOUT,*)" --- (primitive) Internal Coordinates before Update --- "
     if (zmat_coordinates) then
        if (zmat_format) then
           if (print_debug) then
              call print_internal(s_prim,n_primitive,full=.true.)
           else
              call print_internal(s_prim,n_primitive,full=.false.)
           endif
        elseif (free_format) then
           if (print_debug) then
              call print_internal(s,n_internal,full=.true.)
           else
              call print_internal(s,n_internal,full=.false.)
           endif
        endif
     elseif(delocalized_coordinates) then
        if (print_debug) then
           call print_internal(s_prim,n_primitive,full=.true.)
        else
           call print_internal(s_prim,n_primitive,full=.false.)
        endif
     endif
  endif

   ! check if cartesian hessian is available from ``hesse_cartesian.dat'':
   inquire(EXIST=have_hesse_cartesian, FILE=trim(inpfile('hesse_cartesian.dat')))

   if(have_hesse_cartesian)then
     write(OPT_STDOUT,*) "main_opt: found ``hesse_cartesian.dat'' in input directory ..."
     DPRINT              "main_opt: found ``hesse_cartesian.dat'' in input directory ..."
   endif

   if(analitic_hessian_calculated.and..not.have_hesse_cartesian)then
     ABORT("no hesse_cartesian.dat")
   endif

  if (optimization.or.calc_hessian.or.ts_search )  call transform_gradient()

  if(calc_cart_hess) then ! Calculate Cartesian Hessian numerically
     call cart_step(geo_loop,converged)
  end if

  if(calc_cart_grad) then ! Calculate Cartesian Gradients numerically
     call cart_step_g(geo_loop,converged)
  end if

  ! setup an initial hessian: from the following two methods are possible:
  ! 1. setup an approximate Hessian based on force-fiells parameters
  !    This will call only some routines that setup valence coords from
  !    the cartesian geometry, give some estimates for the force constants
  !    and transform this to a cartesian Hessian, which in turn will be
  !    transformed to the actual set of internals by 'hesse_main'.
  ! 2. Make a symmetry-compatible step for each internal variable and
  !    calculate estimates for the Hessian in the following step, yielding
  !    a diagonal Hessian in the actual set of internal coordinates.

  do_step=.true.

  if (( (estimate_hessian.or.estimate_hessian_m) .and.geo_loop==1).and. &
       .not.valence_format) then
     call valence_setup()
  elseif (step_technique.and.geo_loop<=2) then
     call hesse_eval(geo_loop,do_step)
     !depending on whether this is the first or second geometry loop
     !a real newton- or linear step has to be made
  endif
  DPRINT 'main_opt: do_step finaly defined',do_step

  ! either initialize hesse or read the (updated) hesse matrix
  ! from file and update it ---------------------------------

  DPRINT 'main_opt: either initialize hesse or read the (updated) hesse matrix'
  if(.not.calc_hessian.and.frequency_calculation)then
     DPRINT 'hesse_main for frequency_calculation'
     call hesse_main(geo_loop,hesse_complete)       !(1)
  endif

  hsm1:if (do_step .and. (optimization.or.calc_hessian.or.ts_search) ) then
     DPRINT 'hsm1 block'
     if (optimization .and. line_search2) then
        if (mod(geo_loop,2) /= 0) then
           DPRINT 'hesse_main 1'
           call hesse_main(geo_loop,hesse_complete)   !(2) do_step
           do_newton=.true.
        else
           hesse_complete=.true.
        end if
     else
        DPRINT 'hesse_main 2'
        call hesse_main(geo_loop,hesse_complete)      !(3) do_step
        do_newton=.true.
     end if
     DPRINT 'hesse_complete defined',hesse_complete

     if (optimization.and.(line_search.or.line_search2).and.hesse_complete) then
     DPRINT 'do line_search_main'
        if (line_search2) then
           call line_search_main2(do_newton,geo_loop)
        else
           call line_search_main(do_newton,geo_loop)
        endif
     end if

     if (optimization.and.do_newton.and.hesse_complete) then
     DPRINT 'do a quasi-newton step'
        call newton_step(geo_loop)
     endif
     if (hesse_complete.and.ts_search) then
     DPRINT 'do a  step to TS'
        call ts_main(geo_loop)
     endif

  endif hsm1
  DPRINT ' hsm1 passed'


  if(allocated(hesse_sphere)) then
   deallocate(hesse_sphere,stat=allocopt_stat(15))
   ASSERT(allocopt_stat(15).eq.0)
   allocopt_stat(15)=1
  endif

  if(allocated(grad_sphere)) then
           deallocate(grad_sphere,stat=allocopt_stat(13))
           ASSERT(allocopt_stat(13).eq.0)
           allocopt_stat(13)=1
  endif

  if (convert_internals)  call fake_internals(new_internal, sym_int, q)

  DPRINT 'update cartesian geometry'
  if(cart_coordinates) then
     if (optimization) call cart2cart()
  else
     if (optimization.or.convert_internals.or.calc_hessian.or.ts_search) &
          call internal_to_cart(q_old,q,bmat,bmat_inv) !(1)
  end if
  DPRINT 'done'

  if(.not.cart_coordinates) then
     if (print_internals .or. convert_internals .or.&
          optimization.or.ts_search.or.calc_hessian.or.frequency_calculation) &
          call free_bmat()
  end if

  if (print_xmol.and.&
       (optimization.or.ts_search.or.convert_internals)) &
       call xmol_output("xmol_update")

  if (print_internals.and.(convert_internals.or.optimization.or.&
       calc_hessian.or.ts_search)) then
     write(OPT_STDOUT,*)" --- (primitive) Internal Coordinates after Update --- "
     if (zmat_coordinates) then
        if (zmat_format) then
           if (print_debug) then
              call print_internal(s_prim,n_primitive,full=.true.)
           else
              call print_internal(s_prim,n_primitive,full=.false.)
           endif
        elseif (free_format) then
           call print_internal(s,n_internal,full=.true.)
        endif
     elseif (delocalized_coordinates) then
        call print_internal(s_prim,n_primitive,full=.true.)
     endif
  endif
  DPRINT 'print_internal passed'

  !PRELIMINARY!!!!
  if (optimization.or.convert_internals.or.calc_hessian.or.ts_search.or. &
       calc_cart_hess.or.calc_cart_grad) then
     iwork=iwork-one
     ! write a new gxfile
     io_gx=openget_iounit(status='replace',form='formatted',&
          file=trim(opt_data_dir)//'/gxfile')
   allocopt_stat(27)=0
!     if (open_stat/=0) then
!        stop 'main_opt : could not re-open gxfile'
!     endif

     call coordinates_write(io_gx,iwork)

     call returnclose_iounit(io_gx)
     allocopt_stat(27)=1

     write(char_geo_loop,'(i4)') geo_loop
     io_gx=openget_iounit(status='replace',form='formatted',&
          file=trim(opt_dir)//'/gx.out.'//adjustl(char_geo_loop))
   allocopt_stat(28)=0
         call coordinates_write(io_gx,iwork)
     call returnclose_iounit(io_gx)
   allocopt_stat(28)=1

  endif


  if (optimization.or.calc_hessian.or.ts_search) then
     !
     ! branch 1: optimization, transition state search or numerical hessian calculation ...
     !
     DPRINT 'optimizer: proceed as in optimization, transition state search or numerical hessian calculation ...'
     if (ts_search.and.hesse_complete) then
        !
        ! branch 1.1: WRITE ME!
        !
        write(OPT_STDOUT,1400)grad_max_comp,grad_mean_square
        write(OPT_STDOUT,1500)step_max_comp,step_mean_square
     else
        if (line_search2) then
           !
           ! branch 1.2.1: WRITE ME!
           !
           if(mod(geo_loop,2) /= 0) then
              write(OPT_STDOUT,1400)grad_max_comp,grad_mean_square
              write(OPT_STDOUT,1500)step_max_comp,step_mean_square
           end if
        else
           !
           ! branch 1.2.2: WRITE ME!
           !
           write(OPT_STDOUT,1400)grad_max_comp,grad_mean_square
           write(OPT_STDOUT,1500)step_max_comp,step_mean_square
        end if
     endif

     DPRINT 'hesse_complete at convergence_check blocl',hesse_complete,calc_hessian,optimization,ts_search
     if (optimization.and.hesse_complete) then
        if (line_search2) then
           if(mod(geo_loop,2) /= 0) then
              call convergence_check(converged)
           end if
        else
           call convergence_check(converged)
           DPRINT 'converged ', converged
        end if
     elseif ( calc_hessian.and. .not.optimization .and. .not.ts_search  ) then
        if (hesse_complete) then
           !
           ! branch 1.4.1: numerical hessian is complete, do the frequencies:
           !
           DPRINT 'in frequency_main block',sum(hesse)
!          write(OPT_STDOUT,*)"               --- Final Hesse Matrix --- "
!          call print_matrix(hesse,n_internal,n_internal,10)
!          write(OPT_STDOUT,*)"               --------------------------"
           converged = .true.
           DPRINT 'main_opt: now we can start frequency calculation if desired'
           if(frequency_calculation) call frequency_main(hesse)
           DPRINT 'main_opt: frequency_main passed (1)'
        endif
     elseif ( ts_search.and.hesse_complete ) then
        !
        ! branch 1.5: WRITE ME!
        !
        call convergence_check(converged)
     elseif (( ts_search.or.optimization).and. .not.hesse_complete ) then
        write(OPT_STDOUT,*)" Calculation of Hessian not yet finished"
     endif

  elseif(.not.calc_hessian.and.frequency_calculation) then
    !
    ! branch 2: analytical hessian/frequency calculation:
    !
    DPRINT 'main_opt: calculate the frequencies with available hessian'
    call frequency_main(hesse) ! (2)
    DPRINT 'main_opt: frequency_main passed, tell them we are finished'

    ! tell paragauss that there is no more action required,
    ! the task "freq_analyt" is finished:
    converged = .true.
  endif

  if(allocated(hesse)) then
   deallocate(hesse,hesse_inv,stat=allocopt_stat(14))
   ASSERT(allocopt_stat(14).eq.0)
   allocopt_stat(14)=1
  endif

  if(noline) then
   number_nolines=number_nolines+1
   if(number_nolines.eq.3) converged = .true.
  else
   number_nolines=0
  endif

  if(present(convert_internal)) then
     if(convert_internal) converged=.false.
     if(ts_scan.and.converged) then
        write(io_flepo,*) ' main_opt converged energy_prev energy', converged,energy_prev,energy
        convert_internal=energy.ge.energy_prev
        energy_prev=energy
        step_counter=0
     endif
  end if

DPRINT 'deallocate section'
  if(allocated(grad_cartes)) then
   deallocate(grad_cartes,stat=allocopt_stat(3))
   ASSERT(allocopt_stat(3).eq.0)
   allocopt_stat(3)=1
  endif
DPRINT 'deallocate grad_cartes done'

  if(allocated(dervs_cartes)) then
    deallocate(dervs_cartes,stat=allocopt_stat(4))
    ASSERT(allocopt_stat(4).eq.0)
    allocopt_stat(4)=1
  endif

   if(associated(atom)) deallocate(atom)
   if(associated(atom_reactant)) deallocate(atom_reactant)
   if(associated(atom_product)) deallocate(atom_product)


   if(.not.calc_cart_hess .and..not.cart_coordinates .and..not.calc_cart_grad) then
DPRINT 'dealloc_intcoor'
      call dealloc_intcoor(s_prim,q_prim,s,sym_type,q,q_old)
      if(tsscan_sphere) call dealloc_intcoor(s_prim_reactant,q_prim_reactant,s_reactant)
      if(tsscan_sphere) call dealloc_intcoor(s_prim_pointonmep,q_prim_pointonmep,s_pointonmep)
      if(qst_step.or.(tsscan_sphere.and.exist_product)) call dealloc_intcoor(s=s_product)
      if(tsscan_sphere.and.exist_product) call dealloc_intcoor(s_prim_product,q_prim_product)
DPRINT 'dealloc_reduce'
      call dealloc_reduce()
DPRINT 'dealloc_constraint'
      call dealloc_constraint()
   else
      call dealloc_cartcoor()
      call dealloc_constraint()
   end if
   DPRINT 'dealloc_grad_cart_to_internal'
   call dealloc_grad_cart_to_internal()
   if(allocated(epe))  then
     deallocate(epe,stat=allocopt_stat(16))
     ASSERT(allocopt_stat(16).eq.0)
     allocopt_stat(16)=1
   endif
   if(allocated(gx)) deallocate(gx)
   DPRINT 'main_opt: io_par=',io_par
   if(io_par/=-1) call returnclose_iounit(io_par)
   if(io_par_ffh/=-1) call returnclose_iounit(io_par_ffh)
   allocopt_stat(20)=1
   DPRINT 'main_opt: io_flepo=',io_flepo

   write(OPT_STDOUT,*) "optimizer: the rest of output goes to output ..."
   call returnclose_iounit(io_flepo)
   allocopt_stat(25)=1

   ! re-set the global from opt_data_module:
   OPT_STDOUT = old_stdout
   write(OPT_STDOUT,*) "optimizer: i am back again ..."

   call close_slspar()
#ifdef NEW_EPE
   call free_ff_arrays()
#endif
   if(calc_epeff_hessian.and.allocated(ewa)) then
    deallocate(ewa,stat=allocopt_stat(17))
    ASSERT(allocopt_stat(17).eq.0)
    allocopt_stat(17)=1
   endif

   call hesse_shutdown()
   call print_allocopt()
 if(calc_epeff_hessian.and..not. hesse_complete) cycle
 if(.not.tsscan_mix.or.converged) exit
 enddo internal_loop

9999 CONTINUE ! finalize and exit (should we move it higher?)

     ! duplicate persistent values on disk:
     call persistent_state('save')
     DPRINT 'done optimizer'

     RETURN

!1100 format(3(2x,f13.7))
!1200 format(20(2x,f8.4))
!1300 format('Koordinate No.',i2,2x,'    :',f13.7,5x,A20)
1400 format("max. comp. of gradient   ",F11.8,"      square mean gradient   ",F11.8)
1500 format("max. comp. of step       ",F11.8,"      square mean step       ",F11.8)

contains

  subroutine transform_gradient()
    DPRINT 'transform cartesian gradient to internals'
    if(optimization .and. line_search2) then                  !!!!!!!!!!!
       if(mod(geo_loop,2) /= 0) call grad_cart_to_internal()  !!!!!!!!!!! (1)
    else                                                      !!!!!!!!!!!
       if(cart_coordinates) then
          call grad_cart_for_opt()
       else
          call grad_cart_to_internal()                           !(2)
          if(tsscan_sphere) then
            allocate(grad_sphere(size(grad_intern)-1),stat=allocopt_stat(13))
            ASSERT(allocopt_stat(13).eq.0)
            if(exist_product) then
             grad_sphere= &
              rp_grads(tsscan_rp_var,sphere_dependent_var,grad_max_sphere,grad_mean_sphere,dEdR_sphere)
            else
             grad_sphere=sphere_grads(sphere_dependent_var,grad_max_sphere,grad_mean_sphere)
            endif
          endif
       endif
    end if
  end subroutine transform_gradient

  subroutine setup_intensities()
     ! read dipole moment from the gxfile
     ! it is later used to calculate intensities
     ! *** end of interface ***

     dipole_moment = 0.0_r8_kind
     read(io_gx,*,IOSTAT=istat) dipole_moment
     if ( istat /= 0 ) then
       print *,'OPTIMIZER: There seem to be no dipole moments in gxfile!'
       print *,'           I am setting them to zero, you will get no intensities,'
       print *,'           unless you specify OPERATIONS_DIPOLE next time.'
       WARN('WARNING: no dipole moments!')
       dipole_moment = 0.0_r8_kind
     else
       ! and write it to the file dipole.dat
       if(geo_loop==1) then
          io_dip=openget_iounit(status='unknown',form='formatted',&
           file="dipole.dat")
   allocopt_stat(29)=0
       else
          io_dip=openget_iounit(status='unknown',form='formatted',position='append',&
               file="dipole.dat")
   allocopt_stat(29)=0
       end if
       write(io_dip,'(3F20.10)') dipole_moment
       call returnclose_iounit(io_dip)
       allocopt_stat(29)=1
     endif
  end subroutine setup_intensities

 subroutine print_cartesian_gradients()
 integer(kind=i4_kind):: i,n_grads
#ifdef WITH_EFP
 integer(kind=i4_kind):: j,k
#endif

   n_grads=n_atoms+n_dummy
#ifdef WITH_EFP
   if(efp .and. n_efp > 0) n_grads=n_grads+2*n_efp
#endif
   write(opt_stdout,*) 'complete Cartesian gradients ',n_grads

   do i=1,n_atoms+n_dummy
      write(opt_stdout,'(i4,3f20.12,a19)') i,grad_cartes(i,:)," atomic   gradients"
!      write(opt_stdout,'(i4,3f20.12,a19)') i,grad_cartes(i,:)/sqrt(sum(grad_cartes(i,:)**2) )
   end do

#ifdef WITH_EFP
   if(efp .and. n_efp > 0) then
      j=n_atoms+n_dummy+1
      k=j
      do i=1,n_efp
         write(opt_stdout,'(i4,3f20.12,a19)') k,grad_cartes(j,:)," efp tran gradients"
         j=j+1
         write(opt_stdout,'(i4,3f20.12,a19)') k,grad_cartes(j,:)," efp torq gradients"
         j=j+1
         k=k+1
      enddo
   end if
#endif

  end subroutine print_cartesian_gradients

     subroutine fix_bmat_inv()
     integer(kind=i4_kind):: ibmat_inv,i
             ibmat_inv=openget_iounit(status='unknown',form='formatted',&
                                  file=trim(inpfile('bmat_inv.dat')))
             do i=1,size(bmat_inv,1)
!              write(ibmat_inv,*) bmat_inv(i,:)
              read(ibmat_inv,*) bmat_inv(i,:)
             enddo
             call returnclose_iounit(ibmat_inv)
      end subroutine fix_bmat_inv

     subroutine fix_gradients()
     integer(kind=i4_kind):: ifix_g,i
             ifix_g=openget_iounit(status='unknown',form='formatted',&
                                  file=trim(inpfile('grad.dat')))
             do i=1,size(grad_cartes,1)
!              write(ifix_g,*) grad_cartes(i,:)
              read(ifix_g,*) grad_cartes(i,:)
             enddo
             call returnclose_iounit(ifix_g)
      end subroutine fix_gradients

  subroutine setup_calc_ff()

#ifdef NEW_EPE
    real(r8_kind), allocatable :: atnm(:)
    integer(i4_kind) :: n_epe
#endif

  ff: if(epe_forces)  then

     print*,'epe_forces started',adjustl(trim(inpfile('epe.pcs')))

     inquire (file=adjustl(trim(inpfile('epe.pcs'))), exist=epe_pc_shells)

     if(epe_pc_shells) then
        print*,'file of EPE shells is found'
        io_epe_pcs=openget_iounit(file=adjustl(trim(inpfile('epe.pcs'))), &
             status='old',form='formatted')
        read (io_epe_pcs,*,err=100) epe_nucen
        print*,'number of EPE shells in epe.pcs is equal to ' &
             &        ,epe_nucen
        allocate (epe(epe_nucen+n_atoms+n_dummy),stat=allocopt_stat(16))
        ASSERT(allocopt_stat(16)==0)
        epe(:)%k=0
! first read epe%s, they are fixed for cluster geometry optimization
        do i=1, epe_nucen
           read (io_epe_pcs,*,err=100) epe(i+n_atoms+n_dummy)%s, &
                z_dummy,ind_deck(1:4) , &
                epe(i+n_atoms+n_dummy)%k, &  ! this type is define to 0 in epe.pcs presently
                epe(i+n_atoms+n_dummy)%ant   ! thus it should be set via ant as done for
        enddo ! i=1, nucen                   ! atoms of cluster
        call returnclose_iounit(io_epe_pcs)
     else
        if(.not.calc_epeff_hessian) &
        call error_handler("Main_opt: setup_calc_ff: no epe.pcs found")
     endif

        inquire (file=adjustl(trim(inpfile('epe.r'))), exist=gx_pos_regular)
        if(gx_pos_regular) then
           print*, 'file with gx regular positions found'
           io_epe_r=openget_iounit(file=adjustl(trim(inpfile('epe.r'))), &
                                               status='old',form='formatted')
        else
        if(.not.calc_epeff_hessian) &
           call error_handler("Main_opt: setup_calc_ff: no epe.r found")
        end if

           allocate(gx(n_atoms+n_dummy),stat=istat)
           ASSERT(istat==0)
           gx(:)%ieq=index_unique(1:n_atoms+n_dummy)
           gx(:)%r(1)=0.0_r8_kind
           gx(:)%r(2)=0.0_r8_kind
           gx(:)%r(3)=0.0_r8_kind
           gx(:)%q=0.0_r8_kind
           gx(:)%s(1)=x(:n_atoms+n_dummy)
           gx(:)%s(2)=y(:n_atoms+n_dummy)
           gx(:)%s(3)=z(:n_atoms+n_dummy)
           if(gx_pos_regular) print*, 'reading epe.r file'
           do i=1,n_atoms+n_dummy
            if(.not.gx_pos_regular) exit
            if(impu(i).ne.0) read(io_epe_r,*,err=101) gx(i)%r
 !         print*,'gx(i)%r', i, gx(i)%r,impu(i)
           enddo ! i=1,kolat

           if(calc_epeff_hessian) then
            io_eq=openget_iounit(file=adjustl(trim(inpfile('epe.eq'))), &
                                               status='old',form='formatted')
           do i=1,n_atoms+n_dummy
            if(impu(i).ne.0.or. &
               (calc_epeff_hessian.and.gx(i)%ieq.ne.0)) read(io_eq,*,err=103) gx(i)%eq
           enddo ! i=1,kolat
           call returnclose_iounit(io_eq)
           endif

        if(gx_pos_regular) then
           print*, 'done reading epe.r file'
           call returnclose_iounit(io_epe_r)
        endif

     epe_ff: if(epe_pc_shells) then

       epeffh: if(calc_epeff_hessian) then
           inquire (file=adjustl(trim(inpfile('ewald.pcr'))), exist=ewald_PCs)
           ewPCs: if(ewald_PCs) then
              print*, 'file with Ewald PCs is found'
              io_ewald_PCs=openget_iounit(file=adjustl(trim(inpfile('ewald.pcr'))), &
                   status='old',form='formatted')

              print*,'reading ewald.pcr'
              read(io_ewald_PCs,*) nucen_ewa
              allocate (ewa(nucen_ewa),stat=allocopt_stat(17))
              ASSERT(allocopt_stat(17).eq.0)
              do i=1,nucen_ewa
                 read(io_ewald_PCs,*) ewa(i)%r,ewa(i)%q
              enddo
              print*,'done ewald.pcr'
              call returnclose_iounit(io_ewald_PCs)

              kl=0
              n_counter=1
              do while(n_counter.le.nucen_ewa)
                 do k_counter=1,n_atoms+n_dummy
                    if(gx(k_counter)%ieq.eq.0.or.impu(k_counter).eq.0) cycle
                    if(dot_product(gx(k_counter)%r-ewa(n_counter)%r, &
                         & gx(k_counter)%r-ewa(n_counter)%r).lt.small_dist) then
                       kl=kl+1
                       pos_temp(:)=ewa(kl)%r
                       ant=ewa(kl)%q
                       ewa(kl)%r=ewa(n_counter)%r(:)
                       ewa(kl)%q=ewa(n_counter)%q
                       ewa(n_counter)%r=pos_temp(:)
                       ewa(n_counter)%q=ant
                       exit
                    endif ! small_dist
                 enddo !kolat
                 n_counter=n_counter+1
              enddo ! i=1,nucen
              ewa_low_limit=kl+1
              print *, 'number of coinciding centers',ewa_low_limit-1
           endif ewPCs

              Coul_ewa=0.0_r8_kind

              do k_counter=1,n_atoms+n_dummy
                 if(charge(k_counter)-aint(charge(k_counter)).lt.0.0099_r8_kind.or. &
                      index_unique(k_counter).eq.0) cycle

                 do n_counter=1,max_type_ions
                    if(abs(charge(k_counter)-name_of_type(n_counter)).lt.0.0099_r8_kind) &
                         gx(k_counter)%q=charge_of_type(n_counter)
                 enddo

                 V_ewa=0.0_r8_kind
                 do kl=1,n_atoms+n_dummy
                    if(charge(kl)-aint(charge(kl)).lt.0.0099_r8_kind.or. &
                         index_unique(kl).eq.0 .or. kl.eq. k_counter) cycle
                    do n_counter=1,max_type_ions
                       if(abs(charge(kl)-name_of_type(n_counter)).lt.0.001_r8_kind) then
                          V_ewa=V_ewa+charge_of_type(n_counter)/&
                               sqrt(dot_product(gx(kl)%s-gx(k_counter)%s,gx(kl)%s-gx(k_counter)%s))
                          Coul_ewa=Coul_ewa+0.5_r8_kind*gx(k_counter)%q*charge_of_type(n_counter)/&
                               sqrt(dot_product(gx(kl)%s-gx(k_counter)%s,gx(kl)%s-gx(k_counter)%s))
                       endif
                    enddo
                 enddo

           if(ewald_PCs) then
                 do n_counter=ewa_low_limit,nucen_ewa
                    V_ewa=V_ewa+ewa(n_counter)%q/&
                         sqrt(dot_product(ewa(n_counter)%r-gx(k_counter)%s, &
                         ewa(n_counter)%r-gx(k_counter)%s))
                    Coul_ewa=Coul_ewa+gx(k_counter)%q*ewa(n_counter)%q/&
                         sqrt(dot_product(ewa(n_counter)%r-gx(k_counter)%s,&
                         ewa(n_counter)%r-gx(k_counter)%s))
                 enddo
           endif !ewald_PCs
        print*,k_counter, V_ewa
        enddo
        print* ,'Coul_ewa ',Coul_ewa

        endif epeffh


! read  epe%pcr in epe%r
        inquire (file=adjustl(trim(inpfile('epe.pcr'))), &
             exist=epe_pos_regular)

        if(epe_pos_regular) then
           print*, 'file with EPE regular positions is found'
           io_epe_pcr=openget_iounit(file=adjustl(trim(inpfile('epe.pcr'))), &
                status='old',form='formatted')
           read(io_epe_pcr,*,err=102) nucen_pcr
           if(nucen_pcr.ne.epe_nucen)  &
                stop 'simol_inp: nucen_pcr.ne.epe_nucen'
           do i=1,nucen_pcr
              read(io_epe_pcr,*,err=102) epe(i+n_atoms+n_dummy)%r
           enddo ! nucen_pcr
           print*, ' epe.pcr read'
           call returnclose_iounit(io_epe_pcr) !
        else
           call error_handler("Main_opt: setup_calc_ff: no epe.pcr found")
        endif ! epe_pos_regular

        epe(1:n_atoms+n_dummy)%ant=charge(1:n_atoms+n_dummy)
        epe(1:n_atoms+n_dummy)%s(1)=x(1:n_atoms+n_dummy)
        epe(1:n_atoms+n_dummy)%s(2)=y(1:n_atoms+n_dummy)
        epe(1:n_atoms+n_dummy)%s(3)=z(1:n_atoms+n_dummy)
 !      print*,'first epe positions are set to gxfile values'



        n_counter=1+n_atoms+n_dummy
        ! counter in epe which contains in optimizer
        ! first all atoms of gxfile and then all atoms
        ! of epe environment
        kl=0
        do while(n_counter.le.epe_nucen+n_atoms+n_dummy)
           kl1=kl+1
           do k_counter=kl1,n_atoms+n_dummy
              if(gx(k_counter)%ieq.eq.0.or.impu(k_counter).eq.0) cycle
                ! only atoms of lattice are collected in gx first
              if(dot_product(gx(k_counter)%r-epe(n_counter)%r, &
                   & gx(k_counter)%r-epe(n_counter)%r).lt.small_dist) then
                ! in standart run epe-environment does not contain
                ! coinciding with gx-centers atoms therefore this reorder
                ! is actually dummy
                 kl=kl+1
                 pos_temp(:)=gx(kl)%r
                 ieqt=gx(kl)%ieq
                 gx(kl)%r=gx(k_counter)%r
                 epe(kl)%r=gx(kl)%r
                 epe(kl)%ant=charge(k_counter)
                 gx(kl)%ieq=gx(k_counter)%ieq
                 gx(k_counter)%r=pos_temp(:)
                 gx(k_counter)%ieq=ieqt

                 ieqt=epe(kl+n_atoms+n_dummy)%k
                 pos_temp(:)=epe(kl+n_atoms+n_dummy)%r(:)
                 ant=epe(kl+n_atoms+n_dummy)%ant
                 epe(kl+n_atoms+n_dummy)%r(:)=epe(n_counter)%r(:)
                 epe(kl+n_atoms+n_dummy)%ant=epe(n_counter)%ant
                 epe(kl+n_atoms+n_dummy)%k=epe(n_counter)%k
                 epe(n_counter)%r(:)=pos_temp(:)
                 epe(n_counter)%ant=ant
                 epe(n_counter)%k=ieqt

                 pos_temp(:)=epe(kl+n_atoms+n_dummy)%s(:)
                 epe(kl+n_atoms+n_dummy)%s(:)=epe(n_counter)%s(:)
                 epe(n_counter)%s(:)=pos_temp(:)
                 exit
              endif ! small_dist
           enddo !kolat
           n_counter=n_counter+1
        enddo ! i=1,nucen

        epe_kl=kl+1+n_atoms+n_dummy ! indicate first atom of environment
        DPRINT 'epe centers go from number' ,epe_kl,kl


     if(.not.tsscan_mix) then
#ifdef NEW_EPE
        n_epe=size(epe)
        allocate(atnm(n_epe),stat=istat)
        ASSERT(istat==0)
        do i=1,n_epe
           atnm(i)=epe(i)%ant
        end do
        call read_interface_param(n_epe,atnm,n_atoms+n_dummy)
#endif

        if(n_types_central_atoms_3body.ne.0) then
           ks: do k_counter=1,n_atoms+n_dummy + epe_nucen
                                              ! recalculated instead read
                                              ! why they eq 0 in epe.pcs
              do n_counter=1,max_type_ions
                 if(abs(epe(k_counter)%ant-name_of_type(n_counter))&
                      .lt.0.00001_r8_kind)  then
                    epe(k_counter)%k=n_counter
                    exit
                 end if
              end do
           end do ks

           DPRINT 'call building_tet'
           call building_tet(n_atoms+n_dummy)
        endif
!           print*, 'number of EPE coinciding centers in cluster',kl

        if(calc_epeff_hessian) then
 !      print*, 'add_epe_f'
           if(epe_pc_shells) call add_epe_f(n_atoms+n_dummy,'ffh')
 !      print*, 'done'
        else
           if(epe_pc_shells) call add_epe_f(n_atoms+n_dummy)
        endif

       if(n_types_central_atoms_3body > 0) then
           energy=energy+energy_3_body(n_atoms+n_dummy)
           call gradients_3_body(n_atoms+n_dummy)
       endif

     DPRINT 'done forces due to external links'
     else
        print*, 'no file with epe shells'
     endif


#ifdef NEW_EPE
     deallocate(atnm,stat=istat)
     ASSERT(istat==0)
#endif

  endif  epe_ff

    if(.not.tsscan_mix) then
     if(calc_epeff_hessian) then
 !      print*,'add_sls_f'
        call add_sls_f(n_atoms+n_dummy,'ffh')
 !      print*,'done add_sls_f'
        call gradients_bending3b(n_atoms+n_dummy)
     else
        call add_sls_f(n_atoms+n_dummy)
     endif
    endif

     print*, energy,' EPE corrected total energy'
     write(opt_stdout,*) energy,' EPE corrected total energy'
     write(output_unit,*)  energy,' EPE corrected total energy' !!!!!!for processing by Viewmol

    if(fixed_orientation) write(opt_stdout,*) 'Fixed orientation'

    if(allocated(grad_cartes)) then

     do i=1,n_atoms+n_dummy
        DPRINT i,grad_cartes(i,:)
             write(output_unit,'(3f20.12,a14)') grad_cartes(i,:)," epe gradients"
     enddo ! i=1,kolat
     endif ! grad_cartes allocated

  endif  ff

  return

100 call error_handler("Main_opt: setup_calc_ff: error while reading in epe.pcs")
101 call error_handler("Main_opt: setup_calc_ff: error while reading in epe.r")
102 call error_handler("Main_opt: setup_calc_ff: error while reading in epe.pcr")
103 call error_handler("Main_opt: setup_calc_ff: error while reading in epe.eq")

  end subroutine setup_calc_ff

  subroutine make_gx_test()
     ! check if the gxfile contains
     ! redundant coordinates
     ! and if it preserves symmetry
     if(single_step) then
        call coordinates_setup(geo_loop-1)
     else
        call coordinates_setup(geo_loop/2)
     end if

     call gradient_read(io_gx,zero=.true.) !!! (1) gx_test
     call returnclose_iounit(io_gx)

     call grad_cart_to_internal()          !(3) gx_test
     call free_bmat()
     call frequency_main(hesse)
     call alloc_bmat(n_primitive)
     call hesse_main(geo_loop,hesse_complete) !(0) gx_test
     call internal_to_cart(q_old,q,bmat,bmat_inv) ! (2) test
     iwork=iwork-one
     ! write a new gxfile
       call returnclose_iounit(io_gx)
     io_gx=openget_iounit(status='replace',form='formatted',&
          file=trim(opt_data_dir)//'/gxfile')

!!$          file=adjustl(trim(opt_data_dir))//'/gxfile',iostat=open_stat)
!     if (open_stat/=0) then
!        stop 'main_opt : could not re-open gxfile'
!     endif

     call coordinates_write(io_gx,iwork)

     call returnclose_iounit(io_gx)
!     if (open_stat/=0) then
!        print*,'main_opt : warning - the gxfile could not be closed again'
!     endif! --

     io_gx=openget_iounit(status='replace',form='formatted',&
          file=trim(opt_dir)//'/gx.out.'//adjustl(char_geo_loop))
         call coordinates_write(io_gx,iwork)
     call returnclose_iounit(io_gx)

     if(hesse_complete) then
        converged = .true.
        write(OPT_STDOUT,*) 'Congratulations for choosing this gxfile. All tests have been passed.'
     endif
!    ! FIXME:
!    print *,'main_opt: stop!'
!    stop
   end subroutine make_gx_test

  subroutine add_sls_f(kolat,ffh)
    !       purpose: adds SLS forces for intra-cluster interactions
        use coortype_module
    implicit none
    integer(kind=i4_kind), intent(in):: kolat
    character(len=3),optional,intent(in):: ffh
    ! *** end of interface ***

    real(r8_kind), dimension(3):: rr
    real(r8_kind)::sls_e,sls_g,k_par,r0
    real(r8_kind)::sls_d
    integer(kind=i4_kind)::k
    integer(kind=i4_kind)::n1,n2
    real(kind=r8_kind)::rss1,rss2,rss4,rss6,rss8,rss10
    real(kind=r8_kind), dimension(3):: rr_n1,rr_n2
#ifdef NEW_EPE
    integer(i4_kind) :: it,jt
    real(r8_kind) :: b,ro,c,d,cut
#else
    real(kind=r8_kind)::mn,mx
    integer(kind=i4_kind)::nr
#endif


    sls_e=0.0_r8_kind
    do n1=2,kolat
     rr_n1(1)=x(n1)
     rr_n1(2)=y(n1)
     rr_n1(3)=z(n1)


       if(charge(n1)-aint(charge(n1)).lt.0.0009_r8_kind &
                .or. index_unique(n1).eq.0) cycle
#ifdef NEW_EPE
       it=at_type(n1)
#endif

       do n2=1,n1-1
        rr_n2(1)=x(n2)
        rr_n2(2)=y(n2)
        rr_n2(3)=z(n2)
          if(charge(n2)-aint(charge(n2)).lt.0.000999_r8_kind.or. &
               index_unique(n2).eq.0) cycle
#ifdef NEW_EPE
          jt=at_type(n2)
          b=interface_ff(it,jt)%b
          ro=interface_ff(it,jt)%r
          c=interface_ff(it,jt)%c
          d=interface_ff(it,jt)%d
          cut=interface_ff(it,jt)%cut
          if(abs(cut) > 0.1_r8_kind) then
#else
          mn=min(charge(n1),charge(n2))
          mx=max(charge(n1),charge(n2))

          if(present(ffh)) then
 !           print*,'get_ffhp'
             call get_ffhp(mn,mx,0,nr)
             sls_e=sls_e+gx(n1)%q*gx(n2)%q/&
                  sqrt(dot_product(gx(n1)%s-gx(n2)%s,gx(n1)%s-gx(n2)%s))
             sls_g=gx(n1)%q*gx(n2)%q/&
                  sqrt(dot_product(gx(n1)%s-gx(n2)%s,gx(n1)%s-gx(n2)%s))**3
             grad_cartes(n1,:)=grad_cartes(n1,:)-(gx(n1)%s-gx(n2)%s)*sls_g
             grad_cartes(n2,:)=grad_cartes(n2,:)+(gx(n1)%s-gx(n2)%s)*sls_g
 !           print*,'charge contrib added'
          else
             call get_slsp(mn,mx,0,nr)
          endif
 !        print*,'nr=',nr
          if(nr.ne.0) then
#endif
             sls_g=0.0
             sls_d=0.0
             rss1=sqrt(dot_product(rr_n1-rr_n2,rr_n1-rr_n2))*angsau
             rss2=rss1**2
#ifdef NEW_EPE
           if(rss2.le.cut**2) then
#else
           if(rss2.le.par%item%cutoff**2) then
#endif
                rss4=rss2**2
                rss6=rss4*rss2
                rss8=rss4**2
                rss10=rss8*rss2
#ifdef NEW_EPE
                sls_e=sls_e-evau*(c/rss6-d/rss8)
                sls_g=sls_g+gepe_const*(6.0_r8_kind*c/rss8+8.0_r8_kind*d/rss10)
                sls_e=sls_e+evau*b*exp(-rss1/ro)
                sls_g=sls_g-gepe_const*b/(ro*rss1)*exp(-rss1/ro)
                if(analitic_hessian_calculated) then
                   sls_d=sls_d-gepe_const*angsau**2* &
                        (48.0_r8_kind*c/rss10+80.0_r8_kind*d/rss10/rss2)
                   sls_d=sls_d+gepe_const*angsau**2*b/(ro*rss1)**2*exp(-rss1/ro)
                   sls_d=sls_d+gepe_const*angsau**2*b/ro/rss2/rss1*exp(-rss1/ro)
                endif
#else
                sls_e=sls_e-evau*(par%item%c/rss6-par%item%d/rss8)
                sls_g=sls_g+gepe_const*(6.0_r8_kind*par%item%c/rss8 &
                     &                   +8.0_r8_kind*par%item%d/rss10)
                sls_e=sls_e+evau*par%item%b*exp(-rss1/par%item%r)
                sls_g=sls_g- &
                &  gepe_const*par%item%b/(par%item%r*rss1)*exp(-rss1/par%item%r)
            if(analitic_hessian_calculated.or.update_fromcartessian) then
              sls_d=sls_d-gepe_const*angsau**2* &
                (48.0_r8_kind*par%item%c/rss10+80.0_r8_kind*par%item%d/rss10/rss2)
              sls_d=sls_d+gepe_const*angsau**2*par%item%b/(par%item%r*rss1)**2*exp(-rss1/par%item%r)
              sls_d=sls_d+gepe_const*angsau**2*par%item%b/par%item%r/rss2/rss1*exp(-rss1/par%item%r)
            endif
#endif
            endif

!    if(n1.eq.24) print*,mn,mx,nr,par%item%k,rss1,par%item%r0/angsau
!    now harmonic pair potential terms

#ifndef NEW_EPE
          if(par%item%k.ne.0.0 .and. abs(rss1-par%item%r0).lt.0.1) then
              sls_e=sls_e+par%item%k*(rss1/angsau-par%item%r0/angsau)**2
              sls_g=sls_g+2.0_r8_kind*par%item%k*(rss1-par%item%r0)/rss1
           endif

           if(par%item%k1.ne.0.0 .and. abs(rss1-par%item%r1).lt.0.1) then
              sls_e=sls_e+par%item%k1*(rss1/angsau-par%item%r1/angsau)**2
              sls_g=sls_g+2.0_r8_kind*par%item%k1*(rss1-par%item%r1)/rss1
           endif
#endif

             rr(1)=(x(n1)-x(n2))
             rr(2)=(y(n1)-y(n2))
             rr(3)=(z(n1)-z(n2))


        if(allocated(grad_cartes))  then
             grad_cartes(n1,:)=grad_cartes(n1,:)+rr(:)*sls_g
             grad_cartes(n2,:)=grad_cartes(n2,:)-rr(:)*sls_g
        endif

        if(.false..and.analitic_hessian_calculated) then
         forall(k=1:3)
          dervs_cartes(n1,n1,k,k)=dervs_cartes(n1,n1,k,k)+sls_g
          dervs_cartes(n2,n2,k,k)=dervs_cartes(n2,n2,k,k)+sls_g
          dervs_cartes(n1,n2,k,k)=dervs_cartes(n1,n2,k,k)-sls_g
          dervs_cartes(n2,n1,k,k)=dervs_cartes(n2,n1,k,k)-sls_g
         end forall
          dervs_cartes(n1,n1,:,:)=dervs_cartes(n1,n1,:,:) &
                    + sls_d*spread(rr,1,3)*spread(rr,2,3)
          dervs_cartes(n2,n2,:,:)=dervs_cartes(n2,n2,:,:) &
                    + sls_d*spread(rr,1,3)*spread(rr,2,3)
          dervs_cartes(n1,n2,:,:)=dervs_cartes(n1,n2,:,:) &
                    - sls_d*spread(rr,1,3)*spread(rr,2,3)
          dervs_cartes(n2,n1,:,:)=dervs_cartes(n2,n1,:,:) &
                    - sls_d*spread(rr,1,3)*spread(rr,2,3)
        endif
          elseif(present(ffh)) then

             sls_g=0.0
             sls_d=0.0
             rss1=sqrt(dot_product(rr_n1-rr_n2,rr_n1-rr_n2))
             rr(:)=rr_n1-rr_n2

 !        print*,'k_parameter'
          k_par=k_parameter(r0,n1,n2,type=b_length)
          sls_e=sls_e+K_par*(rss1-r0)**2
              sls_g=sls_g+2*k_par*(rss1-r0)/rss1
          grad_cartes(n1,:)=grad_cartes(n1,:)+rr(:)*sls_g
          grad_cartes(n2,:)=grad_cartes(n2,:)-rr(:)*sls_g
          endif ! nr.ne.0

       enddo ! n=2,kolat
    enddo ! n=1,kolat
    energy=energy+sls_e
  end subroutine add_sls_f


        subroutine add_epe_f(kolat,ffh)
         ! **** add energy contributions and forces due to liks outside cluster
        use coortype_module

        integer(kind=I4_kind), intent(in):: kolat
        character(len=3),optional,intent(in):: ffh

        real(kind=r8_kind), dimension(3):: rr,rr_n1
        integer(kind=I4_kind)::n1,n2
        real(kind=r8_kind)::sls_e,sls_g,rss1,rss2,rss4,rss6,rss8,rss10
        real(kind=r8_kind)::k_par,r0
#ifdef NEW_EPE
        integer(kind=I4_kind) :: it,jt
        real(r8_kind) :: b,ro,c,d,cut
#else
        real(kind=r8_kind)::mn,mx
        integer(kind=I4_kind)::nr
#endif

        integer(kind=I4_kind)::k
        real(kind=r8_kind)::sls_d

      sls_e=0.0_r8_kind

      atoms:  do n1=1,kolat
 !    DPRINT 'n1=',n1
        rr_n1(1)=x(n1)
        rr_n1(2)=y(n1)
        rr_n1(3)=z(n1)

        if(charge(n1)-aint(charge(n1)).lt.0.00099_r8_kind.or. &
             index_unique(n1).eq.0) cycle

#ifdef NEW_EPE
        it=at_type(n1)
#endif

        if(present(ffh).and.ewald_PCs) then
          DPRINT 'add coulomb contribution',ewa_low_limit,nucen_ewa
          grad_cartes(n1,:)=0.0_r8_kind
          do n2=ewa_low_limit,nucen_ewa
           rss1=sqrt(dot_product(gx(n1)%s-ewa(n2)%r,gx(n1)%s-ewa(n2)%r))
           if(rss1.lt.0.01_r8_kind) cycle
           sls_e=sls_e+gx(n1)%q*ewa(n2)%q/rss1
           sls_g=gx(n1)%q*ewa(n2)%q/rss1**3
           grad_cartes(n1,:)=grad_cartes(n1,:)-(gx(n1)%s-ewa(n2)%r)*sls_g
          enddo

        endif

        DPRINT 'calculate FF contribs for atom',n1
        do n2=epe_kl,epe_nucen+kolat ! loop over epe-environment

        if(epe(n2)%ant-aint(epe(n2)%ant).lt.0.009999_r8_kind) cycle
#ifdef NEW_EPE
        jt=at_type(n2)
        b=interface_ff(it,jt)%b
        ro=interface_ff(it,jt)%r
        c=interface_ff(it,jt)%c
        d=interface_ff(it,jt)%d
        cut=interface_ff(it,jt)%cut
        if(abs(cut) >= 0.1_r8_kind) then
#else
        mn=min(charge(n1),epe(n2)%ant)
        mx=max(charge(n1),epe(n2)%ant)

       if(present(ffh)) then
         call get_ffhp(mn,mx,0,nr)
       else
         call get_slsp(mn,mx,0,nr)
       endif

        if(nr.ne.0) then
#endif

           sls_g=0.0_r8_kind
           sls_d=0.0_r8_kind

           rr=rr_n1-epe(n2)%s(:)
           rss2=dot_product(rr,rr)
           rss1=sqrt(rss2)
#ifndef NEW_EPE
           if(rss1*angsau.le.par%item%rqq) then
              sls_e=sls_e- par%item%qq/rss1
              sls_g=par%item%qq/rss2/rss1
            if(analitic_hessian_calculated.or.update_fromcartessian)  &
               sls_d=sls_d-3.0_r8_kind*par%item%qq/rss2**2/rss1
           endif
!           DPRINT 'dist to',n2,  rss1, sls_g, 'force on n1'

           grad_cartes(n1,:)=grad_cartes(n1,:)+rr(:)*sls_g
#endif

           ! from now on r is in angsterms
           rss1=rss1*angsau
           rss2=rss1**2


           sls_g=0.0
#ifdef NEW_EPE
           if(rss2.le.cut**2) then
#else
           if(rss2.le.par%item%cutoff**2) then
#endif
              rss4=rss2**2
              rss6=rss4*rss2
              rss8=rss4**2
              rss10=rss8*rss2
#ifdef NEW_EPE
              sls_e=sls_e-evau*(c/rss6-d/rss8)
              sls_g=sls_g+gepe_const*(6.0_r8_kind*c/rss8+8.0_r8_kind*d/rss10)
              sls_e=sls_e+evau*b*exp(-rss1/ro)
              sls_g=sls_g-gepe_const*b/(ro*rss1)*exp(-rss1/ro)
              if(analitic_hessian_calculated) then
                 sls_d=sls_d-gepe_const*angsau**2* &
                      (48.0_r8_kind*c/rss10+80.0_r8_kind*d/rss10/rss2)
                 sls_d=sls_d+gepe_const*angsau**2*b/(ro*rss1)**2*exp(-rss1/ro)
                 sls_d=sls_d+gepe_const*angsau**2*b/ro/rss2/rss1*exp(-rss1/ro)
              endif
#else
              sls_e=sls_e-evau*(par%item%c/rss6-par%item%d/rss8)
              sls_g=sls_g+gepe_const*(6.0_r8_kind*par%item%c/rss8 &
                   &          +8.0_r8_kind*par%item%d/rss10)
              sls_e=sls_e+evau*par%item%b*exp(-rss1/par%item%r)
              sls_g=sls_g- gepe_const*par%item%b/(par%item%r*rss1)*exp(-rss1/par%item%r)
            if(analitic_hessian_calculated.or.update_fromcartessian) then
              sls_d=sls_d-gepe_const*angsau**2* &
                (48.0_r8_kind*par%item%c/rss10+80.0_r8_kind*par%item%d/rss10/rss2)
              sls_d=sls_d+gepe_const*angsau**2*par%item%b/(par%item%r*rss1)**2*exp(-rss1/par%item%r)
              sls_d=sls_d+gepe_const*angsau**2*par%item%b/par%item%r/rss2/rss1*exp(-rss1/par%item%r)
            endif
#endif
           endif! rss2.le.rcuts

#ifndef NEW_EPE
          if(par%item%k.ne.0.0 .and. abs(rss1-par%item%r0).lt.0.1) then
           DPRINT 'harmonic_bond_fix', rss1
          sls_e=sls_e+par%item%k*(rss1/angsau-par%item%r0/angsau)**2
              sls_g=sls_g+2.0_r8_kind*par%item%k*(rss1-par%item%r0)/rss1
          endif

         if(par%item%k1.ne.0.0 .and. abs(rss1-par%item%r1).lt.0.1) then
         sls_e=sls_e+par%item%k1*(rss1/angsau-par%item%r1/angsau)**2
              sls_g=sls_g+2.0_r8_kind*par%item%k1*(rss1-par%item%r1)/rss1
         endif
#endif

        if(allocated(grad_cartes)) &
           grad_cartes(n1,:)=grad_cartes(n1,:)+rr(:)*sls_g

        if(analitic_hessian_calculated.or.update_fromcartessian) then
         do k=1,3
          dervs_cartes(n1,n1,k,k)=dervs_cartes(n1,n1,k,k)+sls_g
          dervs_cartes(n1,n1,:,:)=dervs_cartes(n1,n1,:,:) &
                    + sls_d*spread(rr,1,3)*spread(rr,2,3)
         enddo
        endif

        else
         if(.false..and.present(ffh)) then
!          print*,'FF parameters are not found',charge(n1),epe(n2)%ant
           sls_g=0.0_r8_kind
           sls_d=0.0_r8_kind

           rr=rr_n1-epe(n2)%s(:)
           rss2=dot_product(rr,rr)
           rss1=sqrt(rss2)

          k_par=k_parameter(r0,n1,n2,type=b_length_toepe)
          sls_e=sls_e+K_par*(rss1-r0)**2
              sls_g=sls_g+2*k_par*(rss1-r0)/rss1
          grad_cartes(n1,:)=grad_cartes(n1,:)+rr(:)*sls_g

         endif
!!$              stop
        endif                   ! nr.ne.0
        enddo ! n=2,kolat

#if 0
      if(analitic_hessian_calculated) then
      print*,'grad_cartes dervs_cartes'
      print*,grad_cartes(n1,:)
      print*,' '
      do k=1,3
      print*,dervs_cartes(n1,n1,k,:)
      enddo
      endif
#endif

        enddo atoms

        energy=energy+sls_e
      end subroutine add_epe_f

      subroutine gradients_bending3b(kolat)
        use coortype_module
        use gradient_module, only:grad_cartes

        real(kind=r8_kind):: rss1,rss2,scal,cos_th,sin_th,g3b,theta
        integer(kind=i4_kind),intent(in)::kolat
        integer(kind=i4_kind)::i,ia_1,ia_2,ia_3,ia1,ia2
        real(kind=r8_kind) :: e1(3),e2(3)
        real(kind=r8_kind) :: theta_0,minrss1,k_par

 !         print*,'kolat', kolat, size(gx)
           do ia_2=1,kolat
            if(impu(ia_2).ne.0.or.gx(ia_2)%ieq.eq.0) cycle
            if(charge(ia_2)-aint(charge(ia_2)).lt.00001) cycle
            minrss1=100000.0_r8_kind
            do i=1,kolat
             if(impu(i).eq.0.or.gx(i)%ieq.eq.0) cycle
            if(charge(i)-aint(charge(i)).lt.00001) cycle
             rss1=dot_product(gx(ia_2)%eq-gx(i)%eq, &
                              gx(ia_2)%eq-gx(i)%eq)
             if(rss1.ge.minrss1) cycle
             minrss1=rss1
             ia_1=i
            enddo
 !          print*,ia_1,minrss1,ia_2,impu(ia_2)
            rss1=sqrt(dot_product(gx(ia_2)%s-gx(ia_1)%s, &
                             gx(ia_2)%s-gx(ia_1)%s))
            do ia_3=1,kolat
            if(ia_3.eq.ia_1.or.impu(ia_3).eq.0.or.gx(ia_3)%ieq.eq.0) cycle
            if(charge(ia_3)-aint(charge(ia_3)).lt.00001) cycle

                    rss2=dot_product(gx(ia_3)%s-gx(ia_1)%s, &
                                     gx(ia_3)%s-gx(ia_1)%s)
                    scal=dot_product(gx(ia_2)%s-gx(ia_1)%s, &
                                     gx(ia_3)%s-gx(ia_1)%s)

                    rss2=sqrt(rss2)

                    k_par=k_parameter(theta_0,ia_2,ia_1,ia_3,type=b_angle)
                    if(k_par.lt.0.000001) cycle

                    cos_th=scal/(rss1*rss2)
                    if(cos_th.gt.0.9999) cos_th=0.9999
                    if(cos_th.lt.-0.9999) cos_th=-0.9999
                    theta=acos(cos_th)
                    sin_th=sin(theta)
                    g3b=k_par*(theta-theta_0)
!                    print*,'theta theta_0 g3b',theta,theta_0,g3b
                    e1=(gx(ia_2)%s-gx(ia_1)%s)/rss1
                    e2=(gx(ia_3)%s-gx(ia_1)%s)/rss2
                         grad_cartes(ia_3,:)=grad_cartes(ia_3,:)+ &
                            g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)

                         grad_cartes(ia_2,:)=grad_cartes(ia_2,:)+ &
                         g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)

                         grad_cartes(ia_1,:)=grad_cartes(ia_1,:)+ &
                         g3b*((rss1-rss2*cos_th)*e1(:) &
                         +(rss2-rss1*cos_th)*e2(:))/ &
                         (rss2*rss1*sin_th)
            enddo
            ia1=ia_2
            ia2=ia_1
            do ia_3=1,kolat
            if(ia_3.eq.ia_2.or.impu(ia_3).ne.0.or.gx(ia_3)%ieq.eq.0) cycle
            if(charge(ia_3)-aint(charge(ia_3)).lt.00001) cycle

                    rss2=dot_product(gx(ia_3)%s-gx(ia1)%s, &
                                     gx(ia_3)%s-gx(ia1)%s)
                    scal=dot_product(gx(ia2)%s-gx(ia1)%s, &
                                     gx(ia_3)%s-gx(ia1)%s)

                    rss2=sqrt(rss2)

                    k_par=k_parameter(theta_0,ia2,ia1,ia_3,type=b_angle)
                    if(k_par.lt.0.000001) cycle

                    cos_th=scal/(rss1*rss2)
                    if(cos_th.gt.0.9999) cos_th=0.9999
                    if(cos_th.lt.-0.9999) cos_th=-0.9999
                    theta=acos(cos_th)
                    sin_th=sin(theta)

                    g3b=k_par*(theta-theta_0)
!                    print*,'theta theta_0 g3b',theta,theta_0,g3b
                    e1=(gx(ia2)%s-gx(ia1)%s)/rss1
                    e2=(gx(ia_3)%s-gx(ia1)%s)/rss2
                         grad_cartes(ia_3,:)=grad_cartes(ia_3,:)+ &
                            g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)

                         grad_cartes(ia2,:)=grad_cartes(ia2,:)+ &
                         g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)

                         grad_cartes(ia1,:)=grad_cartes(ia1,:)+ &
                         g3b*((rss1-rss2*cos_th)*e1(:) &
                         +(rss2-rss1*cos_th)*e2(:))/ &
                         (rss2*rss1*sin_th)
            enddo
#if 1
            ia1=ia_2
            do ia_1=1,kolat
            ia2=ia_1
            if(charge(ia_1)-aint(charge(ia_1)).lt.00001) cycle
            if(ia_1.eq.ia_2.or.impu(ia_1).ne.0.or.gx(ia_1)%ieq.eq.0) cycle
            rss1=sqrt(dot_product(gx(ia_2)%s-gx(ia_1)%s, &
                             gx(ia_2)%s-gx(ia_1)%s))
            do ia_3=ia_1+1,kolat
            if(charge(ia_3)-aint(charge(ia_3)).lt.00001) cycle
            if(ia_3.eq.ia_2.or.impu(ia_3).ne.0.or.gx(ia_3)%ieq.eq.0) cycle
                    rss2=dot_product(gx(ia_3)%s-gx(ia1)%s, &
                                     gx(ia_3)%s-gx(ia1)%s)
                    scal=dot_product(gx(ia2)%s-gx(ia1)%s, &
                                     gx(ia_3)%s-gx(ia1)%s)

                    rss2=sqrt(rss2)

                    k_par=k_parameter(theta_0,ia2,ia1,ia_3,type=b_angle)
                    if(k_par.lt.0.000001) cycle

                    cos_th=scal/(rss1*rss2)
                    if(cos_th.gt.0.9999) cos_th=0.9999
                    if(cos_th.lt.-0.9999) cos_th=-0.9999
                    theta=acos(cos_th)
                    sin_th=sin(theta)
                    g3b=k_par*(theta-theta_0)
!                    print*,'theta theta_0 g3b',theta,theta_0,g3b
                    e1=(gx(ia2)%s-gx(ia1)%s)/rss1
                    e2=(gx(ia_3)%s-gx(ia1)%s)/rss2
                         grad_cartes(ia_3,:)=grad_cartes(ia_3,:)+ &
                            g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)

                         grad_cartes(ia2,:)=grad_cartes(ia2,:)+ &
                         g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)

                         grad_cartes(ia1,:)=grad_cartes(ia1,:)+ &
                         g3b*((rss1-rss2*cos_th)*e1(:) &
                         +(rss2-rss1*cos_th)*e2(:))/ &
                         (rss2*rss1*sin_th)
            enddo
            enddo
#endif

           enddo

      end subroutine gradients_bending3b

#ifdef NEW_EPE
      subroutine gradients_3_body(kolat)
        use gradient_module, only:grad_cartes
        real(kind=r8_kind) :: deg2rad
        real(kind=r8_kind):: rss1,rss2,scal,cos_th,sin_th,g3b,theta
        integer(kind=i4_kind),intent(in)::kolat
        integer(kind=i4_kind)::i,l,ia_1,i1,j,ia_2,k1,j1,k,ia_3
        real(kind=r8_kind) :: e1(3),e2(3)
        integer(kind=i4_kind)::num_3b_links

        if(n_types_central_atoms_3body > 0 ) then
           print*, 'EPE starts from ', epe_kl
           deg2rad=pi/180.0_r8_kind
           num_3b_links=0
           fst: do i=1,n_tetrahedrons
              do l=1,5
                 if (tetra_atoms(i,l) <= epe_nucen+kolat) goto 1
              enddo
              cycle fst
1             ia_1=tetra_atoms(i,1)
              if( ia_1<1 .or. ia_1 > size(epe) )then
                 ABORT('please fix code, remove logic!')
                 cycle fst
              endif
              i1=at_type(ia_1)
              scnd :do j=2,4
                 ia_2=tetra_atoms(i,j)
                 if(ia_2.eq.0) exit
                 j1=at_type(ia_2)
                 thrd: do k=j+1,5
                    ia_3=tetra_atoms(i,k)
                    if(ia_3.eq.0) exit
                    if(ia_1.ge.epe_kl.and.ia_2.ge.epe_kl &
                         .and.ia_3.ge.epe_kl) cycle thrd
                    ! all atoms  outside cluster
                    k1=at_type(ia_3)
                    rss1=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                         epe(ia_2)%s-epe(ia_1)%s)
                    rss2=dot_product(epe(ia_3)%s-epe(ia_1)%s, &
                         epe(ia_3)%s-epe(ia_1)%s)
                    scal=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                         epe(ia_3)%s-epe(ia_1)%s)

                    rss1=sqrt(rss1)
                    rss2=sqrt(rss2)
                    cos_th=scal/(rss1*rss2)
                    theta=acos(cos_th)
                    sin_th=sin(theta)
                    g3b=ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
                    e1=(epe(ia_2)%s-epe(ia_1)%s)/rss1
                    e2=(epe(ia_3)%s-epe(ia_1)%s)/rss2
                    num_3b_links=num_3b_links+1
                    if(ia_3 .lt.epe_kl)  then
                       grad_cartes(ia_3,:)=grad_cartes(ia_3,:)+ &
                            evau*g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
                    endif

                    if(ia_2.lt.epe_kl ) &
                         grad_cartes(ia_2,:)=grad_cartes(ia_2,:)+ &
                         evau*g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)

                    if(ia_1.lt.epe_kl) &
                         grad_cartes(ia_1,:)=grad_cartes(ia_1,:)+ &
                         evau*g3b*((rss1-rss2*cos_th)*e1(:) &
                         +(rss2-rss1*cos_th)*e2(:))/ &
                         (rss2*rss1*sin_th)
                    if(ia_1 .lt.epe_kl) then
                       DPRINT evau*g3b*((rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/(rss2*rss1*sin_th)
                    endif

                 enddo thrd!k=j+1,4
              enddo scnd!j=1,3
           enddo fst !i=1,n_tetra_atoms
        endif
      end subroutine gradients_3_body

      real(kind=r8_kind) function energy_3_body(kolat)
        implicit none
        integer(kind=i4_kind), intent(in):: kolat
        ! *** end of interface ***

        real(kind=r8_kind) :: deg2rad
        real(kind=r8_kind) :: rss1,rss2,scal1
        real(kind=r8_kind) :: theta_a,E3b

        integer(kind=i4_kind) :: i,j,k,i1,j1,k1,l
        integer(kind=i4_kind) :: ia_1,ia_2,ia_3

        deg2rad=pi/180.0_r8_kind
        e3b=zero
        fst:  do i=1,n_tetrahedrons
!!$     print*,'tetrahedron No', i
           do l=1,5
              if (tetra_atoms(i,l)<=epe_nucen+kolat) goto 1
           enddo
           cycle fst
1          ia_1=tetra_atoms(i,1)
           if(ia_1.eq.0) cycle fst
           i1=at_type(ia_1)

           scnd: do j=2,4
              ia_2=tetra_atoms(i,j)
              if(ia_2.eq.0) cycle

              j1=at_type(ia_2)
              thrd: do k=j+1,5
                 ia_3=tetra_atoms(i,k)
                 if(ia_3.eq.0) cycle
                 if(ia_1 >= epe_kl.and.ia_2 >= epe_kl.and.ia_3 >= epe_kl) cycle thrd
                 ! epe_kl number of atoms in cluster+1 , if all 3 atoms are
                 ! outside cycle
                 k1=at_type(ia_3)
                 rss1=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                      epe(ia_2)%s-epe(ia_1)%s)
                 rss2=dot_product(epe(ia_3)%s-epe(ia_1)%s, &
                      epe(ia_3)%s-epe(ia_1)%s)
                 scal1=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                      epe(ia_3)%s-epe(ia_1)%s)
                 rss1=sqrt(rss1)
                 rss2=sqrt(rss2)
                 theta_a=acos(scal1/(rss1*rss2))
                 e3b=e3b+evau*0.5*ki(j1,i1,k1)*((theta_a-theta_0(j1,i1,k1)*deg2rad)**2)
              enddo thrd!k=j+1,5
           enddo scnd!j=2,4
        enddo fst!i=1,n_tetrahedrons
        energy_3_body= e3b
        print*, energy_3_body,  ' energy_3_body'
        write(io_flepo,*) e3b, 'e3b from all contribs'
      end function energy_3_body

      subroutine building_tet(kolat)
        ! **procedure looks for tetrahedrally coordinated atoms
        ! **and their neighbours

        integer(kind=i4_kind) :: i,i1,j,k,k1,l,m,n,ind
        integer(kind=i4_kind) :: status
        integer(kind=i4_kind), intent(in)::kolat
        logical :: exit_cycle
        real(kind=r8_kind) :: dist
        real(kind=r8_kind), dimension(4) :: buf_dist
        integer(kind=i4_kind), dimension(4) :: buf_index
        integer(kind=i4_kind) :: max_ind(1)
        real(kind=r8_kind) :: max_val

        if (allocated(tetra_atoms)) then
           deallocate(tetra_atoms,stat=status)
           ASSERT(status==0)
        end if
        n_tetrahedrons=0

        do i=1,epe_nucen+kolat
           k=at_type(i)
           do j=1,n_types_central_atoms_3body
              if(k == types(j,1)) then
                 n_tetrahedrons=n_tetrahedrons+1
                 exit
              endif
           enddo
        enddo
        DPRINT  'building_tet: tetrahedrons found', n_tetrahedrons

        allocate(tetra_atoms(n_tetrahedrons,5),stat=status)
        ASSERT(status.eq.0)
        tetra_atoms=0

        i1=0
        lab1:do i=1,epe_nucen+kolat
           k=at_type(i)
           lab2: do j=1,n_types_central_atoms_3body
              if(k == types(j,1)) then
                 exit_cycle=.false.
                 ind=j
                 exit lab2
              else
                 exit_cycle=.true.
              endif
           enddo lab2

           if (exit_cycle) cycle lab1

           i1=i1+1
           tetra_atoms(i1,1)=i

           buf_dist=0.0_r8_kind
           buf_index=0_i4_kind
           lab3: do l=1,epe_nucen+kolat
              if (l==i) cycle lab3

              k1=at_type(l)
              if(k1.eq.0) cycle  lab3
              lab4: do m=2,5
                 if (k1 == types(ind,m)) then
                    exit_cycle=.false.
                    exit lab4
                 else
                    exit_cycle=.true.
                 endif
              enddo lab4

              if (exit_cycle) cycle lab3

              dist=sqrt(dot_product(epe(i)%s-epe(l)%s,epe(i)%s-epe(l)%s))
              if(dist*0.529177_r8_kind > r3b(ind)) cycle lab3

              lab5: do n=1,4
                 if(buf_dist(n) == 0.0_r8_kind) then
                    buf_dist(n)=dist
                    buf_index(n)=l
                    cycle lab3
                 endif
              enddo lab5
              max_ind=maxloc(buf_dist)
              max_val=maxval(buf_dist)
              if(dist < max_val) then
                 buf_dist(max_ind)=dist
                 buf_index(max_ind)=l
              endif
           enddo lab3
           tetra_atoms(i1,2:5)=buf_index
           DPRINT  'tetra_atoms',i1, tetra_atoms(i1,1:5)
        enddo lab1
      end subroutine building_tet
#endif

  subroutine convergence_check(converged)
    ! Purpose: performs the various convergence checks and,
    !          if convergence has not yet been reached, the
    !          'conv'-file is deleted.
    ! ----------------------------------------------------
    logical, intent(out):: converged
    logical :: conv
    conv=.true.

    if(tsscan_sphere) then
     conv=grad_max_sphere.lt.max_comp_grad.and.conv
     write(io_flepo,*) grad_max_sphere.lt.max_comp_grad, 'grad_max_sphere conv'
     conv=abs(dEdR_sphere).lt.max_dEdR_sphere.and.conv
     write(io_flepo,*) abs(dEdR_sphere).lt.max_dEdR_sphere, 'max_dEdR_sphere conv', &
     abs(dEdR_sphere), max_dEdR_sphere
    else
    if (grad_max_comp>max_comp_grad) then
       conv=.false.
    endif
    endif

    if(tsscan_sphere) then
    conv=grad_mean_sphere.lt.rms_grad.and.conv
     write(io_flepo,*) grad_mean_sphere.lt.rms_grad, 'grad_mean_sphere conv'
    else
    if (grad_mean_square>rms_grad) then
       conv=.false.
    endif
    endif

    if (step_mean_square>rms_step) then
       conv=.false.
    endif
    if (step_max_comp>max_comp_step) then
       conv=.false.
    endif


    if (conv) then
       converged = .true.
       write(OPT_STDOUT,*)" main_opt: the optimization is converged"
    else
       write(OPT_STDOUT,*)" main_opt: optimization not yet converged"
    endif
  end subroutine convergence_check

  function k_parameter(r0,i,j,k,type) result(k_par)
    !
    ! Purpose : returns the force parameter in accordance with model
    !           stretching-bending potential
    ! Reference.: M.V.Fernandez_Serra, E.Artacho, J.N.Soler
    !             Phys.Rev.B.67,100101(R) (2003)
    !
    !
    !------- Modules used -----------------------------------------
    !------- Declaration of formal parameters ---------------------

    use coortype_module
    use atom_data_module
    use slspar_module

    real(kind=r8_kind)::             k_par
    integer(kind=i4_kind), optional, intent(in)    :: i,j,k,type
    real(kind=r8_kind), intent(out):: r0
    !------- Declaration of local variables -----------------------
    real(kind=r8_kind), parameter :: A = 308.722_r8_kind*evau, &
                                     B = 0.05_r8_kind
    real(kind=r8_kind)            :: r_ij,r_jk, R_i,R_j,R_k,angcos
    !------- executable code --------------------------------------

    K_par=0.0_r8_kind
    r0=0.0_r8_kind
    select case(type)
       case(b_length)
          r_ij  = sqrt( dot_product( gx(i)%eq-gx(j)%eq, gx(i)%eq-gx(j)%eq) )
          if(r_ij*angsau.gt.19.9) return
!          R_i   = covrad(nint(epe(i)%ant))/angsau
!          R_j   = covrad(nint(epe(j)%ant))/angsau
          R_i   = covrad(nint(charge(i)))/angsau
          R_j   = covrad(nint(charge(j)))/angsau
          K_par = A*((R_i+R_j)/r_ij)**8 / 10.0_r8_kind
          r0=r_ij
#if 0
          write(*,*) "K_parameter : i,j    ",i,j
          write(*,*) "K_parameter : cov    ",covrad(nint(epe(i)%ant)),covrad(nint(epe(j)%ant))
          write(*,*) "K_parameter : atom i ", gx(i)%eq
          write(*,*) "K_parameter : atom j ", gx(j)%eq
          write(*,*) "K_parameter : r_ij ", r_ij*angsau
          write(*,*) "K_parameter : K    ", K_par
#endif
       case(b_length_toepe)
          r_ij  = sqrt( dot_product( gx(i)%r-epe(j)%s, gx(i)%r-epe(j)%s) )
          R_i   = covrad(nint(epe(i)%ant))/angsau
          R_j   = covrad(nint(epe(j)%ant))/angsau
          K_par = A*((R_i+R_j)/r_ij)**8 / 10.0_r8_kind
          r0=r_ij
         if(r_ij.lt.5.66) then
          write(*,*) "K_parameter : i,j    ",i,j
          write(*,*) "K_parameter : cov    ",covrad(nint(epe(i)%ant)),covrad(nint(epe(j)%ant))
          write(*,*) "K_parameter : atom i ", gx(i)%r
          write(*,*) "K_parameter : atom j ", epe(j)%s
          write(*,*) "K_parameter : r_ij ", r_ij*0.529
          write(*,*) "K_parameter : K    ", K_par
         endif
       case(b_angle)
          r_ij  = sqrt( dot_product(gx(i)%eq-gx(j)%eq,gx(i)%eq-gx(j)%eq) )
          if(r_ij*angsau.gt.3.1_r8_kind) return
          r_jk  = sqrt( dot_product(gx(j)%eq-gx(k)%eq,gx(j)%eq-gx(k)%eq) )
          if(r_jk*angsau.gt.3.1_r8_kind) return
!          R_i   = covrad(nint(epe(i)%ant))/angsau
!          R_j   = covrad(nint(epe(j)%ant))/angsau
!          R_k   = covrad(nint(epe(k)%ant))/angsau
          R_i   = covrad(nint(charge(i)))/angsau
          R_j   = covrad(nint(charge(j)))/angsau
          R_k   = covrad(nint(charge(k)))/angsau
          angcos=dot_product(gx(j)%eq-gx(i)%eq,gx(j)%eq-gx(k)%eq)/r_ij/r_jk
          if(angcos.gt.0.9999) angcos=0.9999
          if(angcos.lt.-0.9999) angcos=-0.9999
          if(abs(angcos).gt.0.96) return
          K_par = B*A * ((R_i+R_j)/r_ij)**4 * ((R_k+R_j)/r_jk)**4 * r_ij*r_jk  / 10.0_r8_kind
          r0=acos(angcos)
#if 0
          write(*,*) "K_parameter : i,j    ",i,j,k
          write(*,*) "K_parameter : atom i ", gx(i)%eq
          write(*,*) "K_parameter : atom j ", gx(j)%eq
          write(*,*) "K_parameter : atom k ", gx(k)%eq
          write(*,*) "K_parameter : K   angle ", K_par,r0
#endif
    end select

  end function k_parameter
end subroutine main_opt

  subroutine persistent_state(act)
    !
    ! optimizer seems to be using some SAVEed variables to store
    ! persistent data over geometry iterations. This state is lost
    ! after crash/restart cycle.
    !
    ! PLEASE DONT (AB)USE PERSISTENT STATE IF POSSIBLE!
    !
    use opt_data_module, only: OPT_STDOUT
    use iounitadmin_module, only: openget_iounit, returnclose_iounit
    use filename_module, only: inpfile, MAXPATH=>filename_namelengthmax

    ! these are the persistant varibales, KEEP THEM TO MINIMUM!:
    use hesse_module, only: step_counter
    use step_module, only: trust_radius => r_curr &
                         , total_energy => e_prev

    ! these are the subroutines that given action and iounit do additional work:
    use step_module, only: step_module_persistent_state
    use ts_module, only: ts_module_persistent_state
    implicit none
    character(len=*), intent(in) :: act ! restore or save
    ! *** end of interface ***

    namelist /optimizer_state/ step_counter &
                             , trust_radius &
                             , total_energy

    character(len=MAXPATH) :: state_file
    integer(i4_kind)       :: iounit,stat

    ! for compatibility reasons dont save/restore persistent state on disk:
    logical, save          :: exists = .false.

    state_file = trim(inpfile('optimizer.state'))

    ! if "optimizer.state" does exist => user requested persistency:
    inquire(exist=exists,file=trim(state_file))

    select case(act)

    case('restore')

      if( .not.exists )then
        write(*,*)" main_opt: no saved state found, relying on initial values ..."
        goto 999 ! finalize and exit (so far only dump the state onto tty)
      endif

      write(OPT_STDOUT,*)" main_opt: found saved state, setting persistent variables ..."

      iounit = openget_iounit(status='old',form='formatted',file=trim(state_file))

      ! this will fail with IOSTAT/=0 if the file is empty:
      read(iounit,nml=optimizer_state,IOSTAT=stat)

      if( stat == -1 )then
        ! an immediate EOF => user "touched" an empty file to initiate persistency
        write(OPT_STDOUT,*)" main_opt: saved state is empty, relying on initial values ..."
      else
        ! continue reading ...
        call step_module_persistent_state('restore',iounit)
        call ts_module_persistent_state('restore',iounit)
      endif

      call returnclose_iounit(iounit)

    case('save')

      if( .not.exists )then
#ifndef FPP_OPTIMIZER
        ! ignore saving state if "optimizer.state" does not exist
        write(*,*)" main_opt: no saved state found, not saving either"
        goto 999 ! finalize and exit (so far only dump the state onto tty)
#else
        ! if optimizer is compiled as a standalone prog it makes sense
        ! to create an empty "optimizer.state" by default:
        write(*,*)" main_opt: no saved state found, creating an empty file"
        iounit = openget_iounit(status='new',form='formatted',file=trim(state_file))
        call returnclose_iounit(iounit)
#endif
      endif

      write(OPT_STDOUT,*)" main_opt: overwriting on-disk persistent state ..."

      iounit = openget_iounit(status='old',form='formatted',file=trim(state_file))

      write(iounit,nml=optimizer_state)

      call step_module_persistent_state('save',iounit)
      call ts_module_persistent_state('save',iounit)

      call returnclose_iounit(iounit)
    case default
      ABORT('no such action')
    end select

999 CONTINUE ! clean up and exit
  end subroutine persistent_state


  !--------------- End of module -------------------------------------
end module optimizer

