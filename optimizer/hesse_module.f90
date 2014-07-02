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
module  hesse_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Contains Hesse-Matrix  (in internal coordinates)
  !           Routines for input and output of the Hesse-Matrix,
  !           Transformation Routines ...
  !
  !  Module called by: main_opt ...
  !
  !  References: ...
  !  Author: FN
  !  Date: 6/97
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: AG
  ! Date:   9/03
  ! Description: CG-method
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  !   Contents of the file 'HESSE.DAT':
  ! - 1 integer (kind=i4_kind) : this is interpreted as the actual
  !   cycle number belong to the data on the file in all places
  !   except for oone. This exception is the exact calculation of the
  !   Hessian. Here, it is interpreted as the symmetry type (->'sym_type')
  !   if the actual internal variable being elongated.
  !
  ! - The coordinate vector q. Note that this coordinate vector is interpreted
  !   as follows:
  !   a) ZMAT_COORDINATES:
  !   Here, it is the coordinate q (already symmetry reduced) as it enters the
  !   optimization procedure.
  !   b) DELOCALIZED COORDINATES:
  !   In this case, the coordinate vector is stored in the PRIMITIVE coordinate
  !   system. The reason is, that when calculating differences between two geometries,
  !   the dihedral angles need special consideration. This can be done ONLY in primitive
  !   variables. Thus the routine 'delta_coordinate' will first calculated
  !   the difference in primitive coordinates, and then transform the difference to
  !   delocalized coordinates.
  !
  ! - The gradient vector: this is always stored within the actual set of coordinates
  !
  ! - OPTIONALLY: for a exact calculation of the Hessian with a two-step finite
  !   difference algorithms, also the coordinate vector and the gradient vector
  !   from the first step are stored. Since an exact calculation of the Hessian
  !   is possible ONLY in ZMAT_COORDINATES, these two items are always interpreted
  !   within the ZMAT_COORDINATE system.
  !
  ! - The Hessian itself in quadratic storage mode (could be improved) within
  !   the actual set of coordinates
  !
  ! - (AG) while CG-method does not require Hessian calculation, it is implemented
  !   using Hesse-matrix formally - this will not require more memory than Hesse-methods,
  !   but the program structure can be kept almost unchanged and reused in case of
  !   CG-method as well;
  !--------------------------------------------------------------------------------
  !
#include <def.h>
  use type_module ! type specification parameters
  use coortype_module, only: int_coor
  use iounitadmin_module
  use allocopt_module
#ifdef WITH_EFP
  use qmmm_interface_module, only: efp
  use efp_module, only: n_efp, efp_fixed, qm_fixed
#endif

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of constants and variables ---------------
  integer(kind=i4_kind), public :: step_counter=0
  integer(kind=i4_kind), public :: mep_point=0

  real(kind=r8_kind),allocatable,public    :: hesse(:,:),hesse_inv(:,:),pmat(:,:)
  real(kind=r8_kind),allocatable,public    :: hesse_eigvec(:,:),&
                                              hesse_eigval(:),hesse_diaval(:,:)
  real(kind=r8_kind),allocatable,public    :: hesse_sphere(:,:)

  logical, public :: reinit_hess=.false.
  !------------ public functions and subroutines ---------------------
  public hesse_main,hesse_update,hesse_write,hesse_eval,rewrite_hesse
  public hesse_shutdown, hess_and_grad_project,write_old_hess
  public cart_step

  !===================================================================
  ! End of public interface of module
  !===================================================================
  real(kind=r8_kind),parameter    :: delta_q=0.01_r8_kind

  integer(kind=i4_kind),parameter          :: unit_matrix=1,&
                                              rough_estimate=2,&
                                              force_field=3
  integer(kind=i4_kind)                    :: io_hesse, io_hesse_t
  real(kind=r8_kind),allocatable           :: q_last(:),g_last(:)
  real(kind=r8_kind),allocatable           :: dq(:),dg(:)
  type(int_coor),allocatable               :: s_local(:)
  integer(kind=i4_kind)                    :: alloc_sta(9),last_loop
  logical                                  :: do_projection

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine hesse_shutdown()
    implicit none
    ! *** end of interface ***

    if(allocated(hesse_eigval)) then
      deallocate(hesse_eigval,hesse_eigvec,stat=allocopt_stat(6))
      ASSERT(allocopt_stat(6).eq.0)
      allocopt_stat(6)=1
    endif

    if(allocated(hesse_diaval)) then
      deallocate(hesse_diaval,stat=allocopt_stat(7))
      ASSERT(allocopt_stat(7).eq.0)
      allocopt_stat(7)=1
    endif
  end subroutine hesse_shutdown

  !*************************************************************

  subroutine cart_step(geo_loop,converged)
    ! Calculate numerical hessian
    ! To activate:
    ! PG input - operations_geo_opt = t
    !            max_geo_iterations - 3*N_atoms+1
    ! optimizer input - within namelist "operations"
    !           calc_cart_hess = t
    use filename_module, only: inpfile
    use gradient_module, only: grad_cartes
    use opt_data_module, only: n_atoms,n_dummy,x,y,z,dummy_list,step_size
    use coordinates_module
    implicit none
    integer(i4_kind),intent(in)   :: geo_loop
    logical         ,intent(out)  :: converged
    ! *** end of interface ***

    real(r8_kind) :: cxyz(3*n_atoms),gxyz(3*n_atoms),gxyz1(3*n_atoms)
    real(r8_kind) :: hess(3*n_atoms,3*n_atoms)
    real(r8_kind) :: step
    integer(i4_kind) :: i_step
    real(r8_kind) :: sstep,dstep
    integer(i4_kind) :: i,j,k
    integer(i4_kind) :: igrads,loop
    logical :: yes

    DPRINT  'cart_step: entered, geo_loop=',geo_loop
    sstep=step_size
    dstep=2.0_r8_kind*sstep
    converged=.false.

    igrads=get_iounit()
    if(geo_loop > 1) then
       inquire(file=trim(inpfile('grads')), exist=yes)
       if(.not.yes) call error_handler('HESSE_MODULE: cart_step: file "grads" does not exist ')
    end if
    open(igrads,file=trim(inpfile('grads')), err=100)

    if(geo_loop==1) then
       ! take current geometry from optimizer globals:
       k=0
       do i=1,n_atoms+n_dummy
          if(.not.dummy_list(i)) then
             k=k+1; cxyz(k)=x(i)
             k=k+1; cxyz(k)=y(i)
             k=k+1; cxyz(k)=z(i)
          end if
       end do
       ! save central geometry to a file:
       write(igrads,*) cxyz(:)
    else
       ! read in the central geometry from the file:
       read(igrads,*) cxyz(:)

       do i=1,geo_loop-2
          read(igrads,*) loop,gxyz(:)
       end do
       k=0
       do i=1,n_atoms+n_dummy
          do j=1,3
             if(.not.dummy_list(i)) then
                k=k+1
                gxyz(k)=grad_cartes(i,j)
             end if
          end do
       end do
       write(igrads,*) geo_loop,gxyz(:)

       if(geo_loop == 6*n_atoms+1) then
          converged=.true.

          rewind igrads
          read(igrads,*) cxyz(:)

          do i=1,3*n_atoms
             read(igrads,*) loop,gxyz(:)
             read(igrads,*) loop,gxyz1(:)
             hess(:,i)=(gxyz-gxyz1)/dstep
          end do

          ! dump cartesian hessian into a file and compute frequencies:
          call write_cart_hess(hess)
       end if
    end if

    call returnclose_iounit(igrads)

    DPRINT  'cart_step: exit, converged=',converged
    if( .not.converged )then
      if (mod(geo_loop,2) /= 0) then
         step=sstep
         i_step=int(geo_loop/2)+1
      else
         step=-sstep
         i_step=int(geo_loop/2)
      end if

      ! make a step into the new direction:
      cxyz(i_step)=cxyz(i_step)+step
      DPRINT  'cart_step: step by',step,'in direction',i_step

      ! feed new geometry to optimizer globals:
      k=0
      do i=1,n_atoms+n_dummy
         if(.not.dummy_list(i)) then
            k=k+1; x(i)=cxyz(k)
            k=k+1; y(i)=cxyz(k)
            k=k+1; z(i)=cxyz(k)
         end if
      end do
    endif
    return
100 call error_handler('HESSE_MODULE: cart_step: error - file "grads"')
  end subroutine cart_step

  !*************************************************************

  subroutine hesse_main(geo_loop,hesse_complete)
    ! Purpose: allocate hesse and hesse_inv. If geo_loop=1,
    !          initialize hesse and hesse_inv. Write this
    !          directly to the hesse.dat together with coordinates q
    !          and gradient grad_intern.
    !          If geo_loop is not 1 then read
    !          in from file:
    !          q_last,g_last,hesse_last
    !          and update it to hesse_act taking q_act and
    !          g_act from gxfile.
    !
    ! calculation of hessian:
    ! a) single_step: The hessian is calculated as
    !    g0-grad_intern/(q0-q). The first call to 'optimizer'
    !    will calculate q0 and g0, allocate and intialize the hessian
    !    and write them to file. A step for the first coordinate is
    !    taken. The next calls will read q0,g0 and hesse from file,
    !    calculate the hessian with the actual values grad_intern and q
    !    and write q0,g0 and hesse to file again.
    ! b) .not.single_step: The hessian is calculated as
    !     1/2 ( g0-g1/q0-q1  + grad_intern-g0/q-q0 ).
    !     The first call to 'optimizer' calculates q0 and g0, allocates
    !     hesse, q1 and g1 , initializes them and writes them to file.
    !     The next steps have the following possibilities:
    !     (I) after reading 'hesse.dat', q1 and g1 are found to be zero.
    !         This means that the actual step is actually q1 and
    !         grad_intern=g1. Set those variables correspondingly,
    !         make the second step for the same coordinate and
    !         write q0,g0,q1,g1 and hesse to file.
    !    (II) after reading 'hesse.dat' q1 and g1 are found to be different
    !         from zero. This means that the actual step corresponds to q2
    !         (and g2). Given q0,g0,q1,g1,q,grad_intern, calculate the
    !         appropriate matrix elements of hesse. Set q1 and g1 to zero,
    !         make the first step for the next variable and write
    !         q0,g0,q1,g1 and hesse to file.
    !
    ! calculation of hessian for delocalized internals:
    !  This is possible ONLY if the coordinates are specified in zmat_format
    !  in order to be able to save the symmetry equivalent steps.
    !  Thus the calculation of the hessian runs along the same lines
    !  as described above, the only difference being that as sson as
    !  the hessian is complete, it is transformed back to delocalized
    !  coordinates. This means that as long as the calculation is still
    !  running, the data on 'hesse.dat' are interpreted as primitive
    !  internal coordinates. It is for this reason that the variable
    !  s_local is used here.
    !
    ! start calculation with already existing
    ! Hessian in 'hesse.dat':
    ! a) usual minimization without initial calculation of hessian
    !    i.e. mo ts_search-mode:
    !    If 'hesse.dat' is present, read in this file. If geo_loop==1
    !    throw away the information g_last and q_last and do NOT
    !    update the Hessian.
    !
    ! b) frequency calculation or ts-mode calculation (minimization or
    !    transition state search):
    !    If geo_loop==1 the Hessian is read in from file and the switch
    !    hesse_complete is set to true. This also requires that
    !    calc_hessian is set to FALSE in order to keep the program
    !    from trying to go on with the calculation of the Hessian
    !    in cycles above 1 bit below n_hesse_step. Here also the information
    !    about last point and gradient are discarded if geo_loop==1.
    !    If 1 < geo_loop < n_hesse_step:  This is interpreted as a continuation
    !    if a Hessian calculation.
    !
    ! Re-calculation of hessian every 'n_ts_step' steps:
    ! This mean ins practive that every n_ts_step+n_hesse_step
    ! the calculation has to be started again. In this case the
    ! existence of the file 'hesse.dat' has to be ignored.
    ! Strategy for handling the existence of 'hesse.dat' and the
    ! switch 'calc_hessian':
    !     'hesse.dat' exists    calc_hessian
    !             T                  T         file will be overwritten
    !             T                  F         contents will be used
    !             F                  T         file will be created
    !             F                  F         usual initial hesse matrix
    !------------ Modules used ---------------------------------
    use coortype_module
    use opt_data_module
    use coordinates_module
    use math_module
    use gradient_module, only: grad_intern
    use valence_coord_module
    use filename_module, only: inpfile
    !------------ Declaration of formal parameters ---------------
    implicit none
    integer(kind=i4_kind),intent(in)   :: geo_loop
    logical,              intent(out)  :: hesse_complete

    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind),allocatable  :: q0(:),g0(:),q1(:),g1(:),&
                                       g_actual(:),q_actual(:),&
                                       help_mat(:,:)
    integer(kind=i4_kind)           :: i,j,alloc_stat,sym_last,sym_act,&
                                       n_var,n_hesse_step,h_calc
    logical                         :: have_hesse_dat,err,back
    logical :: not_freq_only
    logical :: prepare_hesse
    !------------ Executable code ------------------------------------

    ! first find out how many cycles the calculation of the
    ! hessian takes
    if(cart_coordinates) then
       do_projection=.true.
#ifdef WITH_EFP
       if(efp_fixed .or. qm_fixed) do_projection=.false.
#endif
       if(epe_interfaced_mode) do_projection=.false.
    end if

    DPRINT 'hesse_main: single_step=',single_step

    call hesse_alloc() ! 'hesse' and 'hesse_inv' will be allocated if not
                       !  already allocated

    not_freq_only=.not.(.not.calc_hessian.and.frequency_calculation)

    DPRINT 'hesse_main: not_freq_only=',not_freq_only
   not_freq_onl: if(not_freq_only) then
    allocate(s_local(n_internal),q_actual(n_internal),&
             g_actual(n_internal),STAT=alloc_sta(1))
    ASSERT(alloc_sta(1).eq.0)
    alloc_sta(1)=1

    DPRINT 'hesse_main: var allocated',shape(s),shape(q),shape(grad_intern)


    if (zmat_coordinates) then
       s_local = s
       q_actual = q
       if (.not.gx_test) then
          g_actual= grad_intern
       else
          g_actual = zero
       endif

    elseif (delocalized_coordinates.and.zmat_format) then
       ! This option will ONLY work if the primitive internals
       ! have been specified in zmat_format. Thus n_internal=n_primitive
       s_local = s_prim
       q_actual= s_prim(:)%value
       g_actual = matmul(umat,grad_intern)

    elseif (delocalized_coordinates.and..not.zmat_format) then
       s_local = s
       q_actual = q
       g_actual = grad_intern
    else if(cart_coordinates) then
       q_actual = q
       g_actual= grad_intern
    endif

    DPRINT 'hesse_main: now check if hesse.dat exist'

    err=.false.
    inquire(EXIST=have_hesse_dat,FILE=trim(inpfile('hesse.dat'))) ! if so 'optimizer'
    !                                                       ! will try to read the
    !                                                       ! Hessian from this file
    if (.not.have_hesse_dat.and.(calc_hessian.and.geo_loop>1)) then
       write(OPT_STDOUT,*)" hesse_main: no file hesse.dat found although &
            &cycle number indicates that calculation of Hessian &
            &is in progress"
       ABORT("no hesse.dat see flepo")
    endif
    DPRINT 'hesse_main: exist hesse.dat',have_hesse_dat

    calchess: if( calc_hessian )then
      n_var = maxval(sym_type(:))

      if (single_step) then
         n_hesse_step = n_var + 1
      else
         n_hesse_step = 2*n_var + 1
      endif

      h_calc=mod(geo_loop,(n_hesse_step+n_recalc-1))
      if(h_calc==0)  h_calc=n_hesse_step

    DPRINT 'call calculate_hessian now(calc_hessian h_calc n_hesse_step'
    DPRINT calc_hessian,h_calc,n_hesse_step

    if (( h_calc>=1 .and. h_calc <=n_hesse_step)) then
       if (.not.have_hesse_dat .and. geo_loop>1) call error_handler("this should not occur")

       write(OPT_STDOUT,*)"                ----------------               "
       write(OPT_STDOUT,*)"   Calculation of Hessian for cycle ",geo_loop
       write(OPT_STDOUT,*)"                ----------------               "
       write(OPT_STDOUT,*)"  "

       call calculate_hessian()

       if (hesse_complete) then

          if (delocalized_coordinates) then
             call print_matrix(hesse,n_internal,n_internal,10_i4_kind)
             allocate(help_mat(n_internal,n_internal),STAT=alloc_stat)
             if (alloc_stat/=0) call error_handler&
                  (" hesse_main: allocation of help_mat failed")
             help_mat = matmul(hesse,umat)
             hesse=matmul(umat_trans,help_mat)
             deallocate(help_mat,STAT=alloc_stat)
             if (alloc_stat/=0) call error_handler&
                  ("hesse_main: deallocation of help_mat failed")
          endif
          write(OPT_STDOUT,*)"                       ----------------                          "
          write(OPT_STDOUT,*)"  Hesse_Main: calculation of hessian finished in cycle ",geo_loop
          write(OPT_STDOUT,*)"                       ----------------                          "
          write(OPT_STDOUT,*)" "
          if (print_hesse) then
             call print_matrix(hesse,n_internal,n_internal,10_i4_kind)
             write(OPT_STDOUT,*)" ---------------------------------- initial hessian ----"
          endif

             if(tsscan_sphere) then
                if(exist_product) then
                 call hessian_rp_model(hesse_sphere,hesse,sphere_dependent_var)
                else
                 call hessian_sphere(hesse_sphere,hesse,sphere_dependent_var)
                endif
                call eigensolve_hessian(hesse_sphere,n_internal-1)
             else
             call eigensolve_hessian(hesse,n_internal) !(1) numerically calculated
             endif

          call hesse_internal_to_cart(int_flag=geo_loop)

       endif! hesse_complete==.true.
       call dea_main
      return
    endif ! if(h_calc>=1 .and. h_calc <=n_hesse_step)
    endif calchess !if (calc_hessian)
 endif not_freq_onl

    DPRINT 'hesse_main:              have_hesse_dat=',have_hesse_dat
    DPRINT 'hesse_main: analitic_hessian_calculated=',analitic_hessian_calculated
    DPRINT 'hesse_main:       update_fromcartessian=',update_fromcartessian

    prepare_hesse = (.not.have_hesse_dat) &
                     .or. analitic_hessian_calculated &
                     .or. update_fromcartessian
                     ! Absoft does not evaluate this inside of ``if'' statement properly

    DPRINT 'hesse_main: prepare_hesse=',prepare_hesse

    if ( prepare_hesse ) then  ! hesse.dat has to be initialized

     DPRINT 'hesse_main: Hessian has to be initialized, analytic=',analitic_hessian_calculated

       if (geo_loop == 1 .or. analitic_hessian_calculated.or.update_fromcartessian ) then
          err=.true.
          ! This routine checks, if a file 'hesse_cartesian.dat' is present
          ! and converts the content to internal coordinates,
          ! otherwise returns err=.true.:
          call hesse_prepare_internal(hesse, err)
          DPRINT 'hesse_main: hesse_prepare_internal() returned err=',err

          if (err) then
             write(OPT_STDOUT,*)" hesse_main: reading of an initial Hessian from file "
             write(OPT_STDOUT,*)"             HESSE_CARTESIAN.DAT failed. Proceed with "
             write(OPT_STDOUT,*)"             initialization of Hessian"
             write(OPT_STDOUT,*)" "
          endif

          if (.not.err ) then
             write(OPT_STDOUT,*)" hesse_main: reading of initial Hessian from file"
             write(OPT_STDOUT,*)"             HESSE_CARTESIAN.DAT successful."
             write(OPT_STDOUT,*)" "

             hesse_inv = hesse
             if(cart_coordinates) then
                if(do_projection) then
                   call hess_project()
                   hesse=matmul(matmul(pmat,hesse),pmat)
                end if
             else
                call invert_matrix(n_internal,hesse_inv)
             end if

             if (print_hesse) then
                write(OPT_STDOUT,*)" ---- Hessian from HESSE_CART_TO_INTERNAL -----------------"
                call print_matrix(hesse,n_internal,n_internal,11_i4_kind)
                write(OPT_STDOUT,*)" ----------------------------------------------------------"
                write(OPT_STDOUT,*)
             endif

             if(not_freq_only) then
                 ! write the Hessian to file HESSE.DAT
                 if (zmat_coordinates) then
                    call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
                         dimen=n_internal,int_flag=geo_loop,temp_hesse=.false.)
                    call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
                         dimen=n_internal,int_flag=geo_loop,store=.true.,temp_hesse=.false.)
                 elseif (delocalized_coordinates ) then
                    call hesse_write(q_act=q_prim,g_act=grad_intern,h_mat=hesse,&
                         dimen=n_internal,int_flag=geo_loop,deloc=.true., temp_hesse=.false.)
                 elseif (cart_coordinates) then
                    call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
                         dimen=n_internal,int_flag=geo_loop,temp_hesse=.false.)
                    call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
                         dimen=n_internal,int_flag=geo_loop,store=.true.,temp_hesse=.false.)
                 endif
                 if(.not.cart_coordinates) then
                    if(tsscan_sphere) then
                    if(exist_product) then
                     call hessian_rp_model(hesse_sphere,hesse,sphere_dependent_var)
                    else
                     call hessian_sphere(hesse_sphere,hesse,sphere_dependent_var)
                    endif
                       call eigensolve_hessian(hesse_sphere,n_internal-1)
                    else
                       call eigensolve_hessian(hesse,n_internal) !(2) analitic or initialized and updated
                    endif
                 end if
                 call dea_main
             endif

             hesse_complete = .true.
             return !!! RETURN POINT !!!
           endif

       else ! geo_loop /= 1 and 'hesse.dat' does not exist
!!$          write(OPT_STDOUT,*)" hesse_main: The cycle number in your input gxfile  is ",geo_loop
!!$          write(OPT_STDOUT,*)"             but the file HESSE.DAT does not exist."
!!$          write(OPT_STDOUT,*)"             Please either "
!!$          write(OPT_STDOUT,*)"             - copy the desired HESSE.DAT in your input directory "
!!$          write(OPT_STDOUT,*)"             or "
!!$          write(OPT_STDOUT,*)"             - set the cycle number to -1.0 "
!!$          call error_handler(" ")
          err=.true.
       endif

       if (err) then
          if (estimate) then
             if (zmat_coordinates) then
                write(OPT_STDOUT,*)" hesse_main : start with initial Hessian containing the following"
                write(OPT_STDOUT,*)"              estimates:"
                write(OPT_STDOUT,*)"              bond-stretches   0.5 Hartree/au**2"
                write(OPT_STDOUT,*)"              angle bends      0.2 Hartree/rad**2"
                write(OPT_STDOUT,*)"              torsions         0.1 Hartree/rad**2"
             elseif (delocalized_coordinates) then
                write(OPT_STDOUT,*)" hesse_main: start with initial Hessian containing the "
                write(OPT_STDOUT,*)"             following value as diagonal elements :"
                write(OPT_STDOUT,*)"             all delocalized coordinates     0.4 Hartree/au**2"
             endif
             call hesse_init(rough_estimate)
          elseif (estimate_hessian) then
             if (.not.(delocalized_coordinates.and.valence_format))&
                  call error_handler&
                  (" hesse_main: control is lost here")
             write(OPT_STDOUT,*)" hesse_main: set up the initial Hessian with empirical "
             write(OPT_STDOUT,*)"             force constants based on valence-forcefield"
             write(OPT_STDOUT,*)"             approximation for the valence coordinates  "
             write(OPT_STDOUT,*)"             underlying the delocalized coordinates"
             call hesse_init(force_field)
          else
             write(OPT_STDOUT,*)" hesse_main: start with a unit Hessian"
             call hesse_init(unit_matrix)
          endif
       endif

       DPRINT 'main_hesse: call invert_matrix(n_internal,hesse_inv)'
       hesse_inv = hesse
       if(cart_coordinates) then
          call break_degen()
          if(do_projection) then
             call hess_project()
!!$print*,'Init HESSE before P'
!!$do i=1,n_primitive
!!$print*,hesse(i,:)
!!$end do
             hesse=matmul(matmul(pmat,hesse),pmat)
!!$print*,'init HESSE after P'
!!$do i=1,n_primitive
!!$print*,hesse(i,:)
!!$end do
          end if
       else
          call invert_matrix(n_internal,hesse_inv)
       end if
       DPRINT 'done'

       ! transform back to cartesians
       if (zmat_coordinates) then
          call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
               dimen=n_internal,int_flag=geo_loop,temp_hesse=.false.)
          call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
               dimen=n_internal,int_flag=geo_loop,temp_hesse=.false., store=.true.)
       elseif (delocalized_coordinates) then
          call hesse_write(q_act=q_prim,g_act=grad_intern,h_mat=hesse, &
               dimen=n_internal,int_flag=geo_loop,deloc=.true.,temp_hesse=.false.)
       elseif (cart_coordinates) then
          call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
               dimen=n_internal,int_flag=geo_loop,temp_hesse=.false.)
          call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
               dimen=n_internal,int_flag=geo_loop,store=.true.,temp_hesse=.false.)
       endif
       hesse_complete=.true.
       if(.not.cart_coordinates) then
          if(tsscan_sphere) then
                if(exist_product) then
                 call hessian_rp_model(hesse_sphere,hesse,sphere_dependent_var)
                else
                 call hessian_sphere(hesse_sphere,hesse,sphere_dependent_var)
                endif
             call eigensolve_hessian(hesse_sphere,n_internal-1)
          else
             call eigensolve_hessian(hesse,n_internal) !(3)
          endif
       end if
!Added symmetrization of hessian
    elseif ( have_hesse_dat ) then

       DPRINT 'analitic_hessian_calculated',analitic_hessian_calculated
       if(.not.(analitic_hessian_calculated.or.update_fromcartessian )) then
       write(OPT_STDOUT,*)"hesse_main:  Proceed with update procedure of the Hessian ... "
        call update_hessian()
       endif

!               .or.update_fromcartessian    &
!               .or.step_counter.eq.0))  call update_hessian()
!       write(OPT_STDOUT,*)"             ... ok"
       if(.not.cart_coordinates) then
          if(tsscan_sphere) then
                if(exist_product) then
                 call hessian_rp_model(hesse_sphere,hesse,sphere_dependent_var)
                else
                 call hessian_sphere(hesse_sphere,hesse,sphere_dependent_var)
                endif
             call eigensolve_hessian(hesse_sphere,n_internal-1)
          else
             call eigensolve_hessian(hesse,n_internal) !(4)
          endif
       end if

    endif

  call dea_main()

  DPRINT 'main_hessse done'

  contains

    subroutine break_degen()
      ! Purpose: Artificially break degeneracies so that accidentally degenerate
      !          modes are split and therefore step restriction along modes
      !          is well defined
      !------------ Modules used -----------------------------------
      !------------ Declaration of formal parameters ---------------
      !------------ Declaration of local variables -----------------
      integer(i4_kind) :: i
      !------------ Executable code --------------------------------

      do i=1,n_internal
         hesse(i,i)=hesse(i,i)+(i-1)*1.0e-7_r8_kind
      end do

    end subroutine break_degen

    subroutine hess_project()
      ! Purpose: used to calculate Pmat to project out from Hessian linearly independent
      !          set of coordinates. In the cartesian case translation
      !          and rotation modes are project out
      !          Now used only for Cartesian Hessian
      !------------ Modules used -----------------------------------
#ifdef WITH_EFP
      use pointcharge_module, only: rcm
#endif
      !------------ Declaration of formal parameters ---------------
      !------------ Declaration of local variables -----------------
      real(r8_kind) :: center_mass(3),r1(3),r2(3),fx
      real(r8_kind), allocatable :: work_array(:,:)
      real(r8_kind),parameter :: small_hp=1.0e-6_r8_kind,small_hp1=1.0e-8_r8_kind
      integer(i4_kind) :: i,j,k,status,n_atom_frg
#ifdef WITH_EFP
      integer(i4_kind) :: l
#endif
      !------------ Executable code --------------------------------

      !Construct normalized vectors in work_array in the direction
      !of the rotations and translations.

      n_atom_frg=n_atoms
#ifdef WITH_EFP
      if(efp .and. n_efp > 0) n_atom_frg=n_atom_frg+n_efp
#endif
      center_mass=zero
      do i=1,n_atoms
         center_mass=center_mass+atom(i)%x/n_atom_frg
      end do
#ifdef WITH_EFP
      if(efp .and. n_efp > 0) then
         do i=1,n_efp
            center_mass=center_mass+rcm(:,i)/n_atom_frg
         end do
      end if
#endif

      allocate(work_array(n_primitive,6), stat=status)
      !n_primitive==3*n_atoms+6*n_efp
      ASSERT(status==0)
      work_array=zero

      !x,y,z translations
      do i=1,3
         do j=1,n_atoms
            k=3*(j-1)
            work_array(k+i,i)=sqrt(one/n_atom_frg)
         end do
#ifdef WITH_EFP
         if(efp .and. n_efp > 0) then
            l=n_atoms+1
            do j=1,n_efp
               k=3*(l-1)
               work_array(k+i,i)=sqrt(one/n_atom_frg)
               l=l+2
            end do
         end if
#endif
      end do
      !x,y,z rotation
      do i=4,6
         do j=1,n_atoms
            r1=atom(j)%x-center_mass
            select case(i)
            case(4)
               r2(1)= zero
               r2(2)=-r1(3)
               r2(3)= r1(2)
            case(5)
               r2(1)= r1(3)
               r2(2)= zero
               r2(3)=-r1(1)
            case(6)
               r2(1)=-r1(2)
               r2(2)= r1(1)
               r2(3)= zero
            end select
            k=3*(j-1)
            work_array(k+1:k+3,i)=r2
         end do
#ifdef WITH_EFP
         if(efp .and. n_efp > 0) then
            l=n_atoms+1
            do j=1,n_efp
               r1=rcm(:,j)-center_mass
               select case(i)
               case(4)
                  r2(1)= zero
                  r2(2)=-r1(3)
                  r2(3)= r1(2)
               case(5)
                  r2(1)= r1(3)
                  r2(2)= zero
                  r2(3)=-r1(1)
               case(6)
                  r2(1)=-r1(2)
                  r2(2)= r1(1)
                  r2(3)= zero
               end select
               k=3*(l-1)
               work_array(k+1:k+3,i)=r2
               l=l+1

               select case(i)
               case(4)
                  r2(1)=one
                  r2(2)=zero
                  r2(3)=zero
               case(5)
                  r2(1)=zero
                  r2(2)=one
                  r2(3)=zero
               case(6)
                  r2(1)=zero
                  r2(2)=zero
                  r2(3)=one
               end select
               k=3*(l-1)
               work_array(k+1:k+3,i)=r2
               l=l+1
            end do
         end if
#endif
         do j=1,i-1
            fx=dot_product(work_array(:,j),work_array(:,i))
            work_array(:,i)=work_array(:,i)-fx*work_array(:,j)
         end do
         fx=sqrt(dot_product(work_array(:,i),work_array(:,i)))
         if(fx > small_hp) then
            work_array(:,i)=work_array(:,i)/fx
         else
            work_array(:,i)=zero
         end if
      end do
!!$print*,'WORK'
!!$do i=1,n_primitive
!!$print*,work_array(i,:)
!!$end do


      if(allocated(pmat)) then
         deallocate(pmat, stat=status)
         ASSERT(status==0)
      end if
      allocate(pmat(n_primitive,n_primitive), stat=status)
      ASSERT(status==0)
      pmat=zero
      do i=1,n_primitive
         pmat(i,i)=one
      end do
      do i=1,n_primitive
         do j=1,n_primitive
            do k=1,6
               pmat(i,j)=pmat(i,j)-work_array(j,k)*work_array(i,k)
            end do
         end do
      end do
      do i=1,n_primitive
         do j=1,n_primitive
            if(abs(pmat(i,j)) < small_hp1) pmat(i,j)=zero
         end do
      end do

!!$print*,'PMAT',n_primitive
!!$do i=1,n_primitive
!!$print*,pmat(i,:)
!!$end do


      deallocate(work_array, stat=status)
      ASSERT(status==0)
    end subroutine hess_project

    subroutine dea_main()
     if(allocated(s_local)) &
     deallocate(s_local,STAT=alloc_sta(1))
     ASSERT(alloc_sta(1).eq.0)
     if(allocated(q_actual)) &
     deallocate(q_actual,g_actual,STAT=alloc_sta(1))
     ASSERT(alloc_sta(1).eq.0)
   end subroutine dea_main

    subroutine hessian_sphere(hesse_sphere,hesse,kk)

      real(kind=r8_kind),intent(out):: hesse_sphere(:,:)
      real(kind=r8_kind),intent(in):: hesse(:,:)
      integer(kind=i4_kind), intent(in) :: kk

      real(kind=r8_kind):: ssk,ssi,ssj
      real(kind=r8_kind):: dRdT,d2Ed2t
      integer(kind=i4_kind):: i,j,ii,jj

    ssk=s(kk)%value-s_reactant(kk)%value
    ii=0

    do i=1,size(hesse,1)
    if(i.eq.kk) cycle
     ssi=s(i)%value-s_reactant(i)%value
     if(select_sphere_vars.and..not.s(i)%sphere) ssi=zero
     ii=ii+1
     jj=0
     do j=1, size(hesse,1)
     if(j.eq.kk) cycle
     ssj=s(j)%value-s_reactant(j)%value
     if(select_sphere_vars.and..not.s(j)%sphere) ssj=zero
     jj=jj+1
     hesse_sphere(ii,jj)=hesse(i,j)-hesse(i,kk)*ssj/ssk-hesse(kk,j)*ssi/ssk &
          +hesse(kk,kk)*ssj*ssi/ssk**2-grad_intern(kk)*ssj*ssi/ssk**3
     enddo
     if(.not.select_sphere_vars.or.s(i)%sphere) &
     hesse_sphere(ii,ii)=hesse_sphere(ii,ii)-grad_intern(kk)/ssk
    enddo

    dRdT=distance_to_reactant/ssk
    d2Ed2t= grad_intern(kk)*(1-dRdT**2)/ssk+hesse(kk,kk)*dRdT**2
    distance_to_ts=-dRdT*grad_intern(kk)/d2Ed2t

    write(io_flepo,*) 'd2Ed2t distance_to_ts' , d2Ed2t, distance_to_ts

    end subroutine hessian_sphere

    subroutine hessian_rp_model(hesse_rp,hesse,kk)

      real(kind=r8_kind),intent(out):: hesse_rp(:,:)
      real(kind=r8_kind),intent(in):: hesse(:,:)
      integer(kind=i4_kind), intent(in) :: kk

      real(kind=r8_kind):: ssrk,sspk,ssri,sspi,ssrj,sspj,ssk,ssi,ssj,rp_var
      real(kind=r8_kind):: P2,R2,dRdT,d2Ed2t,d2RdT2
      integer(kind=i4_kind):: i,j,ii,jj

    rp_var=tsscan_rp_var
    ssrk=s(kk)%value-s_reactant(kk)%value
    sspk=s(kk)%value-s_product(kk)%value
    ssk=rp_var**2*ssrk-(1-rp_var)**2*sspk
    ii=0

    do i=1,size(hesse,1)
    if(i.eq.kk) cycle
     ssri=s(i)%value-s_reactant(i)%value
     sspi=s(i)%value-s_product(i)%value
     ssi=rp_var**2*ssri-(1-rp_var)**2*sspi
     if(select_sphere_vars.and..not.s(i)%sphere) ssi=zero
     ii=ii+1
     jj=0
     do j=1, size(hesse,1)
     if(j.eq.kk) cycle
     ssrj=s(j)%value-s_reactant(j)%value
     sspj=s(j)%value-s_product(j)%value
     ssj=rp_var**2*ssrj-(1-rp_var)**2*sspj
     if(select_sphere_vars.and..not.s(j)%sphere) ssj=zero
     jj=jj+1
     hesse_rp(ii,jj)=hesse(i,j)-hesse(i,kk)*ssj/ssk-hesse(kk,j)*ssi/ssk &
    +hesse(kk,kk)*ssj*ssi/ssk**2-(2*rp_var-1)*grad_intern(kk)*ssj*ssi/ssk**3
     enddo
     if(.not.select_sphere_vars.or.s(i)%sphere) &
     hesse_rp(ii,ii)=hesse_rp(ii,ii)-grad_intern(kk)*(2*rp_var-1)/ssk
    enddo

    P2=distance_to_product**2
    R2=distance_to_reactant**2
    dRdT=-(P2*(1-rp_var)+R2*rp_var)/ssk
    d2RdT2=( P2-R2 -4*dRdT*(rp_var*ssrk+(1-rp_var)*sspk)-(2*rp_var-1)*dRdT**2 )/ssk
    d2Ed2t= grad_intern(kk)*d2RdT2+hesse(kk,kk)*dRdT**2
    distance_to_ts=-dRdT*grad_intern(kk)/d2Ed2t

    write(io_flepo,*) 'd2Ed2t(g,h) distance_to_ts' , d2Ed2t,grad_intern(kk)*d2RdT2,hesse(kk,kk)*dRdT**2, distance_to_ts
    write(io_flepo,*) grad_intern(kk),'grad_intern(kk)',kk,hesse(kk,kk),'h(k,k)'




    write(io_flepo,*) dRdT, d2RdT2,'dRdT d2R/dt2'




    end subroutine hessian_rp_model

    subroutine calculate_hessian()
      ! Purpose: provide logic for setting up an exact Hessian
      !          in a step-by-step calculation (finite differences).
      !
      !          This can be applied ONLY for ZMAT_COORDINATES !!!
      !
      ! ----------------------------------------------------------
      ! --- declaration of local variables -----------------------
      allocate(g0(n_internal),q0(n_internal),STAT=alloc_stat)
      ASSERT(alloc_stat.eq.0)
      if (.not.single_step) then
         allocate(g1(n_internal),q1(n_internal),STAT=alloc_stat)
         ASSERT(alloc_stat.eq.0)
         g1=zero
         q1=zero
      endif

      if (h_calc==1) then
         write(OPT_STDOUT,*)" calculate_hessian: writing initial hesse matrix to file"
         sym_act=1_i4_kind
         if (single_step) then
            call hesse_write(q_act=q_actual,g_act=g_actual,&
                 h_mat=hesse,dimen=n_internal,int_flag=sym_act,temp_hesse=.false.)
            call hesse_write(q_act=q_actual,g_act=g_actual,&
                 h_mat=hesse,dimen=n_internal,int_flag=sym_act,store=.true.,temp_hesse=.false.)
         else
            call hesse_write(q_actual,g_actual,q1,g1,hesse,n_internal,sym_act,temp_hesse=.false.)
            call hesse_write(q_actual,g_actual,q1,g1,hesse,n_internal,sym_act,temp_hesse=.false.,store=.true.)
         endif
         call hesse_step(q_actual,sym_act,step_size)
         hesse_complete=.false.
      else
         if (single_step) then
            call hesse_read(q_act=q0,g_act=g0,h_mat=hesse,&
                 dimen=n_internal,int_flag=sym_last,temp_hesse=.false.)
         else
            call hesse_read(q0,g0,q1,g1,hesse,n_internal,sym_last,temp_hesse=.false.)
         endif

         ! determine the step
         if (single_step) then
            sym_act=sym_last+1
         else
            if (sum(q1)==zero) then ! we can safely assume that
               sym_act=sym_last     ! step q1 has not yet been taken
               back=.true.
            else
               sym_act=sym_last+1
               back=.false.
            endif
         endif
         if (h_calc<n_hesse_step) then
            if (single_step) then
               call hesse_calc(q0=q0,g0=g0,q_act=q_actual,g_act=g_actual,sym=sym_last)
               call hesse_step(q0,sym_act,step_size)
            else
               if (sum(q1)/=zero ) then
                  call hesse_calc(q0,g0,q1,g1,q_actual,g_actual,sym_last)
                  q1=zero
                  g1=zero
               else
                  q1=q_actual
                  g1=g_actual
               endif
               call hesse_step(q0,sym_act,step_size,back)
            endif
            hesse_complete=.false.
         elseif(h_calc==n_hesse_step) then
            if (single_step) then
               call hesse_calc(q0=q0,g0=g0,q_act=q_actual,g_act=g_actual,sym=sym_last)
            else
               call hesse_calc(q0,g0,q1,g1,q_actual,g_actual,sym_last)
            endif
            DPRINT 'symmetrize the hessian after initial calculation'
            do i=1,n_internal
               if (.not.s_local(i)%var) &
                    hesse(i,i) = one
               do j=1,i
                  hesse(i,j) = half*(hesse(i,j) + hesse(j,i))
                  hesse(j,i) = hesse(i,j)
               enddo
            enddo
            hesse_complete=.true.
            q = q0
         endif

         if (hesse_complete.or.single_step) then
            call hesse_write(q_act=q0,g_act=g0,h_mat=hesse,&
                 dimen=n_internal,int_flag=sym_act,temp_hesse=.false.)
            call hesse_write(q_act=q0,g_act=g0,h_mat=hesse,&
                 dimen=n_internal,int_flag=sym_act,store=.true. ,temp_hesse=.false.)
         else
            call hesse_write(q0,g0,q1,g1,hesse,n_internal,sym_act,temp_hesse=.false.)
            call hesse_write(q0,g0,q1,g1,hesse,n_internal,sym_act,store=.true.,temp_hesse=.false.)
         endif
      endif

      deallocate(q0,g0,STAT=alloc_stat)
      ASSERT(alloc_stat.eq.0)
      if (.not.single_step) then
         deallocate(q1,g1,STAT=alloc_stat)
         ASSERT(alloc_stat.eq.0)
      endif
    end subroutine calculate_hessian

    subroutine update_hessian()
      ! Purpose: contains logic for updating the Hessian
      ! -------------------------------------------------------
      integer(kind=i4_kind)     :: n_dim
      real(kind=r8_kind)        :: gm

      if (geo_loop==1) then
         write(OPT_STDOUT,*)" update hessian: starting with initial hessematrix from file HESSE.DAT "
         write(OPT_STDOUT,*)"             additional information about last point and last gradient"
         write(OPT_STDOUT,*)"             on file will be ignored since geo_loop=1"
      endif

      if ( zmat_coordinates .or. cart_coordinates) then
         n_dim = n_internal
      elseif (delocalized_coordinates ) then
         n_dim = n_primitive
      endif

      allocate(q_last(n_dim),g_last(n_internal),&
           dq(n_internal),dg(n_internal),STAT=alloc_stat)
      if (alloc_stat/=0) call error_handler&
           (" update hessian: allocation (1) failed")

      DPRINT 'update_hessian: hesse_read'
      call hesse_read(q_act=q_last,g_act=g_last,h_mat=hesse,&
           dimen=n_internal,int_flag=last_loop, &
           deloc=delocalized_coordinates,temp_hesse=.false.)

      DPRINT 'update_hessian: invert_matrix'
      hesse_inv = hesse
      if(cart_coordinates) then
         if(do_projection) call hess_project()
      else
         call invert_matrix(n_internal,hesse_inv)
      end if

      if (zmat_coordinates .or. cart_coordinates) then
         call delta_coordinate(q_last,q,dq)
      elseif (delocalized_coordinates ) then
         call delta_coordinate(q_last,q_prim,dq)
      endif
      DPRINT 'delta_coordinate passed'


      dg = grad_intern - g_last
      gm = dot_product(g_last,g_last)

!      if (geo_loop/=1) then
      if (step_counter/=0) then
         call hesse_update(dq,dg,gm)
         if (print_hesse) then
!            write(OPT_STDOUT,*)"     --- Hesse-Matrix after Update ---"
            call print_matrix(hesse,n_internal,n_internal,10_i4_kind)
            write(OPT_STDOUT,*)" "
         endif
      endif

      if(.not.cart_coordinates) then
         deallocate(q_last,g_last,dq,dg,STAT=alloc_stat)
         ASSERT(alloc_stat==0)
      end if

      if (zmat_coordinates) then
         call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
              dimen=n_internal,int_flag=geo_loop,deloc=.false.,temp_hesse=.false.)
         call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
              dimen=n_internal,int_flag=geo_loop,deloc=.false.,store=.true.,temp_hesse=.false.)
      elseif ( delocalized_coordinates ) then
         call hesse_write(q_act=q_prim,g_act=grad_intern,h_mat=hesse, &
              dimen=n_internal,int_flag=geo_loop,deloc=.true., temp_hesse= .true. )
      elseif ( cart_coordinates ) then
         call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
              dimen=n_internal,int_flag=geo_loop,deloc=.false.,temp_hesse=.false.)
         call hesse_write(q_act=q,g_act=grad_intern,h_mat=hesse,&
              dimen=n_internal,int_flag=geo_loop,deloc=.false.,store=.true.,temp_hesse=.false.)
      endif
      DPRINT 'update_hessian hesse_write passed'

      hesse_complete=.true.

    end subroutine update_hessian

  end subroutine hesse_main
  !*************************************************************

  subroutine write_old_hess()
    use opt_data_module, only : n_internal

    integer(i4_kind) :: alloc_stat

    call hesse_write(q_act=q_last,g_act=g_last,h_mat=hesse_inv,&
         dimen=n_internal,int_flag=last_loop,deloc=.false.,temp_hesse=.false.)
    call hesse_write(q_act=q_last,g_act=g_last,h_mat=hesse,&
         dimen=n_internal,int_flag=last_loop,deloc=.false.,store=.true.,temp_hesse=.false.)

    deallocate(q_last,g_last,dq,dg,STAT=alloc_stat)
    ASSERT(alloc_stat==0)

  end subroutine write_old_hess

  !*************************************************************
  subroutine eigensolve_hessian(hesse,n_internal)
    use opt_data_module, only : force_directed, cart_coordinates
    use matrix_eigenval
    USE_DEBUG
    use math_module, only : zero

    integer(kind=i4_kind), intent(in) :: n_internal
    real(kind=r8_kind), intent(in):: hesse(:,:)

    ! do allocations and eigensolving of hessian
    ! ---------------------------------------------------
    logical :: h_eigvec, h_eigval,h_diaval
    integer(i4_kind):: i
    h_eigvec=allocated(hesse_eigvec)
    h_eigval=allocated(hesse_eigval)
    h_diaval=allocated(hesse_diaval)
    if ((.not.h_eigvec).or.(.not.h_eigval)) then
       allocate(hesse_eigvec(n_internal,n_internal),hesse_eigval(n_internal),STAT=allocopt_stat(6))
       ASSERT(allocopt_stat(6).eq.0)
    endif
#if 1
    DPRINT 'force_directed',force_directed
    if(.not.h_diaval.and.force_directed) then
       allocate(hesse_diaval(n_internal,n_internal),STAT=allocopt_stat(7))
       ASSERT(allocopt_stat(7).eq.0)
       hesse_diaval=0.0_r8_kind
       do i=1,n_internal
          hesse_diaval(i,i)=hesse(i,i)
       enddo
    endif
#endif

    hesse_eigvec=zero
    hesse_eigval=zero
!     call eigensolver(hesse,n_internal,hesse_eigval,hesse_eigvec)
!VVP Use Lapack instead eigensolver of Voytjuk.
    if(.not.cart_coordinates) then
       call eigs(hesse,hesse_eigval,hesse_eigvec)
    else
       hesse_eigvec=hesse
       call jacobi(n_internal,hesse_eigvec,n_internal,hesse_eigval)
    end if
    DCALL show("Eigenvectors of Hessian",hesse_eigvec(:,:))
    DCALL show("Eigenvalues of Hessian",hesse_eigval(:))
    if(force_directed) then
       call eigs(hesse_diaval,hesse_eigval,hesse_eigvec)
       deallocate(hesse_diaval,STAT=allocopt_stat(7))
       ASSERT(allocopt_stat(7).eq.0)
       allocopt_stat(7)=1
    endif
  end subroutine eigensolve_hessian

  !*************************************************************

  subroutine hess_and_grad_project(pg,g)
    ! Purpose: Project and shift the Hessian and gradients
    !------------ Modules used -----------------------------------
    use opt_data_module, only : n_internal
    !------------ Declaration of formal parameters ---------------
    real(r8_kind)   , intent(out) :: pg(n_internal)
    real(r8_kind)   , intent(in)  :: g(n_internal)
    !------------ Declaration of local variables ---------------------
    real(r8_kind),parameter :: big=1000000_r8_kind,small_hgp=1.0e-8_r8_kind
    integer(i4_kind) :: i,j,status
    !------------ Executable code ------------------------------------

    if(allocated(q_last)) then
       deallocate(q_last,g_last,dq,dg,STAT=status)
       ASSERT(status==0)
    end if

    do_proj: if(do_projection) then
    ASSERT(size(pg)==size(pmat,1))
    !Gradient projection PMAT*G
!!$print*,'HESSE before projection'
!!$do i=1,n_primitive
!!$print*,hesse(i,:)
!!$end do
!!$print*,'GRAGS before projection'
!!$print*,g
    pg=matmul(pmat,g)
!!$print*,'GRAGS after projection'
!!$print*,pg

    !hessian projection PMAT*HESSE*PMAT+1000*(1-PMAT)
    hesse=matmul(matmul(pmat,hesse),pmat)
!!$print*,'HESSE after projection before shift'
!!$do i=1,n_primitive
!!$print*,hesse(i,:)
!!$end do

    hesse=hesse-big*pmat
    do i=1,n_internal
       hesse(i,i)=hesse(i,i)+big
    end do
!!$print*,'HESSE after projection & shift'
!!$do i=1,n_primitive
!!$print*,hesse(i,:)
!!$end do
    do i=1,n_internal
       do j=1,n_internal
          if(abs(hesse(i,j)) < small_hgp) hesse(i,j)=0.0_r8_kind
       end do
    end do
    else
       pg=g
    end if do_proj

    call eigensolve_hessian(hesse,n_internal)
!!$print*,'HESSE eigenvalues'
!!$print*,hesse_eigval
!!$print*,'HESSE eigenvectors'
!!$do i=1,n_internal
!!$print*,hesse_eigvec(i,:)
!!$end do

    if(do_projection) then
       deallocate(pmat,stat=status)
       ASSERT(status==0)
    end if

  end subroutine hess_and_grad_project
  !*************************************************************

  subroutine hesse_eval(geo_loop,do_step)
    use coortype_module
    use opt_data_module
    use coordinates_module
    use math_module
    use gradient_module, only: grad_intern
    use valence_coord_module
    use filename_module, only: inpfile
    ! Purpose: provide logic for the 'Vladimir'-technique: at the
    !          first geometry loop make a symmetry-compatible step
    !          for each component of the coordinate vector. In
    !          the second step use the gradient to get an estimate
    !          for the diagonal elements of the Hessian.
    ! ---------------------------------------------------------
    integer(kind=i4_kind),intent(inout)   :: geo_loop
    logical,intent(inout)              :: do_step
    real(kind=r8_kind),allocatable     :: q_actual(:),q_actual_p(:)
    real(kind=r8_kind),allocatable     :: q_last(:),g_last(:)
    integer(kind=i4_kind)              :: n_coor,alloc_stat
    ! --------------------------------------------------------
    call hesse_alloc()
    ! this routine assumes that 'ZMAT_FORMAT' is used thus
    ! allowing for symmetry information. The check is performed
    ! in read_input already.
    if(zmat_coordinates) then
       n_coor = n_internal
       allocate(q_actual(n_internal),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler("hesse_eval: allocation (1) failed")
    elseif (delocalized_coordinates) then
       n_coor = n_primitive
       allocate(q_actual(n_internal),q_actual_p(n_primitive),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler("hesse_eval: allocation (2) failed")
    endif
    allocate (s_local(n_coor),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler("hesse_eval: allocation (3) failed")
    allocate(q_last(n_internal),g_last(n_internal),STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    if (zmat_coordinates) then
       s_local = s
    elseif(delocalized_coordinates) then
       s_local = s_prim
       q_actual_p = q_prim
    endif
    q_actual = q

    if (geo_loop==1) then
       call step_one()
       do_step=.false.
    elseif(geo_loop==2) then
       call step_two()
       do_step=.true.
    endif

    call hesse_write(q_act=q_actual,g_act=grad_intern,h_mat=hesse,&
         dimen=n_internal,int_flag=geo_loop,temp_hesse=.false.)
    call hesse_write(q_act=q_actual,g_act=grad_intern,h_mat=hesse,&
         dimen=n_internal,int_flag=geo_loop,store=.true.,temp_hesse=.false.)
    if (geo_loop==2) then
       !restore the old point
       q=q_last
       grad_intern=g_last
    endif
    deallocate(s_local,q_actual,q_last,g_last,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("hesse_eval: deallocation (1) failed")
    if (delocalized_coordinates) then
       deallocate(q_actual_p,STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            ("hesse_eval: deallocation (2) failed")
    endif

  contains
    subroutine step_one()
      real(kind=r8_kind)           :: sign_q
      integer(kind=i4_kind)        :: i
      if (zmat_coordinates) then
         do i=1,n_coor
            sign_q = one
            if (q(i)<zero) sign_q = -one
            q(i)=q(i)+sign_q*delta_q
         enddo
      elseif(delocalized_coordinates) then
         do i=1,n_coor
            if (s_prim(i)%var) then
               sign_q=one
               if (q_actual_p(i)<zero) sign_q=-one
               q_actual_p(i)=q_actual_p(i)+sign_q*delta_q
            endif
         enddo

         q = matmul(umat_trans,q_actual_p)
      endif
      hesse=zero
    end subroutine step_one
    subroutine step_two()
      logical :: exist
      integer(kind=i4_kind)          :: alloc_stat,dummy,i
      real(kind=r8_kind),allocatable :: dq(:)

      inquire(EXIST=exist, FILE=trim(inpfile('hesse.dat')))
      if (.not.exist) then
         call error_handler&
              ("hesse_eval/step_two: file HESSE.DAT does not exist")
      endif
      call hesse_read(q_act=q_last,g_act=g_last,h_mat=hesse,&
           dimen=n_internal,int_flag=dummy,temp_hesse=.false.)
      allocate(dq(n_internal),STAT=alloc_stat)
      if (alloc_stat/=0) call error_handler &
           (" step_two: allocation (1) failed")
      dq = zero

      call delta_coordinate(q_last,q_actual,dq)

      do i=1,n_internal
         if (abs(dq(i))>small) then
            hesse(i,i) = (grad_intern(i)-g_last(i))/&
                 (q_actual(i)-q_last(i))
            if (hesse(i,i)<zero) hesse(i,i)=one
         else
            hesse(i,i)=one
         endif
      enddo
      deallocate(dq,STAT=alloc_stat)
      ASSERT(alloc_stat.eq.0)

    end subroutine step_two
  end subroutine hesse_eval

  subroutine hesse_alloc()
    use opt_data_module, only: n_internal,tsscan_sphere
    use math_module, only: zero
    !  Purpose: allocates space for hesse matrix
    !------------ Modules used ---------------------------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------
    if (.not.allocated(hesse)) then
       allocate(hesse(n_internal,n_internal),&
            hesse_inv(n_internal,n_internal),&
            STAT=allocopt_stat(14))
       ASSERT(allocopt_stat(14).eq.0)
       hesse = zero
       hesse_inv=zero
    end if

    if(tsscan_sphere) then
      allocate(hesse_sphere(n_internal-1,n_internal-1),stat=allocopt_stat(15))
      ASSERT(allocopt_stat(15).eq.0)
    endif
  end subroutine hesse_alloc

  subroutine hesse_init(start_value)
    ! Purpose : initializes the Hesse Matrix either with a
    !           unit matrix (force_field=.FALSE. or with a
    !           force field approximation.
    !           This probalby has to be replaced by something better...
    ! -------------- Modules used ---------------------------
    use opt_data_module
    use coortype_module
    use coordinates_module
    use math_module, only: one,zero
    use valence_coord_module, only: set_hesse
    ! ---------- Declaration of formal parameters -----------
    integer(kind=i4_kind),intent(in) :: start_value
    ! ---------- Declaration of local variables -------------
    integer(kind=i4_kind)    :: i,alloc_stat
    real(kind=r8_kind),allocatable :: help_mat(:),inter(:,:)

    select case (start_value)
    case( unit_matrix)
       hesse=zero
       do i=1,n_internal
          hesse(i,i) = one
          if(cart_coordinates) hesse(i,i)=hesse(i,i)/2.0_r8_kind
       enddo
    case( rough_estimate)
       if (delocalized_coordinates) then
          ! ------- ATTENTION: Alexander Voityuks Values used here -----------
          hesse=zero
          do i=1,n_internal
             hesse(i,i)=0.4_r8_kind
          enddo
          ! ------------------------------------------------------------------
       elseif (zmat_coordinates ) then
          hesse=zero
          do i=1,n_internal
             if  (s(i)%typ == b_length) then
                hesse(i,i)=0.5_r8_kind
             elseif (s(i)%typ == b_angle ) then
                hesse(i,i) = 0.1_r8_kind
             elseif (s(i)%typ == d_angle) then
                hesse(i,i) = 0.05_r8_kind
             else
                call error_handler(" invalid value detected by hesse_init")
             endif
          enddo

       endif
    case (force_field)
       hesse=zero
       allocate(help_mat(n_primitive),inter(n_primitive,n_internal),&
            STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            (" hesse_init: allocation (1) failed")
       help_mat=zero
       inter=zero
       call set_hesse(s_prim,n_primitive,help_mat)
       where (help_mat == zero )
          help_mat = one
       end where
       do i=1,n_primitive
          inter(i,:) = help_mat(i)*umat(i,:)
       enddo
       hesse=matmul(umat_trans,inter)
       deallocate(help_mat,inter,STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            (" hesse_init: deallocation (1) failed")
    case default
       call error_handler(" hesse_init: error ")
    end select
  end subroutine hesse_init
  !*************************************************************

  subroutine hesse_update(dq,dg,gm)
    ! Purpose: updates the hesse matrix and its inverse.
    !          After update, the hesse matrix is written
    !          to file.
    !          Algorithms implemented for update:
    !          -DFP      (Davidon-Fletcher-Powell) 			 [1]
    !          -BFGS     (Broyden-Fletcher-Goldfarb-Shanno)		 [1]
    !          -CG       (Fletcher-Reeves) (Hesse-matrix used formally !)[?]
    !          -BFGS_SR1 (SR1+BFGS)                                      [3]
    ! References:
    !          [1] H.B. Schlegel, Adv.Chem.Phys. 67, 249 (1987)
    !          (Review Article) and references therein
    !          [2] Farkas, O.,Schlegel, H.B. Methods for optimizing
    !          large molecules. II Qudratic Search. // J. Phys. Chem. -
    !          1999. - V. 111. - N 24. - P. 10806-10814
    !
    ! -------------- Modules used ---------------------------
    use opt_data_module
    use math_module
    use gradient_module, only: grad_intern
    USE_DEBUG
    ! ---------- Declaration of formal parameters -----------
    real(kind=r8_kind),intent(in)   :: dq(:),dg(:)
    real(kind=r8_kind),intent(in)   :: gm
    ! ---------- Declaration of local variables -------------
    logical                         :: hesse_methods!,update_bfgs_sr1
    real(kind=r8_kind)   :: bfgs,dfp,gq,ghg,qhq,t1,t2,t3,s1,qq,vh,vv,phi_k,sk,tht
    real(kind=r8_kind),allocatable  :: help(:),hesse_help(:,:),v(:),m_bfgs(:,:)
    real(kind=r8_kind),parameter    :: ac_angle=5_r8_kind,small1=1.0e-8_r8_kind
    integer(kind=i4_kind)           :: alloc_stat,i,k,j
    if (print_debug) then
       write(OPT_STDOUT,*)' Difference in grad_intern :'
       write(OPT_STDOUT,1011)(dg(i),i=1,n_internal)
       write(OPT_STDOUT,*)' '
       write(OPT_STDOUT,*)' Difference in q :'
       write(OPT_STDOUT,1011)(dq(i),i=1,n_internal)
       write(OPT_STDOUT,*)' '
    endif
1011 format(2x,10(f12.8,2x))


    ! determine method of update
    hesse_methods = .not.update_cg ! .AND. (anything else ...) [AG]
    if ( update_dfp ) then
       dfp=one
    else
       dfp=zero
    endif
    if ( update_bfgs ) then
       bfgs = one
    else
       bfgs = zero
    endif

   hes_met: if(hesse_methods) then
    allocate(help(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler(" hesse_update: allocation (3) failed")
    if(.not.cart_coordinates) then
       allocate(hesse_help(n_internal,n_internal),STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    end if

    if (abs_value(dq)/=zero.and.abs_value(dg)/=zero) then
       if ( .not.update_direct ) then
          if( print_debug) then
             write(OPT_STDOUT,*)" Inverse Hessian before Update"
             call print_matrix(hesse_inv,n_internal,n_internal,10_i4_kind)
          endif
          write(io_flepo,*)"hesse_update: largest and smallest elements of Hessian"
          write(io_flepo,*)"               BEFORE update:",&
               maxval(reshape(hesse,(/n_internal*n_internal/))),&
               minval(reshape(hesse,(/n_internal*n_internal/)))
!         write(OPT_STDOUT,*)"               BEFORE update:",&
!              maxval(reshape(hesse_inv,(/n_internal*n_internal/))),&
!              minval(reshape(hesse_inv,(/n_internal*n_internal/)))
!          write(OPT_STDOUT,*)" "
!          write(io_flepo,*)" "
          ! first update inverse
          gq = dot_product(dg,dq)
          if(abs(gq) < small1) gq=sign(small1,gq)

          if(.not.cart_coordinates) then
             help = matmul(hesse_inv,dg)
             ghg = dot_product(dg,help)

             hesse_help = zero
          else
             help = matmul(hesse,dq)
             qhq = dot_product(dq,help)
          end if

          if(.not.cart_coordinates) then
          do i=1,n_internal
             s1 = dq(i)/gq
             do k=1,n_internal
                t1 = s1*dq(k)
                t2 = - sum(dg*hesse_inv(i,:)) * &
                     sum(dg*hesse_inv(:,k))/ghg

                t3 = ghg*(dq(i)/gq-sum(hesse_inv(i,:)*dg)/ghg)* &
                     (dq(k)/gq-sum(hesse_inv(:,k)*dg)/ghg)

                hesse_help(i,k) =  dfp*( t1 + t2 ) + bfgs*t3
             enddo
          enddo
          hesse_inv = hesse_inv + hesse_help
          else
             !update direct hessian
             do i=1,n_internal
                do k=1,n_internal
                   hesse(i,k)=hesse(i,k)+dg(i)*dg(k)/gq-help(i)*help(k)/qhq
                end do
             end do
!!$print*,'HESSE UP',hesse
          end if
          if(.not.cart_coordinates) then
             hesse=hesse_inv
             call invert_matrix(n_internal,hesse)
          end if
          if( print_debug) then
             write(OPT_STDOUT,*)" Inverse Hessian after Update with gq and ghg"
             call print_matrix(hesse_inv,n_internal,n_internal,10_i4_kind)
             write(OPT_STDOUT,*)
          endif
          write(io_flepo,*)"               AFTER  update:",&
               maxval(reshape(hesse,(/n_internal*n_internal/))),&
               minval(reshape(hesse,(/n_internal*n_internal/)))
!VVP gq,ghg - i don't see reason of print this value.
!               minval(reshape(hesse_inv,(/n_internal*n_internal/))), gq,ghg
!          write(io_flepo,*)" "
!          write(OPT_STDOUT,*)"               AFTER update:",&
!               maxval(reshape(hesse_inv,(/n_internal*n_internal/))),&
!               minval(reshape(hesse_inv,(/n_internal*n_internal/)))
!          write(OPT_STDOUT,*)" "

       else ! if 'update_direct' is set
          if (print_debug) then
             write(OPT_STDOUT,*)" "
             write(OPT_STDOUT,*)" Direct Hesse before update :"
             call print_matrix(hesse,n_internal,n_internal,10_i4_kind)
          endif
          hesse_help =zero
          ! now update direct hesse
          ! Construct help matrix:
          help = matmul(hesse,dq)
          ! Calculate scale parameter
          gq = dot_product(dg,dq)
          ! Calculate scale parameter
          qhq = dot_product(dq,help)
          if(update_dfp.or.update_bfgs) then
             do i=1,n_internal
                s1 = dg(i)/gq
                do k=1,n_internal
                   t1 = s1*dg(k)
                   t2 = -sum(dq*hesse(i,:)) * &
                        sum(dq*hesse(:,k))/qhq
                   t3 = qhq* (dg(i)/gq - sum(hesse(i,:)*dq)/qhq) * &
                        (dg(k)/gq-sum(hesse(:,k)*dq)/qhq)

                   hesse_help(i,k) = dfp*( t1 + t2 ) + bfgs*t3

                enddo
             enddo
          end if
!New scheme for update of hessian by BFGS formula:
          if (update_bfgs) then
            write(io_flepo,*)"hesse_update: largest and smallest elements of Hessian"
            write(io_flepo,*)"               BEFORE update:",&
            maxval(reshape(hesse,(/n_internal*n_internal/))),&
            minval(reshape(hesse,(/n_internal*n_internal/)))
           allocate(m_bfgs(n_internal,n_internal),STAT=alloc_stat)
            if (alloc_stat/=0) call error_handler&
            (" hesse_update: allocation m_bfgs failed")
            m_bfgs=zero
            sk=abs(gq**2/dot_product(dg,dg)*dot_product(dq,dq))
!Construct bfgs correction
            tht=acos(-dot_product(g_last,dq)/(dsqrt(dot_product(dq,dq))*&
            dsqrt(dot_product(g_last,g_last))))*convert1
#ifdef WITH_ISNAN
            if(isNaN(tht)) tht=zero
#else
            ! if not optimized away, use the fact that NaN /= NaN
            if(1+tht /= 1+tht) tht=zero
#endif
            DPRINT  "hesse_main: theta=", tht
            if((dsqrt(sk)>small).and.(abs(tht-90_r8_kind).gt.ac_angle)) then
              do i=1,n_internal
               do j=1,n_internal
                m_bfgs(i,j)=dg(i)*dg(j)/gq-help(i)*help(j)/qhq
               enddo
              enddo
             if(print_debug)then
               DCALL show("BFGS Correction",m_bfgs(:,:))
             endif
            else
             print *,"Skip update Hessian"
             write(io_flepo,'(A)')"Skip update hessian"
            endif
!Symmetriezired correction to current hessian:
          hesse_help=(m_bfgs+transpose(m_bfgs))/two
          deallocate(m_bfgs,STAT=alloc_stat)
          if (alloc_stat/=0) call error_handler&
                  (" hesse_update: deallocation m_bfgs failed")
          endif

          if(update_powell.or.update_bofill) then
             allocate(v(n_internal),stat=alloc_stat)
             if (alloc_stat/=0) call error_handler&
                  (" hesse_update: allocation v failed")
             qq=dot_product(dq,dq)
             v=dg-help
             vh=dot_product(v,dq)
             vv=dot_product(v,v)
             do i=1,n_internal
                do k=1,n_internal
                   hesse_help(i,k)=1.0_r8_kind/qq*(v(i)*dq(k)+v(k)*dq(i)-&
                        1.0_r8_kind/qq*dq(i)*dq(k)*vh)
                end do
             end do
             if(update_bofill) then
                phi_k=1.0_r8_kind-vh*vh/qq/vv
             do i=1,n_internal
                do k=1,n_internal
                   hesse_help(i,k)=phi_k*hesse_help(i,k)+(1.0_r8_kind-phi_k)*&
                        v(k)*v(i)/vh
                end do
             end do
             write(OPT_STDOUT,*) "phi_k:",phi_k
             end if
             deallocate(v,stat=alloc_stat)
             if (alloc_stat/=0) call error_handler&
                  (" hesse_update: deallocation v failed")
          end if
          hesse = hesse + hesse_help
          write(io_flepo,*)"               AFTER  update:",&
          maxval(reshape(hesse,(/n_internal*n_internal/))),&
          minval(reshape(hesse,(/n_internal*n_internal/)))

          if (print_debug) then
             write(OPT_STDOUT,*)" "
             write(OPT_STDOUT,*)" Direct Hesse after Update "
             call print_matrix(hesse,n_internal,n_internal,10_i4_kind)
          endif
          hesse_inv=hesse
          call invert_matrix(n_internal,hesse_inv)
       endif! update_direct


       if(.not.cart_coordinates) then
          deallocate(hesse_help,STAT=alloc_stat)
          if (alloc_stat/=0) &
               stop ' hesse_update: deallocation (1) failed'
       end if
       deallocate(help,STAT=alloc_stat)
       if (alloc_stat/=0) &
            stop ' hesse_update: deallocation (2) failed'
    else !if dq is zero
     write(OPT_STDOUT,*)" hesse_update: dq,i.e. dg was found to be zero. Skip hesse_update"
       hesse_inv=hesse
       if(.not.cart_coordinates) then
          call invert_matrix(n_internal,hesse_inv)
       end if
    endif! dq=zero?
   else hes_met! non-Hesse-methods (CG etc.)
        if(update_CG) then
           hesse_inv = zero
           do i=1,n_internal
              do k=1,n_internal
                 if(i==k) hesse_inv(i,k)= one
                 hesse_inv(i,k) = hesse_inv(i,k) - dq(i)*grad_intern(k)/gm
              end do
           end do
           hesse = hesse_inv ! just to save inverse HEssian formally
        DPRINT 'update_CG done for gm=',gm
        DPRINT hesse_inv(1,:)
        DPRINT hesse_inv(2,:)
        DPRINT hesse_inv(3,:)
        else
           call error_handler("hesse_update : What is it ?")
        end if
   endif hes_met

  end subroutine hesse_update
  !*************************************************************
  subroutine hesse_write(q_act,g_act,q1,g1,h_mat,dimen,int_flag,deloc,store,temp_hesse)
    ! Purpose: writes the hesse matrix and the actual position
    !          'q_act' and gradient 'g_act' to file.
    !          The parameter 'int_flag' has been introduced for
    !          the routines that calculation the initial hessian
    !          by a step-by-step elongation of the internal
    !          variables. There 'int_flag' will be the sym_type
    !          of the last variable elongated.
    use opt_data_module
    use filename_module, only: inpfile
    ! ----- Declaration of formal parameters ------------------
    real(kind=r8_kind),intent(in)     :: q_act(:),g_act(:),h_mat(:,:)
    real(kind=r8_kind),optional       :: q1(:),g1(:)
    integer(kind=i4_kind),intent(in)  :: dimen,int_flag
    logical,intent(in),optional       :: deloc,store
    logical,intent(in)                :: temp_hesse
    ! ----- Declaration of local vrariables -------------------
    integer(kind=i4_kind) :: i
    logical               :: deloc_local
    character(len=4) ::  char_int_flag


    if (present(deloc)) then
       if (deloc) then
          deloc_local=.true.
       else
          deloc_local=.false.
       endif
    else
       deloc_local = .false.
    endif

    if (deloc_local) then ! check if q_act has the right dimension
       if (ubound(q_act,1) /= n_primitive ) then
          write(OPT_STDOUT,*)" hesse_write: for writing information q_act does not"
          write(OPT_STDOUT,*)"              have the right dimension."
          write(OPT_STDOUT,*)"              Dimension found on entry to write_hesse:",ubound(q_act,1)
          write(OPT_STDOUT,*)"              Correct dimension                      :",n_primitive
          call error_handler(" ")
       endif
    endif

    if (dimen /= ubound(g_act,1)) then
       write(OPT_STDOUT,*)" hesse_write: trying to write Hessian with wrong dimension"
       write(OPT_STDOUT,*)"              Dimension on entry to this routine :",dimen
       write(OPT_STDOUT,*)"              Actual dimension of Gradient       :",ubound(g_act,1)
       call error_handler(" ")
    endif
    if(present(store) ) then
       if(.not.store) return
       write(char_int_flag,'(i4)') int_flag
       io_hesse=openget_iounit(file=trim(opt_dir)//'/hesse.'//adjustl(char_int_flag), &
            status='replace',form='formatted')
    else
    if(temp_hesse) then
       io_hesse=openget_iounit(file=trim(inpfile('hesse_t.dat')), &
            status='replace',form='formatted')
    else
       io_hesse=openget_iounit(file=trim(inpfile('hesse.dat')), &
            status='replace',form='formatted')
    DPRINT 'io_hesse=',io_hesse
    endif
    endif

    write(io_hesse,*)int_flag
    write(io_hesse,*)q_act
    write(io_hesse,*)g_act
    if (present(q1).and.present(g1)) then
       write(io_hesse,*)q1
       write(io_hesse,*)g1
    endif
    do i=1,dimen
       write(io_hesse,*)h_mat(:,i)
    enddo
       call returnclose_iounit(io_hesse)

  end subroutine hesse_write

  !*************************************************************
  subroutine hesse_read(q_act,g_act,q1,g1,h_mat,dimen,int_flag, &
                        deloc,temp_hesse)
    ! read the hessematrix from file. Musrt be consistent
    ! with 'hesse_write'
    ! --------------------------------------------------------
    use opt_data_module
    use filename_module, only: inpfile
    real(kind=r8_kind),intent(inout)     :: q_act(:),g_act(:),h_mat(:,:)
    real(kind=r8_kind),optional          :: q1(:),g1(:)
    integer(kind=i4_kind),intent(inout)  :: dimen,int_flag
    logical,intent(in),optional          :: deloc
    logical,intent(in)                   :: temp_hesse
    ! ----- Declaration of local vrariables -------------------
    integer(kind=i4_kind) :: i
    integer(kind=i4_kind) :: stat
    logical               :: deloc_local

    if (present(deloc)) then
       if (deloc) then
          deloc_local=.true.
       else
          deloc_local=.false.
       endif
    else
       deloc_local = .false.
    endif

    if (deloc_local) then ! check if q_act has the right dimension
       if (ubound(q_act,1) /= n_primitive ) then
          write(OPT_STDOUT,*)" hesse_read: for reading information q_act does not"
          write(OPT_STDOUT,*)"             have the right dimension."
          write(OPT_STDOUT,*)"             Dimension found on entry to read_hesse:",ubound(q_act,1)
          write(OPT_STDOUT,*)"             Correct dimension                     :",n_primitive
          call error_handler(" ")
       endif
    endif
    if (dimen /= ubound(g_act,1)) then
       write(OPT_STDOUT,*)" hesse_read: trying to read Hessian with wrong dimension"
       print*,'dimen on file ',dimen
       print*,'dimension of gradient :',ubound(g_act,1)
       call error_handler(" ")
    endif

   if(temp_hesse) then
    io_hesse=openget_iounit(form='formatted', file=trim(inpfile('hesse_t.dat')))
   else
    io_hesse=openget_iounit(form='formatted', file=trim(inpfile('hesse.dat')))
    DPRINT 'io_hesse=',io_hesse !(1)
    DPRINT 'hesse.dat read from '//trim(inpfile('hesse.dat'))
   endif

    !
    ! Iteration count:
    !
    read(io_hesse, *, iostat=stat) int_flag
    ASSERT(stat==0)

    !
    ! Previous (?) position:
    !
    read(io_hesse, *, iostat=stat) q_act(:)
    ASSERT(stat==0)

    !
    ! Previous (?) gradient:
    !
    read(io_hesse, *, iostat=stat) g_act(:)
    ASSERT(stat==0)

    !
    ! FIXME: ?
    !
    if (present(g1).and.present(q1)) then
       read(io_hesse, *, iostat=stat) q1(:)
       ASSERT(stat==0)

       read(io_hesse, *, iostat=stat) g1(:)
       ASSERT(stat==0)
    endif

    !
    ! Current state of the hessian:
    !
    do i=1,dimen
       read(io_hesse, *, iostat=stat) h_mat(:,i)
       ASSERT(stat==0)
    enddo
    call returnclose_iounit(io_hesse)
  end subroutine hesse_read

   subroutine rewrite_hesse(geo_loop)
     use opt_data_module
     ! Purpose :
     !----- Declaration of formal paramaters -------------------
     integer(kind=i4_kind), intent(in) :: geo_loop
     !----- Declaration of local variables ---------------------
!    logical    :: exist
     real(kind=r8_kind), allocatable :: q_temp(:),g_temp(:),hesse_temp(:,:)
     integer(kind=i4_kind)           :: n_dim,int_tmp,alloc_stat
     !----- executable code ------------------------------------

     if(geo_loop==1) return
     if( zmat_coordinates        ) n_dim = n_internal
     if( delocalized_coordinates ) n_dim = n_primitive

     allocate(q_temp(n_dim),g_temp(n_dim), hesse_temp(n_dim,n_dim), &
                                                     stat=alloc_stat)
     if (alloc_stat/=0) &
          call error_handler(" rewrite_hesse : allocation  failed")
     write(*,*) "============== rewrite hesse =============="
     write(*,*) "N_dim ", n_dim, "Cycle ", geo_loop
     write(*,*) "==========================================="
     call hesse_read(q_act=q_temp, g_act=g_temp, h_mat=hesse_temp,&
                     dimen=n_internal, int_flag=int_tmp,          &
                     deloc=delocalized_coordinates,temp_hesse=.true.)
     call hesse_write(q_act=q_temp,g_act=g_temp, h_mat=hesse_temp,&
                      dimen=n_internal,int_flag=int_tmp,          &
                      deloc=delocalized_coordinates,temp_hesse=.false.)

     deallocate(q_temp, g_temp, hesse_temp, stat=alloc_stat)
     if (alloc_stat/=0) &
          call error_handler(" rewrite_hesse : deallocation failed")

   end subroutine rewrite_hesse

  subroutine hesse_prepare_internal(hesse,error)
    !
    ! Purpose: read in a Hesse matrix in cartesian coordinates
    !          if it exist, augment it with other contributions
    !          and transform to internal coordinates.
    !
    ! 3. Now also is used to prepare suiteable for optimization
    !    hessian in cartesian coordinates
    !
    !----- Modules used ---------------------------------------
    use opt_data_module, only: OPT_STDOUT, n_atoms, n_dummy, charge, xyz
    use opt_data_module, only: tsscan_mix, epe_forces
    use gradient_module, only: dervs_cartes
    use filename_module, only: inpfile
#ifndef FPP_OPTIMIZER
    use gradient_data_module, only: optimizer_write_cart_hess
#endif
    implicit none
    !----- Declaration of formal paramaters -------------------
    real(r8_kind), intent(out) :: hesse(:,:)
    logical      , intent(out) :: error
    ! *** end of interface ***

    !----- Declaration of local variables ---------------------
    logical    :: exist
    real(kind=r8_kind),allocatable :: hesse_cartesian(:,:),hesse_cartesian_dummy(:,:)
    integer(kind=i4_kind):: alloc_stat,io_stat,io_hesse_cart
    !----- executable code ------------------------------------

    error=.false.
#ifdef WITH_EFP
    if(efp .and. n_efp > 0 .and. .not.efp_fixed) then
       error=.true.
       return
    end if
#endif

    allocate(hesse_cartesian(3*n_atoms,3*n_atoms),&
             hesse_cartesian_dummy(3*(n_atoms+n_dummy),3*(n_atoms+n_dummy)),&
         STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

   if(tsscan_mix) then
     ! Dont get it from file but calculate second dervs:
     call tsscan_hesse_cartesian(hesse_cartesian)
   else
     inquire(EXIST=exist, FILE=trim(inpfile('hesse_cartesian.dat')))
     if (exist) then
        DPRINT 'hesse_cart_to_internal: found hesse_cartesian.dat'
        io_hesse_cart=openget_iounit(status='unknown',form='formatted',file=&
             trim(inpfile('hesse_cartesian.dat')))

     else
        DPRINT 'hesse_cart_to_internal: hesse_cartesian.dat not found'
        write(OPT_STDOUT,*)" hesse_cart_to_internal: cartesian Hesse does not exist"
        error=.true.
        deallocate(hesse_cartesian,hesse_cartesian_dummy,STAT=alloc_stat)
        ASSERT(alloc_stat.eq.0)
        return
     endif

     read(io_hesse_cart,*,IOSTAT=io_stat) hesse_cartesian
     if( io_stat/=0 )then
       print *,'hesse_cart_to_internal: ERROR reading hesse_cartesian.dat, iostat=',io_stat
       print *,'Maybe your hesse_cartesian.dat corresponds to a different job?'
       ABORT('error reading hesse_cartesian.dat, see tty')
     endif
     call returnclose_iounit(io_hesse_cart)! if it is not required futher on?
    endif

    DPRINT 'hesse_prepare_internal: hesse_cartesian read/calculated'

    call add_dummies(hesse_cartesian,hesse_cartesian_dummy)

    if( allocated(dervs_cartes) )then
      call add_dervs(dervs_cartes,hesse_cartesian_dummy)
    endif

    if(epe_forces)then
      call add_epe_contrib(hesse_cartesian_dummy)
    endif

#ifndef FPP_OPTIMIZER
    ! this diagonalizes the hessian and prints the frequencies on stdout,
    ! the name is misleading. To disentangle standalone optimizer from PG
    ! this subroutine from gradient_data_module should be better moved
    ! into a separate file that can be used both in PG and optimizer.
    call optimizer_write_cart_hess(stdout_unit,hesse_cartesian_dummy,charge,xyz,n_atoms+n_dummy)
#else
    print *,'WARNING: skipping diagonalization of hessian for technical reasons'
    print *,'WARNING: may be fixed on request'
#endif
    call diff_to_numhess(hesse_cartesian_dummy)

    print*,' -------------------------- hesse -------------------------'

    ! now convert it to zmat-, delocalized, or other coordinates:
    call hesse_cart_to_internal(hesse_cartesian_dummy,hesse)

    deallocate(hesse_cartesian,hesse_cartesian_dummy,STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)


  contains ! of hesse_prepare_internal

    subroutine diff_to_numhess(hesse)
      real(kind=r8_kind), intent(in) :: hesse(:,:)
      ! *** end of interface ***

      integer(kind=i4_kind):: i,ihess
      logical:: num_cart_hes
      real(kind=r8_kind) cart_hes(size(hesse,1),size(hesse,1))

      inquire(file=trim(inpfile('HESS.cart')), exist=num_cart_hes)
      print*,'file exist '//trim(inpfile('HESS.cart')), num_cart_hes
      num_cart_hes=.false.
      if(num_cart_hes) then
        ihess=openget_iounit(status='old',form='formatted',&
                             file=trim(inpfile('HESS.cart')))
        DPRINT 'ihess=',ihess
        do i=1,size(hesse,1)
         read(ihess,*) cart_hes(i,:)
        enddo
        call returnclose_iounit(ihess)
        print*,'diff hesse cart_hes',sum((hesse-cart_hes)**2)
!       hesse=cart_hes
      endif
    end subroutine diff_to_numhess

    subroutine add_epe_contrib(hesse_cartesian_dummy)
      ! Contributions to Cartesian Hessian due to short-range interactions between
      ! embedded QM cluster and MM environment (EPE method)
#ifdef NEW_EPE
      use qm_epe_interface_module
#endif
      use opt_data_module, only: n_atoms,n_dummy,x,y,z,charge,index_unique
      use slspar_module
      use math_module, only: pi
      use constants, only: angstrom
      implicit none
      real(r8_kind), intent(inout) :: hesse_cartesian_dummy(:,:) ! (3(NA+NX))^2
      ! *** end of interface ***

      real(r8_kind) :: r1(3),r2(3),r3(3),rr(3),dr,dr7,dr8,dr9,dr10
      real(r8_kind) :: rr1(3),rr2(3),dr1,dr2,e1(3),e2(3)
      real(r8_kind) :: dot_prod,cos_th,theta,dE_dt,d2E_dt2,sin_th
      real(r8_kind) :: p(10),dE_dr,d2E_dr2,f1(3),f2(3),f3(3),ff(6,6),ff1(9,9),fstore
      real(r8_kind) :: deg2rad
      integer(i4_kind) :: i,j,k,k2,kk1,kk2,kk3,l,m,m1,l1,l2,l3
      integer(i4_kind) :: ia_1,ia_2,ia_3,i1,j1,k1,num_3b_links
#ifdef NEW_EPE
      integer(i4_kind) :: it,jt
      real(r8_kind) :: cut
#else
      integer(i4_kind) :: nr
      real(r8_kind) :: mini,maxi
#endif

      real(r8_kind), parameter :: one=1.0_r8_kind,two=2.0_r8_kind,six=6.0_r8_kind,eight=8.0_r8_kind
      real(kind=r8_kind), parameter :: small=1.0e-8_r8_kind

!      hesse_cartesian_dummy=0.0_r8_kind

      i_lab: do i=1,n_atoms+n_dummy
         r1(1)=x(i); r1(2)=y(i); r1(3)=z(i)
         if(charge(i)-aint(charge(i)).lt.0.0009_r8_kind &
              .or. index_unique(i).eq.0) cycle i_lab

#ifdef NEW_EPE
         it=at_type(i)
#endif

         if(.true.) then
         j_lab: do j=1,i-1
            r2(1)=x(j); r2(2)=y(j); r2(3)=z(j)
!!!         print*, 'r2 defined'
            if(charge(j)-aint(charge(j)).lt.0.0009_r8_kind &
                 .or. index_unique(j).eq.0) cycle j_lab

#ifdef NEW_EPE
            jt=at_type(j)
            p(1)=interface_ff(it,jt)%b
            p(2)=interface_ff(it,jt)%r
            p(3)=interface_ff(it,jt)%c
            p(4)=interface_ff(it,jt)%d
            cut=interface_ff(it,jt)%cut
            if(abs(cut) > 0.1_r8_kind) then
#else
            mini=min(charge(i),charge(j))
            maxi=max(charge(i),charge(j))
            call get_slsp(mini,maxi,0,nr)
            if(nr /= 0) then
               p(1)=par%item%b
               p(2)=par%item%r
               p(3)=par%item%c
               p(4)=par%item%d
#endif

               dE_dr=0.0_r8_kind; d2E_dr2=0.0_r8_kind

               rr=r1-r2
               dr=sqrt(dot_product(rr,rr)) / angstrom
#ifdef NEW_EPE
               if(dr <= cut) then
#else
               if(dr <= par%item%cutoff) then
#endif
                  dr7=dr*dr*dr*dr*dr*dr*dr
                  dr8=dr7*dr
                  dr9=dr8*dr
                  dr10=dr9*dr

                  dE_dr=-(p(1)/p(2))*exp(-dr/p(2))+six*p(3)/dr7+eight*p(4)/dr9
                  dE_dr=dE_dr*evau / angstrom
                  d2E_dr2=(p(1)/(p(2)*p(2)))*exp(-dr/p(2))-42.0_r8_kind*p(3)/dr8- &
                       72_r8_kind*p(4)/dr10
                  d2E_dr2=d2E_dr2*evau / angstrom**2
               end if

#ifndef NEW_EPE
               p(5)=par%item%k
               p(6)=par%item%r0
               p(7)=par%item%k1
               p(8)=par%item%r1
               if(p(5) /= 0.0 .and. abs(dr-p(6)).lt.0.1_r8_kind) then
                  dE_dr=dE_dr+two*p(5)*(dr-p(6)) * angstrom
                  d2E_dr2=d2E_dr2+p(5)
               end if
               if(p(7) /= 0.0 .and. abs(dr-p(8)).lt.0.1_r8_kind) then
                  dE_dr=dE_dr+two*p(7)*(dr-p(8)) * angstrom
                  d2E_dr2=d2E_dr2+p(7)
               end if
#endif
               dr=dr * angstrom
               f1=rr/dr; f2=-f1

               ff(1,1)=one/dr-rr(1)*rr(1)/dr**3
               ff(1,2)=-rr(1)*rr(2)/dr**3
               ff(1,3)=-rr(1)*rr(3)/dr**3
               ff(1,4)=-ff(1,1)
               ff(1,5)=-ff(1,2)
               ff(1,6)=-ff(1,3)
               ff(2,1)=-rr(2)*rr(1)/dr**3
               ff(2,2)=one/dr-rr(2)*rr(2)/dr**3
               ff(2,3)=-rr(2)*rr(3)/dr**3
               ff(2,4)=-ff(2,1)
               ff(2,5)=-ff(2,2)
               ff(2,6)=-ff(2,3)
               ff(3,1)=-rr(3)*rr(1)/dr**3
               ff(3,2)=-rr(3)*rr(2)/dr**3
               ff(3,3)=one/dr-rr(3)*rr(3)/dr**3
               ff(3,4)=-ff(3,1)
               ff(3,5)=-ff(3,2)
               ff(3,6)=-ff(3,3)
               ff(4,1:3)=ff(1:3,4)
               ff(5,1:3)=ff(1:3,5)
               ff(6,1:3)=ff(1:3,6)
               ff(4,4:6)=ff(1,1:3)
               ff(5,4:6)=ff(2,1:3)
               ff(6,4:6)=ff(3,1:3)

               k1=3*(i-1)
               k2=3*(j-1)
               do l=1,3
                  do m=1,6
                     if(m <= 3) then
                        fstore=f1(m)
                        m1=k1+m
                     elseif(m > 3) then
                        m1=k2+(m-3)
                        fstore=f2(m-3)
                     end if
                     l1=k1+l; l2=k2+l
                     hesse_cartesian_dummy(l1,m1)=hesse_cartesian_dummy(l1,m1)+ &
                          d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
                     hesse_cartesian_dummy(l2,m1)=hesse_cartesian_dummy(l2,m1)+ &
                          d2E_dr2*fstore*f2(l)+dE_dr*ff(l+3,m)
                  end do
               end do
            end if
         end do j_lab
         endif

         if(.false..and.epe_pc_shells) then !!! something wrong in this block
                                            !!! switched out up to AS correction
         jj_lab: do j=epe_kl,epe_nucen+n_atoms+n_dummy
            if(epe(j)%ant-aint(epe(j)%ant).lt.0.009999_r8_kind) cycle jj_lab

#ifdef NEW_EPE
            jt=at_type(j)
            p(1)=interface_ff(it,jt)%b
            p(2)=interface_ff(it,jt)%r
            p(3)=interface_ff(it,jt)%c
            p(4)=interface_ff(it,jt)%d
            cut=interface_ff(it,jt)%cut
            if(abs(cut) > 0.1_r8_kind) then
#else
            mini=min(charge(i),epe(j)%ant)
            maxi=max(charge(i),epe(j)%ant)
            call get_slsp(mini,maxi,0,nr)

            if(nr.ne.0) then
#endif
               r2(:3)=epe(j)%s(:3)
               rr=r1-r2
               dr=sqrt(dot_product(rr,rr)) / angstrom

               dE_dr=0.0_r8_kind; d2E_dr2=0.0_r8_kind

#ifndef  NEW_EPE
               p(9)=par%item%rqq
               p(10)=par%item%qq
               if(dr <= p(9)) then
                  dE_dr=-(p(10)/(dr*dr)) / angstrom**2
                  d2E_dr2=(two*p(10)/(dr*dr*dr)) / angstrom**3
               end if

               p(1)=par%item%b
               p(2)=par%item%r
               p(3)=par%item%c
               p(4)=par%item%d
#endif

#ifdef NEW_EPE
               if(dr <= cut) then
#else
               if(dr <= par%item%cutoff) then
#endif
                  dr7=dr*dr*dr*dr*dr*dr*dr
                  dr8=dr7*dr
                  dr9=dr8*dr
                  dr10=dr9*dr

                  dE_dr=dE_dr-((p(1)/p(2))*exp(-dr/p(2))+six*p(3)/dr7+eight*p(4)/dr9)*evau / angstrom
                  d2E_dr2=((p(1)/(p(2)*p(2)))*exp(-dr/p(2))-42.0_r8_kind*p(3)/dr8- &
                       72_r8_kind*p(4)/dr10)*evau / angstrom**2
               end if

#ifndef NEW_EPE
               p(5)=par%item%k
               p(6)=par%item%r0
               p(7)=par%item%k1
               p(8)=par%item%r1
               if(p(5) /= 0.0 .and. abs(dr-p(6)).lt.0.1_r8_kind) then
                  dE_dr=dE_dr+two*p(5)*(dr-p(6)) * angstrom
                  d2E_dr2=d2E_dr2+p(5)
               end if
               if(p(7) /= 0.0 .and. abs(dr-p(8)).lt.0.1_r8_kind) then
                  dE_dr=dE_dr+two*p(7)*(dr-p(8)) * angstrom
                  d2E_dr2=d2E_dr2+p(7)
               end if
#endif
               dr=dr * angstrom
               f1=rr/dr

               ff(1,1)=one/dr-rr(1)*rr(1)/dr**3
               ff(1,2)=-rr(1)*rr(2)/dr**3
               ff(1,3)=-rr(1)*rr(3)/dr**3
               ff(2,1)=-rr(2)*rr(1)/dr**3
               ff(2,2)=one/dr-rr(2)*rr(2)/dr**3
               ff(2,3)=-rr(2)*rr(3)/dr**3
               ff(3,1)=-rr(3)*rr(1)/dr**3
               ff(3,2)=-rr(3)*rr(2)/dr**3
               ff(3,3)=one/dr-rr(3)*rr(3)/dr**3

               k1=3*(i-1)
               do l=1,3
                  do m=1,3
                     fstore=f1(m)
                     m1=k1+m
                     l1=k1+l
                     hesse_cartesian_dummy(l1,m1)=hesse_cartesian_dummy(l1,m1)+ &
                          d2E_dr2*fstore*f1(l)+dE_dr*ff(l,m)
                  end do
               end do
            end if
        end do jj_lab
        endif
      end do i_lab


!!!      hesse_cartesian_dummy=0.0
      if(n_types_central_atoms_3body > 0 ) then
         deg2rad=pi/180.0_r8_kind
         num_3b_links=0
         fst: do i=1,n_tetrahedrons
            do l=1,5
               if (tetra_atoms(i,l) <= epe_nucen+n_atoms+n_dummy) goto 1
            enddo
            cycle fst
1           ia_1=tetra_atoms(i,1)
#ifdef NEW_EPE
            i1=at_type(ia_1)
#else
            i1=epe(ia_1)%k
#endif
            scnd :do j=2,4
               ia_2=tetra_atoms(i,j)
               if(ia_2.eq.0) exit
#ifdef NEW_EPE
               j1=at_type(ia_2)
#else
               j1=epe(ia_2)%k
#endif
               thrd: do k=j+1,5
                  ia_3=tetra_atoms(i,k)
                  if(ia_3.eq.0) exit
                  ! all atoms  outside cluster
                  if(ia_1.ge.epe_kl.and.ia_2.ge.epe_kl &
                       .and.ia_3.ge.epe_kl) cycle thrd
#ifdef NEW_EPE
                  k1=at_type(ia_3)
#else
                  k1=epe(ia_3)%k
#endif
                  r2=epe(ia_1)%s; r1=epe(ia_2)%s; r3=epe(ia_3)%s
                  rr1=r1-r2; rr2=r3-r2
                  dr1=sqrt(dot_product(rr1,rr1)); dr2=sqrt(dot_product(rr2,rr2));

                  dot_prod=dot_product(rr1,rr2)
                  cos_th=dot_prod/(dr1*dr2)
                  theta=acos(cos_th)

                  dE_dt=ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)*evau
                  d2E_dt2=ki(j1,i1,k1)*evau

                  sin_th=sin(theta)
                  sin_th = sign(max(small,abs(sin_th)),sin_th)

                  e1=rr1/dr1
                  e2=rr2/dr2
                  num_3b_links=num_3b_links+1
                  f1=(cos_th*e1-e2)/(dr1*sin_th)
                  f3=(cos_th*e2-e1)/(dr2*sin_th)
                  f2=((dr1-dr2*cos_th)*e1+(dr2-dr1*cos_th)*e2)/(dr1*dr2*sin_th)

                  ff1(1,1)=(-sin_th*f1(1)*e1(1)+cos_th*(one/dr1-rr1(1)*rr1(1)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(-rr1(1)/(dr1**3*sin_th)-cos_th*f1(1)/(dr1*sin_th**2))

                  ff1(1,2)=(-sin_th*f1(2)*e1(1)+cos_th*(-rr1(1)*rr1(2)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(-rr1(2)/(dr1**3*sin_th)-cos_th*f1(2)/(dr1*sin_th**2))

                  ff1(1,3)=(-sin_th*f1(3)*e1(1)+cos_th*(-rr1(1)*rr1(3)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(-rr1(3)/(dr1**3*sin_th)-cos_th*f1(3)/(dr1*sin_th**2))

                  ff1(1,4)=(-sin_th*f2(1)*e1(1)+cos_th*(-one/dr1+rr1(1)*rr1(1)/dr1**3)+ &
                       one/dr2-rr2(1)*rr2(1)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(rr1(1)/(dr1**3*sin_th)-cos_th*f2(1)/(dr1*sin_th**2))

                  ff1(1,5)=(-sin_th*f2(2)*e1(1)+cos_th*(rr1(1)*rr1(2)/dr1**3)- &
                       rr2(1)*rr2(2)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(rr1(2)/(dr1**3*sin_th)-cos_th*f2(2)/(dr1*sin_th**2))

                  ff1(1,6)=(-sin_th*f2(3)*e1(1)+cos_th*(rr1(1)*rr1(3)/dr1**3)- &
                       rr2(1)*rr2(3)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(rr1(3)/(dr1**3*sin_th)-cos_th*f2(3)/(dr1*sin_th**2))

                  ff1(1,7)=(-sin_th*f3(1)*e1(1)-one/dr2+rr2(1)*rr2(1)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(-cos_th*f3(1)/(dr1*sin_th**2))

                  ff1(1,8)=(-sin_th*f3(2)*e1(1)+rr2(1)*rr2(2)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(-cos_th*f3(2)/(dr1*sin_th**2))

                  ff1(1,9)=(-sin_th*f3(3)*e1(1)+rr2(1)*rr2(3)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(1)-e2(1))*(-cos_th*f3(3)/(dr1*sin_th**2))
                  !........................................................
                  ff1(2,1)=(-sin_th*f1(1)*e1(2)+cos_th*(-rr1(2)*rr1(1)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(-rr1(1)/(dr1**3*sin_th)-cos_th*f1(1)/(dr1*sin_th**2))

                  ff1(2,2)=(-sin_th*f1(2)*e1(2)+cos_th*(one/dr1-rr1(2)*rr1(2)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(-rr1(2)/(dr1**3*sin_th)-cos_th*f1(2)/(dr1*sin_th**2))

                  ff1(2,3)=(-sin_th*f1(3)*e1(2)+cos_th*(-rr1(2)*rr1(3)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(-rr1(3)/(dr1**3*sin_th)-cos_th*f1(3)/(dr1*sin_th**2))

                  ff1(2,4)=(-sin_th*f2(1)*e1(2)+cos_th*(rr1(2)*rr1(1)/dr1**3)- &
                       rr2(2)*rr2(1)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(rr1(1)/(dr1**3*sin_th)-cos_th*f2(1)/(dr1*sin_th**2))

                  ff1(2,5)=(-sin_th*f2(2)*e1(2)+cos_th*(-one/dr1+rr1(2)*rr1(2)/dr1**3)+ &
                       one/dr2-rr2(2)*rr2(2)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(rr1(2)/(dr1**3*sin_th)-cos_th*f2(2)/(dr1*sin_th**2))

                  ff1(2,6)=(-sin_th*f2(3)*e1(2)+cos_th*(rr1(2)*rr1(3)/dr1**3)- &
                       rr2(2)*rr2(3)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(rr1(3)/(dr1**3*sin_th)-cos_th*f2(3)/(dr1*sin_th**2))

                  ff1(2,7)=(-sin_th*f3(1)*e1(2)+rr2(2)*rr2(1)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(-cos_th*f3(1)/(dr1*sin_th**2))

                  ff1(2,8)=(-sin_th*f3(2)*e1(2)-one/dr2+rr2(2)*rr2(2)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(-cos_th*f3(2)/(dr1*sin_th**2))

                  ff1(2,9)=(-sin_th*f3(3)*e1(2)+rr2(2)*rr2(3)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(2)-e2(2))*(-cos_th*f3(3)/(dr1*sin_th**2))
                  !........................................................
                  ff1(3,1)=(-sin_th*f1(1)*e1(3)+cos_th*(-rr1(3)*rr1(1)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(-rr1(1)/(dr1**3*sin_th)-cos_th*f1(1)/(dr1*sin_th**2))

                  ff1(3,2)=(-sin_th*f1(2)*e1(3)+cos_th*(-rr1(3)*rr1(2)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(-rr1(2)/(dr1**3*sin_th)-cos_th*f1(2)/(dr1*sin_th**2))

                  ff1(3,3)=(-sin_th*f1(3)*e1(3)+cos_th*(one/dr1-rr1(3)*rr1(3)/dr1**3))/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(-rr1(3)/(dr1**3*sin_th)-cos_th*f1(3)/(dr1*sin_th**2))

                  ff1(3,4)=(-sin_th*f2(1)*e1(3)+cos_th*(rr1(3)*rr1(1)/dr1**3)- &
                       rr2(3)*rr2(1)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(rr1(1)/(dr1**3*sin_th)-cos_th*f2(1)/(dr1*sin_th**2))

                  ff1(3,5)=(-sin_th*f2(2)*e1(3)+cos_th*(rr1(3)*rr1(2)/dr1**3)- &
                       rr2(3)*rr2(2)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(rr1(2)/(dr1**3*sin_th)-cos_th*f2(2)/(dr1*sin_th**2))

                  ff1(3,6)=(-sin_th*f2(3)*e1(3)+cos_th*(-one/dr1+rr1(3)*rr1(3)/dr1**3)+ &
                       one/dr2-rr2(3)*rr2(3)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(rr1(3)/(dr1**3*sin_th)-cos_th*f2(3)/(dr1*sin_th**2))

                  ff1(3,7)=(-sin_th*f3(1)*e1(3)+rr2(3)*rr2(1)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(-cos_th*f3(1)/(dr1*sin_th**2))

                  ff1(3,8)=(-sin_th*f3(2)*e1(3)+rr2(3)*rr2(2)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(-cos_th*f3(2)/(dr1*sin_th**2))

                  ff1(3,9)=(-sin_th*f3(3)*e1(3)-one/dr2+rr2(3)*rr2(3)/dr2**3)/(dr1*sin_th)+ &
                       (cos_th*e1(3)-e2(3))*(-cos_th*f3(3)/(dr1*sin_th**2))
                  !........................................................
                  ff1(7,1:3)=ff1(1:3,7)

                  ff1(7,4)=(-sin_th*f2(1)*e2(1)+cos_th*(-one/dr2+rr2(1)*rr2(1)/dr2**3)+ &
                       one/dr1-rr1(1)*rr1(1)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(1)-e1(1))*(rr2(1)/(dr2**3*sin_th)-cos_th*f2(1)/(dr2*sin_th**2))

                  ff1(7,5)=(-sin_th*f2(2)*e2(1)+cos_th*(rr2(1)*rr2(2)/dr2**3)- &
                       rr1(1)*rr1(2)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(1)-e1(1))*(rr2(2)/(dr2**3*sin_th)-cos_th*f2(2)/(dr2*sin_th**2))

                  ff1(7,6)=(-sin_th*f2(3)*e2(1)+cos_th*(rr2(1)*rr2(3)/dr2**3)- &
                       rr1(1)*rr1(3)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(1)-e1(1))*(rr2(3)/(dr2**3*sin_th)-cos_th*f2(3)/(dr2*sin_th**2))

                  ff1(7,7)=(-sin_th*f3(1)*e2(1)+cos_th*(one/dr2-rr2(1)*rr2(1)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(1)-e1(1))*(-rr2(1)/(dr2**3*sin_th)-cos_th*f3(1)/(dr2*sin_th**2))

                  ff1(7,8)=(-sin_th*f3(2)*e2(1)+cos_th*(-rr2(1)*rr2(2)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(1)-e1(1))*(-rr2(2)/(dr2**3*sin_th)-cos_th*f3(2)/(dr2*sin_th**2))

                  ff1(7,9)=(-sin_th*f3(3)*e2(1)+cos_th*(-rr2(1)*rr2(3)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(1)-e1(1))*(-rr2(3)/(dr2**3*sin_th)-cos_th*f3(3)/(dr2*sin_th**2))
                  !........................................................
                  ff1(8,1:3)=ff1(1:3,8)

                  ff1(8,4)=(-sin_th*f2(1)*e2(2)+cos_th*(rr2(2)*rr2(1)/dr2**3)- &
                       rr1(2)*rr1(1)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(2)-e1(2))*(rr2(1)/(dr2**3*sin_th)-cos_th*f2(1)/(dr2*sin_th**2))

                  ff1(8,5)=(-sin_th*f2(2)*e2(2)+cos_th*(-one/dr2+rr2(2)*rr2(2)/dr2**3)+ &
                       one/dr1-rr1(2)*rr1(2)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(2)-e1(2))*(rr2(2)/(dr2**3*sin_th)-cos_th*f2(2)/(dr2*sin_th**2))

                  ff1(8,6)=(-sin_th*f2(3)*e2(2)+cos_th*(rr2(2)*rr2(3)/dr2**3)- &
                       rr1(2)*rr1(3)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(2)-e1(2))*(rr2(3)/(dr2**3*sin_th)-cos_th*f2(3)/(dr2*sin_th**2))

                  ff1(8,7)=(-sin_th*f3(1)*e2(2)+cos_th*(-rr2(2)*rr2(1)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(2)-e1(2))*(-rr2(1)/(dr2**3*sin_th)-cos_th*f3(1)/(dr2*sin_th**2))

                  ff1(8,8)=(-sin_th*f3(2)*e2(2)+cos_th*(one/dr2-rr2(2)*rr2(2)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(2)-e1(2))*(-rr2(2)/(dr2**3*sin_th)-cos_th*f3(2)/(dr2*sin_th**2))

                  ff1(8,9)=(-sin_th*f3(3)*e2(2)+cos_th*(-rr2(2)*rr2(3)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(2)-e1(2))*(-rr2(3)/(dr2**3*sin_th)-cos_th*f3(3)/(dr2*sin_th**2))
                  !........................................................
                  ff1(9,1:3)=ff1(1:3,9)

                  ff1(9,4)=(-sin_th*f2(1)*e2(3)+cos_th*(rr2(3)*rr2(1)/dr2**3)- &
                       rr1(3)*rr1(1)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(3)-e1(3))*(rr2(1)/(dr2**3*sin_th)-cos_th*f2(1)/(dr2*sin_th**2))

                  ff1(9,5)=(-sin_th*f2(2)*e2(3)+cos_th*(rr2(3)*rr2(2)/dr2**3)- &
                       rr1(3)*rr1(2)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(3)-e1(3))*(rr2(2)/(dr2**3*sin_th)-cos_th*f2(2)/(dr2*sin_th**2))

                  ff1(9,6)=(-sin_th*f2(3)*e2(3)+cos_th*(-one/dr2+rr2(3)*rr2(3)/dr2**3)+ &
                       one/dr1-rr1(3)*rr1(3)/dr1**3)/(dr2*sin_th)+ &
                       (cos_th*e2(3)-e1(3))*(rr2(3)/(dr2**3*sin_th)-cos_th*f2(3)/(dr2*sin_th**2))

                  ff1(9,7)=(-sin_th*f3(1)*e2(3)+cos_th*(-rr2(3)*rr2(1)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(3)-e1(3))*(-rr2(1)/(dr2**3*sin_th)-cos_th*f3(1)/(dr2*sin_th**2))

                  ff1(9,8)=(-sin_th*f3(2)*e2(3)+cos_th*(-rr2(3)*rr2(2)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(3)-e1(3))*(-rr2(2)/(dr2**3*sin_th)-cos_th*f3(2)/(dr2*sin_th**2))

                  ff1(9,9)=(-sin_th*f3(3)*e2(3)+cos_th*(one/dr2-rr2(3)*rr2(3)/dr2**3))/(dr2*sin_th)+ &
                       (cos_th*e2(3)-e1(3))*(-rr2(3)/(dr2**3*sin_th)-cos_th*f3(3)/(dr2*sin_th**2))
                  !........................................................
                  ff1(4,1:3)=ff1(1:3,4)

                  ff1(4,4)=((-e1(1)+e2(1)*cos_th+dr2*sin_th*f2(1))*e1(1)+ &
                       (dr1-dr2*cos_th)*(-one/dr1+rr1(1)*rr1(1)/dr1**3)+ &
                       (-e2(1)+e1(1)*cos_th+dr1*sin_th*f2(1))*e2(1)+ &
                       (dr2-dr1*cos_th)*(-one/dr2+rr2(1)*rr2(1)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(1)+(dr2-dr1*cos_th)*e2(1))* &
                       (rr1(1)/(dr1**3*dr2*sin_th)+rr2(1)/(dr1*dr2**3*sin_th)-cos_th*f2(1)/(dr1*dr2*sin_th**2))

                  ff1(4,5)=((-e1(2)+e2(2)*cos_th+dr2*sin_th*f2(2))*e1(1)+ &
                       (dr1-dr2*cos_th)*(rr1(1)*rr1(2)/dr1**3)+ &
                       (-e2(2)+e1(2)*cos_th+dr1*sin_th*f2(2))*e2(1)+ &
                       (dr2-dr1*cos_th)*(rr2(1)*rr2(2)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(1)+(dr2-dr1*cos_th)*e2(1))* &
                       (rr1(2)/(dr1**3*dr2*sin_th)+rr2(2)/(dr1*dr2**3*sin_th)-cos_th*f2(2)/(dr1*dr2*sin_th**2))

                  ff1(4,6)=((-e1(3)+e2(3)*cos_th+dr2*sin_th*f2(3))*e1(1)+ &
                       (dr1-dr2*cos_th)*(rr1(1)*rr1(3)/dr1**3)+ &
                       (-e2(3)+e1(3)*cos_th+dr1*sin_th*f2(3))*e2(1)+ &
                       (dr2-dr1*cos_th)*(rr2(1)*rr2(3)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(1)+(dr2-dr1*cos_th)*e2(1))* &
                       (rr1(3)/(dr1**3*dr2*sin_th)+rr2(3)/(dr1*dr2**3*sin_th)-cos_th*f2(3)/(dr1*dr2*sin_th**2))

                  ff1(4,7:9)=ff1(7:9,4)
                  !........................................................
                  ff1(5,1:3)=ff1(1:3,5)

                  ff1(5,4)=((-e1(1)+e2(1)*cos_th+dr2*sin_th*f2(1))*e1(2)+ &
                       (dr1-dr2*cos_th)*(rr1(2)*rr1(1)/dr1**3)+ &
                       (-e2(1)+e1(1)*cos_th+dr1*sin_th*f2(1))*e2(2)+ &
                       (dr2-dr1*cos_th)*(rr2(2)*rr2(1)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(2)+(dr2-dr1*cos_th)*e2(2))* &
                       (rr1(1)/(dr1**3*dr2*sin_th)+rr2(1)/(dr1*dr2**3*sin_th)-cos_th*f2(1)/(dr1*dr2*sin_th**2))

                  ff1(5,5)=((-e1(2)+e2(2)*cos_th+dr2*sin_th*f2(2))*e1(2)+ &
                       (dr1-dr2*cos_th)*(-one/dr1+rr1(2)*rr1(2)/dr1**3)+ &
                       (-e2(2)+e1(2)*cos_th+dr1*sin_th*f2(2))*e2(2)+ &
                       (dr2-dr1*cos_th)*(-one/dr2+rr2(2)*rr2(2)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(2)+(dr2-dr1*cos_th)*e2(2))* &
                       (rr1(2)/(dr1**3*dr2*sin_th)+rr2(2)/(dr1*dr2**3*sin_th)-cos_th*f2(2)/(dr1*dr2*sin_th**2))

                  ff1(5,6)=((-e1(3)+e2(3)*cos_th+dr2*sin_th*f2(3))*e1(2)+ &
                       (dr1-dr2*cos_th)*(rr1(2)*rr1(3)/dr1**3)+ &
                       (-e2(3)+e1(3)*cos_th+dr1*sin_th*f2(3))*e2(2)+ &
                       (dr2-dr1*cos_th)*(rr2(2)*rr2(3)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(2)+(dr2-dr1*cos_th)*e2(2))* &
                       (rr1(3)/(dr1**3*dr2*sin_th)+rr2(3)/(dr1*dr2**3*sin_th)-cos_th*f2(3)/(dr1*dr2*sin_th**2))

                  ff1(5,7:9)=ff1(7:9,5)
                  !........................................................
                  ff1(6,1:3)=ff1(1:3,6)

                  ff1(6,4)=((-e1(1)+e2(1)*cos_th+dr2*sin_th*f2(1))*e1(3)+ &
                       (dr1-dr2*cos_th)*(rr1(3)*rr1(1)/dr1**3)+ &
                       (-e2(1)+e1(1)*cos_th+dr1*sin_th*f2(1))*e2(3)+ &
                       (dr2-dr1*cos_th)*(rr2(3)*rr2(1)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(3)+(dr2-dr1*cos_th)*e2(3))* &
                       (rr1(1)/(dr1**3*dr2*sin_th)+rr2(1)/(dr1*dr2**3*sin_th)-cos_th*f2(1)/(dr1*dr2*sin_th**2))

                  ff1(6,5)=((-e1(2)+e2(2)*cos_th+dr2*sin_th*f2(2))*e1(3)+ &
                       (dr1-dr2*cos_th)*(rr1(3)*rr1(2)/dr1**3)+ &
                       (-e2(2)+e1(2)*cos_th+dr1*sin_th*f2(2))*e2(3)+ &
                       (dr2-dr1*cos_th)*(rr2(3)*rr2(2)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(3)+(dr2-dr1*cos_th)*e2(3))* &
                       (rr1(2)/(dr1**3*dr2*sin_th)+rr2(2)/(dr1*dr2**3*sin_th)-cos_th*f2(2)/(dr1*dr2*sin_th**2))

                  ff1(6,6)=((-e1(3)+e2(3)*cos_th+dr2*sin_th*f2(3))*e1(3)+ &
                       (dr1-dr2*cos_th)*(-one/dr1+rr1(3)*rr1(3)/dr1**3)+ &
                       (-e2(3)+e1(3)*cos_th+dr1*sin_th*f2(3))*e2(3)+ &
                       (dr2-dr1*cos_th)*(-one/dr2+rr2(3)*rr2(3)/dr2**3))/(dr1*dr2*sin_th)+ &
                       ((dr1-dr2*cos_th)*e1(3)+(dr2-dr1*cos_th)*e2(3))* &
                       (rr1(3)/(dr1**3*dr2*sin_th)+rr2(3)/(dr1*dr2**3*sin_th)-cos_th*f2(3)/(dr1*dr2*sin_th**2))


                  ff1(6,7:9)=ff1(7:9,6)
                  !........................................................

                  kk1=3*(ia_2-1)
                  kk2=3*(ia_1-1)
                  kk3=3*(ia_3-1)
                  do l=1,3
                     do m=1,9
                        if(m <= 3) then
                           fstore=f1(m)
                           m1=kk1+m
                        else if(m > 3 .and. m <= 6) then
                           m1=kk2+(m-3)
                           fstore=f2(m-3)
                        else if(m > 6 .and. m <= 9) then
                           m1=kk3+(m-6)
                           fstore=f3(m-6)
                        end if
                        l1=kk1+l; l2=kk2+l; l3=kk3+l
                        if(l1<epe_kl*3-3 .and. m1<epe_kl*3-3) hesse_cartesian_dummy(l1,m1)= &
                             hesse_cartesian_dummy(l1,m1)+d2E_dt2*fstore*f1(l)+dE_dt*ff1(l,m)
 !!!                    if(l1<epe_kl*3 .and. m1<epe_kl*3.and.l1.eq.3.and.m1.eq.3) print*,'3b 1',&
 !!!                       l,m,d2E_dt2*fstore*f1(l),dE_dt*ff1(l,m)
                        if(l2<epe_kl*3-3 .and. m1<epe_kl*3-3) hesse_cartesian_dummy(l2,m1)= &
                             hesse_cartesian_dummy(l2,m1)+d2E_dt2*fstore*f2(l)+dE_dt*ff1(l+3,m)
 !!!                    if(l2<epe_kl*3 .and. m1<epe_kl*3.and.l2.eq.3.and.m1.eq.3) print*,'3b 2',&
 !!!                         l,m,d2E_dt2*fstore*f2(l)+dE_dt,ff1(l+3,m)
                        if(l3<epe_kl*3-3 .and. m1<epe_kl*3-3) hesse_cartesian_dummy(l3,m1)= &
                             hesse_cartesian_dummy(l3,m1)+d2E_dt2*fstore*f3(l)+dE_dt*ff1(l+6,m)
 !!!                    if(l3<epe_kl*3 .and. m1<epe_kl*3.and.l3.eq.3.and.m1.eq.3) print*,'3b 3',&
 !!!                         l,m,d2E_dt2*fstore*f3(l),dE_dt*ff1(l+6,m)

                     end do
                  end do

               end do thrd
            end do scnd
         end do fst
      end if
 !    print*,'hesse_cartesian_dummy'
 !    do l=1,6
 !    print*, hesse_cartesian_dummy(l,:)
 !    enddo
    end subroutine add_epe_contrib

    subroutine add_dummies(hcart,hdumm)
      use opt_data_module, only: n_atoms, n_dummy, atom
      implicit none
      real(r8_kind), intent(in)  :: hcart(:,:) ! 3NA^2
      real(r8_kind), intent(out) :: hdumm(:,:) ! (3NA+3NX)^2
      ! *** end of interface ***

      integer(kind=i4_kind):: i,j,k,l

      i = n_atoms+n_dummy
      j = n_atoms
      ASSERT(size(hcart,1)==3*j)
      ASSERT(size(hcart,2)==3*j)
      ASSERT(size(hdumm,1)==3*i)
      ASSERT(size(hdumm,2)==3*i)

      ! introduce zero entries at the positions which belong to
      hdumm = 0.0_r8_kind

      k=0
      do i=1,n_atoms+n_dummy
         if(atom(i)%dummy) CYCLE
         k=k+1
         l=0
         do j=1,n_atoms+n_dummy
            if(atom(j)%dummy) CYCLE
            l=l+1
            hdumm(i*3-3+1:i*3-3+3,j*3-3+1:j*3-3+3) = &
            hcart(k*3-3+1:k*3-3+3,l*3-3+1:l*3-3+3)
         enddo
      enddo
    end subroutine add_dummies

    subroutine add_dervs(dervs,hdumm)
      use opt_data_module, only: n_atoms, n_dummy, atom
      implicit none
      real(r8_kind), intent(in)    :: dervs(:,:,:,:) !
      real(r8_kind), intent(inout) :: hdumm(:,:) ! (3NA+3NX)^2
      ! *** end of interface ***

      integer(kind=i4_kind):: i,j

      i = n_atoms+n_dummy
      ASSERT(size(hdumm,1)==3*i)
      ASSERT(size(hdumm,2)==3*i)

      do i=1,n_atoms+n_dummy
         if(atom(i)%dummy) CYCLE
         do j=1,n_atoms+n_dummy
            if(atom(j)%dummy) CYCLE

            hdumm(i*3-3+1:i*3-3+3,j*3-3+1:j*3-3+3) = &
            hdumm(i*3-3+1:i*3-3+3,j*3-3+1:j*3-3+3) + dervs(i,j,:3,:3)
         enddo
      enddo
    end subroutine add_dervs
  end subroutine hesse_prepare_internal

  !*************************************************************

  subroutine hesse_cart_to_internal(hesse_cartesian_dummy,hesse)
    ! Purpose: read in a Hesse matrix in cartesian coordinates
    !          if it exist
    !         (usually from the routines in the
    !          'valence_coord_module') and transform it to
    !         internal (Z-matrix) coordinates:
    !           hesse = bmat**(-1)*hesse_cartesian*bmat**(-1)T
    !         where bmat**(-1)T is the transpose of the inverse
    !         B-Matrix.
    !
    ! This has to work for the different sets of coordinates in
    ! the following way:
    ! 1. zmat_coordinates (symmetry-reduced)
    !   Transform the hessematrix from cartesian coordinates to
    !   primitive (non- symmetryreduced) Z-Matrixcoordinates :
    !   hesse_prim(n_primitive,n_primitive). In the next step
    !   reduce it with:
    !   hesse = reduc_mat * hesse_prim * expand_mat
    !
    ! 2. delocalized coordinates (based on non-symmetryreduced
    !    Z-Matrixcoordinates)
    !   Transform hesse_cartesian directly to delocalized
    !   internal coordinates using
    !   hesse = bmat**(-1) * hesse_cartesian * bmat**(-1)T
    !
    !  The two branches are coded in two different subroutines
    !  for greater convenience when extending or changing the
    !  code. The intermingling of the two transformation algorithms
    !  is in my view a serious obstacle to modularity.FN 3/98
    !
    use opt_data_module, only: zmat_coordinates, delocalized_coordinates, cart_coordinates
   !use coordinates_module, only: q_old ! and more below
   !use math_module, only: print_matrix,small,pi
   !use coordinates_module, only: bmat_inv !,expand_mat,reduc_mat
    implicit none
    !----- Declaration of formal paramaters -------------------
    real(kind=r8_kind), intent(in)  :: hesse_cartesian_dummy(:,:)
    real(kind=r8_kind), intent(out) :: hesse(:,:)
    ! *** end of interface ***

    if (zmat_coordinates) then

       call zmat_transform(hesse_cartesian_dummy,hesse)

    elseif (delocalized_coordinates) then

       call deloc_transform(hesse_cartesian_dummy,hesse)

    else if(cart_coordinates) then

       call cart_transform(hesse_cartesian_dummy,hesse)

    else
       ABORT('dont know what to do!')
    endif

  contains ! of hesse_cart_to_internal

    subroutine cart_transform(hesse_cartesian_dummy,hesse)
      implicit none
      real(kind=r8_kind), intent(in)  :: hesse_cartesian_dummy(:,:)
      real(kind=r8_kind), intent(out) :: hesse(:,:)
      ! *** end of interface ***

      hesse = hesse_cartesian_dummy
    end subroutine cart_transform

    subroutine zmat_transform(hesse_cartesian_dummy,hesse)
      ! Purpose: see above
      ! ----------------------------------------------------
      use opt_data_module, only: q, n_primitive, n_internal, n_atoms, n_dummy
      use opt_data_module, only: free_format, zmat_format
      use gradient_module, only: grad_cartes
      use coordinates_module, only: q_old ! and more below
      use coordinates_module, only: bmat_inv,expand_mat,reduc_mat
      use math_module, only: small
!     use debug
      implicit none
      real(kind=r8_kind), intent(in)  :: hesse_cartesian_dummy(:,:)
      real(kind=r8_kind), intent(out) :: hesse(:,:)
      ! *** end of interface ***

      integer(kind=i4_kind)           :: alloc_stat,i,j
      real(kind=r8_kind),allocatable  :: hesse_prim(:,:),help1(:,:),help2(:,:)
      real(kind=r8_kind), allocatable :: xrr(:,:,:)
      ! second derivatives d^2 Internal / d Cartesian^2 = d Bmat(:,Cartesian) / d Cartesian


      allocate(help1(n_primitive,3*(n_atoms+n_dummy)),&
           hesse_prim(n_primitive,n_primitive),&
           help2(n_primitive,n_internal),&
           STAT=alloc_stat)
      ASSERT(alloc_stat.eq.0)
!     call show('Bmat',bmat_inv)
!     call show('Hcart',hesse_cartesian_dummy)

      ! internal cartes dervs
      do i=1,n_primitive
         do j=1,3*(n_atoms+n_dummy)
            help1(i,j) = sum(hesse_cartesian_dummy(j,:)*bmat_inv(:,i))
         enddo
      enddo

      do i=1,n_primitive
         do j=1,n_primitive
            hesse_prim(i,j) = sum(bmat_inv(:,i)*help1(j,:))
         enddo
      enddo


!     call show('Hint',hesse_prim)
      ! in constrast to the above lines the following is not
      ! intendend to gain optimal performance of the program.
      ! Whoever thinks that this may become an issue is kindly
      ! invited to change these lines as well as all the other
      ! places in the program where clarity was retained at the
      ! expense of performance.
      if (zmat_format) then

         help2 = matmul(hesse_prim,expand_mat)
         hesse = matmul(reduc_mat,help2)

         ! outside of extrema derivatives of bmat are needed:
         allocate(xrr(size(bmat_inv,1),size(q),size(q)))

         ! This subroutine calculates derivatives of the B^-1 matrix
         ! with respect to internal coordinates ones.
         ! It is necessary to convert the hessian matrix into internal coordinates.
         ASSERT(associated(q))
         call calc_xrr(q,xrr) ! calculates ``xrr''
         ASSERT(all(q==q_old))

         ! outside of extrema non-zero (cartesian) forces contribute to (internal) hessian:
         hesse = hesse + xrr_grad_cartes(grad_cartes,xrr)

         deallocate(xrr)

      elseif (free_format) then
         hesse = hesse_prim
      endif

      ! check if there are coolumns which consist only of zeros
#if 1
      do i=1,n_internal
         if(sum(abs(hesse(:,i)))<small) print*,i,n_internal
         if(sum(abs(hesse(:,i)))<small) call error_handler(&
              "cart_to_internal: something wrong with your cartesian hesse matrix")
      end do
#endif

      deallocate(help1,help2,hesse_prim,STAT=alloc_stat)
      ASSERT(alloc_stat.eq.0)
    end subroutine zmat_transform

    function xrr_grad_cartes(grad_cartes,xrr)
!     use gradient_module, only: grad_cartes
      real(kind=r8_kind), intent(in) :: grad_cartes(:,:) ! (n_atoms,3)
      real(kind=r8_kind), intent(in) :: xrr(:,:,:)       ! (?,?,?)
      real(kind=r8_kind) :: xrr_grad_cartes(size(xrr,2),size(xrr,2))
      ! *** end if interface ***

      real(kind=r8_kind) :: gradients(size(xrr,1))

      integer(kind=i4_kind):: i,j
      do i=1,3
        do j=1,size(grad_cartes,1)
          gradients(i+3*j-3)=grad_cartes(j,i)
        enddo
      enddo
      forall(i=1:size(xrr,2),j=1:size(xrr,2)) &
        xrr_grad_cartes(i,j)=sum(gradients*xrr(:,i,j))
    end function xrr_grad_cartes

    subroutine calc_xrr(q,xrr)
      ! This subroutine calculates derivatives of the B^-1 matrix
      ! with respect to internal coordinates ones.
      ! It is necessary to convert the hessian matrix into internal coordinates.
      use opt_data_module, only: OPT_STDOUT
      use coordinates_module, only: bmat, bmat_inv, reduc_mat, internal_to_cart
      implicit none
      real(r8_kind), intent(in)  :: q(:)
      real(r8_kind), intent(out) :: xrr(:,:,:)
      ! *** end of interface ***

      integer(kind=i4_kind) :: displaced
      real(r8_kind)         :: bmat_inv_d(size(bmat_inv,1),size(q),2)

      real(r8_kind), parameter          :: delta = 0.01_r8_kind
      real(r8_kind), dimension(size(q)) :: qo,qp,qm ! original, plus, minus

      write(OPT_STDOUT,*) "main_opt: computing derivatives of B^-1 wrt internal coordinates"
      DPRINT              "main_opt: computing derivatives of B^-1 wrt internal coordinates"

      ! save original values:
      qo = q ! or global q_old?

      do displaced=1,size(q)
        ! FIXME: dont we destroy the cartesian coordinates and
        !        the values of bmat and bmat_inv in global vars?
        !        Try to go in circles: original->plus delta->minus delta->original
        !        (qo->qp->qm->qo)

        ! prepare two geometries...
        qp = qo
        qm = qo
        qp(displaced) = qp(displaced) + delta
        qm(displaced) = qm(displaced) - delta
        ! ... and move the system updating the B-matrices:

        ! o->p: update cartesian, update B-matrices:
        call internal_to_cart(qo,qp,bmat,bmat_inv)

        ! reduced version:
        bmat_inv_d(:,:,1) = matmul(bmat_inv,transpose(reduc_mat))

        ! p->m: update cartesian, update B-matrices:
        call internal_to_cart(qp,qm,bmat,bmat_inv)

        ! reduced version:
        bmat_inv_d(:,:,2) = matmul(bmat_inv,transpose(reduc_mat))

        ! finite differencing for the B^-1:
        xrr(:,:,displaced)=(bmat_inv_d(:,:,1)-bmat_inv_d(:,:,2)) / (2*delta)

        ! m->o: restore cartesian coordinates and B-matrices:
        call internal_to_cart(qm,qo,bmat,bmat_inv)
      enddo
    end subroutine calc_xrr

    subroutine deloc_transform(hesse_cartesian_dummy,hesse)
      ! Purpose: see above
      ! ----------------------------------------------------
      use opt_data_module, only: n_internal, n_atoms, n_dummy
      use coordinates_module, only: bmat_inv
      implicit none
      real(kind=r8_kind), intent(in)  :: hesse_cartesian_dummy(:,:)
      real(kind=r8_kind), intent(out) :: hesse(:,:)
      ! *** end of interface ***

      integer(kind=i4_kind)          :: alloc_stat,i,j
      real(kind=r8_kind),allocatable :: help(:,:)
!     real(kind=r8_kind)             :: mini,maxi

      allocate( help(n_internal,3*(n_atoms+n_dummy)), STAT=alloc_stat)
      if (alloc_stat/=0) call error_handler &
           ("deloc_transform: allocation (1) failed")

      do i=1,n_internal
         do j=1,3*(n_atoms+n_dummy)
            help(i,j) = sum(hesse_cartesian_dummy(j,:)*bmat_inv(:,i))
         enddo
      enddo
      do i=1,n_internal
         do j=1,n_internal
            hesse(i,j) = sum(bmat_inv(:,i)*help(j,:))
         enddo
         !mini=minval(hesse(i,:))
         !maxi=maxval(hesse(i,:))
         !if (abs(mini)<=small.and.abs(maxi)<small) then
         !   hesse(i,i)=one
         !endif
         !if(sum(abs(hesse(i,:)))<small) call error_handler(&
         !     "cart_to_internal: something wrong with your cartesian hesse matrix")
      enddo
!!$   print*,' This is the Hessian in delocalized coordinates :'
!!$   call print_matrix(hesse,n_internal,n_internal,10_i4_kind)
!!$   print*,' -------------------------- hesse -------------------------'
      deallocate(help,STAT=alloc_stat)
      if (alloc_stat/=0 ) call error_handler&
           (" deloc_transform: deallocation (1) failed")
    end subroutine deloc_transform
 end subroutine hesse_cart_to_internal

  !*************************************************************

subroutine hesse_internal_to_cart(int_flag)
  ! Prupose: after 'exact' caluation of the Hessian, it is transformed
  !          to cartesian coordinaten for greater flexibility in
  !          re-using the precious thing.
  !          The transformation is done as follows:
  !          1. Expand the matrix 'hesse' to the full coordinate
  !          space using
  !          hesse_prim = expand_mat * hesse * reduc_mat
  !          2. Then the cartesian Hessian is obtaied as usual:
  !          hesse_cartesian = B**T * Hesse_prim * B
  !          After transformation hesse_cartesian is stored in
  !          the file 'hesse_cartesian.dat' which is copied into
  !          the optimizer-directory by the script 'optimizer'.
  ! -----------------------------------------------------------------
  use math_module
  use opt_data_module
  use coordinates_module
  use filename_module, only: inpfile
  real(kind=r8_kind), allocatable :: hesse_cartesian_dummy(:,:),hesse_cartesian(:,:),&
       help_mat1(:,:),hesse_prim(:,:),help_mat2(:,:), reduc_help(:,:)
  integer(kind=i4_kind)           :: alloc_stat,i,j,n_sym,i_counter,j_counter,io_hesse_cart,io_hesse_cart_store
  integer(kind=i4_kind), intent(in), optional :: int_flag
! real(kind=r8_kind)              :: mini,maxi
  character(len=4) char_int_flag

  allocate(hesse_cartesian(3*n_atoms,3*n_atoms),&
       hesse_cartesian_dummy(3*(n_atoms+n_dummy),3*(n_atoms+n_dummy)),&
       help_mat1(n_primitive,3*(n_atoms+n_dummy)),&
       help_mat2(n_internal,n_primitive),&
       hesse_prim(n_primitive,n_primitive),&
       reduc_help(n_internal,n_primitive),&
       STAT=alloc_stat)
  if (alloc_stat/=0) call error_handler&
       ("hesse_internal_to_cart: allocation (1) failed")
  help_mat1=zero
  help_mat2=zero
  hesse_cartesian=zero
  hesse_prim=zero


  do i=1,n_internal
     n_sym = sum(abs(reduc_mat(i,:)))
     reduc_help(i,:)=reduc_mat(i,:)/n_sym
  enddo
  help_mat2=matmul(hesse,reduc_help)
  hesse_prim=matmul(transpose(reduc_help),help_mat2)

  !do i=1,n_primitive
  !   mini=minval(hesse_prim(:,i))
  !   maxi=maxval(hesse_prim(:,i))
  !   if (abs(mini)<small.and.abs(maxi)<small) hesse_prim(i,i)=one
  !enddo


  ASSERT(allocated(bmat))
  help_mat1 = matmul(hesse_prim,bmat)
  hesse_cartesian_dummy = matmul(transpose(bmat),help_mat1)

  ! now we elliminate dummy atoms from the hessian matrix
  i_counter=0
  do i=1,n_atoms+n_dummy
     if(.not.atom(i)%dummy) then
        i_counter=i_counter+1
        j_counter=0
        do j=1,n_atoms+n_dummy
           if(.not.atom(j)%dummy) then
              j_counter=j_counter+1
              hesse_cartesian(i_counter*3-3+1:i_counter*3-3+3,j_counter*3-3+1:j_counter*3-3+3) = &
                   hesse_cartesian_dummy(i*3-3+1:i*3-3+3,j*3-3+1:j*3-3+3)
           end if
        end do
     end if
  end do
  call print_matrix(hesse_cartesian,3*(n_atoms),3*(n_atoms),10)

  io_hesse_cart=openget_iounit(status='unknown',form='formatted',&
       FILE=trim(inpfile('hesse_cartesian.dat')))
  if(present(int_flag)) then
  write(char_int_flag,'(i4)')  int_flag
  io_hesse_cart_store=openget_iounit(status='unknown',form='formatted',&
       FILE=trim(opt_dir)//'/hesse_cartesian.'//adjustl(char_int_flag))
  endif
  if (print_debug) then
     write(OPT_STDOUT,*)"hesse_internal_to_cart: writing HESSE_CARTESIAN to file"
  endif
  write(io_hesse_cart,*,IOSTAT=alloc_stat)hesse_cartesian
  if(present(int_flag)) write(io_hesse_cart_store,*)hesse_cartesian
  if (alloc_stat/=0) call error_handler&
       ("hesse_internal_to_cart: writing HESSE_CARTESIAN to file failed")
!  close(io_hesse)
  call returnclose_iounit(io_hesse_cart)
  if(present(int_flag)) call returnclose_iounit(io_hesse_cart_store)
  deallocate(hesse_cartesian_dummy,hesse_cartesian,help_mat1,help_mat2,hesse_prim,&
       reduc_help,STAT=alloc_stat)
  if (alloc_stat/=0) call error_handler&
       ("hesse_internal_to_cart: deallocation (1) failed")
end subroutine hesse_internal_to_cart

  !*************************************************************

  subroutine tsscan_hesse_cartesian(hesse_cartesian)
    !
    ! Does what?
    !
    use opt_data_module, only: n_atoms,n_dummy,atom,xyz,xyz_reactant,xyz_product
    use opt_data_module, only: rpmix,epe_forces
    use slspar_module
    implicit none
    real(r8_kind), intent(out) :: hesse_cartesian(:,:) ! (3NA,3NA)
    ! *** end of interface ***

    real(kind=r8_kind) :: fdev,zero=0.0_r8_kind,pmix=0.2_r8_kind,r2_r,r2_p
    integer(i4_kind)   :: i,j,i_c,j_c

    DPRINT 'tsscan_hesse_cartesian: calculating hesse_cartesian '
    !SSERT(size(hesse_cartesian,1)==3*n_atoms)
    !SSERT(size(hesse_cartesian,2)==3*n_atoms)

    pmix=rpmix
    hesse_cartesian=zero
     i_c=0
     do i=1,n_atoms+n_dummy
      if(atom(i)%dummy) cycle
      j_c=0
      do j = 1, n_atoms+n_dummy
        if(atom(j)%dummy) cycle
        if(i.eq.j) cycle

!       grad_cartes(i,:)=grad_cartes(i,:)+4*(xyz(:,i)-xyz(:,j)) &
!           *( sum( (xyz(:,i)-xyz(:,j))**2 ) &
!             -sum( (xyz_reactant(:,i)-xyz_reactant(:,j))**2 )/2 &
!             -sum( (xyz_product(:,i)-xyz_product(:,j))**2 )/2)

        r2_r=sum( (xyz_reactant(:,i)-xyz_reactant(:,j))**2 )
        r2_p=sum( (xyz_product(:,i)-xyz_product(:,j))**2 )
        fdev=( sum( (xyz(:,i)-xyz(:,j))**2 )-pmix*r2_r-(1-pmix)*r2_p)

        hesse_cartesian(i_c+1,i_c+1)=hesse_cartesian(i_c+1,i_c+1) &
         + 4*(fdev+2*(xyz(1,i)-xyz(1,j))**2 )/(r2_r+r2_p)**2
        hesse_cartesian(i_c+2,i_c+2)=hesse_cartesian(i_c+2,i_c+2) &
         + 4*(fdev+2*(xyz(2,i)-xyz(2,j))**2 )/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,i_c+3)=hesse_cartesian(i_c+3,i_c+3) &
         + 4*(fdev+2*(xyz(3,i)-xyz(3,j))**2 )/(r2_r+r2_p)**2

        hesse_cartesian(i_c+1,i_c+2)=hesse_cartesian(i_c+1,i_c+2) &
         + 8*(xyz(1,i)-xyz(1,j))*(xyz(2,i)-xyz(2,j))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+2,i_c+1)=hesse_cartesian(i_c+2,i_c+1) &
         + 8*(xyz(2,i)-xyz(2,j))*(xyz(1,i)-xyz(1,j))/(r2_r+r2_p)**2

        hesse_cartesian(i_c+1,i_c+3)=hesse_cartesian(i_c+1,i_c+3) &
         + 8*(xyz(1,i)-xyz(1,j))*(xyz(3,i)-xyz(3,j))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,i_c+1)=hesse_cartesian(i_c+3,i_c+1) &
         + 8*(xyz(3,i)-xyz(3,j))*(xyz(1,i)-xyz(1,j))/(r2_r+r2_p)**2

        hesse_cartesian(i_c+2,i_c+3)=hesse_cartesian(i_c+2,i_c+3) &
         + 8*(xyz(2,i)-xyz(2,j))*(xyz(3,i)-xyz(3,j))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,i_c+2)=hesse_cartesian(i_c+3,i_c+2) &
         + 8*(xyz(3,i)-xyz(3,j))*(xyz(2,i)-xyz(2,j))/(r2_r+r2_p)**2


        hesse_cartesian(i_c+1,j_c+1)=hesse_cartesian(i_c+1,j_c+1) &
         - 4*(fdev+2*(xyz(1,i)-xyz(1,j))**2 )/(r2_r+r2_p)**2
        hesse_cartesian(i_c+2,j_c+2)=hesse_cartesian(i_c+2,j_c+2) &
         - 4*(fdev+2*(xyz(2,i)-xyz(2,j))**2 )/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,j_c+3)=hesse_cartesian(i_c+3,j_c+3) &
         - 4*(fdev+2*(xyz(3,i)-xyz(3,j))**2 )/(r2_r+r2_p)**2

        hesse_cartesian(i_c+1,j_c+2)=hesse_cartesian(i_c+1,j_c+2) &
         - 8*(xyz(1,i)-xyz(1,j))*(xyz(2,i)-xyz(2,j))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+2,j_c+1)=hesse_cartesian(i_c+2,j_c+1) &
         - 8*(xyz(2,i)-xyz(2,j))*(xyz(1,i)-xyz(1,j))/(r2_r+r2_p)**2

        hesse_cartesian(i_c+1,j_c+3)=hesse_cartesian(i_c+1,j_c+3) &
         - 8*(xyz(1,i)-xyz(1,j))*(xyz(3,i)-xyz(3,j))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,j_c+1)=hesse_cartesian(i_c+3,j_c+1) &
         - 8*(xyz(3,i)-xyz(3,j))*(xyz(1,i)-xyz(1,j))/(r2_r+r2_p)**2

        hesse_cartesian(i_c+2,j_c+3)=hesse_cartesian(i_c+2,j_c+3) &
         - 8*(xyz(2,i)-xyz(2,j))*(xyz(3,i)-xyz(3,j))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,j_c+2)=hesse_cartesian(i_c+3,j_c+2) &
         - 8*(xyz(3,i)-xyz(3,j))*(xyz(2,i)-xyz(2,j))/(r2_r+r2_p)**2

       j_c=j_c+3
     enddo
     if(epe_forces) then
      do j = epe_kl,epe_nucen+epe_kl-1

!       grad_cartes(i,:)=grad_cartes(i,:)+4*(xyz(:,i)-xyz(:,j)) &
!           *( sum( (xyz(:,i)-xyz(:,j))**2 ) &
!             -sum( (xyz_reactant(:,i)-xyz_reactant(:,j))**2 )/2 &
!             -sum( (xyz_product(:,i)-xyz_product(:,j))**2 )/2)

        r2_r=sum( (xyz_reactant(:,i)-epe(j)%s)**2 )
        r2_p=sum( (xyz_product(:,i)-epe(j)%s)**2 )
        fdev=( sum( (xyz(:,i)-epe(j)%s)**2 )-pmix*r2_r-(1-pmix)*r2_p)

        hesse_cartesian(i_c+1,i_c+1)=hesse_cartesian(i_c+1,i_c+1) &
         + 4*(fdev+2*(xyz(1,i)-epe(j)%s(1))**2 )/(r2_r+r2_p)**2
        hesse_cartesian(i_c+2,i_c+2)=hesse_cartesian(i_c+2,i_c+2) &
         + 4*(fdev+2*(xyz(2,i)-epe(j)%s(2))**2 )/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,i_c+3)=hesse_cartesian(i_c+3,i_c+3) &
         + 4*(fdev+2*(xyz(3,i)-epe(j)%s(3))**2 )/(r2_r+r2_p)**2

        hesse_cartesian(i_c+1,i_c+2)=hesse_cartesian(i_c+1,i_c+2) &
         + 8*(xyz(1,i)-epe(j)%s(1))*(xyz(2,i)-epe(j)%s(2))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+2,i_c+1)=hesse_cartesian(i_c+2,i_c+1) &
         + 8*(xyz(2,i)-epe(j)%s(2))*(xyz(1,i)-epe(j)%s(1))/(r2_r+r2_p)**2

        hesse_cartesian(i_c+1,i_c+3)=hesse_cartesian(i_c+1,i_c+3) &
         + 8*(xyz(1,i)-epe(j)%s(1))*(xyz(3,i)-epe(j)%s(3))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,i_c+1)=hesse_cartesian(i_c+3,i_c+1) &
         + 8*(xyz(3,i)-epe(j)%s(3))*(xyz(1,i)-epe(j)%s(1))/(r2_r+r2_p)**2

        hesse_cartesian(i_c+2,i_c+3)=hesse_cartesian(i_c+2,i_c+3) &
         + 8*(xyz(2,i)-epe(j)%s(2))*(xyz(3,i)-epe(j)%s(3))/(r2_r+r2_p)**2
        hesse_cartesian(i_c+3,i_c+2)=hesse_cartesian(i_c+3,i_c+2) &
         + 8*(xyz(3,i)-epe(j)%s(3))*(xyz(2,i)-epe(j)%s(2))/(r2_r+r2_p)**2
     enddo
     endif
     i_c=i_c+3
    enddo
  end subroutine tsscan_hesse_cartesian


  subroutine hesse_calc(q0,g0,q1,g1,q_act,g_act,sym)
    ! performs a step by step calculation of the Hessian
    ! around the point 'q0' (read in from 'hesse.dat').
    !---------------------------------------------------------
    use opt_data_module
    use math_module
    use coortype_module
    use coordinates_module
    real(kind=r8_kind),intent(in)               :: q0(:),g0(:)
    real(kind=r8_kind),intent(in)               :: q_act(:),g_act(:)
    real(kind=r8_kind),optional,intent(in)      :: q1(:),g1(:)
    integer(kind=i4_kind),intent(in)            :: sym
    ! --- declaration of local variables ---------------------
    integer(kind=i4_kind)          :: i,j,alloc_stat
    real(kind=r8_kind),allocatable :: dq(:),dg(:),dq1(:),dg1(:)


    allocate(dg(n_internal),dq(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" hesse_calc: allocation (1) failed")
    dg = g0 - g_act
    dq = q0 - q_act
    ! If Z-Matrix coordinates are used, care has to be taken
    ! with diherdral angles
    do i=1,n_internal
       if (s_local(i)%typ == d_angle ) then
          if (dq(i) > pi ) dq(i) = -two*pi+dq(i)
          if (dq(i) < -pi ) dq(i) =  two*pi+dq(i)
       endif
    enddo
    if (present(q1).and.present(g1)) then
       allocate(dg1(n_internal),dq1(n_internal),STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            (" hesse_calc: allocation (2) failed")
       dg1 = g1 - g0
       dq1 = q1 - q0
       do i=1,n_internal
          if (s_local(i)%typ == d_angle ) then
             if (dq1(i) > pi ) dq1(i) = -two*pi+dq1(i)
             if (dq1(i) < -pi ) dq1(i) =  two*pi+dq1(i)
          endif
       enddo
    endif
    ! update all those lines in hesse wich correspond to
    ! internal variables of sym_type 'sym'.
    ! Attention has to be paid:
    ! - Different Coordinates of the same symmetry-type have to be
    !   decoupled
    !   => update only if both coordinates are unique OR
    !      if both are NOT unique, but never if one is
    !      unique and the other not unique.
    if (zmat_coordinates) then
       if (abs(dq(sym)) < small ) call error_handler&
            ("hesse_calc: delta q close to zero")
       do i=1,n_internal
          if (present(q1).and.present(g1)) then
             hesse(i,sym) = half*(dg(i)/dq(sym) + dg1(i)/dq1(sym))
          else
             hesse(i,sym) = dg(i)/dq(sym)
          endif
       enddo
    elseif (delocalized_coordinates) then
       do i=1,n_internal
          if (sym_type(i)/=sym) cycle
          do j=1,n_internal
             if (abs(dq(i))<=small) then
                write(OPT_STDOUT,*)" hesse_calc: delta q close to zero"
                stop 1
             endif
             if ((s(i)%unique.and.s(j)%unique).or.&
                  (.not.s(i)%unique.and..not.s(j)%unique)) then
                if (present(q1).and.present(g1)) then
                   hesse(j,i) = half*(dg(j)/dq(i) + dg1(j)/dq1(i))
                else
                   hesse(j,i) = dg(j)/dq(i)
                endif
             endif
          enddo
       enddo
    endif
    deallocate(dq,dg,STAT=alloc_stat)
    if (alloc_stat/=0) then
       write(OPT_STDOUT,*)" hesse_calc: deallocation (1) failed"
       stop 1
    endif
    if (present(q1).and.present(g1)) then
       deallocate(dq1,dg1,STAT=alloc_stat)
       if (alloc_stat/=0) then
          write(OPT_STDOUT,*)" hesse_calc: deallocation (2) failed"
          stop 1
       endif
    endif
  end subroutine hesse_calc
  !*************************************************************

  subroutine hesse_step(q0,sym,delta_q,back)
    ! Purpose: performs a step for the calculation of the hessian
    !          Since the the smal steps taken must be consistent
    !          with the symmetry of the molecule, variables of the
    !          the same symmetry type (->'sym_type') have to be
    !          elongated simultaneously by the the same amount.
    ! ------------------------------------------------------------
    use opt_data_module
    use coordinates_module
    real(kind=r8_kind),intent(in)    :: q0(:)
    integer(kind=i4_kind),intent(in) :: sym
    real(kind=r8_kind),intent(in)    :: delta_q
    logical,optional                 :: back
    ! --- declaration of local variables -------------------------
    logical                 :: local
    integer(kind=i4_kind)   :: i
    if (present(back)) then
       if (back) then
          local=.true.
       else
          local=.false.
       endif
    else
       local=.false.
    endif
    if (delocalized_coordinates) then ! right symmetry type has to be
       do i=1,n_internal ! found out since the internals are not symmetry_reduced
          if (sym_type(i)==sym) then
             if (local) then
                s_local(i)%value = q0(i) - delta_q
             else
                s_local(i)%value = q0(i) + delta_q
             endif
          else
             s_local(i)%value=q0(i)
          endif
       enddo
       q = matmul(umat_trans,s_local(:)%value) ! to delocalized coords
    elseif (zmat_coordinates) then
       s_local(:)%value = q0
       if (local) then
          s_local(sym)%value = q0(sym) - delta_q
       else
          s_local(sym)%value = q0(sym) + delta_q
       endif
       q = s_local(:)%value
    endif

  end subroutine hesse_step

  !*************************************************************

  subroutine symmetrize()

    use coortype_module
    use opt_data_module
    use coordinates_module
    use math_module
    use valence_coord_module
    ! Purpose: build symmetry-consistent linear combinations
    !          of eigenvectors belonging to degenerate eigevalues.
    !          This is necessary since the symmetry-equivalent
    !          coordinates have been completely decoupled when
    !          setting up the Hessian leading to also completely
    !          decoupled eigenmodes. These in turn are not
    !          symmetry-consistent any more but since any linear
    !          combination of degenerate eigenvectors is an
    !          eigenvector again, we can symmetrize them.
    ! ---------------------------------------------------------
    integer(kind=i4_kind),allocatable    :: equiv(:)
    real(kind=r8_kind),allocatable       :: help_vec(:)
    integer(kind=i4_kind)                :: alloc_stat,i_eig,j,i,&
         eq,symm
    real(kind=r8_kind)                   :: sign_first,sign_add,eig

    allocate(equiv(n_internal),help_vec(n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("symmetrize : allocation (1) failed")
    equiv=0_i4_kind
    help_vec=zero

    do i_eig=1,n_internal
       if(equiv(i_eig)/=0) cycle
       eig = hesse_eigval(i_eig)
       do j=1,n_internal
          if (abs(hesse_eigval(j) - eig)<=small ) then
             equiv(j) = i_eig
          endif
       enddo
    enddo

    ! add those eigenvecs belonging to degenerate eigenvalues
    do i_eig=1,n_internal
       eq=equiv(i_eig)
       help_vec=zero
       symm=-99
       ! pick the sign of first non-zero component,
       ! also save the symmetry-type
       loop_1: do j=1,n_internal
          if (equiv(j)==eq) then
             do i=1,n_internal
                if (hesse_eigvec(i,j)==zero) cycle
                if (hesse_eigvec(i,j)<zero) then
                   sign_first = -one
                else
                   sign_first=one
                endif
                symm=sym_type(i)
             enddo
             exit loop_1
          endif
       enddo loop_1
       do j=1,n_internal ! Loop through eigenvectors and skip those which do not belong
          if (equiv(j)/=eq) cycle ! to the same eigenvalue
          ! now look at the sign of the first component of
          ! sym_type(i) == symm.
          loop_2: do i=1,n_internal
             if (hesse_eigvec(i,j) == zero ) cycle
             if (sym_type(i) == symm) then
                if (hesse_eigvec(i,j) < zero ) then
                   sign_add = -one
                else
                   sign_add = one
                endif
                exit loop_2
             endif
          enddo loop_2
          if (sign_add == sign_first) then
             sign_add = one
          else
             sign_add = -one
          endif
          help_vec=help_vec+sign_add*hesse_eigvec(:,j)
       enddo
       help_vec = help_vec/abs_value(help_vec)
       do j=1,n_internal
          if (equiv(j)/=eq) cycle
          hesse_eigvec(:,j) = help_vec
       enddo
    enddo

    deallocate(equiv,help_vec,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ("symmetrize: deallocation (1) failed")

  end subroutine symmetrize

  !*************************************************************
  subroutine delta_coordinate(q_1,q_2,delta_q)
    use coortype_module
    use opt_data_module
    use coordinates_module
    use math_module
    ! Purpose: build the difference
    !          delta_q = q_2 - q_1
    !          for zmat_coordinates and for delocalized coordinates
    !          with special consideration of dihedral angles.
    !          1. ZMAT_COORDINATES
    !          The difference can be calculated directly with only
    !          one if-branch for dihedral angles, that will
    !          add/subtract 2*pi.
    !          2. DELOCALIZED COORDINATES
    !          Here, the input information q_1 and q_2 is expected
    !          in primitive coordinates, whereas the output information
    !          will be given in delocalized coordinates. This means
    !          that on entry to this routine, delta_q HAS TO HAVE a
    !          different dimension than q_1 and q_2.
    !          The difference in primtives is first build in the same
    !          way as for ZMAT_COORDINATES and then the difference
    !          is trandformed to delocalized coordiantes.
    !
    ! ---------------------------------------------------------
    real(kind=r8_kind),intent(in)         :: q_1(:),q_2(:)
    real(kind=r8_kind)                    :: delta_q(:)
    integer(kind=i4_kind)                 :: dimen
    ! ------ help variables for delocalized coords ---
    real(kind=r8_kind),allocatable        :: delta_qprim(:)
    integer(kind=i4_kind)                 :: alloc_stat,i
    ! ---------------------------------------------------------

    dimen=ubound(q_1,1)

    if ( zmat_coordinates ) then
       ! use the variable 's' directly instead of 's_local'
       do i=1,dimen
          delta_q(i) = q_2(i) - q_1(i)
          if (s(i)%typ == d_angle ) then
             if (delta_q(i) > pi ) delta_q(i) = -two*pi+delta_q(i)
             if (delta_q(i) < -pi ) delta_q(i) =  two*pi+delta_q(i)
          endif
       enddo

    elseif (delocalized_coordinates ) then
       ! here use the variable 's_prim' for typ identification
       allocate(delta_qprim(n_primitive),STAT=alloc_stat)
       if (alloc_stat /=0 ) call error_handler &
            ("delta_coordinate: allocation (1) failed")
       delta_qprim=zero

       do i=1,n_primitive
          delta_qprim(i) = q_2(i) - q_1(i)
          if (s_prim(i)%typ == d_angle ) then
             if (delta_qprim(i) > pi ) then
                delta_qprim(i) = -two*pi+delta_qprim(i)
             endif
             if (delta_qprim(i) < -pi ) then
                delta_qprim(i) =  two*pi+delta_qprim(i)
             endif
          endif
       enddo
       delta_q = matmul(umat_trans,delta_qprim)
       deallocate(delta_qprim,STAT=alloc_stat)
       if (alloc_stat/=0 ) call error_handler &
            ("delta_coordinate: deallocation (1) failed ")

    elseif(cart_coordinates) then
          delta_q = q_2 - q_1
    endif
  end subroutine delta_coordinate

  !*************************************************************

  subroutine write_cart_hess(hess)
    use filename_module, only: inpfile
    use atom_data_module, only: nuc_mass
    use frequency_module, only: freq,freq_print,dipole_der
    use opt_data_module, only: n_atoms,charge,OPT_STDOUT,n_dummy,dummy_list
    implicit none
    real(r8_kind)   ,intent(in)   :: hess(:,:)
    ! *** end of interface ***

    real(r8_kind), parameter     :: mp = 1836.15267261_r8_kind ! Mp/Me
    real(r8_kind)                :: mass(3*n_atoms)
    real(r8_kind)                :: fre(3*n_atoms)
    real(r8_kind)                :: mode(3*n_atoms,3*n_atoms)
    integer(i4_kind)             :: ihess,i,j
    real(kind=r8_kind)           :: dipder(3,3*n_atoms)
    real(r8_kind)                :: intensities(3*n_atoms)
    logical                      :: error

    ihess=get_iounit()
    open(ihess, file=trim(inpfile('HESS.cart')), err=100)

    ! dump cartesian hessian into a file:
    do i=1,3*n_atoms
       write(ihess,'(3F20.10)') ( hess(i,j), j=1,3*n_atoms )
    enddo

    call returnclose_iounit(ihess)

    ! build (diagonal) mass matrix:
    j=0
    do i=1,n_atoms+n_dummy
       if(.not.dummy_list(i)) then
          j=j+1
          mass(3*j-2:3*j) = nuc_mass( NINT(charge(i)) ) * (mp/nuc_mass(1))
       end if
    enddo

    ! compute dipole derivatives wrt cartesian coordinates
    call dipole_der(3*n_atoms,dipder,error)
    if (error) then ! dipole switch is off, so dipole.dat is not present for dipder
     print*,'hesse_module: write_cart_hess: dipole is not present'
      ! compute frequencies:
      call freq(hess,mass,fre,mode)

      ! for printing the frequencies in the flepo file and stdout
      call freq_print(OPT_STDOUT,mass,fre,mode)
      call freq_print(         6,mass,fre,mode) ! FIXME: stdout
    else
     print*,'hesse_module: write_cart_hess: intensity'
      call freq(hess,mass,fre,mode,dipder,intensities)
      call freq_print(OPT_STDOUT,mass,fre,mode,intensities)
      call freq_print(         6,mass,fre,mode,intensities) ! FIXME: stdout
    endif
    return
100 call error_handler('HESSE_MODULE: write_cart_hess: error - file "HESS.cart"')
  end subroutine write_cart_hess

  !--------------- End of module -------------------------------------
end module hesse_module
