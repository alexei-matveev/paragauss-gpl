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
module  gradient_module
  !---------------------------------------------------------------
  !  Purpose: contains gradient  in cartesian and internal
  !           coordinates plus routines to transform, print
  !           out, allocate and deallocate etc. blablah
  !
  !  Module called by: ...
  !
  !  References: ...
  !
  !  Author: FN
  !  Date: 2/98
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
#include "def.h"
  use type_module ! type specification parameters
   use opt_data_module
#ifdef WITH_EFP
  use qmmm_interface_module, only: efp
  use efp_module, only: n_efp, efp_fixed, qm_fixed
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------
  real(kind=r8_kind),allocatable,public,target    :: grad_cartes(:,:)
  real(kind=r8_kind),allocatable,public,target    :: dervs_cartes(:,:,:,:)
  real (r8_kind), allocatable, public, target :: grad_intern(:)
  real(kind=r8_kind),allocatable,public           :: grad_sphere(:)
  real(kind=r8_kind),public                :: energy,energy_ph,&
       grad_mean_square,grad_max_comp
  real(kind=r8_kind),public:: grad_max_sphere,grad_mean_sphere,dEdR_sphere

  !------------ public functions and subroutines ------------------
  public gradient_read,grad_cart_to_internal,grad_cart_for_opt
  public dealloc_grad_cart_to_internal
  public sphere_grads,eg_rpmix,rp_grads
  public cart_step_g

  !================================================================
  ! End of public interface of module
  !================================================================
  integer(i4_kind)                  :: n_cart_grads

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine cart_step_g(geo_loop,converged)
    ! Calculate numerical gradients
    ! To activate:
    ! PG input - operations_geo_opt = t
    !            max_geo_iterations - should be enough
    !                                 to calculate all
    !                                 choosen componets
    ! optimizer input - within namelist "operations"
    !           calc_cart_grad = t
    use filename_module, only: inpfile
    use iounitadmin_module
    use coordinates_module
    implicit none
    integer(i4_kind),intent(in)   :: geo_loop
    logical         ,intent(out)  :: converged
    ! *** end of interface ***

    type ::cart_coor_map
       integer(i4_kind) :: n_atoms
       integer(i4_kind) :: atom_numbers(500)
       integer(i4_kind) :: n_xyz(500)
       integer(i4_kind) :: xyz_numbers(500,3)
       integer(i4_kind) :: current_number(500*3)
    end type cart_coor_map
    type(cart_coor_map) :: xyz_comp
    integer(i4_kind) :: n_calcs

    real(r8_kind) :: cxyz(3*n_atoms),en,en1
    real(r8_kind) :: grad(n_atoms,3)
    real(r8_kind) :: step
    integer(i4_kind) :: i_step
    real(r8_kind) :: sstep,dstep
    integer(i4_kind) :: i,j,k,l,m
    integer(i4_kind) :: ienergies,loop,igrads
    logical :: yes

    DPRINT  'cart_step_g: entered, geo_loop=',geo_loop

    xyz_comp%n_atoms=0; k=0
    do i=1,n_atoms
       if(logic_coor_map(i,1).or.logic_coor_map(i,2).or. &
            logic_coor_map(i,3))then
          xyz_comp%n_atoms=xyz_comp%n_atoms+1
          k=k+1
          xyz_comp%atom_numbers(k)=i
          xyz_comp%n_xyz(k)=0; l=0
          do j=1,3
             if(logic_coor_map(i,j)) then
                xyz_comp%n_xyz(k)=xyz_comp%n_xyz(k)+1
                l=l+1
                xyz_comp%xyz_numbers(k,l)=j
             end if
          end do
       end if
    end do

    n_calcs=0
    do i=1,xyz_comp%n_atoms
       do j=1,xyz_comp%n_xyz(i)
          n_calcs=n_calcs+1
       end do
    end do

    k=0
    do i=1,xyz_comp%n_atoms
       l=xyz_comp%atom_numbers(i)
       do j=1,xyz_comp%n_xyz(i)
          m=xyz_comp%xyz_numbers(i,j)
          k=k+1
          xyz_comp%current_number(k)=3*(l-1)+m
       end do
    end do
    n_calcs=2*n_calcs+1


    sstep=step_size
    dstep=2.0_r8_kind*sstep
    converged=.false.

    ienergies=get_iounit()
    if(geo_loop > 1) then
       inquire(file=trim(inpfile('energies')), exist=yes)
       if(.not.yes) call error_handler('GRADIENT_MODULE: cart_step_g: file "energies" does not exist ')
    end if
    open(ienergies,file=trim(inpfile('energies')), err=100)

    igrads=get_iounit()

    if(geo_loop==1) then

       open(igrads,file=trim(inpfile('grads.cart')), err=101)

       write(igrads,*) 'Analytical Gradients'
       k=0
       do i=1,n_atoms+n_dummy
          if(.not.dummy_list(i)) then
             k=k+1
             write(igrads,'(I5,5X,3F17.12)') k,grad_cartes(i,:)
          end if
       end do

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
       write(ienergies,*) cxyz(:)
    else
       ! read in the central geometry from the file:
       read(ienergies,*) cxyz(:)

       do i=1,geo_loop-2
          read(ienergies,*) loop,en
       end do
       write(ienergies,*) geo_loop,energy

       if(geo_loop == n_calcs) then
          converged=.true.

          rewind ienergies
          read(ienergies,*) cxyz(:)

          grad=0.0_r8_kind

          do i=1,xyz_comp%n_atoms
             k=xyz_comp%atom_numbers(i)
             do j=1,xyz_comp%n_xyz(i)
                l=xyz_comp%xyz_numbers(i,j)
                read(ienergies,*) loop,en
                read(ienergies,*) loop,en1
                grad(k,l)=(en-en1)/dstep
             end do
          end do

          ! dump cartesian gradients into a file:
          open(igrads,file=trim(inpfile('grads.cart')), position="append",err=101)

          write(igrads,*) 'Numerical Gradients'
          do i=1,n_atoms
             write(igrads,'(I5,5X,3F17.12)') i,grad(i,:)
          end do
       end if
    end if

    call returnclose_iounit(igrads)
    call returnclose_iounit(ienergies)

    DPRINT  'cart_step_g: exit, converged=',converged
    if( .not.converged )then
      if (mod(geo_loop,2) /= 0) then
         step=sstep
         i_step=int(geo_loop/2)+1
      else
         step=-sstep
         i_step=int(geo_loop/2)
      end if

      i_step=xyz_comp%current_number(i_step)

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
100 call error_handler('GRADIENT_MODULE: cart_step_g: error - file "energies"')
101 call error_handler('GRADIENT_MODULE: cart_step_g: error - file "grads.cart"')
  end subroutine cart_step_g

  subroutine gradient_alloc(allocate_dervs)
    !  Purpose: allocates space for gradient in cartesian as well
    !           as internal coordinates
    !------------ Modules used ----------------------------------
    use math_module,only: zero
    use allocopt_module, only: allocopt_stat
    implicit none
    logical,optional:: allocate_dervs
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    if (.not.allocated(grad_cartes)) then
       allocate(grad_cartes(n_cart_grads,3),STAT=allocopt_stat(3))
       ASSERT(allocopt_stat(3).eq.0)
       allocopt_stat(3)=1
       grad_cartes=zero

     if(present(allocate_dervs)) then
      if(allocate_dervs) then
       allocate(dervs_cartes(n_atoms+n_dummy,n_atoms+n_dummy,3,3),STAT=allocopt_stat(4))
       ASSERT(allocopt_stat(4).eq.0)
       dervs_cartes=zero
       DPRINT 'dervs_cartes allocated'
      endif
     endif
    else
       write(OPT_STDOUT,*)"You stupid programmer! grad_cartes is already allocated"
    endif

  end subroutine gradient_alloc

  subroutine gradient_read (io_gx, zero, io_gxstore)
    ! Purpose: read in the scf and post-scf energies and
    !          the actual gradient from the gxfile. The gradient
    !          in (in cartesians) includes, of course, also those
    !          componentes which arise from dummy atoms, although
    !          these are set and kept to zero. To fill in only those
    !          which belong to a 'real' atom, the dummy_list is used.
    implicit none
    integer(kind=i4_kind),intent(in)  :: io_gx
    integer(kind=i4_kind),intent(in),optional  :: io_gxstore
    !** End of interface ****************************************

    integer(kind=i4_kind)             :: n,i
    integer(kind=i4_kind)             :: stat
    logical, optional                 :: zero ! do not raed gradients,
    ! but set to zero. Can be used for testing gxfiles.
    logical                           :: zero_loc

    DPRINT "gm::gradient_read: entered",io_gx
    zero_loc=.false.
    if(present(zero)) zero_loc=zero

    n_cart_grads=n_atoms+n_dummy
#ifdef WITH_EFP
    if(efp .and. n_efp > 0) n_cart_grads=n_cart_grads+2*n_efp
#endif

    call gradient_alloc(allocate_dervs=analitic_hessian_calculated.or.update_fromcartessian)

    if(zero_loc) then
       grad_cartes=0.0_r8_kind
       energy=0.0_r8_kind
       energy_ph=0.0_r8_kind
    else
      if(.not.calc_epeff_hessian) then
       read(io_gx,'(2F24.12)',IOSTAT=stat) energy,energy_ph
       if ( stat /= 0 ) then
         print *,'OPTIMIZER: Error while reading DFT energy from gxfile.'
         print *,'           Make sure to run ParaGauss before using optimizer!'
         ABORT('error reading energy')
       endif
       if(present(io_gxstore)) write(io_gxstore,'(2F24.12)') energy,energy_ph
       DPRINT "gm::gradient_read: energy,energy_ph=",energy,energy_ph
       do i=1,n_cart_grads
          DPRINT "gm::gradient_read: i=",i
          if (.not.dummy_list(i)) then
             read(io_gx,'(I5,5X,3F17.12)',IOSTAT=stat) n,grad_cartes(i,:)
             if ( stat /= 0 ) then
               print *,'OPTIMIZER: Error while reading DFT gradients from gxfile.'
               print *,'           Make sure to run ParaGauss before using optimizer!'
               ABORT('error reading gradients')
             endif
             if(present(io_gxstore)) write(io_gxstore,'(I5,5X,3F17.12)')n,grad_cartes(i,:)
             DPRINT "gm::gradient_read: n,grad_cartes(i,:)=",n,grad_cartes(i,:)
          else
             DPRINT "gm::gradient_read: skip dummy"
          endif
       enddo
      else
        grad_cartes=0.0_r8_kind
        energy=0.0_r8_kind
        energy_ph=0.0_r8_kind
      endif
    end if
  end subroutine gradient_read

  function eg_rpmix(epe_kl,epe_nucen,epe) result(dev)
   integer(kind=i4_kind), intent(in):: epe_kl,epe_nucen

  type (epe_shells):: epe(:)


   real(kind=r8_kind) :: dev,r2mean,zero=0.0_r8_kind,dev_int,dev_out
   real(kind=r8_kind) :: pmix=0.2_r8_kind
   integer(kind=i4_kind):: i,j
   real(kind=r8_kind) :: r2_r,r2_p

   pmix=rpmix
    n_cart_grads=n_atoms+n_dummy
    call gradient_alloc(allocate_dervs=.true.)
  dev=0.0_r8_kind
  dev_out=0.0_r8_kind
  dev_int=0.0_r8_kind
   print* ,'epe_kl,epe_nucen+epe_kl-1', epe_kl,epe_nucen+epe_kl-1
  do i=1,n_atoms+n_dummy
   if(atom(i)%dummy) cycle
   do j = i+1, n_atoms+n_dummy
   if(atom(j)%dummy) cycle
   r2_r=sum((xyz_reactant(:,i)-xyz_reactant(:,j))**2)
   r2_p=sum((xyz_product(:,i)-xyz_product(:,j))**2)
   r2mean= ( sum((xyz(:,i)-xyz(:,j))**2)-pmix*r2_r-(1-pmix)*r2_p)**2
   dev=dev+r2mean/(r2_r+r2_p)**2
   dev_int=dev_int+r2mean/(r2_r+r2_p)**2
   enddo

   if(epe_forces) then
   do j=epe_kl,epe_nucen+epe_kl-1
    r2_r=sum((xyz_reactant(:,i)-epe(j)%s)**2)
    r2_p=sum((xyz_product(:,i)-epe(j)%s)**2)
    r2mean=( sum((xyz(:,i)-epe(j)%s)**2)-pmix*r2_r-(1-pmix)*r2_p)**2
    dev=dev+r2mean/(r2_r+r2_p)**2
    dev_out=dev_out+r2mean/(r2_r+r2_p)**2
   enddo
   endif
  enddo
 print*, 'dev calculated',dev, pmix, 'rpmix'
   print*,'deviatrions from internal part', dev_int
   print*,'deviatrions from outer part', dev_out
  grad_cartes = zero
  do i=1,n_atoms+n_dummy
   if(atom(i)%dummy) cycle
   do j = 1, n_atoms+n_dummy
    if(i.eq.j.or.atom(j)%dummy) cycle

    r2_r=sum((xyz_reactant(:,i)-xyz_reactant(:,j))**2)
    r2_p=sum((xyz_product(:,i)-xyz_product(:,j))**2)


    grad_cartes(i,:)=grad_cartes(i,:)+4*(xyz(:,i)-xyz(:,j)) &
      *( sum( (xyz(:,i)-xyz(:,j))**2 )-pmix*r2_r-(1-pmix)*r2_p)/(r2_r+r2_p)**2
   enddo

   if(epe_forces) then
   do j = epe_kl,epe_nucen+epe_kl-1

    r2_r=sum((xyz_reactant(:,i)-epe(j)%s)**2)
    r2_p=sum((xyz_product(:,i)-epe(j)%s)**2)


    grad_cartes(i,:)=grad_cartes(i,:)+4*(xyz(:,i)-epe(j)%s) &
      *( sum( (xyz(:,i)-epe(j)%s)**2 )-pmix*r2_r-(1-pmix)*r2_p)/(r2_r+r2_p)**2
   enddo
   endif
    print*,'grad cart', grad_cartes(i,:),dev
  enddo

  end function eg_rpmix
   !*************************************************************
  subroutine grad_cart_to_internal()
    ! Purpose: transform the gradient from cartesian coordinates
    !          to internal coordinates using the inverse bmat:
    !
    !          grad_intern = B**-1 * grad_cartes
    !          B**-1 = bmat_inv
    use coordinates_module, only: bmat_inv, tmat, reduc_mat
    use math_module, only: zero, round, invert_matrix
    USE_DEBUG
    implicit none
    !** End of interface ****************************************

    integer (i4_kind) :: k, i, start, alloc_stat
    real (r8_kind), allocatable, target :: grad_prim(:)
    real (r8_kind), pointer :: grad(:)

    allocate (grad_intern(n_internal), grad_prim(n_primitive), STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    grad_intern = zero
    grad_prim = zero
    if (zmat_coordinates .and. zmat_format) then
       grad => grad_prim
    else
       grad => grad_intern
    endif

    do i = 1, size (grad)
       do k = 1, n_atoms + n_dummy
          start = (k - 1) * 3 + 1
          grad(i) = grad(i) + &
               bmat_inv(start + 0, i) * grad_cartes(k, 1) + &
               bmat_inv(start + 1, i) * grad_cartes(k, 2) + &
               bmat_inv(start + 2, i) * grad_cartes(k, 3)
       enddo
       ! FIXME: hm, what is this?
       if (abs(grad(i)) <= 1.0e-10_r8_kind) grad(i) = zero
    enddo
    if (zmat_coordinates .and. zmat_format) then
       grad_intern = matmul (reduc_mat, grad_prim)
    endif
    grad_intern = matmul (tmat, grad_intern)

    if (print_debug) then
       write (OPT_STDOUT, *) " full internal gradient is:"
       do i = 1, size (grad)
          write (OPT_STDOUT, '(I4,3x,f9.5)') i, grad(i)
       enddo
       write (OPT_STDOUT, *) " Symmetry reduced gradient is "
       do i = 1, n_internal
          write (OPT_STDOUT, '(i4,3x,f9.5)') i, grad_intern(i)
       enddo
    endif

    grad_max_comp = maxval (abs (grad_intern))
    grad_mean_square = sqrt (sum (grad_intern**2) / n_internal)

    DCALL show("Internal gradient", grad_intern(:))
  end subroutine grad_cart_to_internal

  function sphere_grads(kk,grad_max_sphere,grad_mean_sphere) result(grad_sphere)
    real(kind=r8_kind):: grad_sphere(size(grad_intern))
    real(kind=r8_kind),intent(out):: grad_max_sphere, grad_mean_sphere
    real(kind=r8_kind):: ssi,zero=0.0_r8_kind


    integer(kind=i4_kind), intent(in):: kk
    integer(kind=i4_kind)   :: i,ii

    write(io_flepo,*)'sphere_dependent_var', kk, s(kk)%value-s_reactant(kk)%value

    dEdR_sphere=zero
    ii=0
   do i=1,size(grad_intern)

    ssi=s(i)%value-s_reactant(i)%value

    if(i.eq.kk) cycle
    if(select_sphere_vars.and..not.s(i)%sphere) ssi=zero
    ii=ii+1
    grad_sphere(ii)=grad_intern(i)-grad_intern(kk)*ssi/(s(kk)%value-s_reactant(kk)%value)
    !write(io_flepo,*) ii,i,grad_sphere(ii),grad_intern(i)
   enddo

    grad_max_sphere = maxval(abs(grad_sphere(1:n_internal-1)))
    grad_mean_sphere = sqrt(sum(grad_sphere(1:n_internal-1)**2)/size(grad_sphere))
    dEdR_sphere=grad_intern(kk)*distance_to_reactant/(s(kk)%value-s_reactant(kk)%value)
    write(io_flepo,*) 'grad_max_sphere ',grad_max_sphere,'         grad_mean_sphere ', grad_mean_sphere
    write(io_flepo,*) 'dEdR_sphere' , dEdR_sphere

  end function sphere_grads

  function rp_grads(rp_var,kk,grad_max_rp,grad_mean_rp,dEdR_rp) result(grad_rp)
    real(kind=r8_kind):: grad_rp(size(grad_intern))
    real(kind=r8_kind),intent(out):: grad_max_rp, grad_mean_rp,dEdR_rp
    real(kind=r8_kind),intent(in):: rp_var
    real(kind=r8_kind):: ssr,ssp,zero=0.0_r8_kind,sqt,sqmt,ssk,dkdi,drdt
    real(kind=r8_kind):: R2,P2


    integer(kind=i4_kind), intent(in):: kk
    integer(kind=i4_kind)   :: i,ii

    sqt=rp_var**2
    sqmt=(1.0_r8_kind-rp_var)**2

    write(io_flepo,*)'rp_dependent_var', kk, s(kk)%value-s_reactant(kk)%value, &
                                             s(kk)%value-s_product(kk)%value

    dEdR_rp=zero

    ii=0
    ssk= (s(kk)%value-s_reactant(kk)%value)*sqt-(s(kk)%value-s_product(kk)%value)*sqmt
   do i=1,size(grad_intern)

    ssr=(s(i)%value-s_reactant(i)%value)*sqt
    ssp=(s(i)%value-s_product(i)%value)*sqmt

    dkdi=-(ssr-ssp)/ssk

    if(i.eq.kk) cycle
    if(select_sphere_vars.and..not.s(i)%sphere) dkdi=zero
    ii=ii+1
    grad_rp(ii)=grad_intern(i)+grad_intern(kk)*dkdi
    write(io_flepo,*) ii,i,grad_rp(ii),grad_intern(i),dkdi, 'dkdi'
   enddo

    grad_max_rp = maxval(abs(grad_rp(1:n_internal-1)))
    grad_mean_rp = sqrt(sum(grad_rp(1:n_internal-1)**2)/size(grad_rp))

    P2=distance_to_product**2
    R2=distance_to_reactant**2
    drdt=-(P2*(1.0_r8_kind-rp_var)+R2*rp_var)/ssk
    dEdR_rp=grad_intern(kk)*drdt



    write(io_flepo,*) 'grad_max_rp ',grad_max_rp,'         grad_mean_rp ', grad_mean_rp
    write(io_flepo,*) 'dEdR_rp' , dEdR_rp

  end function rp_grads

  subroutine dealloc_grad_cart_to_internal()
    integer(kind=i4_kind)   :: alloc_stat

    if(allocated(grad_intern)) then
        deallocate(grad_intern, STAT=alloc_stat)
        if (alloc_stat/=0) call error_handler&
         ("gradient_alloc : deallocation (1) failed")
    endif
  end subroutine dealloc_grad_cart_to_internal

  !*************************************************************
  subroutine grad_cart_for_opt ()
    use math_module, only: zero
    implicit none
    !** End of interface ****************************************

    integer(kind=i4_kind)   :: i,alloc_stat
    real (r8_kind), allocatable, target :: grad_prim(:)

    allocate(grad_intern(n_internal),grad_prim(n_primitive), STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)
    grad_intern=zero
    grad_prim =zero

    do i=1,n_cart_grads
       grad_prim(3*i-2)=grad_cartes(i,1)
       grad_prim(3*i-1)=grad_cartes(i,2)
       grad_prim(3*i)  =grad_cartes(i,3)
    end do

    grad_intern=grad_prim(1:n_internal)
#ifdef WITH_EFP
    if(qm_fixed) then
       grad_intern=grad_prim(3*n_atoms+1:3*n_atoms+n_internal)
    end if
#endif

    grad_max_comp = maxval(abs(grad_intern))
    grad_mean_square = sqrt(sum(grad_intern**2)/n_internal)

  end subroutine grad_cart_for_opt
  !*************************************************************

  !--------------- End of module ----------------------------------
end module gradient_module
