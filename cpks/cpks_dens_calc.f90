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
module cpks_dens_calc
  !---------------------------------------------------------------
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

! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:   &
                   i4_kind &
                 , r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: density_calc_ph_v2
  public :: fix_rho_grads

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  integer(i4_kind), parameter, private ::&
       & VALUE=0, X=1, Y=2, Z=3,&
       & UP=1, DN=2,&
       & UPUP=1, DNDN=2, UPDN=3,&
       & RE=1, IM=2


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine density_calc_ph_v2( vl, ispin                       &
                               , rho,gamma,grarho                & ! out: rho, gamma, and rho gradient
                               , nuc_grarho                      & ! out: rho gradient wrt nuc
                               , orbs_ob,orbs_grads,nuc_grads    & ! input
                               , phi_ob                          & ! output
                               , phi_gr_nuc                      & ! output
                               , rho_gr_nuc                      & ! output
                               , nuc_dervsrho                    & ! dervs only
                               , nuc_dervs_grarho                &
                               , nuc_sec_der                     &
                               , atnum)                            ! grid hosting UA

    ! Purpose : Calculation of the density
    !           Also the gradient of the density is calculated
    use eigen_data_module
    use machineparameters_module
    use occupation_module
    use cpks_grid_utils, only: eigMO
    use datatype, only: arrmat4, arrmat5, arrmat6, arrmat7
    use unique_atom_module, only: n_unique_atoms, unique_atoms
    use occupied_levels_module, only: eigvec_occ, occ_num_occ
    use orbitalstore_module, only: orbital_type &
                                 , orbital_spin_type &
                                 , orbital_gradient_type &
                                 , orbital_nuclear_gradient_type &
                                 , orbital_nuclear_sec_der_type
    USE_MEMLOG
    implicit none
    integer(kind=i4_kind),intent(in) :: vl
    integer(kind=i4_kind),intent(in) :: ispin
    real(kind=r8_kind),intent(out)   :: rho(:,:),gamma(:,:),grarho(:,:,:) 

    ! these are derivatives wrt nuclear coordinates:
    type(arrmat4),intent(inout)            :: nuc_grarho(:)         ! OUTPUT: gradient of rho
    type(arrmat6),intent(inout), optional  :: nuc_dervsrho(:,:)     ! OUTPUT: dervs of rho wrt nuclear coordinates
    type(arrmat7),intent(inout), optional  :: nuc_dervs_grarho(:,:) ! OUTPUT: dervs only, partial dervs of rho

    type(orbital_type)     ,intent(in)     :: orbs_ob(:)    ! INPUT : basis orbitals

    type(orbital_spin_type),intent(inout)  :: phi_ob(:)         ! OUTPUT: MO-orbitals (eigenfunctions)
                                            ! phi_ob(irr)%o(:vl,n_orbs,n_part,ispin)

    type(arrmat5)          ,intent(inout)  :: phi_gr_nuc(:)     ! OUTPUT: MO-orbital grads wrt nuc
                                            ! phi_gr_nuc(irr)%m(:vl,n_orbs,n_part,n_modes,ispin)

    real(r8_kind)          ,intent(out)    :: rho_gr_nuc(:,:,:) ! OUTPUT: density grads wrt nuc
    optional :: rho_gr_nuc                  ! rho_gr_nuc(:vl,n_modes,ispin)

    type(orbital_gradient_type),intent(in) :: orbs_grads(:) ! INPUT : gradients of basis orbitals

    type(orbital_nuclear_gradient_type),intent(in) :: nuc_grads(:,:) ! intent(in)
    type(orbital_nuclear_sec_der_type),intent(in),optional :: nuc_sec_der(:,:)

    integer(i4_kind),intent(in)            :: atnum  ! grid hosting UA

    optional :: phi_ob            ! required for CPKS RHS and final second derivatives.
    optional :: phi_gr_nuc, atnum ! required for CPKS RHS and final second derivatives.
    !** End of interface *****************************************

    real(kind=r8_kind)  :: &
         occ_real,occ_real_2,help_vec(vl)
    real(kind=r8_kind),pointer  :: ob(:,:,:),occ(:,:),eigv(:,:,:),&
         grads(:,:,:,:),nuc(:,:,:,:,:),nuc_grarho_p(:,:,:,:)
    integer(i4_kind) :: s,j,l,eig_dim,counter,i_ua,i_ea,&
         full_vl,occ_dim,orb_dim
    integer(i4_kind) :: i_ua2,i_ea2,orb_dim2,counter2,i_irr
    real(kind=r8_kind),parameter:: two=2.0_r8_kind

    real(kind=r8_kind), pointer :: phi_arr(:,:)
    real(kind=r8_kind), allocatable :: graphi_arr(:,:,:), &
                        grad_help(:,:,:), nuc_help(:,:,:,:),  &
                        help_nuc_arr(:,:,:,:),help_sec_der_arr(:,:,:,:), &
                        sec_der_help(:,:,:,:)


   type(arrmat4),allocatable :: help_nuc_arr2(:)

    integer(i4_kind) :: xyz
    integer(i4_kind) :: i_grad

    integer(i4_kind) :: alloc_stat(34)

    integer(i4_kind) :: n_irrep                 ! := size(orbs_ob)
    integer(i4_kind) :: partners(size(orbs_ob)) ! (n_irrep)

    ! for convenient loop ranges:
    n_irrep = size(orbs_ob)
    do i_irr=1,n_irrep
      partners(i_irr) = size(orbs_ob(i_irr)%o,3)
    enddo

    full_vl=machineparameters_veclen
    
        if(present(nuc_dervsrho)) then
         allocate(help_nuc_arr2(n_unique_atoms),stat=alloc_stat(16))
         ASSERT(alloc_stat(16).eq.0)
        endif

    rho(1:vl,:)=0.0_r8_kind
    grarho = 0.0_r8_kind
    do i_ua=1,n_unique_atoms               
       nuc_grarho(i_ua)%m(1:vl,:,:,:)=0.0_r8_kind

    if(present(nuc_dervsrho)) then
     do i_ua2=1,n_unique_atoms
       nuc_dervsrho(i_ua,i_ua2)%m(1:vl,:,:,:,:,:)=0.0_r8_kind
     enddo
    endif

    enddo

#if 1
    if( present(phi_ob) )then
    ![[=== calculate molecular orbitlas ===
    ! FIXME: sometimes only the occupied orbitals are needed...
    do i_irr=1,n_irrep
      ASSERT(i_irr<=size(phi_ob))
      ASSERT(associated(phi_ob(i_irr)%o))
    do s=1,ispin
    do l=1,partners(i_irr)
      phi_ob(i_irr)%o(:,:,l,s) = 0.0_r8_kind
      call eigMO(vl, 0                     &
                , orbs_ob(i_irr)%o(:,:,l)  &
                , eigvec(i_irr)%m(:,:,s)   &
                , phi_ob(i_irr)%o(:,:,l,s) &
                )
    enddo
    enddo
    enddo
    !]]====================================
    endif

    ![[=== calculate molecular orbital grads wrt nuc ===
    if( present(phi_gr_nuc) )then
      ASSERT(present(atnum))
      ! FIXME: sometimes only the occupied orbitals are needed...
      do i_irr=1,n_irrep
         ASSERT(i_irr<=size(phi_gr_nuc))
         !SSERT(associated(phi_gr_nuc(i_irr)%m))
      do s=1,ispin
        call orbNGra( nuc_grads(:,i_irr)             &
                    , eigvec(i_irr)%m(:,:,s)         &
                    , phi_gr_nuc(i_irr)%m(:,:,:,:,s) &
                    , atnum                          &
                    )
      enddo
      enddo
    endif
    !]]=================================================

#if 1
    if( present(phi_ob) )then
    ! if not present -- compute rho the old way down there!
    ![[=== calculate density of occupied orbitals ======
    rho(:,:) = 0.0_r8_kind
    do i_irr=1,n_irrep
      occ     => occ_num_occ(i_irr)%m
      occ_dim =  size(eigvec_occ(i_irr)%m,2) ! FIXME: a better source?
      if( occ_dim == 0 ) cycle
      do s=1,ispin
        do l=1,partners(i_irr)
          do j=1,occ_dim
            occ_real=occ(j,s)/partners(i_irr)
            ! density of the (occupied) orbitals
            rho(:vl,s) = rho(:vl,s)                    &
                       + phi_ob(i_irr)%o(:vl,j,l,s)**2 &
                       * occ_real
          enddo
        enddo
      enddo
    enddo
    !]]=================================================
    endif
#endif

#if 1
    ![[=== calculate density grad wrt nuc ==============
    if( present(rho_gr_nuc) )then
      ASSERT(present(phi_ob))
      ASSERT(present(phi_gr_nuc))
      rho_gr_nuc(:,:,:) = 0.0_r8_kind
      do i_irr=1,size(phi_ob) ! n_irrep
        occ     => occ_num_occ(i_irr)%m
        occ_dim =  size(eigvec_occ(i_irr)%m,2) ! FIXME: a better source?
        if( occ_dim == 0 ) cycle
        do s=1,ispin
          do i_grad=1,size(phi_gr_nuc(i_irr)%m,4) ! n_modes
            do l=1,partners(i_irr)
              do j=1,occ_dim
                occ_real=occ(j,s)/partners(i_irr)
                ! density gradient wrt nuc: SUM[j]{  2 * occ(j) * phi(j) * d/dN phi(j)
                rho_gr_nuc(:vl,i_grad,s) = rho_gr_nuc(:vl,i_grad,s) &
                           + phi_ob(i_irr)%o(:vl,j,l,s)             &
                           * phi_gr_nuc(i_irr)%m(:vl,j,l,i_grad,s)  &
                           * occ_real * 2
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    !]]=================================================
#endif
#endif

    irr: do i_irr=1,n_irrep
       ob      => orbs_ob(i_irr)%o
       occ     => occ_num_occ(i_irr)%m
       eigv    => eigvec_occ(i_irr)%m
       grads   => orbs_grads(i_irr)%o
       eig_dim =  size(eigv,1)
       occ_dim =  size(eigv,2)
       if( occ_dim == 0 ) cycle irr

       allocate(phi_arr(vl,occ_dim),graphi_arr(vl,3,occ_dim),stat=alloc_stat(8))
!      allocate(graphi_arr(vl,3,occ_dim),stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(34)=0 ! graphi_arr
       MEMLOG(size(phi_arr)+size(graphi_arr))

       spin: do s=1,ispin
         part: do l=1,partners(i_irr)

             ! calculation of orbitals and their grads

             ! FIXME: use phi_ob instead:
             call dgemm('n','n',vl,occ_dim,&
                  eig_dim,1.0_r8_kind,ob(1,1,l),full_vl,eigv(1,1,s),&
                  eig_dim,0.0_r8_kind,phi_arr,vl)
!            phi_arr => phi_ob(i_irr)%o(:,:,l,s)

             ! and gradient
             if(vl==full_vl) then	
                call dgemm('n','n',vl*3,occ_dim,&
                     eig_dim,1.0_r8_kind,grads(1,1,1,l),full_vl*3,eigv(1,1,s),&
                     eig_dim,0.0_r8_kind,graphi_arr(1,1,1),vl*3)
             else
                allocate(grad_help(vl,3,eig_dim),stat=alloc_stat(9))
                ASSERT(alloc_stat(9).eq.0)
                grad_help(1:vl,:,:)=grads(1:vl,:,:,l)
                call dgemm('n','n',vl*3,occ_dim,&
                     eig_dim,1.0_r8_kind,grad_help(1,1,1),vl*3,eigv(1,1,s),&
                     eig_dim,0.0_r8_kind,graphi_arr(1,1,1),vl*3)	          
                deallocate(grad_help,stat=alloc_stat(9))
                ASSERT(alloc_stat(9).eq.0)
                alloc_stat(9)=1
             endif

                
             calc_grarho: do j=1,occ_dim
                occ_real=occ(j,s)/partners(i_irr)
                occ_real_2=occ_real*2.0_r8_kind
                if( .not. present(phi_ob) )then
                ! compute rho traditional way:
                ! Calculation of the density of the orbital               
                rho(1:vl,s)=rho(1:vl,s)+phi_arr(1:vl,j)*phi_arr(1:vl,j)*occ_real
                endif

                ! calculation of the gradient of the density
                help_vec=phi_arr(:,j)*occ_real_2
                phi_arr(:,j)=phi_arr(:,j)*occ_real_2  !!! modefied from now on
                do xyz=1,3
                   grarho(1:vl,xyz,s) = grarho(1:vl,xyz,s)&
                        & + help_vec*graphi_arr(:,xyz,j)
                enddo
             enddo calc_grarho
                
             counter=0

               ! now calculate gradient with respect to nuclear displacements

           nuc_grarho_dervsrho_ua:  do i_ua=1,n_unique_atoms
                nuc=> nuc_grads(i_ua,i_irr)%o
                orb_dim=size(nuc,4)
                if (orb_dim == 0) cycle nuc_grarho_dervsrho_ua

                allocate(help_nuc_arr(vl,3,unique_atoms(i_ua)%n_equal_atoms,occ_dim),& 
                         stat=alloc_stat(10))
                         ! to hold orbital nuclear diplacement grads
                ASSERT(alloc_stat(10).eq.0)

                if(present(nuc_dervsrho)) then
                 allocate(help_sec_der_arr(vl,6,unique_atoms(i_ua)%n_equal_atoms,occ_dim),&
                          stat=alloc_stat(18))
                 ASSERT(alloc_stat(18).eq.0)
                endif

                if(vl==full_vl) then
                   call dgemm('n','n',vl*3*unique_atoms(i_ua)%n_equal_atoms,occ_dim,orb_dim, &
                        1.0_r8_kind, &
                        nuc_grads(i_ua,i_irr)%o(1,1,1,1,l),vl*3*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1,s),eig_dim, &
                        0.0_r8_kind,help_nuc_arr(1,1,1,1), &
                        vl*3*unique_atoms(i_ua)%n_equal_atoms)

                  if(present(nuc_dervsrho)) then
                   ! now help_sec_der_arr
                   call dgemm('n','n',vl*6*unique_atoms(i_ua)%n_equal_atoms,&
                                      occ_dim,orb_dim, &
                        1.0_r8_kind,&
                        nuc_sec_der(i_ua,i_irr)%o(1,1,1,1,l), &
                                   vl*6*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1,s),eig_dim, &
                        0.0_r8_kind,help_sec_der_arr(1,1,1,1),&
                                    vl*6*unique_atoms(i_ua)%n_equal_atoms)
                  endif

                else
                   ! copy data first
                   allocate(nuc_help(vl,3,unique_atoms(i_ua)%n_equal_atoms,orb_dim), &
                            stat=alloc_stat(12))
                   ASSERT(alloc_stat(12).eq.0)

                   nuc_help(1:vl,:,:,:)=nuc_grads(i_ua,i_irr)%o(1:vl,:,:,:,l)

                   call dgemm('n','n',vl*3*unique_atoms(i_ua)%n_equal_atoms,&
                                      occ_dim,orb_dim, &
                        1.0_r8_kind,&
                        nuc_help(1,1,1,1),vl*3*unique_atoms(i_ua)%n_equal_atoms,&
                        eigv(1+counter,1,s),eig_dim, &
                        0.0_r8_kind,help_nuc_arr(1,1,1,1),&
                                    vl*3*unique_atoms(i_ua)%n_equal_atoms)

                   if(present(nuc_dervsrho)) then
                    allocate(sec_der_help(vl,6,unique_atoms(i_ua)%n_equal_atoms,orb_dim),&
                              stat=alloc_stat(17))
                     ASSERT(alloc_stat(17).eq.0)
                   
                     sec_der_help(1:vl,:,:,:)=nuc_sec_der(i_ua,i_irr)%o(1:vl,:,:,:,l)
                   
                     ! now help_sec_der_arr
                     call dgemm('n','n',vl*6*unique_atoms(i_ua)%n_equal_atoms,&
                          occ_dim,orb_dim, &
                          1.0_r8_kind,&
                          sec_der_help(1,1,1,1),vl*6*unique_atoms(i_ua)%n_equal_atoms,&
                          eigv(1+counter,1,s),eig_dim,0.0_r8_kind,help_sec_der_arr(1,1,1,1),&
                          vl*6*unique_atoms(i_ua)%n_equal_atoms)
                     deallocate(sec_der_help, stat=alloc_stat(17))
                     ASSERT(alloc_stat(17).eq.0)
                     alloc_stat(17)=1
                   endif

                   deallocate(nuc_help, stat=alloc_stat(12))
                   ASSERT(alloc_stat(12).eq.0)
                   alloc_stat(12)=1
                endif

                nuc_grarho_p=>nuc_grarho(i_ua)%m

                mo: do j=1,occ_dim
                   equals: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      

                    ! here and below phi_arr = phi_arr*2*occ(j), thus no other
                    ! occupation factor is required when phi_arr is contributing:

                      nuc_grarho(i_ua)%m(1:vl,:,i_ea,s)=&
                           nuc_grarho(i_ua)%m(1:vl,:,i_ea,s)+&
                           help_nuc_arr(1:vl,:,i_ea,j)*spread(phi_arr(:,j),2,3)

                     if(present(nuc_dervsrho)) then
                      nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(:vl,1,i_ea,1,i_ea,s)+ &
                              help_sec_der_arr(:vl,1,i_ea,j)*phi_arr(:,j)

                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,1,i_ea,2,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,1,i_ea,2,i_ea,s)+ &
                              help_sec_der_arr(1:vl,2,i_ea,j)*phi_arr(:,j)

                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,2,i_ea,1,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,2,i_ea,1,i_ea,s)+ &
                              help_sec_der_arr(1:vl,2,i_ea,j)*phi_arr(:,j)

                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,1,i_ea,3,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,1,i_ea,3,i_ea,s)+ &
                              help_sec_der_arr(1:vl,3,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,3,i_ea,1,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,3,i_ea,1,i_ea,s)+ &
                              help_sec_der_arr(1:vl,3,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,2,i_ea,2,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,2,i_ea,2,i_ea,s)+ &
                              help_sec_der_arr(1:vl,4,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,2,i_ea,3,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,2,i_ea,3,i_ea,s)+ &
                              help_sec_der_arr(1:vl,5,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,3,i_ea,2,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,3,i_ea,2,i_ea,s)+ &
                              help_sec_der_arr(1:vl,5,i_ea,j)*phi_arr(:,j)
                      nuc_dervsrho(i_ua,i_ua)%m(1:vl,3,i_ea,3,i_ea,s)= &
                          nuc_dervsrho(i_ua,i_ua)%m(1:vl,3,i_ea,3,i_ea,s)+ &
                              help_sec_der_arr(1:vl,6,i_ea,j)*phi_arr(:,j)
                     endif ! present nuc_dervsrho
                enddo  equals
                enddo mo

                 counter2=0
           ![[=== second derivatives wrt nuc coordinates ===
           if(present(nuc_dervsrho)) then
           ua2:  do i_ua2=1,n_unique_atoms
                orb_dim2=size(nuc_grads(i_ua2,i_irr)%o,4)
                if( orb_dim2 == 0) cycle ua2
               
                calc_help_nuc_arr2: if(.not.associated(help_nuc_arr2(i_ua2)%m)) then 

                allocate(help_nuc_arr2(i_ua2)%m(vl,3,unique_atoms(i_ua2)%n_equal_atoms,occ_dim),&
                         stat=alloc_stat(14))
                        ! to hold orbital nuc displacement grads for second unique
                ASSERT(alloc_stat(14).eq.0)

                fullorelse: if(vl==full_vl) then
                   call dgemm('n','n',vl*3*unique_atoms(i_ua2)%n_equal_atoms,occ_dim,orb_dim2, &
                        1.0_r8_kind, &
                        nuc_grads(i_ua2,i_irr)%o(1,1,1,1,l),vl*3*unique_atoms(i_ua2)%n_equal_atoms,&
                        eigv(1+counter2,1,s),eig_dim, &
                        0.0_r8_kind,help_nuc_arr2(i_ua2)%m(1,1,1,1), &
                        vl*3*unique_atoms(i_ua2)%n_equal_atoms)

                else fullorelse
                   ! copy data first
                   allocate(nuc_help(vl,3,unique_atoms(i_ua2)%n_equal_atoms,&
                                                    orb_dim2), &
                            sec_der_help(vl,6,unique_atoms(i_ua2)%n_equal_atoms,&
                                                        orb_dim2),&
                            stat=alloc_stat(15))
                   ASSERT(alloc_stat(15).eq.0)

                   nuc_help(1:vl,:,:,:)=      nuc_grads(i_ua2,i_irr)%o(1:vl,:,:,:,l)
                   sec_der_help(1:vl,:,:,:)=nuc_sec_der(i_ua2,i_irr)%o(1:vl,:,:,:,l)

                   call dgemm('n','n',vl*3*unique_atoms(i_ua2)%n_equal_atoms,&
                                             occ_dim,orb_dim2, &
                        1.0_r8_kind,&
                        nuc_help(1,1,1,1),vl*3*unique_atoms(i_ua2)%n_equal_atoms,&
                        eigv(1+counter2,1,s),eig_dim, &
                        0.0_r8_kind,help_nuc_arr2(i_ua2)%m(1,1,1,1),&
                          vl*3*unique_atoms(i_ua2)%n_equal_atoms)

                   deallocate(nuc_help,sec_der_help, stat=alloc_stat(15))
                   ASSERT(alloc_stat(15).eq.0)
                   alloc_stat(15)=1
                endif fullorelse
               endif calc_help_nuc_arr2

!              if(present(nuc_dervsrho)) then
                
               ! cross contrib to nuc_dervsrho, 2nd of two contribs
           
                mo2: do j=1,occ_dim

                occ_real=occ(j,s)/partners(i_irr)
                occ_real_2=occ_real*2.0_r8_kind

                   equals1: do i_ea=1,unique_atoms(i_ua)%n_equal_atoms                      

                   equals2: do i_ea2=1,unique_atoms(i_ua2)%n_equal_atoms                      

                   ! in expresions below help_nuc_arr and help_nuc_arr2 do not
                   ! depend on occupation, thus occ_real_2 factor is used
                   ! to account for occupation dependence

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,1,i_ea,1,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,1,i_ea,1,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,1,i_ea,j)* &
                           help_nuc_arr2(i_ua2)%m(1:vl,1,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,1,i_ea,2,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,1,i_ea,2,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,1,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(1:vl,2,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,2,i_ea,1,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,2,i_ea,1,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,2,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(1:vl,1,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,2,i_ea,2,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,2,i_ea,2,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,2,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(1:vl,2,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,1,i_ea,3,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,1,i_ea,3,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,1,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(1:vl,3,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,3,i_ea,1,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,3,i_ea,1,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,3,i_ea,j)* &
                            help_nuc_arr2(i_ua2)%m(1:vl,1,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,2,i_ea,3,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,2,i_ea,3,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,2,i_ea,j)* &
                             help_nuc_arr2(i_ua2)%m(1:vl,3,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,3,i_ea,2,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,3,i_ea,2,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,3,i_ea,j)* &
                             help_nuc_arr2(i_ua2)%m(1:vl,2,i_ea2,j)

                      nuc_dervsrho(i_ua,i_ua2)%m(1:vl,3,i_ea,3,i_ea2,s)= &
                          nuc_dervsrho(i_ua,i_ua2)%m(1:vl,3,i_ea,3,i_ea2,s)+ &
                          occ_real_2*help_nuc_arr(1:vl,3,i_ea,j)* &
                             help_nuc_arr2(i_ua2)%m(1:vl,3,i_ea2,j)

                   enddo equals2
                   enddo equals1
                enddo mo2
!              endif ! present  nuc_dervsrho

                counter2=counter2+orb_dim2
           enddo ua2
           endif ! present(nuc_dervsrho)
           !]]==============================================


                counter=counter+orb_dim

                deallocate(help_nuc_arr, stat=alloc_stat(10))
                ASSERT(alloc_stat(10).eq.0)
                alloc_stat(10)=1

                if(present(nuc_dervsrho)) then
                  deallocate(help_sec_der_arr, stat=alloc_stat(18))
                  ASSERT(alloc_stat(18).eq.0)
                  alloc_stat(18)=1
                endif

             enddo nuc_grarho_dervsrho_ua

             if( allocated(help_nuc_arr2) )then
               do i_ua=1,n_unique_atoms
                 if(associated(help_nuc_arr2(i_ua)%m)) then
                  deallocate(help_nuc_arr2(i_ua)%m,stat=alloc_stat(14))
                  ASSERT(alloc_stat(14).eq.0)
                  alloc_stat(14)=1
                 endif
               enddo
             endif

          enddo part


       enddo spin
       MEMLOG(-size(phi_arr)-size(graphi_arr))
       deallocate(phi_arr,graphi_arr,stat=alloc_stat(8))
!      deallocate(graphi_arr,stat=alloc_stat(8))
       ASSERT(alloc_stat(8).eq.0)
       alloc_stat(8)=1
       alloc_stat(34)=1 ! graphi_arr
    enddo irr

         if(present(nuc_dervsrho)) then
          MEMLOG(-size(help_nuc_arr2))
          deallocate(help_nuc_arr2,stat=alloc_stat(16))
          ASSERT(alloc_stat(16).eq.0)
          alloc_stat(16)=1
         endif


    if(ispin==1) then
       gamma(1:vl,UPUP) = SUM(grarho(1:vl,X:Z,UP)**2, DIM=2)
    else
       gamma(1:vl,UPUP) = SUM(grarho(1:vl,X:Z,UP)**2, DIM=2)

       gamma(1:vl,DNDN) = SUM(grarho(1:vl,X:Z,DN)**2, DIM=2)

       gamma(1:vl,UPDN) = SUM(grarho(1:vl,X:Z,UP)*grarho(1:vl,X:Z,DN), DIM=2)
    endif

  end subroutine density_calc_ph_v2

  subroutine orbNGra( nuc_grads, eigvec, gphi, i_m )
    !
    ! Compute gradients of the molecular orbitals
    ! wrt symmetric modes:
    !
    !             gphi = - d/dN phi
    !
    use orbitalstore_module, only: orbital_nuclear_gradient_type
    use unique_atom_module, only: unique_atoms
    implicit none
    type(orbital_nuclear_gradient_type), intent(in)  :: nuc_grads(:) ! (nua)%o(:vl,3,nea,orb_dim,n_partn)
    real(r8_kind)   , intent(in)  :: eigvec(:,:)   ! (eig_dim,eig_dim)
    real(r8_kind)   , intent(out) :: gphi(:,:,:,:) ! (vl,eig_dim,n_part,n_modes)
    integer(i4_kind), intent(in)  :: i_m ! grid-hosting ua
    ! *** end of interface ***

    integer(I4_kind) :: vl,eig_dim,n_part,n_cart
    integer(I4_kind) :: n_modes
    integer(I4_kind) :: istat
    real(r8_kind), allocatable :: xphi(:,:,:,:) ! (vl,eig_dim,n_part,n_cart)

    ! number of cartesian gradients:
    n_cart   = 3 * sum( unique_atoms(:)%n_equal_atoms )

    eig_dim  = size(eigvec,1)

    ASSERT(eig_dim==size(gphi,2))
    vl       = size(gphi,1)
    n_part   = size(gphi,3)
    n_modes  = size(gphi,4)

    allocate( xphi(vl,eig_dim,n_part,n_cart), stat=istat )
    ASSERT(istat==0)

    ! compute cartesian gradients of eigenfunctions:
    call eigNGra( nuc_grads  &
                , eigvec     &
                , xphi       &
                )

    ! MAGIC: fix translational invariance of eigenfunction gradients:
    call fix_nuc_grads(i_m,xphi)

    ! transform cartesian gradients into mode-gradients:
    call symOrbNGra(vl,xphi,gphi)

    deallocate( xphi, stat=istat )
    ASSERT(istat==0)
  end subroutine orbNGra

  subroutine eigNGra(nuc_grads,eigvec,gphi)
    !
    ! transforms gradients wrt nuclear position to MO basis.
    !
    use orbitalstore_module, only: orbital_nuclear_gradient_type
    use cpks_grid_utils, only: eigMO
    implicit none
    type(orbital_nuclear_gradient_type), intent(in)  :: nuc_grads(:) ! (nua)%o(:vl,3,nea,orb_dim,n_partn)
    real(r8_kind)                      , intent(in)  :: eigvec(:,:)  ! (eig_dim,eig_dim)
    real(r8_kind)                      , intent(out) :: gphi(:,:,:,:) ! (=vl,eig_dim,n_partn,n_modes)
    ! *** end of interface ***

    integer(i4_kind) :: vl
    integer(i4_kind) :: i_ua, l, k
    integer(i4_kind) :: n_ea, e
    integer(i4_kind) :: orb_ind, gra_ind
    integer(i4_kind) :: orb_dim
    integer(i4_kind) :: eig_dim, n_modes

    ! first to MO-basis, then symmetrize (because of magic)

    vl = size(gphi,1)
    eig_dim = size(eigvec,2)
    n_modes = size(gphi,4)

    ! cumulative orbital/mode counter:
    orb_ind = 0
    gra_ind = 0
    gphi    = 0.0_r8_kind
    do i_ua=1,size(nuc_grads) !n_unique_atoms

      ! number of orbitlas of symmetry i_irr, localized on i_ua:
      orb_dim = size(nuc_grads(i_ua)%o,4)
      n_ea    = size(nuc_grads(i_ua)%o,3)

      do e=1,n_ea
        do l=1,size(nuc_grads(i_ua)%o,5) ! num of partners

        ! transforms gradients to MO basis:
        ! (basis functions range from orb_ind+1 to orb_ind+orb_dim)
          do k=1,3
            ! FIXME: avoid copy-in-out:
            call eigMO(vl, orb_ind                   &
                      , nuc_grads(i_ua)%o(:,k,e,:,l) &
                      , eigvec                       &
                      , gphi(:,:,l,gra_ind+k)        &
                      )
          enddo
        enddo
        ! increment cumulative counter:
        gra_ind = gra_ind + 3
      enddo

      ! increment cumulative counters:
      orb_ind = orb_ind + orb_dim
    enddo
    ASSERT(gra_ind==n_modes)
    ASSERT(orb_ind==eig_dim)
  end subroutine eigNGra

  subroutine symOrbNGra(vl,x,q)
    !
    ! transforms cartesin grads x(i_cart) to symmetric grads q(i_grad)
    !
    use unique_atom_module, only: unique_atom_grad_info
    implicit none
    integer(i4_kind), intent(in)  :: vl
    real(r8_kind)   , intent(in)  :: x(:,:,:,:) ! (:vl,n_orb,n_part,n_cart)
    real(r8_kind)   , intent(out) :: q(:,:,:,:) ! (:vl,n_orb,n_part,n_mode)
    ! *** end of interface ***

    real(r8_kind), pointer :: rot(:,:,:) ! => (ng,3,n_ea)
    integer(i4_kind) :: ua,ea
    integer(i4_kind) :: i_cart,i_grad,k_cart,n_grad,k,n_ea

    i_cart = 0
    i_grad = 0
    do ua=1,size(unique_atom_grad_info)

      rot => unique_atom_grad_info(ua)%m!(ng,1:3,n_ea)
      n_grad = size(rot,1)
      n_ea   = size(rot,3)

      do k=1,n_grad
        i_grad = i_grad + 1
        q(:vl,:,:,i_grad) = 0.0_r8_kind

        k_cart = i_cart
        do ea=1,size(rot,3)
          q(:vl,:,:,i_grad) = q(:vl,:,:,i_grad)                 &
                            + x(:vl,:,:,k_cart+1) * rot(k,1,ea) &
                            + x(:vl,:,:,k_cart+2) * rot(k,2,ea) &
                            + x(:vl,:,:,k_cart+3) * rot(k,3,ea)
          k_cart = k_cart + 3
        enddo
      enddo
      i_cart = i_cart + 3 * n_ea
    enddo
    ASSERT(i_grad==size(q,4))
    ASSERT(i_cart==size(x,4))
  end subroutine symOrbNGra

  subroutine fix_nuc_grads(i_m,gphi)
    !
    ! fix translational invariance of the eigenfunction gradients
    !
    use unique_atom_module, only: unique_atoms
    implicit none
    integer(i4_kind), intent(in)    :: i_m           ! grid hosting ua
    real(r8_kind)   , intent(inout) :: gphi(:,:,:,:) !(vl,eig_dim,partners(i_irr),n_cart)
    ! *** end of interface ***

    integer(i4_kind) :: o,a,n

    ! number of cartesian coordiantes:
    n = 3 * sum( unique_atoms(:)%n_equal_atoms )
    ASSERT(n==size(gphi,4))

    ! offset in cartesian coordiantes:
    o = 3 * sum( unique_atoms(:i_m-1)%n_equal_atoms )

    ! clear first ea of i_m ua:
    gphi(:,:,:,o+1) = 0.0_r8_kind
    gphi(:,:,:,o+2) = 0.0_r8_kind
    gphi(:,:,:,o+3) = 0.0_r8_kind

    ! loop over cartesian triplets:
    do a=0,n-3,3
        if( a == o ) cycle
        gphi(:,:,:,o+1) = gphi(:,:,:,o+1) - gphi(:,:,:,a+1)
        gphi(:,:,:,o+2) = gphi(:,:,:,o+2) - gphi(:,:,:,a+2)
        gphi(:,:,:,o+3) = gphi(:,:,:,o+3) - gphi(:,:,:,a+3)
    enddo
  end subroutine fix_nuc_grads

#if 1
  subroutine fix_rho_grads(i_m,vl,nuc_grarho,nuc_sec_derrho)
    !!! make nuc_grarho and nuc_sec_derrho correct
    use datatype, only: arrmat4, arrmat5
    use unique_atom_module, only: n_unique_atoms, unique_atoms
    implicit none
    integer(i4_kind), intent(in)    :: i_m ! UA hosting the grid
    integer(i4_kind), intent(in)    :: vl  ! grid batch size
    type(arrmat4)   , intent(inout) :: nuc_grarho(:)     ! (n_ua)%(:vla,3,n_ea,spin)
    type(arrmat5)   , intent(inout) :: nuc_sec_derrho(:) ! (n_ua)%(:vla,3,3,n_ea,spin)
    optional :: nuc_sec_derrho
    ! *** end of interface ***

    integer(i4_kind) :: spin,i_ua,i_ea
    integer(i4_kind) :: ispin

    ispin = size(nuc_grarho(i_m)%m,4)

    do spin=1,ispin
      nuc_grarho(i_m)%m(:vl,:3,1,spin)=0.0_r8_kind
      do i_ua=1,n_unique_atoms
       do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
        if(i_ua.eq.i_m.and.i_ea.eq.1) cycle
        nuc_grarho(i_m )%m(:vl,:,1,spin)= &
        nuc_grarho(i_m )%m(:vl,:,1,spin)  &
       -nuc_grarho(i_ua)%m(:vl,:,i_ea,spin)
       enddo
      enddo
    enddo

    if( .not. present(nuc_sec_derrho) ) RETURN

    do spin=1,ispin
      nuc_sec_derrho(i_m)%m(:vl,:3,:3,1,spin)=0.0_r8_kind
      do i_ua=1,n_unique_atoms
       do i_ea=1,unique_atoms(i_ua)%n_equal_atoms
        if(i_ua.eq.i_m.and.i_ea.eq.1) cycle
        nuc_sec_derrho(i_m )%m(:vl,:3,:3,1,spin)= &
        nuc_sec_derrho(i_m )%m(:vl,:3,:3,1,spin)  &
       -nuc_sec_derrho(i_ua)%m(:vl,:3,:3,i_ea,spin)
       enddo
      enddo
    enddo
  end subroutine fix_rho_grads
#endif

  !--------------- End of module ----------------------------------
end module cpks_dens_calc
