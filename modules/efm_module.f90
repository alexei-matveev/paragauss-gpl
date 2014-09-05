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
module  efm_module
  !-------------------------------------------------------------------
  !
  !  Purpose: contains subroutines which perform the Eigen Function
  !           Method (EFM)
  !
  !  References: Group Representation Theory for Physicists
  !              Jin-Quan Chen, World Scientific Singapore 1989
  !
  !  Author: MM
  !  Date: 10/96
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

#include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind,&
       & CK=>c16_kind  ! type specification parameters
  use efm_decl         ! efm-types, some methods and vars
#ifdef FPP_AIX_XLF
  use matrix_module, only: matmult
# define MATMUL(a,b) matmult(a,b)
!#define MATMUL(a,b) 1.0
#endif
  implicit none
  save            ! keep contents of module
  private         ! by default, all names are private


!== Interrupt end of public interface of module =================

  ! ------------ Declaration of types ------------------------------

  ! make type declarations in efm_decl public:
  public &
       & efm_cscoII_eigenspace,&
       & efm_csco_eigenspace,&
       & efm_cscoII_diag_invspace


  !------------ public functions and subroutines ---------------------
  public efm_trafo_gen,&
       & efm_proj_trafo_gen,&
       & efm_irrep_gen,&
       & efm_proj_irrep_gen,&
       & efm_cg_gen,&
       & efm_proj_cg_gen,&
       & efm_cscoII_solve,&
       & efm_proj_cscoII_solve,&
       & efm_proj_irrep_label,&
       & efm_done

#if FPP_COMPILE_UNUSED_CODE
  public :: efm_test_ylm
  public :: efm_test_ch
  public :: efmm_check_proj_irrep
#endif


  ! make some routines from efm_decl public:
  public &
       & efm_alloc,&
       & efm_free,&
       & efm_direct_product,&
       & efm_show_invspace,&
       & efm_rephase,&
       & efm_gauge,&
       & extract_basis,  &
       & load_basis,     &
       & phase

  !------------ Declaration of constants and variables ---------------
  ! the following factors are used to contruct the efmII and efmIII matrix
  ! they are chosen in order to obtain a optimal conditioned efm-matrix
  real(RK), public                           :: efm_subgroup_factor
  ! factor for constructing the efmII using subgroups
  real(RK), public                           :: efm_group_factor,efm_intrgroup_factor
  ! factor which weighs the efm-matrices of the group and the intrinsic group
  real(RK), public                           :: efm_proj_subgroup_factor
  ! factor for constructing the efmII using subgroups
  real(RK), public                           :: efm_proj_group_factor,efm_proj_intrgroup_factor
  ! factor which weighs the efm-matrices of the group and the intrinsic group
  !------------ Declaration of variables --------------------------
  type(efm_cscoII_diag_invspace),pointer,public          :: efm_cg(:,:)
  ! glebsch gordan coefficients
  ! efm_cg(group_num_ir,group_num_ir)
  type(efm_cscoII_diag_invspace),pointer,public          :: efm_proj_cg(:,:)
  ! clebsch gordan coefficients for the combination of
  ! vector and projective irreps
  ! efm_cg(group_num_ir,group_num_re)


  !===================================================================
  ! End of public interface of module
  !===================================================================


contains

  subroutine efm_done(spin_orbit)
    implicit none
    logical,optional,intent(in) :: spin_orbit
    ! *** end of interface ***

    logical     :: spor
    integer(IK) :: i,j
    integer(IK) :: memstat

    DPRINT 'efm/efm_done: entered'

    spor = present(spin_orbit)
    if(present(spin_orbit)) spor = spin_orbit

    DPRINT 'efm/efm_done: free efm_cg:'
    do i=1,size(efm_cg,1)
       do j=1,size(efm_cg,2)
          call efm_free(efm_cg(i,j))
       enddo
    enddo
    deallocate(efm_cg,STAT=memstat)
    ASSERT(memstat==0)

    DPRINT 'efm/efm_done: free efm_proj_cg:'
    if(spor)then
       do i=1,size(efm_proj_cg,1)
          do j=1,size(efm_proj_cg,2)
             call efm_free(efm_proj_cg(i,j))
          enddo
       enddo
       deallocate(efm_proj_cg,STAT=memstat)
       ASSERT(memstat==0)
    endif

    DPRINT 'efm/efm_done: call efm_irrep_dealloc():'
    call efm_irrep_dealloc()

    DPRINT 'efm/efm_done: exit'
  end subroutine efm_done

  !*************************************************************
  subroutine efm_irrep_gen
    !  Purpose: generates and allocates all irreps
    !           of the group
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    !------------ Declaration of local variables ---------------------
    real(RK),allocatable             :: reg_trafo(:,:,:)
    ! trafo matrix of the regular representation for the  group
    real(RK),allocatable             :: reg_trafo_intr(:,:,:)
    ! trafo matrix of the regular representation for the intrinsic group
    !------------ Executable code ------------------------------------

    ! allocate trafo matrices of the regular representation
    allocate(reg_trafo(group_num_el,group_num_el,group_num_el))
    allocate(reg_trafo_intr(group_num_el,group_num_el,group_num_el))

    ! generate regular representation
    call efm_regrep_gen(reg_trafo,reg_trafo_intr)

    ! generate and allocate irrep matrices
    call efm_cscoIII_solve(reg_trafo,reg_trafo_intr)

    ! deallocate arrays
    deallocate(reg_trafo,reg_trafo_intr)

  end subroutine efm_irrep_gen
  !*************************************************************

  !*************************************************************
  subroutine efm_proj_irrep_gen
    !  Purpose: generates and allocates all projective irreps
    !           of the group
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    !------------ Declaration of local variables ---------------------
    real(RK),allocatable             :: reg_trafo(:,:,:)
    ! trafo matrix of the regular representation for the  group
    real(RK),allocatable             :: reg_trafo_intr(:,:,:)
    ! trafo matrix of the regular representation for the intrinsic group

    !------------ Executable code ------------------------------------

    ! allocate trafo matrices of the regular representation
    allocate(reg_trafo(group_num_el,group_num_el,group_num_el))
    allocate(reg_trafo_intr(group_num_el,group_num_el,group_num_el))

    ! generate regular representation
    call efm_proj_regrep_gen(reg_trafo,reg_trafo_intr)

    ! generate and allocate irrep matrices
    call efm_proj_cscoIII_solve(reg_trafo,reg_trafo_intr)

    ! deallocate arrays
    deallocate(reg_trafo,reg_trafo_intr)

  end subroutine efm_proj_irrep_gen
  !*************************************************************

  !*************************************************************
  subroutine efm_irrep_dealloc
    !  Purpose: deallocate irrep matrices
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    !------------ Declaration of local variables ---------------------
    integer(IK)                      :: i
    ! counter

    !------------ Executable code ------------------------------------

    do i=1,group_num_ir
       deallocate(irrep_can(i)%irrep_gen_coeff)
       deallocate(irrep_can(i)%irrep_gen_oper)
       deallocate(irrep_can(i)%irrep_matrix)
       deallocate(irrep_can(i)%characters)
    end do

    deallocate(irrep_can)

  end subroutine efm_irrep_dealloc
  !*************************************************************


#if FPP_COMPILE_UNUSED_CODE
  !*************************************************************
  subroutine efm_test_ylm
    !  Purpose: try efm_ylmrot
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    !------------ Declaration of local variables ---------------------
    integer(IK)                      :: l,nls
    ! angular momentum
    real(RK),allocatable             :: trafo_matrix(:,:,:)
    ! transformation matrices of all group elements
    ! dimensions: trafo_matrix(group_num_el,2*l+1,2*l+1)
    type(efm_cscoII_diag_invspace)             :: so3_invspace
    ! irreducible eigenspace of SO(3) = invariant space of group

    l = 2
    nls = 2*l+1

    ! allocate transformation matrix
    allocate(trafo_matrix(nls,nls,group_num_el))

    ! allocate eigensolutions
    call efm_alloc(group_num_ir,so3_invspace,DIM=nls)

    ! generate transformation matrices
    call efm_trafo_gen(l,trafo_matrix)

    ! do efm
    call efm_cscoII_solve(trafo_matrix,nls,so3_invspace)

    ! deallocate arrays
    deallocate(trafo_matrix)
    call efm_free(so3_invspace)

  end subroutine efm_test_ylm
  !*************************************************************
#endif


  !*************************************************************
  subroutine efm_cg_gen
    !  Purpose: generates the glebsch gordan coefficients
    !** End of interface *****************************************
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    use group_module, n_irr => group_num_ir
    use clebsch_gordan
    !------------ Declaration of local variables ---------------------
    integer(IK)                      :: i_irrep,j_irrep,k_irrep,i,j
    ! counter
    integer(IK)                      :: irrep_dim_i,irrep_dim_j
    ! dimension of the irreps i and j
    integer(IK)                      :: alloc_stat
    ! allocation status
    integer(IK)                      :: product_dim
    ! dimension of the direct product space
    real(RK),allocatable             :: trafo_matrix(:,:,:)
    ! transformation matrix of the direct product

    !------------ Executable code ------------------------------------

    ! allocate cg-coefficients
    allocate(efm_cg(n_irr,n_irr),stat=alloc_stat)
    ASSERT(alloc_stat==0)

    allocate(cg(n_irr,n_irr,n_irr),stat=alloc_stat)
    ASSERT(alloc_stat==0)

    ! loop over all combinations of irreps
    do i_irrep=1,n_irr
       irrep_dim_i = irrep_can(i_irrep)%dimension
       do j_irrep=1,n_irr
          irrep_dim_j = irrep_can(j_irrep)%dimension

          ! determine dimension of product space
          product_dim = irrep_dim_i*irrep_dim_j

          call efm_alloc(n_irr,efm_cg(i_irrep,j_irrep),dim=product_dim)

          ! allocate transformation matrix
          allocate(trafo_matrix(product_dim,product_dim,group_num_el),&
               STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("efm_cg_gen: temporary allocation failed")

          ! generate the transformation matrices of the direct product
          do i=1,group_num_el
             call efm_direct_product(&
                  & irrep_can(i_irrep)%irrep_matrix(:,:,i),&
                  & irrep_can(j_irrep)%irrep_matrix(:,:,i),&
                  & trafo_matrix(:,:,i)&
                  & )
          end do

          ! symmetrize direct product
          call efm_cscoII_solve(trafo_matrix,product_dim,efm_cg(i_irrep,j_irrep))

          do k_irrep=1,n_irr
             call load_efm_to_cg(&
                  & irrep_dim_i,irrep_dim_j,&
                  & efm_cg(i_irrep,j_irrep)%csco_eigenspaces(k_irrep),&
                  & cg(k_irrep,i_irrep,j_irrep)&
                  & )
          enddo

          ! deallocate transformation matrix
          deallocate(trafo_matrix,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       end do
    end do

    do k_irrep=1,n_irr
    do i_irrep=1,n_irr
    do j_irrep=1,i_irrep-1
      do i=1,cg(k_irrep,i_irrep,j_irrep)%mult
      do j=1,size(cg(k_irrep,i_irrep,j_irrep)%sub(i)%c,1)
             cg(k_irrep,i_irrep,j_irrep)%sub(i)%c(j,:,:) =  &
             transpose(                                     &
             cg(k_irrep,j_irrep,i_irrep)%sub(i)%c(j,:,:)    &
                      )
      enddo
      enddo
    enddo
    enddo
    enddo

!!    call show_cg(cg)
  end subroutine efm_cg_gen
  !*************************************************************

  !*************************************************************
  subroutine efm_proj_cg_gen
    !  Purpose: generates the glebsch gordan coefficients
    !           for projective irreps
    !** End of interface *****************************************
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    use group_module, n_virr => group_num_ir, n_pirr => group_num_re
    use clebsch_gordan
    !------------ Declaration of local variables ---------------------
    integer(IK)                      :: i_irrep,j_irrep,k_irrep,i
    ! counter
    integer(IK)                      :: irrep_dim_i,irrep_dim_j
    ! dimension of the irreps i and j
    integer(IK)                      :: alloc_stat
    ! allocation status
    integer(IK)                      :: product_dim
    ! dimension of the direct product space
    complex(CK),allocatable         :: trafo_matrix(:,:,:)
    ! transformation matrix of the direct product

    complex(CK),allocatable         :: trafo_su2(:,:,:)
    ! transformation matrices of su(2)
    type(efm_cscoII_diag_invspace)  :: invspace
    ! for reduction of product space

    !------------ Executable code ------------------------------------

    ! allocate cg-coefficients
    allocate(efm_proj_cg(n_virr,n_pirr))

    ! allocate vpcg struct:
    allocate(vpcg(n_pirr,n_virr,n_pirr))


    ! loop over all combinations of vector irreps and projective irreps
    do i_irrep=1,n_virr ! loop over vector irreps
       irrep_dim_i = irrep_can(i_irrep)%dimension
       do j_irrep=1,n_pirr ! loop over projective irreps

          DPRINT "    COMBINATION:"
          DPRINT "   ",irrep_can(i_irrep)%label," x ",proj_irrep_can(j_irrep)%label

          irrep_dim_j = proj_irrep_can(j_irrep)%dimension
          ! determine dimension of product space
          product_dim = irrep_dim_i*irrep_dim_j

          call efm_alloc(n_pirr,efm_proj_cg(i_irrep,j_irrep),DIM=product_dim,PROJ=.true.)

          ! allocate transformation matrix
          allocate(trafo_matrix(product_dim,product_dim,group_num_el),&
               STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("efm_cg_gen: temporary allocation failed")

          ! generate the transformation matrices of the direct product
          do i=1,group_num_el
             call efm_direct_product(&
                  & irrep_can(i_irrep)%irrep_matrix(:,:,i),&
                  & proj_irrep_can(j_irrep)%irrep_matrix(:,:,i),&
                  & trafo_matrix(:,:,i))
          end do

          ! symmetrize direct product
          call efm_proj_cscoII_solve(trafo_matrix,product_dim,efm_proj_cg(i_irrep,j_irrep))

          do k_irrep=1,n_pirr
             call load_efm_proj_to_cg(&
                  & irrep_dim_i,irrep_dim_j,&
                  & efm_proj_cg(i_irrep,j_irrep)%csco_eigenspaces(k_irrep),&
                  & vpcg(k_irrep,i_irrep,j_irrep)&
                  & )
          enddo

          ! deallocate transformation matrix
          deallocate(trafo_matrix)
       end do
    end do

    !
    ! now direct products of vector irreps with su(2) irrep:
    !
    allocate(vsu2cg(n_pirr,n_virr)) ! from cg-module
    allocate(trafo_su2(2,2,group_num_el))

    ! generate rotation matrices:
    do i=1,group_num_el
       call efm_so3_irrep(symm_el(i)%quaternion,1,trafo_su2(:,:,i))
    enddo

    do i_irrep=1,n_virr ! loop over vector irreps
       irrep_dim_i = irrep_can(i_irrep)%dimension

       DPRINT "    COMBINATION:"
       DPRINT "   ",irrep_can(i_irrep)%label," x SU(2)"

       product_dim = 2 * irrep_dim_i

       call efm_alloc(n_pirr,invspace,DIM=product_dim,PROJ=.true.)

       ! allocate transformation matrix
       allocate(trafo_matrix(product_dim,product_dim,group_num_el))

       ! generate the transformation matrices of the direct product
       do i=1,group_num_el
          call efm_direct_product(&
               & irrep_can(i_irrep)%irrep_matrix(:,:,i),&
               & trafo_su2(:,:,i),&
               & trafo_matrix(:,:,i))
       end do

       ! symmetrize direct product
       call efm_proj_cscoII_solve(trafo_matrix,product_dim,invspace)

       do k_irrep=1,n_pirr
          call load_efm_proj_to_cg(&
               & irrep_dim_i,2,&
               & invspace%csco_eigenspaces(k_irrep),&
               & vsu2cg(k_irrep,i_irrep)&
               & )
       enddo

       ! deallocate transformation matrix
       deallocate(trafo_matrix)

       ! clean invspace:
       call efm_free(invspace)
    enddo

    deallocate(trafo_su2)

#ifdef FPP_DEBUG
    call show_vsu2cg(vsu2cg)
#endif

  end subroutine efm_proj_cg_gen
  !*************************************************************


  !*************************************************************
  subroutine efm_proj_irrep_label
    !  Purpose: label the projective irreps which are prelabled by
    !           group_proj_irrep_pre_label
    use type_module
    use iounitadmin_module
    use group_module
    !------------ Declaration of formal parameters ---------------
    !------------ Declaration of local variables ---------------------
    type(efm_cscoII_diag_invspace)     :: ylm_adapt
    ! eigenspace of CSCO-II
    type(group_proj_irrep),pointer     :: proj_irrep_can_imed(:)
    integer(IK)              :: two_j,l,s,nls,i,label_point,dim_ind,n_irreps,up_dn
    integer(IK)              :: n_irreps_appeared,j,m
    ! number of irreps appeared in investigated spherical harmonics
    integer(IK),allocatable  :: irrep_appeared_first(:)
    ! 2j of spherical harmonic where irreps appeared at first
    integer(IK),allocatable  :: irrep_sequence(:)
    ! sequence of appearance
    logical                            :: pseudo
    !------------ Executable code ------------------------------------

    ! allocate and initialize indicator if irrep appeared
    allocate(irrep_appeared_first(group_num_pir),irrep_sequence(group_num_pir))
    irrep_appeared_first = 0
    irrep_sequence = 0
    proj_irrep_can_imed => proj_irrep_can

    if (inversion.ne.0) then
       n_irreps = group_num_pir/2
    else
       n_irreps = group_num_pir
    endif

    l = 0
    n_irreps_appeared = 0
    ! loop over gerade sperical harmonics
    do
       do s=1,2
          if ((l.eq.0).and.(s.eq.1)) then
             cycle
          endif
          two_j = 2*l + 2*s - 3
          nls = two_j + 1

          if (output_unit > 0) then
             write(output_unit,*) "l = ",l
             write(output_unit,*) "j = ",two_j,"/2"
             write(output_unit,*) "dimension of space: ",nls
          endif
          ! determine CSCO-II eigenspaces of sperical harmonics
          allocate(ylm_adapt%csco_eigenspaces(group_num_pir))
          call efm_proj_cscoII_solve(ylm_proj_trafos(l,s)%matrix,nls,ylm_adapt)
          do i=1,group_num_pir
             if (irrep_appeared_first(i).eq.0) then
                if (ylm_adapt%csco_eigenspaces(i)%exists) then
                   n_irreps_appeared = n_irreps_appeared + 1
                   if (output_unit > 0) then
                      write(output_unit,*) " identified irrep # ",i
                      write(output_unit,*) "n_irreps_appeared: ",n_irreps_appeared
                   endif
                   irrep_appeared_first(i) = two_j
                   irrep_sequence(n_irreps_appeared) = i
                endif
             endif
          enddo
          call efm_free(ylm_adapt)
       enddo
       if (n_irreps_appeared.eq.n_irreps) then
          exit
       endif
       l = l + 2
    enddo

    if (inversion.ne.0) then
       l = 1
       ! loop over ungerade sperical harmonics
       n_irreps_appeared = 0
       do
          do s=1,2
             two_j = 2*l + 2*s - 3
             nls = two_j + 1
             if (output_unit > 0) then
                write(output_unit,*) "l = ",l
                write(output_unit,*) "j = ",two_j,"/2"
                write(output_unit,*) "dimension of space: ",nls
             endif
             ! determine CSCO-II eigenspaces of sperical harmonics
             allocate(ylm_adapt%csco_eigenspaces(group_num_pir))
             call efm_proj_cscoII_solve(ylm_proj_trafos(l,s)%matrix,nls,ylm_adapt)
             do i=1,group_num_pir
                if (irrep_appeared_first(i).eq.0) then
                   if (ylm_adapt%csco_eigenspaces(i)%exists) then
                      n_irreps_appeared = n_irreps_appeared + 1
                      if (output_unit > 0) then
                         write(output_unit,*) " identified irrep # ",i
                         write(output_unit,*) "n_irreps_appeared: ",n_irreps_appeared
                      endif
                      irrep_appeared_first(i) = two_j
                      irrep_sequence(n_irreps_appeared + n_irreps) = i
                   endif
                endif
             enddo
             call efm_free(ylm_adapt)
          enddo
          if (n_irreps_appeared.eq.n_irreps) then
             exit
          endif
          l = l + 2
       enddo
    endif

    ! allocate new irreps
    allocate(proj_irrep_can(group_num_re))
    do i=1,group_num_pir
       ! FIXME: memory leak?
       allocate(proj_irrep_can(i)%characters(group_num_re))
    enddo

    ! now label irreps
    up_dn = 1
    do i=1,group_num_pir
       if (output_unit > 0) then
          write(output_unit,*) "efm_proj_irrep_label: label irrep ",i
          write(output_unit,*) "efm_proj_irrep_label: irrep_sequence ",irrep_sequence(i)
       endif
       label_point = 1
       dim_ind = proj_irrep_can_imed(irrep_sequence(i))%dimension
       pseudo = proj_irrep_can_imed(irrep_sequence(i))%pseudo

       ! attach irrep properties
       proj_irrep_can(i) = proj_irrep_can_imed(irrep_sequence(i))

       ! FIXME: arent they set by an assignment above?
       proj_irrep_can(i)%time_dimension = dim_ind
       proj_irrep_can(i)%label = "        "

       ! calculate 2*j
       two_j = irrep_appeared_first(irrep_sequence(i))

       ! discriminate "partners" of pseudo-2D irreps
       if (pseudo) then
          dim_ind = 2*dim_ind
          if (proj_irrep_can(i)%sign_of_jz.gt.0.0_rk) then
             proj_irrep_can(i)%label(label_point:label_point) = '1'
          else
             proj_irrep_can(i)%label(label_point:label_point) = '2'
          endif
          label_point = label_point + 1
          proj_irrep_can(i)%sign_of_jz = proj_irrep_can(i)%sign_of_jz*two_j*0.5_rk
       endif

       select case(dim_ind)
       case(1)
          proj_irrep_can(i)%label(label_point:label_point) = 'A'
       case(2)
          proj_irrep_can(i)%label(label_point:label_point) = 'E'
       case(3)
          proj_irrep_can(i)%label(label_point:label_point) = 'T'
       case(4)
          proj_irrep_can(i)%label(label_point:label_point) = 'G'
       case(5)
          proj_irrep_can(i)%label(label_point:label_point) = 'H'
       case(6)
          proj_irrep_can(i)%label(label_point:label_point) = 'I'
       end select
       label_point = label_point + 1
       !am: print*,"two_j",two_j
       if (two_j.ge.10) then
          write(proj_irrep_can(i)%label(label_point:label_point+4),'(I2,a2)') two_j,'/2'
          label_point = label_point + 4
       else
          write(proj_irrep_can(i)%label(label_point:label_point+3),'(I1,a2)') two_j,'/2'
          label_point = label_point + 3
       endif
       if (inversion.ne.0) then
          if (i.gt.n_irreps) then
             proj_irrep_can(i)%label(label_point:label_point) = 'U'
          else
             proj_irrep_can(i)%label(label_point:label_point) = 'G'
          endif
       endif
    enddo
    if (output_unit > 0) then
       write(output_unit,*) ">>> Final Labeling of Irreps <<<"
       write(output_unit,*) "found ",group_num_re," regular classes and Irreps"
       write(output_unit,601) group_name, &
            (symm_el(klass(regular(j)) % conjug_el(1)) % name, j = 1, group_num_re)
       do j = 1,group_num_pir
          write(output_unit,602) proj_irrep_can(j) % label, &
               (proj_irrep_can(j) % characters(m), m = 1, group_num_re)
       end do
    endif

601 format(/1x,'Group ',a,' eigenvalues of class operator ',&
         &   'and primitive characters :'/ 5(10x,4(a8,8x :)/))
602 format(1x,A8,2x,4(f8.3,f8.3:)/4(10x,4(f8.3,f8.3:)/))

    deallocate(irrep_appeared_first)

  end subroutine efm_proj_irrep_label
  !*************************************************************

  !*************************************************************
  subroutine efm_regrep_gen(trafo_matrix,trafo_matrix_intr)
    !  Purpose: generates the transformation matrices of the
    !           regular representation of the group and the
    !           intrinisic group
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    !------------ Declaration of formal parameters ---------------
    real(RK),    intent(out)  :: trafo_matrix(:,:,:)
    ! transformation of the regular representation
    real(RK),    intent(out)  :: trafo_matrix_intr(:,:,:)
    ! transformation of the regular representation for the intrinsic
    ! group
    !** End of interface *****************************************

    !------------ Declaration of local variables ---------------------
    integer(IK)                :: i,j
    ! counter
    !------------ Executable code ------------------------------------

    ! initialize trafo matrices to zero, as most of their elements
    ! will be
    trafo_matrix = 0.0_rk
    trafo_matrix_intr = 0.0_rk
    ! the regular representation matrices are determined
    ! according to the formulas:
    ! D(R)k,l = delta(k,Rl) for the rep of the group and
    ! Dintr(R)k,l = delta(k,lR) for the rep of the intrinsic group
    do i=1,group_num_el
       do j=1,group_num_el
          trafo_matrix(group_multab(i,j),j,i) = 1.0_rk
          trafo_matrix_intr(group_multab(j,i),j,i) = 1.0_rk
       end do
    end do

  end subroutine efm_regrep_gen
  !*************************************************************


  !*************************************************************
  subroutine efm_proj_regrep_gen(trafo_matrix,trafo_matrix_intr)
    !  Purpose: generates the transformation matrices of the
    !           regular projective representation of the group and the
    !           intrinisic group
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    !------------ Declaration of formal parameters ---------------
    real(RK),    intent(out)  :: trafo_matrix(:,:,:)
    ! transformation of the regular representation
    real(RK),    intent(out)  :: trafo_matrix_intr(:,:,:)
    ! transformation of the regular representation for the intrinsic
    ! group
    !** End of interface *****************************************

    !------------ Declaration of local variables ---------------------
    integer(IK)                :: i,j
    ! counter
    !------------ Executable code ------------------------------------

    ! initialize trafo matrices to zero, as most of their elements
    ! will be
    trafo_matrix = 0.0_rk
    trafo_matrix_intr = 0.0_rk
    ! the regular representation matrices are determined
    ! according to the formulas:
    ! D(R)k,l = delta(k,Rl)*[R,l] for the rep of the group and
    ! Dintr(R)k,l = delta(k,lR)*[l,R] for the rep of the intrinsic group
    do i=1,group_num_el
       do j=1,group_num_el
          trafo_matrix(group_multab(i,j),j,i) = group_factab(i,j)
          trafo_matrix_intr(group_multab(j,i),j,i) = group_factab(j,i)
       end do
    end do

  end subroutine efm_proj_regrep_gen
  !*************************************************************



  !*************************************************************
  subroutine efm_trafo_gen(l,trafo_matrix)
    !  Purpose: generates transformation matrices for a given
    !           angular momentum l
    !------------ Modules used ------------------- ---------------
    use iounitadmin_module
    use group_module
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in)   :: l
    ! angular momentum
    real(RK),    intent(out)  :: trafo_matrix(:,:,:)
    ! transformation matrices of all group elements
    ! dimensions: trafo_matrix(group_num_el,2*l+1,2*l+1)
    !** End of interface *****************************************


    !------------ Declaration of local variables ---------------------
    real(RK)                   :: euler(4),quaternion(5)
    ! euler angles and quaternions
    integer(IK)                :: nls
    ! angular momentum
    integer(IK)                :: kl
    ! counter
    real(RK),allocatable       :: dmat(:,:)
    ! rotation matrix to be determined
    real(RK),allocatable       :: dbeta(:,:)
    ! auxiliary array for old version of ylmrot

    !------------ Executable code ------------------------------------

    nls = 2*l+1

    ! allocate rotation matrix
    allocate(dmat(2*l+1,2*l+1))

    ! allocate auxiliary array
    allocate(dbeta(2*l+1,2*l+1))

    ! test of rotation
    euler(1) = 0.6_rk
    euler(2) = 0.0000001_rk
    euler(3) = 0.0_rk
    euler(4) = 1.0_rk
    call efm_ylmrot(l,nls,nls,euler(1),euler(2),euler(3),dmat,&
         &euler(4),dbeta)

    do kl = 1,group_num_el

       ! generate rotation matrices
       quaternion = symm_el(kl)%quaternion

       call efm_trafo_quat_euler(quaternion,euler)

       call efm_ylmrot(l,nls,nls,euler(1),euler(2),euler(3),dmat,&
            &euler(4),dbeta)

       ! read to trafo_matrix
       trafo_matrix(:,:,kl) = dmat

    end do

    ! deallocate arrays
    deallocate(dmat,dbeta)

  end subroutine efm_trafo_gen

  !*************************************************************
  subroutine efm_proj_trafo_gen(l,s,trafo_matrix,complex_to_real_up,complex_to_real_down)
    !  Purpose: generates transformation matrices for a given
    !           total angular momentum j = l + s - 3/2
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
!!$    use options_module, only: options_debug_key !<<<debug
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in)       :: l,s
    ! angular momentum
    complex(CK),    intent(out)  :: trafo_matrix(:,:,:)
    ! transformation matrices of all group elements
    ! dimensions: trafo_matrix(group_num_el,2*j+1,2*j+1)
    complex(CK),    intent(out)  :: complex_to_real_up(:,:),complex_to_real_down(:,:)
    ! transformation matrices of all group elements
    ! dimensions: trafo_matrix(group_num_el,2*l+1,2*l+1)
    !** End of interface *****************************************
    !------------ Declaration of local constants  -----------------
    real(RK),parameter         :: small = 1e-12_rk,zero=0.0_rk,one=1.0_rk,two=2.0_rk
    complex(CK),parameter     :: i = (0.0_rk,1.0_rk),czero = (0.0_rk,0.0_rk)
    ! imaginary one

    !------------ Declaration of local variables ---------------------
    integer(IK)                :: two_j,two_j_plus_one
    ! angular momentum
    integer(IK)                :: kl,m,m1,alloc_stat
    ! counter
    real(RK)                   :: quaternion(5),sqrt_2_inv,minus_exp_m,little_value,sqrt_2l_plus_1
    ! quaternion
    complex(CK),allocatable   :: U_matrix(:,:)
    ! transformation of real to complex spherical harmonics

    ! for debugging >>

    complex(CK),allocatable   :: trafo_matrix_real(:,:,:)
    ! transformation matrices of all group elements
    ! dimensions: trafo_matrix(group_num_el,2*j+1,2*j+1)
    complex(CK),allocatable   :: trafo_matrix_new(:,:,:)
    ! transformation matrices of all group elements
    complex(CK),allocatable   :: trafo_spinor(:,:,:)
    ! transformation matrices of spin=1/2 spinors
    complex(CK),allocatable   :: complex_to_real(:,:)
    ! transformation matrices of spin=1/2 spinors

    ! for debugging <<



    !------------ Executable code ------------------------------------

    allocate(U_matrix(2*l+1,2*l+1)&
         ,stat = alloc_stat)
    if (alloc_stat .ne. 0 ) call error_handler( &
         "efm_proj_trafo_gen: allocation of U_matrix failed")

    !
    ! determine transformation matrix U_matrix
    !
    ! determine inverse square root of two
    sqrt_2_inv = 1/(sqrt(2.0_rk))
    minus_exp_m = one
    U_matrix = czero
    U_matrix(l+1,l+1) = (one,zero) ! i
    !do m=1,l
    !   minus_exp_m = (-one)*minus_exp_m
    !   U_matrix(l+1-m,l+1-m) = cmplx(sqrt_2_inv,zero,RK)
    !   U_matrix(l+1-m,l+1+m) = i*sqrt_2_inv
    !   U_matrix(l+1+m,l+1-m) = cmplx(minus_exp_m*sqrt_2_inv,zero,RK)
    !   U_matrix(l+1+m,l+1+m) = -i*minus_exp_m*sqrt_2_inv
    !end do
    do m=1,l
       minus_exp_m = (-one)*minus_exp_m
       U_matrix(l+1-m,l+1+m) = -i*sqrt_2_inv
       U_matrix(l+1-m,l+1-m) = sqrt_2_inv
       U_matrix(l+1+m,l+1+m) = i*minus_exp_m*sqrt_2_inv
       U_matrix(l+1+m,l+1-m) = minus_exp_m*sqrt_2_inv
    end do
    !U_matrix = conjg(U_matrix)
    ! determine total angular momentum
    two_j = 2*l + 2*s - 3
    two_j_plus_one = two_j + 1
    sqrt_2l_plus_1 = one/sqrt(2*l + one)

    if (output_unit > 0) then
       write(output_unit,*) "----------------------"
       write(output_unit,*) "total angular momentum ",two_j,"/ 2"
       write(output_unit,*) "  angular momentum ",l," spin ",s
       write(output_unit,*) "----------------------"
    endif


    !allocate(intermediate(two_j_plus_one,2*l+1)&
    !     ,stat = alloc_stat)
    !if (alloc_stat .ne. 0 ) call error_handler( &
    !     "efm_proj_trafo_gen: allocation of intermediate failed")

    do kl = 1,group_num_el

       ! generate rotation matrices
       quaternion = symm_el(kl)%quaternion
       call efm_so3_irrep(quaternion,two_j,trafo_matrix(:,:,kl))

       if (NINT(quaternion(5)) == -1) then
          ! write(output_unit,*) "element ",symm_el(kl)%name," contains inversion"
          trafo_matrix(:,:,kl) = trafo_matrix(:,:,kl)*(-1.0_rk)**l
       end if

    end do

!!$    !mdf>>>
!!$    select case(options_debug_key(1))
!!$    case(0)
!!$    case(1)
!!$       print *,'efm/proj_trafo_gen: Ylm gauge changed to -1'
!!$       call efm_gauge(trafo_matrix,(-1.0_rk,0.0_rk))
!!$       ! negated, but still Pauli gauge
!!$    case(2)
!!$       print *,'efm/proj_trafo_gen: Ylm gauge changed to -i'
!!$       call efm_gauge(trafo_matrix,-i);
!!$       ! Inv == -i, Cartan gauge, must not work
!!$    end select
!!$    !<<<mdf

    U_matrix = conjg(U_matrix)

    if (s.eq.2) then
       ! spin up component
       complex_to_real_up(:,1) = zero
       do m1 = 2,two_j_plus_one
          complex_to_real_up(:,m1) = conjg(U_matrix(m1-1,:))*sqrt(m1-one)*sqrt_2l_plus_1
       end do

       !
       ! Out of boundary access was living in the next loop.
       ! Note that
       !
       !      two_j_plus_one == 2L + 2S - 2 == 2L + 2 here!
       !
       ! but U_matrix(2L+2, :) is out of bounds.
       ! "Fortunately", the square root of 2L+2-m1 is zero for
       !
       !      m1 == 2L+2
       !
       ! That is probably why this bug lived that long.
       !
       ASSERT(2*l+two-two_j_plus_one==0)
! BUG: do m1 = 1,two_j_plus_one
       do m1 = 1, two_j_plus_one - 1
          complex_to_real_down(:,m1) = conjg(U_matrix(m1,:))*sqrt(2*l+two-m1)*sqrt_2l_plus_1
       end do
       complex_to_real_down(:, two_j_plus_one) = zero
    else
       ! spin down component
       do m1 = 1,two_j_plus_one
          complex_to_real_up(:,m1) = -conjg(U_matrix(m1,:))*sqrt(2*l+one-m1)*sqrt_2l_plus_1
          complex_to_real_down(:,m1) = conjg(U_matrix(m1+1,:))*sqrt(m1+zero)*sqrt_2l_plus_1
       end do
    end if

    ! now test complex_to_real

    ! transformation matrix of the direct product of spin=1/2 and real sperical
    ! harmonics of l
    allocate(trafo_matrix_real(4*l+2,4*l+2,group_num_el),trafo_spinor(2,2,group_num_el),&
         complex_to_real(4*l+2,two_j_plus_one),trafo_matrix_new(two_j_plus_one,two_j_plus_one,group_num_el)&
         ,stat = alloc_stat)
    if (alloc_stat .ne. 0 ) call error_handler( &
         "efm_proj_trafo_gen: allocation of trafo_matrix_real failed")
    complex_to_real(1:2*l+1,:) = complex_to_real_down
    complex_to_real(2*l+2:4*l+2,:) = complex_to_real_up

    do kl = 1,group_num_el
       ! generate rotation matrices
       quaternion = symm_el(kl)%quaternion
       call efm_so3_irrep(quaternion,1,trafo_spinor(:,:,kl))

!!$       !mdf>>>
!!$       select case(options_debug_key(1))
!!$       case(0)
!!$       case(1)
!!$          call error_handler("efm/efm_proj_trafo_gen: options_debug_key no more supported")
!!$          call efm_gauge(kl,trafo_spinor(:,:,kl),-1.0_rk)
!!$       case(2)
!!$          call error_handler("efm/efm_proj_trafo_gen: options_debug_key no more supported")
!!$          call efm_gauge(kl,trafo_spinor(:,:,kl),-i)
!!$       end select
!!$       !<<<mdf

       call efm_direct_product(&
            & trafo_spinor(:,:,kl),&
            & ylm_trafos(l)%matrix(:,:,kl),&
            & trafo_matrix_real(:,:,kl))
       trafo_matrix_new(:,:,kl) =  matmul(conjg(transpose(complex_to_real)),matmul(trafo_matrix_real(:,:,kl),complex_to_real))

       little_value = maxval(abs(trafo_matrix_new(:,:,kl)-trafo_matrix(:,:,kl)))
       !write(output_unit,*) " difference of tranformation ",kl," :",little_value
       if (little_value.gt.small) then
          if (output_unit > 0) then
             write(output_unit,*) "Error: wrong tranformation between real and complex sherical harmonics)"
             write(output_unit,*) "       for symmetry element:",symm_el(kl)%name
          endif
          call error_handler("efm/efm_proj_trafo_gen: error occured see output_unit")
       end if
    end do



    deallocate(trafo_matrix_real,trafo_spinor,complex_to_real,trafo_matrix_new,stat = alloc_stat)
    if (alloc_stat .ne. 0 ) call error_handler( &
         "efm_proj_trafo_gen: deallocation of trafo_matrix_real failed")

    deallocate(U_matrix&
         ,stat = alloc_stat)
    if (alloc_stat .ne. 0 ) call error_handler( &
         "efm_proj_trafo_gen: deallocation of U_matrix failed")
  end subroutine efm_proj_trafo_gen
  !*************************************************************


  !*************************************************************
  subroutine efm_trafo_quat_euler(quaternion,euler)
    !  Purpose: Transforms quaternion parameters into Euler-
    !           angles
    !  according to Altmann: Rotations, Quaternions and Double
    !                        Groups ; section 12.11
    !------------ Modules used ------------------- ---------------
    use type_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(RK),    intent(in)  :: quaternion(5)
    ! quaternion parameters (including inversion)
    real(RK),    intent(out) :: euler(4)
    ! euler angles (alpha,beta,gamma,inversion)
    !** End of interface *****************************************

    !------------ Declaration of local parameters ----------------
    real(RK),parameter         :: small = 0.000001_rk
    ! euler angles

    !------------ Declaration of local variables ---------------------
    real(RK)                   :: alpha,beta,gamma
    ! euler angles
    real(RK)                   :: sinb2,cosb2
    ! sin and cos of beta/2
    real(RK)                   :: sinapg2,cosapg2,apg
    ! sin and cos of (alpha + gamma)/2
    real(RK)                   :: sinamg2,cosamg2,amg
    ! sin and cos of (alpha - gamma)/2


    !------------ Executable code ------------------------------------

    ! take over inversion symmetry
    euler(4) = quaternion(5)

    ! calculate sin(beta/2) and cos(beta/2)
    sinb2 = sqrt(quaternion(2)**2 + quaternion(3)**2)
    cosb2 = sqrt(quaternion(1)**2 + quaternion(4)**2)

    ! calculate beta
    beta = 2*acos(cosb2)

    if (cosb2.ne.0) then
       ! determine sin and cos (alpha + gamma)/2
       sinapg2 = quaternion(4)/cosb2
       cosapg2 = quaternion(1)/cosb2

       ! determine alpha + gamma
       apg = 2*atan2(sinapg2,cosapg2)
    else
       apg = 0
    endif

    if (sinb2.ne.0) then
       ! determine sin and cos (alpha + gamma)/2
       sinamg2 = -quaternion(2)/sinb2
       cosamg2 = quaternion(3)/sinb2

       ! determine alpha - gamma
       amg = 2*atan2(sinamg2,cosamg2)
    else
       amg = 0
    endif

    ! now calculate alpha and gamma
    alpha = 0.5_rk*(apg+amg)
    gamma = 0.5_rk*(apg-amg)

    ! identify euler angles
    euler(1) = alpha
    euler(2) = beta
    euler(3) = gamma


  end subroutine efm_trafo_quat_euler
  !*************************************************************


  !*************************************************************
  subroutine efm_ylmrot(l,nms,ldim,alpha,beta,gamma,dmat,invers,dbeta)
    !  Purpose: generates rotation matrix of real spherical
    !           harmonics
    !
    !        CONSTRUCTS A ROTATION MATRIX 'DMAT' OVER THE BASIS OF
    !      CARTESIAN SPHERICAL HARMONICS OF ANGMOM 'L'. THE ROTATION
    !      IS DESCRIBED BY EULER ANGLES 'ALPHA,BETA,GAMMA'. INVERSION
    !      IS INCLUDED IF 'INV'=.TRUE.
    !
    !        THE ROWS & COLS ARE INDEXED BY M, M = L --> -L
    !
    !------------ Modules used ----------------------------------
    use type_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK) :: ldim, l, nms
    real(RK) :: invers,alpha,beta,gamma, &
         DMAT(LDIM,LDIM),DBETA(LDIM,LDIM)
    !** End of interface *****************************************
    integer(IK) :: lp1,m1,m2,mc,mr,lpmr,lmmr,lpmc,lmmc, &
         kmin,mdif,isinx,icosx,k,m1neg,kmax,m2neg,isgn1,isgn2
    real(RK) :: SQR2,BO2,COSBO2,SINBO2,SQRTF,COSF,SINF, &
         ANGF,FACTRL,TERM,SGN,COSAR,SINAR,DBPLUS,DBMINS,COSGC,SINGC
    logical :: INV,COSZO,SINZO,ONEZO
    real(RK), parameter :: ZERO=0.0_rk, &
         ONE=1.0_rk, TWO=2.0_rk, THR=0.001_rk
    real(RK) :: fac(nms)
    ! factorials
    integer(IK) :: max_fac
    ! maximum factorial

    ! prepare factorials
    max_fac = nms
!   allocate(fac(max_fac))
    call facset

    if (NINT(invers).eq.1) then
       inv = .false.
    else
       inv = .true.
    endif

    SQR2= SQRT(TWO)
    BO2=BETA/TWO
    COSBO2= COS(BO2)
    SINBO2= SIN(BO2)
    COSZO=.FALSE.
    SINZO=.FALSE.
    IF ( ABS(COSBO2).LT.THR) COSZO=.TRUE.
    IF ( ABS(SINBO2).LT.THR) SINZO=.TRUE.
    ONEZO=COSZO.OR.SINZO
    LP1=L+1
    !
    !-----------------------------------------------------------------------
    !  GENERATE SMALL DLM`M(BETA):  (ROSE EQ. 4.13, II.17)
    !-----------------------------------------------------------------------
    !----- LOOP OVER LOWER TRIANGLE ----------------------------------------
    !
    DO M1=1,NMS
       MR=LP1-M1
       LPMR=L+MR
       LMMR=L-MR
       DO M2=1,M1
          MC=LP1-M2
          MDIF=MC-MR            ! <- new location
          LPMC=L+MC
          LMMC=L-MC
          SQRTF= SQRT(FAC(LPMR+1)*FAC(LMMR+1))* SQRT(FAC(LPMC+1)*&
          &FAC(LMMC+1))
          DBETA(M1,M2)=ZERO
          KMAX=MIN0(LMMR,LPMC)
          KMIN=MAX0(MC-MR,0)
          IF (KMIN.GT.KMAX) GO TO 50
 !        MDIF=MC-MR            ! <- old location
          IF (ONEZO) GO TO 20
          !
          !----- NEITHER OF COS(BETA/TWO), SIN(BETA/TWO) = ZERO : ---------------
          !
          ISINX=KMIN+KMIN-MDIF
          ICOSX=L+L-ISINX
          COSF=ONE
          SINF=ONE
          IF (ICOSX.NE.0) COSF=  COSBO2 **ICOSX
          IF (ISINX.NE.0) SINF=(-SINBO2)**ISINX
          ANGF=SINF*COSF
          IF (MOD(KMIN,2).EQ.1) ANGF=-ANGF
          DO K=KMIN,KMAX
             FACTRL=FAC(K+1)*FAC(LMMR-K+1)*FAC(LPMC-K+1)*FAC(K-MDIF+1)
             TERM=ANGF/FACTRL
             DBETA(M1,M2)=DBETA(M1,M2)+TERM
             ANGF=-ANGF*(SINBO2/COSBO2)**2
          ENDDO
          GO TO 50
20        CONTINUE
          !
          !----- ONE OF THE ANGULAR FACTORS EQUALS ZERO: ------------------------
          !
          IF (MOD(MDIF,2).EQ.1) GO TO 50
          IF (SINZO) GO TO 30
          K=MDIF/2+L
          ISINX=K+K-MDIF
          ANGF=ONE
          IF (ISINX.NE.0) ANGF=(-SINBO2)**(ISINX)
          GO TO 40
30        K=MDIF/2
          ICOSX=L+L+MDIF-K-K
          ANGF=ONE
          IF (ICOSX.NE.0) ANGF=  COSBO2 **(ICOSX)
40        IF (K.LT.KMIN.OR.K.GT.KMAX) GO TO 50
          IF (MOD(K,2).EQ.1) ANGF=-ANGF
          FACTRL=FAC(K+1)*FAC(LMMR-K+1)*FAC(LPMC-K+1)*FAC(K-MDIF+1)
          DBETA(M1,M2)=ANGF/FACTRL
          !
          !----- ADD SQUARE ROOT FACTOR, FILL IN
          !        THIS ELEMENT OF MATRIX AND THE TRANSPOSE EL`T -----------------
          !
50        DBETA(M1,M2)=SQRTF*DBETA(M1,M2)
          SGN=ONE
          IF (MOD(MDIF,2).EQ.1) SGN=-ONE
          IF (M1.NE.M2) DBETA(M2,M1)=SGN*DBETA(M1,M2)
       ENDDO
    ENDDO
    IF (NMS.EQ.1) GO TO 95
    !
    !-----------------------------------------------------------------------
    !  GENERATE CAPITAL DLM`M(ALPHA,BETA,GAMMA) FOR CARTESIAN HARMONICS:
    !    THE FULL MATRIX IS GENERATED FROM POS M, M` VALUES
    !-----------------------------------------------------------------------
    !
    DO M1=1,L
       M1NEG=NMS-M1+1
       MR=LP1-M1
       COSAR= COS(ALPHA* DBLE(MR))
       SINAR= SIN(ALPHA* DBLE(MR))
       DO M2=1,L
          M2NEG=NMS-M2+1
          MC=LP1-M2
          SGN=ONE
          IF (MOD(MC,2).EQ.1) SGN=-ONE
          DBPLUS=DBETA(M1,M2)+SGN*DBETA(M1,M2NEG)
          DBMINS=DBETA(M1,M2)-SGN*DBETA(M1,M2NEG)
          COSGC= COS(GAMMA* DBLE(MC))
          SINGC= SIN(GAMMA* DBLE(MC))
          DMAT(M1,M2)=       COSAR*COSGC*DBPLUS-SINAR*SINGC*DBMINS
          DMAT(M1,M2NEG)=   -COSAR*SINGC*DBPLUS-SINAR*COSGC*DBMINS
          DMAT(M1NEG,M2)=    SINAR*COSGC*DBPLUS+COSAR*SINGC*DBMINS
          DMAT(M1NEG,M2NEG)=-SINAR*SINGC*DBPLUS+COSAR*COSGC*DBMINS
       ENDDO
       DMAT(M1,LP1)=    SQR2*COSAR*                 DBETA(M1,LP1)
       DMAT(M1NEG,LP1)= SQR2*SINAR*                 DBETA(M1,LP1)
       DMAT(LP1,M1)=    SQR2* COS(GAMMA* DBLE(MR))*DBETA(LP1,M1)
       DMAT(LP1,M1NEG)=-SQR2* SIN(GAMMA* DBLE(MR))*DBETA(LP1,M1)
    ENDDO
95  DMAT(LP1,LP1)=DBETA(LP1,LP1)
    !
    !----- CHANGE THE PHASE OF THE REAL SPHERICAL HARMONICS BY (-1)**M
    !
    !
    !  a more straight forward (and correct) way:
    !
    DO M1=1,NMS
       ISGN1 = 2*MOD(L+M1,2) - 1
       DO M2=1,NMS
          ISGN2 = 2*MOD(L+M2,2) - 1
          DMAT(M1,M2) = ISGN1*ISGN2*DMAT(M1,M2)
       ENDDO
    ENDDO
    !

!   deallocate(fac)

    IF (.NOT.INV.OR.MOD(L,2).EQ.0) RETURN
    !
    !----- IF INVERSION OP`N ON AN ODD-L BASIS, CHANGE SIGN OF MATRIX ------
    !
    !  seems to be a more transparent way of doing the same thing:
    !
    DO M1=1,NMS
       DO M2=1,NMS
          DMAT(M1,M2) = -DMAT(M1,M2)
       ENDDO
    ENDDO
    !

  contains

    subroutine facset
      real(RK), parameter   :: one = 1.0_rk
      integer(IK)           :: i
      real(RK)              :: factmp

      fac(1) = one
      factmp = one
      do i = 2,max_fac
         factmp = dble(i-1)*factmp
         fac(i) = factmp
      end do
    end subroutine facset

  end subroutine efm_ylmrot
  !*************************************************************

#ifdef FPP_COMPILE_UNUSED_CODE
  !*************************************************************
  subroutine efm_real_ylm(quaternion,j,irrep_matrix)
    !  Purpose: determines irreducible representation
    !           matrices for a symmetry operation given
    !           in quaternion form for integral j in
    !           a basis of real spherical harmonics
    !------------ Modules used ----------------------------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in) :: j
    ! j
    real(RK),intent(in)  :: quaternion(5)
    ! quaternion parameters
    real(RK),intent(out) :: irrep_matrix(:,:)
    ! irreducible representation matrix of real spherical
    ! harmonics
    !** End of interface *****************************************
    !------------ Declaration of local constants  ----------------
    complex(CK),parameter :: i = (0.0_rk,1.0_rk)
    ! imaginary one
    real(RK),parameter     :: zero = 0.0_rk,one=1.0_rk
    !------------ Declaration of local variables ---------------------
    integer(IK)          :: alloc_stat
    ! help variables
    integer(IK)          :: dim_trafo,m
    ! 2j+1
    real(RK)             :: sqrt_2_inv,minus_exp_m
    ! inverse of square root of two,(-1)**m
    complex(CK),allocatable :: U_matrix(:,:)
    ! transformation matrix from complex to real spherical
    ! harmonics
    complex(CK),allocatable :: irrep_matrix_c(:,:)
    ! irrep matrix of complex spherical harmonics
    !------------ Executable code ------------------------------------
    ! determine size of transformation
    dim_trafo = 2*j+1
    ! determine inverse square root of two
    sqrt_2_inv = 1/(sqrt(2.0_rk))
    ! allocate U_matrix and irrep_matrix_c
    allocate(U_matrix(dim_trafo,dim_trafo),irrep_matrix_c(dim_trafo,dim_trafo),&
         STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("efm_real_ylm: allocation of fac failed")

    !
    ! determine transformation matrix U_matrix
    !
    U_matrix = zero
    U_matrix(j+1,j+1) = (one,zero)
    minus_exp_m = one
    do m=1,j
       minus_exp_m = (-one)*minus_exp_m
       U_matrix(j+1-m,j+1-m) = sqrt_2_inv
       U_matrix(j+1-m,j+1+m) = i*sqrt_2_inv
       U_matrix(j+1+m,j+1-m) = minus_exp_m*sqrt_2_inv
       U_matrix(j+1+m,j+1+m) = -i*minus_exp_m*sqrt_2_inv
    end do

    !
    ! now determine irrep matrix of complex spherical harmonics
    !
    call efm_so3_irrep(quaternion,2*j,irrep_matrix_c)

    !
    ! determine irrep matrix of real spherical harmonics
    !
    irrep_matrix_c = matmul(conjg(transpose(U_matrix)),irrep_matrix_c)
    irrep_matrix_c = matmul(irrep_matrix_c,U_matrix)
    ! obtain irrep matrix as real part of complex irrep matrix
    ! (imaginary part should vanish)
    irrep_matrix = real(irrep_matrix_c,RK)
    if (NINT(quaternion(5)).ne.1) then
       irrep_matrix = irrep_matrix*(-one)**j
    end if


    ! deallocate U_matrix and irrep_matrix_c
    deallocate(U_matrix,irrep_matrix_c,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("efm_real_ylm: deallocation of fac failed")

  end subroutine efm_real_ylm
  !*************************************************************
#endif

  !*************************************************************
  subroutine efm_so3_irrep(quaternion,two_j,irrep_matrix)
    !  Purpose: determines irreducible representation
    !           matrices in SO3 for a symmetry operation given
    !           in quaternion form
    !------------ Modules used ----------------------------------
    use type_module
    use iounitadmin_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in) :: two_j
    ! 2 * j
    real(RK),intent(in)  :: quaternion(5)
    ! quaternion parameters
    complex(CK),intent(out) :: irrep_matrix(:,:)
    ! irreducible representation matrix
    !** End of interface *****************************************
    !------------ Declaration of local constants  ----------------
    complex(CK),parameter :: i = (0.0_rk,1.0_rk)
    ! imaginary one
    real(RK),parameter     :: zero = 0.0_rk,small = 1.0e-8_rk
    !------------ Declaration of local variables ---------------------
    integer(IK)          :: alloc_stat,max_fac
    ! otherwise
    integer(IK)          :: k
    ! counter
    integer(IK)          :: k_start,k_end
    ! range of k
    integer(IK)          :: m1,m2
    ! magnetic quantum numbers
    integer(IK)          :: two_j_plus_1,j_plus_m,j_minus_mp,&
         m_minus_mp,j_plus_mp,j_minus_m
    ! 2j + 1, j+m, j-m`
    complex(CK)         :: a,b,a_c,b_c
    ! powers of this variables will be taken
    real(RK),allocatable :: fac(:)
    ! array of factorials

    !------------ Executable code ------------------------------------
    two_j_plus_1 = two_j + 1
    ! allocate and set factorials
    allocate(fac(two_j_plus_1),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("efm_so3_irrep: allocation of fac failed")
    max_fac = two_j_plus_1
    call facset
    ! determine a,b
    a   =  quaternion(1) - i*quaternion(4)
    a_c =  quaternion(1) + i*quaternion(4)
    b   = -quaternion(3) - i*quaternion(2)
    b_c = -quaternion(3) + i*quaternion(2)

    !
    ! now determine irrep matrix
    !
    ! initialize irrep_matrix
    irrep_matrix = zero
    ! distinguish cases
    if (abs(b).lt.small) then
       !write(output_unit,*) "b equals zero"
       do m1 = 1,two_j_plus_1
          j_plus_m = m1-1               ! j + m
          j_minus_m = two_j_plus_1-m1   ! j - m
          irrep_matrix(m1,m1) = a**j_plus_m*a_c**j_minus_m
       end do
    elseif (abs(a).lt.small) then
       !write(output_unit,*) "a equals zero"
       do m1 = 1,two_j_plus_1
          m2 = two_j_plus_1 + 1 - m1
          j_plus_m = m1-1               ! j + m
          j_minus_m = two_j_plus_1-m1   ! j - m
          irrep_matrix(m2,m1) = b**j_minus_m*(-b_c)**j_plus_m
       end do
    else
       do m1 = 1,two_j_plus_1
          do m2 =1,two_j_plus_1
             j_plus_m = m1-1               ! j + m
             j_minus_m = two_j_plus_1-m1   ! j - m

             j_plus_mp = m2 - 1            ! j + m`
             j_minus_mp = two_j_plus_1-m2  ! j - m`
             m_minus_mp = m1-m2            ! m - m`
             ! determine range for k
             k_start = max(m_minus_mp,0)
             k_end   = min(j_plus_m,j_minus_mp)

             do k=k_start,k_end
                irrep_matrix(m2,m1) = irrep_matrix(m2,m1) +&
                     a**(j_plus_m-k)*(-b_c)**k*&
                     b**(k-m_minus_mp)*a_c**(j_minus_mp-k)/&
                     (fac(j_plus_m-k+1)*fac(k+1)*fac(k-m_minus_mp+1)&
                     *fac(j_minus_mp-k+1))
             end do
             irrep_matrix(m2,m1) = sqrt(fac(j_plus_m+1)*fac(j_minus_m+1)*&
                  fac(j_plus_mp+1)*fac(j_minus_mp+1))*irrep_matrix(m2,m1)
          end do
       end do
    end if

    ! deallocate factorials
    deallocate(fac,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("efm_so3_irrep: deallocation of fac failed")

  contains

    subroutine facset
      real(RK), parameter   :: one = 1.0_RK
      integer(IK)           :: i
      real(RK)              :: factmp

      fac(1) = one
      factmp = one
      do i = 2,max_fac
         factmp = real(i-1,kind=RK)*factmp
         fac(i) = factmp
      end do
    end subroutine facset


  end subroutine efm_so3_irrep
  !*************************************************************




#if FPP_COMPILE_UNUSED_CODE
  !*************************************************************
  subroutine efm_test_ch
    !  Purpose: little test of the complex eigensolver
    !** End of interface *****************************************
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    !------------ Declaration of local variablrs  -----------------

    integer(IK)              :: i,j,matz,ierr
    ! counters
    integer(IK)              :: dim_trafo
    ! dimension of transformation matrix
    real(RK)                 :: factor
    ! work variable
    real(RK),allocatable     :: exact_eigenvalue(:,:)
    ! matrix which contains the eigenvalues
    complex(CK),allocatable :: exact_eigenvectors(:,:)
    ! complex eigenvectors
    complex(CK),allocatable :: matrix(:,:)
    ! matrix of eigenvalue problem
    real(RK),allocatable     :: calc_eigenvalue(:)
    ! matrix which contains the calculated eigenvalues
    real(RK),allocatable     :: calc_eigenvectors_r(:,:),calc_eigenvectors_i(:,:)
    ! real and imaginary parts of the calculated eigenvectors
    complex(CK),allocatable :: calc_eigenvectors(:,:)
    ! complex calculated eigenvectors
    real(RK),allocatable     :: work1(:),work2(:),work3(:,:)
    ! work arrays

    !------------ Executable code ------------------------------------

    matz = 1
    dim_trafo = 4

    ! allocate arrays
    allocate(exact_eigenvalue(dim_trafo,dim_trafo))
    allocate(exact_eigenvectors(dim_trafo,dim_trafo))
    allocate(calc_eigenvalue(dim_trafo))
    allocate(calc_eigenvectors_r(dim_trafo,dim_trafo),calc_eigenvectors_i(dim_trafo,dim_trafo))
    allocate(calc_eigenvectors(dim_trafo,dim_trafo))
    allocate(matrix(dim_trafo,dim_trafo))
    allocate(work1(dim_trafo),work2(dim_trafo),work3(dim_trafo,dim_trafo))

    ! set eigenvectors
    exact_eigenvectors = 0.0_rk
    exact_eigenvectors(1,1) = (1.0_rk,0.0_rk)
    exact_eigenvectors(2,1) = (0.0_rk,1.0_rk)
    exact_eigenvectors(1,2) = (1.0_rk,0.0_rk)
    exact_eigenvectors(2,2) = (0.0_rk,-1.0_rk)
    exact_eigenvectors(3,3) = (1.0_rk,0.0_rk)
    exact_eigenvectors(4,3) = (0.0_rk,1.0_rk)
    exact_eigenvectors(3,4) = (1.0_rk,0.0_rk)
    exact_eigenvectors(4,4) = (0.0_rk,-1.0_rk)

    exact_eigenvectors = exact_eigenvectors/sqrt(2.0_rk)

    ! set eigenvalues
    exact_eigenvalue = 0.0_rk
    exact_eigenvalue(1,1) = 1.0_rk
    exact_eigenvalue(2,2) = 2.0_rk
    exact_eigenvalue(3,3) = 3.0_rk
    exact_eigenvalue(4,4) = 4.0_rk

    ! construct matrix
    matrix = matmul(exact_eigenvalue,conjg(transpose(exact_eigenvectors)))
    matrix = matmul(exact_eigenvectors,matrix)

    ! solve eigenvalue problem
    call ch(dim_trafo,dim_trafo,&                                        ! dimensions of transformation
         real(matrix),aimag(matrix),&                                    ! matrix of eigenvalue problem
         calc_eigenvalue,matz,calc_eigenvectors_r,calc_eigenvectors_i,&  ! eigenvalues and eigenvectors
         work1,work2,work3,ierr)                                         ! work arrays
    if (ierr.ne.0) then
       write(*,*) "efm_test_ch: stopped due to error in eigensolver IERR =",ierr
       call error_handler("efm_test_ch: error in eigensolver")
    end if

    ! construct complex eigenvector
    calc_eigenvectors = cmplx(calc_eigenvectors_r,calc_eigenvectors_i,RK)

    write(output_unit,*) "Test of complex eigen solver "
    write(output_unit,*) "calculated Eigenvectors"
    do i = 1,dim_trafo
       write(output_unit,602)&
            &         (calc_eigenvectors(i,j), j = 1,dim_trafo)
    end do

    do i=1,dim_trafo
       factor = sqrt(dot_product(calc_eigenvectors(:,i),calc_eigenvectors(:,i)))
       write(output_unit,*)" Eigenvector ",i, " normalized to ",factor
       calc_eigenvectors(:,i) = calc_eigenvectors(:,i)/factor
    end do

    write(output_unit,*) "calculated Eigenvectors normalized"
    do i = 1,dim_trafo
       write(output_unit,602)&
            &         (calc_eigenvectors(i,j), j = 1,dim_trafo)
    end do

    write(output_unit,*) "exact Eigenvectors"
    do i = 1,dim_trafo
       write(output_unit,602)&
            &         (exact_eigenvectors(i,j), j = 1,dim_trafo)
    end do
    write(output_unit,*) "Eigenvalues "
    do i = 1,dim_trafo
       write(output_unit,*) calc_eigenvalue(i)
    end do

    ! deallocate arrays
    deallocate(exact_eigenvalue)
    deallocate(exact_eigenvectors)
    deallocate(matrix)
    deallocate(calc_eigenvalue)
    deallocate(calc_eigenvectors_r,calc_eigenvectors_i)
    deallocate(calc_eigenvectors)
    deallocate(work1,work2,work3)

602 format(1x,f17.12,1x,f17.12,2x,4(f17.12,f17.12:)/4(10x,4(f17.12,f17.12:)/))

  end subroutine efm_test_ch
  !*************************************************************
#endif


  !*************************************************************
  subroutine efm_solve_nonsq(col_dim,row_dim,matrix,in_out,perm)
    !  Purpose: solves a nonsquare linear equation
    !           On input is a col_dim x row_dim matrix.
    !           matrix contains row_dim column vectors.
    !           It must hold: row_dim >= col_dim.
    !           The col_dim x col_dim matrix in_out
    !           contains col_dim columnvectors.
    !           The equation matrix * x = b is solved for all column
    !           vectors. The solution is a linear combination of
    !           only col_dim of the matrix column vectors
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    !------------ Declaration of formal parameters  -----------------
    integer(IK), intent(in)           :: col_dim,row_dim
    ! dimensions col_dim x row_dim of matrix
    real(RK),    intent(inout)        :: matrix(col_dim,row_dim)
    ! the matrix, it contains row_dim columnvector of length col_dim
    real(RK),    intent(inout)        :: in_out(:,:)
    ! contains the right sides of the equation on input and
    ! the solution vectors on output
    integer(IK), intent(out)          :: perm(:)
    ! list permutations of the row_dim column vectors of matrix
    !** End of interface *****************************************
    !------------ Declaration of local parameters ---------------
    real(RK),parameter               :: small = 0.000001
    ! very small value

    !------------ Declaration of local variables  ---------------
    integer(IK),allocatable          :: permutation(:)
    ! list permutations of the row_dim column vectors of matrix
    real(RK),allocatable             :: intermediate(:)
    ! help variables used for determation of the pivot element
    integer(IK)                      :: i,j,k,max_j,some_int
    ! counters
    real(RK)                         :: max_value,new_value,some
    ! help variables used for determation of the pivot element

    !------------ Executable code ------------------------------------

    ! allocate help arrays
    allocate(permutation(row_dim))
    allocate(intermediate(col_dim))

    ! check if col_dim > row_dim
    if (col_dim.gt.row_dim) call error_handler("efm_solve_nonsq: col_dim > row_dim")

    ! initialize permutation
    do i=1,row_dim
       permutation(i) = i
    end do

    do i=1,col_dim
       ! determine element of maximal size
       ! in the row i
       max_j = i
       max_value = abs(matrix(i,i))
       do j=i+1,row_dim
          new_value = abs(matrix(i,j))
          if (new_value.gt.max_value) then
             max_j = j
             max_value = new_value
          end if
       end do
       ! if no nonzero element was found then stop execution
       if (max_value.lt.small) call error_handler("efm_solve_nonsq: rank of matrix < col_dim")

       ! interchange rows i and max_j
       if (max_j.ne.i) then
          intermediate(:) = matrix(:,max_j)
          matrix(:,max_j) = matrix(:,i)
          matrix(:,i) = intermediate
          some_int = permutation(max_j)
          permutation(max_j) = permutation(i)
          permutation(i) = some_int
       end if

       ! do Gauss-elimination
       do j=i+1,col_dim
          some = -matrix(j,i)/matrix(i,i)
          matrix(j,i) = 0.0_rk
          do k=i+1,row_dim
             matrix(j,k) = matrix(j,k) + some*matrix(i,k)
          end do
          do k=1,col_dim
             in_out(j,k) = in_out(j,k) + some*in_out(i,k)
          end do
       end do
    end do

    ! do backsubstitution
    do i=col_dim,1,-1
       some = 1/matrix(i,i)
       in_out(i,:) = in_out(i,:)*some
       do j=i-1,1,-1
          some = -matrix(j,i)
          in_out(j,:) = in_out(j,:) + in_out(i,:)*some
       end do
    end do

    perm = permutation(1:col_dim)

    ! allocate help arrays
    deallocate(permutation)
    deallocate(intermediate)

  end subroutine efm_solve_nonsq
  !*************************************************************

  !*************************************************************
  subroutine efm_solve_nonsq_complex(col_dim,row_dim,matrix,in_out,perm)
    !  Purpose: solves a nonsquare linear equation
    !           On input is a col_dim x row_dim matrix.
    !           matrix contains row_dim column vectors.
    !           It must hold: row_dim >= col_dim.
    !           The col_dim x col_dim matrix in_out
    !           contains col_dim columnvectors.
    !           The equation matrix * x = b is solved for all column
    !           vectors. The solution is a linear combination of
    !           only col_dim of the matrix column vectors
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    !------------ Declaration of formal parameters  -----------------
    integer(IK), intent(in)           :: col_dim,row_dim
    ! dimensions col_dim x row_dim of matrix
    complex(CK),    intent(inout)        :: matrix(col_dim,row_dim)
    ! the matrix, it contains row_dim columnvector of length col_dim
    complex(CK),    intent(inout)        :: in_out(:,:)
    ! contains the right sides of the equation on input and
    ! the solution vectors on output
    integer(IK), intent(out)          :: perm(:)
    ! list permutations of the row_dim column vectors of matrix
    !** End of interface *****************************************
    !------------ Declaration of local parameters ---------------
    real(RK),parameter               :: small = 0.000001_rk
    ! very small value

    !------------ Declaration of local variables  ---------------
    integer(IK),allocatable          :: permutation(:)
    ! list permutations of the row_dim column vectors of matrix
    complex(CK),allocatable         :: intermediate(:)
    ! help variables used for determation of the pivot element
    integer(IK)                      :: i,j,k,max_j,some_int
    ! counters
    real(RK)                         :: max_value,new_value
    complex(CK)                     :: some
    ! help variables used for determation of the pivot element

    !------------ Executable code ------------------------------------

    ! allocate help arrays
    allocate(permutation(row_dim))
    allocate(intermediate(col_dim))

    ! check if col_dim > row_dim
    if (col_dim.gt.row_dim) call error_handler("efm_solve_nonsq: col_dim > row_dim")

    ! initialize permutation
    do i=1,row_dim
       permutation(i) = i
    end do

    do i=1,col_dim
       ! determine element of maximal size
       ! in the row i
       max_j = i
       max_value = abs(matrix(i,i))
       do j=i+1,row_dim
          new_value = abs(matrix(i,j))
          if (new_value.gt.max_value) then
             max_j = j
             max_value = new_value
          end if
       end do
       ! if no nonzero element was found then stop execution
       if (max_value.lt.small) call error_handler("efm_solve_nonsq_complex: rank of matrix < col_dim")

       ! interchange rows i and max_j
       if (max_j.ne.i) then
          intermediate(:) = matrix(:,max_j)
          matrix(:,max_j) = matrix(:,i)
          matrix(:,i) = intermediate
          some_int = permutation(max_j)
          permutation(max_j) = permutation(i)
          permutation(i) = some_int
       end if

       ! do Gauss-elimination
       do j=i+1,col_dim
          some = -matrix(j,i)/matrix(i,i)
          matrix(j,i) = 0.0_rk
          do k=i+1,row_dim
             matrix(j,k) = matrix(j,k) + some*matrix(i,k)
          end do
          do k=1,col_dim
             in_out(j,k) = in_out(j,k) + some*in_out(i,k)
          end do
       end do
    end do

    ! do backsubstitution
    do i=col_dim,1,-1
       some = 1/matrix(i,i)
       in_out(i,:) = in_out(i,:)*some
       do j=i-1,1,-1
          some = -matrix(j,i)
          in_out(j,:) = in_out(j,:) + in_out(i,:)*some
       end do
    end do

    perm = permutation(1:col_dim)

    ! allocate help arrays
    deallocate(permutation)
    deallocate(intermediate)

  end subroutine efm_solve_nonsq_complex
  !*************************************************************

#if FPP_COMPILE_UNUSED_CODE
  !*************************************************************
  subroutine efm_cscoIII_solve_r(trafo_matrix,trafo_matrix_intr)
    !  Purpose: solves eigenequation for the CSCO-III of the
    !           regular representation,
    !           thus determines irreducible representations of
    !           the group.
    !           The CSCO-III is solved as a real equation
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    use group_module
    implicit none
    !------------ Declaration of formal parameters  -----------------
    real(RK),    intent(in)           :: trafo_matrix(:,:,:)
    ! transformation matrices of the group operation in the regular representation
    real(RK),    intent(in)           :: trafo_matrix_intr(:,:,:)
    ! transformation matrices of the intrinsic group operation in the regular
    ! representation
    !** End of interface *****************************************

    !------------ Declaration of local parameters ---------------
    integer(IK),parameter                :: matz = 1
    ! parameter for eigensolver
    real(RK),parameter                   :: small = 0.00001_rk
    ! very small value
    real(RK),parameter                   :: sq2 = 0.353553391_rk
    ! prefactor
    !------------ Declaration of local variables  ---------------
    real(RK),allocatable             :: efm_matrix(:,:),efmII_matrix(:,:)
    ! contains matrix of eigenvalue problem for CSCO and CSCO-II
    real(RK),allocatable             :: efmIII_matrix(:,:)
    ! contains matrix of eigenvalue problem for CSCO-III
    real(RK),allocatable             :: csco_sub_intr(:,:)
    ! csco-matrix of the intrinsic group (labels multiple irreps)
    real(RK),allocatable             :: csco_sub(:,:)
    ! csco-matrix of the group (labels irreps)
    real(RK),allocatable             :: projector(:,:)
    ! work array of dimension group_num_el x group_num_el
    real(RK),allocatable             :: intermediate(:,:)
    ! work array of dimension group_num_el x group_num_el
    real(RK),allocatable             :: intermediate_irr(:,:)
    ! work array of dimension group_num_el x irrep_dim
    real(RK),allocatable             :: intermediate_irr2(:,:)
    ! work array of dimension irrep_dim x irrep_dim
    real(RK),allocatable                 :: intermediate_irr3(:,:)
    ! work array of dimension irrep_dim x group_num_el
    real(RK),allocatable                 :: eigen_value(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable                 :: eigen_valueII(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable                 :: eigen_valueIII(:)
    ! eigenvalues of CSCO-III
    real(RK),allocatable                 :: eigen_value_intr(:)
    ! eigenvalues of the CSCO of the intrinsic subgroups
    real(RK),allocatable                 :: eigen_value_sub(:)
    ! eigenvalues of the CSCO of the subgroups
    real(RK),allocatable                :: eigen_vector(:,:)
    ! eigenvectors of CSCO-III and CSCO-II and CSCO
    ! eigen_vector(group_num_el,group_num_el)
    ! eigen_vector(basis_function,eigen_vector_number)
    real(RK),allocatable                 :: irrep_basis(:,:)
    ! basis for one irrep: irrep_basis(group_num_el,irrep_dim)
    real(RK),allocatable                :: eigen_space(:,:,:)
    ! eigen_space of a certain irrep
    real(RK),allocatable                 :: irrep_eigen_value(:)
    ! eigen_values of irrep
    real(RK),allocatable                 :: work1(:),work2(:)
    ! real work arrays
    integer(IK) :: i, j, k, m, n, ierr
    ! counters
    integer(IK)                          :: irrep_dim
    ! help variable for dimension of irrep
    real(RK)                             :: factor
    ! factor for subgroups
    real(RK)                             :: irrep_eigenvalue, irrep_mult_label
    ! irrep and irrep multiplicity labels
    real(RK)                             :: irrep_character
    ! calculated character of some operation of the irrep

    logical :: found

    !------------ Executable code ------------------------------------

    ! allocate local arrays
    allocate(efm_matrix(group_num_el,group_num_el))
    allocate(efmII_matrix(group_num_el,group_num_el))
    allocate(efmIII_matrix(group_num_el,group_num_el))
    allocate(csco_sub_intr(group_num_el,group_num_el),csco_sub(group_num_el,group_num_el))
    allocate(projector(group_num_el,group_num_el))
    allocate(intermediate(group_num_el,group_num_el))
    allocate(eigen_value(group_num_el))
    allocate(eigen_valueII(group_num_el))
    allocate(eigen_valueIII(group_num_el))
    allocate(eigen_value_intr(group_num_el),eigen_value_sub(group_num_el))
    allocate(eigen_vector(group_num_el,group_num_el))
    allocate(irrep_eigen_value(group_num_el))

    ! allocate work arrays
    allocate(work1(group_num_el),work2(group_num_el))

    ! initialize efm-matrix
    efm_matrix = 0

    ! first build efm for CSCO
    do i = 1,3
       write(output_unit,*) "ncsco (",i,") ",ncsco(i)
       write(output_unit,*) "fcsco (",i,") ",fcsco(i)
       if (ncsco(i).ne.0) then
          if (klass(ncsco(i))%ambivalent) then
             do j=1,klass(ncsco(i))%number
                write(output_unit,*) symm_el(klass(ncsco(i))%conjug_el(j))%name
                efm_matrix = efm_matrix +&
                     trafo_matrix(:,:,klass(ncsco(i))%conjug_el(j))*fcsco(i)
             end do
          else
             do j=1,klass(ncsco(i))%number
                write(output_unit,*) symm_el(klass(ncsco(i))%conjug_el(j))%name
                intermediate = trafo_matrix(:,:,klass(ncsco(i))%conjug_el(j))
                efm_matrix = efm_matrix +&
                     3*(intermediate + transpose(intermediate))*fcsco(i)*sq2
             end do
          endif
       endif
    end do

    ! now build efm for CSCO-II

    efmII_matrix = efm_matrix

    !
    ! build CSCO of subgroup chain
    !

    csco_sub = 0.0_rk
    factor = 10.0_rk
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups(m)%ncsco(i).ne.0) then
             if (subgroups(m)%klass(subgroups(m)%ncsco(i))%ambivalent) then
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   csco_sub = csco_sub +&
                        trafo_matrix(:,:,subgroups(m)%klass(subgroups(m)%ncsco(i))%&
                        &conjug_el(j))*subgroups(m)%fcsco(i)*factor
                end do
             else
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   intermediate = trafo_matrix(:,:,subgroups(m)%klass(subgroups(m)&
                        &%ncsco(i))%conjug_el(j))
                   csco_sub = csco_sub + subgroups(m)%fcsco(i)*factor*sq2*&
                        &(3*(intermediate + transpose(intermediate)))
                end do
             endif
          endif
       end do
       factor = factor*10.0_rk
    end do

    ! add projectors for subgroup T(if existent)
    projector = 0.0_rk
    if (subgroups(1)%name.eq.'T   ') then
       do i=1,group_num_el
          projector = projector + E_irrep%matrix(1,1,i)*trafo_matrix(:,:,i)
       end do
       projector = projector*200.0_rk/group_num_el
       csco_sub = csco_sub + projector
    end if

    efmII_matrix = efmII_matrix + csco_sub


    ! now build efm for CSCO-III
    efmIII_matrix = efmII_matrix


    !
    ! build CSCO of intrinsic subgroup chain
    !

    ! initialize csco-matrix of intrinsic subgroups
    csco_sub_intr = 0.0_rk

    factor = 10.0_rk
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups(m)%ncsco(i).ne.0) then
             if (subgroups(m)%klass(subgroups(m)%ncsco(i))%ambivalent) then
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   csco_sub_intr = csco_sub_intr +&
                        trafo_matrix_intr(:,:,subgroups(m)%klass(subgroups(m)%ncsco(i))%&
                        &conjug_el(j))*subgroups(m)%fcsco(i)*factor
                end do
             else
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   intermediate = trafo_matrix_intr(:,:,subgroups(m)%klass(subgroups(m)&
                        &%ncsco(i))%conjug_el(j))
                   csco_sub_intr = csco_sub_intr + subgroups(m)%fcsco(i)*factor*sq2*&
                        &(3*(intermediate + transpose(intermediate)))
                end do
             endif
          endif
       end do
       factor = factor*10.0_rk
    end do

    ! add projectors for subgroup T(if existent)
    projector = 0.0_rk
    if (subgroups(1)%name.eq.'T   ') then
       do i=1,group_num_el
          projector = projector + E_irrep%matrix(1,1,i)*trafo_matrix_intr(:,:,i)
       end do
       projector = projector*200.0_rk/group_num_el
       csco_sub_intr = csco_sub_intr + projector
       write(output_unit,*) "Projector"
       do i = 1,group_num_el
          write(output_unit,602)&
               &         (projector(i,j), j = 1,group_num_el)
       end do

    end if

    ! add up final CSCO-III matrix
    efmIII_matrix = efmIII_matrix + csco_sub_intr*100.0_rk


    ! solve eigenvalue problem for CSCO-III
    ierr = 0
    call rs(group_num_el,group_num_el,&               ! dimensions of transformation
         efmIII_matrix,&                              ! matrix of eigenvalue problem
         eigen_valueIII,matz,eigen_vector,&           ! eigenvalues and eigenvectors
         work1,work2,ierr)                            ! work arrays


    !
    ! --- In the following eigenvalues of the different CSCO`s as CSCO-I, CSCO-II, -----
    ! --- and the CSCO-III are determined. ---------------------------------------------
    !

    ! determine eigenvalues of CSCO (i.e. the irrep labels)

    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(efm_matrix,eigen_vector)

    ! determine the multiple
    eigen_value = 0
    do j=1,group_num_el
       do i=1,group_num_el
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value(j) = -intermediate(i,j)/eigen_vector(i,j)
             exit
          endif
       enddo
    enddo

    ! determine eigenvalues of CSCO of the intrinisc subgroups (i.e. the
    ! multiplicity labels)

    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(csco_sub_intr,eigen_vector)

    ! determine the multiple
    eigen_value_intr = 0
    do j=1,group_num_el
       do i=1,group_num_el
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value_intr(j) = -intermediate(i,j)/eigen_vector(i,j)
             exit
          endif
       enddo
    enddo

    ! determine eigenvalues of CSCO of the subgroups (i.e. the irrep
    ! labels)

    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(csco_sub,eigen_vector)

    ! determine the multiple
    eigen_value_sub = 0
    do j=1,group_num_el
       do i=1,group_num_el
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value_sub(j) = -intermediate(i,j)/eigen_vector(i,j)
             exit
          endif
       enddo
    enddo

    ! show the eigen values
    do j=1,group_num_el
       write(output_unit,*) "Nr. ",j,"irrep ",eigen_value(j),&
            "basis: ",eigen_value_sub(j)," multi ",eigen_value_intr(j)
    end do

    !
    ! -------------- determine the irrep matrices -------------------
    !

    ! loop over all irreps
    do i=1,group_num_ir

       ! read eigenvalue and dimension of irrep
       irrep_eigenvalue = irrep_can(i)%eigen_value
       irrep_dim = irrep_can(i)%dimension

       ! for subsequent use allocate help arrays of the suitable size
       ! they will be reallocated for the next irrep
       allocate(intermediate_irr(group_num_el,irrep_dim))
       allocate(intermediate_irr2(irrep_dim,irrep_dim))
       allocate(intermediate_irr3(irrep_dim,group_num_el))
       allocate(eigen_space(group_num_el,irrep_dim,irrep_dim))
       allocate(irrep_basis(group_num_el,irrep_dim))

       ! look up for one eigenvector of the current irrep and memorize
       ! its multiplicy label
       irrep_mult_label = HUGE (irrep_mult_label) ! sentinel
       found = .false.
       do j = 1, group_num_el
          if (abs(eigen_value(j)-irrep_eigenvalue).lt.small) then
             irrep_mult_label = eigen_value_intr(j)
             found = .true.
             exit
          endif
       end do
       ASSERT(found)

       ! determine the labels of the subspace
       n = 1
       do j=1,group_num_el
          if ((abs(eigen_value(j)-irrep_eigenvalue).lt.small).and.&
               &(abs(eigen_value_intr(j)-irrep_mult_label).lt.small)) then
             irrep_eigen_value(n) = eigen_value_sub(j)
             n = n+1
          end if
       end do

       ! fill up array eigen_space
       do n=1,irrep_dim
          do k=1,irrep_dim
             do j=1,group_num_el
                if ((abs(eigen_value(j)-irrep_eigenvalue).lt.small).and.&
                     &(abs(eigen_value_sub(j)-irrep_eigen_value(n)).lt.small).and.&
                     &(abs(eigen_value_intr(j)-irrep_eigen_value(k)).lt.small)) then
                   eigen_space(:,n,k) = eigen_vector(:,j)
                   exit
                end if
             end do
          end do
       end do

       ! fill irrep basis, the basis from which the transformation properties
       ! will be obtained
       irrep_basis(:,:) = eigen_space(:,:,1)

       ! show irrep basis
       write(output_unit,*) "after first normalization"
       do n = 1,group_num_el
          write(output_unit,602)&
               &         (irrep_basis(n,m), m = 1,irrep_dim)
       end do


       ! determine eigenvalues of CSCO-II for the basis of the irrep

       ! multiply real CSCO-II efm-matrix by the eigenvectors
       intermediate_irr = matmul(efmII_matrix,irrep_basis)

       ! determine the multiple
       eigen_valueII = 0
       do j=1,irrep_dim
          do n=1,group_num_el
             if (abs(intermediate_irr(n,j)).gt.small) then
                eigen_valueII(j) = -(intermediate_irr(n,j)/irrep_basis(n,j))
                exit
             endif
          enddo
       enddo

       !
       ! -----  some tests -------------------------------------------------------------
       !

       ! show irrep basis
       write(output_unit,*) "basis of irrep"
       do n = 1,group_num_el
          write(output_unit,602)&
               &         (irrep_basis(n,m), m = 1,irrep_dim)
       end do

       ! check if irrep basis is orthonormalized
       intermediate_irr2 = matmul(transpose(irrep_basis), irrep_basis)
       write(output_unit,*)
       write(output_unit,*) "======================================================="
       write(output_unit,*) "Irrep ",irrep_can(i)%label
       write(output_unit,*) "overlap matrix of basis:"
       do n = 1,irrep_dim
          write(output_unit,602)&
               &         (intermediate_irr2(n,m), m = 1,irrep_dim)
       end do

       ! show all eigenvectors of the irrep
       write(output_unit,*) "eigen spaces of the irrep",irrep_can(i)%label
       do k = 1,irrep_dim
          write(output_unit,*)
          write(output_unit,*) "Eigenspace Nr. ",k
          do n = 1,group_num_el
             write(output_unit,602)&
                  &         (eigen_space(n,m,k), m = 1,irrep_dim)
          end do
       end do

       !
       ! --------- determine the representations of the irreps ---------------------
       !

       ! allocate representation matrix of the irrep
       allocate(irrep_can(i)%irrep_matrix(irrep_dim,irrep_dim,group_num_el))
       ! allocate irrep basis generators
       allocate(irrep_can(i)%irrep_gen_oper(irrep_dim))
       allocate(irrep_can(i)%irrep_gen_coeff(irrep_dim,irrep_dim))

       ! initialize irrep basis generator
       irrep_can(i)%irrep_gen_coeff = 0.0_rk
       do j=1,irrep_dim
          irrep_can(i)%irrep_gen_coeff(j,j) = 1.0_rk
       end do

       ! now determine the irrep matrix from the formula
       ! irrep matrix = C^T * TRAFO * C, where the columns of C are the
       ! basis vectors of the irrep
       do j=1,group_num_el
          intermediate_irr = matmul(trafo_matrix(:,:,j),irrep_basis)
          intermediate_irr2 = matmul(transpose(irrep_basis), intermediate_irr)
          ! show irrep matrix
          write(output_unit,*) "Irrep ",irrep_can(i)%label
          write(output_unit,*) "Representation of ",symm_el(j)%name
          do n = 1,irrep_dim
             write(output_unit,602)&
                  &         (intermediate_irr2(n,m), m = 1,irrep_dim)
          end do
          ! write irrep matrix to its destination
          irrep_can(i)%irrep_matrix(:,:,j) = intermediate_irr2
          ! determine character of irrep matrix
          irrep_character = 0.0_rk
          do n=1,irrep_dim
             irrep_character = irrep_character + intermediate_irr2(n,n)
          end do
          write(output_unit,*) "Character of the irrep calculated ",irrep_character
          write(output_unit,*) "Character of the irrep exact      ",irrep_can(i)%&
               &characters(symm_el(j)%klass)
          ! check unitarity of basis
          intermediate_irr2 = matmul(transpose(intermediate_irr2), intermediate_irr2)
          write(output_unit,*) "test unitarity of representation"
          do n = 1,irrep_dim
             write(output_unit,602)&
                  &         (intermediate_irr2(n,m), m = 1,irrep_dim)
          end do
       end do

       ! determine eigenvector of the basisvector from which the eigenfunctions
       ! are generated
       irrep_can(i)%gen_eigen_value  = eigen_valueII(1)
       irrep_can(i)%build_projector = .false.


       ! determine generators of the irrep basis

       ! collect the first columns of all irrep matrices, which are the pictures
       ! of the first basis vector for the corresponding group elements
       do j=1,group_num_el
          intermediate_irr3(:,j) = irrep_can(i)%irrep_matrix(:,1,j)
       end do
       call efm_solve_nonsq(irrep_dim,group_num_el,intermediate_irr3(:,:),&
            irrep_can(i)%irrep_gen_coeff,irrep_can(i)%irrep_gen_oper)

       ! test irrep basis generator
       do j=1,irrep_dim
          intermediate_irr3(:,j) = irrep_can(i)%irrep_matrix(:,1,&
               &irrep_can(i)%irrep_gen_oper(j))
       end do
       intermediate_irr3(:,1:irrep_dim) = matmul(intermediate_irr3(:,1:irrep_dim),irrep_can(i)%irrep_gen_coeff)
       write(output_unit,*) "Irrep ",irrep_can(i)%label
       write(output_unit,*) "test of gauss jordan"
       do n = 1,irrep_dim
          write(output_unit,602)&
               &         (intermediate_irr3(n,m), m = 1,irrep_dim)
602       format(1x,f9.4,1x,f9.4,2x,4(f9.4,f9.4:)/4(10x,4(f9.4,f9.4:)/))
       end do

       ! deallocate all help arrays for this irrep
       deallocate(intermediate_irr,intermediate_irr2,intermediate_irr3)
       deallocate(irrep_basis,eigen_space)
    end do


    ! deallocate local arrays
    deallocate(efm_matrix,efmII_matrix,efmIII_matrix)
    deallocate(csco_sub_intr,csco_sub,intermediate)
    deallocate(projector)
    deallocate(eigen_value,eigen_valueII,eigen_valueIII,eigen_value_intr)
    deallocate(eigen_value_sub,eigen_vector)
    deallocate(irrep_eigen_value)

    ! deallocate work arrays
    deallocate(work1,work2)

  end subroutine efm_cscoIII_solve_r
  !*************************************************************
#endif


  !*************************************************************
  subroutine efm_cscoIII_solve(trafo_matrix,trafo_matrix_intr)
    !  Purpose: solves eigenequation for the CSCO-III of the
    !           regular representation,
    !           thus determines irreducible representations of
    !           the group.
    !           the representations are made real
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    use group_module
    implicit none
    !------------ Declaration of formal parameters  -----------------
    real(RK),    intent(in)           :: trafo_matrix(:,:,:)
    ! transformation matrices of the group operation in the regular representation
    real(RK),    intent(in)           :: trafo_matrix_intr(:,:,:)
    ! transformation matrices of the intrinsic group operation in the regular
    ! representation
    !** End of interface *****************************************

    !------------ Declaration of local parameters ---------------
    integer(IK),parameter                :: matz = 1
    ! parameter for eigensolver
    real(RK),parameter                   :: small = 0.00001_rk
    ! very small value

    !------------ Declaration of local variables  ---------------
    complex(CK),allocatable             :: efm_matrix(:,:),efmII_matrix(:,:)
    ! contains matrix of eigenvalue problem for CSCO and CSCO-II
    complex(CK),allocatable             :: efmIII_matrix(:,:)
    ! contains matrix of eigenvalue problem for CSCO-III
    complex(CK),allocatable             :: csco_sub_intr(:,:)
    ! csco-matrix of the intrinsic group (labels multiple irreps)
    complex(CK),allocatable             :: csco_sub(:,:)
    ! csco-matrix of the group (labels irreps)
    complex(CK),allocatable             :: intermediate(:,:)
    ! work array of dimension group_num_el x group_num_el
    complex(CK),allocatable             :: intermediate_irr(:,:)
    ! work array of dimension group_num_el x irrep_dim
    complex(CK),allocatable             :: intermediate_irr2(:,:)
    ! work array of dimension irrep_dim x irrep_dim
    real(RK),allocatable                 :: intermediate_irr3(:,:)
    ! work array of dimension irrep_dim x group_num_el
    real(RK),allocatable                 :: eigen_value(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable                 :: eigen_valueII(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable                 :: eigen_valueIII(:)
    ! eigenvalues of CSCO-III
    real(RK),allocatable                 :: eigen_value_intr(:)
    ! eigenvalues of the CSCO of the intrinsic subgroups
    real(RK),allocatable                 :: eigen_value_sub(:)
    ! eigenvalues of the CSCO of the subgroups
    complex(CK),allocatable             :: eigen_vector(:,:)
    ! eigenvectors of CSCO-III and CSCO-II and CSCO
    ! eigen_vector(group_num_el,group_num_el)
    ! eigen_vector(basis_function,eigen_vector_number)
    real(RK),allocatable                 :: eig_vec_r(:,:),eig_vec_i(:,:)
    ! real and imaginary parts of eigen_vector
    complex(CK),allocatable             :: irrep_basis(:,:)
    ! basis for one irrep: irrep_basis(group_num_el,irrep_dim)
    complex(CK),allocatable             :: eigen_space(:,:,:)
    ! eigen_space of a certain irrep
    real(RK),allocatable                 :: irrep_eigen_value(:)
    ! eigen_values of irrep
    real(RK),allocatable                 :: work1(:),work2(:),work3(:,:)
    ! real work arrays
    complex(CK),allocatable             :: workc1(:),workc2(:)
    ! complex work arrays
    integer(IK) :: i, j, k, l, m, n, ierr
    ! counters
    integer(IK)                          :: irrep_dim
    ! help variable for dimension of irrep
    integer(IK)                          :: mult_label,basis_label
    ! help variables which hold the labels of some conjugated element
    real(RK)                             :: factor,sq2
    ! factor for subgroups
    real(RK)                             :: pi
    ! pi, the circle number
    real(RK)                             :: abs_factor
    ! factor for making real of eigenvectors
    complex(CK)                         :: c_factor
    ! complex factor for making real of eigenvectors
    complex(CK)                         :: phase_factor
    ! complex phase factor which adjusts two dimensional irreps to
    ! a special transformation (only for the groups Td,O,Oh and I,Ih)
    real(RK)                             :: irrep_eigenvalue, irrep_mult_label
    ! irrep and irrep multiplicity labels
    real(RK)                             :: eigen_sub,eigen_intr
    ! eigenvalues of conjugated element
    logical                                        :: first_vector_real
    ! is the first basisvector of an irrepbasis real?

    logical :: found


    !------------ Executable code ------------------------------------

    ! allocate local arrays

    allocate(efm_matrix(group_num_el,group_num_el))
    allocate(efmII_matrix(group_num_el,group_num_el))
    allocate(efmIII_matrix(group_num_el,group_num_el))
    allocate(csco_sub_intr(group_num_el,group_num_el),csco_sub(group_num_el,group_num_el))
    allocate(intermediate(group_num_el,group_num_el))
    allocate(eigen_value(group_num_el))
    allocate(eigen_valueII(group_num_el))
    allocate(eigen_valueIII(group_num_el))
    allocate(eigen_value_intr(group_num_el),eigen_value_sub(group_num_el))
    allocate(eigen_vector(group_num_el,group_num_el))
    allocate(eig_vec_r(group_num_el,group_num_el),eig_vec_i(group_num_el,group_num_el))
    allocate(irrep_eigen_value(group_num_el))

    ! allocate work arrays

    allocate(work1(group_num_el),work2(group_num_el))
    allocate(workc1(group_num_el),workc2(group_num_el))
    allocate(work3(2,group_num_el))


    !
    ! determine phasefactor (only needed for polyhedral irreps)
    !

    select case(group_name)
    case('Td  ','O   ','Oh  ')
       pi = acos(-1.0_rk)*2.0_rk/3.0_rk
       phase_factor = cmplx(cos(pi),sin(pi),RK)
    case('I   ','Ih  ')
       phase_factor = (0.6123724356958248_rk,0.7905694150420700_rk)
    case default
       phase_factor = (1.0_rk,0.0_rk)
    end select

    !
    ! determine factors for EFM-Matrix
    !

    select case(group_name)
    case('Td  ','O   ','Oh  ')
       efm_subgroup_factor = 5.0_rk
       efm_group_factor = 1.0_rk
       efm_intrgroup_factor = 2.0_rk
    case('I   ','Ih  ')
       efm_subgroup_factor = 5.0_rk
       efm_group_factor = 3.0_rk
       efm_intrgroup_factor = 6.0_rk
    case default
       efm_subgroup_factor = 10.0_rk
       efm_group_factor = 4.0_rk
       efm_intrgroup_factor = 12.0_rk
    end select

    !
    ! construct EFM-matrix of the CSCO-III
    !

    ! initialize efm-matrix
    efm_matrix = (0.0_rk,0.0_rk)

    ! prefactor
    sq2 = 1/(2.0_rk*sqrt(2.0_rk))

    ! first build efm for CSCO
    do i = 1,3
       if (ncsco(i).ne.0) then
          if (klass(ncsco(i))%ambivalent) then
             do j=1,klass(ncsco(i))%number
                efm_matrix = efm_matrix +&
                     trafo_matrix(:,:,klass(ncsco(i))%conjug_el(j))*fcsco(i)
             end do
          else
             do j=1,klass(ncsco(i))%number
                intermediate = trafo_matrix(:,:,klass(ncsco(i))%conjug_el(j))
                efm_matrix = efm_matrix +&
                     (3.0_rk*(intermediate + transpose(intermediate))+&
                     cmplx(0.0_rk,1.0_rk,RK)*(intermediate - transpose(intermediate)))*fcsco(i)*sq2
             end do
          endif
       endif
    end do

    ! now build efm for CSCO-II

    efmII_matrix = efm_matrix
    csco_sub = (0.0_rk,0.0_rk)
    factor = 1.0_rk
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups(m)%ncsco(i).ne.0) then
             if (subgroups(m)%klass(subgroups(m)%ncsco(i))%ambivalent) then
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   csco_sub = csco_sub +&
                        trafo_matrix(:,:,subgroups(m)%klass(subgroups(m)%ncsco(i))%&
                        &conjug_el(j))*subgroups(m)%fcsco(i)*factor
                end do
             else
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   intermediate = trafo_matrix(:,:,subgroups(m)%klass(subgroups(m)&
                        &%ncsco(i))%conjug_el(j))
                   csco_sub = csco_sub + subgroups(m)%fcsco(i)*factor*sq2*&
                        &(3.0_rk*(intermediate + transpose(intermediate))+&
                        &10.0_rk*cmplx(0.0_rk,1.0_rk,RK)*(intermediate - transpose(intermediate)))
                end do
             endif
          endif
       end do
       factor = factor*efm_subgroup_factor
    end do

    efmII_matrix = efmII_matrix + csco_sub*efm_group_factor/(efm_subgroup_factor**group_num_sg)

    ! now build efm for CSCO-III
    efmIII_matrix = efmII_matrix
    ! initialize csco-matrix of intrinsic subgroups
    csco_sub_intr = (0.0_rk,0.0_rk)

    factor = 1.0_rk
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups(m)%ncsco(i).ne.0) then
             if (subgroups(m)%klass(subgroups(m)%ncsco(i))%ambivalent) then
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   csco_sub_intr = csco_sub_intr +&
                        trafo_matrix_intr(:,:,subgroups(m)%klass(subgroups(m)%ncsco(i))%&
                        &conjug_el(j))*subgroups(m)%fcsco(i)*factor
                end do
             else
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   intermediate = trafo_matrix_intr(:,:,subgroups(m)%klass(subgroups(m)&
                        &%ncsco(i))%conjug_el(j))
                   csco_sub_intr = csco_sub_intr + subgroups(m)%fcsco(i)*factor*sq2*&
                        &(3.0_rk*(intermediate + transpose(intermediate))+&
                        &10.0_rk*cmplx(0.0_rk,1.0_rk,RK)*(intermediate - transpose(intermediate)))
                end do
             endif
          endif
       end do
       factor = factor*efm_subgroup_factor
    end do

    ! add csco-matrix of the intrinsic subgroups
    efmIII_matrix = efmIII_matrix + csco_sub_intr*efm_intrgroup_factor/(efm_subgroup_factor**(group_num_sg+1))

    ! solve eigenvalue problem for CSCO-III
    ierr = 0
    call ch(group_num_el,group_num_el,&               ! dimensions of transformation
         real(efmIII_matrix),aimag(efmIII_matrix),&   ! matrix of eigenvalue problem
         eigen_valueIII,matz,eig_vec_r,eig_vec_i,&    ! eigenvalues and eigenvectors
         work1,work2,work3,ierr)                      ! work arrays
    if (ierr.ne.0) then
       write(*,*) "efm_cscoIII_solve: stopped due to error in eigensolver IERR =",ierr
       call error_handler("efm_cscoIII_solve: error in eigensolver")
    end if


    ! construct complex eigenvectors
    eigen_vector = cmplx(eig_vec_r,eig_vec_i,RK)

    !
    ! --- In the following eigenvalues of the different CSCO`s as CSCO-I, CSCO-II, -----
    ! --- and the CSCO-III are determined. ---------------------------------------------
    !

    ! determine eigenvalues of CSCO (i.e. the irrep labels)

    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(efm_matrix,eigen_vector)

    ! determine the multiple
    eigen_value = 0
    do j=1,group_num_el
       do i=1,group_num_el
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value(j) = -real(intermediate(i,j)/eigen_vector(i,j))
             exit
          endif
       enddo
    enddo

    ! determine eigenvalues of CSCO of the intrinisc subgroups (i.e. the
    ! multiplicity labels)

    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(csco_sub_intr,eigen_vector)

    ! determine the multiple
    eigen_value_intr = 0.0_rk
    do j=1,group_num_el
       do i=1,group_num_el
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value_intr(j) = -real(intermediate(i,j)/eigen_vector(i,j))
             exit
          endif
       enddo
    enddo

    ! determine eigenvalues of CSCO of the subgroups (i.e. the irrep
    ! labels)

    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(csco_sub,eigen_vector)

    ! determine the multiple
    eigen_value_sub = 0.0_rk
    do j=1,group_num_el
       do i=1,group_num_el
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value_sub(j) = -real(intermediate(i,j)/eigen_vector(i,j))
             exit
          endif
       enddo
    enddo


    !
    ! ---  make the eigenvectors come real (as far as possible) -------
    !

    ! loop over all eigenvectors
    do i=1,group_num_el
       ! Normalize eigenvectors, they  are not exactly normalized from
       ! the eigensolver.  DOT_PRODUCT(X, Y)  does the right thing for
       ! complex vectors, namely the SUM(CONJG(X)  * Y). So that for X
       ! = Y the  output is real, formally. Still  the return value is
       ! complex by type constrains. So the cast to real is safe here:
       abs_factor = sqrt (real (dot_product (eigen_vector(:, i), &
                                             eigen_vector(:, i)), &
                                kind (abs_factor)))
       eigen_vector(:, i) = eigen_vector(:, i) / abs_factor
       !write(output_unit,*) "normalization of eigenvector Nr.",i," : ",abs_factor
       ! loop over basisfunctions
       do j=1,group_num_el
          abs_factor = abs(eigen_vector(j,i))
          if (abs_factor.gt.small) then
             c_factor = conjg(eigen_vector(j,i)/abs_factor)
             eigen_vector(:,i) = eigen_vector(:,i)*c_factor
             exit
          end if
       end do
    end do


    !
    ! -------------- determine the irrep matrices -------------------
    !

    ! loop over all irreps
    do i=1,group_num_ir

       ! assume eigenvectors of first row to be real
       ! This is true for one-dimensional and pseudo-2D irreps
       ! In case of multi-dimensional irreps it might or not
       first_vector_real = .true.

       ! read eigenvalue and dimension of irrep
       irrep_eigenvalue = irrep_can(i)%eigen_value
       irrep_dim = irrep_can(i)%dimension

       ! for subsequent use allocate help arrays of the suitable size
       ! they will be reallocated for the next irrep
       allocate(intermediate_irr(group_num_el,irrep_dim))
       allocate(intermediate_irr2(irrep_dim,irrep_dim))
       allocate(intermediate_irr3(irrep_dim,group_num_el))
       allocate(eigen_space(group_num_el,irrep_dim,irrep_dim))
       allocate(irrep_basis(group_num_el,irrep_dim))

       if (irrep_can(i)%pseudo) then
          ! treatment of pseudo-2D irreps
          do j=1,group_num_el
             if (abs(eigen_value(j)-irrep_eigenvalue).lt.small) then
                irrep_basis(:,1) = real(eigen_vector(:,j))
                irrep_basis(:,2) = aimag(eigen_vector(:,j))
                exit
             endif
          end do
          ! normalize basis
          do j=1,irrep_dim
             factor = real(dot_product(irrep_basis(:,j),irrep_basis(:,j)))
             factor = sqrt(1/factor)
             irrep_basis(:,j) = irrep_basis(:,j)*factor
          end do

          ! determine eigenvalues of CSCO-II for the basis of the irrep

          ! multiply real CSCO-II efm-matrix by the eigenvectors
          intermediate_irr = matmul(real(efmII_matrix),irrep_basis)

          ! determine the multiple
          eigen_valueII = 0.0_rk
          do j=1,irrep_dim
             do n=1,group_num_el
                if (abs(intermediate_irr(n,j)).gt.small) then
                   eigen_valueII(j) = -real(intermediate_irr(n,j)/irrep_basis(n,j))
                   exit
                endif
             enddo
          enddo

       else if (irrep_dim.eq.1) then
          ! treatment for one-dimensional irreps

          ! fill basis of irrep
          do j=1,group_num_el
             if (abs(eigen_value(j)-irrep_eigenvalue).lt.small) then
                irrep_basis(:,1) = eigen_vector(:,j)
                exit
             endif
          end do
          ! determine eigenvalues of CSCO-II for the basis of the irrep

          ! multiply real CSCO-II efm-matrix by the eigenvectors
          intermediate_irr = matmul(real(efmII_matrix),irrep_basis)

          ! determine the multiple
          eigen_valueII = 0.0_rk
          do j=1,irrep_dim
             do n=1,group_num_el
                if (abs(intermediate_irr(n,j)).gt.small) then
                   eigen_valueII(j) = -real(intermediate_irr(n,j)/irrep_basis(n,j))
                   exit
                endif
             enddo
          enddo

       else
          ! treatment for multi-dimensional irreps


          ! look  up for  one  eigenvector of  the  current irrep  and
          ! memorize its multiplicy label
          irrep_mult_label = HUGE (irrep_mult_label) ! sentinel
          found = .false.
          do j = 1, group_num_el
             if (abs(eigen_value(j)-irrep_eigenvalue).lt.small) then
                irrep_mult_label = eigen_value_intr(j)
                found = .true.
                exit
             endif
          end do
          ASSERT(found)

          ! determine the labels of the subspace
          n = 1
          do j=1,group_num_el
             if ((abs(eigen_value(j)-irrep_eigenvalue).lt.small).and.&
                  &(abs((eigen_value_intr(j)-irrep_mult_label)/irrep_mult_label).lt.small)) then
                irrep_eigen_value(n) = eigen_value_sub(j)
                n = n+1
             end if
          end do


          ! fill up array eigen_space
          do n=1,irrep_dim
             do k=1,irrep_dim
                do j=1,group_num_el
                   if ((abs(eigen_value(j)-irrep_eigenvalue).lt.small).and.&
                        &(abs((eigen_value_sub(j)-irrep_eigen_value(n))/irrep_eigen_value(n)).lt.small).and.&
                        &(abs((eigen_value_intr(j)-irrep_eigen_value(k))/irrep_eigen_value(k)).lt.small)) then
                      eigen_space(:,n,k) = eigen_vector(:,j)
                      exit
                   end if
                end do
             end do
          end do

          ! fill irrep basis, the basis from which the transformation properties
          ! will be obtained
          irrep_basis(:,:) = eigen_space(:,:,1)


          ! determine eigenvalues of CSCO-II for the basis of the irrep

          ! multiply real CSCO-II efm-matrix by the eigenvectors
          intermediate_irr = matmul(real(efmII_matrix),irrep_basis)

          ! determine the multiple
          eigen_valueII = 0.0_rk
          do j=1,irrep_dim
             do n=1,group_num_el
                if (abs(intermediate_irr(n,j)).gt.small) then
                   eigen_valueII(j) = -real(intermediate_irr(n,j)/irrep_basis(n,j))
                   exit
                endif
             enddo
          enddo


          !
          !--- make the basis come real --------------------------------------------------
          !

          ! settings if no complex coefficients are found
          basis_label = 0
          mult_label  = 0

          ! check if there are still complex coefficients in the basis
          ! if yes look for the conjugated element
          first_vector_real = .false.
          conjug: do j=1,irrep_dim
             do n=1,group_num_el
                if (aimag(irrep_basis(n,j)).gt.small) then
                   ! nonreal eigen_vector found
                   ! determine labeling of conjugated element
                   workc1 = irrep_basis(:,j)
                   workc2 = conjg(workc1)
                   workc1 = MATMUL(csco_sub,workc2)
                   do k=1,group_num_el
                      if (abs(workc2(k)).gt.small) then
                         eigen_sub = -real(workc1(k)/workc2(k))
                         workc1 = MATMUL(csco_sub_intr,workc2)
                         eigen_intr = -real(workc1(k)/workc2(k))
                         ! determine coordinates of conjugated element
                         labeling:     do l=1,irrep_dim
                            do m=1,irrep_dim
                               if ((abs((eigen_sub-irrep_eigen_value(l))/eigen_sub).lt.small).and.&
                                    &(abs((eigen_intr-irrep_eigen_value(m))/eigen_intr).lt.small)) then
                                  basis_label   = l
                                  mult_label    = m
                                  ! labels found --> exit loops
                                  exit labeling
                               end if
                            end do
                         end do labeling
                         ! eigen values determined and labels found --> exit loop
                         exit
                      end if
                   end do
                   exit conjug
                end if
             end do
             ! the program only can reach this point if the first vector was real. Thus:
             first_vector_real = .true.
          end do conjug



          ! if the conjugated element was not found in the initial irrep space
          ! then add the irrep space where the conjugated elements were found
          ! to the inital space

          if (mult_label.gt.1) then
             irrep_basis(:,:) = irrep_basis(:,:) + eigen_space(:,:,mult_label)
             do j=1,irrep_dim
                factor = real(dot_product(irrep_basis(:,j),irrep_basis(:,j)))
                factor = sqrt(1/factor)
                irrep_basis(:,j) = irrep_basis(:,j)*factor
             end do
          end if


          ! now irrep_basis contains its own conjugated elements,
          ! the basis can be made real

          do j=1,irrep_dim
             do n=1,group_num_el
                if (aimag(irrep_basis(n,j)).gt.small) then
                   ! nonreal eigen_vector found
                   ! determine labeling of conjugated element
                   workc1 = irrep_basis(:,j)
                   workc2 = conjg(workc1)
                   workc1 = MATMUL(csco_sub,workc2)
                   ! determine eigenvalue and label of conjugated element
                   do k=1,group_num_el
                      if (abs(workc2(k)).gt.small) then
                         eigen_sub = -real(workc1(k)/workc2(k))
                         do l=1,irrep_dim
                            if (abs(eigen_sub-irrep_eigen_value(l)).lt.small) then
                               basis_label   = l
                               ! label found --> exit loop
                               exit
                            end if
                         end do
                         ! one element is sufficient for the determination of
                         ! the eigenvalue --> exit loop
                         exit
                      end if
                   end do
                   ! now the two conjugated elements are substituted by
                   ! their real and imaginary parts
                   factor = sqrt(2.0_rk)
                   ! apply phasefactor
                   workc2 = workc2*phase_factor
                   ! now obtain the basis
                   irrep_basis(:,j) = real(workc2)*factor
                   irrep_basis(:,basis_label) = -aimag(workc2)*factor
                   ! one element was enough to see that the element was complex
                   ! --> exit loop
                   exit
                end if
             end do
          end do
       end if! if irrep_can(i)%pseudo


       !
       ! --------- determine the representations of the irreps ---------------------
       !

       ! NO MORE: workaround in order to circumvent compiler error on SP2.
       allocate(irrep_can(i)%irrep_gen_oper(irrep_dim))
       allocate(irrep_can(i)%irrep_gen_coeff(irrep_dim, irrep_dim))
       allocate(irrep_can(i)%irrep_matrix(irrep_dim, irrep_dim, group_num_el))


       ! initialize irrep basis generator
       irrep_can(i)%irrep_gen_coeff = 0.0_rk
       do j=1,irrep_dim
          irrep_can(i)%irrep_gen_coeff(j,j) = 1.0_rk
       end do

       ! now determine the irrep matrix from the formula
       ! irrep matrix = C^T * TRAFO * C, where the columns of C are the
       ! basis vectors of the irrep
       do j=1,group_num_el
          intermediate_irr  = matmul(trafo_matrix(:,:,j),irrep_basis)
          intermediate_irr2 = matmul(conjg(transpose(irrep_basis)), intermediate_irr)

          ! write irrep matrix to its destination
          irrep_can(i)%irrep_matrix(:,:,j) = real(intermediate_irr2)

       end do

       ! determine eigenvector of the basisvector from which the eigenfunctions
       ! are generated
       irrep_can(i)%gen_eigen_value  = eigen_valueII(1)
       irrep_can(i)%build_projector = .false.
       if (.not.first_vector_real) then
          irrep_can(i)%build_projector =.true.
          irrep_can(i)%gen_eigen_value = irrep_can(i)%gen_eigen_value - 0.314_rk
       end if

       ! determine generators of the irrep basis
       ! collect the first columns of all irrep matrices, which are the pictures
       ! of the first basis vector for the corresponding group elements
       do j=1,group_num_el
          intermediate_irr3(:,j) = irrep_can(i)%irrep_matrix(:,1,j)
       end do
       ! determine set of irrep_dim operators, which generate the irrep
       call efm_solve_nonsq(irrep_dim,group_num_el,intermediate_irr3(:,:),&
            irrep_can(i)%irrep_gen_coeff,irrep_can(i)%irrep_gen_oper)

       ! deallocate all help arrays for this irrep
       deallocate(intermediate_irr,intermediate_irr2,intermediate_irr3)
       deallocate(irrep_basis,eigen_space)
    end do


    ! deallocate local arrays
    deallocate(efm_matrix,efmII_matrix,efmIII_matrix)
    deallocate(csco_sub_intr,csco_sub,intermediate)
    deallocate(eigen_value,eigen_valueII,eigen_valueIII,eigen_value_intr)
    deallocate(eigen_value_sub,eigen_vector)
    deallocate(eig_vec_r,eig_vec_i)
    deallocate(irrep_eigen_value)

    ! deallocate work arrays
    deallocate(work1,work2,work3,workc1,workc2)

  end subroutine efm_cscoIII_solve
  !*************************************************************

  !*************************************************************
  subroutine efm_proj_cscoIII_solve(trafo_matrix,trafo_matrix_intr)
    !  Purpose: solves eigenequation for the CSCO-III of the
    !           regular representation of the projective group,
    !           thus determines irreducible representations of
    !           the double group.
    !------------ Modules used -----------------------------------
    use type_module
    use iounitadmin_module
    use group_module
    implicit none
    !------------ Declaration of formal parameters  -----------------
    real(RK),    intent(in)           :: trafo_matrix(:,:,:)
    ! transformation matrices of the group operation in the regular representation
    real(RK),    intent(in)           :: trafo_matrix_intr(:,:,:)
    ! transformation matrices of the intrinsic group operation in the regular
    ! representation
    !** End of interface *****************************************

    !------------ Declaration of local parameters ---------------
    integer(IK),parameter                :: matz = 1
    ! parameter for eigensolver
    real(RK),parameter                   :: small = 0.00001_rk
    ! very small value

    !------------ Declaration of local variables  ---------------
    complex(CK),allocatable             :: efm_matrix(:,:),efmII_matrix(:,:)
    ! contains matrix of eigenvalue problem for CSCO and CSCO-II
    complex(CK),allocatable             :: efmIII_matrix(:,:)
    ! contains matrix of eigenvalue problem for CSCO-III
    complex(CK),allocatable             :: csco_sub_intr(:,:)
    ! csco-matrix of the intrinsic group (labels multiple irreps)
    complex(CK),allocatable             :: csco_sub(:,:)
    ! csco-matrix of the group (labels irreps)
    complex(CK),allocatable             :: intermediate(:,:)
    ! work array of dimension group_num_el x group_num_el
    complex(CK),allocatable             :: intermediate_irr(:,:)
    ! work array of dimension group_num_el x irrep_dim
    complex(CK),allocatable             :: intermediate_irr2(:,:)
    ! work array of dimension irrep_dim x irrep_dim
    complex(CK),allocatable             :: intermediate_irr3(:,:)
    ! work array of dimension irrep_dim x group_num_el
    real(RK),allocatable                 :: eigen_value(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable                 :: eigen_valueII(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable                 :: eigen_valueIII(:)
    ! eigenvalues of CSCO-III
    real(RK),allocatable                 :: eigen_value_intr(:)
    ! eigenvalues of the CSCO of the intrinsic subgroups
    real(RK),allocatable                 :: eigen_value_sub(:)
    ! eigenvalues of the CSCO of the subgroups
    complex(CK),allocatable             :: eigen_vector(:,:)
    ! eigenvectors of CSCO-III and CSCO-II and CSCO
    ! eigen_vector(group_num_el,group_num_el)
    ! eigen_vector(basis_function,eigen_vector_number)
    real(RK),allocatable                 :: eig_vec_r(:,:),eig_vec_i(:,:)
    ! real and imaginary parts of eigen_vector
    complex(CK),allocatable             :: irrep_basis(:,:)
    ! basis for one irrep: irrep_basis(group_num_el,irrep_dim)
    complex(CK),allocatable             :: eigen_space(:,:,:)
    ! eigen_space of a certain irrep
    real(RK),allocatable                 :: irrep_eigen_value(:)
    ! eigen_values of irrep
    real(RK),allocatable                 :: work1(:),work2(:),work3(:,:)
    ! real work arrays
    complex(CK),allocatable             :: workc1(:),workc2(:)
    ! complex work arrays
    integer(IK) :: i, j, k, m, n, ierr
    ! counters
    integer(IK)                          :: irrep_dim
    ! help variable for dimension of irrep
    real(RK)                             :: factor,sq2
    ! factor for subgroups
    real(RK)                             :: abs_factor
    ! factor for making real of eigenvectors
    complex(CK)                         :: c_factor
    ! complex factor for making real of eigenvectors
    complex(CK)                         :: phase_factor
    ! complex phase factor which adjusts two dimensional irreps to
    ! a special transformation (only for the groups Td,O,Oh and I,Ih)
    real(RK)                             :: irrep_eigenvalue, irrep_mult_label
    ! irrep and irrep multiplicity labels
    complex(CK)                         :: eigen_val
    ! for debugging
    logical                                        :: strategy_one_equation
    ! if .true. only one CSCO-III matrix is used
    ! if .false. first Eigenfunctions of CSCO, then of CSCO-II intr and of CSCO-II
    ! are determined
    complex(CK),allocatable            :: eigen_space_CSCO(:,:), &
         CSCOIIintr_eigen_space(:,:),eigen_vector_CII_transf(:,:),eigen_vector_CII(:,:),CSCOII_irrep_sub(:,:)
    ! work arrays for strategy with succesive Diagonalizations

    logical :: found

    !------------ Executable code ------------------------------------

    ! allocate local arrays

    allocate(efm_matrix(group_num_el,group_num_el))
    allocate(efmII_matrix(group_num_el,group_num_el))
    allocate(efmIII_matrix(group_num_el,group_num_el))
    allocate(csco_sub_intr(group_num_el,group_num_el),csco_sub(group_num_el,group_num_el))
    allocate(intermediate(group_num_el,group_num_el))
    allocate(eigen_value(group_num_el))
    allocate(eigen_valueII(group_num_el))
    allocate(eigen_valueIII(group_num_el))
    allocate(eigen_vector(group_num_el,group_num_el))
    allocate(eig_vec_r(group_num_el,group_num_el),eig_vec_i(group_num_el,group_num_el))
    allocate(irrep_eigen_value(group_num_el))

    ! allocate work arrays

    allocate(work1(group_num_el),work2(group_num_el))
    allocate(workc1(group_num_el),workc2(group_num_el))
    allocate(work3(2,group_num_el))

    !
    ! determine phasefactor (only needed for polyhedral irreps)
    !

    !select case(group_name)
    !case('Td  ','O   ','Oh  ')
    !   pi = acos(-1.0_rk)*2.0_rk/3.0_rk
    !   phase_factor = cmplx(cos(pi),sin(pi),RK)
    !case('I   ','Ih  ')
    !   phase_factor = (0.6123724356958248_rk,0.7905694150420700_rk)
    !case default
       phase_factor = (1.0_rk,0.0_rk)
    !end select

    !
    ! determine factors for EFM-Matrix
    !

    select case(group_name)
    case('Td  ','O   ','Oh  ')
       efm_proj_subgroup_factor = 5.0_rk
       efm_proj_group_factor = 1.0_rk
       efm_proj_intrgroup_factor = 2.0_rk
       strategy_one_equation = .true.
    case('I   ','IH  ')
       efm_proj_subgroup_factor = 5.0_rk
       efm_proj_group_factor = 6.0_rk
       efm_proj_intrgroup_factor = 10.0_rk
       strategy_one_equation = .false.
    case default
       efm_proj_subgroup_factor = 10.0_rk
       !efm_group_factor = 4.0_rk
       efm_proj_group_factor = 8.0_rk
       efm_proj_intrgroup_factor = 12.0_rk
       strategy_one_equation = .true.
    end select

    !
    ! construct EFM-matrix of the CSCO-III
    !

    ! initialize efm-matrix
    efm_matrix = (0.0_rk,0.0_rk)

    ! prefactor
    sq2 = 1/(2.0_rk*sqrt(2.0_rk))

    ! first build efm for CSCO
    do i = 1,3
       if (ncsco_proj(i).ne.0) then
          if (klass(ncsco_proj(i))%ambivalent_proj) then
             do j=1,klass(ncsco_proj(i))%number
                efm_matrix = efm_matrix +&
                     trafo_matrix(:,:,klass(ncsco_proj(i))%conjug_el(j))*fcsco_proj(i)
             end do
          else
             do j=1,klass(ncsco_proj(i))%number
                intermediate = trafo_matrix(:,:,klass(ncsco_proj(i))%conjug_el(j))
                efm_matrix = efm_matrix +&
                     (3.0_rk*(intermediate + transpose(intermediate))+&
                     cmplx(0.0_rk,1.0_rk,RK)*(intermediate - transpose(intermediate)))*fcsco_proj(i)*sq2
             end do
          endif
       endif
    end do

    ! now build efm for CSCO-II

    efmII_matrix = efm_matrix
    ! initialize efmII-matrix for debugging
    !efmII_matrix = (0.0_rk,0.0_rk)

    csco_sub = (0.0_rk,0.0_rk)
    factor = 1.0_rk
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups_proj(m)%ncsco_proj(i).ne.0) then
             if (subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%ambivalent_proj) then
                do j=1,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%number
                   csco_sub = csco_sub +&
                        trafo_matrix(:,:,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%&
                        &conjug_el(j))*subgroups_proj(m)%fcsco_proj(i)*factor
                end do
             else
                do j=1,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%number
                   intermediate = trafo_matrix(:,:,subgroups_proj(m)%klass(subgroups_proj(m)&
                        &%ncsco_proj(i))%conjug_el(j))
                   csco_sub = csco_sub + subgroups_proj(m)%fcsco_proj(i)*factor*sq2*&
                        &(3.0_rk*(intermediate + transpose(intermediate))+&
                        &10.0_rk*cmplx(0.0_rk,1.0_rk,RK)*(intermediate - transpose(intermediate)))
                end do
             endif
          endif
       end do
       factor = factor*efm_proj_subgroup_factor
    end do

    ! show csco sub matrix
    !write(output_unit,*) "--------------------------------------------"
    !write(output_unit,*) "         CSCO sub  matrix"
    !write(output_unit,*) "--------------------------------------------"
    !do n=1,group_num_el
    !         write(output_unit,603)&
    !              &         (csco_sub(n,m), m = 1,group_num_el)
    !end do


    ! initialize csco-matrix of intrinsic subgroups
    csco_sub_intr = (0.0_rk,0.0_rk)

    factor = 1.0_rk
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups_proj(m)%ncsco_proj(i).ne.0) then
             if (subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%ambivalent_proj) then
                do j=1,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%number
                   csco_sub_intr = csco_sub_intr +&
                        trafo_matrix_intr(:,:,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%&
                        &conjug_el(j))*subgroups_proj(m)%fcsco_proj(i)*factor
                end do
             else
                do j=1,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%number
                   intermediate = trafo_matrix_intr(:,:,subgroups_proj(m)%klass(subgroups_proj(m)&
                        &%ncsco_proj(i))%conjug_el(j))
                   csco_sub_intr = csco_sub_intr + subgroups_proj(m)%fcsco_proj(i)*factor*sq2*&
                        &(3.0_rk*(intermediate + transpose(intermediate))+&
                        &10.0_rk*cmplx(0.0_rk,1.0_rk,RK)*(intermediate - transpose(intermediate)))
                end do
             endif
          endif
       end do
       factor = factor*efm_proj_subgroup_factor
    end do

    ! show csco sub intr matrix
    !write(output_unit,*) "--------------------------------------------"
    !write(output_unit,*) "         CSCO sub_intr  matrix"
    !write(output_unit,*) "--------------------------------------------"
    !do n=1,group_num_el
    !         write(output_unit,603)&
    !              &         (csco_sub_intr(n,m), m = 1,group_num_el)
    !end do

    ! build up CSCO-II
    efmII_matrix = efm_matrix + csco_sub*efm_proj_group_factor/(efm_proj_subgroup_factor**group_num_sg)

    ! build up CSCO-III as one matrix
    efmIII_matrix = efmII_matrix + csco_sub_intr*efm_proj_intrgroup_factor/(efm_proj_subgroup_factor**(group_num_sg+1))



    ! show efmIII matrix
    !write(output_unit,*) "--------------------------------------------"
    !write(output_unit,*) "         CSCO-III matrix"
    !write(output_unit,*) "--------------------------------------------"
    !do n=1,group_num_el
    !         write(output_unit,603)&
    !              &         (efmIII_matrix(n,m), m = 1,group_num_el)
    !end do

    if (strategy_one_equation) then
       ! allocate eigenvalues of CSCO-II and CSCO-II (intr)
       allocate(eigen_value_intr(group_num_el),eigen_value_sub(group_num_el))

       ! solve eigenvalue problem for CSCO-III
       ierr = 0
       call ch(group_num_el,group_num_el,&               ! dimensions of transformation
            real(efmIII_matrix),aimag(efmIII_matrix),&   ! matrix of eigenvalue problem
            eigen_valueIII,matz,eig_vec_r,eig_vec_i,&    ! eigenvalues and eigenvectors
            work1,work2,work3,ierr)                      ! work arrays
       if (ierr.ne.0) then
          write(*,*) "efm_proj_cscoIII_solve: stopped due to error in eigensolver IERR =",ierr
          call error_handler("efm_proj_cscoIII_solve: error in eigensolver")
       end if


       ! construct complex eigenvectors
       eigen_vector = cmplx(eig_vec_r,eig_vec_i,RK)

       !
       ! --- In the following eigenvalues of the different CSCO`s as CSCO-I, CSCO-II, -----
       ! --- and the CSCO-III are determined. ---------------------------------------------
       !

       ! determine eigenvalues of CSCO (i.e. the irrep labels)

       ! multiply CSCO efm-matrix by the eigenvectors
       intermediate = matmul(efm_matrix,eigen_vector)

       ! determine the multiple
       eigen_value = 0
       do j=1,group_num_el
          do i=1,group_num_el
             if ((abs(intermediate(i,j)).gt.small).and.(abs(eigen_vector(i,j)).gt.small)) then
                eigen_value(j) = -real(intermediate(i,j)/eigen_vector(i,j))
                exit
             endif
          enddo
       enddo

       ! determine eigenvalues of CSCO of the intrinsic subgroups (i.e. the multiplicity
       ! labels)

       ! multiply CSCO efm-matrix by the eigenvectors
       intermediate = matmul(csco_sub_intr,eigen_vector)

       ! determine the multiple
       eigen_value_intr = 0.0_rk
       do j=1,group_num_el
          do i=1,group_num_el
             if ((abs(intermediate(i,j)).gt.small).and.(abs(eigen_vector(i,j)).gt.small)) then
                !      if (eigen_val_old.ne.0_rk) then
                !         if (abs(eigen_val_old-eigen_val).gt.small) then
                !            write(output_unit,*) "eigenvector ",j," is not eigenvector of csco_sub_intr"
                !         else
                !            write(output_unit,*) "eigenvector ",j," is an eigenvector of csco_sub_intr"
                !         endif
                !      end if
                !      eigen_val_old = eigen_val
                eigen_value_intr(j) = -real(intermediate(i,j)/eigen_vector(i,j))
                exit
             endif
          enddo
          !eigen_value_intr(j) = -real(eigen_val)
       enddo

       ! determine eigenvalues of CSCO of the subgroups (i.e. the irrep
       ! labels)

       ! multiply CSCO efm-matrix by the eigenvectors
       intermediate = matmul(csco_sub,eigen_vector)

       ! determine the multiple
       eigen_value_sub = 0.0_rk
       do j=1,group_num_el
          !eigen_val_old = 0.0_rk
          do i=1,group_num_el
             if ((abs(intermediate(i,j)).gt.small).and.(abs(eigen_vector(i,j)).gt.small)) then
                !eigen_val = intermediate(i,j)/eigen_vector(i,j)
                !write(output_unit,*) "eigenvalue ",eigen_val
                !if (eigen_val_old.ne.0_rk) then
                !   if (abs(eigen_val_old-eigen_val).gt.small) then
                !      write(output_unit,*) "eigenvector ",j," is not eigenvector of csco_sub"
                !   else
                !      write(output_unit,*) "eigenvector ",j," is an eigenvector of csco_sub"
                !   endif
                !end if
                !eigen_val_old = eigen_val
                eigen_value_sub(j) = -real(intermediate(i,j)/eigen_vector(i,j))
                exit
             endif
          enddo
          !eigen_value_sub(j) = -real(eigen_val)
       enddo

    else
       !
       ! now determine eigenvectors of CSCO
       !
       ierr = 0
       call ch(group_num_el,group_num_el,&               ! dimensions of transformation
            real(efm_matrix),aimag(efm_matrix),&         ! matrix of eigenvalue problem
            eigen_value,matz,eig_vec_r,eig_vec_i,&    ! eigenvalues and eigenvectors
            work1,work2,work3,ierr)                      ! work arrays
       if (ierr.ne.0) then
          write(*,*) "efm_proj_cscoIII_solve: stopped due to error in eigensolver IERR =",ierr
          call error_handler("efm_proj_cscoIII_solve: error in eigensolver")
       end if


       ! construct complex eigenvectors
       eigen_vector = cmplx(eig_vec_r,eig_vec_i,RK)

       ! due to some convention
       eigen_value = -eigen_value

       ! show solutions
       !write(output_unit,*) "--------------------------------------------"
       !write(output_unit,*) "         eigenvectors of CSCO (rows)"
       !write(output_unit,*) "--------------------------------------------"
       !do n=1,group_num_el
       !   write(output_unit,603)&
       !        &         (eigen_vector(n,m), m = 1,group_num_el)
       !end do


    end if ! strategy_one_equation

    ! deallocate work arrays
    deallocate(work1,work2,work3,workc1,workc2)
    deallocate(eig_vec_r,eig_vec_i)

    !
    ! test if csco_sub and csco_intr commute
    !

    !   intermediate = matmul(efm_matrix,csco_sub_intr)-matmul(csco_sub_intr,efm_matrix)
    !if (maxval(abs(intermediate)).gt.small) then
    !   write(output_unit,*) "csco and csco_sub_intr do not commute"
    !else
    !   write(output_unit,*) "csco and csco_sub_intr do commute"
    !end if

    !intermediate = matmul(efm_matrix,csco_sub)-matmul(csco_sub,efm_matrix)
    ! show commutator
    !write(output_unit,*) "--------------------------------------------"
    !write(output_unit,*) "         [CSCO,CSCO_sub]"
    !write(output_unit,*) "--------------------------------------------"
    !do n=1,group_num_el
    !         write(output_unit,603)&
    !              &         (intermediate(n,m), m = 1,group_num_el)
    !end do
    !if (maxval(abs(intermediate)).gt.small) then
    !   write(output_unit,*) "csco and csco_sub do not commute"
    !else
    !   write(output_unit,*) "csco and csco_sub do commute"
    !end if

    !intermediate = matmul(csco_sub,csco_sub_intr)-matmul(csco_sub_intr,csco_sub)
    !if (maxval(abs(intermediate)).gt.small) then
    !   write(output_unit,*) "csco_sub and csco_sub_intr do not commute"
    !else
    !   write(output_unit,*) "csco_sub and csco_sub_intr do commute"
    !end if

    !do i=1,group_num_el
    !   intermediate = matmul(trafo_matrix(:,:,i),csco_sub_intr)-matmul(csco_sub_intr,trafo_matrix(:,:,i))
    !   if (maxval(abs(intermediate)).gt.small) then
    !      write(output_unit,*) "trafo ",i," and csco_sub_intr do not commute"
    !   else
    !      write(output_unit,*) "trafo ",i," and csco_sub_intr do commute"
    !   end if
    !end do

    !do i=1,group_num_el
    !   intermediate = matmul(trafo_matrix(:,:,i),csco_sub)-matmul(csco_sub,trafo_matrix(:,:,i))
    !   if (maxval(abs(intermediate)).gt.small) then
    !      write(output_unit,*) "trafo ",i," and csco_sub do not commute"
    !   else
    !      write(output_unit,*) "trafo ",i," and csco_sub do commute"
    !   end if
    !end do

    !
    ! ---  make the eigenvectors come real (as far as possible) -------
    !

    ! loop over all eigenvectors
    do i=1,group_num_el
       ! Normalize eigenvectors, they  are not exactly normalized from
       ! the eigensolver.  DOT_PRODUCT(X, Y)  does the right thing for
       ! complex vectors, namely the SUM(CONJG(X)  * Y). So that for X
       ! = Y the  output is real, formally. Still  the return value is
       ! complex by type constrains. So the cast to real is safe here:
       abs_factor = sqrt (real (dot_product (eigen_vector(:, i), &
                                             eigen_vector(:, i)), &
                                kind (abs_factor)))
       eigen_vector(:, i) = eigen_vector(:, i) / abs_factor
       !write(output_unit,*) "normalization of eigenvector Nr.",i," : ",abs_factor
       ! loop over basisfunctions
       do j=1,group_num_el
          abs_factor = abs(eigen_vector(j,i))
          if (abs_factor.gt.small) then
             c_factor = conjg(eigen_vector(j,i)/abs_factor)
             eigen_vector(:,i) = eigen_vector(:,i)*c_factor
             exit
          end if
       end do
    end do
    !
    ! -------------- determine the irrep matrices -------------------
    !

    ! loop over all (projective) irreps
    do i=1,group_num_re


       !am: print*,"processing # irrep ",i
       ! read eigenvalue and dimension of irrep
       irrep_eigenvalue = proj_irrep_can(i)%eigen_value
       irrep_dim = proj_irrep_can(i)%dimension

       ! allocate basis of irreps
       allocate(irrep_basis(group_num_el,irrep_dim))

       !am: print*,"now determine basis of irrep"
       !
       ! Determine basis of irrep
       !

       if (strategy_one_equation) then
          ! allocate eigenspace of CSCO
          allocate(eigen_space(group_num_el,irrep_dim,irrep_dim))
          if (irrep_dim.eq.1) then
             ! treatment for one-dimensional irreps

             !am: print*,"  processing one dimensional irrep"
             ! fill basis of irrep
             do j=1,group_num_el
                if (abs(eigen_value(j)-irrep_eigenvalue).lt.small) then
                   irrep_basis(:,1) = eigen_vector(:,j)
                   exit
                endif
             end do

          else
             ! treatment for multi-dimensional irreps

             ! look up  for one eigenvector  of the current  irrep and
             ! memorize its multiplicy label
             irrep_mult_label = HUGE (irrep_mult_label) ! sentinel
             found = .false.
             do j = 1, group_num_el
                if (abs(eigen_value(j)-irrep_eigenvalue).lt.small) then
                   irrep_mult_label = eigen_value_intr(j)
                   found = .true.
                   exit
                endif
             end do
             ASSERT(found)

             ! determine the labels of the subspace
             n = 1
             do j=1,group_num_el
                if ((abs(eigen_value(j)-irrep_eigenvalue).lt.small).and.&
                     &(abs((eigen_value_intr(j)-irrep_mult_label)/irrep_mult_label).lt.small)) then
                   irrep_eigen_value(n) = eigen_value_sub(j)
                   n = n+1
                end if
             end do

             ! fill up array eigen_space
             do n=1,irrep_dim
                do k=1,irrep_dim
                   do j=1,group_num_el
                      if ((abs(eigen_value(j)-irrep_eigenvalue).lt.small).and.&
                           &(abs((eigen_value_sub(j)-irrep_eigen_value(n))/irrep_eigen_value(n)).lt.small).and.&
                           &(abs((eigen_value_intr(j)-irrep_eigen_value(k))/irrep_eigen_value(k)).lt.small)) then
                         eigen_space(:,n,k) = eigen_vector(:,j)
                         exit
                      end if
                   end do
                end do
             end do

             ! fill irrep basis, the basis from which the transformation properties
             ! will be obtained
             irrep_basis(:,:) = eigen_space(:,:,1)

          end if ! irrep_dim.eq.1

          ! deallocate work arrays
          deallocate(eigen_space)
       else
          !
          ! now determine irrep eigenspace
          !
          !am: print*," succesive equations"
          ! allocate eigen_space of CSCO
          allocate(eigen_space_CSCO(group_num_el,irrep_dim**2))
          n = 1
          do j=1,group_num_el
             if (abs(eigen_value(j)-irrep_eigenvalue).lt.small) then
                eigen_space_CSCO(:,n) = eigen_vector(:,j)
                n = n+1
             endif
          end do
          !am: print*," now diagonalize CSCO-II (intr)"
          ! allocate CSCO-II (transformed to eigenspace of CSCO)
          allocate(CSCOIIintr_eigen_space(irrep_dim**2,irrep_dim**2))
          ! allocate eigenvalues of CSCO-II and CSCO-II (intr)
          allocate(eigen_value_intr(irrep_dim**2),eigen_value_sub(irrep_dim))
          ! allocate help array
          allocate(intermediate_irr(group_num_el,irrep_dim**2))

          !am: print*," now transform CSCO-II (intr) to eigenspace of CSCO"
          ! transform CSCO-II (intr) to eigenspace of CSCO
          !CSCOIIintr_eigen_space = matmul(transpose(conjg(eigen_space_CSCO)),matmul(csco_sub_intr,eigen_space_CSCO))
          intermediate_irr = matmul(csco_sub_intr,eigen_space_CSCO)
          CSCOIIintr_eigen_space = matmul(transpose(conjg(eigen_space_CSCO)),intermediate_irr)

          ! deallocate help array
          deallocate(intermediate_irr)
          !am: print*," now allocate arrays"
          !
          ! now diagonalize CSCO-II (intr) in eigenspace of CSCO
          !
          ! allocate work arrays
          allocate(work1(irrep_dim**2),work2(irrep_dim**2),work3(2,irrep_dim**2))
          ! allocate real and imaginary part of solution
          allocate(eig_vec_r(irrep_dim**2,irrep_dim**2),eig_vec_i(irrep_dim**2,irrep_dim**2))
          ! allocate eigenvectors
          allocate(eigen_vector_CII_transf(irrep_dim**2,irrep_dim**2),eigen_vector_CII(group_num_el,irrep_dim**2))

          !am: print*,"now solve equation"
          ierr = 0
          call ch(irrep_dim**2,irrep_dim**2,&                                 ! dimensions of transformation
               real(CSCOIIintr_eigen_space),aimag(CSCOIIintr_eigen_space),&   ! matrix of eigenvalue problem
               eigen_value_intr,matz,eig_vec_r,eig_vec_i,&                      ! eigenvalues and eigenvectors
               work1,work2,work3,ierr)                                        ! work arrays
          if (ierr.ne.0) then
             write(*,*) "efm_proj_cscoIII_solve: stopped due to error in eigensolver IERR =",ierr
             call error_handler("efm_proj_cscoIII_solve: error in eigensolver")
          end if


          !am: print*,"now build up solution in complex space"
          ! construct complex eigenvectors
          eigen_vector_CII_transf = cmplx(eig_vec_r,eig_vec_i,RK)

          ! deallocate work arrays
          deallocate(work1,work2,work3)
          ! deallocate real and imaginary part of solution
          deallocate(eig_vec_r,eig_vec_i)

          !am: print*,"now transform back to space of symmetry elements"
          ! transform eigenvectors back to space of symmetry elements
          eigen_vector_CII   = matmul(eigen_space_CSCO,eigen_vector_CII_transf)

          ! deallocate eigen_vector_CII_transf
          deallocate(eigen_vector_CII_transf)
          ! deallocate eigen_space of CSCO
          deallocate(eigen_space_CSCO,CSCOIIintr_eigen_space)

          !
          ! Determine One Irreducible Subspace
          !

          ! allocate (irreducible) subspace of CSCO
          allocate(eigen_space_CSCO(group_num_el,irrep_dim))
          ! now assemble one irreducible eigenspace
          irrep_mult_label = eigen_value_intr(1)
          n=1
          do j=1,irrep_dim**2
             if (abs(eigen_value_intr(j)-irrep_mult_label).lt.small) then
                eigen_space_CSCO(:,n) = eigen_vector_CII(:,j)
                n=n+1
             end if
          end do

          ! deallocate eigenvectors
          deallocate(eigen_vector_CII)

          !
          ! Now determine Basis of Irrep
          !
          ! allocate irreducible invariant subspace
          allocate(CSCOII_irrep_sub(irrep_dim,irrep_dim))
          ! allocate help array
          allocate(intermediate_irr(group_num_el,irrep_dim))

          ! transform CSCO-II to irreducible invariant subspace
          !CSCOII_irrep_sub = matmul(transpose(conjg(eigen_space_CSCO)),matmul(csco_sub,eigen_space_CSCO))
          intermediate_irr = matmul(csco_sub,eigen_space_CSCO)
          CSCOII_irrep_sub = matmul(transpose(conjg(eigen_space_CSCO)),intermediate_irr)

          ! deallocate help array
          deallocate(intermediate_irr)

          !
          ! now diagonalize CSCO-II in irreducible invariant eigenspace
          !
          ! allocate eigenvectors
          allocate(eigen_vector_CII_transf(irrep_dim,irrep_dim))
          ! allocate work arrays
          allocate(work1(irrep_dim),work2(irrep_dim),work3(2,irrep_dim))
          ! allocate real and imaginary part of solution
          allocate(eig_vec_r(irrep_dim,irrep_dim),eig_vec_i(irrep_dim,irrep_dim))

          ierr = 0
          call ch(irrep_dim,irrep_dim,&                                 ! dimensions of transformation
               real(CSCOII_irrep_sub),aimag(CSCOII_irrep_sub),&         ! matrix of eigenvalue problem
               eigen_value_sub,matz,eig_vec_r,eig_vec_i,&               ! eigenvalues and eigenvectors
               work1,work2,work3,ierr)                                  ! work arrays
          if (ierr.ne.0) then
             write(*,*) "efm_proj_cscoIII_solve: stopped due to error in eigensolver IERR =",ierr
             call error_handler("efm_proj_cscoIII_solve: error in eigensolver")
          end if


          ! construct complex eigenvectors
          eigen_vector_CII_transf = cmplx(eig_vec_r,eig_vec_i,RK)

          ! transform eigenvectors back to space of symmetry elements
          ! ==> Basis of Irrep results
          irrep_basis   = matmul(eigen_space_CSCO,eigen_vector_CII_transf)


          ! deallocate irreducible invariant subspace
          deallocate(CSCOII_irrep_sub)

          ! deallocate work arrays
          deallocate(work1,work2,work3)

          ! deallocate real and imaginary part of solution
          deallocate(eig_vec_r,eig_vec_i)

          ! deallocate eigenvectors
          deallocate(eigen_vector_CII_transf)

          ! deallocate (irreducible) eigen_space of CSCO
          deallocate(eigen_space_CSCO)

          ! deallocate eigenvalues of CSCO-II and CSCO-II (intr)
          deallocate(eigen_value_intr,eigen_value_sub)
       endif

       ! for subsequent use allocate help arrays of the suitable size
       ! they will be reallocated for the next irrep
       allocate(intermediate_irr(group_num_el,irrep_dim))
       allocate(intermediate_irr2(irrep_dim,irrep_dim))
       allocate(intermediate_irr3(irrep_dim,group_num_el))

       ! determine eigenvalues of CSCO-II for the basis of the irrep

       ! multiply real CSCO-II efm-matrix by the eigenvectors
       intermediate_irr = matmul(efmII_matrix,irrep_basis)

       ! determine the multiple
       eigen_valueII = 0.0_rk
       do j=1,irrep_dim
          do n=1,group_num_el
             if (abs(intermediate_irr(n,j)).gt.small) then
                eigen_val = intermediate_irr(n,j)/irrep_basis(n,j)
                if (output_unit > 0) then
                   write(output_unit,*) "eigenvalue ",j," of CSCO-II ",eigen_val
                endif
                eigen_valueII(j) = -real(intermediate_irr(n,j)/irrep_basis(n,j))
                exit
             endif
          enddo
       enddo

       !
       ! --------- determine the representations of the irreps ---------------------
       !

       ! NO MORE: workaround in order to circumvent compiler error on SP2.
       allocate(proj_irrep_can(i)%irrep_gen_oper(irrep_dim))
       allocate(proj_irrep_can(i)%irrep_gen_coeff(irrep_dim, irrep_dim))
       allocate(proj_irrep_can(i)%irrep_matrix(irrep_dim, irrep_dim, group_num_el))

       ! initialize irrep basis generator
       proj_irrep_can(i)%irrep_gen_coeff = 0.0_rk
       do j=1,irrep_dim
          proj_irrep_can(i)%irrep_gen_coeff(j,j) = 1.0_rk
       end do

       ! now determine the irrep matrix from the formula
       ! irrep matrix = C^T * TRAFO * C, where the columns of C are the
       ! basis vectors of the irrep
       do j=1,group_num_el
          intermediate_irr  = matmul(trafo_matrix(:,:,j),irrep_basis)
          intermediate_irr2 = matmul(conjg(transpose(irrep_basis)),&
               intermediate_irr)

          ! write irrep matrix to its destination
          proj_irrep_can(i)%irrep_matrix(:,:,j) = intermediate_irr2

       end do

       !am: print*,"  determine generator"
       ! determine eigenvector of the basisvector from which the eigenfunctions
       ! are generated
       proj_irrep_can(i)%gen_eigen_value  = eigen_valueII(1)
       proj_irrep_can(i)%build_projector = .false.

       ! determine generators of the irrep basis
       ! collect the first columns of all irrep matrices, which are the pictures
       ! of the first basis vector for the corresponding group elements
       do j=1,group_num_el
          intermediate_irr3(:,j) = proj_irrep_can(i)%irrep_matrix(:,1,j)
       end do
       ! determine set of irrep_dim operators, which generate the irrep
       call efm_solve_nonsq_complex(irrep_dim,group_num_el,intermediate_irr3(:,:),&
            proj_irrep_can(i)%irrep_gen_coeff,proj_irrep_can(i)%irrep_gen_oper)


       ! deallocate all help arrays for this irrep
       deallocate(irrep_basis)
       deallocate(intermediate_irr,intermediate_irr2,intermediate_irr3)

    end do

    if (strategy_one_equation) then
       ! deallocate eigenvalues of CSCO-II and CSCO-II (intr)
       deallocate(eigen_value_intr,eigen_value_sub)
    end if

! 603 format(1x,f8.3,1x,f8.3,2x,8(f8.3,f8.3:)/8(10x,4(f8.3,f8.3:)/))

    ! deallocate local arrays
    deallocate(efm_matrix,efmII_matrix,efmIII_matrix)
    deallocate(csco_sub_intr,csco_sub,intermediate)
    deallocate(eigen_value,eigen_valueII,eigen_valueIII,eigen_vector)
    deallocate(irrep_eigen_value)


  end subroutine efm_proj_cscoIII_solve
  !*************************************************************


  !*************************************************************
  subroutine efm_cscoII_solve(trafo_matrix,dim_trafo,eigen_space_can)
    !  Purpose: solves eigenequation for a special
    !           CSCO-II
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    implicit none
    !------------ Declaration of formal parameters  -----------------
    integer(IK), intent(in)           :: dim_trafo
    ! dimension of the transformation matrix
    real(RK),    intent(in)           :: trafo_matrix(:,:,:)
    ! transformation matrices of the invariant eigenspace
    type(efm_cscoII_diag_invspace), intent(inout) :: eigen_space_can
    ! canonically ordered eigenvectors of cscoII
    !** End of interface *****************************************
    !------------ Declaration of local parameters ---------------
    integer(IK),parameter          :: matz = 1
    ! parameter for eigensolver
    real(RK),parameter             :: small = 0.0000001_rk
    ! very small value
    real(RK),parameter             :: smaller = 1.0e-12_rk
    ! smaller value
    real(RK),parameter             :: even_smaller = 1.0e-13_rk
    ! even smaller value

    !------------ Declaration of local variables  -------------------
    real(RK),allocatable             :: efm_matrix(:,:),efmII_matrix(:,:)
    ! contains matrix of eigenvalue problem for CSCO and CSCO-II
    complex(CK),allocatable         :: efm_matrix_c(:,:)
    ! complex matrix of the eigenvalue problem for the CSCO
    real(RK),allocatable             :: projector(:,:)
    ! projectors, which must be added in order to obtain real coefficients
    real(RK),allocatable             :: intermediate(:,:)
    ! work array
    real(RK),allocatable             :: intermediate_irr(:,:)
    ! work array
    complex(CK),allocatable         :: intermediate_irr_c1(:,:),intermediate_irr_c2(:,:)
    ! work array
    ! complex efm
    real(RK),allocatable             :: eig_vec_r(:,:),eig_vec_i(:,:)
    ! real and imaginary part of the solution vectors
    real(RK),allocatable             :: eigen_value_c(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable             :: eigen_value(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable             :: eigen_valueII(:)
    ! eigenvalues of CSCO-II
    real(RK),allocatable             :: eigen_vector(:,:)
    ! eigenvectors of CSCO-II and CSCO
    real(RK),allocatable             :: gen_matrix(:,:,:)
    ! generator matrices for the irrep basis functions
    ! gen_matrix(5,dim_trafo,dim_trafo) 5 is the maximal dimension of all irreps of
    ! the pointgroups
    real(RK),allocatable                 :: work1(:), work2(:),work3(:,:)
    ! work arrays
    integer(IK)                          :: i,j,m,n,l,ierr,status
    ! counters
    integer(IK)                          :: irrep_dim,multiplicity,counter
    ! help variable for dimension of irrep
    real(RK)                             :: sq2
    ! 1/(2*sqrt(2))
    real(RK)                             :: factor
    ! factor for subgroups
    real(RK)                             :: gen_eigen_value
    ! eigenvalue of basisvector, from which basis is generated

    !------------ Declaration of local handles ----------------------
    type(efm_csco_eigenspace)  , pointer :: eigen_space_can_eigenspace
    type(efm_cscoII_eigenspace), pointer :: eigen_space_can_subspace
    !------------ Executable code ------------------------------------

    ! allocate work arrays
    allocate(efm_matrix(dim_trafo,dim_trafo),efmII_matrix(dim_trafo,dim_trafo), &
         stat=status)
    ASSERT(status==0)
    allocate(efm_matrix_c(dim_trafo,dim_trafo), stat=status)
    ASSERT(status==0)
    allocate(projector(dim_trafo,dim_trafo), stat=status)
    ASSERT(status==0)
    allocate(eigen_value(dim_trafo),eigen_vector(dim_trafo,dim_trafo), stat=status)
    ASSERT(status==0)
    allocate(eigen_valueII(dim_trafo), stat=status)
    ASSERT(status==0)
    allocate(gen_matrix(dim_trafo,dim_trafo,5), stat=status)
    ASSERT(status==0)
    allocate(work1(dim_trafo),work2(dim_trafo),intermediate(dim_trafo,dim_trafo), &
         stat=status)
    ASSERT(status==0)

    ! we determine vector irreps
    eigen_space_can%projective = .false.

    !
    !----- build (real) efm-matrix ------------------------------------------------------------
    !


    ! initialize efm-matrix
    efm_matrix = 0.0_rk
    efm_matrix_c = (0.0_rk,0.0_rk)

    ! prefactor
    sq2 = 1/(2.0_rk*sqrt(2.0_rk))

    ! first build efm for CSCO
    do i = 1,3
       if (ncsco(i).ne.0) then
          if (klass(ncsco(i))%ambivalent) then
             do j=1,klass(ncsco(i))%number
                efm_matrix = efm_matrix +&
                     trafo_matrix(:,:,klass(ncsco(i))%conjug_el(j))*fcsco(i)
             end do
          else
             do j=1,klass(ncsco(i))%number
                intermediate = trafo_matrix(:,:,klass(ncsco(i))%conjug_el(j))
                efm_matrix = efm_matrix +&
                     3*(intermediate + transpose(intermediate))*fcsco(i)*sq2
                efm_matrix_c = efm_matrix_c + &
                     &cmplx(0.0,1.0_rk,RK)*(intermediate - transpose(intermediate))*fcsco(i)*sq2
             end do
          endif
       endif
    end do

    ! construct complex CSCO-II matrix
    efm_matrix_c = cmplx(efm_matrix,0.0_rk,RK) + efm_matrix_c

    ! now build efm for CSCO-II

    efmII_matrix = efm_matrix
    factor = efm_group_factor/(efm_subgroup_factor**group_num_sg)
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups(m)%ncsco(i).ne.0) then
             if (subgroups(m)%klass(subgroups(m)%ncsco(i))%ambivalent) then
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   efmII_matrix = efmII_matrix +&
                        trafo_matrix(:,:,subgroups(m)%klass(subgroups(m)%ncsco(i))%&
                        &conjug_el(j))*subgroups(m)%fcsco(i)*factor
                end do
             else
                do j=1,subgroups(m)%klass(subgroups(m)%ncsco(i))%number
                   intermediate = trafo_matrix(:,:,subgroups(m)%klass(subgroups(m)&
                        &%ncsco(i))%conjug_el(j))
                   efmII_matrix = efmII_matrix + subgroups(m)%fcsco(i)*factor*sq2*&
                        &(3*(intermediate + transpose(intermediate)))
                   ! complex efm (disabled)
                   !cmplx(0.0,1,RK)*(intermediate - transpose(intermediate)))
                end do
             end if
          endif
       end do
       factor = factor*efm_subgroup_factor
    end do

    ! add projection operators to the csco-II matrix

    ! loop over all irreps
    do i=1,group_num_ir
       if (irrep_can(i)%build_projector) then
          projector = 0.0_rk
          do j=1,group_num_el
             projector = projector + trafo_matrix(:,:,j)*&
                  &irrep_can(i)%irrep_matrix(1,1,j)
          end do
          projector = projector*irrep_can(i)%dimension*0.314_rk&
               &/group_num_el
          efmII_matrix = efmII_matrix + projector
       end if
    end do


    ! solve eigenvalue problem for CSCO-II
    ierr = 0
    call rs(dim_trafo,dim_trafo,&                  ! dimensions of transformation
         efmII_matrix,&                            ! matrix of eigenvalue problem
         eigen_valueII,matz,eigen_vector,&         ! eigenvalues and eigenvectors
         work1,work2,ierr)                         ! work arrays
    if (ierr.ne.0) then
       write(*,*) "efm_cscoII_solve: stopped due to error in eigensolver IERR =",ierr
       call error_handler("efm_cscoII_solve: error in eigensolver")
    end if

!   if( group_name == "C1" )then
!      WARN('C1: overwriting eigenvectors')
!
!      eigen_vector = 0.0_rk
!      do i=1,dim_trafo
!         ! check that all eigenvalues are 1:
!         ASSERT(eigen_valueII(i)-1.0_rk<small)
!         eigen_valueII(i)  = 1.0_rk
!         ! force regular eigenvectors:
!         eigen_vector(i,i) = 1.0_rk
!      enddo
!      if( dim_trafo==3 )then
!         ! with this order ROTMAT for gradients is unity:
!         eigen_vector(:,:) = eigen_vector(:,(/1,3,2/))
!      endif
!   endif


    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(efm_matrix,eigen_vector)

    ! determine the multiple
    eigen_value = 0
    do j=1,dim_trafo
       do i=1,dim_trafo
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value(j) = -(intermediate(i,j)/eigen_vector(i,j))
             exit
          endif
       enddo
    enddo

    ! owing to some convention
    eigen_valueII = -eigen_valueII


    !
    ! ---- now generate the irreducible bases of the irreps ----------------
    !

    ! initialize canonically ordered basis of the irreps
    eigen_space_can%dimension = dim_trafo
    eigen_space_can%csco_eigenspaces(:)%exists = .false.
    eigen_space_can%csco_eigenspaces(:)%dim_irrep = 0

    ! count number of identified basis functions
    counter = 0

    ! loop over irreps
    do i=1,group_num_ir
       eigen_space_can_eigenspace => eigen_space_can%csco_eigenspaces(i)

       multiplicity = 0

       irrep_dim = irrep_can(i)%dimension
       gen_eigen_value = irrep_can(i)%gen_eigen_value

       ! in order to prevent division by zero
       if (abs(gen_eigen_value).lt.smaller) then
          gen_eigen_value = even_smaller
       endif

       ! determine multiplicity of the irrep
       do j=1,dim_trafo
          if ((abs((eigen_valueII(j)-gen_eigen_value)/gen_eigen_value).lt.small).or.&
               (abs(eigen_valueII(j)-gen_eigen_value).lt.smaller)) then
             multiplicity = multiplicity + 1
          end if
       end do
       counter = counter + multiplicity*irrep_can(i)%time_dimension

       if (irrep_can(i)%pseudo) then
          ! handle pseudo-2D irreps
          ! (a complex eigenvalue problem must be solved)

          ! allocate data structure of eigen_space_can for current irrep
          if (multiplicity.gt.0) then
             call efm_alloc(irrep_dim,eigen_space_can_eigenspace)
             do n=1,irrep_dim
                eigen_space_can_subspace => &
                     eigen_space_can_eigenspace%cscoII_subspaces(n)
                call efm_alloc(dim_trafo,multiplicity/2,eigen_space_can_subspace)
             end do
          else
             cycle
          end if

          ! allocate work array
          allocate(intermediate_irr(dim_trafo,multiplicity))

          ! fill basis of first eigen value of CSCO-II
          multiplicity = 0
          do j=1,dim_trafo
             if ((abs((eigen_valueII(j)-gen_eigen_value)/gen_eigen_value).lt.small).or.&
                  (abs(eigen_valueII(j)-gen_eigen_value).lt.smaller))then
                multiplicity = multiplicity + 1
                intermediate_irr(:,multiplicity) = eigen_vector(:,j)
             end if
          end do

          ! allocate auxiliary arrays
          allocate(intermediate_irr_c1(multiplicity,multiplicity), stat=status)
          ASSERT(status==0)
          allocate(intermediate_irr_c2(dim_trafo,multiplicity), stat=status)
          ASSERT(status==0)
          allocate(eig_vec_r(multiplicity,multiplicity),eig_vec_i(multiplicity,multiplicity), &
               stat=status)
          ASSERT(status==0)
          allocate(eigen_value_c(multiplicity), stat=status)
          ASSERT(status==0)
          allocate(work3(multiplicity,multiplicity), stat=status)
          ASSERT(status==0)

          ! determine CSCO-II in eigenspace of the irrep
          intermediate_irr_c2 = matmul(efm_matrix_c,intermediate_irr)
          intermediate_irr_c1 = -matmul(transpose(intermediate_irr),intermediate_irr_c2)

          ! complex eigenvalue problem
          call ch(multiplicity,multiplicity,&                           ! dimensions of transformation
               real(intermediate_irr_c1),&                              ! matrix of eigenvalue problem
               aimag(intermediate_irr_c1),&                             ! matrix of eigenvalue problem
               eigen_value_c,matz,eig_vec_r,eig_vec_i,&                 ! eigenvalues and eigenvectors
               work1(1:multiplicity),work2(1:multiplicity),work3,ierr)  ! work arrays
          if (ierr.ne.0) then
             write(*,*) "efm_cscoII_solve: stopped due to error in eigensolver IERR =",ierr
             call error_handler("efm_cscoII_solve: error in eigensolver")
          end if


          ! determine basis
          l=1
          do j=1,multiplicity
             if (abs(eigen_value_c(j)-irrep_can(i)%eigen_value).lt.small) then
                eigen_space_can_subspace => &
                     eigen_space_can_eigenspace%cscoII_subspaces(1)
                factor = dot_product(eig_vec_r(:,j),eig_vec_r(:,j))
                factor = sqrt(1/factor)
                eigen_space_can_subspace%eigenfunction(:,l) = &
                     factor*MATMUL(intermediate_irr,eig_vec_r(:,j))
                eigen_space_can_subspace => &
                     eigen_space_can_eigenspace%cscoII_subspaces(2)
                factor = dot_product(eig_vec_i(:,j),eig_vec_i(:,j))
                factor = sqrt(1/factor)
                eigen_space_can_subspace%eigenfunction(:,l) = &
                     factor*MATMUL(intermediate_irr,eig_vec_i(:,j))
                l = l+1
             end if
          end do

          deallocate(intermediate_irr_c1, stat=status)
          ASSERT(status==0)
          deallocate(intermediate_irr_c2, stat=status)
          ASSERT(status==0)
          deallocate(eig_vec_r,eig_vec_i, stat=status)
          ASSERT(status==0)
          deallocate(eigen_value_c, stat=status)
          ASSERT(status==0)
          deallocate(work3, stat=status)
          ASSERT(status==0)

       else

          ! allocate data structure of eigen_space_can for current irrep
          if (multiplicity.gt.0) then
             call efm_alloc(irrep_dim,eigen_space_can_eigenspace)
             do n=1,irrep_dim
                eigen_space_can_subspace => &
                     eigen_space_can_eigenspace%cscoII_subspaces(n)
                call efm_alloc(dim_trafo,multiplicity,eigen_space_can_subspace)
             end do
          else
             cycle
          end if

          ! allocate work array
          allocate(intermediate_irr(dim_trafo,multiplicity))

          ! fill basis of first eigen value of CSCO-II
          multiplicity = 0
          do j=1,dim_trafo
             if ((abs((eigen_valueII(j)-gen_eigen_value)/gen_eigen_value).lt.small).or.&
                  (abs(eigen_valueII(j)-gen_eigen_value).lt.smaller)) then
                multiplicity = multiplicity + 1
                eigen_space_can_subspace => &
                     eigen_space_can_eigenspace%cscoII_subspaces(1)
                eigen_space_can_subspace%eigenfunction(:,multiplicity) = &
                     eigen_vector(:,j)
             end if
          end do

          !
          ! now generate the other basis vectors from the first one
          !

          ! in order to do so the generator matrices are constructed
          do j=1,irrep_dim
             gen_matrix(:,:,j) = 0.0_rk
             do n=1,irrep_dim
                gen_matrix(:,:,j) = gen_matrix(:,:,j) + trafo_matrix(:,:,irrep_can(i)%&
                     irrep_gen_oper(n))*irrep_can(i)%irrep_gen_coeff(n,j)
             end do
          end do

          ! now generate all basis vectors
          eigen_space_can_subspace => &
               eigen_space_can_eigenspace%cscoII_subspaces(1)
          intermediate_irr = eigen_space_can_subspace%eigenfunction
          do j=2,irrep_dim
             eigen_space_can_subspace => &
                  eigen_space_can_eigenspace%cscoII_subspaces(j)
             eigen_space_can_subspace%eigenfunction = &
                  matmul(gen_matrix(:,:,j),intermediate_irr)
          end do

       end if

       deallocate(intermediate_irr, stat=status)
       ASSERT(status==0)
    end do

    ! check if all basis functions were identified
    if (counter.ne.dim_trafo) then
       if (output_unit > 0) then
          write(output_unit,*) "efm_cscoII_solve: number of identified irreps/rows ",counter,&
               " is lower than dimension ",dim_trafo
          write(output_unit,*) "efm_cscoII_solve: Please contact Markus Mayer mayer@theochem.tu-muenchen.de"
       endif
       stop
    endif

    ! deallocate work arrays
    deallocate(efm_matrix,efmII_matrix, stat=status)
    ASSERT(status==0)
    deallocate(efm_matrix_c, stat=status)
    ASSERT(status==0)
    deallocate(projector, stat=status)
    ASSERT(status==0)
    deallocate(eigen_value,eigen_vector, stat=status)
    ASSERT(status==0)
    deallocate(eigen_valueII, stat=status)
    ASSERT(status==0)
    deallocate(gen_matrix, stat=status)
    ASSERT(status==0)
    deallocate(intermediate, stat=status)
    ASSERT(status==0)
    deallocate(work1,work2, stat=status)
    ASSERT(status==0)

  end subroutine efm_cscoII_solve
  !*************************************************************

  !*************************************************************
  subroutine efm_proj_cscoII_solve(trafo_matrix,dim_trafo,eigen_space_can)
    !  Purpose: solves eigenequation for a special
    !           CSCO-II for projective representations
    !------------ Modules used ------------------- ---------------
    use type_module
    use iounitadmin_module
    use group_module
    implicit none
    !------------ Declaration of formal parameters  -----------------
    integer(IK), intent(in)           :: dim_trafo
    ! dimension of the transformation matrix
    complex(CK),    intent(in)       :: trafo_matrix(:,:,:)
    ! transformation matrices of the invariant eigenspace
    type(efm_cscoII_diag_invspace), intent(out) :: eigen_space_can
    ! canonically ordered eigenvectors of cscoII
    !** End of interface *****************************************
    !------------ Declaration of local parameters ---------------
    integer(IK),parameter          :: matz = 1
    ! parameter for eigensolver
    real(RK),parameter             :: small = 0.0000001_rk
    ! very small value
    real(RK),parameter             :: smaller = 1.0e-12_rk
    ! even smaller value

    !------------ Declaration of local variables  -------------------
    complex(CK),allocatable         :: efm_matrix(:,:)
    ! complex matrix of the eigenvalue problem for the CSCO
    complex(CK),allocatable         :: efmII_matrix(:,:)
    ! complex matrix of the eigenvalue problem for the CSCO-II
    complex(CK),allocatable         :: csco_sub(:,:)
    ! complex matrix of the CSCO of subgroups
    complex(CK),allocatable         :: intermediate(:,:)
    ! work array
    complex(CK),allocatable          :: intermediate_irr(:,:)
    ! work array
    ! complex efm
    real(RK),allocatable             :: eig_vec_r(:,:),eig_vec_i(:,:)
    ! real and imaginary part of the solution vectors
    real(RK),allocatable             :: eigen_value(:),eigen_valueII(:),eigen_value_sub(:)
    ! eigenvalues of CSCO and CSCO-II
    complex(CK),allocatable         :: eigen_vector(:,:)
    ! eigenvectors of CSCO-II and CSCO
    complex(CK),allocatable         :: gen_matrix(:,:,:)
    ! generator matrices for the irrep basis functions
    ! gen_matrix(5,dim_trafo,dim_trafo) 5 is the maximal dimension of all irreps of
    ! the pointgroups
    real(RK),allocatable                 :: work1(:), work2(:),work3(:,:)
    ! work arrays
    integer(IK)                          :: i,j,m,n,ierr
    ! counters
    integer(IK)                          :: irrep_dim,multiplicity,counter
    ! help variable for dimension of irrep
    real(RK)                             :: sq2
    ! 1/(2*sqrt(2))
    real(RK)                             :: factor
    ! factor for subgroups
    real(RK)                             :: gen_eigen_value
    ! eigenvalue of basisvector, from which basis is generated

    !------------ Executable code ------------------------------------

    ! allocate work arrays
    allocate(efm_matrix(dim_trafo,dim_trafo),efmII_matrix(dim_trafo,dim_trafo),csco_sub(dim_trafo,dim_trafo))
    allocate(eig_vec_r(dim_trafo,dim_trafo),eig_vec_i(dim_trafo,dim_trafo))
    allocate(eigen_value(dim_trafo),eigen_vector(dim_trafo,dim_trafo))
    allocate(eigen_valueII(dim_trafo),eigen_value_sub(dim_trafo))
    allocate(gen_matrix(dim_trafo,dim_trafo,6))
    allocate(work1(dim_trafo),work2(dim_trafo),work3(2,dim_trafo),&
         intermediate(dim_trafo,dim_trafo))

    !
    !----- build efm-matrix ------------------------------------------------------------
    !

    ! we determine projective irreps
    eigen_space_can%projective = .true.

    !am: print*," Entered CSCO-II Solver"
    ! initialize efm-matrix
    efm_matrix = 0.0_rk

    ! prefactor
    sq2 = 1/(2.0_rk*sqrt(2.0_rk))

    ! first build efm for CSCO
    do i = 1,3
       if (ncsco_proj(i).ne.0) then
          if (klass(ncsco_proj(i))%ambivalent_proj) then
             do j=1,klass(ncsco_proj(i))%number
                efm_matrix = efm_matrix +&
                     trafo_matrix(:,:,klass(ncsco_proj(i))%conjug_el(j))*fcsco_proj(i)
             end do
          else
             do j=1,klass(ncsco_proj(i))%number
                intermediate = trafo_matrix(:,:,klass(ncsco_proj(i))%conjug_el(j))
                efm_matrix= efm_matrix + &
                     (3.0_rk*(intermediate + transpose(conjg(intermediate)))+&
                     &cmplx(0.0,1.0_rk,RK)*(intermediate - transpose(conjg(intermediate))))*fcsco_proj(i)*sq2
             end do
          endif
       endif
    end do

    ! now build efm for CSCO-II

    ! for debugging only >>
    !efm_matrix = 0.0_rk
    ! for debugging only <<
    csco_sub = 0.0_rk
    !efmII_matrix = efm_matrix
    factor = 1.0_rk
    ! loop over subgroups
    do m = 1,group_num_sg
       do i = 1,3
          if (subgroups_proj(m)%ncsco_proj(i).ne.0) then
             if (subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%ambivalent_proj) then
                do j=1,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%number
                   csco_sub = csco_sub +&
                        trafo_matrix(:,:,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%&
                        &conjug_el(j))*subgroups_proj(m)%fcsco_proj(i)*factor
                end do
             else
                do j=1,subgroups_proj(m)%klass(subgroups_proj(m)%ncsco_proj(i))%number
                   intermediate = trafo_matrix(:,:,subgroups_proj(m)%klass(subgroups_proj(m)&
                        &%ncsco_proj(i))%conjug_el(j))
                   csco_sub = csco_sub + subgroups_proj(m)%fcsco_proj(i)*factor*sq2*&
                        &(3.0_rk*(intermediate + transpose(conjg(intermediate)))+&
                        &10.0_rk*cmplx(0.0_rk,1.0_rk,RK)*(intermediate - transpose(conjg(intermediate))))
                end do
             endif
          endif
       end do
       factor = factor*efm_proj_subgroup_factor
    end do

    ! build CSCO-II matrix
    !am: print*,"group_num_sg",group_num_sg
    efmII_matrix = efm_matrix + csco_sub*efm_proj_group_factor/(efm_proj_subgroup_factor**group_num_sg)

    !am: print*,"solve eigenvalue problem for CSCO-II"
    ! solve eigenvalue problem for CSCO-II
    ierr = 0
    call ch(dim_trafo,dim_trafo,&                     ! dimensions of transformation
         real(efmII_matrix),aimag(efmII_matrix),&     ! matrix of eigenvalue problem
         eigen_valueII,matz,eig_vec_r,eig_vec_i,&     ! eigenvalues and eigenvectors
         work1,work2,work3,ierr)                      ! work arrays
    if (ierr.ne.0) then
       write(*,*) "efm_proj_cscoII_solve: stopped due to error in eigensolver IERR =",ierr
       call error_handler("efm_proj_cscoII_solve: error in eigensolver")
    end if

    !am: print*," construct complex eigenvectors"
    ! construct complex eigenvectors
    eigen_vector = cmplx(eig_vec_r,eig_vec_i,RK)

    ! owing to some convention
    eigen_valueII = -eigen_valueII

    ! multiply CSCO efm-matrix by the eigenvectors
    ! (in order to obtain irrep labels)
    intermediate = matmul(efm_matrix,eigen_vector)

    ! determine the multiple
    eigen_value = 0
    do j=1,dim_trafo
       do i=1,dim_trafo
          if (abs(intermediate(i,j)).gt.small) then
             eigen_value(j) = -real(intermediate(i,j)/eigen_vector(i,j))
             exit
          endif
       enddo
    enddo

    ! multiply CSCO efm-matrix by the eigenvectors
    intermediate = matmul(csco_sub,eigen_vector)

    ! determine the multiple
    eigen_value_sub = 0.0_rk
    do j=1,dim_trafo
       do i=1,dim_trafo
          if ((abs(intermediate(i,j)).gt.small).and.(abs(eigen_vector(i,j)).gt.small)) then
             eigen_value_sub(j) = -real(intermediate(i,j)/eigen_vector(i,j))
             exit
          endif
       enddo
    enddo

    ! now deallocate work arrays
    deallocate(efmII_matrix)
    deallocate(efm_matrix)
    deallocate(csco_sub)
    deallocate(eig_vec_r,eig_vec_i)
    deallocate(work1,work2,work3,intermediate)

    !
    ! ---- now generate the irreducible bases of the irreps ----------------
    !

    ! initialize canonically ordered basis of the irreps
    eigen_space_can%dimension = dim_trafo
    eigen_space_can%csco_eigenspaces%dim_irrep = 0 !am: Will it be ever set ???
    eigen_space_can%csco_eigenspaces%exists = .false.

    ! count number of identified basis functions
    counter = 0
    ! loop over irreps
    do i=1,group_num_re
       multiplicity = 0

       irrep_dim = proj_irrep_can(i)%dimension
       gen_eigen_value = proj_irrep_can(i)%gen_eigen_value

       ! determine multiplicity of the irrep
       do j=1,dim_trafo
          if (abs((eigen_valueII(j)-gen_eigen_value)/gen_eigen_value).lt.small) then
             multiplicity = multiplicity + 1
          end if
       end do
       counter = counter + multiplicity*proj_irrep_can(i)%dimension

       ! allocate data structure of eigen_space_can for current irrep
       if (multiplicity.gt.0) then

!!$          eigen_space_can%csco_eigenspaces(i)%exists = .true.
!!$          !mdf>>>
!!$          eigen_space_can%csco_eigenspaces(i)%dim_irrep = irrep_dim
!!$          !<<<mdf
!!$          allocate(eigen_space_can%csco_eigenspaces(i)%cscoII_subspaces(irrep_dim))

          call efm_alloc(irrep_dim,eigen_space_can%csco_eigenspaces(i),EXISTS=.true.)

          do n=1,irrep_dim
!!$             eigen_space_can%csco_eigenspaces(i)%cscoII_subspaces(n)%multiplicity =&
!!$                  multiplicity
!!$             allocate(eigen_space_can%csco_eigenspaces(i)%&
!!$                  &cscoII_subspaces(n)%eigenfunction_proj(dim_trafo,multiplicity))
             call efm_alloc(&
                  & dim_trafo,multiplicity,&
                  & eigen_space_can%csco_eigenspaces(i)%cscoII_subspaces(n),&
                  & PROJ=.true.)
          end do
       else
          cycle
       end if

       ! allocate work array
       allocate(intermediate_irr(dim_trafo,multiplicity))

       ! fill basis of first eigen value of CSCO-II
       multiplicity = 0
       do j=1,dim_trafo
          if (abs((eigen_valueII(j)-gen_eigen_value)/gen_eigen_value).lt.small) then
             multiplicity = multiplicity + 1
             eigen_space_can%csco_eigenspaces(i)%cscoII_subspaces(1)%&
                  eigenfunction_proj(:,multiplicity) = eigen_vector(:,j)
          end if
       end do

       !
       ! now generate the other basis vectors from the first one
       !

       ! in order to do so the generator matrices are constructed
       do j=1,irrep_dim
          gen_matrix(:,:,j) = 0.0_rk
          do n=1,irrep_dim
             gen_matrix(:,:,j) = gen_matrix(:,:,j) + trafo_matrix(:,:,proj_irrep_can(i)%&
                  irrep_gen_oper(n))*proj_irrep_can(i)%irrep_gen_coeff(n,j)
          end do
       end do

       ! now generate all basis vectors
       intermediate_irr = eigen_space_can%csco_eigenspaces(i)%cscoII_subspaces(1)%&
            eigenfunction_proj
       do j=2,irrep_dim
          eigen_space_can%csco_eigenspaces(i)%cscoII_subspaces(j)%&
               eigenfunction_proj = matmul(gen_matrix(:,:,j),intermediate_irr)
       end do

       deallocate(intermediate_irr)

    end do


    ! check if all basis functions were identified
    if (counter.ne.dim_trafo) then
       write(output_unit,*) "efm_proj_cscoII_solve: number of identified irreps/rows ",counter,&
            " is lower than dimension ",dim_trafo
       write(output_unit,*) "efm_proj_cscoII_solve: Please contact Markus Mayer mayer@theochem.tu-muenchen.de"
       stop
    endif


    !mdf>>>
    call efm_rephase(eigen_space_can)

    deallocate(eigen_value,eigen_valueII,eigen_value_sub)
    deallocate(gen_matrix)
    deallocate(eigen_vector)
  end subroutine efm_proj_cscoII_solve
  !*************************************************************

#if FPP_COMPILE_UNUSED_CODE
  !*************************************************************
  subroutine efm_check_proj_irrep(eigen_space_can,trafo_matrix,dim_trafo)
    ! checks if eigen_space_can is actually
    ! ordered by irreducible invariant subspaces
    ! of the group
    use iounitadmin_module
    use group_module
    implicit none
    !------------ Declaration of formal parameters  -----------------
    type(efm_cscoII_diag_invspace)         :: eigen_space_can
    ! canonically ordered eigenvectors of cscoII
    complex(CK),intent(in)      :: trafo_matrix(:,:,:)
    ! transformation matrix trafo_matrix(:,:,symmetry_element)
    integer(IK),intent(in)       :: dim_trafo
    ! dimension of transformation
    !** End of interface *****************************************
    !------------ Declaration of local constants -----------------
    real(RK),parameter           :: small = 0.000001_rk
    !------------ Declaration of local variables ---------------------
    integer(IK)                  :: i,j,i_irrep,multiplicity,irrep_dim,n
    ! counters
    complex(CK),allocatable     :: product_matrix(:,:),diff_matrix(:,:),irrep_basis(:,:)
    ! product of two rep matrices
    !------------ Executable code ------------------------------------
    !am: print*,"entered efm_check_proj_irrep"
    ! loop over irreps
    outer:do i_irrep=1,group_num_re
       !am: print*," check irrep ",proj_irrep_can(i_irrep)%label
       if (eigen_space_can%csco_eigenspaces(i_irrep)%exists) then
          !am: print*," irrep exists"
          multiplicity = eigen_space_can%csco_eigenspaces(i_irrep)%&
               &cscoII_subspaces(1)%multiplicity
          irrep_dim = proj_irrep_can(i_irrep)%dimension
          allocate(irrep_basis(dim_trafo,irrep_dim))
          allocate(product_matrix(irrep_dim,irrep_dim))
          allocate(diff_matrix(irrep_dim,irrep_dim))
          do i=1,multiplicity
             !am: print*," test multiplicity ",i
             do j=1,irrep_dim
                irrep_basis(:,j) = eigen_space_can%csco_eigenspaces(i_irrep)%&
                     &cscoII_subspaces(j)%eigenfunction_proj(:,i)
             end do
             !am: print*," determined basis of irrep"
             do n=1,group_num_el
                !am: print*, " test trafo ",n
                product_matrix = matmul(transpose(conjg(irrep_basis)),matmul(trafo_matrix(:,:,n),irrep_basis))
                !am: print*," determined product"
                diff_matrix = proj_irrep_can(i_irrep)%irrep_matrix(:,:,n) - product_matrix
                if (maxval(abs(diff_matrix)).gt.small) then
                   write(output_unit,*) "CAUTION: eigenspace is not an irrep"
                   stop
                   exit outer
                end if
             end do
             !am: print*," tested basis"
          end do
          deallocate(irrep_basis,product_matrix,diff_matrix)
       end if
    end do outer
    write(output_unit,*) "    Irreps are alright"
    !am: print*,"left efm_check_proj_irrep"

  end subroutine efm_check_proj_irrep
  !*************************************************************
#endif

#if FPP_COMPILE_UNUSED_CODE
  !*************************************************************
  subroutine efm_check_rep(trafo_matrix,dim_trafo)
    ! checks if transformation
    ! is a projective representation of the group
    use iounitadmin_module
    use group_module
    implicit none
    !------------ Declaration of formal parameters  -----------------
    complex(CK),intent(in)      :: trafo_matrix(:,:,:)
    ! transformation matrix trafo_matrix(:,:,symmetry_element)
    integer(IK),intent(in)       :: dim_trafo
    ! dimension of transformation
    !** End of interface *****************************************
    !------------ Declaration of local constants -----------------
    real(RK),parameter           :: small = 0.000001_rk
    !------------ Declaration of local variables ---------------------
    integer(IK)                  :: i,j
    ! counters
    complex(CK),allocatable     :: product_matrix(:,:),diff_matrix(:,:)
    ! product of two rep matrices
    !------------ Executable code ------------------------------------

    allocate(product_matrix(dim_trafo,dim_trafo))
    allocate(diff_matrix(dim_trafo,dim_trafo))

    check:do i=1,group_num_el
       do j=1,group_num_el
          product_matrix = matmul(trafo_matrix(:,:,i),trafo_matrix(:,:,j))
          diff_matrix = product_matrix - group_factab(i,j)*trafo_matrix(:,:,group_multab(i,j))
          if (maxval(abs(diff_matrix)).gt.small) then
             write(output_unit,*) "product (",i,",",j,") is wrong "
             write(output_unit,*) "CAUTION: NO PROJECTIVE REPRESENTATION"
             stop
             exit check
          else
             !write(output_unit,*) "product (",i,",",j,") is correct"
          endif
       end do
    end do check

    deallocate(product_matrix,diff_matrix)
  end subroutine efm_check_rep
#endif

#if FPP_COMPILE_UNUSED_CODE
  subroutine efm_init_c1_space(efm_type_inv,multi)
    ! Purpose : set invariant eigenspace to original basis
    ! Subroutine called by : symm_symmadapt
    !
    !** End of interface *****************************************
    type(efm_cscoII_diag_invspace)      :: efm_type_inv
    integer(IK)               :: multi
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of local variables ---------------------
    integer(IK)   :: error,i
    external error_handler
    efm_type_inv%dimension = multi
    efm_type_inv%projective = .false.
    allocate(efm_type_inv%csco_eigenspaces(1)%cscoII_subspaces(1),&
         & stat=error)
    if (error/=0) call error_handler &
         ("efm_init_c1_space : allocation (1) failed")
    efm_type_inv%csco_eigenspaces(1)%dim_irrep = 1
    efm_type_inv%csco_eigenspaces(1)%exists = .true.
    allocate(efm_type_inv%csco_eigenspaces(1)%cscoII_subspaces(1)%&
         & eigenfunction(multi,multi),stat=error)
    if (error/=0) call error_handler &
         ("efm_init_c1_space : allocation (1) failed")
    efm_type_inv%csco_eigenspaces(1)%cscoII_subspaces(1)%multiplicity = multi
    efm_type_inv%csco_eigenspaces(1)%cscoII_subspaces(1)%eigenfunction(:,:)= 0.0_rk
    do i=1,multi
       efm_type_inv%csco_eigenspaces(1)%cscoII_subspaces(1)%&
            & eigenfunction(i,i)= 1.0_rk
    enddo

  end subroutine efm_init_c1_space
#endif
  !*************************************************************

  !*************************************************************
  !*************************************************************

end module efm_module
