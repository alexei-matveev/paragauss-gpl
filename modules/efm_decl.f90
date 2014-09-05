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
module  efm_decl
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
!== Interrupt of public interface of module =====================
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

  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind,&
       & CK=>c16_kind ! type specification parameters
  implicit none
  save            ! keep contents of module
  private         ! by default, all names are private

       
!== Interrupt end of public interface of module =================

  ! ------------ Declaration of types ------------------------------

  ! The following types describe the structure of an invariant eigenspace
  ! of the group, for example a SO(3) eigenspace for some angular momentum l,
  ! in terms of eigenfunctions of the CSCO-II.

  type,  public ::  efm_cscoII_eigenspace
     ! eigenspace of csco-II 
     real(RK)                    :: eigenvalue
     ! eigenvalue of csco-II
     integer(IK)                 :: multiplicity
     ! number of eigenfunctions of the eigenspace
     real(RK),pointer            :: eigenfunction(:,:)
     ! eigenfunction of eigenspace, eigenfunction(dim_trafo,multiplicity)
     complex(CK),pointer         :: eigenfunction_proj(:,:)
     ! eigenfunction of eigenspace, eigenfunction(dim_trafo,multiplicity)
     !
     logical                     :: projective ! added
  end type efm_cscoII_eigenspace

  type, public :: efm_csco_eigenspace
     ! eigenspace of csco
     ! contains all representation functions for a special irrep
     logical                     :: exists
     ! does eigenspace for irrep exist ?
     real(RK)                    :: eigenvalue
     ! eigenvalue of csco
     integer(IK)                 :: dim_irrep
     ! dimension of the irrep = number of csco-II subspaces
     type(efm_cscoII_eigenspace), pointer  :: cscoII_subspaces(:)
     ! subspaces of eigenspace of csco, which are eigenspaces of csco-II
  end type efm_csco_eigenspace

  type, public :: efm_cscoII_diag_invspace
     ! an invariant space of the group diagonalized to
     ! the csco-II of the group
     integer(IK)                           :: dimension,n_irr
     ! dimension of invariant space
     logical                               :: projective
     ! is representation projective?
     type(efm_csco_eigenspace),pointer     :: csco_eigenspaces(:)
     ! eigenspaces of the csco = eigenspaces of the irreps
     ! csco_eigenspaces(group_num_ir) 
  end type efm_cscoII_diag_invspace

  !------------ Interface statements ---------------------------

  interface  efm_direct_product
     module procedure direct_product_rr
     module procedure direct_product_rc
     module procedure direct_product_cr
     module procedure direct_product_cc
  end interface

  interface merge_indices
     module procedure merge_indices_real2D
     module procedure merge_indices_real4D
  end interface

  interface restore_indices
     module procedure restore_indices_real2D
     module procedure restore_indices_cmplx2D
  end interface

  interface efm_alloc  ! public alloc
     module procedure alloc_invspace
     module procedure alloc_csco_eigenspace
     module procedure alloc_cscoII_eigenspace
     module procedure alloc_arr_cscoII_eigenspace
  end interface

  interface alloc      ! private alloc
     module procedure alloc_invspace
     module procedure alloc_csco_eigenspace
     module procedure alloc_cscoII_eigenspace
     module procedure alloc_arr_cscoII_eigenspace
  end interface

  interface dealloc     ! private dealloc
     module procedure dealloc_invspace
     module procedure dealloc_csco_eigenspace
     module procedure dealloc_cscoII_eigenspace
  end interface

  interface efm_free    ! public recursive dealloc
     module procedure free_invspace
     module procedure free_csco_eigenspace
  end interface

  interface free ! private recursive dealloc
     module procedure free_invspace
     module procedure free_csco_eigenspace
     module procedure dealloc_cscoII_eigenspace ! no further recursion required
  end interface

  interface extract_basis
     module procedure extract_proj_basis
     module procedure extract_vec_basis
  end interface

  interface load_basis
     module procedure load_proj_basis
     module procedure load_vec_basis 
  end interface

  interface efm_gauge
     module procedure efm_inv_cgauge_io
     module procedure efm_inv_rgauge_io
  end interface
  
  !------------ public functions and subroutines ---------------------
  public :: &
       & efm_direct_product,&
       & efm_alloc,&
       & efm_free,&
       & efm_show_invspace,&
       & efm_gauge,&
       & efm_rephase,&
       & extract_basis,  &
       & load_basis,     &
       & phase,          &
       & multiplicity,   &
       & multiplicities, &
       & dim_vector,&
       & dim_rep,&
       & restore_indices,&
       & efm_find_one


  !===================================================================
  ! End of public interface of module
  !===================================================================


contains


  function efm_find_one(sa) result(ix)
    use group_module, only: group_num_re, group_num_ir
    implicit none
    type(efm_cscoII_diag_invspace),intent(in) :: sa
    integer(IK) :: ix !<<<result
    ! *** end of interface ***
    
    integer(IK) :: i,n

    if(count(sa%csco_eigenspaces(:)%exists) /= 1)&
         & call error_handler('efmd/efm_find_one: not ONE')

    if(sa%projective)then
       n = group_num_re
    else ! i.e. vector
       n = group_num_ir
    endif

    if(count(sa%csco_eigenspaces(:)%exists) /= 1)&
         & call error_handler('efmd/efm_find_one: not ONE')

    ix = 0
    do i=1,n
       if(sa%csco_eigenspaces(i)%exists) ix=i
    enddo
    if(ix==0)call error_handler('efmd/efm_find_one: ERROR, irrep not found')
  end function efm_find_one

  subroutine alloc_invspace(n_irr,inv,proj,dim)
    implicit none
    integer(IK)                                  :: n_irr
    type(efm_cscoII_diag_invspace),intent(inout) :: inv
    logical,optional                             :: proj
    integer(IK),optional                         :: dim
    ! *** end of interface ***

    integer(IK)     :: memstat
    logical         :: projective
    integer(IK)     :: dimension

    projective=.false.
    if(present(proj))projective=proj

    dimension = 0
    if(present(dim))dimension = dim

    inv%n_irr      = n_irr
    inv%projective = projective
    inv%dimension  = dimension

    allocate(inv%csco_eigenspaces(n_irr),STAT=memstat)
    if(memstat/=0)call error_handler("efmd/alloc_invspace: alloc failed")
  end subroutine alloc_invspace

  subroutine dealloc_invspace(inv)
    implicit none
    type(efm_cscoII_diag_invspace),intent(inout) :: inv
    ! *** end of interface ***

    integer(IK)     :: memstat

    deallocate(inv%csco_eigenspaces,STAT=memstat)
    if(memstat/=0)call error_handler("efmd/dealloc_invspace: dealloc failed")
    inv%projective = .false.
    inv%dimension  = -1
    inv%n_irr      = -1
  end subroutine dealloc_invspace

  subroutine free_invspace(inv)
    implicit none
    type(efm_cscoII_diag_invspace),intent(inout) :: inv
    ! *** end of interface ***

    integer(IK) :: irr

    do irr=1,size(inv%csco_eigenspaces)
       if(inv%csco_eigenspaces(irr)%exists)then
          call free(inv%csco_eigenspaces(irr))
       endif
    enddo
    call dealloc(inv)
  end subroutine free_invspace

  subroutine alloc_csco_eigenspace(dim_irrep,eigs,exists)
    use error_module, only: error
    implicit none
    integer(IK),intent(in)                  :: dim_irrep
    type(efm_csco_eigenspace),intent(inout) :: eigs
    logical,optional,intent(in)             :: exists
    ! *** end of interface ***

    integer :: memstat
    logical :: exsists_

    exsists_ = .true.
    if(present(exists))exsists_ = exists

    if( dim_irrep>0 )then
       eigs%exists    = exsists_
       eigs%dim_irrep = dim_irrep
       allocate(eigs%cscoII_subspaces(dim_irrep),STAT=memstat)
       if(memstat/=0)call error("efmd/alloc_csco_eigenspace: alloc failed")
    else
       eigs%exists    = .false.
       eigs%dim_irrep = -1
       nullify(eigs%cscoII_subspaces)
    endif
    eigs%eigenvalue = 0.0_rk
  end subroutine alloc_csco_eigenspace

  subroutine dealloc_csco_eigenspace(eigs)
    use error_module, only: error
    implicit none
    type(efm_csco_eigenspace),intent(inout) :: eigs
    ! *** end of interface ***

    integer(IK) :: memstat

    if(eigs%exists) then
       deallocate(eigs%cscoII_subspaces,STAT=memstat)
       if(memstat/=0)call error("efmd/dealloc_csco_eigenspace: dealloc failed")
    else
       call error("efmd/dealloc_csco_eigenspace: doesnt exist")
       return
    endif
    eigs%dim_irrep  = -1
    eigs%exists     = .false.
    eigs%eigenvalue = -1.0_rk
  end subroutine dealloc_csco_eigenspace

  subroutine free_csco_eigenspace(eigs)
    implicit none
    type(efm_csco_eigenspace),intent(inout) :: eigs
    ! *** end of interface ***

    integer(IK) :: p

    do p=1,size(eigs%cscoII_subspaces)
       call free(eigs%cscoII_subspaces(p))
    enddo
    call dealloc(eigs)
  end subroutine free_csco_eigenspace

  subroutine alloc_cscoII_eigenspace(dim_vec,n_indep,eigs,proj)
    implicit none
    integer(IK),intent(in)                    :: dim_vec,n_indep
    type(efm_cscoII_eigenspace),intent(inout) :: eigs
    logical,optional,intent(in)               :: proj
    ! *** end of interface ***

    integer(IK) :: memstat
    logical     :: projective

    projective = .false.
    if(present(proj))projective = proj

    eigs%projective   = projective
    eigs%multiplicity = n_indep
    eigs%eigenvalue   = 0.0_rk

    if(projective)then
       allocate(eigs%eigenfunction_proj(dim_vec,n_indep),STAT=memstat)
       if(memstat/=0)&
            & call error_handler("efmd/alloc_cscoII_eigenspace: alloc PROJ failed")
    else
       allocate(eigs%eigenfunction(dim_vec,n_indep),STAT=memstat)
       if(memstat/=0)&
            & call error_handler("efmd/alloc_cscoII_eigenspace: alloc VEC  failed")
    endif
  end subroutine alloc_cscoII_eigenspace

  subroutine dealloc_cscoII_eigenspace(eigs)
    implicit none
    type(efm_cscoII_eigenspace),intent(inout) :: eigs
    ! *** end of interface ***

    integer(IK) :: memstat

    if(eigs%projective)then
       deallocate(eigs%eigenfunction_proj,STAT=memstat)
       if(memstat/=0)&
            & call error_handler("efmd/dealloc_cscoII_eigenspace: dealloc PROJ failed")
    else
       deallocate(eigs%eigenfunction,STAT=memstat)
       if(memstat/=0)&
            & call error_handler("efmd/dealloc_cscoII_eigenspace: dealloc VEC  failed")
    endif
    eigs%projective   = .false.
    eigs%multiplicity = -1
    eigs%eigenvalue   = -1.0_rk
  end subroutine dealloc_cscoII_eigenspace

  subroutine alloc_arr_cscoII_eigenspace(dim_vec,n_indep,eigs,proj)
    implicit none
    integer(IK),intent(in)                    :: dim_vec,n_indep
    type(efm_cscoII_eigenspace),intent(inout) :: eigs(:)
    logical,optional,intent(in)               :: proj
    ! *** end of interface ***

    integer(IK) :: p,dim_irrep

    dim_irrep = size(eigs)

    do p=1,dim_irrep
       if(present(proj))then
          call alloc(dim_vec,n_indep,eigs(p),proj=proj)
       else
          call alloc(dim_vec,n_indep,eigs(p))
       endif
    enddo
  end subroutine alloc_arr_cscoII_eigenspace

  !--------------------------------------------
  !
  ! Direct product routines have been
  ! moved from efm_module, changed drastically.
  ! will be more changes -- to make it work as "pack"
  !
  subroutine merge_indices_real2D(a,b)
    ! just as pack but the seconon index faster 
    ! -- historically for direct product
    use error_module
    implicit none
    real(RK),intent(in)         :: a(:,:)
    real(RK),intent(out)        :: b(:)
    ! *** end of interface ***

    integer(IK) :: i,j,ij

    if(size(a).ne.size(b))&
         & call error("efmd/merge_indices_real2D: sizes")

    ij = 0
    do i=1,size(a,1)
       do j=1,size(a,2)
          ij = ij + 1
          b(ij) = a(i,j)
       enddo
    enddo
  end subroutine merge_indices_real2D

  subroutine restore_indices_real2D(a,b)
    ! just as unpack but the seconon index faster 
    ! -- historically for direct product
    use error_module
    implicit none
    real(RK),intent(in)         :: a(:)
    real(RK),intent(out)        :: b(:,:)
    ! *** end of interface ***

    integer(IK) :: i,j,ij

    if(size(a).ne.size(b))&
         & call error("efmd/restore_indices_real2D: sizes?")

    ij = 0
    do i=1,size(b,1)
       do j=1,size(b,2)
          ij = ij + 1
          b(i,j) = a(ij)
       enddo
    enddo
  end subroutine restore_indices_real2D

  subroutine restore_indices_cmplx2D(a,b)
    ! just as unpack but the seconon index faster 
    ! -- historically for direct product
    use error_module
    implicit none
    complex(CK),intent(in)   :: a(:)
    complex(CK),intent(out)  :: b(:,:)
    ! *** end of interface ***

    integer(IK) :: i,j,ij

    if(size(a).ne.size(b))&
         & call error("efmd/restore_indices_cmplx2D: sizes?")

    ij = 0
    do i=1,size(b,1)
       do j=1,size(b,2)
          ij = ij + 1
          b(i,j) = a(ij)
       enddo
    enddo
  end subroutine restore_indices_cmplx2D

  subroutine merge_indices_real4D(a,b)
    ! just as pack but the seconon index faster 
    ! -- historically for direct product
    use error_module
    implicit none
    real(RK),intent(in)  :: a(:,:,:,:)
    real(RK),intent(out) :: b(:,:)
    ! *** end of interface ***

    integer(IK) :: i,j,ij

    if(size(b,1).ne.size(a,1)*size(a,2)&
         & .or. size(b,2).ne.size(a,3)*size(a,4)&
         & ) call error("efmd/merge_indices_real4D: sizes?")

    ij = 0
    do i=1,size(a,1)
       do j=1,size(a,2)
          ij = ij + 1
          call merge_indices(a(i,j,:,:),b(ij,:))
       enddo
    enddo
  end subroutine merge_indices_real4D
  
  !*************************************************************
  subroutine direct_product_rr(a,b,ab)
    !  Purpose: computes the direct product of two matrices
    !
    ! !!! must be in accordance with restore_indices !!!
    !
    !------------ Modules used -----------------------------------
    use error_module
    !------------ Declaration of formal parameters ---------------
    real(RK),    intent(in)  :: a(:,:),b(:,:)
    ! matrices to be multiplied
    real(RK),    intent(out) :: ab(:,:)
    ! product matrix
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------

    integer(IK) :: ia,ja,ib,jb,iab,jab

    !------------ Executable code ------------------------------------

    iab = 0
    ia_:do ia=1,size(a,1)
       ib_:do ib=1,size(b,1)
          iab = iab + 1
          jab = 0
          ja_:do ja=1,size(a,2)
             jb_:do jb=1,size(b,2)
                jab = jab + 1
                ab(iab,jab)      = a(ia,ja) * b(ib,jb)
             enddo jb_
          enddo ja_
       enddo ib_
    enddo ia_
  end subroutine direct_product_rr
  !*************************************************************

  !*************************************************************
  subroutine direct_product_rc(a,b,ab)
    !  Purpose: computes the direct product of two matrices
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(RK),intent(in)     :: a(:,:)
    complex(CK),intent(in)  :: b(:,:)
    ! matrices to be multiplied
    complex(CK),intent(out) :: ab(:,:)
    ! product matrix
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------

    integer(IK) :: ia,ja,ib,jb,iab,jab

    !------------ Executable code ------------------------------------

    iab = 0
    ia_:do ia=1,size(a,1)
       ib_:do ib=1,size(b,1)
          iab = iab + 1
          jab = 0
          ja_:do ja=1,size(a,2)
             jb_:do jb=1,size(b,2)
                jab = jab + 1
                ab(iab,jab)      = a(ia,ja) * b(ib,jb)
             enddo jb_
          enddo ja_
       enddo ib_
    enddo ia_
  end subroutine direct_product_rc
  !*************************************************************

  !*************************************************************
  subroutine direct_product_cr(a,b,ab)
    !  Purpose: computes the direct product of two matrices
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    complex(CK),intent(in)  :: a(:,:)
    real(RK),intent(in)     :: b(:,:)
    ! matrices to be multiplied
    complex(CK),intent(out) :: ab(:,:)
    ! product matrix
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------

    integer(IK) :: ia,ja,ib,jb,iab,jab

    !------------ Executable code ------------------------------------

    iab = 0
    ia_:do ia=1,size(a,1)
       ib_:do ib=1,size(b,1)
          iab = iab + 1
          jab = 0
          ja_:do ja=1,size(a,2)
             jb_:do jb=1,size(b,2)
                jab = jab + 1
                ab(iab,jab)      = a(ia,ja) * b(ib,jb)
             enddo jb_
          enddo ja_
       enddo ib_
    enddo ia_
  end subroutine direct_product_cr
  !*************************************************************

  !*************************************************************
  subroutine direct_product_cc(a,b,ab)
    !  Purpose: computes the direct product of two matrices
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    complex(CK),intent(in)  :: a(:,:)
    complex(CK),intent(in)  :: b(:,:)
    ! matrices to be multiplied
    complex(CK),intent(out) :: ab(:,:)
    ! product matrix
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------

    integer(IK) :: ia,ja,ib,jb,iab,jab

    !------------ Executable code ------------------------------------

    iab = 0
    ia_:do ia=1,size(a,1)
       ib_:do ib=1,size(b,1)
          iab = iab + 1
          jab = 0
          ja_:do ja=1,size(a,2)
             jb_:do jb=1,size(b,2)
                jab = jab + 1
                ab(iab,jab)      = a(ia,ja) * b(ib,jb)
             enddo jb_
          enddo ja_
       enddo ib_
    enddo ia_
  end subroutine direct_product_cc
  !*************************************************************

  !*************************************************************
  subroutine efm_show_invspace(efm_type_inv)
    ! Purpose : Helper to show the datastructure
    !
    use type_module
    use iounitadmin_module
    use group_module
    !** End of interface *****************************************
    type(efm_cscoII_diag_invspace)      :: efm_type_inv
    !------------ Modules used ------------------- ---------------
    !------------ Declaration of local variables ---------------------
    integer(IK)   :: i,j,i_irrep,i_partner
    integer(IK)   :: dim_trafo,dim_irrep,multiplicity
    external error_handler
    !------------ Declaration of local handles -------------------
    type(efm_csco_eigenspace)  , pointer :: efm_type_inv_eigenspace
    type(efm_cscoII_eigenspace), pointer :: efm_type_inv_subspace
    logical :: io
    !------------ Executable code ------------------------------------

    ! Slaves do not always have output_unit open:
    io = output_unit > 0 .and. .not. no_output_unit_output

    dim_trafo = efm_type_inv%dimension
    do i_irrep = 1,group_num_ir
       efm_type_inv_eigenspace => efm_type_inv%csco_eigenspaces(i_irrep)
       if (io) then
           write(output_unit,*)
           write(output_unit,*) " Irrep: ",i_irrep,"  ",irrep_can(i_irrep)%label
       endif
       if (efm_type_inv_eigenspace%exists) then
          efm_type_inv_subspace => efm_type_inv_eigenspace%cscoII_subspaces(1)
          if (io) then
              write(output_unit,*) " exists:"
          endif
          dim_irrep = irrep_can(i_irrep)%dimension
          multiplicity = efm_type_inv_subspace%multiplicity
          do i_partner = 1,dim_irrep
             efm_type_inv_subspace => &
                  efm_type_inv_eigenspace%cscoII_subspaces(i_partner)
             if (io) then
                 write(output_unit,*) "  Partner:",i_partner
                 do i = 1,multiplicity
                    write(output_unit,602)&
                         (efm_type_inv_subspace%eigenfunction(j,i),j=1,dim_trafo)
                 end do
             endif
          end do
       else
          if (io) then
              write(output_unit,*) " does not exist."
          endif
       endif
    end do

602 format(1x,f17.12,1x,f17.12,2x,4(f17.12,f17.12:)/4(10x,4(f17.12,f17.12:)/))
  end subroutine efm_show_invspace
  !*************************************************************


  !------------------------------------------------------------------
  !
  ! GAUGE TRANSFORMATION ROUTINES >>>
  !

  !*************************************************************
  subroutine  efm_inv_cgauge_io(trafos,gauge)
    use group_module, only: group_num_el,contains_inversion
    implicit none
    complex(CK),intent(inout)   :: trafos(:,:,:)
    complex(CK),intent(in)      :: gauge
    ! *** end of interface ***

    integer(IK) :: shp(3)
    integer(IK) :: i

    shp(:)=shape(trafos)
    
    if(shp(1)/=shp(2))&
         & call error_handler("efm/cgauge_io: non square matrix")
    if(shp(3)/=group_num_el)&
         & call error_handler("efm/cgauge_io: wrong num of elements")

    do i=1,group_num_el
       if(contains_inversion(i))then
          trafos(:,:,i) = gauge * trafos(:,:,i)
       endif
    enddo
  end subroutine efm_inv_cgauge_io

  !*************************************************************
  subroutine  efm_inv_rgauge_io(trafos,gauge)
    use group_module, only: group_num_el,contains_inversion
    implicit none
    real(RK),intent(inout)   :: trafos(:,:,:)
    real(RK),intent(in)      :: gauge
    ! *** end of interface ***

    integer(IK) :: shp(3)
    integer(IK) :: i

    shp(:)=shape(trafos)

    if(shp(1)/=shp(2))&
         & call error_handler("efm/rgauge_io: non square matrix")
    if(shp(3)/=group_num_el)&
         & call error_handler("efm/rgauge_io: wrong num of elements")

    do i=1,group_num_el
       if(contains_inversion(i))then
          trafos(:,:,i) = gauge * trafos(:,:,i)
       endif
    enddo
  end subroutine efm_inv_rgauge_io

  !-------------------------------------------------------
  !
  ! SELECT A PHASE FOR EIGENFUNCTIONS>>>
  !

  subroutine efm_rephase(inv)
    use group_module, only: group_num_re
    implicit none
    type(efm_cscoII_diag_invspace),intent(inout) :: inv
    ! *** end of interface ***

    integer     :: memstat
    integer(IK) :: n_preps,dim_trafo,irr,n_indep,ifx,dim_irrep
    complex(CK),pointer               :: F(:,:)
    type(efm_csco_eigenspace),pointer :: eigs

    n_preps = group_num_re

    dim_trafo     = inv%dimension

    do irr=1,n_preps

       eigs => inv%csco_eigenspaces(irr)

       if(.NOT. eigs%exists)then
          cycle
       endif

       dim_irrep       = eigs%dim_irrep

       n_indep         = multiplicity(eigs)

       allocate(F(dim_trafo,dim_irrep),STAT=memstat)
       if(memstat/=0)call error_handler("efm/efm_rephase: alloc failed")

       do ifx=1,n_indep
          call extract_basis(F,eigs,INFN=ifx)
          F = phase(F) * F
          call load_basis(F,eigs,INFN=ifx)
       enddo

       deallocate(F,STAT=memstat)
       if(memstat/=0)call error_handler("efm/efm_rephase: dealloc failed")
    enddo
  end subroutine efm_rephase

  function phase(u) result(ph)
    implicit none
    complex(CK),intent(in) :: u(:,:)
    complex(CK)            :: ph !<<<result
    ! *** end of interface ***

    real(RK)    :: eps
    real(RK)    :: fi(size(u,1),size(u,2))
    integer(IK) :: n
    logical     :: mask(size(u,1),size(u,2))
    complex(CK),parameter :: i=(0.0_ck,1.0_ck)
    complex(CK) :: s

    eps = 10000.0*epsilon(1.0_rk)

    where(abs(u)>eps)
       fi = AIMAG(log(u/abs(u)))
       mask = .true.
    elsewhere
       fi = 0.0_rk
       mask = .false.
    end where

    s = sum(u,MASK=mask)

    if(abs(s)>eps)then
       s = s/abs(s)
       ph = conjg(s)
       return
    else     
       n = count(mask)
       if(n>0)then
          ph = sum(fi)/n
       else
          ph = 0.0_rk
       endif
       ph = -i*exp(-i*ph)
       ! ph = exp(-i*ph)
       return
    endif
  end function phase

  subroutine extract_proj_basis(CU,eigs,infn)
    use group_module
    implicit none
    integer(IK),intent(in)               :: infn
    optional :: infn ! indep function
    type(efm_csco_eigenspace),intent(in) :: eigs
    complex(CK),intent(out)              :: CU(:,:)
    ! *** end of interface ***

    integer(IK) :: n,m,i,jnf
    type(efm_cscoII_eigenspace),pointer :: subs

    jnf=1
    if(present(infn))jnf=infn

    n = size(CU,1)
    m = size(CU,2)

    if(n/=dim_vector(eigs).OR.m/=dim_rep(eigs))&
         & call error_handler("efm/extract_proj_basis: shape conflict")

    do i=1,m
       subs=>eigs%cscoII_subspaces(i)
       CU(1:n,i) = subs%eigenfunction_proj(1:n,jnf)
    enddo
  end subroutine extract_proj_basis

  subroutine extract_vec_basis(RU,eigs,infn)
    use group_module
    implicit none
    integer(IK),intent(in)               :: infn
    optional :: infn ! indep function
    type(efm_csco_eigenspace),intent(in) :: eigs
    real(RK),intent(out)                 :: RU(:,:)
    ! *** end of interface ***

    integer(IK) :: n,m,i,jnf
    type(efm_cscoII_eigenspace),pointer :: subs

    jnf=1
    if(present(infn))jnf=infn

    n = size(RU,1)
    m = size(RU,2)

    if(n/=dim_vector(eigs).OR.m/=dim_rep(eigs))&
         & call error_handler("efm/extract_vec_basis: shape conflict")

    do i=1,m
       subs=>eigs%cscoII_subspaces(i)
       RU(1:n,i) = subs%eigenfunction(1:n,jnf)
    enddo
  end subroutine extract_vec_basis

  subroutine load_proj_basis(CU,eigs,infn)
    use group_module
    implicit none
    integer(IK),intent(in),optional         :: infn ! indep function
    type(efm_csco_eigenspace),intent(inout) :: eigs
    complex(CK),intent(in)                  :: CU(:,:)
    ! *** end of interface ***

    integer(IK) :: i,jnf
    integer(IK) :: n,m
    type(efm_cscoII_eigenspace),pointer :: subs

    n=size(CU,1)
    m=size(CU,2)

    jnf=1
    if(present(infn))jnf=infn

    if(n/=dim_vector(eigs).OR.m/=dim_rep(eigs))&
         & call error_handler("efm/load_proj_basis: shape conflict")

    do i=1,m
       subs=>eigs%cscoII_subspaces(i)
       subs%eigenfunction_proj(1:n,jnf) = CU(1:n,i)
    enddo
  end subroutine load_proj_basis

  subroutine load_vec_basis(RU,eigs,infn)
    use group_module
    implicit none
    integer(IK),intent(in),optional         :: infn ! indep function
    type(efm_csco_eigenspace),intent(inout) :: eigs
    real(RK),intent(in)                     :: RU(:,:)
    ! *** end of interface ***

    integer(IK) :: i,jnf
    integer(IK) :: n,m
    type(efm_cscoII_eigenspace),pointer :: subs

    n=size(RU,1)
    m=size(RU,2)

    jnf=1
    if(present(infn))jnf=infn

    if(n/=dim_vector(eigs).OR.m/=dim_rep(eigs))&
         & call error_handler("efm/load_vec_basis: shape conflict")

    do i=1,m
       subs=>eigs%cscoII_subspaces(i)
       subs%eigenfunction(1:n,jnf) = RU(1:n,i)
    enddo
  end subroutine load_vec_basis

  function dim_rep (eigs) result (d)
    implicit none
    type (efm_csco_eigenspace), intent (in) :: eigs
    integer (IK) :: d           ! result
    ! *** end of interface ***

    d = eigs % dim_irrep

    if (.not. eigs % exists) then
       call error_handler ("efm/dim_rep: doesnt exists ???")
    endif

    if (d .ne. size (eigs % cscoII_subspaces)) then
       call error_handler ("efm/dim_rep: dimension ???")
    endif
  end function dim_rep

  function multiplicity(eigs) result(mult)
    !
    ! Historically if multiplicity is 0
    ! structure is not allocated sometimes
    !
    implicit none
    type(efm_csco_eigenspace),intent(in) :: eigs
    integer(IK)                          :: mult
    ! *** end of interface ***

    if(.not.eigs%exists)then
       mult = 0
!!$         & call error_handler("efm/multiplicity: doesnt exists ???")
    else
       mult  =  eigs%cscoII_subspaces(1)%multiplicity

       if(any(mult /=eigs%cscoII_subspaces(:)%multiplicity ))&
            & call error_handler("efm/multiplicity: different sizes ???")
    endif
  end function multiplicity

  function multiplicities(eigs) result(mult)
    !
    ! Historically if multiplicity is 0
    ! structure is not allocated sometimes
    !
    implicit none
    type(efm_csco_eigenspace),intent(in) :: eigs(:)
    integer(IK)                          :: mult(size(eigs))
    ! *** end of interface ***

    integer(IK) :: irr

    do irr=1,size(eigs)
       mult(irr) = multiplicity(eigs(irr))
    enddo
  end function multiplicities

  function dim_vector(eigs) result(dim)
    implicit none
    type(efm_csco_eigenspace),intent(in) :: eigs
    integer(IK)                          :: dim
    ! *** end of interface ***

    type(efm_cscoII_eigenspace),pointer :: subs

    if(.not.eigs%exists)&
         & call error_handler("efm/dim_vector: doesnt exists ???")

    subs => eigs%cscoII_subspaces(1)
    if(subs%projective)then
       dim  =  size(subs%eigenfunction_proj,1)
    else
       dim  =  size(subs%eigenfunction,1)
    endif
  end function dim_vector

end module efm_decl
