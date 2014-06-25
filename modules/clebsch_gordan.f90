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
module clebsch_gordan
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
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
  use type_module, only:&
       & IK => i4_kind,&
       & RK => r8_kind,&
       & CK => c16_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  type, public :: prod_bas
     !
     ! Expansion coefficients of product basis FA
     ! over factor bases FB and FC               a
     !                     b      c
     !
     ! FA = SUM      CG       *   FB  *  FC
     !   a     b, c    a, b, c      b      c
     !

     !
     ! Number of functions for A, B, and C spaces
     ! can be retrived as
     !
     !  na = size(%c, 1)
     !  nb = size(%c, 2)
     !  nc = size(%c, 3)
     !
     ! in real case of vector irreps and as
     !
     !  na = size(%z, 1)
     !  nb = size(%z, 2)
     !  nc = size(%z, 3)
     !
     ! in complex case of projective irreps.
     !

     !
     ! Actual CG coefficients (real case),
     !
     ! CG        = c(a, b, c)
     !   a, b, c
     !
     real(RK), allocatable :: c(:, :, :) ! (na, nb, nc)

     !
     ! Actual CG coefficients (complex case),
     !
     ! CG        = z(a, b, c) = re(a, b, c) + i * im(a, b, c)
     !   a, b, c
     !
     logical :: cmplx
     real(RK), allocatable :: re(:, :, :) ! (na, nb, nc)
     real(RK), allocatable :: im(:, :, :) ! (na, nb, nc)
     complex(CK), allocatable :: z(:, :, :) ! (na, nb, nc)

     !
     ! if ( cmplx ) then
     !   z(na, nb, nc)
     !   re = re(z), im = im(z)
     ! else
     !   c(na, nb, nc)
     ! endif
     !
  end type prod_bas

  type, public :: sym_prod
     !
     ! multiple subspaces of produc space belonging
     ! to one rep (irrep)
     !
     integer(IK) :: mult
     type(prod_bas), allocatable :: sub(:) ! sub(mult)
  end type sym_prod

  !------------ Declaration of constants and variables ------------

  type(sym_prod), allocatable, public ::   cg(:,:,:) ! Vector * Vector     Irreps
  !                                  cg(n_irr,n_irr,n_irr)
  type(sym_prod), allocatable, public :: vpcg(:,:,:) ! Vector * Projective Irreps
  !                                vpcg(n_pirr,n_virr,n_pirr)
  type(sym_prod), allocatable, public :: vsu2cg(:,:) ! Vector * Projective SU(2) Irrep
  !                                vsu2cg(n_pirr,n_virr)

  ! reordered equivalents, empty irreps moved to the end:
  type(sym_prod), allocatable, target, public ::   cg_reordered(:,:,:) ! Vector * Vector     Irreps
  type(sym_prod), allocatable, target, public :: vpcg_reordered(:,:,:) ! Vector * Projective Irreps
  type(sym_prod), allocatable, target, public :: vsu2cg_reordered(:,:) ! Vector * Projective SU(2) Irrep

  ! shortend equivalents, empty irreps removed, are pointers to parts of "reordered":
  type(sym_prod), pointer, public ::   cg_eliminated(:,:,:) => NULL()   ! Vector * Vector     Irreps
  type(sym_prod), pointer, public :: vpcg_eliminated(:,:,:) => NULL()   ! Vector * Projective Irreps
  type(sym_prod), pointer, public :: vsu2cg_eliminated(:,:) => NULL()   ! Vector * Projective SU(2) Irrep

  !------------ Interface statements ------------------------------

  interface cg_alloc ! public
     module procedure alloc_sym_prod
     module procedure alloc_prod_bas
  end interface

  interface cg_free ! public
     module procedure free_sym_prod
     module procedure free_sym_prod_3D
     module procedure dealloc_prod_bas
  end interface

  !------------ public functions and subroutines ------------------

  public :: clebsch_gordan_new!(mult, na, nb, nc, cmplx) -> sym_prod
  public :: clebsch_gordan_dimensions!(sp, mult, na, nb, nc)
  public :: clebsch_gordan_bcast!(spin_orbit)
  public :: clebsch_gordan_eliminate!(which_vec, which_proj)
  public :: cg_alloc, cg_free, &
       load_efm_to_cg,load_efm_proj_to_cg,&
       show_vsu2cg, &
       show_cg
  public :: clebsch_gordan_close

#ifdef FPP_GFORTRAN_BUGS
  interface assignment(=)
     module procedure assign_sym_prod
  end interface

  public :: assignment(=)
#endif

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !
  ! FIXME: why not representing CG coeffs as a "matrix" or rather an
  ! "array" c(1:mult, 1:na, 1:nb, 1:nc) or similar?
  !

  function clebsch_gordan_new(mult, na, nb, nc, cmplx) result(sp)
    !
    ! Returns a fresh datastruct (aka constructor) for Clebsch-Gordan
    ! coeffs of a reduction (k, A) <- B x C. At the moment, if
    ! multiplicity mult is zero then (na, nb, nc) are unused.
    !
    implicit none
    integer(IK), intent(in) :: mult, na, nb, nc
    logical, optional, intent(in) :: cmplx
    type(sym_prod) :: sp ! result
    ! *** end of interface ***

    integer(IK) :: k

    call alloc_sym_prod(mult, sp)
    do k = 1, mult
       call alloc_prod_bas(na, nb, nc, sp%sub(k), cmplx)
    enddo
  end function clebsch_gordan_new

  subroutine clebsch_gordan_dimensions(sp, mult, na, nb, nc)
    !
    ! Returns multiplicity and dimensions of irreps.  At the moment,
    ! if multiplicity mult appears to be zero, then the dimensions
    ! (na, nb, nc) are unspecified.
    !
    implicit none
    type(sym_prod), intent(in) :: sp
    integer(IK), intent(out) :: mult, na, nb, nc
    ! *** end of interface ***

    mult = size(sp%sub)
    ASSERT(mult==sp%mult)

    if ( mult > 0 ) then
       if ( sp%sub(1)%cmplx ) then
          na = size(sp%sub(1)%z, 1)
          nb = size(sp%sub(1)%z, 2)
          nc = size(sp%sub(1)%z, 3)
       else
          na = size(sp%sub(1)%c, 1)
          nb = size(sp%sub(1)%c, 2)
          nc = size(sp%sub(1)%c, 3)
       endif
    else
       na = -1
       nb = -1
       nc = -1
    endif
  end subroutine clebsch_gordan_dimensions

  subroutine clebsch_gordan_close()
    implicit none
    ! *** end of interface ***

    integer(IK) :: memstat

    !
    ! Recursive deallocations of structures with allocatable
    ! components:
    !

    if ( allocated(cg) ) then
      deallocate(cg, STAT=memstat)
      ASSERT(memstat==0)
    endif

    if ( allocated(vpcg) ) then
      deallocate(vpcg, STAT=memstat)
      ASSERT(memstat==0)
    endif

    if ( allocated(vsu2cg) ) then
      deallocate(vsu2cg, STAT=memstat)
      ASSERT(memstat==0)
    endif

    if ( allocated(cg_reordered) ) then
      deallocate(cg_reordered, STAT=memstat)
      ASSERT(memstat==0)
    endif

    if ( allocated(vpcg_reordered) ) then
      deallocate(vpcg_reordered, STAT=memstat)
      ASSERT(memstat==0)
    endif

    if ( allocated(vsu2cg_reordered) ) then
      deallocate(vsu2cg_reordered, STAT=memstat)
      ASSERT(memstat==0)
    endif

    cg_eliminated => NULL()
    vpcg_eliminated => NULL()
    vsu2cg_eliminated => NULL()
  end subroutine clebsch_gordan_close

  subroutine clebsch_gordan_bcast (spin_orbit)
    implicit none
    logical, intent(in), optional :: spin_orbit
    ! *** end of interface ***

  end subroutine clebsch_gordan_bcast


  subroutine clebsch_gordan_eliminate(which_vec, which_proj)
    use error_module
    implicit none
    logical,intent(in),optional :: which_vec(:), which_proj(:)
    ! *** end of interface ***

    integer(IK) :: memstat
    integer(IK) :: nv,np

    integer(IK)    :: v(20), p(20)

    if(present(which_proj).and.present(which_vec))then
       nv = size(which_vec)
       np = size(which_proj)
       ASSERT(nv<=20)
       ASSERT(np<=20)

       if(np.ne.size(vsu2cg,1).or.nv.ne.size(vsu2cg,2))&
            & call error("cg/clebsch_gordan_eliminate: sizes?")

       allocate( vsu2cg_reordered(np,nv), STAT=memstat)
       ASSERT(memstat==0)

       v(:nv) = neworder(which_vec)
       p(:np) = neworder(which_proj)

       vsu2cg_reordered(:,:) = vsu2cg(p(:np),v(:nv))

       nv = count(.not.which_vec)
       np = count(.not.which_proj)

       vsu2cg_eliminated => vsu2cg_reordered(:np,:nv)

    else if(present(which_vec))then

       nv = size(which_vec)
       ASSERT(nv<=20)

       ASSERT(allocated(cg))
       ASSERT(.not.allocated(cg_reordered))
       allocate(cg_reordered(nv,nv,nv),STAT=memstat)
       ASSERT(memstat==0)

       v(:nv) = neworder(which_vec)
       ! this does not perform a deep copy of pointer components:
       cg_reordered(:,:,:) = cg(v(:nv),v(:nv),v(:nv))

       nv = count(.not.which_vec)
       cg_eliminated => cg_reordered(:nv,:nv,:nv)

    endif

  contains

    function neworder(del) result(order)
      implicit none
      logical, intent(in) :: del(:)
      integer(IK)         :: order(size(del))
      ! *** end of interface ***

      integer(IK) :: i,ii

      ii = 0
      do i=1,size(del)
         if(.not.del(i))then
            ii = ii+1
            order(ii) = i
         endif
      enddo
      do i=1,size(del)
         if(   del(i))then
            ii = ii+1
            order(ii) = i
         endif
      enddo
    end function neworder

  end subroutine clebsch_gordan_eliminate

  subroutine show_vsu2cg(cg)
    implicit none
    type(sym_prod),intent(in) :: cg(:,:)
    ! *** end of interface ***

    integer(IK) :: n_pirr, n_virr, i, j, k

    n_pirr = size(cg, 1)
    n_virr = size(cg, 2)
    print *, 'cg/show_vsu2cg: n_virr=', n_virr
    print *, 'cg/show_vsu2cg: n_pirr=', n_pirr

    do i = 1, n_virr
       print *, 'cg/show_vsu2cg:'
       print *, 'cg/show_vsu2cg: vector irrep ',i,' x SU(2) ='
       do j = 1, n_pirr
          if ( cg(j, i)%mult <= 0 ) cycle
          print *, 'cg/show_vsu2cg: ... projective irrep ', j, ' times ', cg(j, i)%mult
          do k = 1, cg(j, i)%mult
             print *, 'cg/show_vsu2cg: .... subspace ', k, ' :'
             call show_prod_bas(cg(j, i)%sub(k)%z)
          enddo
       enddo
    enddo
  contains
    subroutine show_prod_bas(z)
      implicit none
      complex(CK), intent(in) :: z(:, :, :) ! (na, nb, nc)
      ! *** end of interface ***

      integer(IK) :: pa, pc

      do pa = 1, size(z, 1) ! na
         print *, 'cg/show_vsu2cg: ..... partner ', pa
         do pc = 1, size(z, 3) ! nc == 2
            print *, 'cg/show_vsu2cg: ...... pc=', pc, ' re=', real(z(pa, :, pc))
            print *, 'cg/show_vsu2cg: ...... pc=', pc, ' im=', aimag(z(pa, :, pc))
         enddo
      enddo
    end subroutine show_prod_bas
  end subroutine show_vsu2cg

  subroutine show_cg(cg)
    implicit none
    type(sym_prod),intent(in) :: cg(:,:,:)
    ! *** end of interface ***

    integer(IK)            :: n_airr,n_birr,n_cirr,i,j,k,ij

    ! A <- B x C
    n_airr = size(cg,1)
    n_birr = size(cg,2)
    n_cirr = size(cg,3)
    print *,'cg/show_cg: n_airr=',n_airr
    print *,'cg/show_cg: n_birr=',n_birr
    print *,'cg/show_cg: n_cirr=',n_cirr

    do i=1,n_birr
       print *,'cg/show_cg:'
       do j=1,n_cirr
          print *,'cg/show_cg: (',i,') x (',j,') ='
          do ij=1,n_airr
             if ( cg(ij, i, j)%mult <= 0 ) cycle
             print *,'cg/show_cg: ... (',ij,') times ',cg(ij,i,j)%mult
             do k=1,cg(ij,i,j)%mult
                print *,'cg/show_cg: .... instance ',k,' :'
                call show_prod_bas(cg(ij, i, j)%sub(k)%c)
             enddo
          enddo
       enddo
    enddo
  contains
    subroutine show_prod_bas(c)
      implicit none
      real(RK), intent(in) :: c(:, :, :) ! (na, nb, nc)
      ! *** end of interface ***

      integer(IK) :: pa, pb, i

      ! A <- B x C
      do pa = 1, size(c, 1) ! na
        write(*,*) 'cg/show_cg: ..... partner ',pa
        write(*,'(" cg/show_cg: ..... p1\p2:" ,8I11)') (i, i=1, size(c, 3)) ! nc
        do pb = 1, size(c, 2) ! nb
        write(*,'(" cg/show_cg: ..... ",I7," ",8(F10.5," "))')&
             & pb, ( c(pa, pb, i), i=1, size(c, 3) ) ! nc
        enddo
      enddo
    end subroutine show_prod_bas
  end subroutine show_cg

  subroutine load_efm_to_cg(n_bf,n_cf,eis,sp)
    use error_module
    use efm_decl
    implicit none
    integer(IK),intent(in)               :: n_bf,n_cf
    ! eigenspace eis is a (reduced) subspace of a
    ! product space A = B * C.  
    ! n_bf,n_cf are dimensions of vectors in B and C
    type(efm_csco_eigenspace),intent(in) :: eis
    type(sym_prod),intent(inout)         :: sp
    ! *** end of interface ***

    integer(IK)      :: memstat
    integer(IK)      :: mult,i,veclen,n_af,f
    real(RK),pointer :: bas(:,:)

    mult = multiplicity(eis)
    call alloc_sym_prod(mult, sp)
    if(mult>0)then
       ! otherwise eis structure may be not allocated

       ! dimension of product space
       veclen  = dim_vector(eis)
       n_af    = dim_rep(eis) 

       if(veclen/=n_bf*n_cf)&
            & call error("cg/load_efm_to_cg: dims?")

       allocate(bas(veclen,n_af),STAT=memstat)
       ASSERT(memstat==0)

       do i=1,mult
          call alloc_prod_bas(n_af, n_bf, n_cf, sp%sub(i))
          call extract_basis(bas,eis,infn=i)
          do f=1,n_af
             call restore_indices(bas(:,f),sp%sub(i)%c(f,:,:))
          enddo
       enddo

       deallocate(bas,STAT=memstat)
       ASSERT(memstat==0)
    endif
  end subroutine load_efm_to_cg

  subroutine load_efm_proj_to_cg(n_bf,n_cf,eis,sp)
    use error_module
    use efm_decl
    implicit none
    integer(IK),intent(in)               :: n_bf,n_cf
    ! eigenspace eis is a (reduced) subspace of a
    ! product space A = B * C.  
    ! n_bf,n_cf are dimensions of vectors in B and C
    type(efm_csco_eigenspace),intent(in) :: eis
    type(sym_prod),intent(inout)         :: sp
    ! *** end of interface ***

    integer(IK)      :: memstat
    integer(IK)      :: mult,i,veclen,n_af,f
    complex(CK),pointer :: bas(:,:)

    mult = multiplicity(eis)
    call alloc_sym_prod(mult, sp)
    if(mult>0)then
       ! otherwise eis structure may be not allocated

       ! dimension of product space
       veclen  = dim_vector(eis)
       n_af    = dim_rep(eis) 

       if(veclen/=n_bf*n_cf)&
            & call error("cg/load_efm_proj_to_cg: dims?")

       allocate(bas(veclen,n_af),STAT=memstat)
       ASSERT(memstat==0)

       do i=1,mult
          call alloc_prod_bas(n_af, n_bf, n_cf, sp%sub(i), cmplx=.true.)
          call extract_basis(bas,eis,infn=i)
          do f=1,n_af
             call restore_indices(bas(:,f),sp%sub(i)%z(f,:,:))
          enddo
          sp%sub(i)%re =  REAL(sp%sub(i)%z)
          sp%sub(i)%im = AIMAG(sp%sub(i)%z)
       enddo

       deallocate(bas,STAT=memstat)
       ASSERT(memstat==0)
    endif
  end subroutine load_efm_proj_to_cg

  subroutine alloc_prod_bas(n_af,n_bf,n_cf,pb,cmplx)
    implicit none
    integer(IK),intent(in)       :: n_af,n_bf,n_cf
    type(prod_bas),intent(inout) :: pb
    logical,optional,intent(in)  :: cmplx
    ! *** end of interface ***

    integer(IK) :: memstat

    pb%cmplx = present(cmplx)
    if(present(cmplx)) pb%cmplx = cmplx

    if(pb%cmplx)then
       allocate(&
            & pb%z(n_af,n_bf,n_cf),&
            & pb%re(n_af,n_bf,n_cf),&
            & pb%im(n_af,n_bf,n_cf),&
            & STAT=memstat&
            & )
       ASSERT(memstat==0)
    else
       allocate(pb%c(n_af,n_bf,n_cf),STAT=memstat)
       ASSERT(memstat==0)
    endif
  end subroutine alloc_prod_bas

  subroutine dealloc_prod_bas(pb)
    implicit none
    type(prod_bas),intent(inout) :: pb
    ! *** end of interface ***

    integer(IK) :: memstat

    if(pb%cmplx)then
       deallocate(pb%z,pb%re,pb%im,STAT=memstat)
       ASSERT(memstat==0)
    else
       deallocate(pb%c,STAT=memstat)
       ASSERT(memstat==0)
    endif
    pb%cmplx = .false.
  end subroutine dealloc_prod_bas

  subroutine alloc_sym_prod(mult, sp)
    implicit none
    integer(IK),intent(in)       :: mult
    type(sym_prod), intent(out) :: sp
    ! *** end of interface ***

    integer(IK) :: memstat

    sp%mult = mult
    allocate(sp%sub(mult), STAT=memstat)
    ASSERT(memstat==0)
  end subroutine alloc_sym_prod

  subroutine dealloc_sym_prod(sp)
    implicit none
    type(sym_prod),intent(inout) :: sp
    ! *** end of interface ***

    integer(IK) :: memstat

    deallocate(sp%sub, STAT=memstat)
    ASSERT(memstat==0)
    sp%mult   = -1
  end subroutine dealloc_sym_prod

  subroutine free_sym_prod(sp)
    implicit none
    type(sym_prod),intent(inout) :: sp
    ! *** end of interface ***

    integer(IK) :: i

    do i=1,sp%mult
       call dealloc_prod_bas(sp%sub(i))
    enddo
    call dealloc_sym_prod(sp)
  end subroutine free_sym_prod

  subroutine free_sym_prod_3d(sp)
    implicit none
    type(sym_prod),intent(inout) :: sp(:,:,:)
    ! *** end of interface ***

    integer(IK) :: i,j,k

    do i=1,size(sp,1)
       do j=1,size(sp,2)
          do k=1,size(sp,3)
             call free_sym_prod(sp(i, j, k))
          enddo
       enddo
    enddo
  end subroutine free_sym_prod_3d

#ifdef FPP_GFORTRAN_BUGS
  subroutine assign_sym_prod(b, a)
    !
    ! Gfortran 4.3 (Debian Lenny) has troubles with TR18581 assignment
    ! of zero sized components.
    !
    implicit none
    type(sym_prod), intent(out) :: b
    type(sym_prod), intent(in) :: a
    ! *** end of interface ***

    call alloc_sym_prod(a%mult, b)
    b%sub(:) = a%sub(:)
  end subroutine assign_sym_prod
#endif

  !--------------- End of module ----------------------------------
end module clebsch_gordan
