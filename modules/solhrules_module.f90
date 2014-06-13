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
module  solhrules_module
  !---------------------------------------------------------------
  !
  !  Purpose: Coding of Product and differentiation rule of solid
  !           harmonics of arbitrary l
  !
  !  Prerequesits: solhrules_setup() must be called when 
  !                unique_atom_lmax is available
  !
  !  This module uses the meta-index lm to describe the
  !  indices (l,m). That should allow larger loops and faciliate
  !  loop structure. The maping of (l,m) to lm is done in the
  !  following way:
  !     lm = 1
  !     do l = 0, lmax
  !        do m = 1, 2*l+1
  !           xx(lm) = ??
  !           lm = lm + 1
  !        enddo
  !     enddo
  !  The maping is described by the array solhrules_l_and_m_of_lm(:,:)
  !  and by the subroutine solhrules_lm_of_l_and_m(lm,l,m).
  !
  !  The storage requirements for the product rule are
  !  about proportional to (l_max + 1) ** 3 and for the differential
  !  rule proportional to (l_max + 1) ** 4. Thus, for high l_max
  !  it might become necessary to recalculate differential rules for
  !  selected (lm1,lm2) from product rule as needed instaed of keeping all
  !  differential rules in memory. To allow such a strategy, routines
  !  to calculate and free the differential rules for single (lm1,lm2)
  !  were made public.
  !
  !
  !  This module also includes a number of functions that apply the
  !  product and differential rules in various ways. All of them
  !  perform their operations for an array of values simultaneusly.
  !  This is intended for use with the a,b exponent metaindex of
  !  the primitive integral subroutines.
  !
  !  
  !
  !  Module called by: in the first place by primitive integral subroutines:
  !   ss_calculate, ls_calculate, ll_calculate and derived routines.
  !
  !
  !  References: PhD-thesis A. Goerling
  ! 
  !
  !  Author: TB
  !  Date: 1/95
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   6/96
  ! Description: Two functions for easy implimentation of product- 
  !              and differentialrule have been added 
  !
  !----------------------------------------------------------------
  !
  ! Modification 
  ! Author: AM
  ! Date:   02/1999
  ! Description: scalar version and common interface to 
  !              diff_rule subroutine have been added
  !
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
# include "def.h"
  use type_module, only: IK => i4_kind, RK => r8_kind ! type specification parameters

  implicit none
  private         ! by default, all names are private
  save
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------
  type, public ::  solhrules_product_type
     ! to implement the product rule one given lm_result of solid harmonic
     ! solid harmonic of gradient operator
     ! sh(lm_result) = sum[ coef  * (sh_grad(lm_sh1) f) * (sh_grad(lm_sh2) g) ]
     integer(IK)          :: n_summands
     integer(IK), allocatable :: lm_sh1(:) ! lm_sh1(n_summands)
     ! lm index of first solid harmonic contributing
     integer(IK), allocatable :: lm_sh2(:) ! lm_sh2(n_summands)
     ! lm index of secound solid harmonic contributing
     real(RK), allocatable :: coef(:) ! coef(n_summands)
     ! coefficients to multiply products of sh1 and sh2 with
  end type solhrules_product_type

  type, public ::  solhrules_differential_type
     ! to implement the differential rule of one given (lm_grad,lm_res)
     ! with lm_grad: index of solid harmonic of gradient operator
     ! and lm_res: index of solid harmonic that gradient operator should act upon
     ! sh_grad(lm_grad) sh(lm_res) = sum[ coef * sh(lm_sh) ]
     integer(IK)          :: n_summands
     integer(IK), allocatable :: lm_sh(:) ! lm_sh(n_summands)
     ! lm index of solid harmonic contributing
     real(RK), allocatable :: coef(:) ! coef(n_summands)
     ! coefficients to multiply sh with
  end type solhrules_differential_type

    type,public :: nested2_vars
    SEQUENCE
    integer(IK) ::i
    integer(IK) ::l1
    integer(IK) ::l3
    integer(IK) ::index_p
    integer(IK) ::index0_p
    integer(IK) ::index1_p
    integer(IK) ::index2_p
    integer(IK) ::index3_p
    integer(IK) ::index3_p2
    real(RK):: coef_prod
    end type nested2_vars


  !------------ Declaration of constants and variables ------------
  type(solhrules_product_type), allocatable, target, public :: solhrules_product(:)
  ! solhrules_product(lm_max), lm index of solid harmonic of gradient operator
  type(solhrules_differential_type), allocatable, target, public :: solhrules_differential(:,:)
  ! solhrules_differential(lm_max,lm_max)

!   protected :: solhrules_product, solhrules_differential
  !
  ! FIXME: we cannot set protected attribute here because of
  !        pointer assignments like this
  !
  !             shrd => solhrules_differential(2, lm2)
  !
  !        being used at a few places over the code.
  !

  ! first lm: index of solid harmonic of gradient operator
  ! secound lm: index of solid harmonic that gradient operator should act upon
  integer(IK), allocatable, public, protected :: solhrules_l_and_m_of_lm(:,:)
  ! solhrules_l_and_m_of_lm(2,lm_max)
  ! first dimension: 1 l, 2 m

  !------------ public functions and subroutines ------------------

  public :: solhrules_setup!(lmax)
  public :: solhrules_free!()

! public :: solhrules_lm_of_l_and_m
  public :: prod_rule
  public :: diff_rule

  public :: prod_rule_double
  public :: prod_rule_nested2

  public :: diff_rule_nested

  public :: solhrules_print
! public :: solhrules_calculate_diff
! public :: prod_rule2
! public :: prod_rule3
! public :: prod_rule_nested
! public :: prod_rule_double3
! public :: prod_rule_nested3
! public :: prod_rule_cross
! public :: prod_rule_scalar
! public :: nsum_prod_rule_nested2
! public :: summands_prod_rule_nested2
! public :: opt_prod_rule_nested2
! public :: nsums_prod_rule_nested2
! public :: gsummands_prod_rule_nested2

  interface diff_rule

     module procedure diff_rule_sclr
     !     scalar variant:
     !function diff_rule_sclr(f,lm_min1,lm_max1,lm_2) result(p)
     !  ! p(lm_grad) = sh_grad(lm_grad) f(lm_2) = sum[ coef * f(lm_sh) ]
     !  integer(i4_kind),intent(in)             :: lm_min1,lm_max1,lm_2
     !  real(r8_kind),dimension(:),intent(in)   :: f
     !  real(r8_kind),dimension(lm_min1:lm_max1):: p ! result
     !end function diff_rule_sclr       

     module procedure diff_rule_vec
     !     vector variant:
     !function diff_rule_vec(f,lm_min1,lm_max1,lm_2) result(p)
     !  ! p(lm_grad) = sh_grad(lm_grad) f(lm_2) = sum[ coef * f(lm_sh) ]
     !  integer(i4_kind),intent(in)                        :: lm_min1,lm_max1,lm_2
     !  real(r8_kind),dimension(:,:),intent(in)            :: f
     !  real(r8_kind),dimension(size(f,1),lm_min1:lm_max1) :: p ! result
     !end function diff_rule_vec
  end interface

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of constants and variables ----

  !
  ! This is set in solhrules_setup(lmax), reset to illegal value
  ! upon cleanup:
  !
  integer(IK) :: l_max = -1
  integer(IK) :: lm_max ! lm_max = (l_max + 1)**2
  integer(IK), allocatable :: solhrules_first_lm_of_l(:)
  ! solhrules_first_lm_of_l(l_max)
  real(RK), allocatable    :: solhrules_dfac(:)
  ! solhrules_dfac(lm_max),  (2*l + 1)!!



  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine solhrules_setup(lmax)
    !
    !  Purpose: must be called to initialize module.
    !  Does allocation and calculates coefficints and indices in 
    !  solhrules_product(:), solhrules_differential(:,:) and
    !  solhrules_l_and_m_of_lm(:,:)
    !
    !  NOTE: called with lmax = unique_atom_lmax_all
    !
    !------------ Modules used ---------------------------------
    ! use unique_atom_module, only: unique_atom_lmax_all
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(IK), intent(in)    :: lmax
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(IK) :: status
    !------------ Executable code ------------------------------

!   print *, "solhrules_setup(", lmax, ")"
!   print *, "solhrules_setup: l_max=", l_max

    !
    ! Otherwise we may need to call solhrules_free() first:
    !
    if ( l_max >= 0 ) then
        ASSERT(l_max==lmax)
        WARN('why calling?')
        return
    endif

    !
    ! Set module global vars (also indicates that the module
    ! was initialized):
    !
    l_max = lmax
    lm_max = (lmax + 1)**2

    allocate( solhrules_first_lm_of_l(0:l_max), stat=status )
    ASSERT(status==0)

    allocate( solhrules_l_and_m_of_lm(2,lm_max), stat=status )
    ASSERT(status==0)

    allocate( solhrules_differential(lm_max,lm_max), stat=status )
    ASSERT(status==0)

    allocate( solhrules_dfac(lm_max), stat=status )
    ASSERT(status==0)

    call solhrules_setup_mapping()

    !
    ! Allocate and fill intent(out), allocatable:
    !
    call solhrules_setup_products(lmax, solhrules_product)

    call solhrules_setup_dfac()

    ! Calculate differentials, seems to be always used:
    call solhrules_setup_diffs()
  end subroutine solhrules_setup

  subroutine solhrules_free()
    !
    ! Deallocate module variables.
    !
    ! NOTE: deallocation of types with allocatbale components
    !       is recursove!
    !
    implicit none
    ! *** end of interface ***

    integer :: status

!    ASSERT(l_max>=0)
    if(l_max < 0) then
       WARN('l_max < 0')
    endif

!   print *, "solhrules_free()"

    if(allocated(solhrules_first_lm_of_l)) then
       deallocate( solhrules_first_lm_of_l, stat=status )
       ASSERT(status==0)
    end if

    if(allocated(solhrules_l_and_m_of_lm)) then
       deallocate( solhrules_l_and_m_of_lm, stat=status )
       ASSERT(status==0)
    end if

    if(allocated(solhrules_product)) then
       deallocate( solhrules_product, stat=status )
       ASSERT(status==0)
    end if

    if(allocated(solhrules_differential)) then
       deallocate( solhrules_differential, stat=status )
       ASSERT(status==0)
    end if

    if(allocated(solhrules_dfac)) then
       deallocate( solhrules_dfac, stat=status )
       ASSERT(status==0)
    end if

    !
    ! Rest module global vars (also indicates that the module
    ! was cleaned up):
    !
    l_max = -1
  end subroutine solhrules_free


  !*************************************************************
  subroutine solhrules_setup_mapping()
    !  Purpose: set up mapping between (l,m) and lm
    !** End of interface ***************************************
    implicit none
    !------------ Declaration of local variables ---------------
    integer(IK)                :: l,m,lm
    !------------ Executable code ------------------------------
    lm = 1
    do l= 0, l_max
       solhrules_first_lm_of_l(l) = lm
       do m = 1, 2*l+1
          solhrules_l_and_m_of_lm(1,lm) = l
          solhrules_l_and_m_of_lm(2,lm) = m
          lm = lm + 1
       enddo
    enddo
  end subroutine solhrules_setup_mapping
  !*************************************************************


  !*************************************************************
  subroutine solhrules_setup_products(lmax, rules)
    !  Purpose: calculates coefficients and indices for product
    !           rule as described in PhD thesis of A. Goerling
    ! This subroutine was derived from program PRODUCT.f 
    ! of A. Goerling, special thanks to the author.
    !** End of interface ***************************************
    implicit none
    integer(IK), intent(in) :: lmax
    type(solhrules_product_type), allocatable, intent(out) :: rules(:)
    ! *** end of interface ***

    integer(IK) :: L, L1, L2, NNN, M, MT, M1
    real(RK)    ::    FACT, FCLBSH, FACTL 
    real(RK), allocatable ::    FAC(:)
    real(RK), parameter   :: SMALL=0.000001_RK

    type(solhrules_product_type) :: shrp

    ! max_dim = lm_max * lm_max
    integer(IK) :: lm_sh1((lmax+1)**4) ! lm_max * lm_max
    integer(IK) :: lm_sh2((lmax+1)**4) ! lm_max * lm_max
    real(RK) :: coef((lmax+1)**4) ! lm_max * lm_max


    !
    ! For all harmonics up to including LMAX:
    !
    allocate(rules((lmax+1)**2))

    ! allocates and fills local FAC(n), n = -lmax-1, ... 3*lmax+1:
    call setup_fac()

    DO L = 0, lmax

       DO M = 0, L

          DO MT = 1, 2
             !
             ! FOR M=0 ONLY MT=1 IS EXISTENT BY DEFINITION:
             !
             IF ( M == 0 .and. MT == 2 ) cycle

             !
             ! This counts the real nomber of non-trivial product rule entries,
             ! we will be putting them at the beginning of these arrays:
             !
             !          lm_sh1(:), lm_sh2(:), coef(:)
             !
             NNN=0

             DO L1=0,L

                L2=L-L1
                FACTL=BIN(L,L1)/CLEBSH(L1,0,L2,0,L,0)


                IF ( M .EQ. 0 ) THEN
                   ! (FOR M=0 ONLY MT=1 IS EXISTENT BY DEFINITION.)
                   ASSERT(MT==1)


                   FCLBSH=CLEBSH(L1,0,L2,0,L,0)
                   IF (ABS(FCLBSH).GT.SMALL) THEN
                      NNN=NNN+1
                      FACT=FACTL*FCLBSH
                      coef(NNN)=FACT
                      lm_sh1(NNN)=NUMB(L1,0,1)
                      lm_sh2(NNN)=NUMB(L2,0,1)
                   ENDIF

                   DO M1=1,L1

                      FCLBSH=CLEBSH(L1,M1,L2,M-M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         FACT=FACTL*FCLBSH
                         IF (MOD(M1,2).EQ.1) THEN
                            FACT=-FACT
                         ENDIF
                         NNN=NNN+1
                         coef(NNN)=FACT
                         lm_sh1(NNN)=NUMB(L1,M1,1)
                         lm_sh2(NNN)=NUMB(L2,M1,1)
                         NNN=NNN+1
                         coef(NNN)=FACT
                         lm_sh1(NNN)=NUMB(L1,M1,2)
                         lm_sh2(NNN)=NUMB(L2,M1,2)
                      ENDIF

                   ENDDO! M1=1,L1



                ELSEIF ((M.NE.0).AND.(MT.EQ.1)) THEN


                   FCLBSH=CLEBSH(L1,0,L2,M,L,M)
                   IF (ABS(FCLBSH).GT.SMALL) THEN
                      FACT=FACTL*FCLBSH
                      NNN=NNN+1
                      coef(NNN)=FACT
                      lm_sh1(NNN)=NUMB(L1,0,1)
                      lm_sh2(NNN)=NUMB(L2,M,1)
                   ENDIF

                   DO M1=1,L1

                      FCLBSH=CLEBSH(L1,M1,L2,M-M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         FACT=FACTL*FCLBSH
                         NNN=NNN+1
                         IF (M.GT.M1) THEN
                            FACT=FACT/SQRT(2.0_RK)
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,1)
                            lm_sh2(NNN)=NUMB(L2,M-M1,1)
                            ELSEIF (M.EQ.M1) THEN
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,1)
                            lm_sh2(NNN)=NUMB(L2,0,1)
                         ELSE
                            FACT=FACT/SQRT(2.0_RK)
                            IF (MOD(M+M1,2).EQ.1) THEN
                               FACT=-FACT
                            ENDIF
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,1)
                            lm_sh2(NNN)=NUMB(L2,M1-M,1)
                         ENDIF
                      ENDIF
                      FCLBSH=CLEBSH(L1,-M1,L2,M+M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         FACT=FACTL*FCLBSH/SQRT(2.0_RK)
                         IF (MOD(M1,2).EQ.1) THEN
                            FACT=-FACT
                         ENDIF
                         NNN=NNN+1
                         coef(NNN)=FACT
                         lm_sh1(NNN)=NUMB(L1,M1,1)
                         lm_sh2(NNN)=NUMB(L2,M+M1,1)
                      ENDIF
                      FCLBSH=CLEBSH(L1,M1,L2,M-M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         IF (M.GT.M1) THEN
                            FACT=FACTL*FCLBSH
                            FACT=-FACT/SQRT(2.0_RK)
                            NNN=NNN+1
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,2)
                            lm_sh2(NNN)=NUMB(L2,M-M1,2)
                            ELSEIF (M.LT.M1) THEN
                            FACT=FACTL*FCLBSH
                            FACT=FACT/SQRT(2.0_RK)
                            NNN=NNN+1
                            IF (MOD(M+M1,2).EQ.1) THEN
                               FACT=-FACT
                            ENDIF
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,2)
                            lm_sh2(NNN)=NUMB(L2,M1-M,2)
                         ENDIF
                      ENDIF
                      FCLBSH=CLEBSH(L1,-M1,L2,M+M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         FACT=FACTL*FCLBSH/SQRT(2.0_RK)
                         IF (MOD(M1,2).EQ.1) THEN
                            FACT=-FACT
                         ENDIF
                         NNN=NNN+1
                         coef(NNN)=FACT
                         lm_sh1(NNN)=NUMB(L1,M1,2)
                         lm_sh2(NNN)=NUMB(L2,M+M1,2)
                      ENDIF

                   ENDDO! M1=1,L1


                ELSEIF ((M.NE.0).AND.(MT.EQ.2)) THEN


                   FCLBSH=CLEBSH(L1,0,L2,M,L,M)
                   IF (ABS(FCLBSH).GT.SMALL) THEN
                      FACT=FACTL*FCLBSH
                      NNN=NNN+1
                      coef(NNN)=FACT
                      lm_sh1(NNN)=NUMB(L1,0,1)
                      lm_sh2(NNN)=NUMB(L2,M,2)
                   ENDIF

                   DO M1=1,L1

                      FCLBSH=CLEBSH(L1,M1,L2,M-M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         IF (M.GT.M1) THEN
                            FACT=FACTL*FCLBSH
                            FACT=FACT/SQRT(2.0_RK)
                            NNN=NNN+1
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,1)
                            lm_sh2(NNN)=NUMB(L2,M-M1,2)
                            ELSEIF (M.LT.M1) THEN
                            FACT=FACTL*FCLBSH
                            FACT=-FACT/SQRT(2.0_RK)
                            NNN=NNN+1
                            IF (MOD(M+M1,2).EQ.1) THEN
                               FACT=-FACT
                            ENDIF
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,1)
                            lm_sh2(NNN)=NUMB(L2,M1-M,2)
                         ENDIF
                      ENDIF
                      FCLBSH=CLEBSH(L1,-M1,L2,M+M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         FACT=FACTL*FCLBSH/SQRT(2.0_RK)
                         IF (MOD(M1,2).EQ.1) THEN
                            FACT=-FACT
                         ENDIF
                         NNN=NNN+1
                         coef(NNN)=FACT
                         lm_sh1(NNN)=NUMB(L1,M1,1)
                         lm_sh2(NNN)=NUMB(L2,M+M1,2)
                      ENDIF
                      FCLBSH=CLEBSH(L1,M1,L2,M-M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         FACT=FACTL*FCLBSH
                         NNN=NNN+1
                         IF (M.GT.M1) THEN
                            FACT=FACT/SQRT(2.0_RK)
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,2)
                            lm_sh2(NNN)=NUMB(L2,M-M1,1)
                            ELSEIF (M.EQ.M1) THEN
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,2)
                            lm_sh2(NNN)=NUMB(L2,0,1)
                            ELSEIF (M.LT.M1) THEN
                            FACT=FACT/SQRT(2.0_RK)
                            IF (MOD(M+M1,2).EQ.1) THEN
                               FACT=-FACT
                            ENDIF
                            coef(NNN)=FACT
                            lm_sh1(NNN)=NUMB(L1,M1,2)
                            lm_sh2(NNN)=NUMB(L2,M1-M,1)
                         ENDIF
                      ENDIF
                      FCLBSH=CLEBSH(L1,-M1,L2,M+M1,L,M)
                      IF (ABS(FCLBSH).GT.SMALL) THEN
                         FACT=-FACTL*FCLBSH/SQRT(2.0_RK)
                         IF (MOD(M1,2).EQ.1) THEN
                            FACT=-FACT
                         ENDIF
                         NNN=NNN+1
                         coef(NNN)=FACT
                         lm_sh1(NNN)=NUMB(L1,M1,2)
                         lm_sh2(NNN)=NUMB(L2,M+M1,1)
                      ENDIF

                   ENDDO! M1=1,L1


                ENDIF


             ENDDO! L1=0,L

             !
             ! Make a fresh entry (type constructor here):
             !
             shrp = solhrules_product_type(NNN, lm_sh1(:NNN), lm_sh2(:NNN), coef(:NNN))

             !
             ! Put a fresh entry into output array at proper position:
             !
             rules(NUMB(L,M,MT)) = shrp

          ENDDO! MT=1,2

       ENDDO! M=0,L

    ENDDO! L=0,LMAX

    call free_fac()

  contains

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine setup_fac()
      ! allocates array fac and calculates faculties
      !++ End of interface +++++++++++++++++++++++++++++++++++++
      implicit none
      integer(IK) :: status,min_val,max_val,i
      !------------ Executable code ----------------------------
      min_val = -1 * lmax - 1
      max_val = 3 * lmax + 1
      allocate( fac(min_val:max_val), stat=status )
      ASSERT(status==0)

      do i = min_val, 1
         fac(i) = 1.0_RK
      enddo
      do i = 2, max_val
         fac(i) = real(i,RK) * fac(i-1)
      enddo
    end subroutine setup_fac
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine free_fac()
      ! deallocates array fac
      !++ End of interface +++++++++++++++++++++++++++++++++++++
      implicit none
      integer(IK) :: status
      !------------ Executable code ----------------------------
      deallocate( fac, stat=status )
      ASSERT(status==0)
    end subroutine free_fac
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    INTEGER(IK) FUNCTION NUMB(L,M,MT)
      ! CALCULATE THE ABSOLUTE NUMBER OF A REAL SOLID HARMONIC
      ! FROM THE 'QUANTUM NUMBERS' L,M,MT.
      IMPLICIT NONE
      !------------ Declaration of formal parameters -----------
      INTEGER(IK) :: L,M,MT
      !++ End of interface +++++++++++++++++++++++++++++++++++++
      INTEGER(IK) :: N
      !------------ Executable code ----------------------------
      N=L*L
      IF (M.EQ.0) THEN
         IF (MT.EQ.1) THEN
            NUMB=N+1
         ELSE
            NUMB = 0 ! to make compiler happy
            call error_handler( &
                 "solhrules_setup_products: NUMB" )
         ENDIF
      ELSE
         N=N+2*M
         NUMB=N+MT-1
      ENDIF
    END FUNCTION NUMB
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    REAL(RK) FUNCTION BIN(L,L1)
      ! CALCULATE THE BINOMIAL COEFFICIENT OF L ABOVE L1.
      IMPLICIT NONE
      !------------ Declaration of formal parameters -----------
      INTEGER(IK) :: L,L1
      !++ End of interface +++++++++++++++++++++++++++++++++++++
      INTEGER(IK) :: I
      REAL(RK)    :: S, SL, SL1, SLHELP
      !------------ Executable code ---------------------------
      SL=REAL(L,RK)
      SL1=REAL(L1,RK)
      IF ((L1.GE.L).OR.(L1.LT.0)) THEN
         BIN=1.0_RK
      ELSE
         SLHELP=1.0_RK
         S=1.0_RK
         DO I=L1+1,L
            SL1=SL1+1.0_RK
            S=S*SL1/SLHELP
            SLHELP=SLHELP+1.0_RK
         ENDDO
         BIN=S
      ENDIF
    END FUNCTION BIN
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    REAL(RK) FUNCTION CLEBSH(L1,M1,L2,M2,L,M)
      ! CALCULATES A CLEBSH GORDON COEFFICIENT FOR L1+L2=L.
      IMPLICIT NONE
      !------------ Declaration of formal parameters -------------
      INTEGER(IK) :: L1,M1,L2,M2,L,M
      !++ End of interface +++++++++++++++++++++++++++++++++++++
      REAL(RK)    :: HILF
      !------------ Executable code ------------------------------
      IF ((L1+L2).NE.L) THEN
         CLEBSH = 0.0 ! to make compiler happy
         call error_handler( &
              "solhrules_setup_products: ERROR 1 IN CLEBSH" )
      ELSEIF (ABS(M1).GT.L1) THEN
         CLEBSH = 0.0 ! to make compiler happy
         call error_handler( &
              "solhrules_setup_products: ERROR 2 IN CLEBSH" )
      ELSEIF ((ABS(M2).GT.L2).OR.(ABS(M).GT.L)) THEN
         CLEBSH=0.0_RK
      ELSEIF ((M1+M2).NE.M) THEN
         CLEBSH=0.0_RK
      ELSE
         HILF=FAC(2*L1)*FAC(2*L2)/FAC(2*L)
         HILF=HILF*FAC(L+M)*FAC(L-M)/(FAC(L1+M1)*FAC(L1-M1))
         HILF=HILF/(FAC(L2+M2)*FAC(L2-M2))
         CLEBSH=SQRT(HILF)
      ENDIF
    END FUNCTION CLEBSH
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end subroutine solhrules_setup_products
  !*************************************************************


  !*************************************************************
  subroutine solhrules_setup_dfac()
    !  Purpose: calculates array solhrules_dfac(lm) = (2*l(lm) + 1)!!
    !** End of interface *****************************************
    implicit none
    !------------ Declaration of local variables -----------------
    integer(IK)                :: lm,l,m,f
    real(RK)                   :: dfac  
    !------------ Executable code --------------------------------
    solhrules_dfac(1) =  1.0_RK
    f = 1
    lm = 1
    do l = 1, l_max
       f = f * ( 2*l - 1 )
       dfac = real(f,RK)
       do m = 1, 2*l + 1
          lm = lm + 1
          solhrules_dfac(lm) = dfac
       enddo
    enddo
  end subroutine solhrules_setup_dfac
  !*************************************************************


  !*************************************************************
  subroutine solhrules_setup_diffs()
    !  Purpose: calculates oefficients and indices for differential
    !  rules of all (lm1,lm2) from coefficients and indices 
    !  for product rule
    !** End of interface *****************************************
    implicit none
    !------------ Declaration of local variables -----------------
    integer(IK)                :: lm1, lm2
    !------------ Executable code --------------------------------

    do lm1 = 1, lm_max
       do lm2 = 1, lm_max
          call solhrules_calculate_diff(lm1,lm2)
       enddo
    enddo
  end subroutine solhrules_setup_diffs
  !*************************************************************


  !*************************************************************
  subroutine solhrules_calculate_diff(lm1,lm2)
    !  Purpose: calculates oefficients and indices for differential
    !  rule from coefficients and indices for product rule
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(IK), intent(in)  :: lm1,lm2
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(IK)                      :: status, i, n
    real(RK)                         :: dfac
    type(solhrules_product_type), pointer      :: shrp
    type(solhrules_differential_type), pointer :: shrd
    !------------ Executable code --------------------------------
    shrp => solhrules_product(lm2)
    shrd => solhrules_differential(lm1,lm2)

    n = 0
    do i=1,shrp%n_summands
       if ( shrp%lm_sh1(i) .eq. lm1 ) n = n + 1
    enddo
    shrd%n_summands = n

    if ( n .gt. 0 ) then
       allocate(shrd%lm_sh(n), shrd%coef(n), stat=status )
       ASSERT(status==0)

       dfac = solhrules_dfac(lm1)
       n = 0
       do i=1,shrp%n_summands
          if ( shrp%lm_sh1(i) .eq. lm1 ) then
             n = n + 1
             shrd%lm_sh(n) = shrp%lm_sh2(i)
             shrd%coef(n)  = shrp%coef(i) * dfac
             if ( n .eq. shrd%n_summands ) exit
          endif
       enddo
    endif

  end subroutine solhrules_calculate_diff
  !*************************************************************

  !*************************************************************
  subroutine solhrules_print(io_unit)
    !  Purpose: prints debug output
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(IK),intent(in) :: io_unit
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(IK)                :: lm,lm1,lm2,i,l,m
    type(solhrules_product_type), pointer      :: shrp
    type(solhrules_differential_type), pointer :: shrd
    !------------ Executable code ------------------------------

    !
    ! Otherwise module is not initialized:
    !
    ASSERT(l_max>=0)

    write(io_unit,*) 
    write(io_unit,*) '*********************************************************'
    write(io_unit,*) 
    write(io_unit,*) "MAPPING OF L AND M TO LM AND SOLHRULES_DFAC(LM) = (2*L(LM) + 1)!!"
    write(io_unit,*)
    write(io_unit,*) '*********************************************************'
    write(io_unit,*)
    write(io_unit,'(A4,X,A4,X,A4,2X,A20)') "L","m","Lm","dfac(Lm)"
    do l=0,l_max
       do m = 1, 2*l + 1
          call solhrules_lm_of_l_and_m(lm,l,m)
          write(io_unit,'(I4,X,I4,X,I4,2X,F20.0)') l, m, lm, solhrules_dfac(lm)
       enddo
    enddo


    write(io_unit,*) 
    write(io_unit,*) '*********************************************************'
    write(io_unit,*)
    write(io_unit,*) "PRODUCT RULES FOR SOLID HARMONICS UP TO L_MAX", L_MAX
    write(io_unit,*)
    write(io_unit,*) '*********************************************************'
    write(io_unit,*) 
    do lm = 1,lm_max
       shrp => solhrules_product(lm)
       write(io_unit,*)
       write(io_unit,'("|",I2,",",I2,">/(",I2,") == ",I2," summands")') &
            solhrules_l_and_m_of_lm(1,lm), solhrules_l_and_m_of_lm(2,lm), lm, &
            shrp%n_summands
       do i=1,shrp%n_summands
          write(io_unit,'(F20.10,X,"|",I2,",",I2,">/(",I2,") and |",I2,",",I2,">/(",I2,")")') &
                shrp%coef(i),   &
                solhrules_l_and_m_of_lm( 1, shrp%lm_sh1(i) ),   &
                solhrules_l_and_m_of_lm( 2, shrp%lm_sh1(i) ),   &
                shrp%lm_sh1(i), &
                solhrules_l_and_m_of_lm( 1, shrp%lm_sh2(i) ), &
                solhrules_l_and_m_of_lm( 2, shrp%lm_sh2(i) ),   &
                shrp%lm_sh2(i) 
       enddo
    enddo

    write(io_unit,*) 
    write(io_unit,*) '*********************************************************'
    write(io_unit,*)
    write(io_unit,*) "DIFFERENTIAL RULES FOR SOLID HARMONICS UP TO L_MAX", L_MAX
    write(io_unit,*)
    write(io_unit,*) '*********************************************************'
    write(io_unit,*) 
    do lm1 = 1,lm_max
       do lm2 = 1,lm_max
          shrd => solhrules_differential(lm1,lm2)
          if ( shrd%n_summands .gt. 0 ) then
             write(io_unit,*)
             write(io_unit,'("|",I2,",",I2,">/(",I2,")*|",I2,",",I2,">/(",I2,") == ",I2," summands")')  &
                solhrules_l_and_m_of_lm(1,lm1), &
                solhrules_l_and_m_of_lm(2,lm1), &
                lm1,    &
                solhrules_l_and_m_of_lm(1,lm2), &
                solhrules_l_and_m_of_lm(2,lm2), &
                lm2,    &
                shrd%n_summands
             do i=1,shrd%n_summands
                write(io_unit,'(F20.10," |",I2,",",I2,">/(",I2,")")') &
                        shrd%coef(i),   &
                        solhrules_l_and_m_of_lm( 1, shrd%lm_sh(i) ),    &
                        solhrules_l_and_m_of_lm( 2, shrd%lm_sh(i) ),    &
                        shrd%lm_sh(i)
             enddo
          endif
       enddo
    enddo

    write(io_unit,*) 
    write(io_unit,*)

  end subroutine solhrules_print
  !*************************************************************


  !*************************************************************
  subroutine solhrules_lm_of_l_and_m(lm,l,m)
    !  Purpose: mapping of (l,m) to metaindex lm
    implicit none
    !------------ Declaration of formal parameters -------------
    integer(IK), intent(out) :: lm
    integer(IK), intent(in)  :: l,m
    !** End of interface ***************************************
    lm = solhrules_first_lm_of_l(l) + m - 1
  end subroutine solhrules_lm_of_l_and_m
  !*************************************************************


  !*************************************************************
  function prod_rule(f,g,lm_min,lm_max) result(p)
    ! Purpose: this function performs the product rule for two given 
    !          vectors f and g, where the second index is the lm
    !          metaindex and the first the exponent a,b metaindex.
    !          lm_min,lm_max specify the range of lm_result for
    !          which the result should be calculated.
    ! p(lm_result) = sum[ coef  * (sh_grad(lm_sh1) f) * (sh_grad(lm_sh2) g) ]
    integer(IK),intent(in) :: lm_min,lm_max
    real(RK),dimension(:,:),intent(in) :: f,g
    real(RK),dimension(size(f,1),lm_min:lm_max) :: p
    !** End of interface ***************************************  
    integer(IK) :: i,j
    integer(IK),pointer :: index1_p(:),index2_p(:)
    real(RK),pointer :: coef(:)
    p=0.0_RK
    do i=lm_min,lm_max
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       coef=>solhrules_product(i)%coef
       do j=1,solhrules_product(i)%n_summands
          p(:,i)=p(:,i)+coef(j)*f(:,index1_p(j))*g(:,index2_p(j))

       end do
    enddo
  end function prod_rule
  !*************************************************************

#if 0
  function prod_rule2(f,g,l) result(p)
    ! Purpose: this function performs the product rule for two given 
    !          vectors f and g, where the second index is the lm
    !          metaindex and the first the exponent a,b metaindex.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The differenc to routine prod_rule is that the 
    !          summation over  l` is not performed
    ! p(lm_result,l`) = sum(m1,m2)[ coef  * (sh_grad(lm_sh1) f)
    !                                  * (sh_grad(lm_sh2) g) ]
    integer(IK),intent(in) :: l
    real(RK),dimension(:,:),intent(in) :: f,g
    real(RK),dimension(size(f,1),1:(l+1)**2,0:l) :: p
    !** End of interface ***************************************  
    !------------ Declaration of local variables ---------------
    integer(IK) :: i,j,l_i
    integer(IK) :: lm_min,lm_max,i1,i2
    integer(IK), pointer :: index1_p(:),index2_p(:)
    real(RK) :: coeff
    real(RK), pointer :: coeff_p(:)
    !------------ Executable code ------------------------------
    p=0.0_RK
    lm_min=1
    lm_max=(l+1)**2
    do i=lm_min,lm_max
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       coeff_p=>solhrules_product(i)%coef
       do j=1,solhrules_product(i)%n_summands
          l_i=solhrules_l_and_m_of_lm(1,solhrules_product(i)%lm_sh1(j))
          i1=index1_p(j)
          i2=index2_p(j)
          coeff=coeff_p(j)
          p(:,i,l_i)=&
               p(:,i,l_i)&
               +coeff*f(:,i1)*g(:,i2)
       end do
    enddo
  end function prod_rule2

  function prod_rule_nested(f,g,h,l) result(p)
    ! Purpose: this function performs the product rule for three given 
    !          vectors f and g and h, where the second index is the 
    !          lm metaindex and the first the exponent a,b metaindex.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The summation over l` is not performed
    ! p(lm_result,l`) = sum(m1,m2,m3)[ coef  * (sh_grad(lm_sh1) f)
    !                       * (sh_grad(lm_sh2) g) * (sh_grad(lm_sh3) h) ]
    integer(IK),intent(in) :: l
    real(RK),dimension(:,:),intent(in) :: f,g,h
    real(RK),dimension(size(f,1),1:(l+1)**2,(l+1)*(l+2)/2) :: p
    !** End of interface ***************************************  
    integer(IK) :: i,j,index1,l1,l2,jj
    integer(IK) :: lm_min,lm_max
    p=0.0_RK
    lm_min=1
    lm_max=(l+1)**2
    do i=lm_min,lm_max
       do j=1,solhrules_product(i)%n_summands
          index1=solhrules_product(i)%lm_sh1(j)
          l1=solhrules_l_and_m_of_lm(1,index1)
          do jj=1,solhrules_product(index1)%n_summands
             l2=solhrules_l_and_m_of_lm(1,solhrules_product(index1)%lm_sh1(jj))
             p(:,i,l1*(l1+1)/2+l2+1)=&
                  p(:,i,l1*(l1+1)/2+l2+1)&
                  +solhrules_product(i)%coef(j)*&
                  solhrules_product(index1)%coef(jj)*&
                  f(:,solhrules_product(index1)%lm_sh1(jj))*&
                  g(:,solhrules_product(index1)%lm_sh2(jj))*&
                  h(:,solhrules_product(i)%lm_sh2(j))           
          end do! jj
       end do! j
    enddo! i
  end function prod_rule_nested
#endif

  function prod_rule_nested2(f,g,h,l,lb,mb,fact1,fact2) result(p)
    ! Purpose: this function performs the three-fold nested product rule
    !          for three given vectors f, g and h , where the second index
    !          is the lm metaindex and the first the a,b exponent metaindex
    !          For g and h the third index is an lm metaindex lmb = (lb**2)+mb
    !          that specifies the range of results included.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The summation over lmb is not performed and the result is
    !          indexed with lb, intended for summing up.
    ! p(lm_result,lb) = sum(m1,m2,m3)[ coef  * (sh_grad(lm_sh1) f)
    !             * (sh_grad(lm_sh2) g(lmb)) * (sh_grad(lm_sh3) h(lmb)) ]
    integer(IK),intent(in) :: l,lb,mb
    real(RK),dimension(:,:),intent(in) :: f
    real(RK),dimension(:,:,:),intent(in) :: g
    real(RK),dimension(:,:,:),intent(in) :: h
    real(RK),dimension(:),intent(in) :: fact1,fact2
    real(RK),dimension(size(f,1),1:(l+1)**2,0:lb+l) :: p
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    real(RK) :: help_vec(size(f,1))
    integer(IK) :: i,j_1,j_2,j_3,l1,l3,lm_min,lm_max,lm_in,n_sum1,n_sum2,n_sum3, &
         l_j1,l_j3
    integer(IK),pointer,dimension(:) :: index_p,index0_p,index1_p,&
         index2_p,index3_p,index3_p2
    real(RK),pointer,dimension(:) :: coef1_p,coef2_p,coef3_p
    !------------ Executable code ------------------------------
    p=0.0_RK
    lm_in=(lb**2)+mb
    lm_min=1
    lm_max=(l+1)**2
    index_p=>solhrules_product(lm_in)%lm_sh1
    index0_p=>solhrules_product(lm_in)%lm_sh2
    coef1_p=>solhrules_product(lm_in)%coef
    n_sum1=solhrules_product(lm_in)%n_summands
    do i=lm_min,lm_max
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       n_sum2=solhrules_product(i)%n_summands
       coef2_p=>solhrules_product(i)%coef
       do j_1=1,n_sum1
          l1=solhrules_l_and_m_of_lm(1,index_p(j_1))
          help_vec=fact2**l1
          do j_2=1,n_sum2
             index3_p=>solhrules_product(index1_p(j_2))%lm_sh1
             index3_p2=>solhrules_product(index1_p(j_2))%lm_sh2
             n_sum3=solhrules_product(index1_p(j_2))%n_summands
             coef3_p=>solhrules_product(index1_p(j_2))%coef
             do j_3=1,n_sum3
                l3=solhrules_l_and_m_of_lm(1,index3_p(j_3))
                l_j3=solhrules_l_and_m_of_lm(1,index3_p2(j_3))
                l_j1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                if(l_j3<=l_j1) then
                   p(:,i,l1+l3)=&
                        p(:,i,l1+l3)&
                        +coef2_p(j_2)*coef3_p(j_3)*coef1_p(j_1)*&
                        f(:,index3_p(j_3))*g(:,index3_p2(j_3),index_p(j_1))*&
                        h(:,index2_p(j_2),index0_p(j_1))*fact1**l3*help_vec
                end if
             end do! j_3
          end do! j_2
       enddo! j_1
    end do! i=1
  end function prod_rule_nested2

#if 0
  function nsum_prod_rule_nested2(l,lb,mb) result(p)
    ! Purpose: this function performs the three-fold nested product rule
    !          for three given vectors f, g and h , where the second index
    !          is the lm metaindex and the first the a,b exponent metaindex
    !          For g and h the third index is an lm metaindex lmb = (lb**2)+mb
    !          that specifies the range of results included.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The summation over lmb is not performed and the result is
    !          indexed with lb, intended for summing up.
    ! p(lm_result,lb) = sum(m1,m2,m3)[ coef  * (sh_grad(lm_sh1) f)
    !             * (sh_grad(lm_sh2) g(lmb)) * (sh_grad(lm_sh3) h(lmb)) ]
    integer(IK),intent(in) :: l,lb,mb
!    real(RK),dimension(:,:),intent(in) :: f
!    real(RK),dimension(:,:,:),intent(in) :: g
!    real(RK),dimension(:,:,:),intent(in) :: h
!    real(RK),dimension(:),intent(in) :: fact1,fact2
!    real(RK),dimension(size(f,1),1:(l+1)**2,0:lb+l) :: p
    integer(IK) :: p
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(IK) :: i,j_1,j_2,j_3,l1,l3,lm_min,lm_max,lm_in,n_sum1,n_sum2,n_sum3, &
         l_j1,l_j3
    integer(IK),pointer,dimension(:) :: index_p,index0_p,index1_p,&
         index2_p,index3_p,index3_p2
    real(RK),pointer,dimension(:) :: coef1_p,coef2_p,coef3_p
    !------------ Executable code ------------------------------
    p=0
    lm_in=(lb**2)+mb
    lm_min=1
    lm_max=(l+1)**2
    index_p=>solhrules_product(lm_in)%lm_sh1
    index0_p=>solhrules_product(lm_in)%lm_sh2
    coef1_p=>solhrules_product(lm_in)%coef
    n_sum1=solhrules_product(lm_in)%n_summands
    do i=lm_min,lm_max
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       n_sum2=solhrules_product(i)%n_summands
       coef2_p=>solhrules_product(i)%coef
       do j_1=1,n_sum1
          l1=solhrules_l_and_m_of_lm(1,index_p(j_1))
          do j_2=1,n_sum2
             index3_p=>solhrules_product(index1_p(j_2))%lm_sh1
             index3_p2=>solhrules_product(index1_p(j_2))%lm_sh2
             n_sum3=solhrules_product(index1_p(j_2))%n_summands
             coef3_p=>solhrules_product(index1_p(j_2))%coef
             do j_3=1,n_sum3
                l3=solhrules_l_and_m_of_lm(1,index3_p(j_3))
                l_j3=solhrules_l_and_m_of_lm(1,index3_p2(j_3))
                l_j1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                if(l_j3<=l_j1) then
                p=p+1
                end if
             end do! j_3
          end do! j_2
       enddo! j_1
    end do! i=1
  end function nsum_prod_rule_nested2

  function nsums_prod_rule_nested2(l,lb,mb,nne_summands) result(p)
    ! Purpose: this function performs the three-fold nested product rule
    !          for three given vectors f, g and h , where the second index
    !          is the lm metaindex and the first the a,b exponent metaindex
    !          For g and h the third index is an lm metaindex lmb = (lb**2)+mb
    !          that specifies the range of results included.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The summation over lmb is not performed and the result is
    !          indexed with lb, intended for summing up.
    ! p(lm_result,lb) = sum(m1,m2,m3)[ coef  * (sh_grad(lm_sh1) f)
    !             * (sh_grad(lm_sh2) g(lmb)) * (sh_grad(lm_sh3) h(lmb)) ]
    integer(IK),intent(in) :: l,lb,mb
    integer(IK),optional,intent(out)::nne_summands
    integer(IK) :: p
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(IK) :: i,j_1,j_2,j_3,l1,l3,lm_min,lm_max,lm_in,n_sum1,n_sum2,n_sum3, &
         l_j1,l_j3
    integer(IK),pointer,dimension(:) :: index_p,index0_p,index1_p,&
         index2_p,index3_p,index3_p2
    real(RK),pointer,dimension(:) :: coef1_p,coef2_p,coef3_p
    !------------ Executable code ------------------------------
    p=0
    if(present(nne_summands)) nne_summands=0
    lm_in=(lb**2)+mb
    lm_min=1
    lm_max=(l+1)**2
    index_p=>solhrules_product(lm_in)%lm_sh1
    index0_p=>solhrules_product(lm_in)%lm_sh2
    coef1_p=>solhrules_product(lm_in)%coef
    n_sum1=solhrules_product(lm_in)%n_summands
    do i=lm_min,lm_max
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       n_sum2=solhrules_product(i)%n_summands
       coef2_p=>solhrules_product(i)%coef
       do j_1=1,n_sum1
          l1=solhrules_l_and_m_of_lm(1,index_p(j_1))
          do j_2=1,n_sum2
             index3_p=>solhrules_product(index1_p(j_2))%lm_sh1
             index3_p2=>solhrules_product(index1_p(j_2))%lm_sh2
             n_sum3=solhrules_product(index1_p(j_2))%n_summands
             coef3_p=>solhrules_product(index1_p(j_2))%coef
             do j_3=1,n_sum3
                l3=solhrules_l_and_m_of_lm(1,index3_p(j_3))
                l_j3=solhrules_l_and_m_of_lm(1,index3_p2(j_3))
                l_j1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                if(l_j3<=l_j1) p=p+1
                if(present(nne_summands).and.l_j3<l_j1) nne_summands=nne_summands+1
             end do! j_3
          end do! j_2
       enddo! j_1
    end do! i=1
  end function nsums_prod_rule_nested2
    subroutine gsummands_prod_rule_nested2(l,lb,mb,nested2,l1_max,l3_max,nested2ne)
    ! Purpose: this function performs the three-fold nested product rule
    !          for three given vectors f, g and h , where the second index
    !          is the lm metaindex and the first the a,b exponent metaindex
    !          For g and h the third index is an lm metaindex lmb = (lb**2)+mb
    !          that specifies the range of results included.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The summation over lmb is not performed and the result is
    !          indexed with lb, intended for summing up.
    ! p(lm_result,lb) = sum(m1,m2,m3)[ coef  * (sh_grad(lm_sh1) f)
    !             * (sh_grad(lm_sh2) g(lmb)) * (sh_grad(lm_sh3) h(lmb)) ]
    integer(IK),intent(in) :: l,lb,mb
    integer(IK),intent(out)::l1_max,l3_max
!    integer(IK),intent(inout),optional:: nne_summands
    type(nested2_vars),intent(out),dimension(:) :: nested2ne
    type(nested2_vars),intent(out),dimension(:) :: nested2
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(IK) :: i,j_1,j_2,j_3,l1,l3,lm_min,lm_max,lm_in,n_sum1,n_sum2,n_sum3, &
         l_j1,l_j3
    integer(IK),pointer,dimension(:) :: index_p,index0_p,index1_p,&
         index2_p,index3_p,index3_p2
    real(RK),pointer,dimension(:) :: coef1_p,coef2_p,coef3_p
    !------------ Executable code ------------------------------
    integer(IK) ::p,p_ne
    logical:: ne_summands
    p=0

    p_ne=0
    ne_summands=size(nested2ne).ne.0

    l1_max=0
    l3_max=0
    lm_in=(lb**2)+mb
    lm_min=1
    lm_max=(l+1)**2
    index_p=>solhrules_product(lm_in)%lm_sh1
    index0_p=>solhrules_product(lm_in)%lm_sh2
    coef1_p=>solhrules_product(lm_in)%coef
    n_sum1=solhrules_product(lm_in)%n_summands
    do i=lm_min,lm_max
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       n_sum2=solhrules_product(i)%n_summands
       coef2_p=>solhrules_product(i)%coef
       do j_1=1,n_sum1
          l1=solhrules_l_and_m_of_lm(1,index_p(j_1))
          do j_2=1,n_sum2
             index3_p=>solhrules_product(index1_p(j_2))%lm_sh1
             index3_p2=>solhrules_product(index1_p(j_2))%lm_sh2
             n_sum3=solhrules_product(index1_p(j_2))%n_summands
             coef3_p=>solhrules_product(index1_p(j_2))%coef
             do j_3=1,n_sum3
                l3=solhrules_l_and_m_of_lm(1,index3_p(j_3))
                l_j3=solhrules_l_and_m_of_lm(1,index3_p2(j_3))
                l_j1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                if(l_j3<=l_j1) then
                if(l1.gt.l1_max) l1_max=l1
                if(l3.gt.l3_max) l3_max=l3
                p=p+1
                nested2(p)%i=i
                nested2(p)%coef_prod=coef2_p(j_2)*coef3_p(j_3)*coef1_p(j_1)
                nested2(p)%index3_p=index3_p(j_3)
                nested2(p)%index3_p2=index3_p2(j_3)
                nested2(p)%index_p=index_p(j_1)
                nested2(p)%index2_p=index2_p(j_2)
                nested2(p)%index0_p=index0_p(j_1)
                nested2(p)%l3=l3
                nested2(p)%l1=l1
                 if(ne_summands.and.l_j3.lt.l_j1) then
                  p_ne=p_ne+1
                  nested2ne(p_ne)=nested2(p)
                 endif
                end if
             end do! j_3
          end do! j_2
       enddo! j_1
    end do! i=1
  end subroutine gsummands_prod_rule_nested2

    subroutine summands_prod_rule_nested2(l,lb,mb,nested2,l1_max,l3_max)
    ! Purpose: this function performs the three-fold nested product rule
    !          for three given vectors f, g and h , where the second index
    !          is the lm metaindex and the first the a,b exponent metaindex
    !          For g and h the third index is an lm metaindex lmb = (lb**2)+mb
    !          that specifies the range of results included.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The summation over lmb is not performed and the result is
    !          indexed with lb, intended for summing up.
    ! p(lm_result,lb) = sum(m1,m2,m3)[ coef  * (sh_grad(lm_sh1) f)
    !             * (sh_grad(lm_sh2) g(lmb)) * (sh_grad(lm_sh3) h(lmb)) ]
    integer(IK),intent(in) :: l,lb,mb
    integer(IK),intent(out)::l1_max,l3_max
    type(nested2_vars),intent(out),dimension(:) :: nested2
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(IK) :: i,j_1,j_2,j_3,l1,l3,lm_min,lm_max,lm_in,n_sum1,n_sum2,n_sum3, &
         l_j1,l_j3
    integer(IK),pointer,dimension(:) :: index_p,index0_p,index1_p,&
         index2_p,index3_p,index3_p2
    real(RK),pointer,dimension(:) :: coef1_p,coef2_p,coef3_p
    !------------ Executable code ------------------------------
    integer(IK) ::p
    p=0

    l1_max=0
    l3_max=0
    lm_in=(lb**2)+mb
    lm_min=1
    lm_max=(l+1)**2
    index_p=>solhrules_product(lm_in)%lm_sh1
    index0_p=>solhrules_product(lm_in)%lm_sh2
    coef1_p=>solhrules_product(lm_in)%coef
    n_sum1=solhrules_product(lm_in)%n_summands
    do i=lm_min,lm_max
       index1_p=>solhrules_product(i)%lm_sh1
       index2_p=>solhrules_product(i)%lm_sh2
       n_sum2=solhrules_product(i)%n_summands
       coef2_p=>solhrules_product(i)%coef
       do j_1=1,n_sum1
          l1=solhrules_l_and_m_of_lm(1,index_p(j_1))
          do j_2=1,n_sum2
             index3_p=>solhrules_product(index1_p(j_2))%lm_sh1
             index3_p2=>solhrules_product(index1_p(j_2))%lm_sh2
             n_sum3=solhrules_product(index1_p(j_2))%n_summands
             coef3_p=>solhrules_product(index1_p(j_2))%coef
             do j_3=1,n_sum3
                l3=solhrules_l_and_m_of_lm(1,index3_p(j_3))
                l_j3=solhrules_l_and_m_of_lm(1,index3_p2(j_3))
                l_j1=solhrules_l_and_m_of_lm(1,index_p(j_1))
                if(l_j3<=l_j1) then
                if(l1.gt.l1_max) l1_max=l1
                if(l3.gt.l3_max) l3_max=l3
                p=p+1
                nested2(p)%i=i
                nested2(p)%coef_prod=coef2_p(j_2)*coef3_p(j_3)*coef1_p(j_1)
                nested2(p)%index3_p=index3_p(j_3)
                nested2(p)%index3_p2=index3_p2(j_3)
                nested2(p)%index_p=index_p(j_1)
                nested2(p)%index2_p=index2_p(j_2)
                nested2(p)%index0_p=index0_p(j_1)
                nested2(p)%l3=l3
                nested2(p)%l1=l1
                end if
             end do! j_3
          end do! j_2
       enddo! j_1
    end do! i=1
  end subroutine summands_prod_rule_nested2

  function opt_prod_rule_nested2(f,g,h,l,lb,l1_max,l3_max,nested2,fact12) result(p)
    ! Purpose: this function performs the three-fold nested product rule
    !          for three given vectors f, g and h , where the second index
    !          is the lm metaindex and the first the a,b exponent metaindex
    !          For g and h the third index is an lm metaindex lmb = (lb**2)+mb
    !          that specifies the range of results included.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The summation over lmb is not performed and the result is
    !          indexed with lb, intended for summing up.
    ! p(lm_result,lb) = sum(m1,m2,m3)[ coef  * (sh_grad(lm_sh1) f)
    !             * (sh_grad(lm_sh2) g(lmb)) * (sh_grad(lm_sh3) h(lmb)) ]
    integer(IK),intent(in) :: l,lb,l1_max,l3_max
    real(RK),dimension(:,:),intent(in) :: f
    real(RK),dimension(:,:,:),intent(in) :: g
    real(RK),dimension(:,:,:),intent(in) :: h

!    real(RK),dimension(size(f,1),0:l3_max), intent(in) :: fact1
!    real(RK),dimension(size(f,1),0:l1_max), intent(in) :: fact2

    real(RK),dimension(size(f,1),0:l3_max,0:l1_max), intent(in) :: fact12
    
    real(RK),dimension(size(f,1),1:(l+1)**2,0:lb+l) :: p
    type(nested2_vars),intent(in),dimension(:) :: nested2
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
!!$    real(RK) :: help_vec(size(f,1))
    integer(IK) :: i
    !------------ Executable code ------------------------------
    p=0.0_RK
        do i=1,size(nested2)
                   p(:,nested2(i)%i,nested2(i)%l1+nested2(i)%l3)=&
                        p(:,nested2(i)%i,nested2(i)%l1+nested2(i)%l3)&
                        +nested2(i)%coef_prod*&
                        f(:,nested2(i)%index3_p)*g(:,nested2(i)%index3_p2,nested2(i)%index_p)*&
                        h(:,nested2(i)%index2_p,nested2(i)%index0_p)*fact12(:,nested2(i)%l3,nested2(i)%l1)
!                        h(:,nested2(i)%index2_p,nested2(i)%index0_p)*fact1(:,nested2(i)%l3)*fact2(:,nested2(i)%l1)
    end do! i=1
  end function opt_prod_rule_nested2

  function prod_rule_nested3(f,g,l,lb,mb,l_i_max,l_c) result(p)
    ! Purpose: this function performs the product rule for two given 
    !          vectors f and g, where the first index is the a,b exponent 
    !          metaindex, the secound is lm metaindex and the third
    !          index is an lm metaindex lmb = (lb**2)+mb that specifies
    !          the range of results included.
    !          l specifies the angular momentum up to wich results
    !          have to be calculated: lm_result = 1, (l+1)**2.
    !          The summation over l` is not performed
    ! p(lm_result,lb) = sum(m1,m2)[ coef  * (sh_grad(lm_sh1) f(lmb))
    !                       * (sh_grad(lm_sh2) g(lmb)) ]
    integer(IK),intent(in) :: l,lb,mb,l_i_max,l_c
    real(RK),dimension(:,:,:),intent(in) :: f
    real(RK),dimension(:,:,:,0:),intent(in) :: g
    real(RK),dimension(size(f,1),l**2+1:(l+1)**2,0:l_i_max) :: p
    !** End of interface ***************************************  
    !------------ Declaration of local variables ---------------
    integer(IK) :: i,j_1,j_2,lm_min,lm_max,lm_in,n_sum1,n_sum2, &
        l_i,l1,l2 
    integer(IK),pointer,dimension(:) :: index1_1p,&
         index1_2p,index2_1p,index2_2p
    real(RK),pointer :: coeff_p1(:),coeff_p2(:)
    real(RK) :: help_vec(size(f,1))
    !------------ Executable code ------------------------------
    p=0.0_RK
    lm_in=(lb**2)+mb
    lm_min=l**2+1
    lm_max=(l+1)**2
    index1_1p=>solhrules_product(lm_in)%lm_sh1
    index1_2p=>solhrules_product(lm_in)%lm_sh2
    coeff_p1=>solhrules_product(lm_in)%coef
    n_sum1=solhrules_product(lm_in)%n_summands
    do i=lm_min,lm_max
       n_sum2=solhrules_product(i)%n_summands
       index2_1p=>solhrules_product(i)%lm_sh1
       index2_2p=>solhrules_product(i)%lm_sh2
       coeff_p2=>solhrules_product(i)%coef
       do j_1=1,n_sum1
          do j_2=1,n_sum2
             l1=solhrules_l_and_m_of_lm(1,index2_1p(j_2))
             l2=solhrules_l_and_m_of_lm(1,index1_1p(j_1))
             if(l1+l2<=l_c) then
                help_vec(:)=coeff_p2(j_2)*coeff_p1(j_1)*&
                     f(:,index2_1p(j_2),index1_1p(j_1))
                do l_i=0,l_i_max
                   p(:,i,l_i)=p(:,i,l_i)+help_vec(:)*&
                        g(:,index2_2p(j_2),index1_2p(j_1),l_i)
                end do
             end if
          end do! j_2
       enddo! j_1
    end do! i
  end function prod_rule_nested3
#endif

  !*************************************************************
  function diff_rule_sclr(f,lm_min1,lm_max1,lm_2) result(p)
    ! Purpose: Scalar version of diff_rule. 
    !          This function performs the differentialrule for a given f,
    !          where the second index is the lm (derivative) metaindex
    !          lm_min1,lm_max1 specify the range of lm_grad
    !          for which results have to be calculated.
    ! p(lm_grad) = sh_grad(lm_grad) f(lm_2) = sum[ coef * f(lm_sh) ]
    integer(IK),intent(in)             :: lm_min1,lm_max1,lm_2
    real(RK),dimension(:),intent(in)   :: f
    real(RK),dimension(lm_min1:lm_max1):: p
    ! result
    !** End of interface ***************************************
    integer(IK) :: i,k,n_sum
    integer(IK),pointer ::  index_p(:)
    real(RK),pointer :: coef(:)
    p=0.0_RK
    do i=lm_min1,lm_max1
       n_sum=solhrules_differential(i,lm_2)%n_summands
       coef=>solhrules_differential(i,lm_2)%coef
       index_p=>solhrules_differential(i,lm_2)%lm_sh
       do k=1,n_sum
          p(i)=p(i)+coef(k)*f(index_p(k))
       end do
    enddo
  end function diff_rule_sclr

  function diff_rule_vec(f,lm_min1,lm_max1,lm_2) result(p)
    ! Purpose: this function performs the differentialrule for a given vector f,
    !          where the second index is the lm metaindex and the first the a,b 
    !          exponent metaindex. lm_min1,lm_max1 specify the range of lm_grad
    !          for which results have to be calculated.
    ! p(lm_grad) = sh_grad(lm_grad) f(lm_2) = sum[ coef * f(lm_sh) ]
    integer(IK),intent(in) :: lm_min1,lm_max1,lm_2
    real(RK),dimension(:,:),intent(in) :: f
    real(RK),dimension(size(f,1),lm_min1:lm_max1) :: p
    ! result
    !** End of interface ***************************************
    integer(IK) :: i,k,n_sum
    integer(IK),pointer ::  index_p(:)
    real(RK),pointer :: coef(:)
    p=0.0_RK
    do i=lm_min1,lm_max1
       n_sum=solhrules_differential(i,lm_2)%n_summands
       coef=>solhrules_differential(i,lm_2)%coef
       index_p=>solhrules_differential(i,lm_2)%lm_sh
       do k=1,n_sum
          p(:,i)=p(:,i)+coef(k)*f(:,index_p(k))
       end do
    enddo
  end function diff_rule_vec


  function diff_rule_nested(f,lm_min1,lm_max1,lm_2,lm_3) result(p)
    ! Purpose: this function performs the two fold nested 
    !          differentialrul on a given vector f,where the 
    !          second index is the lm metaindex and the first the a,b 
    !          exponent metaindex. lm_min1,lm_max1 specify the ranges 
    !          for which results have to be calculated.
    ! p(lm_grad) = sh_grad(lm_2) (sh_grad(lm_grad) f(lm_3))
    integer(IK),intent(in) :: lm_min1,lm_max1,lm_2,lm_3
    real(RK),dimension(:,:),intent(in) :: f
    real(RK),dimension(size(f,1),lm_min1:lm_max1) :: p
    ! result
    !** End of interface ***************************************
    integer(IK) :: i,k_1,k_2,index1,index2,n_sum1
    integer(IK),pointer :: index1_p(:),index2_p(:)
    real(RK),pointer :: coef1(:),coef2(:)
    p=0.0_RK
    n_sum1=solhrules_differential(lm_2,lm_3)%n_summands
    coef1=>solhrules_differential(lm_2,lm_3)%coef
    index1_p=>solhrules_differential(lm_2,lm_3)%lm_sh
    do i=lm_min1,lm_max1
       do k_1=1,n_sum1
          index1=index1_p(k_1)
          coef2=>solhrules_differential(i,index1)%coef
          index2_p=>solhrules_differential(i,index1)%lm_sh
          do k_2=1,solhrules_differential(i,index1)%n_summands
             index2=index2_p(k_2)
             p(:,i)=p(:,i)+coef2(k_2)*coef1(k_1)*f(:,index2)
          end do
       end do
    enddo
  end function diff_rule_nested


  function prod_rule_double(f,g,lm_min,lm_max,jn_min,jn_max) result(p)
    ! Purpose: this function performs the product rule simultaneoulsly 
    !          for the second and third coefficents of two given 
    !          vectors f and g, where the second lm and third jn index are
    !          the lm metaindices and the first the a,b exponent metaindex.
    !          lm_min,lm_max and jn_min,jn_max specify the range of lm_r1
    !          and lm_r2 for which the result p(lm_r1,lm_r2) should be calculated.
    integer(IK),intent(in) :: lm_min,lm_max,jn_min,jn_max
    real(RK),dimension(:,:,:),intent(in) :: f,g
    real(RK),dimension(size(f,1),lm_min:lm_max,jn_min:jn_max) :: p
    !** End of interface ***************************************  
    !------------ Declaration of local variables ---------------
    integer(IK) :: i_lm,j_lm,i_jn,j_jn,i_1_lm,i_1_jn,i_2_lm,i_2_jn
    real(RK) :: coef_lm,coef_jn,coef
    !------------ Executable code ------------------------------
    p=0.0_RK
    do i_lm=lm_min,lm_max
       do j_lm=1,solhrules_product(i_lm)%n_summands
          i_1_lm=solhrules_product(i_lm)%lm_sh1(j_lm)
          i_2_lm=solhrules_product(i_lm)%lm_sh2(j_lm)   
          coef_lm=solhrules_product(i_lm)%coef(j_lm)
          do i_jn=jn_min,jn_max
             do j_jn=1,solhrules_product(i_jn)%n_summands
                i_1_jn=solhrules_product(i_jn)%lm_sh1(j_jn)
                i_2_jn=solhrules_product(i_jn)%lm_sh2(j_jn) 
                coef_jn=solhrules_product(i_jn)%coef(j_jn)
                coef=coef_jn*coef_lm
                p(:,i_lm,i_jn)=p(:,i_lm,i_jn)+coef*&
                     f(:,i_1_lm,i_1_jn)*&
                     g(:,i_2_lm,i_2_jn)
             end do
          end do
       end do
    enddo
  end function prod_rule_double

#if 0
  function prod_rule_double3(f1,f2,f3,g1,g2,g3,lm_min,lm_max,jn_min,jn_max) result(p)
    ! Purpose: this function performs the product rule simultaneoulsly 
    !          for the second and third coefficents of three sets of two given 
    !          vectors fi and gi, where the second lm and third jn index are
    !          the lm metaindices and the first the a,b exponent metaindex.
    !          lm_min,lm_max and jn_min,jn_max specify the range of lm_r1
    !          and lm_r2 for which the result p(lm_r1,lm_r2) should be calculated.
    !          The result is the sum of the applications of the double product
    !          rule to the three pairs of fi and gi:
    !    p(lm_r1,lm_r2) = sum(i=1,3) [ prod_rule_double(fi,gi,...) ]
    integer(IK),intent(in) :: lm_min,lm_max,jn_min,jn_max
    real(RK),dimension(:,:,:),intent(in) :: f1,f2,f3,g1,g2,g3
    real(RK),dimension(size(f1,1),lm_min:lm_max,jn_min:jn_max) :: p
    !** End of interface ***************************************  
    !------------ Declaration of local variables ---------------
    integer(IK) :: i_lm,j_lm,i_jn,j_jn,i_1_lm,i_1_jn,i_2_lm,i_2_jn
    real(RK) :: coef_lm,coef_jn,coef
    !------------ Executable code ------------------------------
    p=0.0_RK
    do i_lm=lm_min,lm_max
       do j_lm=1,solhrules_product(i_lm)%n_summands
          i_1_lm=solhrules_product(i_lm)%lm_sh1(j_lm)
          i_2_lm=solhrules_product(i_lm)%lm_sh2(j_lm)
          coef_lm=solhrules_product(i_lm)%coef(j_lm)
          do i_jn=jn_min,jn_max
             do j_jn=1,solhrules_product(i_jn)%n_summands
                i_1_jn=solhrules_product(i_jn)%lm_sh1(j_jn)
                i_2_jn=solhrules_product(i_jn)%lm_sh2(j_jn)
                coef_jn=solhrules_product(i_jn)%coef(j_jn)
                coef=coef_jn*coef_lm
                p(:,i_lm,i_jn)=p(:,i_lm,i_jn)+coef*(&
                     f1(:,i_1_lm,i_1_jn)*&
                     g1(:,i_2_lm,i_2_jn)+&
                     f2(:,i_1_lm,i_1_jn)*&
                     g2(:,i_2_lm,i_2_jn)+&
                     f3(:,i_1_lm,i_1_jn)*&
                     g3(:,i_2_lm,i_2_jn))
             end do
          end do
       end do
    enddo
  end function prod_rule_double3

  function prod_rule3(f,g,l,upto_l) result(p)
    ! Purpose: this function performs the product rule for two given 
    !          vectors f and g, where the second index is the lm
    !          metaindex and the first the exponent a,b metaindex.
    !          l specifies the angular momentum up to wich results
    !          have lo be calculated: lm_result = 1, (l+1)**2.
    !          The difference to routine prod_rule is that the 
    !          summation over l` is not performed.
    !          The difference to prod_rule2 is that upto_l specifies the
    !          the range of results to be stored: l` = 0, upto_l
    ! p(lm_result,l`) = sum(m1,m2)[ coef  * (sh_grad(lm_sh1) f)
    !                                  * (sh_grad(lm_sh2) g) ]
    integer(IK),intent(in) :: l, upto_l
    real(RK),dimension(:,:),intent(in) :: f,g
    real(RK),dimension(size(f,1),1:(l+1)**2,0:upto_l) :: p
    !** End of interface ***************************************  
    !------------ Declaration of local variables ---------------
    integer(IK) :: i,j,l_i,i1,i2
    integer(IK) :: lm_min,lm_max
    real(RK) :: coef
    !------------ Executable code ------------------------------
    p=0.0_RK
    lm_min=1
    lm_max=(l+1)**2
    do i=lm_min,lm_max
       do j=1,solhrules_product(i)%n_summands
          i1=solhrules_product(i)%lm_sh1(j)           
          l_i=solhrules_l_and_m_of_lm(1,i1)
          if (l_i<=upto_l) then
             i2=solhrules_product(i)%lm_sh2(j)
             coef=solhrules_product(i)%coef(j)
             p(:,i,l_i)=&
                  p(:,i,l_i)+coef*&
                  f(:,i1)*&
                  g(:,i2)
          end if
       end do
    enddo
  end function prod_rule3
#endif

#if 0
function prod_rule_scalar(f,g,l,upto_l) result(p)
  ! Purpose: this function performs the product rule for two given 
  !          vectors f and g, where the second index is the lm metaindex and
  !          the first the a,b metaindex.
  !          l specifies the angular momentum
  !          The differenc to routine prod_rule is that the 
  !          summation over 
  !          l` is only performed up to upto_l
integer(IK),intent(in) :: l,upto_l
real(RK),dimension(:,:),intent(in) :: f,g
real(RK),dimension(size(f,1),1:(l+1)**2) :: p
!** End of interface ***************************************  
integer(IK) :: i,j,l_i
integer(IK) :: lm_min,lm_max
p=0.0_RK
lm_min=1
lm_max=(l+1)**2
do i=lm_min,lm_max
 do j=1,solhrules_product(i)%n_summands
    l_i=solhrules_l_and_m_of_lm(1,solhrules_product(i)%lm_sh1(j))
    if (l_i.gt.upto_l) then
       exit
    endif
    !if (i.eq.1) then
    !   write(*,*) "i=1"
    !   write(*,*) "solhrules_product(i)%coef(j)",solhrules_product(i)%coef(j)
    !   write(*,*) "f(1,solhrules_product(i)%lm_sh1(j))",f(1,solhrules_product(i)%lm_sh1(j))
    !   write(*,*) "g(1,solhrules_product(i)%lm_sh2(j))",g(1,solhrules_product(i)%lm_sh2(j))
    !endif
       p(:,i)=&
            p(:,i)&
            +solhrules_product(i)%coef(j)*&
            f(:,solhrules_product(i)%lm_sh1(j))*&
            g(:,solhrules_product(i)%lm_sh2(j))
 end do
enddo
end function prod_rule_scalar

function prod_rule_cross(f,g,l,upto_l) result(p)
  ! Purpose: this function performs the product rule for the cross product of two given 
  !          vectors f and g, where the second index is the lm metaindex and
  !          the first the a,b metaindex.
  !          the third index is the spatial vector
  !          l specifies the angular momentum

  integer(IK),intent(in) :: l,upto_l
real(RK),dimension(:,:,:),intent(in) :: f,g
real(RK),dimension(size(f,1),1:2*l+1,3) :: p
!** End of interface ***************************************  
integer(IK) :: i,j,l_i,i_m
integer(IK) :: lm_min,lm_max
p=0.0_RK
lm_min=l**2 + 1
lm_max=(l+1)**2
do i=lm_min,lm_max
 i_m = i - l**2
 do j=1,solhrules_product(i)%n_summands
    l_i=solhrules_l_and_m_of_lm(1,solhrules_product(i)%lm_sh1(j))
    if (l_i.gt.upto_l) then
       exit
    end if
    p(:,i_m,1)=&
         p(:,i_m,1)&
         +solhrules_product(i)%coef(j)*&
         (f(:,solhrules_product(i)%lm_sh1(j),2)*g(:,solhrules_product(i)%lm_sh2(j),3)&
         -f(:,solhrules_product(i)%lm_sh1(j),3)*g(:,solhrules_product(i)%lm_sh2(j),2))
    p(:,i_m,2)=&
         p(:,i_m,2)&
         +solhrules_product(i)%coef(j)*&
         (f(:,solhrules_product(i)%lm_sh1(j),3)*g(:,solhrules_product(i)%lm_sh2(j),1)&
         -f(:,solhrules_product(i)%lm_sh1(j),1)*g(:,solhrules_product(i)%lm_sh2(j),3))
    p(:,i_m,3)=&
         p(:,i_m,3)&
         +solhrules_product(i)%coef(j)*&
         (f(:,solhrules_product(i)%lm_sh1(j),1)*g(:,solhrules_product(i)%lm_sh2(j),2)&
         -f(:,solhrules_product(i)%lm_sh1(j),2)*g(:,solhrules_product(i)%lm_sh2(j),1))
 end do
enddo
end function prod_rule_cross
#endif
!*************************************************************

end module solhrules_module
