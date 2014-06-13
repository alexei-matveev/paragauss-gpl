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
module  valence_coord_module
  !---------------------------------------------------------------
  !
  !  Purpose: routines to establish the (redundant) valence
  !           coordinates of the molecule. Valence coordinates
  !           are considered to be:
  !           - bond lengths based on a table of covalent radii
  !           - bond angles and
  !           - dihedral angles based on the above bond lengths
  !           - out-of-plane deformations ??
  !
  !  Module called by: ...
  !
  !  References: Program 'NABOR' by T.A. Halgren
  !
  !  Author: ...
  !  Date: ...
  !
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
# include "def.h"
  use type_module ! type specification parameters
  use atom_data_module
  use math_module
  use coortype_module
  use opt_data_module !!!!!!!!!!!!
  use coordinates_module, only: set_internal,val_types,&
                                generate_bmat,print_internal
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================
  !------------ Declaration of constants and variables ------------

  !------------ public functions and subroutines ------------------
  public :: valence_setup
  public :: set_hesse
  !================================================================
  ! End of public interface of module
  !================================================================

  integer(kind=i4_kind), private :: n_valence ! number of valence coords
  
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine valence_setup()
    !  Purpose: calculate the maximum number of possible valence
    !           coordinates according to
    !           binomial_coeff(n,2) possible bond lengths
    !           3*binomial_coeff(n,3) possible bond angles
    !           6*binomial_coeff(n,4) possible dihedral angles
    !
    !           where n is the number of atoms
    USE_DEBUG
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: total,alloc_stat
    integer(kind=i4_kind),parameter :: n_col=10
    real(kind=r8_kind),allocatable  :: hesse_diagonal(:)
    real(kind=r8_kind),allocatable  :: hesse_cartesian(:,:)
    real(kind=r8_kind),allocatable  :: hesse_nodummies(:,:)
    real(kind=r8_kind),allocatable  :: bmat_valence(:,:)
    type(int_coor), allocatable     :: v(:) ! (n_valence)
    integer(i4_kind)                :: numa
    !------------ Executable code --------------------------------

    DPRINT 'valence_setup: entered'

    numa = n_atoms+n_dummy

    total=binomial_coeff(numa,2) + &
         3*binomial_coeff(numa,3) + &
         6* binomial_coeff(numa,4)

    write(OPT_STDOUT,*)" valence_setup: the maximal number of possible valence coords is"
    write(OPT_STDOUT,*)total

    call val_types(n_coor=n_valence,deloc=.false.,calc_dimension=.true.,keep_coor=keep_valence)

    write(OPT_STDOUT,*)"                the actual number of valence coordinates is"
    write(OPT_STDOUT,*)n_valence
    write(OPT_STDOUT,*)" "


    allocate(v(n_valence),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         (" valence_setup: allocation (1) failed")
    v(:)%value=zero
    v(:)%var = .true.

    call val_types(coor=v,n_coor=n_valence,deloc=.false.,calc_dimension=.false.,keep_coor=keep_valence)
    call set_internal(v(:)%value,n_valence,v,atom)
    write(OPT_STDOUT,*)"        Valence Coordinates "
    call print_internal(v,n_valence,n_val=n_valence,full=.true.)
    write(OPT_STDOUT,*)
    
    allocate(bmat_valence(n_valence,3*numa),&
         hesse_diagonal(n_valence),&
         hesse_cartesian(3*numa,3*numa),&
         hesse_nodummies(3*n_atoms,3*n_atoms),&
         STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         (" valence_setup: allocation (2) failed")
    hesse_diagonal = zero
    hesse_cartesian = zero

    call set_hesse(v,n_valence,hesse_diagonal)

    write(OPT_STDOUT,*)" valence_setup : these are the diagonals of the hesse matrix"
    write(OPT_STDOUT,*)hesse_diagonal
    write(OPT_STDOUT,*)" "
    call generate_bmat(v,n_valence,bmat_valence)

    ! diagonal in valence -> square in cartesians:
    call transform_hesse(hesse_diagonal,hesse_cartesian,bmat_valence)
    
    deallocate(hesse_diagonal,bmat_valence,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" valence_setup: deallocation (1) failed")

    ! remove rows and columns in hess_cartesian corresponding to dummy atoms
    DCALL octave("with dummies",hesse_cartesian)
    call rem_dummies(hesse_cartesian,hesse_nodummies)
    DCALL octave("without dummies",hesse_nodummies)

    ! write hesse_nodummies to hesse_cartesian.dat:
    call store_hesse(hesse_nodummies)

    write(OPT_STDOUT,*)" -----  Hesse Matrix in cartesian coordinates -------"
    call print_matrix(hesse_nodummies,3*n_atoms,3*n_atoms,n_col)
    write(OPT_STDOUT,*)" ----------------------------------------------------"
    write(OPT_STDOUT,*)" "

    deallocate(hesse_cartesian,hesse_nodummies,v,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" valence_setup : deallocation (2) failed")
  end subroutine valence_setup

  ! ********************************************************************

  subroutine set_hesse(coor,n_coor,diag)
    ! Purpose: set the diagonal values of the hesse matrix in
    !          valence coordinates.
    !---------------------------------------------------------
    type(int_coor),intent(in)          :: coor(:)
    integer(kind=i4_kind),intent(in)   :: n_coor
    real(kind=r8_kind)                 :: diag(:)
    !-------- DECLARATION OF LOCAL VARIABLES -----------------
    integer(kind=i4_kind) :: p1,p2,b1,b2,row1,row2,&
                             index1,index2,index3,i,ap
    real(kind=r8_kind)    :: a,b,r,r_cov,cov1,cov2
    real(kind=r8_kind),parameter   :: dummy_force=0.0_r8_kind
    !-------- executable code --------------------------------
    
    if (ubound(diag,1)/=n_coor) call error_handler&
         (" set_hesse: dimensions wrong")
    ! this can either be 'n_valence' if the forcefield hessian
    ! is set up for Z-matrix coords (plain or as the underlying set
    ! for delocalized coords)
    ! or n_primitive (which would be equal to n_valence) if
    ! the forcefiled hessian is set up for automatic coordinates
    ! (=delocalized coords with valence coords as the underlying set)
    do i=1,n_coor
       if (coor(i)%typ == b_length) then
          p1 = coor(i)%length%partner1
          p2 = coor(i)%length%partner2
          if (atom(p1)%charge == dummy_charge .or. &
               atom(p2)%charge == dummy_charge ) then
             diag(i) = dummy_force
             cycle
          endif
          index1 = nint(atom(p1)%charge)
          index2 = nint(atom(p2)%charge)
          if(.not.estimate_hessian_m) then
             row1=get_row(index1)
             row2=get_row(index2)
             call force_parameter(row1,row2,a,b,b_length)
             diag(i) = a/(coor(i)%value-b)**3
          else ! model force-field Hessian
             call set_K_parameter(k_par=diag(i), i=p1, j=p2, type=b_length)
          end if

       elseif(coor(i)%typ == b_angle) then

          p1 = coor(i)%angle%partner1
          p2 = coor(i)%angle%partner2
          ap = coor(i)%angle%apex
          if (atom(p1)%charge == dummy_charge.or.&
               atom(p2)%charge == dummy_charge .or. &
               atom(ap)%charge == dummy_charge ) then
             diag(i) = dummy_force
             cycle
          endif
          index1 = nint(atom(p1)%charge)
          index2 = nint(atom(p2)%charge)
          index3 = nint(atom(ap)%charge)
          if(.not.estimate_hessian_m) then         
             row1=get_row(index1)
             row2=get_row(index2)
             call force_parameter(row1,row2,a,b,b_angle)
             diag(i) = a
          else ! model force-field Hessian
             call set_K_parameter(k_par=diag(i), i=p1, j=ap, k=p2, type=b_angle)
          end if
          
       elseif(coor(i)%typ == d_angle) then

          b1 = coor(i)%dihedral%base1
          b2 = coor(i)%dihedral%base2
          p1 = coor(i)%dihedral%partner1
          p2 = coor(i)%dihedral%partner2
          if (atom(p1)%charge == dummy_charge .or.&
               atom(p2)%charge == dummy_charge .or. &
               atom(b1)%charge == dummy_charge .or. &
               atom(b2)%charge == dummy_charge ) then
             diag(i) = dummy_force
             cycle
          endif
          index1 = nint(atom(b1)%charge)
          index2 = nint(atom(b2)%charge)
          if(.not.estimate_hessian_m) then 
             cov1 = covrad(index1)
             cov2 = covrad(index2)
             row1 = get_row(index1)
             row2 = get_row(index2)
             r_cov = cov1+cov2
             r = abs_value(atom(b1)%x - atom(b2)%x)
             call force_parameter(row1,row2,a,b,d_angle)
             if (r<=r_cov+a/b) then
                diag(i) = a - b*(r-r_cov)
             else
                diag(i) = 0.001_r8_kind
             endif
          else ! There is no torsion parameters
             diag(i) = dummy_force
          end if
       else
          write(*,*)" set_hesse: this bond type is not supported yet"
       endif
    enddo

  end subroutine set_hesse

   !*************************************************************

  subroutine transform_hesse(hdiag,hcart,bmat)
    ! Purpose: transforms the (diagonal) hesse matrix in valence
    !          coordinates to cartesian coordinates using the
    !          B-Matrix 'bmat_valence':
    !          hesse_cartesian = bmat_valence*hesse_diag*bmat_valenceT
    !          where bmat_valenceT is the transpose of 'bmat_valence'.
    !------ Declaration of local variables ---------------------
    implicit none
    real(kind=r8_kind), intent(in)  :: hdiag(:)   ! (n_valence)
    real(kind=r8_kind), intent(out) :: hcart(:,:) ! (3*(n_atoms+n_dummy))^2 ?
    real(kind=r8_kind), intent(in)  :: bmat(:,:)  ! (n_valence,3*(n_atoms+n_dummy))
    ! *** end of interface ***

    real(kind=r8_kind)     :: help(n_valence,3*(n_atoms+n_dummy))
    integer(kind=i4_kind)  :: i
    !------ executable code ------------------------------------

ASSERT(all(shape(help)==shape(bmat)))
ASSERT(size(hcart,1)==size(hcart,2))
ASSERT(size(hcart,1)==size(bmat,2))

    do i=1,n_valence
       help(i,:) = hdiag(i)*bmat(i,:)
    enddo
    hcart = matmul(transpose(bmat),help)
  end subroutine transform_hesse

  subroutine store_hesse(hesse_cartesian)
    use filename_module, only: inpfile
    use iounitadmin_module
    ! Purpose: store the (cartesian) Hesse matrix in a file
    !          called '$DATADIR/hesse_cartesian.dat'. From this
    !          file the hesse matrix can be read in by other
    !          routines.
    !------ Modules used ---------------------------------------
    implicit none
    real(r8_kind), intent(in) :: hesse_cartesian(:,:)
    ! *** end of interface ***

    !----- Declaration of local variables ----------------------
    integer(kind=i4_kind)    :: io_hesse_cart
    !----- executable code -------------------------------------

    DPRINT 'store_hesse: writing cartesian hessian to hesse_cartesian.dat'
    io_hesse_cart=openget_iounit(status='unknown',form='formatted',file=&
         trim(inpfile('hesse_cartesian.dat')))
    write(io_hesse_cart,*) hesse_cartesian
    call returnclose_iounit(io_hesse_cart)
  end subroutine store_hesse
  
  subroutine force_parameter(row1,row2,a,b,typ)
    !----------------------------------------------------------------
    ! Purpose : returns the appropriate force parameters for
    !           calculation the force constants (diagonals of
    !           the hesse matrix in valence coordinates).
    !
    ! Subroutine called by:
    !----------------------------------------------------------------
    !  Author: FN
    !  Date  : -- 
    !  References:  H.B.Schlegel, Theor.Chim.Acta 66,333-340 (1984)
    !
    !----------------------------------------------------------------
    ! Modifications
    !----------------------------------------------------------------
    !
    ! Modification 
    ! Author: VVP 
    ! Date:   06/05
    ! References: Wittbrodt, J.M., Schlegel, H.B. Estimating stretching force constants
    ! for geometry optimization // Theochem. - 1998. - V. 398-399. - P. 55-61
    !----------------------------------------------------------------
    !
    !------- Modules used -----------------------------------------
    use vff_hessian
    !------- Declaration of formal parameters ---------------------
    integer(kind=i4_kind),intent(in)  :: row1,row2
    integer(kind=i4_kind),intent(in)  :: typ
    real(kind=r8_kind),intent(out)    :: a,b
    !------- Declaration of local variables -----------------------
    integer(kind=i4_kind)             :: r1,r2
    !------- executable code --------------------------------------
     select case(typ)
      case(b_length)
       r1=row1
       r2=row2
       if(row1>6) r1=6
       if(row2>6) r2=6
       a=1.734_r8_kind
       b=stretch_par(r1,r2)
      case(b_angle)
       r1=row1
       r2=row2
       if(row1>2) r1=2
       if(row2>2) r2=2
       b=99.0_r8_kind
       a=angle_par(r1,r2)
       case(d_angle)
       a=torsion_a
       b=torsion_b
      case default
       write(OPT_STDOUT,*)" force_parameter: this bond type is not supported yet"
     end select      
  end subroutine force_parameter
  
  subroutine set_K_parameter(k_par,i,j,k,type)
    !
    ! Purpose : returns the force parameter in accordance with model
    !           stretching-bending potential
    ! Reference.: M.V.Fernandez_Serra, E.Artacho, J.N.Soler
    !             Phys.Rev.B.67,100101(R) (2003)
    !
    ! Subroutine called by: set_hesse
    !
    !------- Modules used -----------------------------------------
    !------- Declaration of formal parameters ---------------------
    real(kind=r8_kind),              intent(inout) :: k_par
    integer(kind=i4_kind), optional, intent(in)    :: i,j,k,type
    !------- Declaration of local variables -----------------------
    real(kind=r8_kind), parameter :: A = 308.722_r8_kind, &
                                     B = 0.1_r8_kind
    real(kind=r8_kind)            :: r_ij,r_jk, R_i,R_j,R_k
    !------- executable code --------------------------------------

    select case(type)
       case(b_length)
          write(*,*) "K_parameter : i,j    ",i,j
          write(*,*) "K_parameter : cov    ",covrad(nint(atom(i)%charge)),covrad(nint(atom(j)%charge))
          write(*,*) "K_parameter : atom i ", atom(i)%x
          write(*,*) "K_parameter : atom j ", atom(j)%x
          r_ij  = sqrt( dot_product( atom(i)%x-atom(j)%x, atom(i)%x-atom(j)%x) ) 
          write(*,*) "K_parameter : r_ij ", r_ij
          R_i   = covrad(nint(atom(i)%charge))
          R_j   = covrad(nint(atom(j)%charge))
          K_par = A*((R_i+R_j)/r_ij)**8 / 10.0_r8_kind
          write(*,*) "K_parameter : K    ", K_par
       case(b_angle)
          r_ij  = sqrt( dot_product(atom(i)%x-atom(j)%x,atom(i)%x-atom(j)%x) )
          r_jk  = sqrt( dot_product(atom(j)%x-atom(k)%x,atom(j)%x-atom(k)%x) )
          write(*,*) "K_parameter : i,j    ",i,j,k
          write(*,*) "K_parameter : cov    ",covrad(i),covrad(j),covrad(k)
          write(*,*) "K_parameter : atom i ", atom(i)%x
          write(*,*) "K_parameter : atom j ", atom(j)%x
          write(*,*) "K_parameter : atom k ", atom(k)%x
          R_i   = covrad(nint(atom(i)%charge))
          R_j   = covrad(nint(atom(j)%charge)) 
          R_k   = covrad(nint(atom(k)%charge));
          K_par = B*A * ((R_i+R_j)/r_ij)**4 * ((R_i+R_j)/r_ij)**4 * r_ij*r_jk  / 10.0_r8_kind
          write(*,*) "K_parameter : K    ", K_par
    end select

  end subroutine set_K_parameter

  subroutine rem_dummies(hdumm,hcart)
     ! Purpose: Removes the rows and columns of hesse_cartesian
     !          corresponding to dummy atoms
     !
     ! subroutine called by: valence_setup
     ! 
     !------- Modules used -----------------------------------------
     use opt_data_module, only: n_atoms, n_dummy
     !------- Declaration of formal parameters ---------------------
     implicit none
     real(r8_kind), intent(in)  :: hdumm(3*(n_atoms+n_dummy),3*(n_atoms+n_dummy))
     real(r8_kind), intent(out) :: hcart(3*n_atoms,3*n_atoms)
     !------- Declaration of local variables -----------------------
     integer(kind=i4_kind):: i, j, k, l
     !------- executable code --------------------------------------
     i = 3*n_atoms
     j = 3*(n_atoms+n_dummy)
     ASSERT(size(hcart,1)==i)
     ASSERT(size(hcart,2)==i)
     ASSERT(size(hdumm,1)==j)
     ASSERT(size(hdumm,2)==j)

     ! remove zero entries at the positions which belong to dummy atoms

     k=0
     do i=0,n_atoms+n_dummy-1
        if(atom(i+1)%dummy) CYCLE
        l=0
        do j=0,n_atoms+n_dummy-1
           if(atom(j+1)%dummy) CYCLE
           hcart(k*3+1:k*3+3,l*3+1:l*3+3) = &
           hdumm(i*3+1:i*3+3,j*3+1:j*3+3)
           l=l+1
        end do
        k=k+1
     end do

  end subroutine rem_dummies

!--------------- End of module ----------------------------------
end module valence_coord_module
