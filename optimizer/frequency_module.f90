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
module frequency_module
  !-------------------------------------------------------------------
  !  Purpose: module for calculation of the vibrational frequencies
  !
  !
  !  Module called by: main_opt
  !
  !  References: optimizer documentation
  !
  ! Author: M. Staufer
  ! Date:   3/98
  !
  !-------------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use iounitadmin_module
  use opt_data_module, only: OPT_STDOUT
  !
  ! I am getting crazy resolving dependencies !!!
  ! Dont _USE_ modules in that way, I beg you !!!
  ! AM
  !
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private

  !------------ Declaration of constants and variables ---------------
  real(kind=r8_kind), parameter :: freq_coeff=5143.05_r8_kind

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: frequency_main
  public :: freq ! (k,m,w,x), m -- diagonal
  public :: freq_print !(iou,m,w,x), w--frequency, x--mode
  public :: dipole_der !(n_mode,dipder)
  !------------ Subroutines ------------------------------------------
contains
  !*************************************************************
  subroutine frequency_main(hesse)
    ! Purpose: main routine for calculation of frequencies
    ! Subroutine called by: main_opt
    use constants , only             : ZERO                  &
                                     , HUNDRED               &
                                     , u_Atom_Mass           &    !   atomic mass constant
                                     , a_Bohr                &    !   bohr-radius -> meter
                                     , E_h2Jmol              &    ! conversion Eh -> j/mol
                                     , c_speedoflight             !  vacuum speed of light
    use opt_data_module, only        : n_internal            &
                                     , n_primitive           &
                                     , n_tot_atoms           &
                                     , atom                  &
                                     , s_prim                &
                                     , calculate_intensities &
                                     , gx_test               &
                                     , io_flepo              &
                                     , step_size             &
                                     , temperature           &
                                     , Delta_temperature     &
                                     , N_temperature_steps   &
                                     , pressure       &
                                     , symmetry_index      ! number of identical positions
                                                                ! accessible via rotations
    use coordinates_module, only:&
         & invert_bmat,bmat,bmat_inv,expand_mat,reduc_mat,&
         & alloc_bmat,free_bmat,generate_bmat
    use math_module, only: small,print_matrix
    use matrix_eigenval, only        : eigs         ! Eigensolver (for moments of inertia)
    use geo_operations_module , only : minerttensor &   ! Calculation of tensor of moments
                                     , centerofmass        ! of inertia and center of mass
    use gradient_module , only       : energy                    ! total electronic energy
    use thermodyn_prop_module, only  : thermodynamic_properties
    implicit none
    !----------------------------------------------------
    ! --- Declaration of formal parameters --------------
    real(r8_kind), intent(inout) :: hesse(:,:) ! Hessian matrix, INTENT(INOUT) because it is modified in-place
                                               ! when transformed into mass-weighted coordinates
    ! *** end of interface ***

    ! --- Declaration of local variables ----------------
    integer(kind=i4_kind) :: alloc_stat, i, j, i_atom, k, ierr, work_dim, matz
    real(kind=r8_kind)    :: mineig
    real(r8_kind)         :: velocity(3) ! 3-vector
    real(kind=r8_kind), allocatable :: mass_matrix(:,:),help_mat(:,:),&
         mass_mat_red(:,:),help_mat2(:,:),help_mat3(:,:),mass_matrix_heavy(:,:),&
         eig_real(:), eig_vec(:,:), mmhalf(:,:), fv1(:), fv2(:), kin_energy(:), &
         eig_xyz(:,:),intensities(:), mode_pop(:), eig_vec_int(:,:)
    real(kind=r8_kind), allocatable :: freq_prep_vec(:)

    !*** For thermodynamic module ********************************************************
    real(kind=r8_kind), dimension(3)        :: masscenter = ZERO &        ! center of mass
                                             , Mom_Inert = ZERO    ! mom of inertia vector

    real(kind=r8_kind), dimension(3,3)      :: Rot_Tensor = ZERO & ! mom of inertia tensor
                                             , EV_Rot_Tensor = ZERO       !rotational axis
                                                                              ! (not used)

    ! --- Executable code -------------------------------
    ! first we build the inverse of the B-Matrix
    call alloc_bmat(n_primitive)
    call generate_bmat(s_prim,n_primitive,bmat)
    call invert_bmat(bmat,bmat_inv,spectroscopic=.true.)
    ! call free_bmat()
    ! then we build the mass matrix
    ! now allocate
    work_dim=10*n_primitive
    allocate(mass_matrix(n_primitive,n_primitive), &
         mass_matrix_heavy(n_primitive,n_primitive), &
         mass_mat_red(n_internal,n_internal), &
         help_mat(n_internal,n_internal), &
         help_mat2(n_internal,n_primitive), &
         help_mat3(n_primitive,n_internal), &
         mmhalf(n_internal,n_internal), &
         eig_real(n_internal), fv1(n_internal), fv2(n_internal), &
         eig_vec(n_internal,n_internal), &
         eig_vec_int(n_internal,n_internal), &
         kin_energy(n_tot_atoms), &
         eig_xyz(3*n_tot_atoms,n_internal), &
         intensities(n_internal), &
         mode_pop(n_internal), &
         stat=alloc_stat)
    if(alloc_stat/=0) &
         stop 'frequency_main: allocation of mass_matrix_failed'
    mass_matrix=0.0_r8_kind
    mass_matrix_heavy=0.0_r8_kind

    print*, 'frequency_main: atomic masses'
    do i_atom=1,n_tot_atoms
      print* , atom(i_atom)%mass
    enddo

    ! now build mass_matrix
    do i=1,n_primitive
       do j=1,n_primitive
          do i_atom=1,n_tot_atoms
             if(.not.atom(i_atom)%heavy) then
                do k=1,3
                   mass_matrix(i,j)=mass_matrix(i,j)+atom(i_atom)%mass*bmat_inv((i_atom-1)*3+k,j)*&
                        bmat_inv((i_atom-1)*3+k,i)
                end do
             else
                do k=1,3
                   mass_matrix_heavy(i,j)=mass_matrix_heavy(i,j)+bmat_inv((i_atom-1)*3+k,j)*&
                        bmat_inv((i_atom-1)*3+k,i)
                end do
             end if
          end do
       end do
    end do

    DPRINT  'frequency_main: now check if infinite heavy atoms are moved'

    help_mat2= matmul(reduc_mat,mass_matrix_heavy)
    mass_mat_red= matmul(help_mat2,transpose(reduc_mat))
    do i=1,n_internal
       do j=1,n_internal
          if( abs(mass_mat_red(i,j))>small) then
             write(OPT_STDOUT,*) 'Error: you tried to move atoms with infinite masses'
             write(OPT_STDOUT,*) 'm_mat(',i,',',j,')=',mass_mat_red(i,j)
             call error_handler('frequency_main: see above')
          end if
       end do
    end do
    ! we have to reduce mass_matrix to the free and symmetry non equivalent coordinates
    help_mat2= matmul(reduc_mat,mass_matrix)
    mass_mat_red= matmul(help_mat2,transpose(reduc_mat))

    write(OPT_STDOUT,*) '   ---------- mass matrix ----------   '
    call print_matrix(mass_mat_red,n_internal,n_internal,min(n_internal,10))

    ! now we build m-1/2
    ! to do this we diagonalize mass_mat_red
    mmhalf=mass_mat_red
    matz=1 ! eigenvalues and eigenvetors are required
    call  rs(n_internal,n_internal,mmhalf,eig_real,matz,eig_vec,fv1,fv2,ierr)
    DPRINT 'mass_mat_red diagonalized'

    if(ierr/=0)  stop 'frequency_main: diagonalisation failed'
    mineig=minval(eig_real)

    if(mineig<small) then
       ! we have a internal mode, which is not connected to
       ! a movement of masses
       write(OPT_STDOUT,*) "Error: The mass matrix has zero eigenvalues. Please check if&
            & your definition of internal variables contains a movement of dummy &
            & atoms"
       help_mat3=matmul(expand_mat,eig_vec)
       eig_xyz=matmul(bmat_inv,help_mat3)

       do i=1,n_internal
          if(eig_real(i)<small) then
             write(OPT_STDOUT,*) 'Eigenvector:'
             write(OPT_STDOUT,'(10F8.4)') eig_vec(:,i)
             write(OPT_STDOUT,*) 'has eigenvalue:',eig_real(i)
             ! to make debug easier we also give the eigenvector in cartesian coordinates
             write(OPT_STDOUT,*) 'Eigenvector in cartesian coordinate:'
             write(OPT_STDOUT,'(10F8.4)') eig_xyz(:,i)
          end if
       end do
       call error_handler('frequency_main: mass_matrix has zero or negative eigenvectors')
    end if
    If(gx_test) then
       write(OPT_STDOUT,*) 'There are no redundant coordinates in your gxfile'
       goto 999 ! leave this subroutine, but do deallocation before
    end if
    DPRINT 'frequency_main: mass_matrix controlled'

    eig_real=1.0_r8_kind/sqrt(eig_real)
    mmhalf=eig_vec
    do i=1,n_internal
       mmhalf(:,i)=mmhalf(:,i)*eig_real(i)
    end do
    mmhalf=matmul(mmhalf,transpose(eig_vec))

    DPRINT 'frequency_main: now transform hessian into mass-weighted coordinates, shape(hesse)=',shape(hesse)

    ! FIXME: dont modify the input argument in-place, use a temp copy:
    hesse=matmul(hesse,mmhalf)
    hesse=matmul(mmhalf,hesse)

    DPRINT 'frequency_main: now we do the final diagonalization'
    call  rs(n_internal,n_internal,hesse,eig_real,matz,eig_vec,fv1,fv2,ierr)
    DPRINT 'frequency_main: ierr=',ierr

    ! now we determine the eigenvectors in cartesian coordinates
    eig_vec_int=matmul(mmhalf,eig_vec)   ! eigenvector in internal coordinates
    help_mat3=matmul(expand_mat,eig_vec_int)
    eig_xyz=matmul(bmat_inv,help_mat3) ! eigenvector in cartesian coordinates.
    ! with the above quantity we can calculate the kinetic energy for every mode
    ! and every atom
    intensities=0.0_r8_kind


    if(calculate_intensities) &
         call intensity_calc(n_internal,eig_vec_int,step_size,intensities)

    !
    ! now we can allready print the frequencies
    !
    do i=n_internal,1,-1 ! output high-freq modes first
       ! mode index counted from high-freq end, for output only:
       k = n_internal - i + 1

       do i_atom=1,n_tot_atoms
          if(.not.(atom(i_atom)%dummy.or.(atom(i_atom)%heavy))) then

             ! "velocity" of an atom when mode "i" is active extracted from mode vector:
             velocity(:) = eig_xyz((i_atom-1)*3+1:(i_atom-1)*3+3,i)

             ! (twice the) "kinetic energy" of an atom when mode "i" is active:
             kin_energy(i_atom) = atom(i_atom)%mass * sum(velocity**2)

             ! eigenvectors are normalized so that
             !
             !   v' M v = 1
             !
             ! therefore the "kinetic energy" is just a contribuiton to
             ! this kind of norm of eigenvector that sum to one over all atoms.
          else
             kin_energy(i_atom)=0.0_r8_kind
          end if
       end do
       if(calculate_intensities) then
          if( k == 1 )then ! print table header:
             write(io_flepo,*)
             write(io_flepo,*) "Vibrational frequencies, intensities and partition of kinetic energy [in %] over atoms:"
             write(io_flepo,*)
             write(io_flepo,1010) 'I','FREQUENCY','INTENSITY',(j,j=1,size(kin_energy))
          endif

1010      format(A4,A10 ,4X,A9  ,6X,(10I7  ),/:(33X,10I7))
1020      format(I4,F9.3,A5,F9.3,A6,(10F7.1),/:(33X,10F7.1))

          if(eig_real(i)>0) then
             write(io_flepo,1020) k, sqrt(eig_real(i)) * freq_coeff, 'cm-1', intensities(i), 'km/mol', &
                  kin_energy * 100
          else
             write(io_flepo,1020) k, sqrt(-eig_real(i)) * freq_coeff,'icm-1', intensities(i), 'km/mol', &
                  kin_energy * 100
          endif
       else
          if( k == 1 )then ! print table header:
             write(io_flepo,*)
             write(io_flepo,*) "Vibrational frequencies and partition of kinetic energy [in %] over atoms:"
             write(io_flepo,*)
             write(io_flepo,1030) 'I','FREQUENCY',(j,j=1,size(kin_energy))
          endif

1030      format(A4,A10 ,4X,(10I7  ),/:(18X,10I7))
1040      format(I4,F9.3,A5,(10F7.1),/:(18X,10F7.1))

          if(eig_real(i)>0) then
             write(io_flepo,1040) k, sqrt(eig_real(i)) * freq_coeff, 'cm-1', kin_energy * 100
          else
             write(io_flepo,1040) k, sqrt(-eig_real(i)) * freq_coeff, 'icm-1', kin_energy * 100
          endif
       end if
    end do

    !
    ! Print detailed "population" analysis of mode vectors in internal coordinates
    !
    write(io_flepo,*)
    write(io_flepo,*) 'Vibrational frequencies and contributions of internal modes [in %]:'
    write(io_flepo,*)
    write(io_flepo,1050) 'I','FREQUENCY',(j,j=1,size(mode_pop))

1050 format(A4,A10 ,4X,(10I7  ),/:(18X,10I7  ))
1060 format(I4,F9.3,A5,(10F7.1),/:(18X,10F7.1))

    do i=n_internal,1,-1
       ! mode index counted from high-freq end, for output only:
       k = n_internal - i + 1

       ! FIXME: quadratic in "eig_vec_int", why not just printing mode vectors?
       mode_pop = matmul(mass_mat_red, eig_vec_int(:,i))
       do j=1,n_internal
          mode_pop(j) = mode_pop(j) * eig_vec_int(j,i)
       end do

       if(eig_real(i)>0) then
          write(io_flepo,1060) k, sqrt(eig_real(i)) * freq_coeff, 'cm-1', mode_pop * 100
       else
          write(io_flepo,1060) k, sqrt(-eig_real(i)) * freq_coeff, 'icm-1', mode_pop * 100
       endif
    end do
    call movie_maker((eig_xyz), eig_real)

    !*** Start of the thermodynamic properties section ***********************************
    allocate(freq_prep_vec(size(eig_real)))

    !*** Write processed eigenvalues into frequency preparation vector "freq_prep_vec" ***
    WHERE(eig_real > 0)
      freq_prep_vec = sqrt(eig_real)
    ELSEWHERE
      freq_prep_vec = -sqrt(-eig_real)
    END WHERE

    !*** Conversion to SI-units [Hz] *****************************************************
    freq_prep_vec = freq_prep_vec*freq_coeff*HUNDRED*c_speedoflight

    !*** Calculate center of mass ********************************************************
    masscenter = centerofmass(atom(:)%x(1) &
                             ,atom(:)%x(2) &
                             ,atom(:)%x(3) &
                             ,atom(:)%mass)

    !*** Calculate tensor of moments of inertia ******************************************
    Rot_Tensor =  minerttensor(atom(:)%x(1) &
                              ,atom(:)%x(2) &
                              ,atom(:)%x(3) &
                              ,atom(:)%mass &
                              ,masscenter)

    !*** Calculate moments of inertia ****************************************************
    call eigs(Rot_Tensor,Mom_Inert,EV_Rot_Tensor)

    !*** Conversion to SI-units [kg*m^2] *************************************************
    Mom_Inert = Mom_Inert * u_Atom_Mass * a_Bohr**2

    !FIXME: include g

    !*** calculate and print thermodynamic corrections ***********************************
    !*** adapted for temperature intervals ***********************************************
    if((Delta_temperature==0.0_r8_kind .and. N_temperature_steps/=1).or. &
       (Delta_temperature/=0.0_r8_kind .and. N_temperature_steps==1).or. &
       (Delta_temperature < 0.0_r8_kind .or. N_temperature_steps < 1 )) then
      write(      06,*) "ERROR: no proper addressation of thermodynamic properties module in temperature-interval modus"
      write(      06,*)
      write(io_flepo,*) "ERROR: no proper addressation of thermodynamic properties module in temperature-interval modus"
      write(io_flepo,*)
    else
      call thermodynamic_properties(      06      ,     freq_prep_vec     ,     Mom_inert &
                                   ,sum(atom(:)%mass)*u_Atom_Mass    ,    energy*E_h2Jmol &
                                   ,symmetry_index,temperature,pressure,Delta_temperature &
                                   ,N_temperature_steps)
      call thermodynamic_properties(io_flepo      ,     freq_prep_vec     ,     Mom_inert &
                                   ,sum(atom(:)%mass)*u_Atom_Mass    ,    energy*E_h2Jmol &
                                   ,symmetry_index,temperature,pressure,Delta_temperature &
                                   ,N_temperature_steps)
    endif
    deallocate(freq_prep_vec)
    !*** End of thermodynamic properties part ********************************************

999 deallocate(mass_matrix,eig_real,help_mat,eig_vec,eig_vec_int,mmhalf,&
         help_mat2,help_mat3,mass_mat_red,mass_matrix_heavy,fv1,fv2,kin_energy,&
         eig_xyz, intensities, mode_pop, stat=alloc_stat)
    if(alloc_stat/=0) &
         stop 'frequency_main: deallocation_failed'
    call free_bmat()
  end subroutine frequency_main
  !*************************************************************

  subroutine freq(k,m,w,x,dipder,intensities)
    ! THIS IS ONLY A COPY FROM gradient_data_module
    ! IF EDIT DO IT AT TWO PLACES!
    ! YET BETTER -- FIND A SUTABLE PLACE FOR IT!
    !
    ! Solves the eigenvalue equation:
    !
    !     ( K - w2 * M ) X = 0
    !
    ! with force matrix K and mass matrix M
    !
    ! Outputs result on "iou"
    !
    use matrix_eigenval, only: eigs
    use opt_data_module, only: n_atoms
    implicit none
    real(r8_kind)   , intent(in)  :: m(:)   ! (3*NA)
    real(r8_kind)   , intent(in)  :: k(:,:) ! (3*NA,3*NA)
    real(r8_kind)   , intent(out) :: w(:)
    real(r8_kind)   , intent(out) :: x(:,:)
    real(r8_kind)   , intent(in)  :: dipder(:,:) !(3, 3*NA)
    real(r8_kind)   , intent(out) :: intensities(:) !(3*NA)
    optional                      :: intensities,dipder
    ! *** end of interface ***

    integer(i4_kind) :: i, i_mode
    real(r8_kind)    :: h(size(m),size(m))
    real(r8_kind)    :: dipderq(3,size(m)) !(3, 3*NA)

    ! ``orthonormalize'' generalized eigenvalue problem:
    do i=1,size(m)
         h(:,i) = k(:,i) / sqrt(m(i))
    enddo
    do i=1,size(m)
         h(i,:) = h(i,:) / sqrt(m(i))
    enddo

    ! solve eigenvalue problem:
    call eigs(h,w,x)

    ! square roots of w^2, negative==imaginary:
    where( w > 0.0_r8_kind )
       w = sqrt(w)
    elsewhere
       w = - sqrt(-w)
    endwhere

    ! back-transform eigenvectors:
    do i=1,size(m)
      x(i,:) = x(i,:) / sqrt(m(i))
    enddo

    ! now calculate intensities
    if( present(intensities)) then
      do i=1,3
         dipderq(i,:)=matmul(dipder(i,:),x)
      end do

      do i_mode=1,3*n_atoms
         intensities(i_mode)=sum(dipderq(:,i_mode)**2)*1000
      end do
    endif

  end subroutine freq
  !****************************************************

  subroutine freq_print(iou,m,w,x,intensities)
    use atom_data_module, only: nuc_mass!(1&6)

    implicit none
    integer(i4_kind), intent(in) :: iou
    real(r8_kind)   , intent(in) :: m(:)
    real(r8_kind)   , intent(in) :: w(:)
    real(r8_kind)   , intent(in) :: x(:,:)
    real(r8_kind)   , intent(in) :: intensities(:)
    optional                     :: intensities
    ! *** end of interface ***

    integer(i4_kind) :: i
    real(r8_kind), parameter   :: mp = 1836.15267261_r8_kind ! Mp/Me
    real(r8_kind)              :: m1

    real(r8_kind), parameter   :: au2ev = 27.211658_r8_kind ! eV/au
    real(r8_kind), parameter   :: ev2cm = 8065.5410_r8_kind ! cm-1/eV
    real(r8_kind), parameter   :: au2cm = au2ev * ev2cm     ! cm-1/au
    real(r8_kind), parameter   :: au2D  = 2.541584_r8_kind ! debye/au
    real(r8_kind), parameter   :: au2A  = 0.529177_r8_kind ! angstram/aueV
    real(r8_kind), parameter   :: D_A   = 23.071130228_r8_kind ! (au2D/au2A)^2
    ! for conversion between different unit of Intensity, please follow JCP 84, 2262 (1986)

    ! atomic mass unit in (electronic) a.u.:
    m1 = (nuc_mass(6)/12) * (mp/nuc_mass(1))

    write(iou,'("SD:")')
    write(iou,'("SD: MASS MATRIX:")')
    write(iou,'("SD:",A4,A20,A20)') "I","MASS(I) [C12]","MASS(I) [Me]"
    do i=1,size(m)/3
    write(iou,'("SD:",I4,2F20.12)') i, m(3*i-2)/m1, m(3*i-2)
    enddo

    if( present(intensities)) then
      write(iou,'("SD:")')
      write(iou,'("SD: FREQUENCIES AND INTENSITIES (NUMERICAL CARTESIAN HESSIAN):")')
      write(iou,'("SD:",A4,5A20)') "I","OMEGA(I) [au]","OMEGA(I) [cm-1]", &
                      "[(D/A)^2amu^-1]", "[Km/mol]", "[cm^-2atm^-1]"
      do i=1,size(w)
      write(iou,'("SD:",I4,5F20.12)') i, w(i), w(i)*au2cm, &
                      intensities(i)*D_A, intensities(i)*974.8706077, intensities(i)*3960.159
      enddo
      ! conversion factor from (D/A)^2amu^-1 to Km/mol, 1(D/A)^2amu^-1 = 42.255 Km/mol
      ! conversion factor from (D/A)^2amu^-1 to cm^-2atm^-1, 1(D/A)^2amu^-1 = 171.65 cm^-2atm^-1
    else
      write(iou,'("SD:")')
      write(iou,'("SD: FREQUENCIES (NUMERICAL CARTESIAN HESSIAN):")')
      write(iou,'("SD:",A4,2A20)') "I","OMEGA(I) [au]","OMEGA(I) [cm-1]"
      do i=1,size(w)
      write(iou,'("SD:",I4,5F20.12)') i, w(i), w(i)*au2cm
      enddo
    endif

    write(iou,'("SD:")')
    write(iou,'("SD: MODE VECTORS (NUMERICAL CARTESIAN HESSIAN):")')
    write(iou,'("SD:",A4,A20,A20)') "I","OMEGA(I) [cm-1]","X(:,I) [arb.u.]"
    do i= 1, size(w)
      write(iou,'("SD:",I4,F20.12,3F20.12,/:("SD:",24X,3F20.12))') &
           i, w(i)*au2cm, x(:,i)
    enddo
  end subroutine freq_print

  !*************************************************************
  subroutine movie_maker(eig_xyz, eig_val)
    use opt_data_module, only:&
         & n_internal,n_atoms,n_dummy,n_tot_atoms,&
         & atom
    use math_module, only: pi
    use constants, only: angstrom
    use atom_data_module, only: symbol
    implicit none
    ! Purpose: make a movie of the normal modes for the xmol programm
    ! Subroutine called by: frequency_main
    !----------------------------------------------------
    ! --- Declaration of formal parameters --------------
    real(kind=r8_kind) :: eig_xyz(:,:) ! eigenvectors of vibrational modes in
    real(kind=r8_kind), intent(in) :: eig_val(:)
    ! cartesian coordinates
    ! --- Declaration of local variables ----------------
    real(kind=r8_kind) :: factor
    real(kind=r8_kind) :: freq
    integer(kind=i4_kind) :: i_mode, i_atom, i,io_xmol=11
    character(len=7) :: extension
    character(len=3)      :: sym
    character(len=6)      :: freq_unit
    ! --- Executable code -------------------------------
    do i_mode=1,n_internal
       ! first we open the file
       write(extension,'(i4)') i_mode
       extension=adjustl(extension)
       io_xmol=openget_iounit(status='unknown',form='formatted',&
            file="xmol_mode."//trim(extension))

       ! Print also the frequency in the xmol files
       if(eig_val(i_mode)>0) then
          freq =  sqrt(eig_val(i_mode)) * 5143.05_r8_kind
          freq_unit = 'cm-1'
       else
          freq =  sqrt(-eig_val(i_mode)) * 5143.05_r8_kind
          freq_unit = 'i cm-1'
       endif

       do i=0,20
          factor=sin(pi/10.0_r8_kind*real(i,r8_kind))
          write(io_xmol,*)n_atoms+n_dummy
          write(io_xmol,*)"interesting molecule, Freq : ", freq, freq_unit
          do i_atom=1,n_tot_atoms
             sym = symbol(nint(atom(i_atom)%charge))
             write(io_xmol, 1000) sym, &
                  (atom(i_atom)%x(1) + factor * eig_xyz((i_atom - 1) * 3 + 1, i_mode)) / angstrom, &
                  (atom(i_atom)%x(2) + factor * eig_xyz((i_atom - 1) * 3 + 2, i_mode)) / angstrom, &
                  (atom(i_atom)%x(3) + factor * eig_xyz((i_atom - 1) * 3 + 3, i_mode)) / angstrom
          enddo
       end do
1000   format(a3,3(1x,f10.6))
       call returnclose_iounit(io_xmol)
       end do
  end subroutine movie_maker
  !*************************************************************


  !*************************************************************
  subroutine intensity_calc(n_internal,eig_int,step_size,intensities)
    ! Purpose: calculate the intensities of the normal modes
    !          I = prefactor * (dmu/dz)**2, where the derivative has
    !          to be taken along the normal modes
    !
    !----------------------------------------------------
   use opt_data_module , only:&
         & single_step
    implicit none
    ! --- Declaration of formal parameters --------------
    integer(kind=i4_kind), intent(in)  :: n_internal ! number of modes
    real(kind=r8_kind), intent(in)     :: step_size
    real(kind=r8_kind), intent(in)     :: eig_int(n_internal,n_internal) ! normal
    ! modes in internal coordinates
    real(kind=r8_kind), intent(out)     :: intensities(n_internal)
    ! *** end of interface ***

    ! --- Declaration of local variables ----------------
    real(kind=r8_kind) :: dmudq(3,n_internal), dipole_eq(3), &
         dipole_dis(3), dipole_dis2(3), &
         dq ! internal displacements
    integer(kind=i4_kind) :: i, i_mode,io_dip=11
    ! --- Executable code -------------------------------
    ! first we have to calculate the derivative of the dipolmoment
    ! with respect to the internal coordinates. This is done via finite
    ! differences. The dipole moments for the displacements are stored
    ! in the file dipole.dat
    ! first open file
    io_dip=openget_iounit(status='old',form='formatted',&
            file="dipole.dat")
    ! now we read the dipole moment for the equilibrium
    read(io_dip,*) dipole_eq
    ! now we read the dipolmoment for all displacements and calculate the
    ! derivatives with respect to internal coordinates dmudq
    dq = step_size ! fixed number from opt_data_module
    do i_mode=1,n_internal
       read(io_dip,*) dipole_dis
       if(single_step) then
          ! only one displacement per coordinate
          dmudq(:,i_mode)=(dipole_dis-dipole_eq)/dq
       else
          ! two displacements per coordinate
          read(io_dip,*) dipole_dis2
          dmudq(:,i_mode)=(dipole_dis2-dipole_dis)/(2*dq)
       end if
    end do
    ! to get the derivatives with respect to the normal modes we have to multiply
    ! with the eigenvectors , which contain the composition of the normal modes
    ! from the internal coordinates
    do i=1,3
       dmudq(i,:)=matmul(transpose(eig_int),dmudq(i,:))
    end do

    do i_mode=1,n_internal
       intensities(i_mode)=sum(dmudq(:,i_mode)**2)*975.46
    end do
    call returnclose_iounit(io_dip)
  end subroutine intensity_calc
  !*****************************************************************************

  !*****************************************************************************
  subroutine dipole_der(n_mode,dipder,error)
   use opt_data_module , only:&
         & single_step,step_size
    implicit none
    ! --- Declaration of formal parameters --------------
    integer(kind=i4_kind), intent(in)  :: n_mode ! number of modes
    real(kind=r8_kind), intent(out)    :: dipder(3,3*n_mode)
    logical, intent(out)               :: error
    logical                            :: exist
    ! *** end of interface ***

    ! --- Declaration of local variables ----------------
    real(kind=r8_kind) :: dipole_eq(3), &
         dipole_dis(3), dipole_dis2(3), &
         dq ! internal displacements
    integer(kind=i4_kind) :: i_mode,io_dip=11
    ! --- Executable code -------------------------------

    inquire (EXIST = exist,file="dipole.dat")
    if (.not.exist) then
      print*,'hesse_module: write_cart_hess: dipole_der is not present'
      error = .true.
      RETURN
    endif
    error = .false.

    io_dip=openget_iounit(status='old',form='formatted',&
            file="dipole.dat")
    ! now we read the dipole moment for the equilibrium
    read(io_dip,*) dipole_eq
    ! now we read the dipolmoment for all displacements and calculate the
    ! derivatives with respect to internal coordinates dmudq
    dq = step_size ! fixed number from opt_data_module
    do i_mode=1,n_mode
       read(io_dip,*) dipole_dis
       if(single_step) then
          ! only one displacement per coordinate
          dipder(:,i_mode)=(dipole_dis-dipole_eq)/dq
       else
          ! two displacements per coordinate
          read(io_dip,*) dipole_dis2
          dipder(:,i_mode)=(dipole_dis2-dipole_dis)/(2*dq)
       end if
    end do
    call returnclose_iounit(io_dip)
  end subroutine dipole_der

end module frequency_module
