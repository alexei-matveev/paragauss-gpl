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
module  spectrum_module
!---------------------------------------------------------------
!
!  Purpose: calculates spectra of hermitean and symmetric matrices
!           serves only for debugging reasons
!
!
!  Module called by: ...
!
!
!  References: ...
! 
!
!  Author: MM
!  Date: 10/97
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

use type_module        ! type specification parameters
use iounitadmin_module ! i/o routines
implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================


!------------ public functions and subroutines ------------------
public spectrum_complex,spectrum_real,spectrum_complex_g,spectrum_real_g,&
     spectrum_coeff_complex,spectrum_coeff_real

!================================================================
! End of public interface of module
!================================================================

  interface rsg
     SUBROUTINE RSG(NM,N,A,B,W,MATZ,Z,IERR)
       use type_module
       INTEGER(kind=i4_kind) N,NM,IERR,MATZ
       REAL(kind=r8_kind) A(:,:),B(:,:),W(:),Z(:,:)
     end SUBROUTINE RSG
  end interface




!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


   !*************************************************************
   subroutine spectrum_complex(matrix_real,matrix_imag,dim_matrix)
     ! purpose: determine spectrum of given matrix
     integer(kind=i4_kind),intent(in) :: dim_matrix
     real(kind=r8_kind),intent(in )   :: matrix_real(dim_matrix,dim_matrix),matrix_imag(dim_matrix,dim_matrix)
     !** End of interface *****************************************
     real(kind=r8_kind),allocatable   :: spectrum(:)
     real(kind=r8_kind),allocatable   :: matrix_real_save(:,:),matrix_imag_save(:,:)
     real(kind=r8_kind),allocatable   :: U_real(:,:),U_imag(:,:),work1(:),work2(:),work3(:,:)
     integer(kind=i4_kind)            :: ierr,i
     !------------ Executable code ------------------------------ 
     allocate(matrix_real_save(dim_matrix,dim_matrix),matrix_imag_save(dim_matrix,dim_matrix))
     allocate(work1(dim_matrix),work2(dim_matrix),work3(2,dim_matrix))
     allocate(U_real(dim_matrix,dim_matrix),U_imag(dim_matrix,dim_matrix),spectrum(dim_matrix))
     matrix_real_save = matrix_real
     matrix_imag_save = matrix_imag

       !
       ! now diagonalize matrix
       ! U(+)*T*U = t
       !
       !? print*,"spectrum_complex: dimension of Matrix: ",dim_matrix
       call ch(dim_matrix,dim_matrix,&                ! dimensions of transformation
            matrix_real_save,matrix_imag_save,&       ! matrix of eigenvalue problem
            spectrum,1,U_real,U_imag,&                ! eigenvalues and eigenvectors
            work1,work2,work3,ierr)                   ! work arrays
       if(ierr.ne.0) call error_handler &
            ("relativistic_trafo_master: error in eigensolver for T-matrix")

       ! show obtained spectrum
       do i=1,dim_matrix
          write(output_unit,*) " eigenvalue (",i,"):",spectrum(i)
       enddo
     
       deallocate(matrix_real_save,matrix_imag_save,U_real,U_imag,spectrum)
       deallocate(work1,work2,work3)
     end subroutine spectrum_complex
   !*************************************************************

   !*************************************************************
   subroutine spectrum_real(matrix,dim_matrix)
     ! purpose: determine spectrum of given matrix
     integer(kind=i4_kind),intent(in) :: dim_matrix
     real(kind=r8_kind),intent(in )   :: matrix(:,:)
     !** End of interface *****************************************
     real(kind=r8_kind),allocatable   :: spectrum(:)
     real(kind=r8_kind),allocatable   :: matrix_save(:,:)
     real(kind=r8_kind),allocatable   :: U_matrix(:,:),work1(:),work2(:)
     integer(kind=i4_kind)            :: ierr,i
     !------------ Executable code ------------------------------ 
     allocate(matrix_save(dim_matrix,dim_matrix))
     allocate(work1(dim_matrix),work2(dim_matrix))
     allocate(U_matrix(dim_matrix,dim_matrix),spectrum(dim_matrix))
     matrix_save = matrix
       !
       ! now diagonalize matrix
       ! U(+)*T*U = t
       !
       call rs(dim_matrix,dim_matrix,&                ! dimensions of transformation
            matrix_save,&                             ! matrix of eigenvalue problem
            spectrum,1,U_matrix,&                     ! eigenvalues and eigenvectors
            work1,work2,ierr)                         ! work arrays
       if(ierr.ne.0) call error_handler &
            ("relativistic_trafo_master: error in eigensolver for T-matrix")

       ! show obtained spectrum
       do i=1,dim_matrix
          write(output_unit,*) " eigenvalue (",i,"):",spectrum(i)
       enddo
     
       deallocate(matrix_save,U_matrix,spectrum)
       deallocate(work1,work2)
     end subroutine spectrum_real
   !*************************************************************

   !*************************************************************
   subroutine spectrum_real_g(matrix,S_matrix,dim_matrix)
     ! purpose: determine spectrum of given matrix
     integer(kind=i4_kind),intent(in) :: dim_matrix
     real(kind=r8_kind),intent(in )   :: matrix(:,:)
     real(kind=r8_kind),intent(in )   :: S_matrix(:,:)
     !** End of interface *****************************************
     real(kind=r8_kind),allocatable   :: spectrum(:)
     real(kind=r8_kind),allocatable   :: matrix_save(:,:)
     real(kind=r8_kind),allocatable   :: U_matrix(:,:)
     integer(kind=i4_kind)            :: ierr,i
     !------------ Executable code ------------------------------ 
     allocate(matrix_save(dim_matrix,dim_matrix))
     allocate(U_matrix(dim_matrix,dim_matrix),spectrum(dim_matrix))
     matrix_save = matrix
       !
       ! now diagonalize matrix
       ! U(+)*T*U = t
       !
       call rsg(dim_matrix,dim_matrix,&               ! dimensions of transformation
            matrix_save,&                             ! matrix of eigenvalue problem
            S_matrix,&                                ! Overlap Matrix
            spectrum,1,U_matrix,ierr)                 ! eigenvalues and eigenvectors

       if(ierr.ne.0) call error_handler &
            ("relativistic_trafo_master: error in eigensolver for T-matrix")

       ! show obtained spectrum
       do i=1,dim_matrix
          write(output_unit,*) " eigenvalue (",i,"):",spectrum(i)
       enddo
     
       deallocate(matrix_save,U_matrix,spectrum)
     end subroutine spectrum_real_g
   !*************************************************************

   !*************************************************************
     subroutine spectrum_complex_g(n,H_real,H_imag,S_real,S_imag)
    !  Purpose: Solves generalized hermitean definite Eigenvalueproblem
    !           H*X = S*X*eigen_value
    !           where H,S are hermitean and S additionally is positive definite
    !  
    !
    !  Author: MM
    !  Date: 10/97
    !
    !------------ Modules --------------------------------------
    use iounitadmin_module
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in)     :: n
    real(kind=r8_kind),intent(in)        :: H_real(n,n),H_imag(n,n)
    real(kind=r8_kind),intent(in)        :: S_real(n,n),S_imag(n,n)
    !** End of interface ***************************************
    !------------ Declaration of local parameters ---------------
    real(kind=r8_kind),parameter                   :: small = 0.00001
    ! very small value
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind)                 :: ierr, i
    real(kind=r8_kind),allocatable        :: eigen_value(:)
    real(kind=r8_kind),allocatable        :: X_real(:,:),X_imag(:,:)
    real(kind=r8_kind),allocatable        :: work1(:),work2(:),work3(:,:)
    real(kind=r8_kind),allocatable        :: s_diag(:),U_real(:,:),U_imag(:,:)
    real(kind=r8_kind),allocatable        :: H_red_real(:,:),H_red_imag(:,:)
    real(kind=r8_kind),allocatable        :: A_real(:,:),A_imag(:,:)
    !------------ Executable code ------------------------------

    !? print*,">>> spectrum_complex_g: entered"
    ! allocate local arrays
    allocate(U_real(n,n),U_imag(n,n),s_diag(n),H_red_real(n,n),H_red_imag(n,n),eigen_value(n),&
         X_real(n,n),X_imag(n,n),A_real(n,n),A_imag(n,n))
    
    ! Diagonalize S
    !? print*,">>> spectrum_complex_g: Diagonalize S"
    ierr = 0
    allocate(work1(n),work2(n),work3(2,n))
    call ch(n,n,&                                     ! dimensions of eigenvalueproblem
         S_real,S_imag,&                              ! matrix of eigenvalue problem
         s_diag,1,U_real,U_imag,&                     ! eigenvalues and eigenvectors
         work1,work2,work3,ierr)                      ! work arrays
    if (ierr.ne.0) then
       write(*, *) "eigen_data_chg: stopped due to error in eigensolver IERR =",ierr
       call error_handler("eigen_data_chg: error in eigensolver")
    end if
    deallocate(work1,work2,work3)
    
    ! Construct transformation matrix Q
    ! Q = U * s^(-1/2)
    do i=1,n
       U_real(:,i) = U_real(:,i)*1/sqrt(s_diag(i))
       U_imag(:,i) = U_imag(:,i)*1/sqrt(s_diag(i))
    end do

    ! determine reduced Hamiltonian H_red
    ! H_red = Q(+)*H*Q
    H_red_real = H_real
    H_red_imag = H_imag
    A_real = U_real
    A_imag = U_imag
    call sim_trafo(U_real,U_imag,H_red_real,H_red_imag,n)
    U_real = A_real
    U_imag = A_imag

    ! Diagonalize H_red
    !? print*,">>> spectrum_complex_g: Diagonalize H_red"
    ierr = 0
    A_real = H_red_real
    A_imag = H_red_imag
    allocate(work1(n),work2(n),work3(2,n))
    call ch(n,n,&                                     ! dimensions of eigenvalueproblem
         H_red_real,H_red_imag,&                      ! matrix of eigenvalue problem
         eigen_value,1,X_real,X_imag,&                ! eigenvalues and eigenvectors
         work1,work2,work3,ierr)                      ! work arrays
    if (ierr.ne.0) then
       write(*, *) "eigen_data_chg: stopped due to error in eigensolver IERR =",ierr
       call error_handler("eigen_data_chg: error in eigensolver")
    end if
    !? print*,ierr
    deallocate(work1,work2,work3)

    ! show spectrum 
    do i = 1,n
       write(output_unit,*) "eigenvalue (",i,"):",eigen_value(i)
    enddo
    H_red_real = A_real
    H_red_imag = A_imag

    ! test eigenvectors
    !A_real = matmul(H_red_real,X_real) - matmul(H_red_imag,X_imag)
    !A_imag = matmul(H_red_imag,X_real) + matmul(H_red_real,X_imag)
    !write(output_unit,*) ">> Spectrum of solutions (before back trafo) <<"
    ! show spectrum 
    !do i = 1,n
    !   do j =1,n
    !      if (A_real(j,i).gt.small) then
    !         write(output_unit,*) "eigenvalue (",i,"):",A_real(j,i)/X_real(j,i)
    !     endif
    !   enddo
    !enddo
    
    ! now transform back Eigenvectors
    !H_red_real = X_real
    !H_red_imag = X_imag
    !X_real = matmul(U_real,H_red_real) - matmul(U_imag,H_red_imag)
    !X_imag = matmul(U_imag,H_red_real) + matmul(U_real,H_red_imag)
    
    ! test eigenvectors
    !H_red_real = matmul(H_real,X_real) - matmul(H_imag,X_imag)
    !H_red_imag = matmul(H_imag,X_real) + matmul(H_real,X_imag)
    !U_real = matmul(S_real,X_real) - matmul(S_imag,X_imag)
    !U_imag = matmul(S_imag,X_real) + matmul(S_real,X_imag)
    !write(output_unit,*) ">> Spectrum of solutions (after backtrafo)<<"
    ! show spectrum 
    !do i = 1,n
    !   do j =1,n
    !      if (U_real(j,i).gt.small) then
    !         write(output_unit,*) "eigenvalue real(",i,"):",H_red_real(j,i)/U_real(j,i)
    !      endif
    !      if (U_imag(j,i).gt.small) then
    !         write(output_unit,*) "eigenvalue imag(",i,"):",H_red_imag(j,i)/U_imag(j,i)
    !      endif
    !   enddo
    !enddo
    ! deallocate local arrays
    deallocate(U_real,U_imag,s_diag,H_red_real,H_red_imag,eigen_value,X_real,X_imag,A_real,A_imag)
    
  contains
    
    !*************************************************************
    subroutine sim_trafo(U_real,U_imag,V_real,V_imag,dimension)
      ! purpose: determine similarity transformation U(+)VU
      integer(kind=i4_kind) :: dimension
      real(kind=r8_kind) :: U_real(dimension,dimension),U_imag(dimension,dimension),&
           V_real(dimension,dimension),V_imag(dimension,dimension),&
           W_real(dimension,dimension),W_imag(dimension,dimension)
      !------------ Executable code ------------------------------ 
      W_real = matmul(transpose(U_real),matmul(V_real,U_real)) +&
           matmul(transpose(U_imag),matmul(V_real,U_imag)) +&
           matmul(transpose(U_imag),matmul(V_imag,U_real)) -&
           matmul(transpose(U_real),matmul(V_imag,U_imag)) 
      W_imag = matmul(transpose(U_real),matmul(V_imag,U_real)) +&
           matmul(transpose(U_imag),matmul(V_imag,U_imag)) -&
           matmul(transpose(U_imag),matmul(V_real,U_real)) +&
           matmul(transpose(U_real),matmul(V_real,U_imag))
      V_real = W_real
      V_imag = W_imag
    end subroutine sim_trafo
    !*************************************************************

  end subroutine spectrum_complex_g
  !*************************************************************

  !*************************************************************
  subroutine spectrum_coeff_complex(orbital_real,orbital_imag,irrep)
    !------------ Modules --------------------------------------
    use iounitadmin_module
    use unique_atom_module
    use symmetry_data_module
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(in)        :: orbital_real(:),orbital_imag(:)
    integer,intent(in)                   :: irrep
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind)                :: i_sa, n
    integer(kind=i4_kind)                :: unique, l, if, cf
    integer(kind=i4_kind)                :: equal,i_spinor,i_pa
    type(unique_atom_type),          pointer :: ua
    type(unique_atom_basis_type),    pointer :: uab
    type(unique_atom_partner_type),  pointer :: sp
    type(unique_atom_sa_int_type),   pointer :: sa_int
    real(kind=r8_kind),allocatable           :: coefficient_real(:,:,:,:),coefficient_imag(:,:,:,:)
!!$    type(unique_atom_renorm_indfct_type), pointer :: r
    real(kind=r8_kind)                        :: sum_spinor !!$ norm,sum_spinor
    ! coefficient(equal,m,contraction,i_spinor)
    !------------ Executable code ------------------------------

    i_pa = 1
    do i_pa =1,symmetry_data_n_partners_proj(irrep)
       i_sa = 1
       write(output_unit,*) "PARTNER: ",i_pa
       do unique =1,N_unique_atoms
          ua=>unique_atoms(unique)
          do l =0,ua%lmax_ob ! loop angular momenta
             uab=>ua%l_ob(l)
             ! allocate coefficients
             allocate(coefficient_real(ua%N_equal_atoms,2*l+1,uab%N_exponents,2),&
                  coefficient_imag(ua%N_equal_atoms,2*l+1,uab%N_exponents,2))
             coefficient_real = 0.0_r8_kind
             coefficient_imag = 0.0_r8_kind
                write(output_unit,*) "*****************************************************"
                write(output_unit,*) "*        l=",l
                write(output_unit,*) "*****************************************************"
                sp=>ua%symadapt_spor_partner(irrep,l)
!!$                r=>ua%renormaliation_spor_partner_ob(irrep,l)
                do if=1,sp%N_independent_fcts
                   write(output_unit,*) "*** INDEPENDENT ***",if
                   do cf=1,uab%N_uncontracted_fcts ! uncontracted exponents
!!$                      norm = r%renorm(if)%c_exp(cf)
                      ! now add up contributions
                      do i_spinor =1,2
                         write(output_unit,*) " spinor: ",i_spinor
                         do equal=1,ua%N_equal_atoms
                            sa_int=>sp%sa_spor_int(equal,i_spinor,if,i_pa)
                            write(output_unit,*) " equal: ",equal
                            do n=1,sa_int%N_fcts
                               coefficient_real(equal,sa_int%m(n),cf,i_spinor)=&
                                    coefficient_real(equal,sa_int%m(n),cf,i_spinor)&
                                    +(sa_int%re(n)*orbital_real(i_sa) - sa_int%im(n)*orbital_imag(i_sa)) !!$*norm 
                               coefficient_imag(equal,sa_int%m(n),cf,i_spinor)=&
                                    coefficient_imag(equal,sa_int%m(n),cf,i_spinor)&
                                    +(sa_int%im(n)*orbital_real(i_sa)+sa_int%re(n)*orbital_imag(i_sa)) !!$*norm
                               write(output_unit,*) "symmadapt m:",sa_int%m(n)," ",sa_int%re(n),sa_int%im(n)
                            enddo
                         enddo
                      enddo
                      i_sa = i_sa + 1
                   enddo
                enddo
             ! show coefficients
             write(output_unit,*) "**************************************************"
             write(output_unit,*) "*         unique:",unique
             write(output_unit,*) "*              l:",l
             write(output_unit,*) "**************************************************"
             do equal=1,ua%N_equal_atoms
                do cf=1,uab%N_exponents
                   do n=1,2*l+1
                      write(output_unit,*) " equal: ",equal,"m: ",n
                      sum_spinor = 0.0_r8_kind
                      do i_spinor=1,2
                         write(output_unit,*)"  spinor(",i_spinor,"): ",coefficient_real(equal,n,cf,i_spinor),&
                              coefficient_imag(equal,n,cf,i_spinor)
                         sum_spinor = sum_spinor + coefficient_real(equal,n,cf,i_spinor)**2 +&
                              coefficient_imag(equal,n,cf,i_spinor)**2
                      enddo
                      write(output_unit,*) " total:",sqrt(sum_spinor)
                   enddo
                enddo
             enddo
             ! deallocate coefficients
             deallocate(coefficient_real,coefficient_imag)
          enddo
       enddo
    enddo
    
    
  end subroutine spectrum_coeff_complex
  !*************************************************************

  !*************************************************************
  subroutine spectrum_coeff_real(orbital,irrep)
    !------------ Modules --------------------------------------
    use iounitadmin_module
    use unique_atom_module
    use orbitalprojection_module
    use symmetry_data_module
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(in)        :: orbital(:)
    integer,intent(in)                   :: irrep
    !** End of interface ***************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind)                :: i_sa, n
    integer(kind=i4_kind)                :: unique, l, if, cf
    integer(kind=i4_kind)                :: equal,i_pa
    type(unique_atom_type),          pointer :: ua
    type(unique_atom_basis_type),    pointer :: uab
    type(unique_atom_partner_type),  pointer :: sp
    type(unique_atom_symadapt_type),  pointer :: sadapt
!!$    type(unique_atom_renorm_indfct_type), pointer :: r
    real(kind=r8_kind),allocatable            :: coefficient(:,:,:)
!!$    real(kind=r8_kind)                        :: norm
    ! coefficient(equal,m,contraction)
    !------------ Executable code ------------------------------

    i_sa = 1
    i_pa = 1
    do i_pa = 1,symmetry_data_n_partners(irrep)
       write(output_unit,*) "PARTNER: ",i_pa
       do unique =1,N_unique_atoms
          ua=>unique_atoms(unique)
          do l =0,ua%lmax_ob ! loop angular momenta
             uab=>ua%l_ob(l)
             i_sa = orbitalprojection_ob(irrep,l,unique)
             sp=>ua%symadapt_partner(irrep,l)
!!$             r=>ua%renormaliation_partner_ob(irrep,l)
             ! allocate coefficients
             allocate(coefficient(ua%N_equal_atoms,2*l+1,uab%N_exponents))
             coefficient = 0.0_r8_kind
             do if=1,sp%N_independent_fcts
                write(output_unit,*) "independent function #: ",if
                do cf=1,uab%N_uncontracted_fcts ! uncontracted exponents
                   write(output_unit,*) "exponent #",cf
!!$                   norm = r%renorm(if)%c_exp(cf)
                   ! now add up contributions
                   sadapt=>sp%symadapt(if,i_pa)
                   do n=1,sadapt%N_fcts
                      coefficient(sadapt%I_equal_atom(n),sadapt%m(n),cf)=&
                           coefficient(sadapt%I_equal_atom(n),sadapt%m(n),cf)+sadapt%c(n)*orbital(i_sa) !!$*norm
                      write(output_unit,*) "m: ",sadapt%m(n),"equal: ",sadapt%I_equal_atom(n)
                      write(output_unit,*) "symmadapt (",n,")",sadapt%c(n)
                   enddo
                   write(output_unit,*) "orbital(",i_sa,"): ",orbital(i_sa)
!!$                   write(output_unit,*) "norm:              ",norm
                   i_sa = i_sa + 1
                enddo! exponents
             enddo
             ! show coefficients
             write(output_unit,*) "**************************************************"
             write(output_unit,*) "*         unique:",unique
             write(output_unit,*) "*              l:",l
             write(output_unit,*) "**************************************************"
             do equal=1,ua%N_equal_atoms
                do cf=1,uab%N_exponents
                   do n=1,2*l+1
                      write(output_unit,*) " equal: ",equal,"m: ",n," coefficient: ",coefficient(equal,n,cf)
                   enddo
                enddo
             enddo
             ! deallocate coefficients
             deallocate(coefficient)
           enddo
       enddo
    enddo
    
    
  end subroutine spectrum_coeff_real
  !*************************************************************

!--------------- End of module ----------------------------------
end module spectrum_module
