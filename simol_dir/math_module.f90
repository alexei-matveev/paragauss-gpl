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
module math_module
  !---------------------------------------------------------------
  !
  !  Purpose: module contains subroutines for various 
  !           mathematical tasks
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  use type_module ! type specification parameters
  implicit none
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public cross_product,abs_value,rotmat,invert_matrix, &
       dot_prod,binomial_coeff,schmidt,round,eigvec_sign_check,&
       print_matrix,equal_vector,ortho2
  !     print_matrix,ortho2

  !================================================================
  ! End of public interface of module
  !================================================================


  real(kind=r8_kind),parameter,public :: one=1.0_r8_kind,zero=0.0_r8_kind, &
       two=2.0_r8_kind,three=3.0_r8_kind,&
       four=4.0_r8_kind,five=5.0_r8_kind,&
       half=0.5_r8_kind
  real(kind=r8_kind),parameter,public :: small=1.0e-10_r8_kind
  real(kind=r8_kind),parameter,public :: pi=3.14159265358979324_r8_kind, &
       convert1=180.0_r8_kind/pi, &
       convert2=0.52917706

contains

  function cross_product(vec1,vec2)
    real(kind=r8_kind),intent(in)   :: vec1(3),vec2(3)
    real(kind=r8_kind)              :: cross_product(3)
    !** End of interface *****************************************

    cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  end function cross_product

  function dot_prod(vec1,vec2)
    real(kind=r8_kind),intent(in)   :: vec1(3),vec2(3)
    real(kind=r8_kind)              :: dot_prod
    !** End of interface *****************************************

    dot_prod = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
  end function dot_prod


  function abs_value(vec1)
    real(kind=r8_kind)             :: abs_value
    real(kind=r8_kind),intent(in)  :: vec1(:)
    !** End of interface *****************************************

    abs_value = sqrt(sum(vec1**2))
  end function abs_value

  function binomial_coeff(n,k)
    real(kind=r8_kind)     :: binomial_coeff
    integer(kind=i4_kind)  :: n,k
    integer(kind=i4_kind)  :: ik,k_fac,nom
    !** End of interface *****************************************

    k_fac=1
    do ik=1,k
       k_fac=k_fac*ik
    enddo
    nom = 1
    do ik=0,k-1
       nom = nom*(n-ik)
    enddo
    binomial_coeff = nom/k_fac
  end function binomial_coeff

  function rotmat(axis,angle)
    integer(kind=i4_kind)  :: axis ! 1=x,2=y,3=z
    real(kind=r8_kind)     :: angle
    real(kind=r8_kind)     :: rotmat(3,3)
    !** End of interface *****************************************

    real(kind=r8_kind)     :: c_a,s_a

    c_a = cos(angle)
    s_a = sin(angle)

    if (axis == 1 ) then
       rotmat(1,1) = one
       rotmat(1,2) = zero
       rotmat(1,3) = zero
       rotmat(2,1) = zero
       rotmat(2,2) = c_a
       rotmat(2,3) = -s_a
       rotmat(3,1) = zero
       rotmat(3,2) = s_a
       rotmat(3,3) = c_a
       elseif (axis == 2 ) then
       rotmat(1,1) = c_a
       rotmat(1,2) = zero
       rotmat(1,3) = -s_a
       rotmat(2,1) = zero
       rotmat(2,2) = one
       rotmat(2,3) = zero
       rotmat(3,1) = c_a
       rotmat(3,2) = zero
       rotmat(3,3) = s_a
       elseif ( axis == 3) then
       rotmat(1,1) = c_a
       rotmat(1,2) = -s_a
       rotmat(1,3) = zero
       rotmat(2,1) = s_a
       rotmat(2,2) = c_a
       rotmat(2,3) = zero
       rotmat(3,1) = zero
       rotmat(3,2) = zero
       rotmat(3,3) = one
    else
       call error_handler(' rotmat : stupid axis specified')
    endif

  end function rotmat

  function round(number,digit)
    real(kind=r8_kind)     :: number,round
    integer(kind=i4_kind)  :: digit
    !** End of interface *****************************************
    !-----------------------
    real(kind=r8_kind) :: hilf
    integer(kind=i4_kind)  :: hilf_int
    if (number<=zero) then
       hilf = number*10.0_r8_kind**(real(digit))-half
    else
       hilf = number*10.0_r8_kind**(real(digit))+half
    endif
    hilf_int = int(hilf,kind=i4_kind)
    round=real(hilf_int,kind=r8_kind)/10.0_r8_kind**(real(digit))
  end function round

  function equal_vector(vec1,vec2)
    logical            :: equal_vector
    real(kind=r8_kind) :: vec1(:),vec2(:)
    !** End of interface *****************************************
    ! ----------------------------------
    integer(kind=i4_kind)  :: dim1,dim2,summ,i

    equal_vector=.false.
    dim1=ubound(vec1,1)
    dim2=ubound(vec2,1)
    if (dim1/=dim2) return
    summ=0
    do i=1,dim1
       if ( abs(vec1(i)-vec2(i))<small ) summ=summ+1
    enddo
    if (summ==dim1) then
       equal_vector=.true.
    endif
    return
  end function equal_vector

  subroutine invert_matrix(dimen,matrix)
    ! Purpose : wrapper for the dge*-LAPACK Routines for
    !          the special case of quadratic matrices (i.e. dim1=dim2).
    !          Up to now only the double-precision general-matrix
    !          routines are implemented, although  for the Hessian
    !          which is supposed to be symmetric (and positive definite
    !          if a minimum is searched) other routines would
    !          be appropriate. This has to be chagned lateron.
    !-------------------------------------------------------
    integer(kind=i4_kind),intent(in)  :: dimen
    real(kind=r8_kind),intent(inout)  :: matrix(:,:)
    !** End of interface *****************************************
    ! ------- Declaration of local variables ---------------
    integer(kind=i4_kind)    :: alloc_stat,info
    integer(kind=i4_kind),allocatable :: ipiv(:)
    real(kind=r8_kind),allocatable    :: work(:)
    ! ------ Executable code -------------------------------
    allocate(ipiv(dimen),STAT=alloc_stat)
    if (alloc_stat/=0) &
         call error_handler(' invert_matrix : allocation (1) failed')

    allocate(work(dimen),STAT=alloc_stat)
    if (alloc_stat/=0) &
         call error_handler(' invert_matrix : allocation (2) failed')
    info=0
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(dimen,dimen,matrix,dimen,ipiv,info)
    if (info < 0 ) then
       call error_handler&
            ("invert_matrix: dgetrf exited with an illegal value as argument")
    elseif (info > 0) then
       write(*,*)" dgetrf:  element ",info," is exactly zero in the LU-factorization"
       write(*,*)"          input matrix is singular"
       call error_handler("invert_matrix: dgetrf exited")
    endif
    info=0

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call dgetri(dimen,matrix,dimen,ipiv,work,dimen,info)
    if (info < 0 ) then
       write(*,*)' invert_bmat : warning - dgetri exited with info < 0'
       write(*,*)" argument ",-info," had an illegal value"
       call error_handler(" invert_matrix: exit")
    elseif (info > 0 ) then
       write(*,*)" invert_matrix: the matrix is singular and cannot be inverted"
       call error_handler(" invert_matrix: exit")
    endif
    deallocate(work,ipiv,STAT=alloc_stat)
    if (alloc_stat/=0) &
         call error_handler('invert_matrix : deallocation (1) failed')

  end subroutine invert_matrix

  !*************************************************************
  subroutine eigvec_sign_check(dimen,eigvec,eigvec_prev)
    ! Purpose: compare the signs of the previous eigenvectors
    !          to the actual ones. If  ALL signs on an eigenvector
    !          have changed, reset them to those of the previous
    !          eigenvector. 
    !          This is necessary for the backtransformation from
    !          delocalized internal coordinates to cartesians
    !          as well as for the comparison of internals and
    !          gradients between two geometry steps.
    ! 
    ! Routine called by: generate_bmat_reduced
    ! --------------------------------------------------------
    integer(kind=i4_kind),intent(in)    :: dimen
    real(kind=r8_kind),intent(inout)    :: eigvec(:,:)
    real(kind=r8_kind),intent(in)       :: eigvec_prev(:,:)
    !** End of interface *****************************************
    ! --- Declaration of local variables ---------------------
    real(kind=r8_kind)    :: dir_1(dimen),dir_2(dimen),c_phi
    integer(kind=i4_kind) :: i

    dir_1=zero
    dir_2=zero
    do i=1,dimen
       where (eigvec(:,i)>zero) dir_1 = one
       where ( eigvec(:,i)<zero) dir_1 = -one
       where ( abs(eigvec(:,i))<small) dir_1 = zero
       where (eigvec_prev(:,i)>zero) dir_2 = one
       where ( eigvec_prev(:,i)<zero) dir_2 = -one
       where ( abs(eigvec_prev(:,i))<small) dir_2 = zero
       
       dir_1 = dir_1/abs_value(dir_1)
       dir_2 = dir_2/abs_value(dir_2)
       c_phi = dot_product(dir_1,dir_2)
       if (abs(abs(c_phi)-one)<small .and. c_phi<zero ) then
          eigvec(:,i) = -eigvec(:,i)
       endif
    enddo
  end subroutine eigvec_sign_check

  !*************************************************************

  subroutine schmidt(n_cons,c_cons,vec_in,vec_out)
    ! Purpose: Perform a Schmidt-Orthogonalization using
    !          <n_cons> -already orthogonal- starting vectors
    !          <c_cons>.
    !          Care has to be taken for the case that one of the 
    !          projected constraint vectors is equal to one
    !          of the vectors in 'vec_in'. In that case the
    !          vector has to be skipped in order to prevent its
    !          orthogonalization to zero.
    ! --------------------------------------------------
    integer(kind=i4_kind),intent(in) :: n_cons
    real(kind=r8_kind),intent(in)    :: c_cons(:,:)
    real(kind=r8_kind),intent(in)    :: vec_in(:,:)
    real(kind=r8_kind),intent(out)   :: vec_out(:,:)
    !** End of interface *****************************************
    ! ---------------------------------------------------------
    real(kind=r8_kind),parameter   :: small_ortho=1.0e-6_r8_kind
    integer(kind=i4_kind)          :: n_prim,n_int,alloc_stat,&
         i_vec,i_help,i_cons,counter,ii,n_reject
    real(kind=r8_kind),allocatable :: help_vec(:)
    real(kind=r8_kind)             :: test_prod
    logical,allocatable            :: list(:)

    n_prim = ubound(vec_in,1)
    n_int = ubound(vec_in,2)
    if (ubound(c_cons,1)/=n_prim) call error_handler &
         (" Schmidt: error in dimensions")

    vec_out=zero
    allocate(help_vec(n_prim),list(n_int),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("Schmidt: allocation (1) failed")
    ! first generate the projected constraint vectors.
    do i_cons=1,n_cons 
       help_vec=zero
       do i_help=1,n_int
          help_vec = help_vec + & 
               dot_product(c_cons(:,i_cons),vec_in(:,i_help))*&
               vec_in(:,i_help)
       enddo
       vec_out(:,i_cons) = help_vec/abs_value(help_vec)
    enddo

    
    write(*,*)" The projected constraint vectors are:"
    do i_vec=1,n_prim
       write(*,'(6(3x,f9. 5))')(vec_out(i_vec,i_cons),i_cons=1,n_cons)
    enddo
    write(*,*)" "
    do i_vec=1,n_cons
       do ii=1,n_cons
          if (i_vec==ii) cycle
          test_prod = dot_product(vec_out(:,i_vec),vec_out(:,ii))
          if (abs(test_prod)>=small_ortho) then
             write(*,*)" schmidt: after projection onto the active subspace"
             write(*,*)"          the constraint vectors are no longer orthogonal."
             write(*,*)"          This means that your constraints imply further "
             write(*,*)"          constraints, that are not specified.            "
             write(*,*)"          Constraint vectors ",i_vec," and ",ii," are no longer orthogonal"
             !call error_handler(" Please revise your constraints")
          endif
       enddo
    enddo
    ! set a logical mask to true for those vectors of vec_in which are
    ! NOT equal to any of the (projected) constraint vectors.
    list = .true.
    do i_cons =1,n_cons
       do i_vec=1,n_int
          help_vec = vec_out(:,i_cons) - vec_in(:,i_vec)
          if (abs_value(help_vec)<small) list(i_vec) = .false.
       enddo
    enddo
    if (any(.not.list)) then
       write(*,*)" schmidt: one of the original vectors in UMAT equals "
       write(*,*)"          one of the constraint vectors. Skip this   "
       write(*,*)"          one in the Schmidt-Orthogonalization"
       write(*,*)" "
    endif

    counter=1
!    counter=n_cons
    n_reject=0
    vectors: do i_vec = 1,n_int
       if (list(i_vec)) then
          counter=counter+1
       else
          cycle vectors
       endif
       if (counter>n_int) exit vectors
       help_vec = zero
       do i_help = 1,counter-1
          help_vec = help_vec + &
               dot_product(vec_in(:,i_vec),vec_out(:,i_help))*&
               vec_out(:,i_help)
       enddo
       vec_out(:,counter) = (vec_in(:,i_vec) - help_vec)
       if ( abs_value(vec_out(:,counter))<=small ) then
          write(*,*)" schmidt: a vector dropped out of the "
          write(*,*)"          orthogonalization scheme"
          write(*,*)" "
          counter=counter-1
          cycle vectors
       endif
       vec_out(:,counter) = vec_out(:,counter)/abs_value(vec_out(:,counter))
       do ii=1,counter-1
          ! test if the latest vector is orthogonal to the
          ! previous ones
          test_prod = dot_product(vec_out(:,counter),vec_out(:,ii))
          if (test_prod>=small_ortho) then
             write(*,*)" schmidt: a vector was found not to be "
             write(*,*)"          orthogonal to the previous ones."
             write(*,*)"          Reject vector ",counter
             counter=counter-1
             n_reject=n_reject+1
             cycle vectors
          endif
       enddo
    enddo vectors
    print*,' Counter is now                   :',counter
    print*,' Total number of vectors required :',n_int
    print*,' Number of constraint vectors     :',n_cons
    print*,' Number of rejected vectors       :',n_reject

!!$    ! Finally replace the projected constraint vectors by the
!!$    ! original ones.
!!$    do i_cons=1,n_cons
!!$       vec_out(:,i_cons) = c_cons(:,i_cons)
!!$    enddo
    deallocate(help_vec,list,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" Schmidt: deallocation (1) failed")
  end subroutine schmidt
    
  !*************************************************************
       
  !*************************************************************
  subroutine ortho2(c_cons,vec_in,vec_out)
    real(kind=r8_kind)   :: c_cons(:,:),vec_in(:,:),vec_out(:,:)
    ! ----------------------------------------------------------
    integer(kind=i4_kind)               :: n_constraint,n_internal,n_prim
    real(kind=r8_kind),allocatable      :: matrix(:,:),mat(:,:),&
                                           coeff(:,:),rhs(:),&
                                           work(:),help_vec(:),&
                                           eigval(:),eigvec(:,:),help_mat(:,:)
    integer(kind=i4_kind)               :: start,dimi,alloc_stat,&
                                           counter_j,counter_k,&
                                           i,j,k,counter_test,kk,&
                                           n_int,n_cons,counter_equal,&
                                           counter_cons
    integer(kind=i4_kind),allocatable   :: ind(:),ind_cons(:),ind_sort(:)
    integer                             :: m,n,lda,ldb,info,lwork,&
                                           nrhs
    real(kind=r8_kind)                  :: abs_val,minni
    real(kind=r8_kind),parameter        :: mat_test=1.0e5_r8_kind
    external dgels
    
    n_constraint = ubound(c_cons,2) ! Anzahl constraints
    n_internal = ubound(vec_in,2)  ! Anzahl Vektoren im active set
    n_prim = ubound(vec_in,1) ! Anzahl primitiver Koordinate

    allocate(ind_sort(n_constraint),ind_cons(n_constraint),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler("ortho2: allocation (1) failed")
    ind_cons=0_i4_kind
    ind_sort=0_i4_kind
    
    ! Stelle zuerst fest, ob es im active set Vektoren gibt, die mit einem
    ! der constraint vektoren identisch sind.
    ! 'ind_sort' : enthaelt Indizes der Eigenvektoren, fuer die es identische
    !              unter den Constraint Vektoren gibt.
    ! 'ind_cons' : enthaelt Indizes der Constraint Vektoren, zu denen es
    !              KEINEN identischen unter den Eigenvektoren gibt.
    counter_equal=0
    counter_cons=0
    do k=1,n_constraint
       do j=1,n_internal
          if (equal_vector(c_cons(:,k),vec_in(:,j))) then
             counter_equal = counter_equal + 1
             ind_sort(counter_equal)=j
          else
             counter_cons=counter_cons+1
             ind_cons(counter_cons)=k
          endif
       enddo
    enddo
    ! Dies ist die Anzahl von Basisvektoren, die zur Konstruktion
    ! des neuen active sets benoetigt werden. Gehoert einer der Vektoren des alten
    ! active set gleichzeitig zu den Constraints, so kann er vom Raum, aus dem der Rest
    ! des neuen active sets konstruiert werden soll, ausgeschlossen werden.
    n_int = n_internal-counter_equal
    ! Gleichzeitig muss er auch von den restlichen constraint vektoren ausgeschlossen
    ! werden.
    n_cons=n_constraint-counter_equal
    if (n_cons/=n_constraint) then
       write(*,*)" ortho2: Among the Constraint Vectors ",counter_equal," were found "
       write(*,*)"         to be identical to one of the vectors of the active set.   "
       write(*,*)" ortho2: "
       write(*,*)" This is a list of indices of those eigenvectors of the (original) "
       write(*,*)" active set for which an identical constraint vector has been found:"
       write(*,*)" ----  ",ind_sort(1:counter_equal)," ---- "
       write(*,*)" "
    endif

    allocate(ind(n_int),coeff(n_int,n_int-n_cons),help_vec(n_prim),STAT=alloc_stat)
    if(alloc_stat/=0)call error_handler("ortho2: allocation (2) failed")
    coeff=zero
    ind=0_i4_kind

    if (n_cons==0) then ! das heisst es gab nur einen constraint vektor
       !                 und der war auch noch identisch mit einem der eigenvektoren
       vec_out=vec_in
       return
    endif

    i=0 ! Zaehlindex der Vektoren v
    vector_v: do 
       i=i+1
       if (i > (n_int-n_cons) )exit vector_v  ! Wir brauchen n_int-n_cons
       !                                     ! Vektoren v
       ! Wie gross wird das zu loesende GLS sein?
       dimi = n_cons+i-1
       ! Die Koeffiziente von 1 ... start sollen frei gewaehlt werden,
       ! die von start+1 ... n_int sollen bestimmt werden:
       start = n_int-dimi

       allocate(matrix(dimi,n_int),mat(dimi,dimi),rhs(dimi),STAT=alloc_stat)
       if(alloc_stat/=0)call error_handler(" ortho2: allocation (3) failed")
       matrix=zero 
       
       ! Hier wird 'ind_sort' lediglich dazu benutzt, um aus
       ! dem vollen active set die Eigenvektoren herauszusuchen,
       ! die mit keinem der Constraint vektoren gleich sind.
       k=n_internal+1
       counter_j=0
       full_active: do j=1,n_internal
          k=k-1
          do kk=1,counter_equal
             if (k==ind_sort(kk)) then
                cycle full_active
             endif
          enddo
          counter_j=counter_j+1
          ind(counter_j)=k
       enddo full_active

       ! Um groesstmoegliche Aehnlichkeit mit den input Vektoren
       ! zu erreichen, setze nur einen Koeffizienten auf eins und
       ! alle anderen frei waehlbaren auf null:
       counter_test=0
       matrix_test: do
          counter_test=counter_test+1
          if (start/=0 .and. counter_test>start+1) then
             write(*,*)" ortho2: No satisfactory solution to the LES could be found"
             write(*,*)"         Try to re-formulate your constraints OR change the"
             write(*,*)"         variable MAT_TEST in routine ORTHO2 to a higher   "
             write(*,*)"         value"
             call error_handler(" ")
          endif
          coeff(:,i) = zero
          ! This is how it worked allus fine         
          if (start-counter_test+1>0) then
             coeff(ind(start-counter_test+1),i) = one
          endif
          if (start==0) coeff(ind(n_int),i)=one
          ! Zunaechst (wird spaeter weggelassen) wird die gesamt Matrix
          ! allociert und besetzt (matrix). Daraus wird die Matrix 'mat'
          ! isoliert und die rechte Seite ausgerechnet:'rhs'
          matrix=zero
          ! erst die Zeilen ,die der Bedingung v*c_cons=0 entsprechen
          ! Was kann hier Schlimmes passieren? Es koennte constraint vektoren
          ! geben, die vollstaendig ausserhalb des active set liegen!
          ! Dann waere eine Zeile von 'matrix' = 0
          do j=1,n_int
             do k=1,n_cons
                matrix(k,j) = dot_product(vec_in(:,ind(j)),c_cons(:,ind_cons(k)))
             enddo
          enddo
          if (i==1) then
             do k=1,n_cons
                if (abs(maxval(matrix(k,:)))<small .and. &
                     abs(minval(matrix(k,:)))<small ) then
                   call error_handler(" ortho_constraint: One of the constraint vectors&
                        & lies entirely out of the active set")
                endif
             enddo
          endif
          ! jetzt die Zeilen ,die der Bedingung entsprechen, das der aktuelle
          ! Vektor auf allen zuvor berechneten orthogonal ist
          do j=1,n_int
             counter_k=0
             do k=n_cons+1,dimi
                counter_k=counter_k+1
                matrix(k,j) = coeff(ind(j),counter_k)
             enddo
          enddo
          
          ! Jetzt berechne die rechte Seite:
          rhs=zero
          do k=1,dimi
             do j=1,start
                rhs(k) = rhs(k) + &
                     matrix(k,j)*coeff(ind(j),i)
             enddo
             rhs(k)=-rhs(k)
          enddo
          ! und jetzt setze 'mat' auf:
          mat=zero
          do k=1,dimi
             counter_j=0
             do j=start+1,n_int
                counter_j=counter_j+1
                if (counter_j>dimi) call error_handler&
                     ("ortho2: sth. fishy (1). Seek technical assistance")
                mat(k,counter_j) = matrix(k,j)
             enddo
          enddo
          
          m=int(dimi)
          n=m
          lda=m
          ldb=m
          info=0
          lwork=2*m
          allocate(work(lwork),STAT=alloc_stat)
          if(alloc_stat/=0) call error_handler(" ortho2: allocation (4) failed")
          work=zero
          nrhs=1
          call dgels('N',m,n,nrhs,mat,lda,rhs,ldb,work,lwork,info)
          if (info<0) then
             write(*,*)" ortho_constraint: on entry to the routine DGELS"
             write(*,*)"                   the ",-info," th argument    "
             write(*,*)"                   had an illegal value         "
             call error_handler(" ")
          endif
          deallocate(work,STAT=alloc_stat)
          if (alloc_stat/=0) call error_handler(" ortho2: deallocation (1) failed")
          
          !Wenn die Loesung zu grosse Komponenten hat, d.h. , wenn der
          !anfaenglich gewaehlte Vektor (durch Wahl von Coeff(1:start,i)
          !bestimmt) zu sehr *verdreht* werden musste, waehle
          !die STartkoeffizienten anders und mache den Mist noch mal
          if (abs(maxval(rhs))<mat_test) then
             exit matrix_test
          endif
          write(*,*)"ortho2: The solution to the LES just found is dropped due to "
          write(*,*)"        too high values in the solution vector. Try again with"
          write(*,*)"        a different initial choice of the arbitrary coefficients"

       enddo matrix_test
       
       ! Die Loesung des GLS steckt jetzt in rhs - umschaufeln
       ! nach coeff(start+1:n_int,i)
       counter_j=start
       do j=1,dimi
          counter_j=counter_j+1
          coeff(ind(counter_j),i)=rhs(j)
       enddo
       
       deallocate(matrix,mat,rhs,STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler(" ortho2: deallocation (2) failed")

       ! normiere den so erhaltenen Vektor
       help_vec = zero
       do j=1,n_int
          help_vec=help_vec+coeff(ind(j),i)*vec_in(:,ind(j))
       enddo
       abs_val = abs_value(help_vec)
       if (abs_val<=small) then
          write(*,*)"ortho2: Sth. fishy (2). The computed vector has length zero"
          call error_handler(" ")
       endif
       coeff(:,i)=coeff(:,i)/abs_val

    enddo vector_v
    vec_out = zero
    do k=1,n_cons
       vec_out(:,k) = c_cons(:,k)
    enddo
    counter_k=0
    do k=n_cons+1,n_int
       counter_k=counter_k+1
       do j=1,n_int
          vec_out(:,k) = vec_out(:,k) + coeff(j,counter_k)*vec_in(:,j)
       enddo
       abs_val = abs_value(vec_out(:,k))
       if (abs(abs_val-one)>small) call error_handler &
            (" ortho2: Sth. fishy (3). The computed vector does not have length one")
       do kk=1,k-1
          abs_val=dot_product(vec_out(:,k),vec_out(:,kk))
          if (abs(abs_val)>1.0e-8) then
             write(*,*)" Ortho2: Two vectors are found to be non-orthogonal"
             write(*,*)"         The calculated dotproduct is :",abs_val
             call error_handler(" ")
          endif
       enddo
    enddo

    ! In order to ensure the the projection of the new active set onto
    ! the original active set does not introduce new linear dependencies
    ! calculate the eigenvalues of the new set
    allocate(eigval(n_internal),eigvec(n_internal,n_internal),&
         help_mat(n_internal,n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler(" ortho2: allocation (5) failed")
    eigval=zero
    eigvec=zero
    help_mat=zero
    help_mat=matmul(transpose(vec_out),vec_in)
    call eigensolver(help_mat,n_internal,eigval,eigvec)
    minni=minval(abs(eigval))
    write(*,*)"ortho2: the smallest eigenvalue the projection of the new active"
    write(*,*)"        set onto the old active set ",minni
    write(*,*)"        If this values is too close to zero, this means that your"
    write(*,*)"        constraints contain linear dependencies "
    if (minni < 1.0e-8_r8_kind ) then
       call error_handler("ortho2: Linear dependencies detected in the constraints")
    endif
    
    deallocate (eigval,eigvec,help_mat,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler(" ortho2: deallocation (3) failed")
    
    
  end subroutine ortho2
    !** End of interface *****************************************
  !*************************************************************

  subroutine print_matrix(matrix,n,m,n_col)
    ! Purpose: print out a matrix if dimension n x m in a pretty
    !          format with n_col columns.
    !          n : numbers rows
    !          m : number of columns
    ! ----------------------------------------------------------
    real(kind=r8_kind),intent(in)      :: matrix(:,:)
    integer(kind=i4_kind),intent(in)   :: n,m,n_col
    !** End of interface *****************************************
    ! --- declaration of local variables -----------------------
    integer(kind=i4_kind) :: n_blocks,n_rest,counter,i_start,i_end,&
         j,i,k
    
    n_blocks=ceiling(real(m,kind=r8_kind)/real(n_col,kind=r8_kind))
    n_rest = mod(m,n_col)
    counter=1_i4_kind
    do j=1,n_blocks
       write(*,1202)(i,i=counter,counter+n_col-1)
       counter=counter+n_col
       i_start=(j-1)*n_col+1
       i_end=min(i_start+n_col-1,(n_blocks-1)*n_col+n_rest)
       do i=1,n
          write(*,1201)i,(matrix(i,k),k=i_start,i_end)
       enddo
       write(*,*)
    enddo

1201 format(i4,20(2x,f10.7))   
1202 format(4x,20(8x,i3)) 
    
  end subroutine print_matrix
!!$  !*************************************************************
!!$  
!!$  !subroutine eigensolver(matrix,dimen,eigval,eigvec,destroy)
!!$    ! Purpose: wrapper for 'evvrsp'-Routine from EISPACK.
!!$    !          This routine is used for compatibilty with
!!$    !          A.Voityuk. Do not change it without having a good 
!!$    !          reason for doing so.
!!$    ! 
!!$    ! ----------------------------------------------------------
!!$    !integer(kind=i4_kind),intent(in)   :: dimen      
!!$    !real(kind=r8_kind),intent(inout)   :: matrix(dimen,dimen)
!!$    !real(kind=r8_kind),intent(inout)   :: eigval(dimen)
!!$    !real(kind=r8_kind),intent(inout)   :: eigvec(dimen,dimen)
!!$    !logical,optional                   :: destroy
!!$    !** End of interface *****************************************
!!$    ! -----------------------------------------------------------
!!$    !real(kind=r8_kind)             :: b(dimen,9)
!!$    !real(kind=r8_kind),allocatable :: help_mat(:,:)
!!$    integer(kind=i4_kind)          :: iwork(dimen)
!!$    integer(kind=i4_kind)          :: ierr,alloc_stat
!!$    logical                        :: local
!!$    
!!$    if (present(destroy)) then
!!$       if (destroy) then
!!$          local = .true.
!!$       else
!!$          local= .false.
!!$       endif
!!$    else
!!$       local = .false.
!!$    endif
!!$
!!$    if (.not.local) then
!!$       allocate(help_mat(dimen,dimen),STAT=alloc_stat)
!!$       if (alloc_stat/=0) then
!!$          write(*,*)"eigensolver: allocation (1) failed"
!!$          stop 1
!!$       endif
!!$       help_mat = matrix
!!$    endif
!!$
!!$    eigval = zero
!!$    eigvec=zero
!!$    iwork = 0_i4_kind
!!$    b = zero
!!$    call evvrsp(dimen,dimen,dimen,matrix,b,iwork,eigval,eigvec,ierr)
!!$
!!$    if (.not.local) then
!!$       matrix=help_mat
!!$       deallocate(help_mat,STAT=alloc_stat)
!!$       if (alloc_stat/=0) then
!!$          write(*,*)" eigensolver: deallocation (1) failed"
!!$          stop 1
!!$       endif
!!$    endif
!!$    if(ierr/=0) then
!!$       if (ierr<0) then
!!$          write(*,*)"eigensolver: Iteration for eigevector ",ierr," failed"
!!$          stop
!!$       else
!!$          write(*,*)"eigensolver: Iteration for eigenvalue ",ierr," failed"
!!$       endif
!!$    endif
!!$
!!$  end subroutine eigensolver

  !*************************************************************

  !--------------- End of module ----------------------------------
subroutine error_handler(message)
  ! Purpose: substitute for the larger, pvm-infected subroutine
  !          'error_handler' used by ParaGAU.
  !------------------------------------------------------------
  use type_module
  implicit none
  !------------ Declaration of formal parameters --------------
  character(LEN=*) :: message
  ! --- executable code ---------------------------------------
  write(*,*)" error_handler: ",message
  stop 1
end subroutine error_handler
end module math_module
