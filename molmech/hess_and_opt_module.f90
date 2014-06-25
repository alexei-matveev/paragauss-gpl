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
module hess_and_opt_module
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use tasks_main_options_module
  use inp_out_module
  use mm_timer_module
  use slab_module
  use species_module
  use n_body_lists_module
  use energy_and_forces_module
  use calc_energy_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  integer(kind=i4_kind), public :: n_iterations
  !------------ public functions and subroutines ------------------
  public read_opt_options, read_hess_options, init_hess, calc_hess, minimize, &
       write_opt_options_to_output, write_hess_options_to_output, shutdown_hess, &
       alloc_H_on_slaves,shutdown_H_on_slaves
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  real(kind=r8_kind), allocatable :: Hess(:)
  logical :: parab_3

  integer(kind=i4_kind) :: hess_update_cycles
  character(len=20) :: method, ls_method
  real(kind=r8_kind) :: gnorm, gmax, rnorm, rmax
  real(kind=r8_kind) :: max_step

  integer(kind=i4_kind) :: df_n_iterations=1000
  integer(kind=i4_kind) :: df_hess_update_cycles=10
  character(len=20) :: df_method="NEWTON"
  character(len=20) :: df_ls_method="3_PARAB" ! "STEP"
  real(kind=r8_kind) :: df_gnorm=default_value
  real(kind=r8_kind) :: df_gmax=default_value1
  real(kind=r8_kind) :: df_rnorm=default_value
  real(kind=r8_kind) :: df_rmax=default_value1
  real(kind=r8_kind) :: df_max_step=one

  namelist /opt_options/ n_iterations, method, gnorm, gmax, rnorm, rmax, &
       hess_update_cycles, max_step, ls_method

  character(len=20) :: hess_type
  character(len=6) :: update_type
  logical :: save_hessian

  character(len=20) :: df_hess_type="UNIT" ! or "NUMERICAL" or "ANALYTICAL"
  character(len=6) :: df_update_type="BFGS"
  logical :: df_save_hessian=.false.
  real(kind=r8_kind) :: df_weight_hess=one

  namelist /hess_options/ hess_type, update_type, save_hessian
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  function read_opt_options()

    logical :: read_opt_options
    integer(i4_kind) :: ii

    n_iterations=df_n_iterations
    method=df_method
    ls_method=df_ls_method
    gnorm=df_gnorm
    gmax=df_gmax
    rnorm=df_rnorm
    rmax=df_rmax
    hess_update_cycles=df_hess_update_cycles
    max_step=df_max_step

    call go_to_first_input_line
    read_opt_options=find_namelist("&OPT_OPTIONS",ii)
    if(.not.read_opt_options) goto 100
    read(input_device,nml=opt_options, end=100, err=200)
    call upcase(method)
    if(trim(method) /= "NEWTON") call error_handler &
         ("MolMech: Optimization method has to be NEWTON")
    call upcase(ls_method)
    if(trim(ls_method) /= "3_PARAB" .and.  trim(ls_method) /= "STEP") call error_handler &
         ("MolMech: Definition of line search method is not correct")
    if(gmax < zero .or. rmax < zero) call error_handler &
         ("MolMech: GMAX and RMAX have to be positive")

    if(n_iterations < 1) n_iterations=1

    read_opt_options=.true.
    return

100 read_opt_options=.false.
    return

200 call input_nm_error(0,"OPT_OPTIONS")

  end function read_opt_options
  !****************************************************************

  !****************************************************************
  subroutine write_opt_options_to_output()

    write(output_device,'(a25,/)') " STRUCTURAL OPTIMIZATION:"

    if(trim(method)=="NEWTON")  write(output_device,'(a22)') " Metod: Newton-Raphson"

    write(output_device,'(a11)') " Tolerance:"
    write(output_device,'(a8,e9.3)') " rnorm: ", rnorm
    write(output_device,'(a8,e9.3)') " rmax:  ", rmax
    write(output_device,'(a8,e9.3)') " gnorm: ", gnorm
    write(output_device,'(a8,e9.3,/)') " gmax:  ", gmax
    write(output_device,'(a33,i5)') " Total number of geometry steps: ", n_iterations
    write(output_device,'(a32,f9.7)') " Maximum size of geometry step: ", max_step
    write(output_device,'(a24,i5,a6)') " Updating hessian every ", hess_update_cycles," steps"
    if(trim(ls_method) == "3_PARAB") then
       parab_3=.true.
       write(output_device,'(a35,i5,a6)') " Line search method: cubic parabola"
    else if(trim(ls_method) == "STEP") then
       parab_3=.false.
       write(output_device,'(a30,i5,a6)') " Line search method: line step"
    end if
!!$    if(update_nb_list >= n_iterations) then
!!$       write(output_device,'(a28,/)') " No updating nonbonding list"
!!$    else
!!$       write(output_device,'(a32,i5,a6)') " Updating nonbonding list every ",update_nb_list," steps"
!!$    end if

  end subroutine write_opt_options_to_output
  !****************************************************************

  !****************************************************************
  function read_hess_options()
 
    logical :: read_hess_options
    integer(i4_kind) :: ii

    hess_type=df_hess_type
    update_type=df_update_type
    save_hessian=df_save_hessian
    
    call go_to_first_input_line
    read_hess_options=find_namelist("&HESS_OPTIONS",ii)
    if(.not.read_hess_options) goto 100
    read(input_device,nml=hess_options, end=100, err=200)
    call upcase(hess_type)
    call upcase(update_type)

    if(solvent) hess_type=df_hess_type
    
    if(trim(hess_type) == "NUMERICAL" .and. calc_strain) then
       call error_handler("MolMech: Nunerical Hessian and strain calculation - impossible")
    end if

    if(trim(hess_type) == 'UNIT') hess_update_cycles=n_iterations

    read_hess_options=.true.
    return

100 if(trim(hess_type) == 'UNIT') hess_update_cycles=n_iterations

    read_hess_options=.false.
    return

200 call input_nm_error(0,"HESS_OPTIONS")

  end function read_hess_options
  !****************************************************************

  !****************************************************************
  subroutine write_hess_options_to_output()

    if(calc_optimization) then
       if(trim(hess_type) == "UNIT") & 
            write(output_device,'(a29)') " Initial hessian: unit matrix"
       if(trim(hess_type) == "NUMERICAL") & 
            write(output_device,'(a27)') " Initial hessian: numerical"
       if(trim(hess_type) == "ANALYTICAL") & 
            write(output_device,'(a28)') " Initial hessian: analytical"
       if(trim(update_type) == "BFGS") &
            write(output_device,'(a27,/)') " BFGS hessian update scheme"
    end if

  end subroutine write_hess_options_to_output
  !****************************************************************

  !****************************************************************
  subroutine init_hess()

    integer(kind=i4_kind) :: hess_length
    integer(kind=i4_kind) :: status,ns

    ns=0
    if(calc_strain) then
       if(lattice_calc) ns=6
       if(slab_calc) ns=3
    end if

    hess_length=(3*n_species)*(3*n_species+1)/2
    if(lattice_calc .or. slab_calc) hess_length= &
         (3*(n_species-1)+ns)*(3*(n_species-1)+ns+1)/2

    allocate(Hess(hess_length),stat=status)
    if(status /= 0) call error_handler( &
         "MolMech: failed HESS allocation")
    Hess=zero

  end subroutine init_hess
  !****************************************************************

  !****************************************************************
  subroutine alloc_H_on_slaves()

    integer(i4_kind) :: status,ns

    ns=0
    if(calc_strain) then
       if(lattice_calc) ns=6
       if(slab_calc) ns=3
    end if

    allocate(h(3*n_species+ns,3*n_species+ns),stat=status)
    if(status /= 0) call error_handler("MolMech: failed H allocation on slaves")
    h=zero
       
  end subroutine alloc_H_on_slaves
  !****************************************************************

  !****************************************************************
  subroutine calc_hess()

    use f77_lapack, only : dsptrf,dsptri,dsyev

    integer(i4_kind) :: i,j,k,jh,status,info,ns,length,fi
    real(r8_kind), allocatable :: he(:,:),r0(:),work(:),eval(:)
    real(r8_kind), allocatable :: g0(:),g1(:),g2(:)
    integer(i4_kind), allocatable :: ipiv(:)
    logical :: store
    integer(i4_kind) :: file_hess,h_size

    write(output_device,'(a19)') " Hessian calculated"
    write(6,'(a19)') " Hessian calculated"

    Hess=zero

    ns=0
    if(calc_strain) then
       if(lattice_calc) ns=6
       if(slab_calc) ns=3
    end if

    if(calc_optimization .and. trim(hess_type) == 'UNIT') then
       fi=0 
       if(lattice_calc) fi=1
       k=0
       do i=1,3*(n_species-fi)+ns
          k=k+i
          Hess(k)=0.001_r8_kind
       end do
    else

       allocate(h(3*n_species+ns,3*n_species+ns),stat=status)
       if(status /= 0) call error_handler("MolMech: failed H allocation")
       h=zero

       if(trim(hess_type) == 'NUMERICAL') then
          allocate(r0(3*n_species),g0(3*n_species),g1(3*n_species), &
               g2(3*n_species),stat=status)
          if(status /= 0) call error_handler("MolMech: failed G1 or G2 allocation")

          store=calc_hessian
          calc_hessian=.false.

          !saving current data
          do i=1,n_species
             do j=1,3
                k=3*(i-1)+j
                r0(k)=atoms_cart(i)%r(j)
                g0(k)=Grad(j,i)
             end do
          end do

          !Hessian calculation
          do k=1,3*n_species
             i=int((k-1)/3)+1
             j=k-3*(i-1)
          
             atoms_cart(i)%r(j)=r0(k)+0.0001_r8_kind
             if(lattice_calc) call cart2frac(n_species,vect,atoms_cart,atoms_frac)

             call total_energy_and_grad(do_corr=.true.)
             g1=reshape(Grad,(/3*n_species/))

             atoms_cart(i)%r(j)=r0(k)-0.0001_r8_kind
             if(lattice_calc) call cart2frac(n_species,vect,atoms_cart,atoms_frac)

             call total_energy_and_grad(do_corr=.true.)
             g2=reshape(Grad,(/3*n_species/))

             atoms_cart(i)%r(j)=r0(k)
             if(lattice_calc) call cart2frac(n_species,vect,atoms_cart,atoms_frac)

             do jh=k,3*n_species
                h(k,jh)=(g1(jh)-g2(jh))/0.0002_r8_kind
                if(k /= jh) h(jh,k)=h(k,jh)
             end do
          end do

!!$          do i=1,3*n_species
!!$             do j=i+1,3*n_species
!!$                h(i,j)=(h(i,j)+h(j,i))/two
!!$                h(j,i)=h(i,j)
!!$             end do
!!$          end do

!!$          do i=1,3*n_species
!!$             do j=1,3*n_species
!!$print*,i,j,H(i,j) !*kjm2ev
!!$             end do
!!$          end do

          if(save_hessian .and. .not. calc_optimization) then
             call get_file_device(file_hess,'molmech.hess','out')

             do i=1,3*n_species+ns
                do j=1,3*n_species+ns
                   if(i <= 3*n_species .and. j <= 3*n_species) &
                        write(file_hess,*)' r',i,'r',j,H(i,j)*kjm2ev
                   if(i <= 3*n_species .and. j >  3*n_species) &
                        write(file_hess,*)' r',i,'s',j-3*n_species,H(i,j)*kjm2ev
                   if(i >  3*n_species .and. j <= 3*n_species) &
                        write(file_hess,*)' s',i-3*n_species,'r',j,H(i,j)*kjm2ev
                   if(i >  3*n_species .and. j >  3*n_species) &
                        write(file_hess,*)' s',i-3*n_species,'s',j-3*n_species,H(i,j)*kjm2ev
                end do
             end do
             
             call close_file_device(file_hess)
          end if

          !restoring current data
          do i=1,n_species
             do j=1,3
                k=3*(i-1)+j
                atoms_cart(i)%r(j)=r0(k)
                Grad(j,i)=g0(k)
             end do
          end do

          calc_hessian=store

          deallocate(r0,g0,g1,g2,stat=status)
          if(status /= 0) call error_handler("MolMech: failed G1 or G2 deallocation")

          if(to_pg_hess .and. .not. calc_optimization) then
             call del_shell_contrib(h_size)
             call get_file_device(file_hess,'hesse_cartesian.dat','inp')
             write(file_hess,*) H(1:h_size,1:h_size)/(h2kJm*a2b*a2b)
             call close_file_device(file_hess)
          end if

          if(.not. calc_optimization) return
       end if

       if(trim(hess_type) == 'ANALYTICAL') then
          store=calc_hessian
          calc_hessian=.true.
          call total_energy_and_grad(do_corr=.true.)
          calc_hessian=store

!          allocate(work(3*(3*n_species+ns)-1),eval(3*n_species+ns),stat=status)
!          if(status /= 0) call error_handler("MolMech: failed WORK allocation(0)")
!
!          length=3*n_species+ns
!          call dsyev('V','L',length,h,3*n_species+ns,eval,work,3*length-1,info)
!          if(info /= 0) call error_handler("MolMech: failed calculation of eigenvalues(0)")
!
!          deallocate(work,stat=status)
!          if(status /= 0) call error_handler("MolMech: failed WORK deallocation(0)")
!
!          do i=1,3*n_species+ns
!             if(eval(i) < zero) then
!                write(output_device,'(a22,i6,e13.5)') " Negative eigenvalue: ",i,eval(i)
!                write(6,'(a22,i6,e13.5)') " Negative eigenvalue: ",i,eval(i)
!             end if
!          end do
!          deallocate(eval,stat=status)
!          if(status /= 0) call error_handler("MolMech: failed EVAL deallocation(0)")

          if(save_hessian .and. .not. calc_optimization) then
             call get_file_device(file_hess,'molmech.hess','out')

             do i=1,3*n_species+ns
                do j=1,3*n_species+ns
                   if(i <= 3*n_species .and. j <= 3*n_species) &
                        write(file_hess,*)' r',i,'r',j,H(i,j)*kjm2ev
                   if(i <= 3*n_species .and. j >  3*n_species) &
                        write(file_hess,*)' r',i,'s',j-3*n_species,H(i,j)*kjm2ev
                   if(i >  3*n_species .and. j <= 3*n_species) &
                        write(file_hess,*)' s',i-3*n_species,'r',j,H(i,j)*kjm2ev
                   if(i >  3*n_species .and. j >  3*n_species) &
                        write(file_hess,*)' s',i-3*n_species,'s',j-3*n_species,H(i,j)*kjm2ev
                end do
             end do
             
             call close_file_device(file_hess)
          end if


          if(to_pg_hess .and. .not. calc_optimization) then
             call del_shell_contrib(h_size)
             call get_file_device(file_hess,'hesse_cartesian.dat','inp')
             write(file_hess,*) H(1:h_size,1:h_size)/(h2kJm*a2b*a2b)
             call close_file_device(file_hess)
          end if

          if(.not. calc_optimization) return
       
          if(lattice_calc .or. slab_calc) call hess2frac()
          
!!$          do i=1,3*n_species+ns
!!$             do j=1,3*n_species+ns
!!$print*,i,j,H(i,j)*kjm2ev
!!$             end do
!!$          end do

       end if

       call start_mm_timer(hinv_time)
       if(lattice_calc.or.slab_calc) then
          !pack 2-d hessian into 1-d array
          k=0
          do i=4,3*n_species+ns
             do j=4,i
                k=k+1
                hess(k)=h(j,i)
!!$print*,k,j-3,i-3,hess(k)*kjm2ev
             end do
          end do

          deallocate(h,stat=status)
          if(status /= 0) call error_handler("MolMech: failed H deallocation(1)")
!!$call cpu_time(tt)
!!$print*,tt
          length=3*(n_species-1)+ns
          allocate(ipiv(length),stat=status)
          if(status /= 0) call error_handler("MolMech: failed ipiv allocation(1)")
          call dsptrf('U',length,hess,ipiv,info)
       
          if(info /= 0) then 
!!$             print*,info
             call error_handler("MolMech: dsptrf: info /= 0")
          else
             allocate(work(length),stat=status)
             if(status /= 0) call error_handler("MolMech: failed work allocation(1)")
             call dsptri('U',length,hess,ipiv,work,info)
             if(info /= 0) call error_handler("MolMech: dsptri: info /= 0")
             deallocate(ipiv,work,stat=status)
             if(status /= 0) call error_handler("MolMech: failed ipiv deallocation(1)")
          end if
!!$call cpu_time(tt)
!!$print*,tt
       else
          ! generalized inversion of hessian
          allocate(work(3*(3*n_species+ns)-1),eval(3*n_species+ns),stat=status)
          if(status /= 0) call error_handler("MolMech: failed WORK allocation")

          length=3*n_species+ns
          call dsyev('V','L',length,h,3*n_species+ns,eval,work,3*length-1,info)
          if(info /= 0) call error_handler("MolMech: failed calculation of eigenvalues(1)")
!!$call cpu_time(tt)
!!$print*,tt

          deallocate(work,stat=status)
          if(status /= 0) call error_handler("MolMech: failed WORK deallocation")

          do i=1,3*n_species+ns
             if(abs(eval(i)) <= small*hundred) then
                eval(i)=zero
             else
                eval(i)=one/eval(i)
             end if
          end do
!!$call cpu_time(tt)
!!$print*,tt

          allocate(he(3*n_species+ns,3*n_species+ns),stat=status)
          if(status /= 0) call error_handler("MolMech: failed HE allocation")

          do i=1,3*n_species+ns
             do j=1,3*n_species+ns
                he(i,j)=zero
                do k=1,3*n_species+ns
                   he(i,j)=he(i,j)+h(i,k)*h(j,k)*eval(k)
                end do
             end do
          end do
!!$call cpu_time(tt)
!!$print*,tt

          deallocate(eval,stat=status)
          if(status /= 0) call error_handler("MolMech: failed EVAL deallocation")
          deallocate(h,stat=status)
          if(status /= 0) call error_handler("MolMech: failed H deallocation")

          k=0
          do i=1,3*n_species+ns
             do j=1,i
                k=k+1
                hess(k)=he(j,i)
!!$print*,k,hess(k) 
             end do
          end do
!!$call cpu_time(tt)
!!$print*,tt

          deallocate(he,stat=status)
          if(status /= 0) call error_handler("MolMech: failed HE deallocation")
       end if
       call stop_mm_timer(hinv_time)
    end if

  contains
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    subroutine hess2frac()

      real(r8_kind) :: transmat(3,3),transmat_s(2,2)
      real(r8_kind) :: tmp1(3,3),tmp2(6,3),tmp(6,3)
      real(r8_kind) :: tmp1_s(2,2),tmp2_s(3,2),tmp_s(3,2)
      integer(i4_kind) :: i1,j1,ki,kj

      if(lattice_calc) then
         transmat(:,1)=vect%v1
         transmat(:,2)=vect%v2
         transmat(:,3)=vect%v3

         do i1=1,n_species
            ki=3*(i1-1)
            !converting coordinate derivatives
            do j1=i1,n_species
               kj=3*(j1-1)
               tmp1=matmul(H(ki+1:ki+3,kj+1:kj+3),transmat)
               H(ki+1:ki+3,kj+1:kj+3)=matmul(transpose(transmat),tmp1)
               H(kj+1:kj+3,ki+1:ki+3)=transpose(H(ki+1:ki+3,kj+1:kj+3))
            end do
         end do

         if(calc_strain) then
            do i1=1,n_species
               ki=3*(i1-1)
               !converting mixed derivatives
               tmp=transpose(H(ki+1:ki+3,3*n_species+1:3*n_species+6))
               tmp2=matmul(tmp,transmat)
               H(ki+1:ki+3,3*n_species+1:3*n_species+6)=transpose(tmp2)
               H(3*n_species+1:3*n_species+6,ki+1:ki+3)=tmp2
            end do
         end if
      else if(slab_calc) then
         transmat_s(:,1)=vect_s%v1
         transmat_s(:,2)=vect_s%v2

         do i1=1,n_species
            ki=3*(i1-1)
            !converting coordinate derivatives
            do j1=i1,n_species
               kj=3*(j1-1)
               tmp1_s=matmul(H(ki+1:ki+2,kj+1:kj+2),transmat_s)
               H(ki+1:ki+2,kj+1:kj+2)=matmul(transpose(transmat_s),tmp1_s)
               H(kj+1:kj+2,ki+1:ki+2)=transpose(H(ki+1:ki+2,kj+1:kj+2))
            end do
         end do

         if(calc_strain) then
            do i1=1,n_species
               ki=3*(i1-1)
               !converting mixed derivatives
               tmp_s=transpose(H(ki+1:ki+2,3*n_species+1:3*n_species+3))
               tmp2_s=matmul(tmp_s,transmat_s)
               H(ki+1:ki+2,3*n_species+1:3*n_species+3)=transpose(tmp2_s)
               H(3*n_species+1:3*n_species+3,ki+1:ki+2)=tmp2_s
            end do
         end if
      end if

    end subroutine hess2frac

  end subroutine calc_hess
  !****************************************************************

  !****************************************************************
  subroutine del_shell_contrib(h_size)
    !This routine is used to eliminate shell contributions to the 
    !total Hessian of system to save it in format appropriate for
    !reading in by Optimizer
    !Be careful, after correction the Hessian cannot be used to continue
    !calculations
    integer(i4_kind) :: h_size
    integer(i4_kind) :: type,n_sh,i_buf,n_cr
    integer(i4_kind) :: i,j,k,l,ll,ll1,kk,ns,status
    real(r8_kind), allocatable :: buffer(:,:),buffer1(:,:),H_ss(:,:)
    integer(i4_kind), allocatable :: i_buffer(:)
    integer(i4_kind), allocatable :: ipiv(:)
    real(r8_kind), allocatable :: work(:)
    integer(i4_kind) :: info
    real(r8_kind) :: aaaa

    n_sh=0
    do i=1,n_species
       type=atoms_cart(i)%type
       if(get_c_s(atoms(type)%c_s) == "S") n_sh=n_sh+1
    end do
    h_size=3*n_species !-3*n_sh
return
    n_cr=n_species-n_sh
    if(n_sh == 0) return

    allocate(buffer(3,3*n_species),buffer1(3*n_species,3), &
         i_buffer(n_species),stat=status)
    if(status /= 0) call error_handler &
         ("MolMech: Hess_and_Opt: failed BUFFER allocation")

    !Rearraging Hessian to store all shell elements on right and
    !down sides
    do i=1,n_species
       i_buffer(:)=i
    end do

    ns=0
    i=1
    do
       k=i_buffer(i)
       type=atoms_cart(k)%type
       if(get_c_s(atoms(type)%c_s) == "S") then
          ns=ns+1
          i_buf=k
          kk=3*(k-1)
          buffer=H(kk+1:kk+3,:)
          do l=k,n_species-1
             ll=3*(l-1); ll1=3*l
             i_buffer(l)=i_buffer(l+1)             
             H(ll+1:ll+3,:)=H(ll1+1:ll1+3,:)
          end do
          H(3*(n_species-1):3*n_species,:)=buffer

          buffer1=H(:,kk+1:kk+3)
          do l=k,n_species-1
             ll=3*(l-1); ll1=3*l
             H(:,ll+1:ll+3)=H(:,ll1+1:ll1+3)
          end do
          H(3*(n_species-1):3*n_species,:)=buffer
          i_buffer(n_species)=i_buf
          if(ns == n_sh) exit
       else
          i=i+1
       end if
    end do

    deallocate(buffer,buffer1,i_buffer,stat=status)
    if(status /= 0) call error_handler &
         ("MolMech: Hess_and_Opt: failed BUFFER deallocation")

    !Inverting shell-shell Hessian 
    allocate(H_ss(3*n_sh,3*n_sh),ipiv(3*n_sh),work(3*n_sh),stat=status)
    if(status /= 0) call error_handler &
         ("MolMech: Hess_and_Opt: failed H_ss allocation")
    H_ss=H(3*n_cr+1:3*n_species,3*n_cr+1:3*n_species)

    call dsytrf('U',3*n_sh,H_ss,3*n_sh,ipiv,work,3*n_sh,info)
    if(info == 0) then
       call dsytri('U',3*n_sh,H_ss,3*n_sh,ipiv,work,info)
       if(info /= 0) call error_handler("MolMech: Hess_and_Opt: dsytri: info /= 0")
    else
       call error_handler("MolMech: Hess_and_Opt: dsytrf: info /= 0")
    endif

    deallocate(ipiv,work,stat=status)
    if(status /= 0) call error_handler &
         ("MolMech: Hess_and_Opt: failed IPIV allocation")

    H(3*n_cr+1:3*n_species,3*n_cr+1:3*n_species)=H_ss

    deallocate(H_ss,stat=status)
    if(status /= 0) call error_handler &
         ("MolMech: Hess_and_Opt: failed H_ss deallocation")

    !correct Hessian h=r-t*s-1*t
    do i=1,3*n_cr
       do j=1,3*n_sh
          aaaa=zero
          do l=1,3*n_sh
             aaaa=aaaa+H(3*n_cr+j,3*n_cr+l)*H(i,3*n_cr+l)
          end do
          H(3*n_cr+j,i)=aaaa
       end do
    end do

    do i=1,3*n_cr
       do j=1,3*n_cr
          aaaa=zero
          do l=1,3*n_sh
             aaaa=aaaa-H(j,3*n_cr+l)*H(3*n_cr+l,i)
          end do
          H(i,j)=H(i,j)+aaaa
       end do
    end do

  end subroutine del_shell_contrib
  !****************************************************************

  !****************************************************************
  subroutine minimize()

    real(kind=r8_kind), allocatable :: rstep(:),r_cur(:),r_old(:)
    real(kind=r8_kind), allocatable :: gstep(:),g_cur(:),g_old(:)
    real(kind=r8_kind) :: e_cur, e_old
    real(kind=r8_kind) :: rs_limit,str(6),str_s(3)
    real(kind=r8_kind) :: gnor, gabs, rnor, rabs, cosgr
    logical :: calc,store,corr_sys,ok
    integer(kind=i4_kind) :: i,j,k,nb,i1,j1,status,ns,nvar,fi

    fi=1
    if(lattice_calc.or.slab_calc) fi=2

    ns=0; nvar=3*(n_species)
    if(calc_strain) then 
       if(lattice_calc) ns=6
       if(slab_calc) ns=3
    end if
    if(lattice_calc.or.slab_calc) nvar=3*(n_species-1)

    gnor=zero; gabs=zero; rnor=zero; rabs=zero
    if(lattice_calc) then
       rs_limit=one/volume**(one/three)
    else if(slab_calc) then
       rs_limit=one/sqrt(area)
    else
       rs_limit=one
    end if

    if(trim(method) == "NEWTON") then

       allocate(rstep(nvar+ns),r_cur(nvar+ns),r_old(nvar+ns), &
            g_cur(nvar+ns),g_old(nvar+ns), stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed RSTEP or GSTEP allocation")

       call start_mm_timer(optim_time)
       calc=.true.
       j=0; nb=0
       call output_minim_run(0)

       n_iter: do i=1,n_iterations

          if(calc) then
             store=calc_hessian
             calc_hessian=.false.
             if(i == 0 .and. trim(hess_type) /= 'NUMERICAL') then
                corr_sys=.false.
             else
                corr_sys=.true.
             end if
             call total_energy_and_grad(do_corr=corr_sys)
             calc_hessian=store

             e_cur=E_total
             if(lattice_calc.or.slab_calc) then
                call grad2frac()
                k=0
                do i1=fi,n_species
                   do j1=1,3
                      k=k+1
                      r_cur(k)=atoms_frac(i1)%r(j1)
                      g_cur(k)=Grad(j1,i1)
                   end do
                end do
                if(calc_strain) then
                   do i1=1,ns
                      if(lattice_calc) r_cur(nvar+i1)=strain(i1)
                      if(slab_calc) r_cur(nvar+i1)=strain_s(i1)
                      g_cur(nvar+i1)=Grad_s(i1)
                   end do
                end if
             else
                do i1=fi,n_species
                   do j1=1,3
                      k=3*(i1-1)+j1
                      r_cur(k)=atoms_cart(i1)%r(j1)
                      g_cur(k)=Grad(j1,i1)
                   end do
                end do
             end if
          end if

          if(i==1 .or. j==hess_update_cycles) then
             call calc_hess()
             j=0
          else
             call update_hess()
          end if

          call stop_mm_timer(optim_time)
          call output_minim_run(i)
          if(i==n_iterations) exit n_iter
!!$          if(converged()) exit n_iter
          call start_mm_timer(optim_time)

          call next_step()
          if(converged()) exit n_iter

          r_old=r_cur
          r_cur=r_cur+rstep
          if(lattice_calc.or.slab_calc) then
             do i1=1,nvar
                if(slab_calc .and. mod(i1,3) == 0) cycle
                if(r_cur(i1) < zero) r_cur(i1)=r_cur(i1)+one
                if(r_cur(i1) >= one) r_cur(i1)=r_cur(i1)-one
             end do
          end if
          g_old=g_cur
          e_old=e_cur

          if(lattice_calc.or.slab_calc) then
             k=0
             do i1=fi,n_species
                do j1=1,3
                   k=k+1
                   atoms_frac(i1)%r(j1)=r_cur(k)
                end do
             end do
             if(calc_strain) then
                if(lattice_calc) then
                   strain=r_cur(nvar+1:nvar+6)
                   str=rstep(nvar+1:nvar+6)
                   call store_vect()
                   call renew_lattice(str)
                else if(slab_calc) then
                   strain_s=r_cur(nvar+1:nvar+3)
                   str_s=rstep(nvar+1:nvar+3)
                   call store_vect_s()
                   call renew_slab(str_s)
                end if
             end if
          else
             do i1=fi,n_species
                do j1=1,3
                   k=3*(i1-1)+j1
                   atoms_cart(i1)%r(j1)=r_cur(k)
                end do
             end do
          end if

          if(parab_3) then
             call line_search()
             ok=.true.
          else
             call line_search_per(ok)
          end if
          
          j=j+1
          if(.not.ok) j=hess_update_cycles

          if(lattice_calc) then
             if(cosgr <= 0.01_r8_kind) j=hess_update_cycles
          else
             if(cosgr <= 0.01_r8_kind) j=hess_update_cycles
          end if

       end do n_iter
       call stop_mm_timer(optim_time)

       deallocate(rstep,r_cur,r_old,g_cur,g_old, stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed RSTEP or GSTEP deallocation")

       store=calc_hessian
       calc_hessian=.false.
       call total_energy_and_grad(do_corr=.true.)
       calc_hessian=store

    end if
    
  contains
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    subroutine update_hess()

      integer(kind=i4_kind) :: iu,ju,ku
      real(kind=r8_kind) :: rs_gs,h_gs_gs
      real(kind=r8_kind), allocatable :: h_gs(:)

      allocate(h_gs(nvar+ns),gstep(nvar+ns),stat=status)
      if(status /= 0) call error_handler( &
           "MolMech: failed H_GS allocation")

      gstep=g_cur-g_old

      rs_gs=dot_product(rstep,gstep)
      if(rs_gs < small) rs_gs=small
      
!!$      ku=0
!!$      do iu=1,nvar+ns
!!$         h_gs(iu)=zero
!!$         do ju=1,iu
!!$            h_gs(iu)=h_gs(iu)+gstep(ju)*hess(ku+ju)
!!$         enddo
!!$         ku=ku+iu
!!$      enddo
!!$
!!$      ku=1
!!$      do iu=2,nvar+ns
!!$         do ju=1,iu-1
!!$            h_gs(ju)=h_gs(ju)+hess(ku+ju)*gstep(iu)
!!$         enddo
!!$         ku=ku+iu
!!$      enddo
      call dspmv('U',nvar+ns,one,hess,gstep,1,zero,h_gs,1)
 
      h_gs_gs=dot_product(h_gs,gstep)
      
      if(trim(update_type) == "BFGS") then
         ku=0
         do iu=1,nvar+ns
            do ju=1,iu
               ku=ku+1
               hess(ku)=hess(ku)-h_gs(ju)*rstep(iu)/rs_gs- &
                    rstep(ju)*h_gs(iu)/rs_gs+rstep(ju)*rstep(iu)/rs_gs+ &
                    rstep(ju)*rstep(iu)*h_gs_gs/(rs_gs*rs_gs)
            end do
         end do
      end if

      deallocate(h_gs,gstep,stat=status)
      if(status /= 0) call error_handler( &
           "MolMech: failed H_GS deallocation")

    end subroutine update_hess
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    subroutine next_step()

      real(r8_kind) :: gc,rr

!!$      ku=0
!!$      do iu=1,nvar+ns
!!$         rstep(iu)=zero
!!$         do ju=1,iu
!!$            rstep(iu)=rstep(iu)+g_cur(ju)*hess(ku+ju)
!!$         enddo
!!$         ku=ku+iu
!!$      enddo
!!$
!!$      ku=1
!!$      do iu=2,nvar+ns
!!$         do ju=1,iu-1
!!$            rstep(ju)=rstep(ju)+hess(ku+ju)*g_cur(iu)
!!$         enddo
!!$         ku=ku+iu
!!$      enddo
      call dspmv('U',nvar+ns,one,hess,g_cur,1,zero,rstep,1)
      rstep=-rstep*max_step

      gc=dot_product(g_cur,rstep)
      rr=sqrt(dot_product(rstep,rstep))
print*,gc,rr,rs_limit*1.5_r8_kind
      if(rr > 1.5_r8_kind*rs_limit) then
         rstep=rstep*1.5_r8_kind*rs_limit/rr
         rr=sqrt(dot_product(rstep,rstep))
      end if
      if(gc > zero) rstep=-rstep

    end subroutine  next_step
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    subroutine line_search()

      integer(i4_kind) :: iu,ju,ku
      real(r8_kind) :: gc,go,rr,rs
      real(r8_kind) :: a3,a2,a1,a0,b1,b2

      store=calc_hessian
      calc_hessian=.false.
      call total_energy_and_grad(do_corr=.true.)
      calc_hessian=store

      if(lattice_calc.or.slab_calc) call grad2frac()
      ku=0
      do iu=fi,n_species
         do ju=1,3
            ku=ku+1
            g_cur(ku)=Grad(ju,iu)
         end do
      end do
      if(calc_strain) then
         do i1=1,ns
            g_cur(nvar+i1)=Grad_s(i1)
         end do
      end if
      e_cur=E_total

      rr=sqrt(dot_product(rstep,rstep))
      go=dot_product(g_old,rstep)/rr
      gc=dot_product(g_cur,rstep)/rr

print*,'go,gc,eo,ec',go,gc,e_old,e_cur
      if(go*gc < zero .or. (go*gc > zero .and. e_cur > e_old) ) then
!!$      if(.true.) then
         calc= .true.
         a0=e_old
         a1=go
         a3=-two*(e_cur-a0-a1*rr)/rr**3+(gc-a1)/rr**2
         if(abs(a3) < small) then
            a2=(gc-a1)/(two*rr)
            if(abs(a2) < small) then
               calc = .false.
               rs_limit=rr
               return
            else
               rs=-a1/(two*a2)
            end if
         else
            a2=three*(e_cur-a0-a1*rr)/(rr*rr)-(gc-a1)/rr
            b1=-two*a2; b2=four*a2*a2-twelve*a3*a1
            if(b2 >= zero) then
               b2=sqrt(b2)
               rs=(b1-b2)/(six*a3)
               if(six*a3*rs+two*a2 < zero) rs=(b1+b2)/(six*a3)
            else
               calc = .false.
               rs_limit=rr
               return
            end if
         end if
         rs_limit=rr 
         rstep=rstep*rs/rr
print*,rs,rr,rs/rr
         if(lattice_calc.or.slab_calc) then
            ku=0
            do iu=fi,n_species
               do ju=1,3
                  ku=ku+1
                  atoms_frac(iu)%r(ju)=r_old(ku)+rstep(ku)
                  if(slab_calc .and. ju==3) cycle
                  if(atoms_frac(iu)%r(ju) < zero) atoms_frac(iu)%r(ju)=atoms_frac(iu)%r(ju)+one
                  if(atoms_frac(iu)%r(ju) >= one) atoms_frac(iu)%r(ju)=atoms_frac(iu)%r(ju)-one
               end do
            end do
            if(calc_strain) then
               if(lattice_calc) then
                  strain=r_old(nvar+1:nvar+6)+rstep(nvar+1:nvar+6)
                  str=rstep(nvar+1:nvar+6)
                  call use_vect()
                  call renew_lattice(str)
               else if(slab_calc) then
                  strain_s=r_old(nvar+1:nvar+3)+rstep(nvar+1:nvar+3)
                  str=rstep(nvar+1:nvar+3)
                  call use_vect_s()
                  call renew_slab(str)
               end if
            end if
         else
            do iu=fi,n_species
               do ju=1,3
                  ku=3*(iu-1)+ju
                  atoms_cart(iu)%r(ju)=r_old(ku)+rstep(ku)
               end do
            end do
         end if

         rs_limit=abs(rs)
      else
         calc = .false.
         rs_limit=rr*1.5_r8_kind
      end if

    end subroutine line_search
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    subroutine line_search_per(search_ok)

      logical :: search_ok
      integer(i4_kind) :: iu,ju,ku
      logical :: store_g,store_h
      real(r8_kind) :: e1,e2,e3,x1,x2,x3
      real(r8_kind) :: step

      store_g=calc_gradients; store_h=calc_hessian
      calc_gradients=.false.; calc_hessian=.false.

      x1=zero; e1=e_cur

      call total_energy_and_grad(do_corr=.true.)
      x2=sqrt(dot_product(rstep,rstep))
      e2=E_total
print*,e1,e2
      if(e2 >= e1) then 
         step=half
         main1:do 
            if(lattice_calc.or.slab_calc) then
               ku=0
               do iu=fi,n_species
                  do ju=1,3
                     ku=ku+1
                     atoms_frac(iu)%r(ju)=r_old(ku)+step*rstep(ku)
                     if(slab_calc .and. ju==3) cycle
                     if(atoms_frac(iu)%r(ju) < zero) atoms_frac(iu)%r(ju)=atoms_frac(iu)%r(ju)+one
                     if(atoms_frac(iu)%r(ju) >= one) atoms_frac(iu)%r(ju)=atoms_frac(iu)%r(ju)-one
                  end do
               end do
               if(calc_strain) then
                  if(lattice_calc) then
                     strain=r_old(nvar+1:nvar+6)+step*rstep(nvar+1:nvar+6)
                     str=step*rstep(nvar+1:nvar+6)
                     call use_vect()
                     call renew_lattice(str)
                  else if(slab_calc) then
                     strain_s=r_old(nvar+1:nvar+3)+rstep(nvar+1:nvar+3)
                     str_s=step*rstep(nvar+1:nvar+3)
                     call use_vect_s()
                     call renew_slab(str_s)
                  end if
               end if
            else
               do iu=fi,n_species
                  do ju=1,3
                     ku=3*(iu-1)+ju
                     atoms_cart(iu)%r(ju)=r_old(ku)+step*rstep(ku)
                  end do
               end do
            end if
            call total_energy_and_grad(do_corr=.true.)
            x3=sqrt(dot_product(step*rstep,step*rstep))
            e3=E_total
print*,e3,'e3'
            if(e3 < e1) then
               search_ok=.true.
               rstep=step*rstep
               calc_gradients=store_g; calc_hessian=store_h
               rs_limit=x3
               calc=.true.
               return
            else
               if(x3 < small) then
                  search_ok=.false.
                  calc_gradients=store_g; calc_hessian=store_h
                  calc=.true.
                  if(lattice_calc.or.slab_calc) then
                     ku=0
                     do iu=fi,n_species
                        do ju=1,3
                           ku=ku+1
                           atoms_frac(iu)%r(ju)=r_cur(ku)
                        end do
                     end do
                     if(calc_strain) then
                        if(lattice_calc) then
                           strain=r_cur(nvar+1:nvar+6)
                           str=rstep(nvar+1:nvar+6)
                           call use_vect()
                           call renew_lattice(str)
                        else if(slab_calc) then
                           strain_s=r_cur(nvar+1:nvar+3)
                           str_s=rstep(nvar+1:nvar+3)
                           call use_vect_s()
                           call renew_slab(str_s)
                        end if
                     end if
                  else
                     do iu=fi,n_species
                        do ju=1,3
                           ku=3*(iu-1)+ju
                           atoms_cart(iu)%r(ju)=r_cur(ku)
                        end do
                     end do
                  end if
                  return
               end if
               step=step/two
            end if
         end do main1
      else
         search_ok=.true.
         calc_gradients=store_g; calc_hessian=store_h
         calc=.true.
         rs_limit=x2*1.5_r8_kind
         return
      end if
      
      calc_gradients=store_g; calc_hessian=store_h

    end subroutine line_search_per
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    function converged()

      logical :: converged
      
      real(kind=r8_kind) :: a1,a2

      converged=.false.

      rnor=sqrt(dot_product(rstep,rstep))
      a1=minval(rstep)
      a2=maxval(rstep)
      rabs=max(abs(a1),abs(a2))

      gnor=sqrt(dot_product(g_cur,g_cur))/(nvar+ns)
      a1=minval(g_cur)
      a2=maxval(g_cur)
      gabs=max(abs(a1),abs(a2))

      cosgr=(dot_product(-rstep,g_cur))/(rnor*gnor*(nvar+ns))
write(6,*)'COS',cosgr

      if( rnor <= rnorm .and. rabs <= rmax .and. &
          gnor <= gnorm .and. gabs <= gmax ) converged=.true.
      if( rnor <= small .and. i > 1 ) converged=.true.

    end function converged
    !--------------------------------------------------------------

    !--------------------------------------------------------------
    subroutine output_minim_run(iu)

      integer(i4_kind),intent(in) :: iu

      if(iu==0) then
         write(output_device,'(80("-"))')
         write(output_device,'(a74)') &
              " Iteration       Energy        Rnorm      Rmax     Gnorm      Gmax    Time"
         write(output_device,'(80("-"))')
      else
         write(output_device,'(3x,i5,3x,f17.8,4(1x,e9.3),1x,f7.2)') &
              iu,e_cur,rnor,rabs,gnor,gabs,optim_time%rt
         write(6,'(3x,i5,3x,f17.8,4(1x,e9.3),1x,f7.2)') &
              iu,e_cur,rnor,rabs,gnor,gabs,optim_time%rt
      end if

    end subroutine output_minim_run

  end subroutine minimize
  !****************************************************************

  !****************************************************************
  subroutine shutdown_hess()

    integer(i4_kind) :: status

    if(allocated(Hess)) then
       deallocate(Hess,stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed HESS deallocation")
    end if

  end subroutine shutdown_hess
  !****************************************************************

  !****************************************************************
  subroutine shutdown_H_on_slaves()

    integer(i4_kind) :: status

    deallocate(h,stat=status)
    if(status /= 0) call error_handler("MolMech: failed H deallocation on slaves")
       
  end subroutine shutdown_H_on_slaves
  !****************************************************************

  !****************************************************************
end module hess_and_opt_module




