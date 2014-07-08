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
module  population_module
  !-------------------------------------------------------------------
  !
  !  Purpose: contains various routines for the
  !           Mullikan population analysis.
  !
  !  Module called by: main_scf
  !
  !  References: See e.g. 'Einfuehrung in die theoretische Chemie'
  !              by W. Kutzelnigg, p. 553
  !
  !  Author: FN
  !  Date: 5/97
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   1/98
  ! Description: The posibility of a user defined population analysis
  !              has been added.  If n_columns>0 in the input n_columns
  !              namelists popcolumn have to follow. In this namelist the number
  !              of contributions n_contrib is specified. Then a list of contributions
  !              follows in the format <number of the atom> <angular momentum>
  !              (s=0,p=1,d=2,...)
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   6/98
  ! Description: Extension to Spin Orbit
  !              added population_spor_mulliken
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use occupation_module
  use density_data_module
  use symmetry_data_module
  use unique_atom_module
  use overlap_module
  use eigen_data_module
  use unique_atom_module
  use iounitadmin_module
  use filename_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================



  !------------ public functions and subroutines ---------------------

  public :: population_read
  public :: population_write
  public :: population_mulliken
  public :: population_spor_mulliken
  public :: population_close!(), clean up the module

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of constants and variables ---------------

  integer(kind=i4_kind),parameter     :: max_len=100
  type(intmat2), allocatable          :: column(:)

  !
  ! This holds atomic charges (as in charges and not in number of electrons).
  ! Allocated and set in population_mulliken() modified in main_gradient()
  ! and in EPE code. Thus cannot be protected.
  !
  real(kind=r8_kind), allocatable, public :: m_charge(:) ! (n_unique_atoms)

  integer(kind=i4_kind), allocatable  :: n_con(:)

  integer(kind=i4_kind), private      :: df_population_level  = 2, & !
                                         df_n_columns         = 0, & ! number of columns
                                         ! for a userdefined population analysis
                                         max_population_level = 3, &
                                         df_irrep = 0, &
                                         df_spin = 0

  logical, private                    :: df_dos_plot=.false., &
                                         df_occupation=.false.

  real(kind=r8_kind)                  :: df_eig_min=-100.0_r8_kind , df_eig_max=10.0_r8_kind

  real(kind=r8_kind)                  :: eig_min, & ! minimal
                                         eig_max !  and maximal eigenvalue for dos plots

  integer(kind=i4_kind)               :: population_level, n_columns, &
                                         irrep, & ! only used for dos plot
                                         spin, &  ! only used for spin
                                         n_contrib

  logical                             :: dos_plot, & ! make dos plot or not
                                         occupation  ! use occupation in the dos plots or not

  namelist /population/ population_level, n_columns, dos_plot, irrep, spin, &
                        occupation, eig_min, eig_max
  namelist /popcolumn/ n_contrib
!---------------------------------------------------------------------

  real(r8_kind), parameter :: zero = 0.0_r8_kind
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine population_read()
    ! read in the output level of the Mulliken population
    ! analysis.
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    use input_module
    integer :: unit, status, i, ii, alloc_stat
    !------------ Executable code ------------------------------------
    population_level = df_population_level
    n_columns = df_n_columns
    dos_plot = df_dos_plot
    irrep = df_irrep
    spin = df_spin
    eig_min = df_eig_min
    eig_max = df_eig_max
    occupation = df_occupation
    if ( input_line_is_namelist("population") ) then
       call input_read_to_intermediate
       unit = input_intermediate_unit()
       read(unit, nml=population, iostat=status)
       if (status .gt. 0) call input_error( &
            "population_read: namelist population.")
    endif
    if(population_level.lt.0.or.population_level.gt.max_population_level)then
       call input_error &
            ("population_read : sorry, this value is not supported")
       population_level = df_population_level
    endif
    if(eig_min>=eig_max) call input_error(&
         'population_read: sensless value for eig_min or eig_max')
    eig_min=eig_min/27.211652_r8_kind
    eig_max=eig_max/27.211652_r8_kind
    if(n_columns<0) call input_error(&
         'Sorry, n_columns must not be negative')
    if (n_columns>0) then
       allocate(column(n_columns),n_con(n_columns),stat=alloc_stat)
       if(alloc_stat/=0) &
            call error_handler('population_read: allocating column')
       do i=1,n_columns
          if ( input_line_is_namelist("popcolumn") ) then
             call input_read_to_intermediate
             unit = input_intermediate_unit()
             read(unit, nml=popcolumn, iostat=status)
             if (status .gt. 0) call input_error( &
                  "population_read: namelist population.")
             if(n_contrib==0) call input_error(&
                  'Sorry, n_contrib=0 does not make sense')
             allocate(column(i)%m(2,n_contrib),stat=alloc_stat)
             if(alloc_stat/=0) &
                  call error_handler('population_read: allocating contrib')
             n_con(i)=n_contrib
             call input_read_to_intermediate
             unit = input_intermediate_unit()
             read(unit,*,iostat=status) column(i)%m
             ! check input
             do ii=1,n_con(i)
                if(column(i)%m(1,ii) > N_unique_atoms .or. &
                     column(i)%m(1,ii) < 1) call input_error(&
                     'Sorry, wrong unique atom')
                if(column(i)%m(2,ii)>&
                     unique_atoms(column(i)%m(1,ii))%lmax_ob.or. &
                     column(i)%m(1,ii) < 0 ) call input_error(&
                     'Sorry, l is out of range')
             end do
          else
             call input_error &
                  ("population_read : sorry, for every column namelist popcolumn is&
                  & required")
          end if
       end do
    end if
  end subroutine population_read

  !*************************************************************

  subroutine population_write(iounit)
    !
    ! Write namelist population to iounit in input format.
    !
    use echo_input_module
    use operations_module, only: operations_echo_input_level
    implicit none
    integer, intent(in) :: iounit
    !** End of interface *****************************************

    integer             :: i

    call start("POPULATION","POPULATION_WRITE", &
         iounit,operations_echo_input_level)
    call intg("POPULATION_LEVEL",population_level,df_population_level)
    call intg("N_COLUMNS       ",n_columns,df_n_columns)
    call flag("DOS_PLOT        ",dos_plot ,df_dos_plot)
    call intg("IRREP           ",irrep,df_irrep)
    call intg("SPIN            ",spin ,df_spin)
    call flag("OCCUPATION      ",occupation ,df_occupation)
    call real("EIG_MIN         ",eig_min*27.211652_r8_kind,df_eig_min,3)
    call real("EIG_MAX         ",eig_max*27.211652_r8_kind,df_eig_max,3)

    call stop()
    do i=1,n_columns
       call start("POPCOLUMN","POPULATION_WRITE", &
            iounit,operations_echo_input_level)
       call intg("N_CONTRIB",n_con(i),0)
       write(iounit,'(30I5)') column(i)%m(:,:)
       call stop()
    end do
    if ( operations_echo_input_level == echo_level_full .or. &
         population_level .ne. df_population_level ) &
         write(iounit,'(a30)') "    # legal values are 0, 1, 2"

  end subroutine population_write

  subroutine population_close()
    implicit none
    ! *** end of interface ***

    integer :: memstat

    ! FIXME: anything else? Look at the many global vars ...
    if ( allocated(m_charge) ) then
        deallocate(m_charge, stat=memstat)
        ASSERT(memstat==0)
    endif
  end subroutine population_close

  subroutine population_mulliken()
    !  Purpose: ..
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use comm, only: comm_rank
    use orbitalprojection_module, only: orbitalprojection_ob, &
         orbitalprojection_print
    use iounitadmin_module, only: output_unit
    !------------ Declaration of local variables ---------------------

    integer(kind=i4_kind)                :: i,alloc_stat,il,is,&
         up,lo,ua,m_help,m,num,n,nu,nl,i_col,i_cont
    integer(kind=i4_kind)                :: i_ir,mu
    integer(kind=i4_kind),allocatable    :: index_arr(:)
    real(kind=r8_kind)                   :: summ
    character                  :: symbol(0:11)
    character(len=max_len)     :: ua_char,char_help
    external error_handler
    integer :: iounit_dos ! iounit for dos plots
    integer(kind=i4_kind)               :: n_tot,n_sum
    type(arrmat2),allocatable           :: contribute_charge(:)
    type(arrmat3),allocatable           :: pop_col(:)

    character(len=max_len), allocatable :: title_cont(:)
    character(len=max_len), allocatable :: title_l(:)
    character(len=max_len), allocatable :: title_sum(:)

    real(kind=r8_kind), allocatable :: contribute(:,:)
    real(kind=r8_kind), allocatable :: contribute_l(:,:)
    real(kind=r8_kind), allocatable :: sum_l(:)
    real(kind=r8_kind), allocatable :: sum_ocup_l(:)
    real(kind=r8_kind), allocatable :: sum_cont(:)
    real(kind=r8_kind), allocatable :: sum_ocup_cont(:)
    real(kind=r8_kind), allocatable :: sum_total(:)
    real(kind=r8_kind), allocatable :: sum_total_updn(:,:)
    real(kind=r8_kind), allocatable :: atomic_charge(:)


    ! FIXME: overlap not availble?
    if (comm_rank() /= 0) return

    symbol = (/'S','P','D','F','G','H','I','J','K','L','M','N'/)


    if (population_level == 0) return

    if(dos_plot) then ! open file for output of dos_plots
       iounit_dos = openget_iounit(trim(outfile("dos.dat")), &
            form="FORMATTED",status="REPLACE",action="WRITE")
    end if

    ! allocate and read overlap from file
    ! MOVED BEFORE SCF LOOP: call read_overlap()
    ASSERT(allocated(overlap))

    ! it is assumed that 'densmat' is still allocated
    ! -> see 'main_scf'
    ASSERT(allocated(densmat))


    ! First calculate the atomic charges -------------------------
    ! calculate diagonal elements of P*S
    ! P = densmat, S = overlap
    allocate(contribute_charge(ssym%n_irrep), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    allocate(atomic_charge(N_unique_atoms), m_charge(N_unique_atoms), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    do i=1,ssym%n_irrep
       allocate(contribute_charge(i)%m(ssym%dim(i), ssym%n_spin), STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    enddo

    atomic_charge = 0.0_r8_kind
    do i=1,ssym%n_irrep
       contribute_charge(i)%m=0.0_r8_kind
       do is=1,ssym%n_spin
          do m=1,ssym%dim(i)
             summ=0.0_r8_kind
             do m_help=1,ssym%dim(i)
                summ = summ + densmat(i)%m(m,m_help,is)*overlap(i)%m(m_help,m)
             enddo
             contribute_charge(i)%m(m,is) = summ

          enddo
       enddo
    enddo

    do ua=1,N_unique_atoms
       do i = 1,ssym%n_irrep
          lo=orbitalprojection_ob(i,0,ua)
          if (ua.eq.N_unique_atoms) then
             up=huge(1_i4_kind)
          else
             up=orbitalprojection_ob(i,0,ua+1)
          endif
          do is=1,ssym%n_spin
             do m=1,ssym%dim(i)
                if (m.ge.lo.and.m.lt.up) then
                   atomic_charge(ua) = atomic_charge(ua)+ &
                        contribute_charge(i)%m(m,is)
                endif
             enddo
          enddo
       enddo
       atomic_charge(ua) = atomic_charge(ua)/unique_atoms(ua)%N_equal_atoms
    enddo


    ! prepare the array 'sum_total' which has to contain the sums over Irreps and Spin
    ! of all contributions -------------------------------------------------
    n_sum = 0
    do ua=1,N_unique_atoms
       n_sum = n_sum + unique_atoms(ua)%lmax_ob + 1
    enddo
    allocate(sum_total(n_sum), sum_total_updn(n_sum, ssym%n_spin), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    allocate(index_arr(n_sum), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    allocate(title_sum(n_sum), STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    sum_total = 0.0_r8_kind
    sum_total_updn = 0.0_r8_kind
    index_arr = 0_i4_kind
    is_: do is=1,ssym%n_spin
       i_ir_: do i_ir=1,ssym%n_irrep

          n = 0
          do ua=1,N_unique_atoms
             do il=0,unique_atoms(ua)%lmax_ob
                n = n + 1
                index_arr(n) = orbitalprojection_ob(i_ir,il,ua)
             enddo
          enddo

          mu_: do mu=1,ssym%dim(i_ir)
             n=1
             ua_: do ua=1,N_unique_atoms
                il_: do il=0,unique_atoms(ua)%lmax_ob
                   lo = index_arr(n)
                   if (n==n_sum) then
                      up = huge(1_i4_kind)
                   else
                      up = index_arr(n+1)
                   endif

                   if ((mu.ge.lo).and.(mu.lt.up)) then
                      sum_total(n) = sum_total(n) + &
                           contribute_charge(i_ir)%m(mu,is)
                      sum_total_updn(n,is) = sum_total_updn(n,is) + &
                           contribute_charge(i_ir)%m(mu,is)
                   endif
                   write(ua_char,*)ua
                   title_sum(n)  = ' '//&
                        & trim(adjustl(ua_char))//' '//&
                        & trim(adjustl(symbol(il)))
                   n = n+1
                enddo il_
             enddo ua_
          enddo mu_
       enddo i_ir_
    enddo is_
    ! ---------------------------------------------------------------------

    write(output_unit,*)
    write(output_unit,*)
    write(output_unit,fmt='(46X,"*** POPULATION ANALYSIS***")')
    write(output_unit,*)
    write(output_unit,*)

    if (population_level == 1) then
       call popout(1, 1, .false., ssym, &
                title_cont, title_l, title_sum,&
                contribute, contribute_l, sum_l, sum_ocup_l, &
                sum_cont, sum_ocup_cont, sum_total, sum_total_updn, atomic_charge)
       goto 999 ! finilize and exit
    endif

    if(n_columns>0) then ! do some allocations
       allocate(pop_col(ssym%n_irrep),stat=alloc_stat)
       if(alloc_stat/=0) call error_handler(&
            'POPULATION_MULLIKEN: allocating pop_col')
       do i_ir=1,ssym%n_irrep
          allocate(pop_col(i_ir)%m(ssym%dim(i_ir),ssym%n_spin,n_columns), &
               stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               'POPULATION_MULLIKEN: allocating pop_col%m')
          pop_col(i_ir)%m=0.0_r8_kind
       end do
    end if

    do is=1,ssym%n_spin
       ! This part prints out information on the contribution of
       ! each CONTRACTED basis function
       ! to the respective MO`s.
       do i_ir=1,ssym%n_irrep

          ! First prepare arrays containing the labels for each contracted
          ! basis function : 'title_cont' and 'title_l'
          allocate(title_cont(ssym%dim(i_ir)),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          n=0
          nl=1
          n_tot = 0  ! n_tot determines the number of columns for contributions
          !          ! sorted by angular momenta
          do ua=1,N_unique_atoms
             il_1 : do il=0,unique_atoms(ua)%lmax_ob
                num = (unique_atoms(ua)%l_ob(il)%N_uncontracted_fcts + &
                     unique_atoms(ua)%l_ob(il)%N_contracted_fcts) * &
                     unique_atoms(ua)%symadapt_partner(i_ir,il)%N_independent_fcts
                write(ua_char,*)ua ! use an internal file to convert int to char
                if (n /= ssym%dim(i_ir) ) then
                   title_cont(n+1:n+num) = trim(adjustl(ua_char))//' '//&
                        & trim(adjustl(symbol(il)))
                   n = n + num
                endif
                if (num /= 0 ) n_tot=n_tot+1
             enddo il_1
          enddo

          ! determine 'title_l': every time the contents of 'title_cont' changes
          ! 'title_l' is set to the respective value. Improvements welcome.
          allocate(title_l(n_tot),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          char_help = ' '
          nl=1
          do n=1,ssym%dim(i_ir)
             if ( verify(title_cont(n),char_help) /= 0 ) then
                title_l(nl) = title_cont(n)
                nl = nl + 1
             endif
             char_help = title_cont(n)
          enddo
          ! -------------------------------------------------



          allocate(contribute(ssym%dim(i_ir),ssym%dim(i_ir)),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(contribute_l(ssym%dim(i_ir),n_tot),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(sum_l(n_tot),sum_ocup_l(n_tot),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(sum_cont(ssym%dim(i_ir)),sum_ocup_cont(ssym%dim(i_ir)),STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          contribute = 0.0_r8_kind
          contribute_l = 0.0_r8_kind
          sum_l = 0.0_r8_kind
          ! 'contribute' contains the contributions of the contracted basis
          ! functions excluding the occupation number
          ! 'contribute_l' contains the contributions sorted by angular
          ! momenta
          WARN('O(N^3) serial step')
          do mu=1,ssym%dim(i_ir)
             do nu=1,ssym%dim(i_ir)
                do n=1,ssym%dim(i_ir)
                   contribute(n,mu) = contribute(n,mu) + eigvec(i_ir)%m(mu,n,is)*&
                        eigvec(i_ir)%m(nu,n,is)*overlap(i_ir)%m(nu,mu)
                enddo
             enddo
             sum_ocup_cont(mu) = contribute_charge(i_ir)%m(mu,is)
          enddo
          if(n_columns>0) call sumup_populations(ssym%dim(i_ir), i_ir, is, &
                contribute, pop_col)

![TS:
      !    if(dos_plot.and.((irrep==0.and.spin==0).or.(irrep==i_ir.and.spin==is))) &
      !         call make_dos_plot(ssym%dim(i_ir),i_ir,is)

          if(dos_plot) then
             if(irrep==0 .and. spin==0) then
                call make_dos_plot(ssym%dim(i_ir), i_ir, is, ssym, &
                        pop_col, iounit_dos)
             elseif(irrep==i_ir .and. spin==is) then
                call make_dos_plot(ssym%dim(i_ir), i_ir, is, ssym, &
                        pop_col, iounit_dos)
             elseif(irrep==0 .and. spin==is) then
                call make_dos_plot(ssym%dim(i_ir), i_ir, is, ssym, &
                        pop_col, iounit_dos)
             elseif(irrep==i_ir .and. spin==0) then
                call make_dos_plot(ssym%dim(i_ir), i_ir, is, ssym, &
                        pop_col, iounit_dos)
             endif
          endif

          do mu=1,ssym%dim(i_ir)
             n=0
             nl=1
             do ua=1,N_unique_atoms
                il_2 : do il=0,unique_atoms(ua)%lmax_ob
                   num = (unique_atoms(ua)%l_ob(il)%N_uncontracted_fcts + &
                        unique_atoms(ua)%l_ob(il)%N_contracted_fcts) * &
                        unique_atoms(ua)%symadapt_partner(i_ir,il)%N_independent_fcts
                   if (num/=0) then
                      contribute_l(mu,nl) = sum(contribute(mu,n+1:n+num))
                      sum_ocup_l(nl) = sum(contribute_charge(i_ir)%m(n+1:n+num,is))
                      nl=nl+1
                      n=n+num
                   endif
                   if (n==ssym%dim(i_ir) ) exit il_2
                enddo il_2
             enddo
          enddo

          do n=1,n_tot
             sum_l(n) = sum(contribute_l(:,n))
          enddo
          do mu=1,ssym%dim(i_ir)
             sum_cont(mu) = sum(contribute(mu,:))
          enddo


          call popout(is, i_ir, .true., ssym, &
                title_cont, title_l, title_sum,&
                contribute, contribute_l, sum_l, sum_ocup_l, &
                sum_cont, sum_ocup_cont, sum_total, sum_total_updn, atomic_charge)

          deallocate(contribute,title_cont,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(contribute_l,title_l,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(sum_l,sum_ocup_l,STAT=alloc_stat)
          ASSERT(alloc_stat==0)
          deallocate(sum_cont,sum_ocup_cont,STAT=alloc_stat)
          ASSERT(alloc_stat==0)

       enddo
    enddo
    ! ------------------------------------------------------------------------
    call popout(is, 1, .false., ssym, &
                title_cont, title_l, title_sum,&
                contribute, contribute_l, sum_l, sum_ocup_l, &
                sum_cont, sum_ocup_cont, sum_total, sum_total_updn, atomic_charge)
    if(n_columns>0) then
       write(output_unit,'(A)') 'Printing user defined population analysis'
       write(output_unit,'(A20)') 'Number of columns:',n_columns
       write(output_unit,'(A)') 'Contributions to every column:'
       do i_col=1,n_columns
          write(output_unit,'(A7,I3)') 'column:',i_col
          do i_cont=1,n_con(i_col)
             if(column(i_col)%m(2,i_cont)==0) &
                  write(output_unit,2000) column(i_col)%m(1,i_cont),'s'
             if(column(i_col)%m(2,i_cont)==1) &
                  write(output_unit,2000) column(i_col)%m(1,i_cont),'p'
             if(column(i_col)%m(2,i_cont)==2) &
                  write(output_unit,2000) column(i_col)%m(1,i_cont),'d'
             if(column(i_col)%m(2,i_cont)==3) &
                  write(output_unit,2000) column(i_col)%m(1,i_cont),'f'
             if(column(i_col)%m(2,i_cont)==4) &
                  write(output_unit,2000) column(i_col)%m(1,i_cont),'g'
             if(column(i_col)%m(2,i_cont)==5) &
                  write(output_unit,2000) column(i_col)%m(1,i_cont),'h'
             if(column(i_col)%m(2,i_cont)==6) &
                  write(output_unit,2000) column(i_col)%m(1,i_cont),'i'
             if(column(i_col)%m(2,i_cont)>6) &
                  write(output_unit,'(2I3)') &
                  column(i_col)%m(1,i_cont),column(i_col)%m(2,i_cont)
          end do
       end do
       call occupation_print_popspectrum(pop_col)
    end if
2000 format(('Atom and angular momentum:',I3,A))

    ! dont forget to invalidate overlap before next geo-loop:
    ! MOVED TO AFTER SCF LOOP: call dealloc_overlap()

    if(n_columns>0) then ! do some deallocations
       ! FIXME: make arrmat3 components allocatable:
       do i_ir=1,ssym%n_irrep
          deallocate(pop_col(i_ir)%m, &
               stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               'POPULATION_MULLIKEN: deallocating pop_col%m')
       end do
       deallocate(pop_col,stat=alloc_stat)
       if(alloc_stat/=0) call error_handler(&
            'POPULATION_MULLIKEN: deallocating pop_col')

       ! FIXME: move this to population_close():
       do i=1,n_columns
          deallocate(column(i)%m,stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               'POPULATION_MULLIKEN: deallocating column%m')
       end do
       deallocate(column,n_con,stat=alloc_stat)
       if(alloc_stat/=0) call error_handler(&
            'POPULATION_MULLIKEN: deallocating column')
    end if

    if(dos_plot) then
       call returnclose_iounit(iounit_dos,status="KEEP")
    endif

999 CONTINUE ! clean up and exit
  end subroutine population_mulliken

  !*************************************************************
  subroutine population_spor_mulliken()
    !  Purpose: ..
    !** End of interface *****************************************
    !------------ Modules used ------------------- ---------------
    use comm, only: comm_rank
    use dimensions, only: ShellDim => BasDimSpor, UADim => UABasDimSpor, &
         IrrDim => IrrBasDimSpor, uaL_proj_dims
    use error_module, only: warn, error
    use datatype, only: arrmat1_c
    use symmetry_data_module, only: ssym
    use density_data_module,  only: densmat_real, densmat_imag
    use unique_atom_module,   only: unique_atoms, N_unique_atoms
    use overlap_module,       only: overlap_real, overlap_imag
    use eigen_data_module,    only: eigvec_real, eigvec_imag
    !------------ Declaration of formal parameters ---------------
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)               ::&
         & irr,m,n,mu,nu,ibas,isub,ua,L,uaL,iuaL,uaL_dim
    integer(kind=i4_kind)               :: alloc_stat
    real(kind=r8_kind)                  :: summ_real,summ_imag
    character(len=4)                    :: symbol(0:11)
    character(len=max_len)      :: ua_char
#if FPP_DEBUG > 1
    use error_module, only: MyID
    character(len=*), parameter :: ME="pm/population_spor_mulliken"
#else
    character(len=*), parameter :: ME="pm"
#endif
    external error_handler
!--
    type(arrmat1_c),allocatable         :: contribute_charge(:)
    character(len=max_len),allocatable  :: title_l(:), title_sum(:)
    complex(c16_kind),allocatable       :: contribute(:,:), contribute_l(:,:),&
         & sum_l(:), sum_ocup_l(:), sum_total(:)
    complex(kind=c16_kind),allocatable  :: atomic_charge(:)
    integer(kind=i4_kind)               :: n_tot,n_sum


    ! FIXME: overlap not availble?
    if (comm_rank() /= 0) return

!!$    symbol = (/'S1/2','P1/2','P3/2','D3/2','D5/2','F5/2','F7/2','G7/2','G9/2'/)
    symbol = (/'S','P','D','F','G','H','I','J','K','L','M','N'/)


    if (population_level == 0 ) return

    if ( ssym%n_proj_irrep .ne. size(IrrDim) )&
         & call error(ME//": n_irr?")

    ! allocate and read overlap from file
    ! MOVED BEFORE SCF LOOP: call read_overlap()
    ASSERT(allocated(overlap_real))
    ASSERT(allocated(overlap_imag))

    ! it is assumed that 'densmat' is still allocated
    ! -> see 'main_scf'
    ASSERT(allocated(densmat_real))
    ASSERT(allocated(densmat_imag))


    ! First calculate the atomic charges -------------------------
    ! calculate diagonal elements of P*S
    ! P = densmat, S = overlap
    allocate(&
         & contribute_charge(ssym%n_proj_irrep),&
         & atomic_charge(N_unique_atoms),&
         & STAT=alloc_stat)
    if (alloc_stat/=0)&
         & call error(ME//": allocation (1) failed")

    do irr=1,ssym%n_proj_irrep
       allocate(&
            & contribute_charge(irr)%m(ssym%dim_proj(irr)),&
            & STAT=alloc_stat)
       if (alloc_stat/=0)&
            & call error(ME//": allocation (2) failed")
    enddo

    do irr = 1, size(IrrDim)
       DPRINT ME,': irrep ',irr
       contribute_charge(irr)%m = zero
       DPRINT ME,': now cycle m=1,',IrrDim(irr),'(==',ssym%dim_proj(irr),')'
       do m = 1, IrrDim(irr)

          summ_real = zero
          summ_imag = zero
          do n = 1, IrrDim(irr)
             summ_real = summ_real&
                  & + densmat_real(irr)%m(m,n) * overlap_real(irr)%m(m,n)&
                  & - densmat_imag(irr)%m(m,n) * overlap_imag(irr)%m(m,n)
             summ_imag = summ_imag&
                  & + densmat_real(irr)%m(m,n) * overlap_imag(irr)%m(m,n)&
                  & + densmat_imag(irr)%m(m,n) * overlap_real(irr)%m(m,n)
          enddo
          contribute_charge(irr)%m(m) = cmplx(summ_real,summ_imag,r8_kind)
       enddo
#if FPP_DEBUG >= 5
       DPRINT ME,': charge contribution CHARGE(irr=',irr,')=',sum(contribute_charge(irr)%m)
#endif
    enddo

    DPRINT ME,': now atomic charges ...'
    atomic_charge = zero
    do irr = 1, size(IrrDim)
       DPRINT ME,': irr=',irr
       ibas = 0
       do ua = 1, size(UADim, 2) ! n_ua
          DPRINT ME,': ua=',ua,' dim(irr,ua)=', UADim(irr, ua)
          do isub = 1, UADim(irr, ua) ! run over UA_space
             ibas = ibas + 1
             atomic_charge(ua) = atomic_charge(ua)&
                  & + contribute_charge(irr)%m(ibas)
          enddo
          DPRINT ME,': ibas=',ibas
       enddo
    enddo
    do ua = 1, size(UADim, 2)
       atomic_charge(ua) = atomic_charge(ua)/unique_atoms(ua)%N_equal_atoms
       DPRINT ME,': atomic_charge(',ua,') =',atomic_charge(ua)
    enddo
    DPRINT ME,': done atomic charges'

    ! prepare the array 'sum_total' which has to contain the sums over Irreps and Spin
    ! of all contributions -------------------------------------------------
    DPRINT ME,': now continue to UA_L shells ...'
    n_sum = 0
    do ua=1,N_unique_atoms
       n_sum = n_sum + ( unique_atoms(ua)%lmax_ob + 1 )
    enddo
    DPRINT ME,': n_sum = ',n_sum,'; do alloc ..'
    allocate(&
         & sum_total(n_sum),&
         & title_sum(n_sum),&
         & STAT=alloc_stat)
    if ( alloc_stat .ne. 0 )&
         & call error(ME//": alloc sum_total ... failed")
    DPRINT ME,': diving into loops ...'
    sum_total = zero
    do irr = 1, size(IrrDim)
       uaL = 0 ! runs over UA_L = 1_S, 1_P ... 2_S, 2_P, ...
       ibas = 0
       do ua = 1, size(UADim, 2) ! == N_unique_atoms
          do L=0,ShellDim(ua)%lmax
             uaL = uaL + 1
             do isub=1,ShellDim(ua)%LM(L,irr) ! runs over UA_L_subspace
                ibas = ibas + 1
                sum_total(uaL) = sum_total(uaL)&
                     & + contribute_charge(irr)%m(ibas)
             enddo
          enddo
       enddo
    enddo
    ! type symbols:
    uaL = 0
    do ua=1,size(ShellDim)
       do L =0,ShellDim(ua)%lmax
          uaL = uaL + 1
          write(ua_char,*) ua
          title_sum(uaL)  = ' '//trim(ua_char)//' '//trim(symbol(L))
          DPRINT ME,': ',trim(title_sum(uaL)),' = ',sum_total(uaL)
       enddo
    enddo
    DPRINT ME,': done counting ...'
    ! ---------------------------------------------------------------------

    if(population_level.gt.1)then
       ! This part prints out information on the contributions of
       ! each UA_L shell to MO`s
       DPRINT ME,': now start "MOs in rows / L values in columns" ...'
       do irr = 1, size(IrrDim)
          DPRINT ME,': irrep ',irr

          DPRINT ME,': calculate contributions ...'
          allocate(contribute(IrrDim(irr), IrrDim(irr)), STAT=alloc_stat)
          if(alloc_stat.ne.0) call error(ME//": alloc contribute failed")
          contribute = zero
          ! 'contribute' contains the contributions of the contracted basis
          ! functions excluding the occupation number
          ! 'contribute_l' contains the contributions sorted by angular
          ! momenta
          do mu = 1, IrrDim(irr) ! == ssym%dim_proj(irr)
             do nu = 1, IrrDim(irr) ! == ssym%dim_proj(irr)
                do n = 1, IrrDim(irr) ! == ssym%dim_proj(irr)
                   contribute(n,mu) = contribute(n,mu) +&
                        cmplx( eigvec_real(irr)%m(mu,n), - eigvec_imag(irr)%m(mu,n),r8_kind)*&
                        cmplx( eigvec_real(irr)%m(nu,n),   eigvec_imag(irr)%m(nu,n),r8_kind)*&
                        cmplx(overlap_real(irr)%m(mu,nu), overlap_imag(irr)%m(mu,nu),r8_kind)
                enddo
             enddo
          enddo



          DPRINT ME,': count how many are there ...'
          n_tot = count(uaL_proj_dims(:,irr).gt.0) ! non-empty UA_L_subspaces
          DPRINT ME,': ',n_tot,' uaLs are not empty in irrep ',irr

          DPRINT ME,': determine their names ...'

          allocate(title_l(n_tot),STAT=alloc_stat)
          call error(alloc_stat,ME//': alloc title_l failed')

          uaL  = 0 ! counts all UA_L`s
          iuaL = 0 ! counts non-empty UA_L`s
          do ua=1,size(ShellDim) ! n_ua
             do L=0,ShellDim(ua)%lmax
                uaL = uaL + 1
                if(uaL_proj_dims(uaL,irr).eq.0) cycle
                iuaL = iuaL + 1
                write(ua_char,*) ua
                title_l(iuaL) = trim(adjustl(ua_char))//' '//trim(symbol(L))
             enddo
          enddo

          DPRINT ME,': sort contributions over uaL`s ...'

          allocate(&
               & contribute_l(IrrDim(irr),n_tot),&
               & sum_l(n_tot),sum_ocup_l(n_tot),&
               & STAT=alloc_stat)
          if(alloc_stat.ne.0) call error(ME//": alloc contribute_l, ... failed")

          do m = 1, IrrDim(irr)
             n    = 0 ! as ibas above, runs over 1..IrrDim(irr)
             iuaL = 0 ! counts non-empty UA_L`s
             do uaL=1,size(uaL_proj_dims,1) ! == uaL_max from dimensions.f90
                uaL_dim = uaL_proj_dims(uaL,irr)
                if(uaL_dim.eq.0) cycle
                iuaL = iuaL + 1
                contribute_l(m,iuaL) = sum(contribute(m,n+1:n+uaL_dim))
                n = n + uaL_dim
             enddo
          enddo

          DPRINT ME,': compute sum_ocup_l ...'

          n    = 0 ! runs over basis: 1..IrrDim(irr)
          iuaL = 0 ! non-empty UA_L''s
          do uaL=1,size(uaL_proj_dims,1) ! == uaL_max from dimensions.f90
             uaL_dim = uaL_proj_dims(uaL,irr)
             if(uaL_dim.eq.0) cycle
             iuaL = iuaL + 1
             sum_ocup_l(iuaL) = sum(contribute_charge(irr)%m(n+1:n+uaL_dim))
             n = n + uaL_dim
          enddo

          DPRINT ME,': compute sum_l ...'

          do iuaL=1,n_tot
             sum_l(iuaL) = sum(contribute_l(:,iuaL))
          enddo

          DPRINT ME,': call popout_spor(irr,.true.) ...'

          call popout_spor(irr,.true.)

          DPRINT ME,': clean up ...'

          deallocate(&
               & contribute,contribute_l,&
               & title_l,sum_l,sum_ocup_l,&
               & STAT=alloc_stat)
          if(alloc_stat.ne.0) call error(ME//": dealloc contribute, ... failed ")

       enddo
    endif! population_level.gt.1
    ! ------------------------------------------------------------------------
    DPRINT ME,': call popout_spor(1,.false.)'

    call popout_spor(1,.false.)

    DPRINT ME,': deallocate ...'

    deallocate(sum_total,STAT=alloc_stat)
    if (alloc_stat/=0)&
         & call error("POPULATION_MULLIKEN: deallocation (1) failed")

    deallocate(title_sum,STAT=alloc_stat)
    if (alloc_stat/=0)&
         & call error("POPULATION_MULLIKEN: deallocation (3) failed")

    do irr=1,ssym%n_proj_irrep
       deallocate(contribute_charge(irr)%m,STAT=alloc_stat)
       if (alloc_stat/=0)&
            & call error ("POPAN : deallocation (2) failed")
    enddo
    deallocate(contribute_charge,STAT=alloc_stat)
    if (alloc_stat/=0)&
         & call error ("POPAN : deallocation (3) failed")

    deallocate(atomic_charge,STAT=alloc_stat)
    if (alloc_stat/=0)&
         & call error ("POPAN: deallocation (1) failed")

    ! dont forget to invalidate overlap before next geo-loop:
    ! MOVED TO AFTER SCF LOOP: call dealloc_overlap()

    DPRINT ME,': exit'

  contains
    subroutine popout_spor(i_ir,flag)
      ! Purpose : output interface to the Mulliken Population
      !           analysis - 'popan' should not be overcrowded
      !           with output stuff.
      ! Input :  is        actual spin
      !          i_ir      actual irrep
      !          flag      true: output of contributions of
      !                          contracted basis fct.s as well
      !                          ass sorted by angular momenta
      !                    false: only the sum contributions
      !                          over irreps and spins and
      !                          the atomic charges are printed
      !
      ! subroutine called by: 'popan'
      ! declaration of formal parameters ------------------------
      use iounitadmin_module, only: output_unit
      use eigen_data_module,  only: eigval
      integer, intent(in) :: i_ir
      logical, intent(in) :: flag  ! determines output
      !** End of interface *****************************************
      ! --- declaration of local variables ----------------------
      integer(kind=i4_kind)    :: mu,n,nl
      real(r8_kind), parameter :: SMALL = 1.0E-08_r8_kind
      !------------ Executable code --------------------------------

      if (flag) then
         ! --- angular momenta ---
         write(output_unit,10000)
         write(output_unit,*)' Interesting Molecule' ! This should be the title
         write(output_unit,13000)
         write(output_unit,15000)i_ir,adjustl(trim(ssym%name_proj(i_ir)))
         write(output_unit,13000)
         write(output_unit,50000)
         write(output_unit,48000)
         write(output_unit,13000)

         write(output_unit,20000)(title_l(n),n=1,n_tot)
#ifdef FPP_POPAN_CMPLX
         do mu=1,ssym%dim_proj(i_ir)
            write(output_unit,25000)mu,eigval(i_ir)%m(mu,1)*27.211652_r8_kind,&
                 (contribute_l(mu,nl),nl=1,n_tot)
         enddo
         write(output_unit,30000)(sum_l(nl),nl=1,n_tot)
         write(output_unit,32000)(sum_ocup_l(nl),nl=1,n_tot)
#else
#ifndef FPP_NODEBUG
         if(any( abs(aimag(contribute_l)) .gt. SMALL ) )then
WARN("complex charges, recompile with FPP_POPAN_CMPLX")
            print *,ME//': maxval(abs(aimag(contribute_l)))=',maxval(abs(aimag(contribute_l)))
         endif
#endif
         do mu=1,ssym%dim_proj(i_ir)
            write(output_unit,25000)mu,eigval(i_ir)%m(mu,1)*27.211652_r8_kind,&
                 (real(contribute_l(mu,nl)),nl=1,n_tot)
         enddo
         write(output_unit,30000)(real(sum_l(nl)),nl=1,n_tot)
         write(output_unit,32000)(real(sum_ocup_l(nl)),nl=1,n_tot)
#endif
         write(output_unit,12000)
      else
         write(output_unit,10000)
         write(output_unit,*)' Interesting Molecule' ! This should be the title
         write(output_unit,13000)
         write(output_unit,90000)
         write(output_unit,92000)
         write(output_unit,13000)
         write(output_unit,94000)(title_sum(n),n=1,n_sum)
#ifdef FPP_POPAN_CMPLX
         write(output_unit,96000)(sum_total(n),n=1,n_sum)
         write(output_unit,80000)
         write(output_unit,81000)(n,n=1,N_unique_atoms)
         write(output_unit,82000)(atomic_charge(n),n=1,N_unique_atoms)
         write(output_unit,83000)(unique_atoms(n)%Z,n=1,N_unique_atoms)
         write(output_unit,84000)(unique_atoms(n)%Z-atomic_charge(n),n=1,N_unique_atoms)
#else
         write(output_unit,96000)(real(sum_total(n)),n=1,n_sum)
         write(output_unit,80000)
         write(output_unit,81000)(n,n=1,N_unique_atoms)
         write(output_unit,82000)(real(atomic_charge(n)),n=1,N_unique_atoms)
         write(output_unit,83000)(unique_atoms(n)%Z,n=1,N_unique_atoms)
         write(output_unit,84000)(unique_atoms(n)%Z-real(atomic_charge(n)),n=1,N_unique_atoms)
#endif

      endif

      ! formats ------------------------
10000 FORMAT(46X,'*** POPULATION ANALYSIS***')
15000 FORMAT(43X,'IRREDUCIBLE REPRESENTATION:',I2,1X,A6)
50000 FORMAT(40X,'CONTRIBUTIONS TO THE MOLECULAR ORBITALS')
48000 FORMAT(43X,'MOS IN ROWS / L VALUES IN COLUMNS')
12000 FORMAT('1')
13000 FORMAT('0')
#ifdef FPP_COMPAT
20000 FORMAT('     EIGENVALUE',5X,14(4X,A8,4X),/: (20X,14(4X,A8,4X)))
25000 FORMAT(1X,I3,1X,F10.3,4X,14(1X,F7.3),/:(19X,14(1X,F7.3)))
30000 FORMAT('0  SUM        ',5X,14(1X,F7.3),/:(17X,14(1X,F7.3)))
32000 FORMAT('0 WITH OCCUP. ',5X,14(1X,F7.3),/:(17X,14(1X,F7.3)))
#else
#ifdef FPP_POPAN_CMPLX
20000 FORMAT('        EIGENVALUE',5X,14(4X,A8,4X),/: (20X,14(4X,A8,4X)))
25000 FORMAT(1X,I3,1X,F13.3,4X,14(1X,F7.3),/:(19X,14(1X,F7.3)))
30000 FORMAT('0 SUM            ',5X,14(1X,F7.3),/:(17X,14(1X,F7.3)))
32000 FORMAT('0 WITH OCCUP.    ',5X,14(1X,F7.3),/:(17X,14(1X,F7.3)))
#else
20000 FORMAT('        EIGENVALUE',7X,14(A8),/: (20X,14(A8,1X)))
25000 FORMAT(1X,I3,1X,F13.3,4X,14(1X,F7.3),/:(19X,14(1X,F7.3)))
30000 FORMAT('0 SUM            ',5X,14(1X,F7.3),/:(17X,14(1X,F7.3)))
32000 FORMAT('0 WITH OCCUP.    ',5X,14(1X,F7.3),/:(17X,14(1X,F7.3)))
#endif
#endif
90000 FORMAT(40X,'         SUMMARY OVER IRREPS     ')
92000 FORMAT(40X,' CONTRIBUTIONS OF THE ATOMIC ORBITALS  ')
#ifdef FPP_POPAN_CMPLX
94000 FORMAT('               ',6X,14(6X,A8,6X),/: (21X,14(6X,A8,6X)))
96000 FORMAT(19X,14(1X,F9.3),/:(19X,14(1X,F9.3)))
80000 FORMAT('1',50X,'ATOMIC CHARGES'//)
81000 FORMAT(' NUCLEUS:          ',14(3X,I2,3X),/:(20X,14(3X,I2,3X)))
82000 FORMAT(' ELECTRONIC CHARGE:',14F8.3,/:(20X,14F8.3))
83000 FORMAT(' NUCLEAR CHARGE:   ',14(1X,F7.3),/:(20X,14(1X,F7.3)))
84000 FORMAT(' ATOMIC CHARGE:    ',14(1X,F12.8),/:(20X,14(1X,F12.8)))
#else
94000 FORMAT(19X,14(2X,A8),/: (21X,14(2X,A8)))
96000 FORMAT(19X,14(1X,F9.3),/:(19X,14(1X,F9.3)))
80000 FORMAT('1',50X,'ATOMIC CHARGES'//)
81000 FORMAT(' NUCLEUS:          ',14(3X,I2,4X),/:(20X,14(3X,I2,4X)))
82000 FORMAT(' ELECTRONIC CHARGE:',14F8.3,/:(20X,14F8.3))
83000 FORMAT(' NUCLEAR CHARGE:   ',14F8.3,/:(20X,14F8.3))
84000 FORMAT(' ATOMIC CHARGE:    ',14F8.3,/:(20X,14F8.3))
#endif

    end subroutine popout_spor

  end subroutine population_spor_mulliken

  !*************************************************************

  !*************************************************************
  subroutine popout (&
       is, &
       i_ir, &
       flag, &
       ssym, &
       title_cont, &            ! ???
       title_l, &               ! ???
       title_sum, &
       contribute, &            ! ???
       contribute_l, &          ! ???
       sum_l, &                 ! ???
       sum_ocup_l, &            ! ???
       sum_cont, &              ! ???
       sum_ocup_cont, &         ! ???
       sum_total, &
       sum_total_updn, &
       atomic_charge)
    ! Purpose : output interface to the Mulliken Population
    !           analysis - 'POPULATION_MULLIKEN' should not be
    !           overcrowded with output stuff.
    ! Input :  is        actual spin
    !          i_ir      actual irrep
    !          flag      true: output of contributions of
    !                          contracted basis fct.s as well
    !                          ass sorted by angular momenta
    !                    false: only the sum contributions
    !                          over irreps and spins and
    !                          the atomic charges are printed
    !
    ! subroutine called by: 'POPULATION_MULLIKEN'
    ! declaration of formal parameters ------------------------
    use iounitadmin_module, only: output_unit
    use operations_module, only: operations_core_density
    implicit none
    integer(kind=i4_kind),intent(in)  :: is, i_ir
    logical,intent(in)  :: flag  ! determines output
    type(sym), intent(in) :: ssym

    character(len=max_len), intent(in) :: title_cont(:)
    character(len=max_len), intent(in) :: title_l(:)
    character(len=max_len), intent(in) :: title_sum(:)

    real(kind=r8_kind), intent(in) :: contribute(:,:)
    real(kind=r8_kind), intent(in) :: contribute_l(:,:)
    real(kind=r8_kind), intent(in) :: sum_l(:)
    real(kind=r8_kind), intent(in) :: sum_ocup_l(:)
    real(kind=r8_kind), intent(in) :: sum_cont(:)
    real(kind=r8_kind), intent(in) :: sum_ocup_cont(:)
    real(kind=r8_kind), intent(in) :: sum_total(:)
    real(kind=r8_kind), intent(in) :: sum_total_updn(:,:)
    real(kind=r8_kind), intent(in) :: atomic_charge(:)
    !** End of interface *****************************************

    ! --- declaration of local variables ----------------------
    integer(kind=i4_kind)   :: mu,n,nl
    integer(kind=i4_kind)   :: n_sum, n_tot
    !------------ Executable code ------------------------------------

    n_sum = size(title_sum)
    n_tot = size(title_l)

    if (flag) then
       ! --- contractions ---
       if (population_level > 2 ) then
          if (ssym%n_spin == 2) then
             if (is == 1) then
                write(output_unit,11000)
             else
                write(output_unit,11500)
             endif
          endif
          write(output_unit,15000)i_ir,adjustl(trim(ssym%name(i_ir)))
          write(output_unit,*)
          write(output_unit,50000)
          write(output_unit,47000)
          write(output_unit,*)
          write(output_unit,20000)(title_cont(n),n=1,ssym%dim(i_ir))
          do mu=1,ssym%dim(i_ir)
             write(output_unit,25000)mu,eigval(i_ir)%m(mu,is)*27.211652_r8_kind,&
                  (contribute(mu,n),n=1,ssym%dim(i_ir))
          enddo
          write(output_unit,30000)(sum_cont(mu),mu=1,ssym%dim(i_ir))
          write(output_unit,32000)(sum_ocup_cont(mu),mu=1,ssym%dim(i_ir))
          write(output_unit,*)
          write(output_unit,*)
       end if
       ! --- angular momenta ---
       if (ssym%n_spin == 2) then
          if (is == 1) then
             write(output_unit,11000)
          else
             write(output_unit,11500)
          endif
       endif
       write(output_unit,15000)i_ir,adjustl(trim(ssym%name(i_ir)))
       write(output_unit,*)
       write(output_unit,50000)
       write(output_unit,48000)
       write(output_unit,*)
       write(output_unit,20000)(title_l(n),n=1,n_tot)
       do mu=1,ssym%dim(i_ir)
          write(output_unit,25000)mu,eigval(i_ir)%m(mu,is)*27.211652_r8_kind,&
               (contribute_l(mu,nl),nl=1,n_tot)
       enddo
       write(output_unit,30000)(sum_l(nl),nl=1,n_tot)
       write(output_unit,32000)(sum_ocup_l(nl),nl=1,n_tot)
       write(output_unit,*)
    else
       write(output_unit,90000)
       write(output_unit,92000)
       write(output_unit,*)
       write(output_unit,94000)(title_sum(n),n=1,n_sum)
       write(output_unit,96000)(sum_total(n),n=1,n_sum)
       if(ssym%n_spin==2) then
          write(output_unit,97000)(sum_total_updn(n,1),n=1,n_sum)
          write(output_unit,98000)(sum_total_updn(n,2),n=1,n_sum)
       end if
       write(output_unit,*)
       write(output_unit,*)
       write(output_unit,80000)
       write(output_unit,*)
       write(output_unit,81000)(n,n=1,N_unique_atoms)
       write(output_unit,82000)(atomic_charge(n),n=1,N_unique_atoms)
       if (operations_core_density) then
          write(output_unit,83000)(unique_atoms(n)%Z,n=1,N_unique_atoms)
          write(output_unit,84000)(unique_atoms(n)%Z-atomic_charge(n), &
               n=1,N_unique_atoms)
       else
          write(output_unit,83000)(unique_atoms(n)%Z-unique_atoms(n)%ZC, &
               n=1,N_unique_atoms)
          write(output_unit,84000)(unique_atoms(n)%Z-unique_atoms(n)%ZC- &
               atomic_charge(n),n=1,N_unique_atoms)
       endif
       write(output_unit,*)
       write(output_unit,*)
    endif

    do n=1,N_unique_atoms
       if (operations_core_density) then
          m_charge(n)=unique_atoms(n)%Z-atomic_charge(n)
       else
          m_charge(n)=unique_atoms(n)%Z-unique_atoms(n)%ZC-atomic_charge(n)
       endif
    enddo

    ! formats ------------------------
11000  FORMAT(52X,'MAJORITY SPIN')
11500  FORMAT(52X,'MINORITY SPIN')
15000  FORMAT(43X,'IRREDUCIBLE REPRESENTATION:',I2,1X,A3)
50000  FORMAT(40X,'CONTRIBUTIONS TO THE MOLECULAR ORBITALS')
47000  FORMAT(41X,'MOS IN ROWS / CONTRACTIONS IN COLUMNS')
48000  FORMAT(43X,'MOS IN ROWS / L VALUES IN COLUMNS')
20000  FORMAT('     EIGENVALUE',5X,14(2X,A5,1X),/:&
            (20X,14(2X,A5,1X)))
25000  FORMAT(1X,I3,1X,F10.3,4X,14(1X,F7.3),/:(19X,14(1X,F7.3)))
30000  FORMAT(' SUM           ',4X,11(1X,F7.3),/:(19X,14(1X,F7.3)))
32000  FORMAT(' WITH OCCUP.   ',4X,11(1X,F7.3),/:(19X,14(1X,F7.3)))
90000  FORMAT(40X,'     SUMMARY OVER SPINS AND IRREPS     ')
92000  FORMAT(40X,' CONTRIBUTIONS OF THE ATOMIC ORBITALS  ')
94000  FORMAT('             ',6X,14(2X,A7,1X),/:&
            (19X,14(2X,A7,1X)))
96000  FORMAT(17X,14(1X,F9.3),/:(17X,14(1X,F9.3)))
97000  FORMAT('Minority Spin    ',14(1X,F9.3),/:(17X,14(1X,F9.3)))
98000  FORMAT('Majority Spin    ',14(1X,F9.3),/:(17X,14(1X,F9.3)))
80000  FORMAT(50X,'ATOMIC CHARGES')
81000  FORMAT(' NUCLEUS:          ',14(4X,I2,3X),/:(20X,14(4X,I2,3X)))
82000  FORMAT(' ELECTRONIC CHARGE:',14F9.4,/:(20X,14F9.4))
83000  FORMAT(' NUCLEAR CHARGE:   ',14(1X,F8.4),/:(20X,14(1X,F8.4)))
84000  FORMAT(' ATOMIC CHARGE:    ',14(1X,F8.4),/:(20X,14(1X,F8.4)))

  end subroutine popout
  !*************************************************************


  !*************************************************************

  !*************************************************************
  subroutine sumup_populations(dim, i_ir, is, &
                contribute, pop_col)
    ! Purpose: sumup contributions to a userdefined population analysis
    !          as defined in column
    ! declaration of formal parameters ------------------------
    use orbitalprojection_module
    integer(kind=i4_kind),intent(in) :: dim,  & ! dim of irrep
                                        i_ir, & ! irrep
                                        is      ! spin
    real(kind=r8_kind), intent(in) :: contribute(:,:)
    type(arrmat3), intent(inout) :: pop_col(:) ! really inout?
    !** End of interface *****************************************
    ! --- declaration of local variables ----------------------
    integer(kind=i4_kind) :: i_col, i, i_ind, i_eig, &
         i_cont, ua_cont, l_cont, index
    !------------ Executable code ------------------------------------
    do i_col=1,n_columns
       do i_cont=1,n_con(i_col)
          l_cont=column(i_col)%m(2,i_cont)
          ua_cont=column(i_col)%m(1,i_cont)
          index=orbitalprojection_ob(i_ir,l_cont,ua_cont)
          do i_ind=1,&
               unique_atoms(ua_cont)%symadapt_partner(i_ir,l_cont)%N_independent_fcts
             do i=1,unique_atoms(ua_cont)%l_ob(l_cont)%N_uncontracted_fcts + &
                  unique_atoms(ua_cont)%l_ob(l_cont)%N_contracted_fcts
                do i_eig=1,dim
                   pop_col(i_ir)%m(i_eig,is,i_col)=pop_col(i_ir)%m(i_eig,is,i_col)+&
                        contribute(i_eig,index)
                end do
                index=index+1
             end do
          end do
       end do
    end do
  end subroutine sumup_populations
  !*************************************************************

  !*************************************************************
  subroutine make_dos_plot(dim, i_ir, is, ssym, &
                pop_col, iounit_dos)
    ! Purpose: write a output to a file, which subsequently can be
    !          used for generating dos plots.
    !          For every orbital, that schould be included, the energy and
    !          the intensity have to be written. Normally the intensity is only
    !          the occupation, but it can also be a population, if the user
    !          does only want to include certain contributions
    implicit none
    integer(kind=i4_kind),intent(in) :: dim,  & ! dim of irrep
         i_ir, & ! irrep
         is      ! spin
    type(sym), intent(in) :: ssym
    type(arrmat3), intent(inout) :: pop_col(:) ! really inout?
    integer, intent(in) :: iounit_dos
    !** End of interface *****************************************

    ! --- declaration of local variables -------------------------
    integer(kind=i4_kind) :: i_eig
    real(kind=r8_kind)    :: occo

    !------------ Executable code ------------------------------------
    if(n_columns>0) then
       do i_eig=1,dim
          if(eigval(i_ir)%m(i_eig,is)>eig_min.and.&
               eigval(i_ir)%m(i_eig,is)<eig_max) then
             if(occupation) then
                occo=occ_num(i_ir)%m(i_eig,is)
             else
                occo=ssym%partner(i_ir)
             end if
             write(iounit_dos,'(10F10.5)') eigval(i_ir)%m(i_eig,is), &
                  occo*pop_col(i_ir)%m(i_eig,is,:)
          end if
       end do
    else
       do i_eig=1,dim
          if(eigval(i_ir)%m(i_eig,is)>eig_min.and.&
               eigval(i_ir)%m(i_eig,is)<eig_max) then
             if(occupation) then
                occo=occ_num(i_ir)%m(i_eig,is)
             else
                occo=ssym%partner(i_ir)
             end if
             write(iounit_dos,'(2F10.5)') eigval(i_ir)%m(i_eig,is), occo
          end if
       end do
    end if
  end subroutine make_dos_plot
  !*************************************************************
  !--------------- End of module -------------------------------------
end module population_module
