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
#include "def.h"
program plot_main
  use readwriteblocked_module
  use type_module
  use machineparameters_module, only : lenblk 
  implicit none
  save
  
  ! central variables :
  ! n_orbitals       : Number of Orbitals to be plotted. This variable is specified
  !                    by the user in the separate input for this plotting interface
  ! n_plots          : If SPLIT_MODE is false, this is equal to n_orbitals
  !                    else it is the total number of plots calculated
  !                    from the contributions of each splitted orbital
  ! n_orb_tot        : Total number of orbitals contained in the tape
  !                    according to its angular momenta and atomar contributions
  ! index_orbs(:)    : List of indices of the orbitals to be plotted
  ! irrep_input(:)   : List of Irreps from the input file
  ! index_input(:)   : List of indices within an irrep from the input file
  ! spin_input(:)    : List of spin-indices from the input file
  ! partner_input(:) : List of the partner indices from the input file
  ! split(:)         : TRUE if there are splitted contributions for this
  !                    orbital
  ! contributions(:) : If an orbital is splitted, this list contains the number of
  !                    contributions (angular momenta and atoms) belonging
  !                    to that orbital
  ! title            : Name of calculation, Title on Plots , character with length=20
  ! orb_actual       : actual index within the data tape 
  ! density          : same as 'density_plot' the only difference being that density
  !                    is used as an input switch which causes the following input
  !                    item to be read in differently. ( If density is set, we can skip
  !                    all input lines describing orbital indices)
  ! density_plot     : TRUE if the data file decribes a density on the grid
  !                    FALSE if this is an orbital plot
  ! phase_swich     : switch to turn dashed lines into full lines and vice versa
  
  logical                           :: density_plot=.false.
  integer(kind=i4_kind)             :: n_grid_x, & ! number of points in x dir
       n_grid_y, & ! number of points in x dir
       n_cont      ! number of contour lines
  integer(kind=i4_kind)             :: n_orb_tot,n_plots, &
       n_densities ! number of densities in a density plot
  integer(kind=i4_kind),allocatable :: index_orbs(:),contributions(:),&
       contribution_list(:)
  logical,allocatable               :: split(:), phase_switch_arr(:)
  real(kind=r8_kind)                ::   step_x,step_y,length_x,length_y
  
  real(kind=r8_kind), allocatable :: x_grid(:), & ! grid points in x direction
                                     y_grid(:) ! grid points in y direction
  real(kind=r8_kind), allocatable :: density_coeff(:)
  ! density_coeff(n_densities) coefficients with wich the densities are multiplied
  ! before summing them up
  real, allocatable :: z_mat(:,:,:)
  integer(kind=i4_kind)           :: orb_actual=0_i4_kind,tape_counter=0_i4_kind,rec_actual=0_i4_kind
       
  ! z_mat(n_grid_x, n_grid_y, 4)
  ! z_mat(:,:,1) function values on grid
  ! z_mat(:,:,2:4) derivatives as calculated from bicub
  real(kind=r8_kind), allocatable :: cont_arr(:) ! value of countour lines

  real(kind=r8_kind) :: x0(3), & ! origin
                        x1(3), & ! x direction
                        x2(3)    ! y direction
  character(len=20)  :: title
  type(readwriteblocked_tapehandle) :: th   ! tape_handle for input file

  ! input variables
  integer(kind=i4_kind)             :: n_orbitals
  logical                           :: plot=.true.,header=.true.,density=.false.
  character(len=3),allocatable      :: irrep_input(:)
  integer(kind=i4_kind),allocatable :: index_input(:)
  integer(kind=i4_kind),allocatable :: spin_input(:)
  integer(kind=i4_kind),allocatable :: partner_input(:)

  integer(kind=i4_kind)             :: n_contour=0
  logical                           :: equidistant=.true., &
                                       linear=.false.,&
                                       exponential=.false.
  real(kind=r8_kind)                :: lower=-1.0_r8_kind,upper=1.0_r8_kind,exp_factor=5.0_r8_kind

  ! variables for  Paper size, units etc.
  real(kind=r8_kind)                :: unit=-99.0_r8_kind, &
                                       sheet_x=190.00_r8_kind,&
                                       sheet_y=270.00_r8_kind
  integer(kind=i4_kind)             :: n_pictures_x=1_i4_kind,&
                                       n_pictures_y=1_i4_kind
  real(kind=r8_kind)                :: margin_x=0.0_r8_kind,&
                                       margin_y=0.0_r8_kind,&
                                       pic_size_x=0.0_r8_kind, &
                                       pic_size_y=0.0_r8_kind
  ! flags for drawing mode, i.e. should several pictures be on one page and how many
  ! etc. ...
  logical                 :: single_sheet=.false., caption=.true.
  integer(kind=i4_kind)   :: n_pics_per_sheet=-99_i4_kind,n_pages=-99_i4_kind,rest_pic=-99_i4_kind, &
       line_thickness=1 ! thickness of contour lines.

  ! colours of background and lines
  real(kind=r8_kind)  :: bac_red=1.0_r8_kind,&
                         bac_green=1.0_r8_kind,&
                         bac_blue=1.0_r8_kind   ! white background
  real(kind=r8_kind)  :: line_red=0.0_r8_kind,&
                         line_green=0.0_r8_kind,&
                         line_blue=0.0_r8_kind, &  ! black lines
                         neg_line_red=0.0_r8_kind,&
                         neg_line_green=0.0_r8_kind,&
                         neg_line_blue=0.0_r8_kind  ! black lines for negative phases
  character(len=10)     :: dash_string=" -F"
  real                  :: dash_width=0.05

  namelist /mode/ header,plot,density
  namelist /orbitals/ n_orbitals
  namelist /contour/ n_contour,equidistant,linear,exponential,lower,upper,exp_factor
  namelist /paper_size/ unit,sheet_x,sheet_y
  namelist /drawing_mode/ caption,single_sheet,n_pics_per_sheet,bac_red,bac_green,bac_blue,&
       line_red, line_green, line_blue, line_thickness, neg_line_red, neg_line_green, neg_line_blue, &
       dash_string, dash_width
  
  integer(kind=i4_kind) :: i,counter,orb_counter,i_pic,pics,i_orb, &
       alloc_stat
  logical               :: phase_switch

  
! init
 lenblk= 8*1024
! end init
  call read_input()
  if (header) then
     if(density) then
        call read_header(1)
     else
        call read_header()
     end if
  end if
  if (.not.plot ) stop

  call set_unit()
  call make_grid()

  counter=0_i4_kind
  orb_counter=0_i4_kind
  print*,' Printing ',n_orbitals,' Orbitals on ',n_pages,' Pages with ',&
       n_pics_per_sheet,' pictures per page ...'
  ! rest_pic is the number of pictures for the last page
  rest_pic=modulo(n_plots,n_pictures_x*n_pictures_y)
  if (rest_pic==0) rest_pic=n_pics_per_sheet

  if(density_plot) then
     allocate(phase_switch_arr(1), stat=alloc_stat)
     if(alloc_stat/=0) &
          call error_handler(&
          'allocation of phase_switch_arr_failed')
     phase_switch_arr(1)=.false.
     call prepare_page(unit,sheet_x,sheet_y)
     call read_orbital(1,1)
     ! note: z_mat(:,:,2) is temporarily used to store the sum
     ! or difference of all densities
     z_mat(:,:,2)=density_coeff(1)*z_mat(:,:,1)
     do i=2,n_densities
        if(header) & 
             call read_header(i)
        call read_orbital(1,1)
        z_mat(:,:,2)=density_coeff(i)*z_mat(:,:,1)+z_mat(:,:,2)
     end do
     z_mat(:,:,1)=z_mat(:,:,2)
     z_mat(:,:,2:4)=0.0
     call plot_orbital(1,1)
  else 
     do i=1,n_pages
        print*,'preparing page for page no. ',i
        call prepare_page(unit,sheet_x,sheet_y)
        pics=n_pics_per_sheet
        if (i==n_pages) pics=rest_pic
        do i_pic=1,pics
           counter=counter+1 ! counter counts all plots including splitted contributions
           if (counter>sum(contributions(1:orb_counter))) then
              orb_counter=orb_counter+1 ! orb_counter counts ORBITALS, 
              ! not splitted contributions
           endif
           if (counter<=n_plots) then
              call read_orbital(orb_counter,counter)
              if (.not.single_sheet) &
                   call new_frame(i_pic)
              call plot_orbital(orb_counter,counter)
           endif
        enddo
     enddo
  end if

  call close()


contains
  
  subroutine make_grid()
    ! Purpose : set up the variables 'x_grid' and 'y_grid'
    ! -----------------------------------------------------
    integer(kind=i4_kind)   :: alloc_stat,i
    real(kind=r8_kind)      :: delta_x,delta_y

    allocate(x_grid(n_grid_x),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("make_grid : allocation (1) failed")
    allocate(y_grid(n_grid_y),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("make_grid : allocation (2) failed")
    
    delta_x=length_x/n_grid_x
    delta_y=length_y/n_grid_y
    
    do i=1,n_grid_x
       x_grid(i) = delta_x*(i-1)
    enddo
    do i=1,n_grid_y
       y_grid(i) = delta_y*(i-1)
    enddo
  end subroutine make_grid

  ! ***************************************************************

  subroutine set_unit()
    ! Purpose: after determining how the n_pics_per_sheet should
    !          should be arragned on the sheet (setting of variables
    !          'n_pictures_x' and 'n_pictures_y'), the units are set
    !          such that each pictures -or the whole sheet in case
    !          the 'single_sheet'-flag is set to TRUE- is filled
    !          completely.
    ! --------------------------------------------------------------
    real(kind=r8_kind) :: u1,u2

    if(single_sheet) then
       pic_size_x=(sheet_x-2.0_r8_kind*margin_x)
       pic_size_y=(sheet_y-2.0_r8_kind*margin_y)
       u1=sheet_x/length_x
       u2=sheet_y/length_y
       unit=min(u1,u2)
       sheet_x=sheet_x/unit
       sheet_y=sheet_y/unit
       n_pages=n_orbitals
       pic_size_x=pic_size_x/unit
       pic_size_y=pic_size_y/unit       
       print*,'set_unit: -----------------------------------------------------'
       write(*,'("unit= ",f8.4)')unit
       write(*,'("Size of Sheet in units ",f8.4," x ",f8.4)')sheet_x,sheet_y
       write(*,'("Size of Sheet in mm    ",f8.4," x ",f8.4)')sheet_x*unit,sheet_y*unit
       write(*,*) " "
       write(*,'("Size of single picture in units ",f8.4," x ",f8.4)')pic_size_x,pic_size_y
       write(*,'("Size of single picture in mm    ",f8.4," x ",f8.4)')pic_size_x*unit,pic_size_y*unit
       print*,'---------------------------------------------------------------'
       print*,' '
       return
    endif

    if (n_pics_per_sheet>1) then
       if (n_pics_per_sheet<=2) then
          n_pictures_x=1
          n_pictures_y=2
       elseif (n_pics_per_sheet <= 4) then
          n_pictures_x=2
          n_pictures_y=2
       elseif ( n_pics_per_sheet <= 6 ) then
          n_pictures_x=2
          n_pictures_y=3
       elseif (n_pics_per_sheet <= 9 ) then
          n_pictures_x=3
          n_pictures_y=3
       elseif (n_pics_per_sheet <= 12 ) then
          n_pictures_x=3
          n_pictures_y=4
       endif
    endif

    if (n_plots>1) then
       margin_x=5.0_r8_kind ! unit is mm
       margin_y=margin_x
       pic_size_x=(sheet_x-n_pictures_x*2.0_r8_kind*margin_x)/n_pictures_x
       pic_size_y=(sheet_y-n_pictures_y*2.0_r8_kind*margin_y)/n_pictures_y
       u1=pic_size_x/length_x
       u2=pic_size_y/length_y
       unit=min(u1,u2)
       sheet_x=n_pictures_x*(pic_size_x+2.0_r8_kind*margin_x)/unit
       sheet_y=n_pictures_y*(pic_size_y+2.0_r8_kind*margin_y)/unit
       pic_size_x=pic_size_x/unit
       pic_size_y=pic_size_y/unit
       margin_x=margin_x/unit
       margin_y=margin_y/unit
    else
       pic_size_x=(sheet_x-n_pictures_x*2.0_r8_kind*margin_x)/n_pictures_x
       pic_size_y=(sheet_y-n_pictures_y*2.0_r8_kind*margin_y)/n_pictures_y
       u1=sheet_x/length_x
       u2=sheet_y/length_y
       unit=min(u1,u2)
       pic_size_x=pic_size_x/unit
       pic_size_y=pic_size_y/unit
       sheet_x=sheet_x/unit
       sheet_y=sheet_y/unit
    endif

    print*,'set_unit: -----------------------------------------------------'
    write(*,'("unit= ",f8.4)')unit
    write(*,'("Size of Sheet in units ",f8.4," x ",f8.4)')sheet_x,sheet_y
    write(*,'("Size of Sheet in mm    ",f8.4," x ",f8.4)')sheet_x*unit,sheet_y*unit
    write(*,*) " "
    write(*,'("Size of single picture in units ",f8.4," x ",f8.4)')pic_size_x,pic_size_y
    write(*,'("Size of single picture in mm    ",f8.4," x ",f8.4)')pic_size_x*unit,pic_size_y*unit
    print*,'---------------------------------------------------------------'
    print*,' '

       

  end subroutine set_unit

  ! ***************************************************************
  subroutine new_frame(pic_no)
    ! Purpose: sets a new window frame for each picture on
    !          a single sheet.
    !
    ! ------------------------------------------------------------
    integer(kind=i4_kind),intent(in)  :: pic_no
    integer(kind=i4_kind)             :: line,column
    real(kind=r8_kind)                :: a,b,c,d

    line=int(ceiling(real(pic_no)/real(n_pictures_x)),kind=i4_kind)
    column=mod(pic_no,n_pictures_x)
    if(column==0) column=n_pictures_x

    a=(column-1)*(2.0*margin_x+pic_size_x)+margin_x
    b=a+pic_size_x
    c=(line-1)*(2.0*margin_y+pic_size_y)+margin_y
    d=c+pic_size_y

    write(*,*)"new_frame: ----------------------------------------------------------------"
    write(*,'(" Lower,left corner of picture ",i2," :  x=",f8.4," y=",f8.4)')pic_no,a*unit,c*unit
    write(*,'("Upper,right corner of picture ",i2," :  x=",f8.4," y=",f8.4)')pic_no,b*unit,d*unit
    write(*,*)"---------------------------------------------------------------------------"
    write(*,*)" "

    call WINDOW(real(a),real(b),real(c),real(d))
    call uscale(real(0.0),real(pic_size_x),real(0.0),real(0.0),real(pic_size_y),real(0.0))
    call lintyp(real(0.01),'-F')
    call frame()

  end subroutine new_frame
    
  ! ***************************************************************

  subroutine read_input()
    ! Purpose: read in the separate input for the plot interface.
    !          The User has to specify:
    !          - number of orbitals to be plotted
    !          - for each orbital to be plotted a NAMELIST 'orbital_list'
    !            consisting of:
    !            Name of Irrep( character(len=3),
    !            Index of orbital within irrep (integer(kind=i4_kind))
    !            Spin of orbital (integer(kind=i4_kind))
    !           
    !         The contents of this namelist is then mapped to the lists
    !         irrep_input(:),index_input(:),spin_input(:),partner_input.
    !
    !         Next, there is the namelist 'contour', which contains
    !         the number of contour lines to be drawn, a switch
    !         which determines if the contours should equidistant
    !         with a range to be specified, or if the user wants
    !         an array 'cont_arr(n_contour)' to be read in which
    !         contains user-specified contour lines.
    ! -------------------------------------------------------------
    character(len=3) :: irrep_name
    integer(kind=i4_kind) :: index_irrep,spin_irrep,partner
    integer(kind=i4_kind),parameter :: io_input=20,max_orbitals=100,&
                                       max_pic=12
    integer(kind=i4_kind) :: status,i
    real(kind=r8_kind)    :: delta,sgn,alpha,arg
    namelist /orbital_list/ irrep_name,index_irrep,spin_irrep, phase_switch, partner
    namelist /densities/ n_densities
    ! default values
    irrep_name = '   '
    index_irrep = 0
    spin_irrep = 0
    DPRINT 'start read input' 
    open(io_input,status='old',form='formatted',&
         file='orb_plot.input',iostat=status)
    if (status/=0) then
       stop 'read_input : could not open orb_plot.input '
    endif

    read(io_input,nml=mode,iostat=status)
    if (status>0) &
         stop 'reading of namelist MODE failed'
    if (.not.header .and. plot ) then
       call error_handler(" For plotting READ_HEADER must be set to TRUE")
    elseif ( header .and. .not.plot ) then
       write(*,*)" read_input: stop reading the input file after namelist MODE "
       write(*,*)"             and proceed to read only the header of plot.dat"
       n_orbitals=0_i4_kind
       return
    endif
    
    if (.not.density) then
       read (io_input,nml=orbitals,iostat=status)
       if (status>0) &
            stop 'reading of namelist ORBITALS failed'
       if (n_orbitals>max_orbitals .or. n_orbitals<0 ) then
          write(*,*)" The specified number of orbitals to be plotted "
          write(*,*)" is either too high or negative. Exit"
          stop
       endif

       allocate(irrep_input(n_orbitals),index_input(n_orbitals),spin_input(n_orbitals),&
            partner_input(n_orbitals), &
            contributions(n_orbitals),split(n_orbitals), &
            phase_switch_arr(n_orbitals), STAT=status)
       if( status/=0 ) &
            stop 'read_input : allocation (1) failed'

       irrep_input = '   '
       index_input = 0
       spin_input = 0
       partner_input=1
       contributions = 1
       split = .false.
       ! Contributions will be set in 'read_header', where this information is read from
       ! the tape, dito for split
       do i=1,n_orbitals
          phase_switch=.false.
          read(io_input,nml=orbital_list,iostat=status)
          if (status>0) then
             write(*,*)"Error:  Could not read in ",i,"th namelist ORBITAL_LIST"
             stop
          endif
          phase_switch_arr(i)=phase_switch
          irrep_input(i)=irrep_name
          index_input(i)=index_irrep
          spin_input(i)=spin_irrep
          partner_input(i)=partner
          print *, 'Reading partner:',partner
       enddo
    else       
       read (io_input,nml=densities,iostat=status)
       if (status>0) &
            stop 'reading of namelist DENSITIES failed'
       if(n_densities<1.or.n_densities>20) then
          write(*,*) 'Wrong value for n_densities'
       end if
       allocate(density_coeff(n_densities), stat=alloc_stat)
       if(alloc_stat/=0) &
            call error_handler('allocating n_densities')
       read(io_input,*) density_coeff
    endif! if (density)

    read(io_input,nml=contour,iostat=status)
    if (status>0) &
         stop 'reading of namelist CONTOUR failed'
    if (n_contour<=0) &
         stop 'Please specify a number of contour lines greater than zero'
    if (equidistant) then ! if contour values are to be computed check these switches
       if ( .not.linear .and. .not.exponential ) &
            stop 'read_input: Please select either LINEAR or EXPONENTIAL'
       if (linear .and. exponential ) &
            stop 'read_input: Please select either LINEAR or EXPONENTIAL'
       if (exp_factor <= 0.0_r8_kind ) &
            stop 'read_input: exponential factor for contour lines must be greater than zero'
    endif! otherwise they do not matter


    allocate(cont_arr(n_contour),STAT=status)
    if( status/=0 ) &    
         stop 'read_input : allocation (2) failed'
    if (equidistant) then
       if (linear) then
          delta=(upper-lower)/n_contour
          do i=1,n_contour
             cont_arr(i)=lower+(i-1)*delta
          enddo
          write(*,'("Contour Values using linear scaling ")')
          do i=1,n_contour
             write(*,'("Line ",i3," Contour value ",es10.3)')i,cont_arr(i)
          enddo
       elseif (exponential) then
          upper=max(abs(upper),abs(lower))
          alpha=upper/(exp(exp_factor*upper)-1.0_r8_kind)
          delta=(upper-lower)/(n_contour-1_i4_kind)
          do i=1,n_contour
             arg=lower+(i-1)*delta
             if (arg>0.0_r8_kind) then
                cont_arr(i)=alpha*(exp(arg*exp_factor) - 1.0_r8_kind)
             else
                cont_arr(i)=-alpha*(exp(abs(arg)*exp_factor) - 1.0_r8_kind)
             endif
          enddo
          write(*,'("Contour Values using function ",f5.2,"*(exp(",f5.2,"*x) -1 )  -----")')&
               alpha,exp_factor
          do i=1,n_contour
             write(*,'("Line ",i3," Contour value ",es10.3)')i,cont_arr(i)
          enddo
       else
          call error_handler("read_input: error when calculating contour lines")
       endif
    else
       do i=1,n_contour
          read(io_input,*)cont_arr(i)
       enddo
    endif

    

    read(io_input,nml=paper_size,iostat=status)
    if (status>0) &
        stop 'reading of namelist PAPER_SIZE failed' 
    if (unit==-99.0_r8_kind) then
       write(*,*)" Take default units"
    else
       if (sheet_x<=0.0 .or. sheet_y<=0.0) call error_handler &
            ("read_input: sheet_x and sheet_y must be greater than zero")
    endif
    
    read(io_input,nml=drawing_mode,iostat=status)
    if (status>0) &
         stop 'reading of namelist DRAWING_MODE failed' 
    if(line_red<0.0_r8_kind.or.line_red>1.0_r8_kind) &
         call error_handler('wrong value for line_red')
    if(line_blue<0.0_r8_kind.or.line_blue>1.0_r8_kind) &
         call error_handler('wrong value for line_blue')
    if(line_green<0.0_r8_kind.or.line_green>1.0_r8_kind) &
         call error_handler('wrong value for line_green')
    if(bac_red<0.0_r8_kind.or.bac_red>1.0_r8_kind) &
         call error_handler('wrong value for bac_red')
    if(bac_blue<0.0_r8_kind.or.bac_blue>1.0_r8_kind) &
         call error_handler('wrong value for bac_blue')
    if(bac_green<0.0_r8_kind.or.bac_green>1.0_r8_kind) &
         call error_handler('wrong value for bac_green')
    if(neg_line_red<0.0_r8_kind.or.neg_line_red>1.0_r8_kind) &
         call error_handler('wrong value for neg_line_red')
    if(neg_line_blue<0.0_r8_kind.or.neg_line_blue>1.0_r8_kind) &
         call error_handler('wrong value for neg_line_blue')
    if(neg_line_green<0.0_r8_kind.or.neg_line_green>1.0_r8_kind) &
         call error_handler('wrong value for neg_line_green')
    if(line_thickness<1) call error_handler('wrong value for level')
    if (single_sheet) then
       n_pics_per_sheet=1_i4_kind
       n_pages=1_i4_kind
    else
       if (n_orbitals>max_pic) then
          if(n_pics_per_sheet==-99) call error_handler &
               ("read_input: Please specify the number of pictures per page n_pics_per_page")
          if(n_pics_per_sheet>max_pic) call error_handler &
               ("read_input: n_pics_per_sheet should not be greater than 12")
       else
          n_pics_per_sheet=min(max_pic,n_pics_per_sheet)
       endif
    endif
       

    close(io_input)
   DPRINT ' end of input read'
  end subroutine read_input
  ! **************************************************************
  subroutine read_header(i_density_dummy)
    ! Purpose: open file with plot data and read some 
    !          general information
    !
    !          The variables
    !          irrep_list,index_list,spin_list,split_list and
    !          contribution_list
    !          which are local to this subroutine are mapped
    !          to the global variables
    !          irrep_input, index_input,spin_input, split and
    !          contributions via the list 'index_orbs' which
    !          contains those indices of orbitals contained in the
    !          tape which are specified in the input for plotting.
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), optional :: i_density_dummy
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: alloc_stat,i,j,i_start, i_density
    real(kind=r8_kind)    :: help_real(1),help_arr(5),help_irrep(3),&
         help_name(20)
    integer(kind=i4_kind),allocatable  :: index_list(:),spin_list(:),partner_list(:)
    logical,allocatable   :: split_list(:)
    character(len=3),allocatable       :: irrep_list(:)
    real(kind=r8_kind)                 :: help_step(4)
    character(len=7)                   :: inter
    !------------ Executable code ------------------------------------
    if(present(i_density_dummy)) then
       i_density=i_density_dummy
    else
       i_density=0
    end if
    ! open file
    help_name=99.0_r8_kind
    if(i_density==0) then
       call readwriteblocked_startread('plot.dat',th)
    else
       write(inter,'(i3)') i_density
       inter=adjustl(inter)
       call readwriteblocked_startread('plot.dat'//trim(inter),th)
    end if
    ! read name of calculation ( not yet implemented )
    call readwriteblocked_read(help_name,th)
    call convert_name(20,help_name,title)
    !
    ! read origin
    call readwriteblocked_read(x0,th)
    ! x-direction
    call readwriteblocked_read(x1,th)
    ! y-direction
    call readwriteblocked_read(x2,th)
    ! step_x step_y length_x length_y
    call readwriteblocked_read(help_step,th)
    step_x=help_step(1)
    step_y=help_step(2)
    length_x=help_step(3)
    length_y=help_step(4)

    ! will this be a density plot or an orbital plot?
    help_real=99.0_r8_kind
    call readwriteblocked_read(help_real,th)
    if (help_real(1)==1.0_r8_kind) then
       density_plot=.true.
       write(*,*)" Header information -------------------------------------"
       write(*,*)" Title of calculation :",title
       write(*,'(" Lower,left corner  ",3(1x,f8.3))')(x0(i),i=1,3)
       write(*,'(" X-direction of Plot",3(1x,f8.3)," Y-direction of Plot",3(1x,f8.3))')&
            (x1(i),i=1,3),(x2(i),i=1,3)
       write(*,*)"     **** DENSITY PLOT ****"
       write(*,*) " -------------------------------------------------------"
       write(*,*)
       n_orb_tot=1_i4_kind
       n_orbitals=1_i4_kind
    else
       if(density) call error_handler &
            ("In the input density is specified, but tape contains orbitals")
       
       density_plot=.false.
    endif

    if (.not.density_plot) then
       ! read total number of orbitals in file
       help_real=99.0_r8_kind
       call readwriteblocked_read(help_real,th)
       n_orb_tot=int(help_real(1),kind=i4_kind)
    endif
    if(i_density<=1) then
       ! let us allocate space for following data
       allocate(irrep_list(n_orb_tot),index_list(n_orb_tot),spin_list(n_orb_tot),&
            partner_list(n_orb_tot),&
            contribution_list(n_orb_tot),split_list(n_orb_tot),STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler &
            ("read_header : allocation (1) failed")
       allocate(index_orbs(n_orbitals),STAT=alloc_stat)
       if(alloc_stat/=0) call error_handler &
            ("read_header : allocation (2) failed")
    end if
    index_orbs=-99_i4_kind
    
    if (.not.density_plot) then
       ! now for every orbital we read irrep spin and number    
       do i=1,n_orb_tot
          call readwriteblocked_read(help_irrep,th)
          call readwriteblocked_read(help_arr,th)
          call convert_name(3,help_irrep,irrep_list(i))
          partner_list(i)=int(help_arr(1),kind=i4_kind)
          spin_list(i)=int(help_arr(2),kind=i4_kind)
          index_list(i)=int(help_arr(3),kind=i4_kind)
          if (help_arr(4)==1.0_r8_kind ) then
             split_list(i)=.true.
          elseif( help_arr(4)==0.0_r8_kind) then
             split_list(i)=.false.
          else
             call error_handler &
                  ("read_header: Error occurred when trying to read split_mode")
          endif
          contribution_list(i)=int(help_arr(5),kind=i4_kind)
       enddo

       ! now identify the indices of the orbitals to be plotted
       i_start=1
       plot: do j=1,n_orbitals   !total number of sets with the SAME set of orbital indices
          orbs :do i=i_start,n_orb_tot !total number of plots sets on the tape (including splitted orbs)
             if (irrep_input(j)/=irrep_list(i) ) cycle orbs
             print *,'i:',i,partner_list(i),partner_input(j)
             if (index_input(j)==index_list(i).and.&
                  spin_input(j)==spin_list(i).and.partner_input(j)==partner_list(i) ) then
                index_orbs(j) = i
                if (split_list(i)) then
                   contributions(j)=contribution_list(i)
                   split(j)=.true. 
                else
                   contributions(j)=1
                   split_list(j)=.false.
                endif
                i_start=i+1
                exit orbs
             endif
          enddo orbs
       enddo plot
    else ! this is a density_plot
       index_orbs=1_i4_kind
       contribution_list=1_i4_kind
       contributions=1_i4_kind
       split=.false.
    endif
       

    if(minval(index_orbs)<0) then
       write(*,*)"error: check the Orbital indices and irreps selected in the input"
       stop
    endif
    if(i_density<=1) then
       ! calculate number of points in x- and y-direction
       n_grid_x=int(floor(length_x/step_x),kind=i4_kind)+1_i4_kind
       n_grid_y=int(floor(length_y/step_y),kind=i4_kind)+1_i4_kind
    end if
    if (i_density>0) then
       n_plots=1_i4_kind
       n_pages=1_i4_kind
    else
       n_plots=0_i4_kind
       do i=1,n_orbitals
          if (split(i)) then
             n_plots=n_plots+contributions(i)  ! if the orbital is split up we have 
          else!                        ! contributions(i) data sets
             n_plots=n_plots+1                 ! if not, only the orbitals itself contributes
          endif
       enddo
       n_pages=int(ceiling(real(n_plots)/real(n_pics_per_sheet)),kind=i4_kind)
    endif



    if (.not.density_plot) then
       ! Output the header information
       write(*,*)" Header information -------------------------------------"
       write(*,*)" Title of calculation :",title
       write(*,'(" Lower,left corner  ",3(1x,f8.3))')(x0(i),i=1,3)
       write(*,'(" X-direction of Plot",3(1x,f8.3)," Y-direction of Plot",3(1x,f8.3))')&
            (x1(i),i=1,3),(x2(i),i=1,3)
       write(*,*)"     **** ORBITAL PLOT ****"
       write(*,*)" Number of orbitals contained in the tape :",n_orb_tot
       write(*,*)" Irrep  Partner Index  Spin  Splitted  Contributions"
       do i=1,n_orb_tot
          if (split_list(i)) then
             write(*,'(4xa3,2x,i4,3x,i4,2x,i3,5x,"True"6x,i3)')irrep_list(i),partner_list(i),index_list(i),spin_list(i),&
                  contribution_list(i)
          else
             write(*,'(4xa3,2x,i4,3x,i4,2x,i3,5x,"False",6x,i3)') irrep_list(i),partner_list(i),index_list(i),spin_list(i),& 
                  1_i4_kind
          endif
       enddo
       write(*,*) " -------------------------------------------------------"
       write(*,*)
    endif
  end subroutine read_header

  subroutine prepare_page(u,w,h)
    ! Purpose: do some initialization stuff connected with 
    !          LRZ-graphic
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind)    :: u,     & ! scaling factor
                             w,     & ! length in x -direction
                             h        ! length in y_direction
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    !------------ Executable code ------------------------------------
    ! open graphic file
    open(67,file='ZZZG67')
    call sheet(real(u),real(w),real(h))
    call setbac(real(bac_red),real(bac_green),real(bac_blue))
    call setrgb(real(0.0),real(0.0),real(0.0))
    call frame()

  end subroutine prepare_page

  subroutine close()
    ! Purpose: do shutdown work
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: alloc_stat 
    !------------ Executable code ------------------------------------
    ! close graphics file
    ! do deallocations
    deallocate(z_mat, cont_arr, &
         stat=alloc_stat)
    if(alloc_stat/=0) call error_handler(&
         'close: deallocation (1) failed')

    if (.not.density) then
       deallocate(irrep_input,index_input,spin_input,partner_input,STAT=alloc_stat)
       if (alloc_stat/=0 ) call error_handler &
            ("close: deallocation (2) failed")
       deallocate(index_orbs,STAT=alloc_stat)
       if (alloc_stat/=0 ) call error_handler &
            ("close: deallocation (3) failed")
       
       deallocate(x_grid,y_grid,STAT=alloc_stat)
       if (alloc_stat/=0 ) call error_handler &
            ("close: deallocation (4) failed")
    endif

    call readwriteblocked_stopread(th)
    close(67)
  end subroutine close

  subroutine read_orbital(i_orbital,i_plot)
    ! Purpose: subroutine reads orbital i_orbital from file
    !          before reading the subroutine skips the correct
    !          amount of data
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind) :: i_orbital ! number of orbital to read
    integer(kind=i4_kind) :: i_plot    ! number of plot, i.e. this
    !                                  ! runs over all splitted contribs.
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: i_orb, orb_to_skip,&
         alloc_stat,i,j
    real(kind=r8_kind),pointer    :: z_pointer(:)
    !------------ Executable code ------------------------------------


    if (.not.allocated(z_mat)) then
       allocate(z_mat(n_grid_x,n_grid_y,4),STAT=alloc_stat)
       if(alloc_stat/=0) &
            stop 'read_orbital : allocation (1) failed'
       z_mat=0.0_r8_kind
    endif
    allocate(z_pointer(n_grid_x*n_grid_y),STAT=alloc_stat)
    if(alloc_stat/=0) &
         stop 'read_orbital : allocation (2) failed'
    z_pointer=0.0_r8_kind


    if(orb_actual==0 .and. index_orbs(i_orbital)==1 ) then
       orb_to_skip=0
       rec_actual=1
    elseif (index_orbs(i_orbital)==index_orbs(orb_actual)) then
       ! we are in the same record consisting of set of split up orbitals
       orb_to_skip=0
       rec_actual=rec_actual+1
    elseif (index_orbs(i_orbital)>orb_actual) then
!!$       print*,'neuer Splitdatensatz'
!!$       print*,'ubound(contribution_list)',ubound(contribution_list,1)
!!$       print*,'orb_actual               ',orb_actual
!!$       print*,'indices :',index_orbs(i_orbital)
!!$       print*,'         ',index_orbs(orb_actual)
!!$       print*,'Anteil 1',sum(contribution_list(1:index_orbs(i_orbital)-1))
!!$       print*,'Anteil 2',sum(contribution_list(1:index_orbs(orb_actual)))
       rec_actual=1
       if(orb_actual==0) then
          orb_to_skip=sum(contribution_list(1:index_orbs(i_orbital)-1))
       else
          orb_to_skip=sum(contribution_list(1:index_orbs(i_orbital)-1))-&
               sum(contribution_list(1:index_orbs(orb_actual)))
       endif
       
    endif
    orb_actual=i_orbital
    
    print*,'read_orbital   i_orb',i_orbital
    print*,'          orb_actual',orb_actual,'  rec_actual ',rec_actual
    print*,'ORB_TO_SKIP  ',orb_to_skip
    
    if(orb_to_skip<0) call error_handler(&
         'read_orbital: orb_to_skip negativ')
    ! first skip records we do not need
    do i_orb=1,orb_to_skip
       call readwriteblocked_read&
            (z_pointer,th)
      tape_counter=tape_counter+1
    end do
    ! now we read actual orbital
    call readwriteblocked_read&
         (z_pointer,th) 
    tape_counter=tape_counter+1
    z_mat(:,:,1)=reshape(real(z_pointer),(/n_grid_x,n_grid_y/))
    deallocate(z_pointer)
  end subroutine read_orbital
    
  subroutine plot_orbital(i_orb,i_plot)
    ! Purpose: subroutine performs actual plotting of orbitals
    !------------ Modules used -----------------------------------    
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind),intent(in)  :: i_orb,i_plot
    integer(kind=i4_kind) :: irb, mue, nue ! parameter for bicub
    integer(kind=i4_kind) :: ier, i_cont
    character(len=10)     :: string
    character(len=27)     :: text_string
    real(kind=r8_kind)    :: aa,ab,ac,ad
    !------------ Executable code ------------------------------------

    if (density_plot) then
       write(text_string,'("Density Plot")')
    else   
       ! This is the caption of the plot
       write(text_string,'(i3,a3,"Pt.",i2,1x,"Sp.",i1)')&
            index_input(i_orb),&
            irrep_input(i_orb),&
            partner_input(i_orb),&
            spin_input(i_orb)
    endif

    aa=0.1*pic_size_x ! This is the left upper corner
    ab=aa+0.4*70.0    ! of the picture
    ac=pic_size_y-1.0 ! 
    ad=ac+0.5         !
    if (caption) then ! blank the area for the caption such that
       !              ! no contour lines are drawn in this area
       call blank(int(i_orb),real(aa),real(ab),real(ac),real(ad))
    endif

    ! first call bicub to get derivatives
    irb=0_i4_kind ! no derivativers are given by the user
    z_mat(:,:,2:4)=0.0_r8_kind
    call bicub(int(n_grid_x),&
         int(n_grid_y),&
         int(n_grid_x),&
         int(n_grid_y),&
         real(x_grid),real(y_grid),&
         z_mat,int(irb),int(ier))
    if(ier==1) call error_handler(&
         'plot_orbital: wrong parameter for bicub')
    ! now we can start plotting
    mue=1.0_r8_kind
    nue=1.0_r8_kind
    write(*,'("Minimal and maximal values for Plot No.",i3,"  ",f11.6,1x,f11.6)')&
         i_plot,minval(z_mat(:,:,1)),maxval(z_mat(:,:,1))
    ! fuer die Beschriftung
    do i_cont=1,n_contour ! loop over contour lines
       if (phase_switch_arr(i_orb)) then
          if(cont_arr(i_cont)<0.0) then
             string='-F'
             call setrgb(real(line_red),real(line_green),real(line_blue))
             call level(int(line_thickness))
          else
             call setrgb(real(neg_line_red),real(neg_line_green),real(neg_line_blue))
             call level(int(line_thickness))
             !string=' -F'
             string=dash_string
          endif
       else
          if(cont_arr(i_cont)<0.0) then
             call setrgb(real(neg_line_red),real(neg_line_green),real(neg_line_blue))
             string=dash_string
             call level(int(line_thickness))
          else
             call setrgb(real(line_red),real(line_green),real(line_blue))
             call level(int(line_thickness))
             string='-F'
          endif
       endif
       call lintyp(real(dash_width),string) 
       call bicont(int(n_grid_x),&
            int(n_grid_y),&
            int(n_grid_x),&
            int(mue),int(n_grid_y),int(nue), &
            real(x_grid),&
            real(y_grid),&
            z_mat,real(cont_arr(i_cont)),int(ier))
!       call licont(int(n_grid_x),&
!            int(n_grid_x),&
!            int(n_grid_y), &
!            real(x_grid),&
!            real(y_grid),&
!            real(z_mat(:,:,1)),real(cont_arr(i_cont)),int(ier))
       if(ier==1) call error_handler&
            ('plot_orbital: wrong parameter for bicont')
       if(ier==2) call error_handler&
            ('plot_orbital: buffer overflow in bicont')
    end do! loop over contour lines
   
    if (caption) then
       call setrgb(real(0.0),real(0.0),real(0.0))
       call moveto(real(aa),real(ac))
       call unblnk(int(i_orb))   ! enable writing in the area of the caption again
       call lintyp(real(0.01),'-F') ! select a line
       call font('CR',real(0.0))    ! and a font
       call size(real(0.8))         ! and a font size
       call text(adjustl(trim(text_string)))
    endif
    ! everything is done
  end subroutine plot_orbital

  subroutine error_handler(message)
    character(len=*) :: message
    print *,message
    stop
  end subroutine error_handler

  subroutine convert_name(n,real_arr,name)
    integer(kind=i4_kind),intent(in)    :: n
    real(kind=r8_kind),intent(in)       :: real_arr(n)
    character(len=n),intent(out)        :: name
    integer                             :: inter,i
    character(len=1)                    :: inter_char(n)
    ! ------------------------------------------
    inter_char=' '
    name = repeat(' ',n)
    do i=1,n
       inter=int(real_arr(i))
       inter_char(i) = char(inter)
    enddo
    do i=1,n
       name(i:i)=inter_char(i)
    enddo
  end subroutine convert_name

end program plot_main
