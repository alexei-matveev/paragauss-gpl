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
module slspar_module
        use smlcom
        use coortype_module
	implicit none	
        logical,public       :: fixed_orientation=.false.
        logical,public       :: extended_format=.true.
        integer(kind=i4_kind),public       :: fixed_atom_1=0_i4_kind,&
                                              fixed_atom_2=0_i4_kind,&
                                              fixed_atom_3=0_i4_kind 

	integer(kind=i4_kind)::n_userdefined_atoms =0
	type atom_defined
	integer(kind=i4_kind) ::no
	real(kind=r8_kind)::mass
	end type atom_defined
	type(atom_defined), allocatable:: userdefined_atom(:)
 ! cartesian coordinates -----------------------------------------------
  type(atom_type),allocatable,public       :: atom(:) 
    real(kind=r8_kind),allocatable,public    :: bmat(:,:)
    real(kind=r8_kind),allocatable,public    :: bmat_inv(:,:)
  real(kind=r8_kind),allocatable,public    :: bmat_trans(:,:) 
  real(kind=r8_kind),allocatable,public    ::constraint_mat(:,:)

!!$ 	real(kind=r8_kind)::rcuts	=100.0_r8_kind
        real(kind=r8_kind),parameter:: evau=1.0_r8_kind/27.213_r8_kind
	integer,parameter:: num_known=4 ! number of known pairs
        integer,parameter::io_par=14 ! unit for epe par file
	integer(kind=i4_kind) :: n_types_central_atoms_3body = 0
	integer(kind=i4_kind) :: max_type_ions
	integer(kind=i4_kind),parameter :: lim_type_ions=20
	real(kind=r8_kind), allocatable :: ki(:,:,:),r3b(:)
	real(kind=r8_kind) :: k_i(6),r_3b,atm_nm(6,3)
	real(kind=r8_kind) :: name_of_type(lim_type_ions)
        real(kind=r8_kind) :: theta_0a,theta_a,theta_b,E3b
        real(kind=r8_kind), parameter :: theta_0=109.47_r8_kind
        real(kind=r8_kind), private  :: e1(3),e2(3)
        real(kind=r8_kind),parameter:: pi=3.14159265358979324_r8_kind, &
             zero=0.0_r8_kind
  type epe_shells
     real(kind=r8_kind),dimension(3)::r
     real(kind=r8_kind),dimension(3)::s
     integer(kind=i4_kind)::k
     real(kind=r8_kind)::ant
  end type epe_shells
	type gx_reg
	 real(r8_kind),dimension(3)::r
	 integer::  ieq
	end type gx_reg
  type(epe_shells),dimension(:),allocatable::epe
	  integer(kind=i4_kind):: epe_kl
	integer(kind=i4_kind) :: n_tetrahedrons
	integer(kind=i4_kind), allocatable :: types(:,:)
	integer(kind=i4_kind), allocatable :: n_3b(:),index(:,:,:)
	integer(kind=i4_kind), allocatable :: tetra_atoms(:,:)
        integer(kind=i4_kind):: epe_nucen !number of shells
	logical:: sls_forces=.false.    ! do not calculate SLS forces
	logical:: ex_epe_parameters=.false.
	type ppd	
	real(r8_kind) l,u,b,r,c,d,qq,rqq,cutoff
	integer n
	end type ppd
	type ppdate	
	type(ppd) item
	type(ppdate),pointer:: next
	end type ppdate

	type (ppdate),pointer::par,top,tree
	type(ppd) ::item=ppd(8.01,13.01,0.0,0.3012,0.0,0.0, 0.0,0.0,2.4, 0)

	contains
          subroutine gradients_3_body(kolat)
            use xyzsml, only:gra
           real(kind=r8_kind):: rss1,rss2,scal,cos_th,sin_th,g3b,theta
           integer(kind=i4_kind),intent(in)::kolat  
          integer(kind=i4_kind)::i,l,ia_1,i1,j,ia_2,k1,j1,k,ia_3
! ** contribution to gradients due to 3-body interaction into I region
     if(n_types_central_atoms_3body > 0 ) then
        theta_0a=theta_0*pi/180.0_r8_kind
        fst: do i=1,n_tetrahedrons
           do l=1,5
              if (tetra_atoms(i,l) <= epe_nucen+kolat) goto 1
           enddo
           cycle fst
1          ia_1=tetra_atoms(i,1)
           i1=epe(ia_1)%k
           scnd :do j=2,4
              ia_2=tetra_atoms(i,j)
              if(ia_2.eq.0) exit
              j1=epe(ia_2)%k
              thrd: do k=j+1,5
                 ia_3=tetra_atoms(i,k)
                 if(ia_3.eq.0) exit
                 if(ia_1.ge.epe_kl.and.ia_2.ge.epe_kl &
                      .and.ia_3.ge.epe_kl) cycle thrd
                 k1=epe(ia_3)%k
                 rss1=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                      epe(ia_2)%s-epe(ia_1)%s)
                 rss2=dot_product(epe(ia_3)%s-epe(ia_1)%s, &
                      epe(ia_3)%s-epe(ia_1)%s)
                 scal=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                      epe(ia_3)%s-epe(ia_1)%s)

                 rss1=sqrt(rss1)
                 rss2=sqrt(rss2)
                 cos_th=scal/(rss1*rss2)
                 theta=acos(cos_th)
                 sin_th=sin(theta)
                 g3b=ki(j1,i1,k1)*(theta-theta_0a)
                 e1=(epe(ia_2)%s-epe(ia_1)%s)/rss1
                 e2=(epe(ia_3)%s-epe(ia_1)%s)/rss2

                 if(ia_3 .lt.epe_kl) &
                      gra(:,ia_3)=gra(:,ia_3)+ & 
                      evau*g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
                 if(ia_2.lt.epe_kl ) &
                 gra(:,ia_2)=gra(:,ia_2)+ & 
                      evau*g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)
                 if(ia_1.lt.epe_kl) &
                 gra(:,ia_1)=gra(:,ia_1)+ & 
                      evau*g3b*((rss1-rss2*cos_th)*e1(:) &
                      +(rss2-rss1*cos_th)*e2(:))/ &
                      (rss2*rss1*sin_th)
              enddo thrd!k=j+1,4
           enddo scnd!j=1,3
        enddo fst !i=1,n_tetra_atoms
     endif
   end subroutine gradients_3_body
   
real(kind=r8_kind) function energy_3_body(kolat)
	
! **calc. of the 3-body interaction part of the relaxation energy of
! **lattice
  integer(kind=i4_kind), intent(in):: kolat
  real(kind=r8_kind) :: rss1,rss2,rsr1,rsr2,scal1,scal2

  integer(kind=i4_kind) :: i,j,k,i1,j1,k1,l
  integer(kind=i4_kind) :: ia_1,ia_2,ia_3,k3
	  theta_0a=theta_0*pi/180.0
  e3b=zero
  fst:  do i=1,n_tetrahedrons
!!$     print*,'tetrahedron No', i
     do l=1,5
        if (tetra_atoms(i,l)<=epe_nucen+kolat) goto 1
     enddo
     cycle fst
1    ia_1=tetra_atoms(i,1)
     if(ia_1.eq.0) cycle fst
     i1=epe(ia_1)%k
  
     scnd: do j=2,4
        ia_2=tetra_atoms(i,j)
        if(ia_2.eq.0) cycle
    
        j1=epe(ia_2)%k
        thrd: do k=j+1,5
           ia_3=tetra_atoms(i,k)
           if(ia_3.eq.0) cycle
           if(ia_1 >= epe_kl.and.ia_2 >= epe_kl.and.ia_3 >= epe_kl) cycle thrd
           ! epe_kl number of atoms in cluster+1 , if all 3 atoms are
           ! outside cycle
           k1=epe(ia_3)%k
           rss1=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                epe(ia_2)%s-epe(ia_1)%s)
           rss2=dot_product(epe(ia_3)%s-epe(ia_1)%s, &
                epe(ia_3)%s-epe(ia_1)%s)
           scal1=dot_product(epe(ia_2)%s-epe(ia_1)%s, &
                epe(ia_3)%s-epe(ia_1)%s)
           rss1=sqrt(rss1)
	rss2=sqrt(rss2)
	           theta_a=acos(scal1/(rss1*rss2))
           e3b=e3b+evau*0.5*ki(j1,i1,k1)*((theta_a-theta_0a)**2)
        enddo thrd!k=j+1,5
     enddo scnd!j=2,4
  enddo fst!i=1,n_tetrahedrons

  energy_3_body= e3b
end function energy_3_body
    subroutine building_tet(kolat)
! **procedure looks for tetrahedrally coordinated atoms
! **and their neighbours

      integer(kind=i4_kind) :: i,i1,j,k,k1,l,m,n,ind
      integer(kind=i4_kind) :: status,n2
      integer(kind=i4_kind), intent(in)::kolat
      logical :: exit_cycle
      real(kind=r8_kind) :: dist
      real(kind=r8_kind), dimension(4) :: buf_dist
      integer(kind=i4_kind), dimension(4) :: buf_index
      integer(kind=i4_kind) :: max_ind(1)
      real(kind=r8_kind) :: max_val

      if (allocated(tetra_atoms)) deallocate(tetra_atoms)
      n_tetrahedrons=0

      do i=1,epe_nucen+kolat
	do j=1,n_types_central_atoms_3body
	if(atm_nm(j,2)-epe(i)%ant.lt.0.00001_r8_kind) then
	n_tetrahedrons=n_tetrahedrons+1
	exit
	endif
	enddo
	enddo
        print*,n_tetrahedrons ,' tetrahedrons found',epe_kl,epe_nucen,&
             n_types_central_atoms_3body
      allocate(tetra_atoms(n_tetrahedrons,5),stat=status)
      tetra_atoms=0
      if(status/=0) then
         call error_handler('building_tet:allocation TETRA_ATOMS is failed')
      endif
      i1=0
      lab1:do i=1,epe_nucen+kolat
      k=epe(i)%k
    
	lab2: do j=1,n_types_central_atoms_3body
              if(k == types(j,1)) then
               exit_cycle=.false.
               ind=j
               exit lab2
            else
               exit_cycle=.true.
            endif
         enddo lab2

         if (exit_cycle) cycle lab1

         i1=i1+1
         tetra_atoms(i1,1)=i

         buf_dist=0.0_r8_kind
         buf_index=0_i4_kind
         lab3: do l=1,epe_nucen+kolat
            if (l==i) cycle lab3

            k1=epe(l)%k
            if(k1.eq.0) cycle  lab3
            lab4: do m=2,5
               if (k1 == types(ind,m)) then
                  exit_cycle=.false.
                  exit lab4
               else
                  exit_cycle=.true.
               endif
            enddo lab4

            if (exit_cycle) cycle lab3
            dist=sqrt(dot_product(epe(i)%s-epe(l)%s,epe(i)%s-epe(l)%s))
            if(dist* 0.529177_r8_kind > r3b(ind)) cycle lab3

            lab5: do n=1,4
               if(buf_dist(n) == 0.0_r8_kind) then
                  buf_dist(n)=dist
                  buf_index(n)=l
                  cycle lab3
               endif
            enddo lab5
            max_ind=maxloc(buf_dist)
            max_val=maxval(buf_dist)
            if(dist < max_val) then
               buf_dist(max_ind)=dist
               buf_index(max_ind)=l
            endif
         enddo lab3
         tetra_atoms(i1,2:5)=buf_index
         print*,tetra_atoms(i1,1:5)
      enddo lab1
      print*, ' done building tetrahedron'
	end subroutine building_tet
	
	subroutine init_slsp()
!	initialize values of EPE  pair potential parameters
!       to be called in simol master
	
	nullify(tree)
	allocate(top)
        
	top=ppdate(ppd(8.01,13.01,1275.2,0.3012,0.0,0.0, 0.0,0.0,2.4, 1),tree)
	tree=>top
!  *Al-O*(1)sls_forces

	allocate(top)
	top=ppdate(ppd(8.01,8.01,22764.3,0.149,20.37,0.0, 0.0,0.0,2.4, 2),tree)
	tree=>top
!  *O-O*(2)
	allocate(top)
	top=ppdate(ppd(13.01,13.01,0.0,1.0,0.0,0.0, 0.0,0.0,2.4, 3),tree)
 	tree=>top


	end subroutine init_slsp

	subroutine 	get_slsp(mn,mx,lev,nr)
!	locate pair potential parameters for given mn mx pair
	real(kind=r8_kind), intent(in)::mn,mx
        integer, intent(out)::nr
        integer, intent(in)::lev
        nr=0
	if(lev.gt.0) then
	par=>top	
	 do while(associated(par))
	  if(lev.eq.par%item%n) then
	   item%b=par%item%b
	   item%r=par%item%r
	   item%c=par%item%c
	   item%d=par%item%d
	   item%qq=par%item%qq
           item%rqq=par%item%rqq
           item%cutoff=par%item%cutoff
	print*,'known parameters assigned to new entry'
	   goto 10
	  endif
	par=>par%next
	 enddo ! while
	stop 'make_epe_namelist: entry is not found'
	endif ! lev.gt.0)
   10	par=>top		
	do while(associated(par))
	if(abs(par%item%l-mn).le.0.0001_r8_kind.and. &
     &     abs(par%item%u-mx).le.0.0001_r8_kind) then
	nr=par%item%n
	return
	endif 
	par=>par%next
        enddo ! i=1,num_known
        end subroutine get_slsp

	subroutine list_epe_par
	par=>top
	    do while(associated(par))
	     print*,par%item%n,par%item%l,par%item%u
		par=>par%next
	    enddo ! while
	end subroutine list_epe_par

	subroutine make_epe_namelist
!!$	real(kind=r8_kind)::rcurs
	integer(kind=i4_kind) :: n_3_body,i1,i2,i3,i,j,k
	namelist /tree_body_interaction/ n_3_body,atm_nm,k_i,r_3b,&
	max_type_ions,name_of_type
        character*1:: data_dir='.'
          integer(kind=i4_kind):: nr,next,status
          integer(kind=i4_kind),parameter::io_par=14 

	
          namelist/epe_parameters/item,n_types_central_atoms_3body &
               ,max_type_ions,fixed_orientation,n_userdefined_atoms, &
               name_of_type,extended_format

          namelist /orientation/ fixed_atom_1,fixed_atom_2,fixed_atom_3

 	  inquire (file=adjustl(trim(data_dir))//'/epe_simol_parameters', &
               exist=ex_epe_parameters)
	  if(ex_epe_parameters) then
	print*,' file epe_simol_parameters exist'
             open(io_par, &
                  file=adjustl(trim(data_dir))//'/epe_simol_parameters', &
                  form='formatted',status='old')

	item%n=1
	do while(item%n.ne.0)
           item%qq=0.0_r8_kind
           item%rqq=0.0_r8_kind
           item%c=0.0_r8_kind
           item%d=0.0_r8_kind
           item%cutoff=2.4_r8_kind
	read (io_par,nml=epe_parameters)
	if(n_userdefined_atoms.ne.0) then
	print*, 'userdefined_atoms'
	allocate(userdefined_atom(n_userdefined_atoms))
	do i=1, n_userdefined_atoms
	read(io_par,*)  userdefined_atom(i)%no, userdefined_atom(i)%mass
	print *,  userdefined_atom(i)%no, userdefined_atom(i)%mass  
	enddo
	endif

           if(fixed_orientation) read (io_par,nml=orientation)
           if(fixed_orientation) write (6,nml=orientation)
           print*,"next stop", fixed_atom_1,fixed_atom_2,fixed_atom_3
           if (n_types_central_atoms_3body /= 0) then
              atm_nm=0.0_r8_kind
              allocate(ki(max_type_ions,max_type_ions,max_type_ions),stat=status)
              if(status.ne.0) call error_handler("allocate ki failed")
              allocate(types(n_types_central_atoms_3body,5),r3b(n_types_central_atoms_3body), &
                   n_3b(n_types_central_atoms_3body),index(n_types_central_atoms_3body,6,3), &
                   stat=status)
              if(status.ne.0) call error_handler("allocate types and r3b failed")
              types=0
              ki=0.0_r8_kind
              r3b=0.0_r8_kind
              n_3b=0
              index=0

     do i=1,n_types_central_atoms_3body
        read (io_par,nml=tree_body_interaction)
        write (6,nml=tree_body_interaction)
	n_3b(i)=n_3_body
        
        i1=0
        i2=0
        i3=0
        do j=1,n_3_body
           do k=1,max_type_ions
              if(abs(atm_nm(j,1) - name_of_type(k)).lt.0.00001) i1=k
              if(abs(atm_nm(j,2) - name_of_type(k)).lt.0.00001) i2=k
              if(abs(atm_nm(j,3) - name_of_type(k)).lt.0.00001) i3=k
           enddo
           print*,'index search is done',i1,i2,i3
           if(i1.eq.0.or.i2.eq.0.or.i3.eq.0) &
                call error_handler('index search failed')
           ki(i1,i2,i3)=k_i(j)
           ki(i3,i2,i1)=k_i(j)
           index(i,j,1)=i1
           index(i,j,2)=i2
           index(i,j,3)=i3
           print*,'indexes are filled',i2
           if(j==1) then
              types(i,1)=i2
           endif
           do k=2,5
              if(types(i,k)==i1) exit
              if(types(i,k)==0) then
                 types(i,k)=i1
                 exit
              end if
           enddo
           do k=2,5
              if(types(i,k)==i3) exit
              if(types(i,k)==0) then
                 types(i,k)=i3
                 exit
              end if
              
           enddo
        enddo
        r3b(i)=r_3b
     enddo
        print*,' done nml=tree_body_interaction'
	endif

	print*,'mamelist', item%l,item%u
	call get_slsp(item%l,item%u,item%n,nr)
	if(nr.ne.0) then
	par%item=item
	par%item%n=nr
	print*,nr
	else
	allocate(top)
	top=ppdate(item,tree)
	top%item%n=top%next%item%n+1
	tree=>top
	print*,top%item%n
	endif ! nr.ne.0
	enddo ! while(item%n.ne.0)
!!$	rcurs=rcuts/angsau**2
  if (n_types_central_atoms_3body /= 0) then
     write(6,*) '-------------------------------------------------'
     write(6,*) '              Three-body parameters'
     write(6,*) '-------------------------------------------------'
     write(6,*) '                  ki        alpha       r3b'
     do i=1,n_types_central_atoms_3body
        do j=1,n_3b(i)
           write(6,'(3f6.2,3f12.5)') name_of_type(index(i,j,1)), &
                name_of_type(index(i,j,2)),name_of_type(index(i,j,3)), &
                ki(index(i,j,1),index(i,j,2),index(i,j,3)), &
                109.47, r3b(i)
        enddo
     enddo
     write(6,*) '-------------------------------------------------'
     deallocate(n_3b,index)
  endif
	else
	  open(io_par,file=adjustl(trim(data_dir))//'/epe_simol_parameters', &
               form='formatted',status='new')
	  write(io_par,nml=epe_parameters)
	  close (io_par)
	stop 'file epe_simol_parameters is created'
	endif ! ex_epe_parameters
        
	end subroutine make_epe_namelist
!===============================================================
! Public interface of module
!===============================================================
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
  subroutine alloc_bmat(n_coor,kolat)
    ! Purpose: just a brief wrapper for the few lines necessary to allocate
    !          the B-Matrix 'bmat'
    ! ----------------------------------------------------------------------
    integer(kind=i4_kind),intent(in)   :: n_coor,kolat
    integer(kind=i4_kind) :: alloc_stat


    allocate(bmat(n_coor,3*(kolat)),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" alloc_bmat : allocation (1) failed")
    allocate(bmat_trans(3*(kolat),n_coor),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" alloc_bmat : allocation (2) failed")
    allocate(bmat_inv(3*(kolat),n_coor),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" alloc_bmat : allocation (3) failed")
    bmat=zero
    bmat_trans = zero
    bmat_inv = zero               
    allocate(constraint_mat(3*(kolat),3*(kolat)),&
         STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler &
            (" invert_bmat : allocation (4) failed")
       constraint_mat=zero
  end subroutine alloc_bmat 
end module slspar_module



