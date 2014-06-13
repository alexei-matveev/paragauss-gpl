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
#         include "def.h"
          use type_module
          use opt_data_module, only:opt_data_dir, epe_shells   !!!!!
          use math_module
          use iounitadmin_module
          use allocopt_module
        implicit none
        private
!        logical,private       :: fixed_orientation=.false.
!        integer(kind=i4_kind),private       :: fixed_atom_1=0_i4_kind,&
!                                               fixed_atom_2=0_i4_kind,&
!                                               fixed_atom_3=0_i4_kind


        real(kind=r8_kind),parameter, public:: evau=1.0_r8_kind/27.211652_r8_kind
        integer,parameter:: num_known=4 ! number of known pairs
        integer, public::io_par      = -1 ! unit for epe par file
        integer, public::io_par_ffh  = -1 ! unit for ffh parameters
        integer(kind=i4_kind), public :: n_types_central_atoms_3body  = 0
        integer(kind=i4_kind), public :: max_type_ions
        integer(kind=i4_kind),parameter :: lim_type_ions=20

        real(kind=r8_kind), allocatable, public :: ki(:,:,:),theta_0(:,:,:),r3b(:)
        real(kind=r8_kind) :: k_i(6),r_3b,theta_i(6),atm_nm(6,3)
        real(kind=r8_kind), public :: name_of_type(lim_type_ions)
        real(kind=r8_kind), public :: charge_of_type(lim_type_ions)
        real(kind=r8_kind) :: theta_0a,theta_a,theta_b,E3b
        real(kind=r8_kind) :: e1(3),e2(3)
!!!        real(kind=r8_kind),parameter:: zero=0.0_r8_kind  pi=3.14159265358979324_r8_kind

!  type ,public :: epe_shells
!     real(kind=r8_kind),dimension(3)::r
!     real(kind=r8_kind),dimension(3)::s
!     integer(kind=i4_kind)::k
!     real(kind=r8_kind)::ant
!  end type epe_shells
  type ,public :: ewa_centers
    real(kind=r8_kind),dimension(3)::r
    real(kind=r8_kind)::q
  end type ewa_centers

  type (ewa_centers),dimension(:),allocatable,public ::ewa

  type(epe_shells),dimension(:),allocatable,public ::epe

          integer(kind=i4_kind),public:: epe_kl
        integer(kind=i4_kind),public :: n_tetrahedrons
        integer(kind=i4_kind), allocatable,public  :: types(:,:)
        integer(kind=i4_kind), allocatable :: n_3b(:),index(:,:,:)
        integer(kind=i4_kind), allocatable,public :: tetra_atoms(:,:)
        integer(kind=i4_kind),public:: epe_nucen !number of shells

        logical:: sls_forces=.false.    ! do not calculate inter gx file atoms forces
        logical:: ex_epe_parameters=.false.
        logical:: ex_ffh_parameters=.false.

        type ,public :: ppd
         real(r8_kind) l, u,k, k1,r1, r0,b,r,c,d,qq,rqq,cutoff
         integer n
        end type ppd

        type ,public :: ppdate
        type(ppd) item
        type(ppdate),pointer:: next
        end type ppdate

        type (ppdate),pointer,public :: par,top,tree &
         ,top_ffh,tree_ffh

        type(ppd),public  ::item=ppd(8.01,13.01, 0.0,0.0, 0.0,0.0,  0.0,0.3012,0.0,0.0, 0.0,0.0,2.4,  0)

        public get_ffh_parameters,list_ffh_par,get_ffhp,get_slsp,close_slspar,init_slsp_ffh
#ifndef NEW_EPE
        public init_slsp,list_epe_par,make_epe_namelist,gradients_3_body,energy_3_body,building_tet
#endif
        contains

#ifndef NEW_EPE
          subroutine gradients_3_body(kolat)
            use gradient_module, only:grad_cartes
            real(kind=r8_kind) :: deg2rad
           real(kind=r8_kind):: rss1,rss2,scal,cos_th,sin_th,g3b,theta
           integer(kind=i4_kind),intent(in)::kolat
          integer(kind=i4_kind)::i,l,ia_1,i1,j,ia_2,k1,j1,k,ia_3
        integer(kind=i4_kind)::num_3b_links

     if(n_types_central_atoms_3body > 0 ) then
     print*, 'EPE starts from ', epe_kl
        deg2rad=pi/180.0_r8_kind
        num_3b_links=0
        fst: do i=1,n_tetrahedrons
           do l=1,5
              if (tetra_atoms(i,l) <= epe_nucen+kolat) goto 1
           enddo
           cycle fst
1          ia_1=tetra_atoms(i,1)
           if( ia_1<1 .or. ia_1 > size(epe) )then
             ABORT('please fix code, remove logic!')
             cycle fst
           endif
           i1=epe(ia_1)%k
           scnd :do j=2,4
              ia_2=tetra_atoms(i,j)
              if(ia_2.eq.0) exit
              j1=epe(ia_2)%k
              thrd: do k=j+1,5
                 ia_3=tetra_atoms(i,k)
                 if(ia_3.eq.0) exit
                 if(ia_1.ge.epe_kl.and.ia_2.ge.epe_kl.and.ia_3.ge.epe_kl) cycle thrd
                 ! all atoms  outside cluster
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
                 g3b=ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)
                 e1=(epe(ia_2)%s-epe(ia_1)%s)/rss1
                 e2=(epe(ia_3)%s-epe(ia_1)%s)/rss2
                  num_3b_links=num_3b_links+1
!                  print*, '3b(1) tet i atom',num_3b_links,theta/deg2rad,theta_0(j1,i1,k1),ia_1,ia_2,ia_3

!!!           e3b=e3b+evau*0.5*ki(j1,i1,k1)*((theta_a-theta_0(j1,i1,k1)*deg2rad)**2)

                 if(ia_3 .lt.epe_kl)  then
                      grad_cartes(ia_3,:)=grad_cartes(ia_3,:)+ &
                      evau*g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
                  print*, '3b(2) tet i atom',num_3b_links,ia_3,theta/deg2rad,theta_0(j1,i1,k1),j1,i1,k1
                  DPRINT evau*g3b*(cos_th*e2(:)-e1(:))/(rss2*sin_th)
                 endif


                 if(ia_2.lt.epe_kl ) then
                 grad_cartes(ia_2,:)=grad_cartes(ia_2,:)+evau*g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)
!                  print*, '3b(3) tet i atom',num_3b_links,ia_2,theta/deg2rad,theta_0(j1,i1,k1),j1,i1,k1
                  DPRINT evau*g3b*(cos_th*e1(:)-e2(:))/(rss1*sin_th)
                 endif

                 if(ia_1.lt.epe_kl) then
                 grad_cartes(ia_1,:)=grad_cartes(ia_1,:)+evau*g3b*( &
                   (rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/ (rss2*rss1*sin_th)
!                  print*, '3b(4) tet i atom',num_3b_links,ia_1,theta/deg2rad,theta_0(j1,i1,k1),j1,i1,k1
!!                  print*, '3b(4) tet i atom',num_3b_links,ia_1,theta/deg2rad,theta_0(j1,i1,k1),ia_1,ia_2,ia_3
                  DPRINT evau*g3b*((rss1-rss2*cos_th)*e1(:)+(rss2-rss1*cos_th)*e2(:))/(rss2*rss1*sin_th)
                 endif

              enddo thrd!k=j+1,4
           enddo scnd!j=1,3
        enddo fst !i=1,n_tetra_atoms
     endif
   end subroutine gradients_3_body

real(kind=r8_kind) function energy_3_body(kolat)
use opt_data_module, only: io_flepo,crossboundary_3b

! **calc. of the 3-body interaction part of the relaxation energy of
! **lattice
  integer(kind=i4_kind), intent(in):: kolat
  real(kind=r8_kind) :: deg2rad,diff3b
  real(kind=r8_kind) :: rss1,rss2,scal1

  integer(kind=i4_kind) :: i,j,k,i1,j1,k1,l
  integer(kind=i4_kind) :: ia_1,ia_2,ia_3

  deg2rad=pi/180.0_r8_kind
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
           e3b=e3b+evau*0.5*ki(j1,i1,k1)*((theta_a-theta_0(j1,i1,k1)*deg2rad)**2)
        enddo thrd!k=j+1,5
     enddo scnd!j=2,4
  enddo fst!i=1,n_tetrahedrons
  energy_3_body= e3b
  print*, energy_3_body,  ' energy_3_body'
  if(crossboundary_3b.ne.0.0_r8_kind) then
  write(io_flepo,*)  crossboundary_3b*evau,energy_3_body,'epe_3b opt_3b'
  diff3b=abs(crossboundary_3b*evau-energy_3_body)
  ASSERT(diff3b.lt.0.0001)
  endif


end function energy_3_body

    subroutine building_tet(kolat)
! **procedure looks for tetrahedrally coordinated atoms
! **and their neighbours

      integer(kind=i4_kind) :: i,i1,j,k,k1,l,m,n,ind
      integer(kind=i4_kind) :: status
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
       k=epe(i)%k
       do j=1,n_types_central_atoms_3body
!          if(atm_nm(j,2)-epe(i)%ant.lt.0.00001_r8_kind) then
           if(k == types(j,1)) then
              DPRINT i,atm_nm(j,2),epe(i)%ant
              n_tetrahedrons=n_tetrahedrons+1
              exit
           endif
        enddo
     enddo
     DPRINT  'building_tet: tetrahedrons found', n_tetrahedrons

!!$        print*,n_tetrahedrons ,' tetrahedrons found',epe_kl,epe_nucen,&
!!$             n_types_central_atoms_3body

      allocate(tetra_atoms(n_tetrahedrons,5),stat=status)
      ASSERT(status.eq.0)
      tetra_atoms=0

      i1=0
      lab1:do i=1,epe_nucen+kolat
      k=epe(i)%k

        lab2: do j=1,n_types_central_atoms_3body
!              DPRINT 'i epe(i)%k types(1,1)', i,epe(i)%k,types(j,1)
              if(k == types(j,1)) then
               exit_cycle=.false.
               ind=j
!               DPRINT 'atom i is of true type',i,types(j,1)
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
            if(k1.eq.0) then
!             DPRINT 'type of l eq 0', l
            endif
            if(k1.eq.0) cycle  lab3
            lab4: do m=2,5
!              DPRINT 'atom l is of true type',l,k1
               if (k1 == types(ind,m)) then
                  exit_cycle=.false.
                  exit lab4
               else
                  exit_cycle=.true.
               endif
            enddo lab4

            if (exit_cycle) cycle lab3

             dist=sqrt(dot_product(epe(i)%s-epe(l)%s,epe(i)%s-epe(l)%s))
             if(dist*0.529177_r8_kind > r3b(ind)) cycle lab3

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
         DPRINT  'tetra_atoms',i1, tetra_atoms(i1,1:5)
!         print*, epe(tetra_atoms(i1,1))%ant,epe(tetra_atoms(i1,2))%ant, &
!                 epe(tetra_atoms(i1,3))%ant,epe(tetra_atoms(i1,4))%ant, epe(tetra_atoms(i1,5))%ant
      enddo lab1
        end subroutine building_tet
#endif

        subroutine init_slsp_ffh()

!       initialize values of EPE FFH  pair potential parameters
!       to be called in master

        nullify(tree_ffh)
        allocate(top_ffh)

        top_ffh=ppdate(ppd(8.01,13.01, 0.0,0.0, 0.0,0.0,  1275.2,0.3012,0.0,0.0, 0.0,0.0, 2.4, 1),tree_ffh)
        tree_ffh=>top_ffh
!  *Al-O*(1)

        allocate(top_ffh)
        top_ffh=ppdate(ppd(8.01,8.01, 0.0,0.0, 0.0,0.0, 22764.3,0.149,20.37,0.0, 0.0,0.0, 2.4, 2),tree_ffh)
        tree_ffh=>top_ffh
!  *O-O*(2)
        allocate(top_ffh)
        top_ffh=ppdate(ppd(13.01,13.01, 0.0,0.0, 0.0,0.0,  0.0,1.0,0.0,0.0, 0.0,0.0, 2.4, 3),tree_ffh)
        tree_ffh=>top_ffh

        print *,'done init_slsp_ffh'
        end subroutine init_slsp_ffh

        subroutine init_slsp()
!       initialize values of EPE  pair potential parameters
!       to be called in master

        nullify(tree)
        allocate(top)

        top=ppdate(ppd(8.01,13.01, 0.0,0.0, 0.0,0.0,  1275.2,0.3012,0.0,0.0, 0.0,0.0, 2.4, 1),tree)
        tree=>top
!  *Al-O*(1)sls_forces

        allocate(top)
        top=ppdate(ppd(8.01,8.01, 0.0,0.0, 0.0,0.0,  22764.3,0.149,20.37,0.0, 0.0,0.0, 2.4, 2),tree)
        tree=>top
!  *O-O*(2)
        allocate(top)
        top=ppdate(ppd(13.01,13.01, 0.0,0.0, 0.0,0.0,  0.0,1.0,0.0,0.0, 0.0,0.0, 2.4, 3),tree)
        tree=>top


        end subroutine init_slsp

      subroutine get_slsp(mn,mx,lev,nr)
!       locate pair potential parameters for given mn mx pair
        real(kind=r8_kind), intent(in)::mn,mx
        integer, intent(out)::nr
        integer, intent(in)::lev
        nr=0
        if(lev.gt.0) then
        par=>top
         do while(associated(par))
          if(lev.eq.par%item%n) then
           item%k=par%item%k
           item%r0=par%item%r0
           item%k1=par%item%k1
           item%r1=par%item%r1
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
   10   par=>top
        do while(associated(par))
        if(abs(par%item%l-mn).le.0.0001_r8_kind.and. &
     &     abs(par%item%u-mx).le.0.0001_r8_kind) then
        nr=par%item%n
        return
        endif
        par=>par%next
        enddo ! i=1,num_known
        end subroutine get_slsp

        subroutine      get_ffhp(mn,mx,lev,nr)
!       locate pair potential parameters for given mn mx pair
        real(kind=r8_kind), intent(in)::mn,mx
        integer, intent(out)::nr
        integer, intent(in)::lev
        nr=0
        if(lev.gt.0) then
        par=>top_ffh
         do while(associated(par))
          if(lev.eq.par%item%n) then
           item%k=par%item%k
           item%r0=par%item%r0
           item%k1=par%item%k1
           item%r1=par%item%r1
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
        stop 'get_ffhp: entry is not found'
        endif ! lev.gt.0)
   10   par=>top_ffh
        do while(associated(par))

        if(abs(par%item%l-mn).le.0.0001_r8_kind.and. &
     &     abs(par%item%u-mx).le.0.0001_r8_kind) then
        nr=par%item%n
        return
        endif

        par=>par%next
        enddo ! i=1,num_known
        end subroutine get_ffhp

        subroutine list_ffh_par
        par=>top_ffh
            do while(associated(par))
!            print*,par%item%n,par%item%l,par%item%u,par%item%k,par%item%r0
                par=>par%next
            enddo ! while
        end subroutine list_ffh_par

        subroutine list_epe_par
        par=>top
            do while(associated(par))
                par=>par%next
            enddo ! while
       end subroutine list_epe_par

       subroutine make_epe_namelist
        use opt_data_module
        use filename_module, only: inpfile
        integer(kind=i4_kind) :: n_3_body,i1,i2,i3,i,j,k

        namelist /three_body_interaction/ n_3_body,atm_nm,k_i,theta_i,r_3b

          integer(kind=i4_kind):: nr

          namelist/epe_parameters/item,n_types_central_atoms_3body &
               ,max_type_ions,fixed_orientation,name_of_type,charge_of_type

          inquire (file=adjustl(trim(inpfile('epe_simol_parameters'))), &
               exist=ex_epe_parameters)
          if(ex_epe_parameters) then
             DPRINT 'slspar::file epe_simol_parameters exist'
             DPRINT 'slspar:: io_par=openget_iounit(...)'
             io_par=openget_iounit( &
                  file=adjustl(trim(inpfile('epe_simol_parameters'))), &
                  form='formatted',status='old')
             allocopt_stat(18)=0
             DPRINT 'slspar:: io_par=',io_par

        item%n=1
        do while(item%n.ne.0)
           item%k=0.0_r8_kind
           item%r0=0.0_r8_kind
           item%k1=0.0_r8_kind
           item%r1=0.0_r8_kind
           item%qq=0.0_r8_kind
           item%rqq=0.0_r8_kind
           item%c=0.0_r8_kind
           item%d=0.0_r8_kind
           item%cutoff=2.4_r8_kind
           n_types_central_atoms_3body=0

        read (io_par,nml=epe_parameters)

        b3b: if (n_types_central_atoms_3body /= 0) then
              atm_nm=0.0_r8_kind
              allocate(ki(max_type_ions,max_type_ions,max_type_ions), &   ! 1
                   theta_0(max_type_ions,max_type_ions,max_type_ions), &  ! 2
                   types(n_types_central_atoms_3body,5), &                ! 3
                   r3b(n_types_central_atoms_3body), &                    ! 4
                   n_3b(n_types_central_atoms_3body), &                   ! 5
                   index(n_types_central_atoms_3body,6,3), &              ! 6
                                                  stat=allocopt_stat(1))
              if(allocopt_stat(1).ne.0) call error_handler("allocate ki failed")
              allocopt_stat(1)=0 ! ki theta_0 types r3b
              allocopt_stat(2)=0 ! index n_3b
              types=0
              ki=0.0_r8_kind
              theta_0=0.0_r8_kind
              r3b=0.0_r8_kind
              n_3b=0
              index=0

     print*, 'make_epe_namelist read three_body_interaction',n_types_central_atoms_3body
     do i=1,n_types_central_atoms_3body
        read (io_par,nml=three_body_interaction)
        write (6,nml=three_body_interaction)
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
!           print*,'index search is done',i1,i2,i3
           if(i1.eq.0.or.i2.eq.0.or.i3.eq.0) &
                call error_handler('index search failed')
           ki(i1,i2,i3)=k_i(j)
           ki(i3,i2,i1)=k_i(j)
           theta_0(i1,i2,i3)=theta_i(j)
           theta_0(i3,i2,i1)=theta_i(j)
           DPRINT 'theta_0', theta_i(j),i1,i2,i3
           index(i,j,1)=i1
           index(i,j,2)=i2
           index(i,j,3)=i3
!           print*,'indexes are filled',i2
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

        DPRINT ' done nml=tree_body_interaction'
     endif b3b

!       print*,'namelist', item%l,item%u
        call get_slsp(item%l,item%u,item%n,nr)
        if(nr.ne.0) then
!       if(nr.eq.6) print*,item%k
        par%item=item
        par%item%n=nr
!       print*,nr
        else
        allocate(top)
        top=ppdate(item,tree)
        top%item%n=top%next%item%n+1
        tree=>top
!       print*,top%item%n
        endif ! nr.ne.0
        enddo ! while(item%n.ne.0)
  if (n_types_central_atoms_3body /= 0) then
     write(OPT_STDOUT,*) '-------------------------------------------------'
     write(OPT_STDOUT,*) '              Three-body parameters'
     write(OPT_STDOUT,*) '-------------------------------------------------'
     write(OPT_STDOUT,*) '                  ki        theta       r3b'
     do i=1,n_types_central_atoms_3body
        do j=1,n_3b(i)
           write(OPT_STDOUT,'(3f6.2,3f12.5)') name_of_type(index(i,j,1)), &
                name_of_type(index(i,j,2)),name_of_type(index(i,j,3)), &
                ki(index(i,j,1),index(i,j,2),index(i,j,3)), &
                theta_0(index(i,j,1),index(i,j,2),index(i,j,3)), r3b(i)
        enddo
     enddo
     write(OPT_STDOUT,*) '-------------------------------------------------'
     deallocate(n_3b,index,stat=allocopt_stat(2))
     ASSERT(allocopt_stat(2).eq.0)
  endif
      else
           DPRINT 'slspar:: io_par=openget_iounit(..2..)'
           io_par=openget_iounit(file=adjustl(trim(inpfile('epe_simol_parameters'))), &
                form='formatted',status='new')
           allocopt_stat(19)=0
           DPRINT 'slspar:: io_par=',io_par
          write(io_par,nml=epe_parameters)
          call returnclose_iounit(io_par)
          allocopt_stat(19)=1
          stop 'file epe_simol_parameters is created'
       endif ! ex_epe_parameters
       DPRINT 'done make_epe_namelist'

    end subroutine make_epe_namelist







        subroutine get_ffh_parameters

        use opt_data_module
        use filename_module, only: inpfile

        integer(kind=i4_kind) :: n_3_body,i1,i2,i3,i,j,k
        integer(kind=i4_kind):: nr

        namelist /three_body_interaction/ n_3_body,atm_nm,k_i,theta_i,r_3b

          namelist/epe_parameters/item,n_types_central_atoms_3body &
               ,max_type_ions,fixed_orientation,name_of_type,charge_of_type

          inquire (file=adjustl(trim(opt_data_dir))//'/ffh_parameters', &
               exist=ex_ffh_parameters)

     if(ex_ffh_parameters) then
           print*,' file ffh_parameters exist'
!             open(io_par_ffh, &
             io_par_ffh=openget_iounit( &
                  file=adjustl(trim(inpfile('ffh_parameters'))), &
                  form='formatted',status='old')
           allocopt_stat(20)=0

        item%n=1

 !      print*,'n_types_central_atoms_3body',n_types_central_atoms_3body
        do while(item%n.ne.0)

           item%k=0.0_r8_kind
           item%r0=0.0_r8_kind
           item%k1=0.0_r8_kind
           item%r1=0.0_r8_kind
           item%qq=0.0_r8_kind
           item%rqq=0.0_r8_kind
           item%c=0.0_r8_kind
           item%d=0.0_r8_kind
           item%cutoff=2.4_r8_kind

        read (io_par_ffh,nml=epe_parameters)
!       print* , ' nml=epe_parameters read',item%l,item%u

!!! to be corrected after pair potentials go
      b3if: if (n_types_central_atoms_3body /= 0) then
              atm_nm=0.0_r8_kind
              allocate(ki(max_type_ions,max_type_ions,max_type_ions), &
                   theta_0(max_type_ions,max_type_ions,max_type_ions), &
                   types(n_types_central_atoms_3body,5), &
                   r3b(n_types_central_atoms_3body), &
                   n_3b(n_types_central_atoms_3body), &
                   index(n_types_central_atoms_3body,6,3), &
                                       stat=allocopt_stat(1))
              if(allocopt_stat(1).ne.0) call error_handler("allocate ki_ffh failed")
                   allocopt_stat(1)=1 ! ki theta_0 types r3b
                   allocopt_stat(2)=1 ! index n_3b
              types=0
              ki=0.0_r8_kind
              theta_0=0.0_r8_kind
              r3b=0.0_r8_kind
              n_3b=0
              index=0

 !       print*, 'get_ffh_parameters read three_body_interaction unit',io_par_ffh,n_types_central_atoms_3body
     b3: do i=1,n_types_central_atoms_3body
        read (io_par_ffh,nml=three_body_interaction)
        write (6,nml=three_body_interaction)
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
!           print*,'index search is done',i1,i2,i3
           if(i1.eq.0.or.i2.eq.0.or.i3.eq.0) &
                call error_handler('index search failed')
           ki(i1,i2,i3)=k_i(j)
           ki(i3,i2,i1)=k_i(j)
           theta_0(i1,i2,i3)=theta_i(j)
           theta_0(i3,i2,i1)=theta_i(j)
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
     enddo b3
 !      print*,' done nml=tree_body_interaction'
     endif b3if

!       print*,'namelist', item%l,item%u, item%k
 !      print*,'get_ffhp'
        call get_ffhp(item%l,item%u,item%n,nr)
 !      print*,'done', item%n,'item%n'
!       if(nr.eq.6) print*,item%k, item%n

        if(nr.ne.0) then

        par%item=item
        par%item%n=nr
!       print*,nr, 'nr'

        else

        allocate(top_ffh)

        top_ffh=ppdate(item,tree_ffh)
        top_ffh%item%n=top_ffh%next%item%n+1
        tree_ffh=>top_ffh

!       print*,top_ffh%item%n, 'allocated',item%l,item%u,item%k

        endif ! nr.ne.0

        enddo ! while(item%n.ne.0)

  if (n_types_central_atoms_3body /= 0) then
     write(6,*) '-------------------------------------------------'
     write(6,*) '              Three-body parameters'
     write(6,*) '-------------------------------------------------'
     write(6,*) '                  ki        theta       r3b'
     do i=1,n_types_central_atoms_3body
        do j=1,n_3b(i)
           write(6,'(3f6.2,3f12.5)') name_of_type(index(i,j,1)), &
                name_of_type(index(i,j,2)),name_of_type(index(i,j,3)), &
                ki(index(i,j,1),index(i,j,2),index(i,j,3)), &
                theta_0(index(i,j,1),index(i,j,2),index(i,j,3)), r3b(i)
        enddo
     enddo
     write(6,*) '-------------------------------------------------'
     deallocate(n_3b,index,stat=allocopt_stat(2))
     ASSERT(allocopt_stat(2).eq.0)
  endif

!         call returnclose_iounit (io_par_ffh)
!          allocopt_stat(20)=1
  else
           DPRINT 'slspar:: io_par=openget_iounit(..3..)'
           io_par=openget_iounit(file=adjustl(trim(inpfile('ffh_parameters'))), &
                form='formatted',status='new')
           allocopt_stat(21)=0

           DPRINT 'slspar:: io_par=',io_par
          write(io_par,nml=epe_parameters)
          call returnclose_iounit (io_par_ffh)
          allocopt_stat(20)=1
        stop 'file ffh_parameters is created'

  endif ! ex_ffh_parameters

        end subroutine get_ffh_parameters

        subroutine close_slspar()
         if(allocated(theta_0)) then
         deallocate(theta_0,ki,types,r3b, stat=allocopt_stat(1))
         ASSERT(allocopt_stat(1).eq.0)
         allocopt_stat(1)=1
         endif
        end subroutine close_slspar

        end module slspar_module



