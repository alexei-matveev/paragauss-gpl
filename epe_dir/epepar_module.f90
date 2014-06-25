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
module epepar_module
 use type_module
 use iounitadmin_module
 use epecom_module, only: epe_input_dir
 implicit none
        integer(kind=i4_kind), parameter::lim_type_ions=40
        real(kind=r8_kind) :: name_of_type(lim_type_ions)
        character(len=3) cname_of_type(lim_type_ions)
        logical,public       :: fixed_orientation=.false.
        integer(kind=i4_kind)::n_userdefined_atoms =0
        integer(kind=i4_kind)::max_type_ions
        integer(kind=i4_kind)::n_types_central_atoms_3body=0
        real(kind=r8_kind),parameter:: evau=1.0_r8_kind/27.211652_r8_kind
        integer(kind=i4_kind)::io_par,nr
        type ppd
        real(r8_kind) l,u, k,r0, k1, r1, b,r,c,d,qq,rqq,cutoff
        integer n
        end type ppd
        type ppdate
        type(ppd) item
        type(ppdate),pointer:: next
        end type ppdate

        type (ppdate),pointer::par,top,tree
        type(ppd) ::item=ppd(8.01,13.01, 0.0,0.0, 0.0,0.0, 0.0,0.3012,0.0,0.0, 0.0,0.0,2.5, 0)
        public treat_epepar_namelist
        contains

        subroutine treat_epepar_namelist
         namelist/epe_parameters/item,n_types_central_atoms_3body &
               ,max_type_ions,fixed_orientation,n_userdefined_atoms &
               ,name_of_type

        nullify(tree)

        allocate(top)
        top=ppdate(ppd(8.01,8.01, 0.0,0.0, 0.0,0.0, 22764.3,0.149,20.37,0.0, 0.0,0.0,2.5, 1),tree)
        tree=>top

        allocate(top)
        top=ppdate(ppd(8.01,12.01, 0.0,0.0, 0.0,0.0, 1275.2,0.3012,0.0,0.0, 0.0,0.0,2.5, 2),tree)
        tree=>top

        allocate(top)
        top=ppdate(ppd(12.01,12.01, 0.0,0.0, 0.0,0.0, 0.0,0.1,0.0,0.0, 0.0,0.0,2.5, 3),tree)
        tree=>top
        
        !        ae O - bare PP correction 
        allocate(top)
        top=ppdate(ppd(8.02,12.03, 0.0,0.0, 0.0,0.0, 61.64367,0.540419,65.80304,0.0, 0.0,0.0,2.5, 4),tree)
        tree=>top

        !  ae Mg - EPE(Lewis)
        allocate(top)
        top=ppdate(ppd(8.01,12.02, 0.0,0.0, 0.0,0.0, 1501.0,0.3028,30.78,0.0, 0.0,0.0,2.5, 5),tree)
        tree=>top

        ! EPE(Lewis) -  bare PP
        allocate(top)
        top=ppdate(ppd(8.01,12.03, 0.0,0.0, 0.0,0.0, 1505.033102,0.308545,30.816343,0.0, 0.0,0.0,2.5, 6),tree)
        tree=>top

        io_par = openget_iounit( &
        trim(epe_input_dir)//'/epe_simol_parameters', &
                                 form='formatted', status='old')
        item%n=1
        do while(item%n.ne.0)
        item%b=0.0_r8_kind
        item%k=0.0_r8_kind
        item%r0=0.0_r8_kind
        item%k1=0.0_r8_kind
        item%r1=0.0_r8_kind
        item%r=1.0_r8_kind
        item%rqq=0.0_r8_kind
        item%c=0.0_r8_kind
        item%d=0.0_r8_kind
        item%cutoff=2.5_r8_kind
        read (io_par,nml=epe_parameters)
        call get_slsp(item%l,item%u,item%n,nr)
        if(nr.ne.0) then
!!$     print*,nr, ' pair function substituted'
        par%item=item
        par%item%n=nr
        else
        allocate(top)
        top=ppdate(item,tree)
        top%item%n=top%next%item%n+1
        tree=>top
        endif ! nr.ne.0
        enddo ! while(item%n.ne.0)              
        call returnclose_iounit(io_par)
        end  subroutine treat_epepar_namelist

        subroutine      get_slsp(mn,mx,lev,nr)
        use epecom_module, only:  output_epe
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
        write(output_epe,*) 'known parameters assigned to new entry'
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
                subroutine list_epe_par
        use epecom_module, only:  output_epe
        par=>top
            do while(associated(par))
             write(output_epe,*) par%item%n,par%item%l,par%item%u,par%item%cutoff,par%item%rqq
             write(output_epe,*) par%item%b,par%item%r,par%item%c,par%item%d,par%item%qq
             write(output_epe,*) par%item%k,par%item%r0,par%item%k1,par%item%r1
                par=>par%next
            enddo ! while
        end subroutine list_epe_par
        
end module epepar_module
