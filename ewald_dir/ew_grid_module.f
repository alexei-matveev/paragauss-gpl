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
!C      program triangl_128

      module ew_grid_module
  ! 
  !  Purpose:  Bulding of a triangle  segmented  surface
  !            for for incorporating madelung field effects 
  !            via surface charge destribution
  !
  !  To be used with Ewald master program
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: Zacharenko A.A...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
	integer, parameter :: r8_kind = selected_real_kind(15)
	 real(kind=r8_kind), dimension(3)::rr_x,rr_y,rr_z
	 real(kind=r8_kind)::zero=0.0_r8_kind

	 type vertex
	 real(kind=selected_real_kind(15)), dimension(3)::r
	 end type vertex

	 type triangle
	  type(vertex), dimension(3):: v
	  type(vertex):: xb,xn	! centers of triangele
	 end type triangle

	 type(triangle)::adiv,adv1,adv
         type(triangle),dimension(4)::bdiv,bdv1,bdv
         type(triangle),dimension(4,4)::cdv
	 type(vertex),dimension(8,4)::sec
         type(vertex),dimension(4,8,4)::sec128
	 type div
	 type(vertex)::sec_to_rp
	 real(kind=selected_real_kind(15))::vel
	 real(kind=selected_real_kind(15))::vch
	 real(kind=selected_real_kind(15))::q
	 end type div
	 type (div),dimension(32)::db
	 real*8, dimension(32,32)::db_a
	 integer, dimension(32,3)::db_ipv

         type (div),dimension(128), target ::db_128
         real*8, dimension(128,128)::db_a_128
         integer, dimension(128,3)::db_ipv_128

         real*8:: rc,r_screep=4.0,db_qtot

         type(vertex) :: xc,xm,screep_xyz
	logical:: action_screep=.false.
	logical:: action_screep_c3=.false.
        logical:: action_multipole=.false.
        logical:: action_screep_c4=.false.
	integer ind_scr	! index of first screep point in the mad array
	 



	 contains
	 subroutine scrot_vertex(v,nn_v)
        	 type(vertex),dimension(nn_v),intent(inout)::v
	
	 integer*4, intent(in)::nn_v
	 integer*4 j
	 real*8,dimension(3)::rw
	 rw=0.0
	 do j=1,nn_v
	    rw(1)=rw(1)+dot_product(v(j)%r,rr_x)
	    rw(2)=rw(2)+dot_product(v(j)%r,rr_y)
	    rw(3)=rw(3)+dot_product(v(j)%r,rr_z)
	    v(j)%r(:)=rw(:)
	 enddo !j=1,3  
	 end subroutine scrot_vertex

                        
                        
      subroutine   triangle_centrum(bdv,adv,sec,r_screep,screep_xyz)
	
	 type(triangle),intent(inout)::adv
                       type(triangle),dimension(4),intent(inout)::bdv
                       type(vertex),dimension(8,4),intent(inout)::sec
	real*8,intent(inout):: r_screep
                       type(vertex),intent(inout) :: screep_xyz
	integer::n1,i,j,k
	real*8:: rc,rr2
                       type(vertex) :: xc,xm
	
	
	
	
	 bdv(1)%v(1)=adv%v(1)
        do i=2,3
         xc%r=(adv%v(1)%r+adv%v(i)%r)*0.5d0
         
	 rc= sqrt(dot_product(xc%r,xc%r))
	 xm%r=xc%r/rc
	 bdv(1)%v(i)=xm
	enddo ! i=2,3
!	done first triangle
	xc%r=(adv%v(2)%r+adv%v(3)%r)*0.5d0
	rc= sqrt(dot_product(xc%r,xc%r))
	xm%r=xc%r/rc
	 bdv(2)%v(1)=bdv(1)%v(2)
	 bdv(2)%v(2)=adv%v(2)
	 bdv(2)%v(3)=xm
!	done second triangle
	bdv(3)%v(1)=bdv(2)%v(1)
	bdv(3)%v(2)=bdv(1)%v(3)
	bdv(3)%v(3)=xm
!	done third triangle
	bdv(4)%v(1)=bdv(1)%v(3)
	bdv(4)%v(2)=xm
	bdv(4)%v(3)=adv%v(3)
	do i=1,4
	bdv(i)%xb%r= (bdv(i)%v(1)%r+bdv(i)%v(2)%r+bdv(i)%v(3)%r)/3.d0
	bdv(i)%xn%r= bdv(i)%xb%r/sqrt(dot_product(bdv(i)%xb%r,bdv(i)%xb%r))

!                          write(1,*) 'bdv%xn ', bdv(i)%xn

	enddo ! i=1,4
	
	
!	do j=1,4
!	do i=1,4
!	cdv(i,j)%xb%r= (cdv(i,j)%v(1)%r+cdv(i,j)%v(2)%r+cdv(i,j)%v(3)%r)/3.d0
!	cdv(i,j)%xn%r= cdv(i,j)%xb%r/sqrt(dot_product(cdv(i,j)%xb%r,cdv(i,j)%xb%r))
!
 !                         write(1,*) 'cdv ', cdv(i,j)
!
!	enddo ! i=1,4
!	enddo ! j=1,4
	
	
	

	sec(1,:)=bdv(:)%xn

	sec(2,:)%r(1)=-bdv(:)%xn%r(1)
	sec(2,:)%r(2)=-bdv(:)%xn%r(2)
	sec(2,:)%r(3)= bdv(:)%xn%r(3)

	sec(3,:)%r(1)= bdv(:)%xn%r(2)
	sec(3,:)%r(2)=-bdv(:)%xn%r(1)
	sec(3,:)%r(3)= bdv(:)%xn%r(3)

	sec(4,:)%r(1)=-bdv(:)%xn%r(2)
	sec(4,:)%r(2)= bdv(:)%xn%r(1)
	sec(4,:)%r(3)= bdv(:)%xn%r(3)
 
	sec(5,:)%r(1)=-bdv(:)%xn%r(1)
	sec(5,:)%r(2)=-bdv(:)%xn%r(2)
	sec(5,:)%r(3)=-bdv(:)%xn%r(3)

	sec(6,:)%r(1)=-sec(5,:)%r(1)
	sec(6,:)%r(2)=-sec(5,:)%r(2)
	sec(6,:)%r(3)= sec(5,:)%r(3)

	sec(7,:)%r(1)= sec(5,:)%r(2)
	sec(7,:)%r(2)=-sec(5,:)%r(1)
	sec(7,:)%r(3)= sec(5,:)%r(3)

	sec(8,:)%r(1)=-sec(5,:)%r(2)
	sec(8,:)%r(2)= sec(5,:)%r(1)
	sec(8,:)%r(3)= sec(5,:)%r(3)

	! if(action_screep_c3) then
	!    rr_x = (/ 1./sqrt(2.),-1./sqrt(2.),0. /)
	!    rr_y = (/ 1./sqrt(2.), 1./sqrt(2.),0. /)
	!    rr_z = (/ 0., 0., 1. /)
	!do i=1,8
	!do j=1,4
	!    call scrot_vertex(sec(i,j),1)
	!enddo
	!enddo
	 !   rr_x = (/ 1.0, 0.0, 0.0 /)
	 !   rr_y = (/ zero, adv%xn%r(3), sqrt(1.-adv%xn%r(3)**2) /)
	 !   rr_z = (/ zero, -sqrt(1.-adv%xn%r(3)**2), adv%xn%r(3) /)
	!do i=1,8
	!do j=1,4
	 !   call scrot_vertex(sec(i,j),1)
	!enddo
	!enddo

!	    endif ! action_screep_c3

	!sec%r(1)=sec%r(1)*r_screep+screep_xyz%r(1)
	!sec%r(2)=sec%r(2)*r_screep+screep_xyz%r(2)
	!sec%r(3)=sec%r(3)*r_screep+screep_xyz%r(3)
	!write(6,*) 128
	!write(6,*) '128 point segmentation'
	!do i=1,8
	!do j=1,4
	 !   write(6,*) 'C ', sec(i,j)%r !, ' ',rr2
	!enddo
!	enddo

                   end  subroutine   triangle_centrum




                    subroutine first_second_div

	    use scr_module, only: rr,ds,scrot_ev
	        real*8, dimension(3,3)::rr_save
	        real*8, dimension(3):: ds_save
                      integer::i,j,k


          ! open(1, file='tri_0.txt', form='formatted',status='unknown')

 !          write(1,*) 'screep sphere coordinates in input units'
  !         write(*,*) 'screep sphere coordinates in input units'



! *** first three vertixes
	 adiv%v(1)%r(1)=1.d0
	 adiv%v(1)%r(2)=0.d0
	 adiv%v(1)%r(3)=0.d0
	 adiv%v(2)%r(1)=0.d0
	 adiv%v(2)%r(2)=1.d0
	 adiv%v(2)%r(3)=0.d0
	 adiv%v(3)%r(1)=0.d0
	 adiv%v(3)%r(2)=0.d0
	 adiv%v(3)%r(3)=1.d0


! 	find coordinates of the a triangle centrum
	adiv%xb%r= (adiv%v(1)%r+adiv%v(2)%r+adiv%v(3)%r)/3.d0
	adiv%xn%r= adiv%xb%r/sqrt(dot_product(adiv%xb%r,adiv%xb%r))

        adv%xn%r(3)=adiv%xn%r(3)


      call triangle_centrum(bdv1,adiv,sec,r_screep,screep_xyz)

                 do i=1,4
                 
                   adiv=bdv1(i) 
        call triangle_centrum(bdiv,adiv,sec,r_screep,screep_xyz)
                   sec128(i,:,:)=sec(:,:)
                   cdv(i,:)=bdiv(:)
                      
                 enddo
                 

         if(action_screep_c3) then
            rr_x = (/ 1./sqrt(2.),-1./sqrt(2.),0. /)
            rr_y = (/ 1./sqrt(2.), 1./sqrt(2.),0. /)
            rr_z = (/ 0., 0., 1. /)
        do i=1,4
        do j=1,8
        do k=1,4
            call scrot_vertex(sec128(i,j,k),1)
        enddo
        enddo
        enddo 
            rr_x = (/ 1.0, 0.0, 0.0 /)
            rr_y = (/ zero, adv%xn%r(3), sqrt(1.-adv%xn%r(3)**2) /)
            rr_z = (/ zero, -sqrt(1.-adv%xn%r(3)**2), adv%xn%r(3) /)
        do i=1,4
        do j=1,8
        do k=1,4
            call scrot_vertex(sec128(i,j,k),1)
        enddo
        enddo
        enddo
            endif ! action_screep_c3

        sec128%r(1)=sec128%r(1)*r_screep+screep_xyz%r(1)
        sec128%r(2)=sec128%r(2)*r_screep+screep_xyz%r(2)
        sec128%r(3)=sec128%r(3)*r_screep+screep_xyz%r(3)
        write(6,*) 128
        write(6,*) '128 point segmentation'
        do i=1,4
        do j=1,8
        do k=1,4
            !rr2=sqrt(sec(i,j)%r(1)**2+sec(i,j)%r(2)**2+sec(i,j)%r(3)**2)
            write(6,*) 'C ', sec128(i,j,k)%r !, ' ',rr2
            !write(1,*) 'C sec%xn ', sec(i,j)%r
           !print *, sec(i,j)%r
        enddo
        enddo
        enddo


                 
                     
        end subroutine first_second_div
	 
	 
	 
	 
	  subroutine first_div

	    use scr_module, only: rr,ds,scrot_ev
	        real*8, dimension(3,3)::rr_save
	        real*8, dimension(3):: ds_save





! *** first three vertixes
	 adv%v(1)%r(1)=1.d0
	 adv%v(1)%r(2)=0.d0
	 adv%v(1)%r(3)=0.d0
	 adv%v(2)%r(1)=0.d0
	 adv%v(2)%r(2)=1.d0
	 adv%v(2)%r(3)=0.d0
	 adv%v(3)%r(1)=0.d0
	 adv%v(3)%r(2)=0.d0
	 adv%v(3)%r(3)=1.d0
! 	find coordinates of the a triangle centrum
	adv%xb%r= (adv%v(1)%r+adv%v(2)%r+adv%v(3)%r)/3.d0
	adv%xn%r= adv%xb%r/sqrt(dot_product(adv%xb%r,adv%xb%r))

!	 if(action_screep_c3) then
!	    rr_x = (/ 1./sqrt(2.),-1./sqrt(2.),0. /)
!	    rr_y = (/ 1./sqrt(2.), 1./sqrt(2.),0. /)
!	    rr_z = (/ 0., 0., 1. /)
!	    call scrot_vertex(adv%v(1),1)
!	    call scrot_vertex(adv%v(2),1)
!	    call scrot_vertex(adv%v(3),1)
!	    write(*,*) 'v1 ',adv%v(1)%r(:)
!	    write(*,*) 'v2 ',adv%v(2)%r(:)
!	    write(*,*) 'v3 ',adv%v(3)%r(:)
!	    rr_x = (/ 1.0, 0.0, 0.0 /)
!	    rr_y = (/ zero, adv%xn%r(3), sqrt(1.-adv%xn%r(3)**2) /)
!	    rr_z = (/ zero, -sqrt(1.-adv%xn%r(3)**2), adv%xn%r(3) /)
!            call scrot_vertex(adv%v(1),1)
!            call scrot_vertex(adv%v(2),1)
!            call scrot_vertex(adv%v(3),1)
!            write(1,*) 'v1 ',adv%v(1)%r(:)
!            write(1,*) 'v2 ',adv%v(2)%r(:)
!            write(1,*) 'v3 ',adv%v(3)%r(:)
!	    endif ! action_screep_c3
	 
	 bdv(1)%v(1)=adv%v(1)
        do i=2,3
         xc%r=(adv%v(1)%r+adv%v(i)%r)*0.5d0
	 rc= sqrt(dot_product(xc%r,xc%r))
	 xm%r=xc%r/rc
	 bdv(1)%v(i)=xm
	enddo ! i=2,3
!	done first triangle
	xc%r=(adv%v(2)%r+adv%v(3)%r)*0.5d0
	rc= sqrt(dot_product(xc%r,xc%r))
	xm%r=xc%r/rc
	 bdv(2)%v(1)=bdv(1)%v(2)
	 bdv(2)%v(2)=adv%v(2)
	 bdv(2)%v(3)=xm
!	done second triangle
	bdv(3)%v(1)=bdv(2)%v(1)
	bdv(3)%v(2)=bdv(1)%v(3)
	bdv(3)%v(3)=xm
!	done third triangle
	bdv(4)%v(1)=bdv(1)%v(3)
	bdv(4)%v(2)=xm
	bdv(4)%v(3)=adv%v(3)
	do i=1,4
	bdv(i)%xb%r= (bdv(i)%v(1)%r+bdv(i)%v(2)%r+bdv(i)%v(3)%r)/3.d0
	bdv(i)%xn%r= bdv(i)%xb%r/sqrt(dot_product(bdv(i)%xb%r,bdv(i)%xb%r))
	enddo ! i=1,4

	sec(1,:)=bdv(:)%xn

	sec(2,:)%r(1)=-bdv(:)%xn%r(1)
	sec(2,:)%r(2)=-bdv(:)%xn%r(2)
	sec(2,:)%r(3)= bdv(:)%xn%r(3)

	sec(3,:)%r(1)= bdv(:)%xn%r(2)
	sec(3,:)%r(2)=-bdv(:)%xn%r(1)
	sec(3,:)%r(3)= bdv(:)%xn%r(3)

	sec(4,:)%r(1)=-bdv(:)%xn%r(2)
	sec(4,:)%r(2)= bdv(:)%xn%r(1)
	sec(4,:)%r(3)= bdv(:)%xn%r(3)
 
	sec(5,:)%r(1)=-bdv(:)%xn%r(1)
	sec(5,:)%r(2)=-bdv(:)%xn%r(2)
	sec(5,:)%r(3)=-bdv(:)%xn%r(3)

	sec(6,:)%r(1)=-sec(5,:)%r(1)
	sec(6,:)%r(2)=-sec(5,:)%r(2)
	sec(6,:)%r(3)= sec(5,:)%r(3)

	sec(7,:)%r(1)= sec(5,:)%r(2)
	sec(7,:)%r(2)=-sec(5,:)%r(1)
	sec(7,:)%r(3)= sec(5,:)%r(3)

	sec(8,:)%r(1)=-sec(5,:)%r(2)
	sec(8,:)%r(2)= sec(5,:)%r(1)
	sec(8,:)%r(3)= sec(5,:)%r(3)

	 if(action_screep_c3) then
	    rr_x = (/ 1./sqrt(2.),-1./sqrt(2.),0. /)
	    rr_y = (/ 1./sqrt(2.), 1./sqrt(2.),0. /)
	    rr_z = (/ 0., 0., 1. /)
	do i=1,8
	do j=1,4
	    call scrot_vertex(sec(i,j),1)
	enddo
	enddo
	    rr_x = (/ 1.0, 0.0, 0.0 /)
	    rr_y = (/ zero, adv%xn%r(3), sqrt(1.-adv%xn%r(3)**2) /)
	    rr_z = (/ zero, -sqrt(1.-adv%xn%r(3)**2), adv%xn%r(3) /)
	do i=1,8
	do j=1,4
	    call scrot_vertex(sec(i,j),1)
	enddo
	enddo
	    endif ! action_screep_c3

	sec%r(1)=sec%r(1)*r_screep+screep_xyz%r(1)
	sec%r(2)=sec%r(2)*r_screep+screep_xyz%r(2)
	sec%r(3)=sec%r(3)*r_screep+screep_xyz%r(3)
	write(6,*) 32
	write(6,*) '32 point segmentation'
	do i=1,8
	do j=1,4
	    !write(1,*) 'C ', sec(i,j)%r
	    write(6,*) 'C ', sec(i,j)%r
	   !print *, sec(i,j)%r	
	enddo
	enddo

	         end subroutine first_div
	 
	 
	 
	 
	 
	 
	end module ew_grid_module

       !   end !program triangl_128

