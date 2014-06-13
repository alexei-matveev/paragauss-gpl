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
module cavity_image_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
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

  use type_module ! type specification parameters
  use iounitadmin_module, only: openget_iounit,returnclose_iounit

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public save_cavity_image, save_tess_image
  public save_points_image, save_cavity_and_points_image

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----
  real(kind=r8_kind) , parameter :: ang2au = 0.529177249_r8_kind

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine save_cavity_image(gepol,N_spheres,xyz_sphere,r_sphere,zero_area)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use filename_module, only: outfile
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: gepol,N_spheres
    real(r8_kind), intent(in)    :: xyz_sphere(:,:),r_sphere(:)
    logical, intent(in)          :: zero_area(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: file_id
    integer(i4_kind) :: i
    !------------ Executable code --------------------------------

!!$    file_id=openget_iounit(trim(outfile('cavity.wrl')), &
!!$            form='formatted', status='unknown')
!!$    write(file_id,'(a16)') '#VRML V1.0 ascii'
!!$    write(file_id,'(a11)') 'Separator {'
!!$
!!$    do i=1,N_spheres
!!$       if(gepol == 93) then
!!$          if(zero_area(i)) cycle
!!$       end if
!!$       write(file_id,'(a11,a16,i5)') 'Separator {','       # sphere ',i
!!$       write(file_id,'(a18)') '        Material {'
!!$       write(file_id,'(a34)') '                ambientColor 0 0 0'
!!$       write(file_id,'(a34)') '                diffuseColor 0 1 1'
!!$       write(file_id,'(a35)') '                specularColor 0 0 0'
!!$       write(file_id,'(a35)') '                emissiveColor 0 0 0'
!!$       write(file_id,'(a27)') '                shininess 0'
!!$       write(file_id,'(a30)') '                transparency 0'
!!$       write(file_id,'(a9)' ) '        }'
!!$       write(file_id,'(a19)') '        Transform {'
!!$       write(file_id,'(a28,3f10.6)') '                translation ',xyz_sphere(i,:)*ang2au
!!$       write(file_id,'(a9)' ) '        }'
!!$       write(file_id,'(a16)') '        Sphere {'
!!$       write(file_id,'(a23,f10.6)') '                radius ',r_sphere(i)*ang2au
!!$       write(file_id,'(a9)' ) '        }'
!!$       write(file_id,'(a1)' ) '}'
!!$    end do
!!$
!!$    write(file_id,'(a1)' ) '}'
!!$
!!$    call returnclose_iounit(file_id)

    file_id=openget_iounit(trim(outfile('cavity1.wrl')), &
            form='formatted', status='unknown')
    write(file_id,'(a15)')              '#VRML V2.0 utf8'
    write(file_id,'(a11)')              'Transform {'
    write(file_id,'(a12)')              '  children ['
    do i=1,N_spheres
       if(gepol == 93) then
          if(zero_area(i)) cycle
          write(file_id,'(a15)')        '    Transform {'
          write(file_id,'(a18,3f10.6)') '      translation ',xyz_sphere(i,:)*ang2au
          write(file_id,'(a16)')        '      children ['
          write(file_id,'(a15)')        '        Shape {'
          write(file_id,'(a27)')        '          geometry Sphere {'
          write(file_id,'(a19,f10.6)')  '            radius ',r_sphere(i)*ang2au
          write(file_id,'(a11)')        '          }'
          write(file_id,'(a33)')        '          appearance Appearance {'
          write(file_id,'(a31)')        '            material Material {'
          write(file_id,'(a34)')        '                diffuseColor 0 1 1'
          write(file_id,'(a35)')        '                specularColor 0 0 0'
          write(file_id,'(a35)')        '                emissiveColor 0 0 0'
          write(file_id,'(a27)')        '                shininess 0'
          write(file_id,'(a30)')        '                transparency 0'
          write(file_id,'(a13)')        '            }'
          write(file_id,'(a11)')        '          }'
          write(file_id,'(a9)' )        '        }'
          write(file_id,'(a7)' )        '      ]'
          write(file_id,'(a5)' )        '    }'
       end if
    end do
    write(file_id,'(a3)' )              '  ]'
    write(file_id,'(a1)' )              '}'

    call returnclose_iounit(file_id)

  end subroutine save_cavity_image
  !*************************************************************

  !*************************************************************
  subroutine save_tess_image(N_total,tes_exp_n,tess_export,vert_bouns)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use filename_module, only: outfile
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: N_total
    integer(i4_kind), intent(in) :: tes_exp_n(:)
    real(r8_kind), intent(in)    :: tess_export(:,:,:)
    integer(i4_kind), intent(in) :: vert_bouns(:,:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: buffer(30)
    integer(i4_kind) :: file_id
    integer(i4_kind) :: i,j,k,n
    !------------ Executable code --------------------------------

    file_id=openget_iounit(trim(outfile('tessera.wrl')), &
            form='formatted', status='unknown')
    write(file_id,'(a16)') '#VRML V1.0 ascii'
    write(file_id,'(a11)') 'Separator {'

    write(file_id,'(a18)') '        Material {'
    write(file_id,'(a34)') '                ambientColor 0 0 0'
    write(file_id,'(a34)') '                diffuseColor 0 1 1'
    write(file_id,'(a35)') '                specularColor 0 0 0'
    write(file_id,'(a35)') '                emissiveColor 0 0 0'
    write(file_id,'(a27)') '                shininess 0'
    write(file_id,'(a30)') '                transparency 0'
    write(file_id,'(a9)' ) '        }'

    do i=1,N_total
       n=tes_exp_n(i)
       write(file_id,'(a40,i5)')       'Separator {                  # tesserea ',i
       write(file_id,'(a20)')          '       Coordinate3 {'
       write(file_id,'(a30)')          '                      point [ '
       do j=1,n
          write(file_id,'(a22,3f10.6,a1)')'                      ',tess_export(i,j,:)*ang2au,','
       end do

       call def_index_order(buffer,n,vert_bouns(i,:,:))

       write(file_id,'(a36)')          '                      ]      # point'
       write(file_id,'(a42)')          '       }                     # Coordinate3'
       write(file_id,'(a23)')          '       IndexedFaceSet {'
       write(file_id,'(a30)')          '                  coordIndex[ '
       write(file_id,'(a22)', &
            ADVANCE='NO')              '                      '
       do k=1,n+1
          write(file_id,'(i20,a1)', ADVANCE='NO') buffer(k)-1, ','
       end do
       write(file_id,'(i3,a1)') -1,','
       write(file_id,'(a41)')          '                  ]          # coordIndex'
       write(file_id,'(a45)')          '        }                    # IndexedFaceSet'
       write(file_id,'(a1)' )          '}'
    end do

    write(file_id,'(a1)' ) '}'

    call returnclose_iounit(file_id)

  end subroutine save_tess_image
  !*************************************************************

  !*************************************************************
  subroutine def_index_order(out_dim,size,in_dim)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(out) :: out_dim(30)
    integer(i4_kind)              :: size
    integer(i4_kind), intent(in)  :: in_dim(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i,j
    !------------ Executable code --------------------------------

    out_dim(1)=1
    out_dim(2)=in_dim(1,1)
    do i=2,size
       j=out_dim(i)
       if(in_dim(j,1) == out_dim(i-1)) then
          out_dim(i+1)=in_dim(j,2)
       else
          out_dim(i+1)=in_dim(j,1)
       end if
    end do

  end subroutine def_index_order
  !*************************************************************

  !*************************************************************
  subroutine save_points_image (N_points, xyz_points)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use filename_module, only: outfile
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (i4_kind), intent (in) :: N_points
    real (r8_kind), intent (in) :: xyz_points(:,:) ! (N_points, 1:3)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: file_id
    integer(i4_kind) :: i
    !------------ Executable code --------------------------------

    file_id=openget_iounit(trim(outfile('points.wrl')), &
            form='formatted', status='unknown')

    write(file_id,'(a15)')           '#VRML V2.0 utf8'
    write(file_id,'(a11)')           'Transform {'
    write(file_id,'(a12)')           '  children ['
    do i=1,N_points
       write(file_id,'(a15)')        '    Transform {'
       write(file_id,'(a18,3f10.6)') '      translation ',xyz_points(i,:)*ang2au
       write(file_id,'(a16)')        '      children ['
       write(file_id,'(a15)')        '        Shape {'
       write(file_id,'(a27)')        '          geometry Sphere {'
       write(file_id,'(a23)')        '            radius 0.05'
       write(file_id,'(a11)')        '          }'
       write(file_id,'(a33)')        '          appearance Appearance {'
       write(file_id,'(a31)')        '            material Material {'
       write(file_id,'(a34)')        '                diffuseColor 1 0 0'
       write(file_id,'(a35)')        '                specularColor 0 0 0'
       write(file_id,'(a35)')        '                emissiveColor 0 0 0'
       write(file_id,'(a27)')        '                shininess 0'
       write(file_id,'(a30)')        '                transparency 0'
       write(file_id,'(a13)')        '            }'
       write(file_id,'(a11)')        '          }'
       write(file_id,'(a9)' )        '        }'
       write(file_id,'(a7)' )        '      ]'
       write(file_id,'(a5)' )        '    }'
    end do
    write(file_id,'(a3)' ) '  ]'
    write(file_id,'(a1)' ) '}'

    call returnclose_iounit(file_id)

  end subroutine save_points_image
  !*************************************************************

  !*************************************************************
  subroutine save_cavity_and_points_image &
       (gepol,N_spheres,xyz_sphere,r_sphere,zero_area,N_points,xyz_points,N_spheres1)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use filename_module, only: outfile
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: gepol,N_spheres
    real(r8_kind), intent(in)    :: xyz_sphere(:,:),r_sphere(:)
    logical, intent(in)          :: zero_area(:)
    integer(i4_kind), intent(in) :: N_points
    real(r8_kind), intent(in)    :: xyz_points(:,:)
    integer(i4_kind), intent(in), optional :: N_spheres1
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: file_id
    integer(i4_kind) :: i
    !------------ Executable code --------------------------------

    file_id=openget_iounit(trim(outfile('cavity_points.wrl')), &
            form='formatted', status='unknown')

    write(file_id,'(a15)')           '#VRML V2.0 utf8'
    write(file_id,'(a11)')           'Transform {'
    write(file_id,'(a12)')           '  children ['
    do i=1,N_spheres
       if(gepol == 93) then
          if(present(N_spheres1) .and. i > N_spheres1) cycle
          if(zero_area(i)) cycle
          write(file_id,'(a15)')        '    Transform {'
          write(file_id,'(a18,3f10.6)') '      translation ',xyz_sphere(i,:)*ang2au
          write(file_id,'(a16)')        '      children ['
          write(file_id,'(a15)')        '        Shape {'
          write(file_id,'(a27)')        '          geometry Sphere {'
          write(file_id,'(a19,f10.6)')  '            radius ',r_sphere(i)*ang2au
          write(file_id,'(a11)')        '          }'
          write(file_id,'(a33)')        '          appearance Appearance {'
          write(file_id,'(a31)')        '            material Material {'
          write(file_id,'(a34)')        '                diffuseColor 0 1 1'
          write(file_id,'(a35)')        '                specularColor 0 0 0'
          write(file_id,'(a35)')        '                emissiveColor 0 0 0'
          write(file_id,'(a27)')        '                shininess 0'
          write(file_id,'(a30)')        '                transparency 0'
          write(file_id,'(a13)')        '            }'
          write(file_id,'(a11)')        '          }'
          write(file_id,'(a9)' )        '        }'
          write(file_id,'(a7)' )        '      ]'
          write(file_id,'(a5)' )        '    }'
       end if
    end do
    do i=1,N_points
       write(file_id,'(a15)')        '    Transform {'
       write(file_id,'(a18,3f10.6)') '      translation ',xyz_points(i,:)*ang2au
       write(file_id,'(a16)')        '      children ['
       write(file_id,'(a15)')        '        Shape {'
       write(file_id,'(a27)')        '          geometry Sphere {'
       write(file_id,'(a23)')        '            radius 0.05'
       write(file_id,'(a11)')        '          }'
       write(file_id,'(a33)')        '          appearance Appearance {'
       write(file_id,'(a31)')        '            material Material {'
       write(file_id,'(a34)')        '                diffuseColor 1 0 0'
       write(file_id,'(a35)')        '                specularColor 0 0 0'
       write(file_id,'(a35)')        '                emissiveColor 0 0 0'
       write(file_id,'(a27)')        '                shininess 0'
       write(file_id,'(a30)')        '                transparency 0'
       write(file_id,'(a13)')        '            }'
       write(file_id,'(a11)')        '          }'
       write(file_id,'(a9)' )        '        }'
       write(file_id,'(a7)' )        '      ]'
       write(file_id,'(a5)' )        '    }'
    end do
    write(file_id,'(a3)' ) '  ]'
    write(file_id,'(a1)' ) '}'

    call returnclose_iounit(file_id)

  end subroutine save_cavity_and_points_image
  !--------------- End of module ----------------------------------
end module cavity_image_module
