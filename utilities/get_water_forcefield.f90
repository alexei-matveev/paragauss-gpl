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
program get_water_forcefield

  implicit none

  character*256 :: buffer

  write(6,'(a6)')  ' &task'
  write(6,'(a23)') '   job_task="gradients"'
  write(6,'(a6)')  ' /task'
  write(6,'(a1)')  ' '
  write(6,'(a14)') ' &main_options'
  write(6,'(a25)') '   energy_unit="kcal/mol"'
  write(6,'(a14)') ' /main_options'
  write(6,'(a1)')  ' '
  write(6,'(a18)') ' &nb_lists_options'
  write(6,'(a26)') '   automatic_nb_lists=true'
  write(6,'(a18)') ' /nb_lists_options'

  do 
     buffer=''
     if(eof(5)) exit
     read(5,'(a256)') buffer
     if(get_word(buffer,"&DESCRIPTION")) then
        call read_write_description()
     else if(get_word(buffer,"&TYPES")) then
        call read_write_types()
     else if(get_word(buffer,"&TIP4P")) then
        call read_write_tip4p()
     else if(get_word(buffer,"&LENNARD_JONES")) then
        call read_write_lj()
     end if
  end do

contains

  subroutine upcase(string)

    character(len=*) :: string
    integer*4 :: i,ln,ich

    ln = len(trim(string))
    do i=1,ln
       ich=iachar(string(i:i))
       if(ich>=97 .and. ich<=122) then
          ich=ich-32
          string(i:i)=achar(ich)
       end if
    end do

  end subroutine upcase

  function get_word(buf,word)
    logical :: get_word

    character(len=*) :: buf,word

    call upcase(buf)

    if(index(buf,trim(word)) /= 0) then
       get_word=.true.
    else
       get_word=.false.
    end if
    
  end function get_word

  subroutine read_write_description()

    character*256 :: buf

    write(6,*) '' 
    do
       read(5,'(a256)') buf
       if(get_word(buf,"&END")) then 
          exit
       else
          write(6,*) '#'//trim(buf)
       end if
    end do
    write(6,*) '' 

  end subroutine read_write_description

  subroutine read_write_types()

    character*256 :: buf
    character*6 :: name,internal_name
    real*8 :: charge

    do
       read(5,'(a256)') buf
       if(get_word(buf,"&END")) then 
          exit
       else
          read(buf,*) name,internal_name,charge
          write(6,'(a9)')       ' &species'
          write(6,*)            '  name="'//trim(name)//'"'
          write(6,*)            '  main_name="'//trim(internal_name)//'"'
          write(6,'(a10,f5.2)') '   charge=',charge
          write(6,'(a9)')       ' /species' 
          write(6,*) '' 
       end if
    end do

  end subroutine read_write_types

  subroutine read_write_tip4p()

    character*256 :: buf
    character*10 :: pot_name
    character*6 :: name1,name2,name3
    real*8 :: k0,r0,theta0

    do
       read(5,'(a256)') buf
       if(get_word(buf,"&END")) then 
          exit
       else
          read(buf,*,err=100) name1,name2,pot_name,k0,r0
          write(6,'(a11)')                    ' &potential'
          write(6,*)                          '  pot_name="'//trim(pot_name)//'"'
          write(6,*)                          '  atom_name="'//trim(name1)//'","'//trim(name2)//'"'
          write(6,'(a16,f3.1,a1,f7.4)')       '   ff_parameter=',k0,',',r0
          write(6,'(a11)')                    ' /potential'
          write(6,*)                          '' 
          cycle
100       read(buf,*,err=100) name1,name2,name3,pot_name,k0,theta0
          write(6,'(a11)')                    ' &potential'
          write(6,*)                          '  pot_name="'//trim(pot_name)//'"'
          write(6,*)                          '  atom_name="'//trim(name1)//'","'//trim(name2)//&
               '","'//trim(name3)//'"'
          write(6,'(a16,f3.1,a1,f7.2)')       '   ff_parameter=',k0,',',theta0
          write(6,'(a11)')                    ' /potential'
          write(6,*)                          '' 
       end if
    end do

  end subroutine read_write_tip4p

  subroutine read_write_lj()

    character*256 :: buf
    character*6 :: name(1000)
    character*1 :: sr(1000)
    real*8 :: s_r(1000),eps(1000)
    real*8 :: Ai,Aj,Aij,Ci,Cj,Cij,sigma
    integer*4 :: ir,jr

    ir=0
    do
       read(5,'(a256)') buf
       if(get_word(buf,"&END")) then 
          exit
       else
          ir=ir+1
          read(buf,*) name(ir),s_r(ir),sr(ir),eps(ir)
       end if
    end do

    if(sr(1)=='R') then
       sigma=2.0d0**5*s_r(1)**6
    else
       sigma=s_r(1)**6
    end if
    Ai=4.0d0*eps(1)*sigma**2
    Ci=4.0d0*eps(1)*sigma
    do jr=1,ir
       if(sr(jr)=='R') then
          sigma=2.0d0**5*s_r(jr)**6
       else
          sigma=s_r(jr)**6
       end if
       Aj=4.0d0*eps(jr)*sigma**2
       Cj=4.0d0*eps(jr)*sigma
       Aij=sqrt(Ai*Aj); Cij=sqrt(Ci*Cj)
          write(6,'(a11)')                    ' &potential'
          write(6,'(a17)')                    '   pot_name="L_J"'
          write(6,*)                          '  atom_name="'//trim(name(1))//'","'//trim(name(jr))//'"'
          write(6,'(a16,f9.1,a1,f6.1)')       '   ff_parameter=',Aij,',',Cij
          write(6,'(a11)')                    ' /potential'
          write(6,*)                          '' 
    end do

  end subroutine read_write_lj

end program get_water_forcefield
