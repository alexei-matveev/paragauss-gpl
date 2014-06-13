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
module strings
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
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

  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  integer(IK),parameter,public ::&
       & stringlength_word       = 16,&
       & stringlength_string     = 64,&
       & stringlength_longstring = 128
!!! word:
!!!wwwwwwwwwwwwwwww
!!! string:
!!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: itoa
  public :: rtoa
  public :: tolower
  public :: toupper
  public :: isalpha
  public :: isdigit
  public :: isalnum
  public :: ispunct
  public :: isspace
  public :: isblank
  public :: spresent
  public :: wpresent
  public :: takeword!(str,pos)
  public :: sindex

!!$       & wboundary,&
!!$       & sindex

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  logical,parameter :: IgnoreCaseDefault = .true.

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function itoa(i) result(a)
     implicit none
     integer(IK), intent(in) :: i
     character(len=16)       :: a
     ! *** end of interface ***

     character(len=16) :: buf

     write(buf,'(I16)') i
     a = trim(adjustl(buf))
   end function itoa

  function rtoa(r,l) result(a)
     implicit none
     real(RK), intent(in)    :: r
     integer(IK), optional, intent(in) :: l
     character(len=20)       :: a
     ! *** end of interface ***

     integer(IK)       :: length
     character(len=20) :: buf
     length = 20
     if ( present(l) ) length = l
     write(buf,'(F20.16)') r
     a = trim(adjustl(buf(:length)))
   end function rtoa

  function wpresent(str,wrd) result(yes)
    ! Word-Present:
    !  returns TRUE if a word
    !  is present in the string:
    !  wpresent("logical function w_present","function") == TRUE
    !  wpresent("logical function w_present","present") == FALSE
    !  for "wrd" other then any word behavior undefined
    implicit none
    character(len=*),intent(in) :: wrd
    character(len=*),intent(in) :: str
    logical                     :: yes !<<< result
    ! *** end of interface ***

    integer(IK) :: n

    yes = .false.

    n = sindex(str,wrd)
    if(n.eq.0)return

    yes = wboundary(n,str)
    if(.not.yes) return

    n = n + LEN(wrd) - 1
    yes = wboundary(n,str)
  end function wpresent

  function spresent(str,wrd) result(yes)
    ! String-Present:
    !  returns TRUE if the string "wrd"
    !  is present in the string "str"
    implicit none
    character(len=*),intent(in) :: wrd
    character(len=*),intent(in) :: str
    logical                     :: yes !<<< result
    ! *** end of interface ***

    yes = (sindex(str,wrd).ne.0)
  end function spresent

  function isalpha(c) result(yes)
    implicit none
    character(len=1),intent(in) :: c
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    yes =  ( LGE(c,"a").and.LLE(c,"z") ) .or.&
         & ( LGE(c,"A").and.LLE(c,"Z") )
  end function isalpha

  function isdigit(c) result(yes)
    implicit none
    character(len=1),intent(in) :: c
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    yes =  ( LGE(c,"0").and.LLE(c,"9") )
  end function isdigit

  function isalnum(c) result(yes)
    implicit none
    character(len=1),intent(in) :: c
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    yes = isalpha(c) .or. isdigit(c)
  end function isalnum

  function isspace(c) result(yes)
    implicit none
    character(len=1),intent(in) :: c
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    yes = isblank(c)
  end function isspace

  function isblank(c) result(yes)
    implicit none
    character(len=1),intent(in) :: c
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    ! space or tab:
    yes = c == ' ' .or. ichar(c) == 9
  end function isblank

  function ispunct(c) result(yes)
    ! punctuation == everything but alphanumeric chars and space
    implicit none
    character(len=1),intent(in) :: c
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    yes = .not. ( isalnum(c) .or. isspace(c) )
  end function ispunct

  function issepar(c) result(yes)
    ! underscore does not count as separator
    implicit none
    character(len=1),intent(in) :: c
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    yes = isspace(c) .or. ispunct(c) .and. c /= '_'
  end function issepar

  function wboundary(n,s) result(yes)
    !
    ! Word Boundary:
    !  returns TRUE if s(n:n) is the first or the last
    !  letter in some word. 
    !  Words defined as assembled of alpha-nums and underscores
    !
    implicit none
    integer(IK),     intent(in) :: n
    character(len=*),intent(in) :: s
    logical                     :: yes !<<<result
    ! *** end of iterface ***

    character   :: c
    integer(IK) :: i,nn

    yes = ( LEN(s) .gt. 0 )
    if(.not.yes) return

    yes = ( n.gt.0 ).and.( n.le.LEN(s) )
    if(.not.yes) return

    c = s(n:n)

    yes = isalnum(c) .or. ( c .eq. "_" )
    if(.not.yes) return

    if( (n.eq.1).or.(n.eq.LEN(s)) ) return

    yes = .false.
    do i = -1, 1, 2
       nn = n + i
       c = s(nn:nn)
       yes = yes.or..not.(isalnum(c).or.c.eq."_")
    enddo
  end function wboundary

  function tolower(in) result(out)
    implicit none
    character(len=*),intent(in) :: in
    character(len=LEN(in))      :: out !<<< result
    ! *** end of interface ***

    integer(IK) :: shift,p
    character   :: c

    shift = ichar('a') - ichar('A')

    ! convert string to lower case letter >>
    do p=1,LEN(in)
       c = in(p:p)
       if ('A' <= c .and. c <= 'Z') then
          out(p:p) = char(ichar(c) + shift)
       else
          out(p:p) = c
       endif
    enddo
  end function tolower

  function toupper(in) result(out)
    implicit none
    character(len=*),intent(in) :: in
    character(len=LEN(in))      :: out !<<< result
    ! *** end of interface ***

    integer(IK) :: shift,p
    character   :: c

    shift = ichar('a') - ichar('A')

    ! convert string to upper case letter >>
    do p=1,LEN(in)
       c = in(p:p)
       if ('a' <= c .and. c <= 'z') then
          out(p:p) = char(ichar(c) - shift)
       else
          out(p:p) = c
       endif
    enddo
  end function toupper

  function sindex(str,wrd,IgnoreCase) result(p)
    implicit none
    character(len=*),intent(in) :: str
    character(len=*),intent(in) :: wrd
    logical,optional,intent(in) :: IgnoreCase
    integer(IK)                 :: p !<<< result
    ! *** end of interface ***

    logical                 :: icase
    character(len=LEN(wrd)) :: cwrd
    character(len=LEN(str)) :: cstr

    p = 0

    if(LEN(wrd).gt.LEN(str)) return

    if(present(IgnoreCase))then
       icase = IgnoreCase
    else
       icase = IgnoreCaseDefault
    endif

    if(icase)then
       cwrd = tolower(wrd)
       cstr = tolower(str)
    else
       cwrd = wrd
       cstr = str
    endif

    p = index(cstr,cwrd)
  end function sindex

  function takeword(str,pos) result(wrd)
     ! returns the first word-token from the string str(pos:)
     ! and increments pos
     implicit none
     character(len=*), intent(in)     :: str
     integer         , intent(inout)  :: pos
     character(len=len(str))          :: wrd
     ! *** end of interface ***

     integer     :: s,e,i

     wrd = repeat(" ",len(wrd))
     if(pos > len(str)) RETURN

     pos = max(pos,1)

     s   = pos ! start at this position

     ! skip word separators:
     do i=pos,len(str)
      if( .not. issepar(str(i:i)) ) exit ! loop
      s = s + 1
     enddo

     e = s-1 ! end at this position

     ! consume next word:
     do i=s,len(str)
      if( issepar(str(i:i)) ) exit ! loop
      e = e + 1
     enddo

     wrd = str(s:e)
     pos = e + 1
   end function takeword

  !--------------- End of module ----------------------------------
end module strings
