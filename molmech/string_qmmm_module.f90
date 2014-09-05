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
module string_qmmm_module
 
  private ! defaults are private use

  public    :: freeunit,gettext,getword

contains

  function freeunit ()

    !     Written by Jay William Ponder (c) 1990
    !
    !     Modified for F90: TK (07.10.1999)
    !
    !     "freeunit" finds an unopened Fortran I/O unit and returns
    !     its numerical value from 1 to 99; the units already assigned
    !     to "input" and "iout" (usually 5 and 6) are skipped since
    !     they have special meaning as the default I/O units
    !

    implicit none
    integer,parameter  :: iout=6,input=5
    integer            :: freeunit
    logical            :: used

    !
    !     try each logical unit until an unopened one is found
    !
    freeunit = 0
    used = .true.
    loop1: do while (used)
       freeunit = freeunit + 1
       if (freeunit.ne.input .and. freeunit.ne.iout) then
          if (freeunit .gt. 99) then
             write (*,*)' FREEUNIT  --  No Available Fortran I/O Units'
             stop 
          end if
          inquire (unit=freeunit,opened=used)
       end if
    end do loop1
 
  end function freeunit

  subroutine gettext (string,text,next)

    !     Written by Jay William Ponder (c) 1990
    !
    !     Modefied for F90: TK (07.10.1999)
    !
    !     ##########################################################
    !     ##                                                      ##
    !     ##  subroutine gettext  --  extract text from a string  ##
    !     ##                                                      ##
    !     ##########################################################
    !
    !     "gettext" searchs an input string for the first string of
    !     non-blank characters; the region from a non-blank character
    !     to the first blank space is returned as "text"; if the
    !     actual text is too long, only the first part is returned
    !
    !     string    input character string to be searched
    !     text      output with the first text string found
    !     next      input with first position of search string;
    !                  output with position following text
    !
    !

    implicit none
    integer           ::  i,j,length,sizetext
    integer           ::  first,last,extent,initial,final
    integer,intent (inout)       ::   next
    character(len=*),intent(in)  ::   string
    character(len=*),intent(out) ::   text
 
    !
    !     get the length of input string and output text
    !
    length = len(string(next:))
    sizetext = len(text)

    !     move through the string one character at a time,
    !     searching for the first non-blank character
    
    first = 1
 
    last = 0
    initial = next
    final = next + length - 1
    loop1: do i = initial, final
       if (string(i:i) .gt. ' ') then
          first = i
          do j = i+1, final
             if (string(j:j) .le. ' ') then
                last = j - 1
                next = j
                exit loop1
             end if
          end do
       end if
    end do loop1

    !     trim the actual text if it is too long to return

    extent = next - first
    final = first + sizetext - 1
    if (extent .gt. sizetext)  last = final

    !     transfer the text into the return string

    j = 0
    do i = first, last
       j = j + 1
       text(j:j) = string(i:i)
    end do
    do i = next, final
       j = j + 1
       text(j:j) = ' '
    end do
    return
  end subroutine gettext

  subroutine getword (string,word,next)


    !     Written by Jay William Ponder (c) 1990
    !
    !     Modefied for F90: TK (07.10.1999)
    !
    !     ################################################################
    !     ##                                                            ##
    !     ##  subroutine getword  --  extract first word from a string  ##
    !     ##                                                            ##
    !     ################################################################
    !
    !
    !     "getword" searchs an input string for the first alphabetic
    !     character (A-Z or a-z); the region from this first character
    !     to the first blank space or comma is returned as a "word"; if
    !     the actual word is too long, only the first part is returned
    !
    !     string    input character string to be searched
    !     word      output with the first word in the string
    !     next      input with first position of search string;
    !                  output with position following word
    !
    !
    
    implicit none
    integer             :: i,j,length,sizetext,next
    integer             :: first,last,extent,initial,final
    character (len=1)   :: letter
    character (len=*)   :: string,word

    !     get the length of input string and output word
    !

    length = len(string(next:))
    sizetext = len(word)

    !     move through the string one character at a time,
    !     searching for the first alphabetic character

    first = 1
    last = 0
    initial = next
    final = next + length - 1
    loop1: do i = initial, final
       letter = string(i:i)
       if ((letter.ge.'A' .and. letter.le.'Z') .or. &
            (letter.ge.'a' .and. letter.le.'z')) then
          first = i
          do j = i+1, final
             if (string(j:j).le.' ' .or. string(j:j).eq.',') then
                last = j - 1
                next = j
                exit loop1
             end if
          end do
       end if
    end do loop1

    !     trim the actual word if it is too long to return

    extent = next - first
    final = first + sizetext - 1
    if (extent .gt. sizetext)  last = final

    !     transfer the word into the return string

    j = 0
    do i = first, last
       j = j + 1
       word(j:j) = string(i:i)
    end do
    do i = next, final
       j = j + 1
       word(j:j) = ' '
    end do

    !     skip over the next character when it is a comma
 
    if (string(next:next) .eq. ',')  next = next + 1

  end subroutine getword

end module string_qmmm_module
