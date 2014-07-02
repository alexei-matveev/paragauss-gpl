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
subroutine gengrp(group,igroup,iorder,elem,namele,ncsco,fcsco)
  !-------------------------------------------------------------------
  !
  !  Purpose:
  !     identifies the group by number in the canonical ordering used in
  !       the tables of Altmann and Herzig and generates the quaternionic
  !       representation of the group elements and their name
  !
  !  author:    N. Roesch, Theoretische Chemie, TU Muenchen
  !  version:   07.05.95
  !
!== Interrupt of public interface of module =====================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

  !------------ Modules used -----------------------------------------
  use type_module, only: i4_kind,r8_kind
  implicit none
!== Interrupt end of public interface of module =================
  !------------ Declaration of formal parameters ---------------------
  character(len=4), intent(in)           :: group
  !       group   name of the group
  integer(kind=i4_kind),intent(in)  :: iorder,igroup
  !       iorder  order of the group
  real(kind=r8_kind),intent(out)    :: elem(5,iorder)
  !       elem    quaternionic parameters of the group elements
  character(len=8),intent(out)           :: namele(iorder)
  !       namele  conventional names of the group elements
  integer(kind=i4_kind),intent(out) :: ncsco(3), fcsco(3)
  !       ncsco   indicates the class operator used to generated the
  !                 characters; for the numbering, see subroutine gencls
  !                 a negative value indicates an non-ambivalent class,
  !                 whence real and imaginary parts will be used as
  !                 separate operators
  !       fcsco   coefficient for the linear combination of class
  !                 operators
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of local constants --------------------
  integer(kind=i4_kind), parameter ::  modnum=1, nptgrp=75, ngrpar=304


  !------------ Declaration of local variables --------------------
  integer(kind=i4_kind) :: naxis, k, koffst, l, i, &
       ncscop(3,nptgrp), fcscop(3,nptgrp)
  character(len=36) :: grppar(ngrpar)
  logical ::appinv

  !
  !     elements of selected groups, angle/axis/parity coded
  !       parameter 1   : m/n -> angle = pi*m/n
  !       parameter 2-4 : unit vector along axis
  !                       coded c -> cos, s -> sin, angle as above,
  !                       special angles in group Oh:
  !                       ct1/2 -> cos of (tetrahedral angle /2),
  !                       sig, rho -> special angles of icosahedral group
  !                       see Altmann & Herzig, Table T75.1,
  !                       but beware of two typos !!!
  !                       0, 1, 1/2 just what you read
  !       parameter 5   : + pure rotation, - apply inversion
  !       parameter 6   : label of group element
  !

  data (grppar(k), k = 1, 2)                     /& ! Ci
       &       '0     0      0      1     + E      ', &
       &       '0     0      0      1     - i      '   /
  data (grppar(k), k = 3, 4)                     /& ! Cs
       &       '0     0      0      1     + E       ', &
       &       '1     0      0      1     - sigh    '  /
  data (grppar(k), k = 5, 12)                    / &! D2d
       &       '0     0      0      1     + E       ', &
       &       '1     0      0      1     + C2      ', &
       &       '1     1      0      0     + C2,1''  ', &
       &       '1     0      1      0     + C2,2''  ',  &
       &       '1/2   0      0      1     - S4-     ', &
       &       '1/2   0      0     -1     - S4+     ', &
       &       '1     c1/4   s1/4   0     - sigd1   ', &
       &       '1    -s1/4   c1/4   0     - sigd2   '  /
  data (grppar(k), k = 13, 24)                   /& ! D3h
       &       '0     0      0      1     + E       ', &
       &       '2/3   0      0      1     + C3+     ', &
       &       '2/3   0      0     -1     + C3-     ', &
       &       '1     1      0      0     + C2,1''  ', &
       &       '1    -c1/3   s1/3   0     + C2,2''  ', &
       &       '1    -c1/3  -s1/3   0     + C2,3''  ', &
       &       '1     0      0      1     - sigh    ',&
       &       '1/3   0      0      1     - S3-     ', &
       &       '1/3   0      0     -1     - S3+     ', &
       &       '1     0      1      0     - sigv1   ', &
       &       '1    -s1/3  -c1/3   0     - sigv2   ', &
       &       '1     s1/3  -c1/3   0     - sigv3   '  /
  data (grppar(k), k = 25, 40)                   / &! D4d
       &       '0     0      0      1     + E', &
       &       '1/2   0      0      1     + C4+     ',&
       &       '1/2   0      0     -1     + C4-     ', &
       &       '1     0      0      1     + C2      ', &
       &       '1     1      0      0     + C2,1''  ', &
       &       '1     0      1      0     + C2,2''  ', &
       &       '1     c1/4   s1/4   0     + C2,3''  ',&
       &       '1    -s1/4   c1/4   0     + C2,4''  ',&
       &       '1/4   0      0      1     - S8^3-   ',&
       &       '1/4   0      0     -1     - S8^3+   ',&
       &       '3/4   0      0      1     - S8-     ',&
       &       '3/4   0      0     -1     - S8+     ', &
       &       '1     c1/8   s1/8   0     - sigd1   ', &
       &       '1    -s1/8   c1/8   0     - sigd2   ', &
       &       '1     s1/8   c1/8   0     - sigd3   ', &
       &       '1    -c1/8   s1/8   0     - sigd4   '  /
  data (grppar(k), k = 41, 60)                   /& ! D5h
       &       '0     0      0      1     + E       ',  &
       &       '2/5   0      0      1     + C5+     ',  &
       &       '2/5   0      0     -1     + C5-     ',  &
       &       '4/5   0      0      1     + C5^2+   ',  &
       &       '4/5   0      0     -1     + C5^2-   ',  &
       &       '1     1      0      0     + C2,1''  ',  &
       &       '1     s1/10  c1/10  0     + C2,2''  ',  &
       &       '1    -c1/5   s1/5   0     + C2,3''  ',  &
       &       '1    -c1/5  -s1/5   0     + C2,4''  ',  &
       &       '1     s1/10 -c1/10  0     + C2,5''  ',  &
       &       '1/5   0      0     -1     - S5^2-   ',  &
       &       '1/5   0      0      1     - S5^2+   ',  &
       &       '3/5   0      0     -1     - S5-     ',  &
       &       '3/5   0      0      1     - S5+     ',  &
       &       '1     0      0      1     - sigh    ',  &
       &       '1     0      1      0     - sigv1   ',  &
       &       '1    -c1/10  s1/10  0     - sigv2   ',  &
       &       '1    -s1/5  -c1/5   0     - sigv3   ',  &
       &       '1     s1/5  -c1/5   0     - sigv4   ',  &
       &       '1     c1/10  s1/10  0     - sigv5   '  /
  data (grppar(k), k = 61, 84)                   /& ! D6d
       &       '0     0      0      1     + E       ',  &
       &       '1/3   0      0      1     + C6+     ',  &
       &       '1/3   0      0     -1     + C6-     ',  &
       &       '2/3   0      0      1     + C3+     ',  &
       &       '2/3   0      0     -1     + C3-     ',  &
       &       '1     0      0      1     + C2      ',  &
       &       '1     1      0      0     + C2,1''  ',  &
       &       '1    -c1/3   s1/3   0     + C2,2''  ',  &
       &       '1    -c1/3  -s1/3   0     + C2,3''  ',  &
       &       '1     0      1      0     + C2,4''  ',  &
       &       '1    -s1/3  -c1/3   0     + C2,5''  ',  &
       &       '1     s1/3  -c1/3   0     + C2,6''  ',  &
       &       '1/6   0      0      1     - S12^5-  ',  &
       &       '1/6   0      0     -1     - S12^5+  ',  &
       &       '1/2   0      0      1     - S4-     ',  &
       &       '1/2   0      0     -1     - S4+     ',  &
       &       '5/6   0      0      1     - S12-    ',  &
       &       '5/6   0      0     -1     - S12+    ',  &
       &       '1     c1/4   s1/4   0     - sigd1   ',  &
       &       '1    -c1/12  s1/12  0     - sigd2   ',  &
       &       '1     s1/12 -c1/12  0     - sigd3   ',  &
       &       '1    -s1/4   c1/4   0     - sigd4   ',  &
       &       '1    -s1/12 -c1/12  0     - sigd5   ',  &
       &       '1     c1/12  s1/12  0     - sigd6   '  /
  data (grppar(k), k = 85, 112)                  /& ! D7h
       &       '0     0      0      1     + E       ',  &
       &       '2/7   0      0      1     + C7+     ',  &
       &       '2/7   0      0     -1     + C7-     ',  &
       &       '4/7   0      0      1     + C7^2+   ',  &
       &       '4/7   0      0     -1     + C7^2-   ',  &
       &       '6/7   0      0      1     + C7^3+   ',  &
       &       '6/7   0      0     -1     + C7^3-   ',  &
       &       '1     1      0      0     + C2,1''  ',  &
       &       '1     c2/7   s2/7   0     + C2,2''  ',  &
       &       '1    -c3/7   s3/7   0     + C2,3''  ',  &
       &       '1    -c1/7   s1/7   0     + C2,4''  ',  &
       &       '1    -c1/7  -s1/7   0     + C2,5''  ',  &
       &       '1    -c3/7  -s3/7   0     + C2,6''  ',  &
       &       '1     c2/7  -s2/7   0     + C27''   ',  &
       &       '1     0      0      1     - sigh    ',  &
       &       '5/7   0      0     -1     - S7+     ',  &
       &       '5/7   0      0      1     - S7-     ',  &
       &       '3/7   0      0     -1     - S7^2+   ',  &
       &       '3/7   0      0      1     - S7^2-   ',  &
       &       '1/7   0      0     -1     - S7^3+   ',  &
       &       '1/7   0      0      1     - S7^3-   ',  &
       &       '1     0      1      0     - sigv1   ',  &
       &       '1    -s2/7   c2/7   0     - sigv2   ',  &
       &       '1    -s3/7  -c3/7   0     - sigv3   ',  &
       &       '1    -s1/7  -c1/7   0     - sigv4   ',  &
       &       '1     s1/7  -c1/7   0     - sigv5   ',  &
       &       '1     s3/7  -c3/7   0     - sigv6   ',  &
       &       '1     s2/7   c2/7   0     - sigv7   '  /
  data (grppar(k), k = 113, 144)                 /& ! D8d &
       &       '0     0      0      1     + E       ',  &
       &       '1/4   0      0      1     + C8+     ',  &
       &       '1/4   0      0     -1     + C8-     ',  &
       &       '1/2   0      0      1     + C4+     ',  &
       &       '1/2   0      0     -1     + C4-     ',  &
       &       '3/4   0      0      1     + C8^3+   ',  &
       &       '3/4   0      0     -1     + C8^3-   ',  &
       &       '1     0      0      1     + C2      ',  &
       &       '1     1      0      0     + C2,1''  ',  &
       &       '1     0      1      0     + C2,2''  ',  &
       &       '1     c1/4   s1/4   0     + C2,3''  ',  &
       &       '1    -s1/4   c1/4   0     + C2,4''  ',  &
       &       '1     c1/8   s1/8   0     + C2,5''  ',  &
       &       '1    -s1/8   c1/8   0     + C2,6''  ',  &
       &       '1     s1/8   c1/8   0     + C2,7''  ',  &
       &       '1    -c1/8   s1/8   0     + C2,8''  ', &
       &       '1/8   0      0      1     - S16^7-  ',  &
       &       '1/8   0      0     -1     - S16^7+  ',  &
       &       '3/8   0      0      1     - S16^5-  ',  &
       &       '3/8   0      0     -1     - S16^5+  ',  &
       &       '5/8   0      0      1     - S16^3-  ',  &
       &       '5/8   0      0     -1     - S16^3+  ',  &
       &       '7/8   0      0      1     - S16-    ',  &
       &       '7/8   0      0     -1     - S16+    ',  &
       &       '1     c1/16  s1/16  0     - sigd1   ',  &
       &       '1    -s1/16  c1/16  0     - sigd2   ',  &
       &       '1     s3/16  c3/16  0     - sigd3   ',  &
       &       '1    -c3/16  s3/16  0     - sigd4   ',  &
       &       '1     c3/16  s3/16  0     - sigd5   ',  &
       &       '1    -s3/16  c3/16  0     - sigd6   ',  &
       &       '1     s1/16  c1/16  0     - sigd7   ',  &
       &       '1    -c1/16  s1/16  0     - sigd8   '  /
  data (grppar(k), k = 145, 180)                 /& ! D9h &
       &       '0     0      0      1     + E       ',  &
       &       '2/9   0      0      1     + C9+     ',  &
       &       '2/9   0      0     -1     + C9-     ',  &
       &       '4/9   0      0      1     + C9^2+   ',  &
       &       '4/9   0      0     -1     + C9^2-   ',  &
       &       '2/3   0      0      1     + C3+     ',  &
       &       '2/3   0      0     -1     + C3-     ',  &
       &       '8/9   0      0      1     + C9^4+   ',  &
       &       '8/9   0      0     -1     + C9^4-   ',  &
       &       '1     1      0      0     + C2,1''  ',  &
       &       '1    -c1/3   s1/3   0     + C2,2''  ',  &
       &       '1    -c1/3  -s1/3   0     + C2,3''  ',  &
       &       '1     c2/9   s2/9   0     + C2,4''  ',  &
       &       '1    -c1/9   s1/9   0     + C2,5''  ',  &
       &       '1     c4/9  -s4/9   0     + C2,6''  ',  &
       &       '1     c4/9   s4/9   0     + C2,7''  ',  &
       &       '1    -c1/9  -s1/9   0     + C2,8''  ',  &
       &       '1     c2/9  -s2/9   0     + C2,9''  ',  &
       &       '1     0      0      1     - sigh    ',  &
       &       '7/9   0      0     -1     - S9+     ',  &
       &       '7/9   0      0      1     - S9-     ',  &
       &       '5/9   0      0     -1     - S9^2+   ',  &
       &       '5/9   0      0      1     - S9^2-   ',  &
       &       '1/3   0      0     -1     - S3+     ',  &
       &       '1/3   0      0      1     - S3-     ',  &
       &       '1/9   0      0     -1     - S9^4+   ',  &
       &       '1/9   0      0      1     - S9^4-   ',  &
       &       '1     0      1      0     - sigv1   ',  &
       &       '1    -s1/3  -c1/3   0     - sigv2   ',  &
       &       '1     s1/3  -c1/3   0     - sigv3   ',  &
       &       '1    -s2/9   c2/9   0     - sigv4   ',  &
       &       '1    -s1/9  -c1/9   0     - sigv5   ',  &
       &       '1     s4/9   c4/9   0     - sigv6   ',  &
       &       '1    -s4/9   c4/9   0     - sigv7   ',  &
       &       '1     s1/9  -c1/9   0     - sigv8   ',  &
       &       '1     s2/9   c2/9   0     - sigv9   '  /
  data (grppar(k), k = 181, 216)                 /& ! D10d
       &       '0     0      0      1     + E       ',  &
       &       '1/5   0      0      1     + C10+    ',  &
       &       '1/5   0      0     -1     + C10-    ',  &
       &       '2/5   0      0      1     + C5+     ',  &
       &       '2/5   0      0     -1     + C5-     ',  &
       &       '3/5   0      0      1     + C10^3+  ',  &
       &       '3/5   0      0     -1     + C10^3-  ',  &
       &       '4/5   0      0      1     + C5^2+   ',  &
       &       '4/5   0      0     -1     + C5^2-   ',  &
       &       '1     0      0      1     + C2      ',  &
       &       '1     1      0      0     + C2,1''  ',  &
       &       '1     s1/10  c1/10  0     + C2,2''  ',  &
       &       '1    -c1/5   s1/5   0     + C2,3''  ',  &
       &       '1    -c1/5  -s1/5   0     + C2,4''  ',  &
       &       '1     s1/10 -c1/10  0     + C2,5''  ',  &
       &       '1     0      1      0     + C2,6''  ',  &
       &       '1    -c1/10  s1/10  0     + C2,7''  ',  &
       &       '1    -s1/5  -c1/5   0     + C2,8''  ',  &
       &       '1     s1/5  -c1/5   0     + C2,9''  ',  &
       &       '1     c1/10  s1/10  0     + C2,10'' ',  &
       &       '1/10  0      0      1     - S20^9-  ',  &
       &       '1/10  0      0     -1     - S20^9+  ',  &
       &       '3/10  0      0      1     - S20^7-  ',  &
       &       '3/10  0      0     -1     - S20^7+  ',  &
       &       '1/2   0      0      1     - S4-     ',  &
       &       '1/2   0      0     -1     - S4+     ',  &
       &       '7/10  0      0      1     - S20^3-  ',  &
       &       '7/10  0      0     -1     - S20^3+  ',  &
       &       '9/10  0      0      1     - S20-    ',  &
       &       '9/10  0      0     -1     - S20+    ',  &
       &       '1     c1/4   s1/4   0     - sigd1   ',  &
       &       '1    -s3/20  c3/20  0     - sigd2   ',  &
       &       '1    -c1/20 -s1/20  0     - sigd3   ',  &
       &       '1    -s1/20 -c1/20  0     - sigd4   ',  &
       &       '1     c3/20 -s3/20  0     - sigd5   ',  &
       &       '1    -s1/4   c1/4   0     - sigd6   ' /
  data (grppar(k), k= 217,220) / &
       &       '1    -c3/20 -s3/20  0     - sigd7   ', &
       &       '1     s1/20 -c1/20  0     - sigd8   ', &
       &       '1     c1/20 -s1/20  0     - sigd9   ', &
       &       '1     s3/20  c3/20  0     - sigd10  '  /
  data (grppar(k), k = 221, 244)                 /& ! O
       &       '0     0      0      1     + E       ', &
       &       '1     1      0      0     + C2x     ', &
       &       '1     0      1      0     + C2y     ', &
       &       '1     0      0      1     + C2z     ', &
       &       '2/3   ct/2   ct/2   ct/2  + C3,1+   ', &
       &       '2/3  -ct/2  -ct/2   ct/2  + C3,2+   ', &
       &       '2/3   ct/2  -ct/2  -ct/2  + C3,3+   ', &
       &       '2/3  -ct/2   ct/2  -ct/2  + C3,4+   ', &
       &       '2/3  -ct/2  -ct/2  -ct/2  + C3,1-   ', &
       &       '2/3   ct/2   ct/2  -ct/2  + C3,2-   ', &
       &       '2/3  -ct/2   ct/2   ct/2  + C3,3-   ', &
       &       '2/3   ct/2  -ct/2   ct/2  + C3,4-   ', &
       &       '1/2   1      0      0     + C4x+    ', &
       &       '1/2   0      1      0     + C4y+    ', &
       &       '1/2   0      0      1     + C4z+    ', &
       &       '1/2  -1      0      0     + C4x-    ', &
       &       '1/2   0     -1      0     + C4y-    ', &
       &       '1/2   0      0     -1     + C4z-    ', &
       &       '1     c1/4   s1/4   0     + C2a''   ', &
       &       '1    -s1/4   c1/4   0     + C2b''   ', &
       &       '1     s1/4   0      c1/4  + C2c''   ', &
       &       '1     0     -c1/4  -s1/4  + C2d''   ', &
       &       '1     c1/4   0     -s1/4  + C2e''   ', &
       &       '1     0     -s1/4   c1/4  + C2f''   '  /
  data (grppar(k), k = 245, 274)                 /& ! I
       &       '0     0      0      1     + E       ', &
       &       '2/5   csig   0      ssig  + C5,1+   ', &
       &       '2/5  -csig   0      ssig  + C5,2+   ', &
       &       '2/5   0      ssig   csig  + C5,3+   ', &
       &       '2/5   0     -ssig   csig  + C5,4+   ', &
       &       '2/5   ssig   csig   0     + C5,5+   ', &
       &       '2/5   ssig  -csig   0     + C5,6+   ', &
       &       '2/5  -csig   0     -ssig  + C5,1-   ', &
       &       '2/5   csig   0     -ssig  + C5,2-   ', &
       &       '2/5   0     -ssig  -csig  + C5,3-   ', &
       &       '2/5   0      ssig  -csig  + C5,4-   ', &
       &       '2/5  -ssig  -csig   0     + C5,5-   ', &
       &       '2/5  -ssig   csig   0     + C5,6-   ', &
       &       '4/5   csig   0      ssig  + C5,1^2+ ', &
       &       '4/5  -csig   0      ssig  + C5,2^2+ ', &
       &       '4/5   0      ssig   csig  + C5,3^2+ ', &
       &       '4/5   0     -ssig   csig  + C5,4^2+ ', &
       &       '4/5   ssig   csig   0     + C5,5^2+ ', &
       &       '4/5   ssig  -csig   0     + C5,6^2+ ', &
       &       '4/5  -csig   0     -ssig  + C5,1^2- ', &
       &       '4/5   csig   0     -ssig  + C5,2^2- ', &
       &       '4/5   0     -ssig  -csig  + C5,3^2- ', &
       &       '4/5   0      ssig  -csig  + C5,4^2- ', &
       &       '4/5  -ssig  -csig   0     + C5,5^2- ', &
       &       '4/5  -ssig   csig   0     + C5,6^2- ', &
       &       '2/3   ct/2   ct/2   ct/2  + C3,1+   ', &
       &       '2/3  -ct/2  -ct/2   ct/2  + C3,2+   ', &
       &       '2/3   ct/2  -ct/2  -ct/2  + C3,3+   ', &
       &       '2/3  -ct/2   ct/2  -ct/2  + C3,4+   ', &
       &       '2/3  -srho   0      crho  + C3,5+   '  /
  data (grppar(k) , k=275,304)  / &
       &       '2/3   0      crho  -srho  + C3,6+   ', &
       &       '2/3   crho  -srho   0     + C3,7+   ', &
       &       '2/3   srho   0      crho  + C3,8+   ', &
       &       '2/3   0      crho   srho  + C3,9+   ', &
       &       '2/3   crho   srho   0     + C3,10+  ', &
       &       '2/3  -ct/2  -ct/2  -ct/2  + C3,1-   ', &
       &       '2/3   ct/2   ct/2  -ct/2  + C3,2-   ', &
       &       '2/3  -ct/2   ct/2   ct/2  + C3,3-   ', &
       &       '2/3   ct/2  -ct/2   ct/2  + C3,4-   ', &
       &       '2/3   srho   0     -crho  + C3,5-   ', &
       &       '2/3   0     -crho   srho  + C3,6-   ', &
       &       '2/3  -crho   srho   0     + C3,7-   ', &
       &       '2/3  -srho   0     -crho  + C3,8-   ', &
       &       '2/3   0     -crho  -srho  + C3,9-   ', &
       &       '2/3  -crho  -srho   0     + C3,10-  ', &
       &       '1     1      0      0     + C2a     ', &
       &       '1     0      1      0     + C2b     ', &
       &       '1     0      0      1     + C2c     ', &
       &       '1     c1/5   s1/10  1/2   + C2d     ', &
       &       '1     1/2   -c1/5  -s1/10 + C2e     ', &
       &       '1     s1/10  1/2   -c1/5  + C2f     ', &
       &       '1     s1/10  1/2    c1/5  + C2g     ', &
       &       '1     1/2   -c1/5   s1/10 + C2h     ', &
       &       '1     c1/5   s1/10 -1/2   + C2i     ', &
       &       '1    -s1/10  1/2    c1/5  + C2j     ', &
       &       '1    -c1/5   s1/10 -1/2   + C2k     ', &
       &       '1    -1/2   -c1/5   s1/10 + C2l     ', &
       &       '1    -c1/5   s1/10  1/2   + C2m     ', &
       &       '1    -s1/10  1/2   -c1/5  + C2n     ', &
       &       '1    -1/2   -c1/5  -s1/10 + C2o     '  /
  !
  !    For each group, ncscop denotes the classes making up a CSCO and
  !      fcscop provides the coefficient in the linear combination of the
  !      class operators. A negative number indicates a non-ambivalent
  !      class whence the real and the imaginary part will be used as
  !      separate class operators. The list of the classes will be
  !      generated in subroutine gencls.
  !
  data ((ncscop(i,k), i = 1,3), k = 1, 18) &
       & /  1, 0, 0,  2, 0, 0, -2, 0, 0, -2, 0, 0, -2, 0, 0, -2, 0, 0,&
       &   -2, 0, 0, -2, 0, 0, -2, 0, 0, -2, 0, 0, -2, 0, 0, -2, 0, 0, &
       &   -4, 0, 0, -6, 0, 0, -8, 0, 0,-10, 0, 0,-12, 0, 0,-14, 0, 0 /
  data ((fcscop(i,k), i = 1,3), k = 1, 18) &
       & /  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,&
       &    1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0, &
       &    1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0 /
  data ((ncscop(i,k), i = 1,3), k = 19, 36) &
       & /-16, 0, 0,-18, 0, 0,-20, 0, 0,  2, 3, 0,  3, 0, 0,  2, 4, 0, &
       &    2, 4, 0,  2, 5, 0,  2, 5, 0,  2, 6, 0,  2, 6, 0,  2, 7, 0, &
       &    2, 3, 5,  3, 4, 0,  2, 4, 6,  2, 4, 7,  2, 5, 7,  2, 5, 6 /
  data ((fcscop(i,k), i = 1,3), k = 19, 36) &
       & /  1, 0, 0,  1, 0, 0,  1, 0, 0,  2, 1, 0,  1, 0, 0,  2, 1, 0,&
       &    2, 1, 0,  2, 1, 0,  2, 1, 0,  2, 2, 0,  2, 1, 0,  2, 1, 0, &
       &    6, 2, 1,  1, 1, 0,  4, 1, 1,  4, 1, 1,  4, 1, 1,  4, 1, 1 /
  data ((ncscop(i,k), i = 1,3), k = 37, 54) &
       & /  2, 6, 8,  2, 6, 7,  2, 7, 9,  2, 4, 6,  4, 5, 0,  4, 6, 0, &
       &    6, 7, 0,  7, 8, 0,  8, 9, 0,  9,10, 0, 10,11, 0, 11,12, 0, &
       &   12,13, 0,  2, 3, 0,  3, 0, 0,  2, 4, 0,  2, 4, 0,  2, 5, 0 /
  data ((fcscop(i,k), i = 1,3), k = 37, 54) &
       & /  2, 2, 1,  2, 1, 1,  4, 2, 1,  2, 1, 1,  2, 1, 0,  2, 1, 0,&
       &    3, 1, 0,  2, 1, 0,  4, 1, 0,  2, 1, 0,  2, 1, 0,  2, 1, 0, &
       &    2, 1, 0,  2, 1, 0,  1, 0, 0,  2, 1, 0,  2, 1, 0,  2, 1, 0 /
  data ((ncscop(i,k), i = 1,3), k = 55, 72) &
       & /  2, 5, 0,  2, 6, 0,  2, 6, 0,  2, 7, 0,  2, 4, 0,  2, 3, 0, &
       &   -2, 4, 0, -2, 5, 0, -2,10, 0, -2, 7, 0, -2, 8, 0, -2, 9, 0, &
       &   -2,10, 0, -2,11, 0,  5, 0, 0, -3, 0, 0,  5, 6, 0, -3, 5, 0 /
  data ((fcscop(i,k), i = 1,3), k = 55, 72) &
       & /  2, 1, 0,  2, 2, 0,  2, 1, 0,  2, 1, 0,  2, 1, 0,  2, 1, 0,&
       &    1, 5, 0,  2, 5, 0,  2, 5, 0,  2, 5, 0,  2, 5, 0,  2, 5, 0, &
       &    2, 5, 0,  2, 5, 0,  1, 0, 0,  1, 0, 0,  2, 1, 0,  1, 5, 0 /
  data ((ncscop(i,k), i = 1,3), k = 73, 75) &
       & /  5, 0, 0,  2, 5, 0,  2, 5, 6 /
  data ((fcscop(i,k), i = 1,3), k = 73, 75) &
       & /  1, 0, 0,  2, 1, 0,  2, 1, 1 /

  !
  ! FIXME: for some reason the subroutine suplab() is called
  ! with uninitialized names of elements from intent(out)
  ! array namele(:) as input, this seems to lead to some
  ! artifacts. Initialize them to empty strings as a workaround:
  !
  do i = 1, size(namele)
    namele(i) = "        "
  enddo

  appinv = .false.
  !
  !     selects proper subgroup or supplement with inversion etc.
  !       then decodes and forms quaternionic parameters elem as well as
  !       label of element from grppar. Where appropriate, supnam is
  !       called for labels of derived elements (i.e. those generated by
  !       applying the inversion)
  !
  if ( igroup.ge.1 .and. igroup.le.10 ) then         ! Cn
     naxis   = iorder              ! order of main rotation axis
     koffst  = 2*naxis*(naxis-1)   ! offset in list grppar
     do k = 1, iorder
        call decdel (grppar(koffst+k),elem(1,k))
        namele(k) = grppar(koffst+k) (29:)
     end do
  else if ( igroup.eq.11 ) then                      ! Ci
     call decdel (grppar(1),elem(1,1))
     namele(1) = grppar(1) (29:)
     call decdel (grppar(2),elem(1,2))
     namele(2) = grppar(2) (29:)
  else if ( igroup.eq.12 ) then                      ! Cs
     call decdel (grppar(3),elem(1,1))
     namele(1) = grppar(3) (29:)
     call decdel (grppar(4),elem(1,2))
     namele(2) = grppar(4) (29:)
  else if ( igroup.ge.13 .and. igroup.le.21 ) then   ! S2n
     naxis   = iorder/2
     koffst  = 2*naxis*(naxis-1)
     do k = 1, naxis
        call decdel (grppar(koffst+k),elem(1,k))
        namele(k) = grppar(koffst+k) (29:)
     end do
     if ( mod(naxis,2) .eq. 0 ) then
        koffst = koffst + iorder
        do k = 1, naxis
           call decdel (grppar(koffst+k),elem(1,naxis+k))
           namele(naxis+k) = grppar(koffst+k) (29:)
        end do
     else
        do k = 1, naxis
           do l = 1, 4
              elem(l,naxis+k) = elem(l,k)
           end do
           elem(5,naxis+k) = - elem(5,k)
           appinv = .true.
           call suplab (namele(k),namele(naxis+k),group,naxis, &
                &                      appinv)
        end do
     end if
  else if ( igroup.ge.22 .and. igroup.le.30 ) then   ! Dn
     naxis   = iorder/2
     koffst  = 2*naxis*(naxis-1)
     do k = 1, iorder
        call decdel (grppar(koffst+k),elem(1,k))
        namele(k) = grppar(koffst+k) (29:)
     end do
     if ( naxis .eq. 2 ) then
        call suplab (namele(2),namele(2),group,naxis,appinv)
        call suplab (namele(3),namele(3),group,naxis,appinv)
        call suplab (namele(4),namele(4),group,naxis,appinv)
     else if ( mod(naxis,2) .eq. 0 .and. naxis .gt. 2 ) then
        do k = 3*naxis/2+1, iorder
           call suplab (namele(k),namele(k),group,naxis,appinv)
        end do
     end if
  else if ( igroup.ge.31 .and. igroup.le.39 ) then   ! Dnh
     naxis   = iorder/4
     koffst  = 2*naxis*(naxis-1)
     if ( mod(naxis,2) .eq. 0 ) then
        do k = 1, iorder/2
           call decdel (grppar(koffst+k),elem(1,k))
           namele(k) = grppar(koffst+k) (29:)
        end do
        appinv = .false.
        if ( naxis .eq. 2 ) then
           call suplab (namele(2),namele(2),group,naxis,appinv)
           call suplab (namele(3),namele(3),group,naxis,appinv)
           call suplab (namele(4),namele(4),group,naxis,appinv)
        else
           do k = 3*naxis/2+1, iorder
              call suplab (namele(k),namele(k),group,naxis,appinv)
           end do
        end if
        appinv = .true.
        do k = 1, iorder/2
           do l = 1, 4
              elem(l,iorder/2+k) = elem(l,k)
           end do
           elem(5,iorder/2+k) = - elem(5,k)
           call suplab (namele(k),namele(iorder/2+k),group, &
                &                      naxis,appinv)
        end do
     else
        do k = 1, iorder
           call decdel (grppar(koffst+k),elem(1,k))
           namele(k) = grppar(koffst+k) (29:)
        end do
     end if
  else if ( igroup.ge.41 .and. igroup.le.49 ) then   ! Dnd
     naxis   = iorder/4
     koffst  = 2*naxis*(naxis-1)
     if ( mod(naxis,2).eq.0 ) then
        do k = 1, iorder
           call decdel (grppar(koffst+k),elem(1,k))
           namele(k) = grppar(koffst+k) (29:)
        end do
     else
        do k = 1, iorder/2
           call decdel (grppar(koffst+k),elem(1,k))
           namele(k) = grppar(koffst+k) (29:)
           do l = 1, 4
              elem(l,iorder/2+k) = elem(l,k)
           end do
           elem(5,iorder/2+k) = - elem(5,k)
           appinv = .false.
           call suplab (namele(k),namele(iorder/2+k),group,naxis, &
                &                     appinv)
        end do
     end if
  else if ( igroup.ge.50 .and. igroup.le.58 ) then   ! Cnv
     naxis   = iorder/2
     koffst  = 2*naxis*(naxis-1)
     do k = 1, naxis
        call decdel (grppar(koffst+k),elem(1,k))
        namele(k) = grppar(koffst+k) (29:)
     end do
     if ( mod(naxis,2).eq.0 ) then
        koffst  = koffst + naxis
        appinv = .true.
        do k = 1, naxis
           call decdel (grppar(koffst+k),elem(1,naxis+k))
           elem(5,naxis+k) = - elem(5,k)
           namele(naxis+k) = grppar(koffst+k) (29:)
           call suplab (namele(naxis+k),namele(naxis+k),group, &
                &                      naxis,appinv)
        end do
     else
        koffst  = koffst + 3*naxis
        do k = 1, naxis
           call decdel (grppar(koffst+k),elem(1,naxis+k))
           namele(naxis+k) = grppar(koffst+k) (29:)
        end do
     end if
  else if ( igroup.ge.60 .and. igroup.le.68 ) then   ! Cnh
     naxis   = iorder/2
     koffst  = 2*naxis*(naxis-1)
     do k = 1, naxis
        call decdel (grppar(koffst+k),elem(1,k))
        namele(k) = grppar(koffst+k) (29:)
     end do
     if ( mod(naxis,2).eq.0 ) then
        koffst  = koffst + naxis
        appinv = .true.
        do k = 1, naxis
           do l = 1, 4
              elem(l,naxis+k) = elem(l,k)
           end do
           elem(5,naxis+k) = - elem(5,k)
           call suplab (namele(k),namele(naxis+k),group,naxis, &
                &                     appinv)
        end do
     else
        koffst  = koffst + 2*naxis
        do k = 1, naxis
           call decdel (grppar(koffst+k),elem(1,naxis+k))
           namele(naxis+k) = grppar(koffst+k) (29:)
        end do
     end if
  else if ( igroup.ge.69 .and. igroup.le.73 ) then   ! octahedral
     koffst = 220
     if ( igroup.eq.69 .or. igroup.eq.70 ) then
        naxis = iorder                             ! O or T
     else
        naxis = iorder/2                           ! Oh, Th, Td
     end if
     do k = 1, naxis
        call decdel (grppar(koffst+k),elem(1,k))
        namele(k) = grppar(koffst+k) (29:)
     end do
     if ( igroup.eq.71 .or. igroup.eq.72 ) then    ! Oh or Th
        appinv = .true.
        do k = 1, naxis
           do l = 1, 4
              elem(l,naxis+k) = elem(l,k)
           end do
           elem(5,naxis+k) = - elem(5,k)
           call suplab (namele(k),namele(naxis+k),group,naxis,&
                &                     appinv)
        end do
     else if ( igroup.eq. 73 ) then                ! Td
        koffst = koffst + naxis
        appinv = .true.
        do k = 1, naxis
           call decdel (grppar(koffst+k),elem(1,naxis+k))
           elem(5,naxis+k) = - elem(5,k)
           namele(naxis+k) = grppar(koffst+k) (29:)
           call suplab (namele(naxis+k),namele(naxis+k),group, &
                &                     naxis,appinv)
        end do
     end if
  else if ( igroup.eq.74 .or.  igroup.eq.75 ) then   ! icosahedral
     koffst = 244
     if ( igroup.eq.74 ) then
        naxis = iorder                             ! I
     else
        naxis = iorder/2                           ! Ih
     end if
     do k = 1, naxis
        call decdel (grppar(koffst+k),elem(1,k))
        namele(k) = grppar(koffst+k) (29:)
     end do
     if ( igroup.eq.75 ) then                      ! Ih
        appinv = .true.
        do k = 1, naxis
           do l = 1, 4
              elem(l,naxis+k) = elem(l,k)
           end do
           elem(5,naxis+k) = - elem(5,k)
           call suplab (namele(k),namele(naxis+k),group,naxis,&
                &                     appinv)
        end do
     end if
  end if
  !
  do l = 1, 3
     ncsco(l) = ncscop (l,igroup)
     fcsco(l) = fcscop (l,igroup)
  end do

contains

  subroutine addinv (namold, namnew)
    implicit none
    !
    !        author:    N. Roesch, Theoretische Chemie, TU Muenchen
    !        version:   30.04.95
    !
    !     ... subroutine parameters
    character(len=8), intent(in)  :: namold
    character(len=8), intent(out) :: namnew
    !
    !     renames a rotation as required when the inversion is applied
    !       based on the representation Cn+, Cn-, Cn^m+, Cn^m- where
    !       the rotation angle is +-2*pi*m/n. The first two forms are
    !       used in case of m = 1. E and C2 are treated explicitly
    !       in subroutine supnam and thus are excluded here. The
    !       subroutine relies on the ASCII collating sequence.
    !
    !       namold  old name
    !       namnew  new name
    !
    !     ... integer
    integer       possig, poscar, m, n, posnup, next
    !     ... character
    character(len=8)   namcha
    !
    possig = index(namold,'+') + index(namold,'-')
    if ( possig .eq. 0 ) then
       if (namold(1:1) .eq. 'E' ) then
          namcha = 'i'
       else if ( namold(1:3) .eq. 'C2 ' ) then
          namcha = 'sigh'
       end if
    else
       poscar = index(namold,'^')
       if (poscar .eq. 0) then
          m = 1
          posnup = possig - 1
       else
          m = ichar(namold(poscar+1:poscar+1)) - 48
          posnup = poscar - 1
       end if
       n = ichar(namold(2:2)) - 48
       if ( posnup .eq. 3 ) then
          n = n*10 + ichar(namold(3:3)) - 48
       end if
       !...   generate improper rotation and eliminate factors if necessary
       m = n - 2*m
       n = n*2
       if (mod(m,2) .eq. 0) then
          m = m/2
          n = n/2
       end if
       if (mod(m,2) .eq. 0) then
          m = m/2
          n = n/2
       end if
       namcha = 'S'
       if (n .gt. 9) then
          namcha(2:2) = char(48+n/10)
          namcha(3:3) = char(48+n-(n/10)*10)
          next = 4
       else
          namcha(2:2) = char(48+n)
          next = 3
       end if
       if (m .gt. 1) then
          namcha(next:next) = '^'
          namcha(next+1:next+1) = char(48+m)
          next = next + 2
       end if
       !...  reverse sense of rotation
       if (namold(possig:possig) .eq. '+') then
          namcha(next:next) = '-'
       else
          namcha(next:next) = '+'
       end if
    end if
    !
    namnew = namcha
    !
  end subroutine addinv


  subroutine decdel (grppar,elem)
    implicit none
    !
    !        author:    N. Roesch, Theoretische Chemie, TU Muenchen
    !        version:   30.04.95
    !
    !     ... subroutine parameters
!!$    integer, parameter :: r8_kind = selected_real_kind(15)
    character(len=36)  grppar
    real(r8_kind) ::   elem(5)
    !
    !     decodes the group element grppar (given in the data statements
    !       of subroutine gengrp) and converts angle/axis parameters into
    !       quaternionic parameters elem:
    !       cos(phi/2), sin(phi/2)*n1, sin(phi/2)*n2, sin(phi/2)*n3), parity
    !       The subroutine relies on the ASCII collating sequence.
    !
    !     ... integer
    integer       k, m, n
    !     ... real
    real(r8_kind) ::        angle, phi2, pi
    real(r8_kind),parameter ::&
         & zero  = 0.0_r8_kind,&
         & one   = 1.0_r8_kind,&
         & two   = 2.0_r8_kind,&
         & three = 3.0_r8_kind,&
         & five  = 5.0_r8_kind
    !     ... character
    character(len=4)   chphi
    character(len=5)   chang
    !
    pi = 4.0_r8_kind*atan(1.0_r8_kind)
    !
    chphi = grppar(1:4)
    m = iachar(chphi(1:1)) - 48
    if ( chphi(2:2).eq.'/' ) then
       if ( chphi(4:4).ne.' ' ) then
          n  = 10*(iachar(chphi(3:3))-48) + (iachar(chphi(4:4))-48)
       else
          n  = iachar(chphi(3:3)) - 48
       endif
    else
       n = 1
    endif
    phi2 = ((pi/two) * m) / n
    !
    do k = 1, 3
       chang = grppar(7*k:7*k+4)
       if ( chang.eq.'1/2  ' ) then
          elem(k+1) = sin(phi2) / two
       else if ( chang.eq.'ct/2 ' ) then
          elem(k+1) = sin(phi2) / sqrt(three)
       else if ( chang.eq.'1    ' ) then
          elem(k+1) = sin(phi2)
       else if ( chang.eq.'0    ' ) then
          elem(k+1) = zero
       else
          if ( chang(2:4).eq.'rho' ) then
             angle = atan((sqrt(five) + three) / two)
          else if ( chang(2:4).eq.'sig' ) then
             angle = atan((sqrt(five) + one  ) / two)
          else
             m = iachar(chang(2:2)) - 48
             if ( chang(5:5).ne.' ' ) then
                n = 10*(iachar(chang(4:4))-48) &
                     &                + (iachar(chang(5:5))-48)
             else
                n = iachar(chang(4:4)) - 48
             endif
             angle = (pi * m) / n
          end if
          if ( chang(1:1).eq.'c' ) then
             elem(k+1) =  sin(phi2) * cos(angle)
          else
             elem(k+1) =  sin(phi2) * sin(angle)
          end if
       end if
       if ( grppar(7*k-1:7*k-1).eq.'-' ) then
          elem(k+1) = -elem(k+1)
       end if
    end do
    !
    elem(1) = cos(phi2)
    if (grppar(27:27).eq.'-' ) then
       elem(5) = -one
    else
       elem(5) =  one
    end if
    !
  end subroutine decdel


  subroutine suplab (namold,namnew,group,naxis,appinv)
    implicit none
    !
    !     part of the module efm_mo
    !        author:    N. Roesch, Theoretische Chemie, TU Muenchen
    !        version:   01.05.95
    !
    !     ... subroutine parameters
    character(len=8), intent(in)  :: namold
    character(len=8), intent(out) :: namnew
    character(len=4), intent(in)  :: group
    integer    , intent(in)  :: naxis
    logical    , intent(in)  :: appinv
    !
    !     supplies additinoal names for group elements not explicitly
    !       given in the data statements of subroutine gengrp,
    !       e.g. for all improper rotations. In most cases, the inversion
    !       is applied automically by subroutine addinv, but E and C2
    !       are handled explictly. In all cases, the names have been
    !       chosen to closely follow the tables of Altmann and Herzig.
    !       The subroutine relies on the ASCII collating sequence.
    !
    !       namold  old name
    !       namnew  new name
    !       group   name of the point group
    !       naxis   order of the main axis
    !       appinv  apply inversion (evaluated only for certain groups !)
    !
    !     ... subroutines called
    !     addinv
    !     ... integer
    integer       possig, poscar, n, poskup
    !     ... character
    character(len=8)   namcha
    character(len=3)   kvalue
    !
    if (group(1:1) .eq. 'S') then                         ! S2n
       !...  groups S2n, n even
       !...  apply inversion to proper rotations
       if (namold(1:1) .eq. 'E') then
          namcha = 'i       '
       else
          call addinv (namold,namcha)
       end if
    else if (group(1:2) .eq. 'D2') then             ! D2, D2h
       !...  groups D2, D2h
       !...  apply inversion to proper rotations, use special names
       if (namold(1:3) .eq. 'C2 ') then
          namcha = 'C2z'
       else if (namold(1:5) .eq. 'C2,1''') then
          namcha = 'C2x'
       else if (namold(1:5) .eq. 'C2,2''') then
          namcha = 'C2y'
       else if (namold(1:1) .eq. 'E') then
          namcha = 'i'
       else if (namold(1:3) .eq. 'C2z') then
          namcha = 'sigz'
       else if (namold(1:3) .eq. 'C2x') then
          namcha = 'sigx'
       else if (namold(1:3) .eq. 'C2y') then
          namcha = 'sigy'
       end if
    else if ( group(1:1).eq.'D' .and. naxis.gt.2 .and. &
         &          mod(naxis,2).eq.0 .and. .not.appinv ) then  ! Dn, n>2
       !...  groups Dn, Dnh and Dnd, n even and n > 2
       !...  generates the labels for the second class C2,k" of C2 operations
       !...  besides C2,k'
       poscar = index(namold,',')
       possig = index(namold,'''')
       n = ichar(namold(poscar+1:poscar+1)) - 48
       if ( poscar+1 .lt. possig-1 ) then
          n = n*10 + ichar(namold(poscar+2:poscar+2)) - 48
       end if
       n = n - naxis/2
       namcha = namold(:poscar)//char(48+n)//"''"
    else if ( group(1:1).eq.'D' .and. naxis.gt.2 .and. &
         &          index(group,'H') .gt. 0 ) then              ! Dnh, n>2
       !...  groups Dnh, n even and n > 2
       !...  apply inversion to proper rotations, generate labels for the
       !...  two sets of vertical mirror planes, sigv and sigd
       poscar = index(namold,',')
       if (poscar .eq. 0 ) then
          call addinv (namold,namcha)
       else
          possig = index(namold,'''')
          namcha = 'sig '//namold(poscar+1:poscar+1)
          if ( mod(naxis,4) .eq. 0 ) then
             if ( possig .gt. 0 ) then
                namcha(4:4) = 'V'
             else
                namcha(4:4) = 'D'
             end if
          else
             if ( possig .gt. 0 ) then
                namcha(4:4) = 'D'
             else
                namcha(4:4) = 'V'
             end if
          end if
       end if
    else if ( group(1:1).eq.'D' .and. mod(naxis,2).ne.0 &
         &          .and. index(group(2:),'D') .gt. 0 ) then        ! Dnd, n odd
       !...  groups Dnd, n odd and n > 2
       !...  apply inversion to proper rotations, generate labels for the
       !...  two class of vertical mirror planes, sigd.
       poscar = index(namold,',')
       if (poscar .eq. 0 ) then
          call addinv (namold,namcha)
       else
          possig = index(namold,'''')
          namcha = 'sigd'//namold(poscar+1:poscar+1)
       end if
    else if ( group(1:1).eq.'C' .and. mod(naxis,2).eq.0 &
         &          .and. index(group,'V') .gt. 0 ) then        ! Cnv, n even
       !...  groups Cnv, n even
       !...  generate labels for the two sets of vertical mirror planes,
       !...  sigvk' and sigvk", C2v is handled explicitly.
       if ( naxis .eq. 2 ) then
          if ( namold(1:5) .eq. 'C2,1''') then
             namcha = 'sigx'
          else if ( namold(1:5) .eq. 'C2,2''') then
             namcha = 'sigy'
          end if
       else if ( namold(1:3) .eq. 'C2,' ) then
          possig = index(namold,'''')
          n = ichar(namold(4:4)) - 48
          if ( possig .gt. 5 ) then
             n = n*10 + ichar(namold(5:5)) - 48
          end if
          if ( n .le. naxis/2 ) then
             namcha = 'sigv'//char(48+n)//''''
          else
             namcha = 'sigv'//char(48+n-naxis/2)//'"'
          end if
       end if
    else if ( group(1:1).eq.'C' .and. mod(naxis,2).eq.0 &
         !...  groups Cnh, n even
         !...  generate labels for the two sets of vertical mirror planes,
         !...  sigvk` and sigvk``, C2v is handled explicitly.
       &          .and. index(group,'H') .gt. 0 ) then        ! Chv, n even
       call addinv (namold,namcha)
    else if ( group(1:2).eq.'OH' .or. group(1:2).eq.'TH' &
         &          .or. group(1:2).eq.'TD') then               ! Oh, Th, Td
       !...  groups Oh, Th, Td
       !...  generate labels for all improper rotations
       !...  by applying the inversion
       possig = index(namold,'''')
       if ( namold(1:2) .eq. 'C2' .and. possig .eq. 0 ) then
          !...  main C2 operations, C2x, C2y, C2z
          namcha = 'sig'//namold(3:3)
       else if ( namold(1:3) .eq. 'C3,' ) then
          namcha = 'S6,'//namold(4:4)
          if ( namold(5:5) .eq. '+' ) then
             namcha(5:5) = '-'
          else
             namcha(5:5) = '+'
          end if
       else if ( namold(1:2) .eq. 'C4' ) then
          namcha = 'S4'//namold(3:3)
          if ( namold(4:4) .eq. '+' ) then
             namcha(4:4) = '-'
          else
             namcha(4:4) = '+'
          end if
       else if ( (namold(1:2) .eq. 'C2') .and. (possig.ne. 0) ) then
          !...  dihedral C2 operations to give sigdk
          namcha = 'sigd'//char(ichar(namold(3:3))-48)
       else if ( namold(1:2) .eq. 'E ' ) then
          namcha = 'i'
       end if
    else if ( group(1:2).eq.'IH' ) then                   ! Ih
       !...  group Ih
       !...  generate labels for all improper rotations
       !...  by applying the inversion
       if ( namold(1:2) .eq. 'E ' ) then
          namcha = 'i'
       else if ( namold(1:2).eq.'C3' .or. namold(1:2).eq.'C5' ) then
          possig = index(namold,'+') + index(namold,'-')
          poscar = index(namold,'^')
          if ( poscar .eq. 0 ) then
             poskup = possig - 1
          else
             poskup = poscar - 1
          end if
          kvalue = namold(3:poskup)
          namcha = namold(1:2)//namold(poskup+1:)
          !
          ! Extra parens  around the intent(in)  argument are intended
          ! to  avoid illegal aliasing  of intent(in)  and intent(out)
          ! arguments:
          !
          call addinv ((namcha), namcha)
          !...  use numbering of proper rotation
          possig = index(namcha,'+') + index(namcha,'-')
          poscar = index(namcha,'^')
          if ( poscar .eq. 0 ) then
             poskup = possig - 1
          else
             poskup = poscar - 1
          end if
          namcha = namcha(1:poskup)//trim(kvalue)//namcha(poskup+1:)
       else if ( namold(1:2) .eq. 'C2' ) then
          !...  use numbering of proper rotation
          namcha = 'sig'//namold(3:3)
       end if
    else
       !...  default case, no change of name
       namcha = namold
    end if
    !
    namnew = namcha
    !
  end subroutine suplab


end subroutine gengrp
