C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
      subroutine fast_erd_scratch_mem_calc(nshells, shtype,ncfps,npfps, 
     *                                atom, coords, alpha, pcoeff,
     *                                ixalpha, ixpcoeff, ccbeg, ccend,
     *                                indx_cc, spherical,
     *                                nuclei, imax, zmax)
c-------------------------------------------------------------------------
c   Calculates the maximum amount (in words) of scratch memory required
c   by the Flocke integral package.
c
c   Arguments:
c	nshells			Number of shells 
c       shtype			Array of shell types.
c	ncfps			Number of contracted functions per shell.
c	npfps			Number of primitive functions per shell.
c       atom                    Array of atomic labels (one for each shell).
c       coords                  Array of coords per shell.
c	imax			Returned number of integer words required.
c	zmax 			Returned number of fl. pt. words required.
c--------------------------------------------------------------------------

      implicit none
      include 'machine_types.h'

      integer nshells 
      integer shtype(nshells), ncfps(nshells), npfps(nshells)
      integer atom(nshells)
      integer ixalpha(nshells), ixpcoeff(nshells)
      integer ccbeg(*), ccend(*), indx_cc(*)
      integer imax, zmax
      integer nuclei, nalpha, ncoeff

      double precision coords(3,nshells)
      double precision alpha(*), pcoeff(*)

      integer max_dim_coeff
      parameter (max_dim_coeff = 5000)
      double precision alpha_pack(max_dim_coeff), 
     *                 pcoeff_pack(max_dim_coeff)
      integer ccbeg_pack(max_dim_coeff), ccend_pack(max_dim_coeff)

      integer i,j
      integer m, n, r, s
      integer iblk, zblk
      integer imin, zmin
      integer iopt, zopt
      integer ixderc
      integer der1x, der1y, der1z
      integer der2x, der2y, der2z
      integer dercx, dercy, dercz
      integer der_flags(12), iflag
      integer lmax, ncmax, npmax

      logical*8 l8true
      logical*8 l8var
      logical spherical
      logical*8 l8spherical

      imax = 0
      zmax = 0
      l8true = .true.

c---------------------------------------------------------------------------
c   Get largest shell values.
c---------------------------------------------------------------------------

      lmax = 0
      ncmax = 0
      npmax = 0
      m = 0
      do i = 1, nshells
         if (shtype(i) .gt. lmax) then
            lmax = shtype(i)
            m    = i
         endif
         ncmax = max(ncmax, ncfps(i))
         npmax = max(npmax, npfps(i))
      enddo

      n = m
      r = m
      s = m

         call pack_coeffs(alpha, ixalpha, pcoeff, ixpcoeff,
     *                        ncfps, npfps, m, n, r, s,
     *                       alpha_pack, nalpha, pcoeff_pack,
     *                       ncoeff, ccbeg, ccend, indx_cc,
     *                       ccbeg_pack, ccend_pack)

c-------------------------------------------------------------------------
c   ERD version 2
c-------------------------------------------------------------------------

            call ERD__MEMORY_ERI_BATCH(
     *                   nalpha, ncoeff,                 
     *                   ncmax, ncmax, ncmax, ncmax,
     *                   npmax, npmax, npmax, npmax,
     *                   shtype(m), shtype(n), shtype(r), shtype(s),
     *                   coords(1,m),coords(2,m),coords(3,m),
     *                   coords(1,n),coords(2,n),coords(3,n), 
     *                   coords(1,r),coords(2,r),coords(3,r), 
     *                   coords(1,s),coords(2,s),coords(3,s), 
     *                   alpha_pack, pcoeff_pack, spherical,
     *                   imin, iblk, zmin, zblk)
         imax = max0(imax, imin)
         zmax = max0(zmax, zmin)

c-------------------------------------------------------------------------
c   Derivative integrals.   Perform loops the same as they will be run
c   at execution time.
c-------------------------------------------------------------------------

         do iflag = 1, 12

            do j = 1, 12
               der_flags(j) = 0
            enddo
            der_flags(iflag) = 1

               if (atom(m) .eq. atom(n)) then
                  if (der_flags(1) .eq. 1) der_flags(4) = 1
                  if (der_flags(4) .eq. 1) der_flags(1) = 1
                  if (der_flags(2) .eq. 1) der_flags(5) = 1
                  if (der_flags(5) .eq. 1) der_flags(2) = 1
                  if (der_flags(3) .eq. 1) der_flags(6) = 1
                  if (der_flags(6) .eq. 1) der_flags(3) = 1
               endif
                                                                                
               if (atom(m) .eq. atom(r)) then
                  if (der_flags(1) .eq. 1) der_flags(7) = 1
                  if (der_flags(7) .eq. 1) der_flags(1) = 1
                  if (der_flags(2) .eq. 1) der_flags(8) = 1
                  if (der_flags(8) .eq. 1) der_flags(2) = 1
                  if (der_flags(3) .eq. 1) der_flags(9) = 1
                  if (der_flags(9) .eq. 1) der_flags(3) = 1
               endif
                                                                                
               if (atom(m) .eq. atom(s)) then
                  if (der_flags(1) .eq. 1) der_flags(10) = 1
                  if (der_flags(10).eq. 1) der_flags(1)  = 1
                  if (der_flags(2) .eq. 1) der_flags(11) = 1
                  if (der_flags(11).eq. 1) der_flags(2)  = 1
                  if (der_flags(3 ).eq. 1) der_flags(12) = 1
                  if (der_flags(12).eq. 1) der_flags(3)  = 1
               endif
                   
               if (atom(n) .eq. atom(r)) then
                  if (der_flags(4) .eq. 1) der_flags(7) = 1
                  if (der_flags(7) .eq. 1) der_flags(4) = 1
                  if (der_flags(5) .eq. 1) der_flags(8) = 1
                  if (der_flags(8) .eq. 1) der_flags(5) = 1
                  if (der_flags(6) .eq. 1) der_flags(9) = 1
                  if (der_flags(9) .eq. 1) der_flags(6) = 1
               endif

               if (atom(n) .eq. atom(s)) then
                  if (der_flags(4) .eq. 1) der_flags(10) = 1
                  if (der_flags(10) .eq. 1) der_flags(4) = 1
                  if (der_flags(5) .eq. 1) der_flags(11) = 1
                  if (der_flags(11) .eq. 1) der_flags(5) = 1
                  if (der_flags(6) .eq. 1) der_flags(12) = 1
                  if (der_flags(12) .eq. 1) der_flags(6) = 1
               endif
                   
               if (atom(r) .eq. atom(s)) then
                  if (der_flags(7) .eq. 1) der_flags(10) = 1
                  if (der_flags(10) .eq. 1) der_flags(7) = 1
                  if (der_flags(8) .eq. 1) der_flags(11) = 1
                  if (der_flags(11) .eq. 1) der_flags(8) = 1
                  if (der_flags(9) .eq. 1) der_flags(12) = 1
                  if (der_flags(12) .eq. 1) der_flags(9) = 1
               endif

c--------------------------------------------------------------------------
c   Set the derivative flag arguments.
c--------------------------------------------------------------------------
 
               call erd__memory_eri_derv_batch(nalpha, ncoeff, 
     *                 ncmax, ncmax, ncmax, ncmax,
     *                 npmax, npmax, npmax, npmax,
     *                 shtype(m), shtype(n), shtype(r), shtype(s),
     *                   coords(1,m),coords(2,m),coords(3,m),
     *                   coords(1,n),coords(2,n),coords(3,n),
     *                   coords(1,r),coords(2,r),coords(3,r),
     *                   coords(1,s),coords(2,s),coords(3,s),
     *                 der_flags(1), der_flags(2), der_flags(3),
     *                 der_flags(4), der_flags(5), der_flags(6),
     *                 der_flags(7), der_flags(8), der_flags(9),
     *                 der_flags(10),der_flags(11), der_flags(12),    
     *                   alpha_pack, pcoeff_pack, spherical,
     *                   imin, iblk, zmin, zblk)   

            imax = max0(imax, imin)
            zmax = max0(zmax, zmin)
         enddo
     
c------------------------------------------------------------------------
c   Now calculate the memory for the one-electron integrals.
c------------------------------------------------------------------------

      do m = 1, nshells
      do n = 1, nshells
         call pack_coeffs_oed(alpha, ixalpha, pcoeff, ixpcoeff,
     *                        ncfps, npfps, m, n, 
     *                       alpha_pack, nalpha, pcoeff_pack,
     *                       ncoeff, ccbeg, ccend, indx_cc,
     *                       ccbeg_pack, ccend_pack,
     *                       max_dim_coeff)
            call oed__memory_nai_batch(nalpha, ncoeff,
     *                   ncfps(m), ncfps(n), 
     *                   npfps(m), npfps(n), 
     *                   shtype(m), shtype(n),
     *                   coords(1,m),coords(2,m),coords(3,m),
     *                   coords(1,n),coords(2,n),coords(3,n),
     *                   nuclei, alpha_pack, pcoeff_pack, spherical, 
     *                   imin, iopt, zmin, zopt)

         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

            call oed__memory_kin_batch(nalpha, ncoeff,
     *                    ncfps(m), ncfps(n), 
     *                    npfps(m), npfps(n), 
     *                    shtype(m), shtype(n),
     *                    coords(1,m),coords(2,m),coords(3,m),
     *                    coords(1,n),coords(2,n),coords(3,n),
     *                    alpha_pack, pcoeff_pack, spherical, 
     *                    imin, iopt, zmin, zopt)

         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

            call oed__memory_ovl_batch(nalpha, ncoeff,
     *                      ncfps(m), ncfps(n), 
     *                      npfps(m), npfps(n), 
     *                      shtype(m), shtype(n),
     *                      coords(1,m),coords(2,m),coords(3,m),
     *                      coords(1,n),coords(2,n),coords(3,n),
     *                      alpha_pack, pcoeff_pack, spherical, 
     *                      imin, iopt, zmin, zopt)

         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

         go to 1000
 
            call oed__memory_nai_derv_batch(nalpha, ncoeff,
     *                   ncfps(m), ncfps(n), 
     *                   npfps(m), npfps(n), 
     *                   shtype(m), shtype(n),
     *                   coords(1,m),coords(2,m),coords(3,m),
     *                   coords(1,n),coords(2,n),coords(3,n),
     *                   nuclei, ixderc, 
     *                   der1x, der1y, der1z,
     *                   der2x, der2y, der2z,
     *                   dercx, dercy, dercz,
     *                   alpha_pack, pcoeff_pack, spherical, 
     *                   imin, iopt, zmin, zopt)

         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

            call oed__memory_kin_derv_batch(nalpha, ncoeff,
     *                   ncfps(m), ncfps(n), 
     *                   npfps(m), npfps(n), 
     *                   shtype(m), shtype(n),
     *                   coords(1,m),coords(2,m),coords(3,m),
     *                   coords(1,n),coords(2,n),coords(3,n),
     *                   der1x, der1y, der1z,
     *                   der2x, der2y, der2z,
     *                   alpha_pack, pcoeff_pack, spherical, 
     *                   imin, iopt, zmin, zopt)

         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

            call oed__memory_ovl_derv_batch(nalpha, ncoeff,
     *                   ncfps(m), ncfps(n), 
     *                   npfps(m), npfps(n), 
     *                   shtype(m), shtype(n),
     *                   coords(1,m),coords(2,m),coords(3,m),
     *                   coords(1,n),coords(2,n),coords(3,n),
     *                   der1x, der1y, der1z,
     *                   der2x, der2y, der2z,
     *                   alpha_pack, pcoeff_pack, spherical, 
     *                   imin, iopt, zmin, zopt)

         imax = max0(iopt, imax)
         zmax = max0(zopt, zmax)

 1000    continue
      enddo
      enddo

      return
      end
