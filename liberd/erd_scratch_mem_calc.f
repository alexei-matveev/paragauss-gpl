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
      subroutine erd_scratch_mem_calc(nshells, shtype, ncfps, npfps, 
     *                                atom, coords, alpha, pcoeff,
     *                                ixalpha, ixpcoeff, ccbeg, ccend,
     *                                indx_cc, spherical,
     *                                nuclei, calc_2der, imax, zmax,
     *                                me, nprocs)
c-------------------------------------------------------------------------
c   Calculates the maximum amount (in words) of scratch memory required
c   by the Flocke integral package.  The calculation is carried out in 
c   parallel over the number of processors in the job. 
c
c   Arguments:
c	nshells			Number of shells 
c       shtype			Array of shell types.
c	ncfps			Number of contracted functions per shell.
c	npfps			Number of primitive functions per shell.
c       atom                    Array of atomic labels (one for each shell).
c       coords                  Array of coords per shell.
c       calc_2der               Logical variable set .true. if we need the
c                               2nd derivative integrals in this job..  
c	imax			Returned number of integer words required.
c	zmax 			Returned number of fl. pt. words required.
c       me                      Rank of current processor.
c       nprocs                  Number of processors.
c--------------------------------------------------------------------------

      implicit none
      include 'machine_types.h'
      include 'dbugcom.h'

      integer nshells 
      integer shtype(nshells), ncfps(nshells), npfps(nshells)
      integer atom(nshells)
      integer ixalpha(nshells), ixpcoeff(nshells)
      integer ccbeg(*), ccend(*), indx_cc(*)
      integer imax, zmax
      integer nuclei, nalpha, ncoeff
      integer me, nprocs

      logical calc_2der

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
      integer der_flags(12), iflag, jflag 
      integer lmax, ncmax, npmax
      integer m1, m2, mrange, mleft

      logical spherical

      imax = 0
      zmax = 0

      mrange = nshells / nprocs
      mleft  = nshells - mrange * nprocs
      if (me .lt. mleft) then
         mrange = mrange + 1
         m1 = me * mrange + 1
      else
         m1 = me*mrange + mleft + 1
      endif  
      m2 = min(m1 + mrange - 1, nshells)
      if (m1 .gt. nshells) return   ! no work to do

      do m = m1, m2
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
     *                   ncfps(m), ncfps(n), ncfps(r), ncfps(s),
     *                   npfps(m), npfps(n), npfps(r), npfps(s),
     *                   shtype(m), shtype(n), shtype(r), shtype(s),
     *                   coords(1,m),coords(2,m),coords(3,m),
     *                   coords(1,n),coords(2,n),coords(3,n), 
     *                   coords(1,r),coords(2,r),coords(3,r), 
     *                   coords(1,s),coords(2,s),coords(3,s), 
     *                   alpha_pack, pcoeff_pack, spherical,
     *                   imin, iblk, zmin, zblk)
         imax = max0(imax, imin)
         zmax = max0(zmax, zmin)
      enddo


      do m = m1, m2
      do n = m, nshells
      do r = 1, nshells
      do s = r, nshells

c-------------------------------------------------------------------------
c   Derivative integrals.   Perform loops the same as they will be run
c   at execution time.
c-------------------------------------------------------------------------

c         do iflag = 1, 12
          do iflag = 1, 3

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
     *                 ncfps(m), ncfps(n), ncfps(r), ncfps(s),
     *                 npfps(m), npfps(n), npfps(r), npfps(s),
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
 
c----------------------------------------------------------------------------
c   Second-derivative integrals. Loop through all possible flag combinations.  
c----------------------------------------------------------------------------

         if (.not. calc_2der) go to 2000
         do iflag = 1, 12
         do jflag = iflag, 12

            do j = 1, 12
               der_flags(j) = 0
            enddo

            der_flags(iflag) = 1
            der_flags(jflag) = der_flags(jflag) + 1

               if (atom(m) .eq. atom(n)) then
                  if (der_flags(1) .ne. 0) der_flags(4) = der_flags(1)
                  if (der_flags(4) .ne. 0) der_flags(1) = der_flags(4)
                  if (der_flags(2) .ne. 0) der_flags(5) = der_flags(2)
                  if (der_flags(5) .ne. 0) der_flags(2) = der_flags(5)
                  if (der_flags(3) .ne. 0) der_flags(6) = der_flags(3)
                  if (der_flags(6) .ne. 0) der_flags(3) = der_flags(6)
               endif

               if (atom(m) .eq. atom(r)) then
                  if (der_flags(1) .ne. 0) der_flags(7) = der_flags(1)
                  if (der_flags(7) .ne. 0) der_flags(1) = der_flags(7)
                  if (der_flags(2) .ne. 0) der_flags(8) = der_flags(2)
                  if (der_flags(8) .ne. 0) der_flags(2) = der_flags(8)
                  if (der_flags(3) .ne. 0) der_flags(9) = der_flags(3)
                  if (der_flags(9) .ne. 0) der_flags(3) = der_flags(9)
               endif

               if (atom(m) .eq. atom(s)) then
                  if (der_flags(1) .ne. 0) der_flags(10) = der_flags(1)
                  if (der_flags(10).ne. 0) der_flags(1)  = der_flags(10)
                  if (der_flags(2) .ne. 0) der_flags(11) = der_flags(2)
                  if (der_flags(11).ne. 0) der_flags(2)  = der_flags(11)
                  if (der_flags(3 ).ne. 0) der_flags(12) = der_flags(3)
                  if (der_flags(12).ne. 0) der_flags(3)  = der_flags(12)
               endif

               if (atom(n) .eq. atom(r)) then
                  if (der_flags(4) .ne. 0) der_flags(7) = der_flags(4)
                  if (der_flags(7) .ne. 0) der_flags(4) = der_flags(7)
                  if (der_flags(5) .ne. 0) der_flags(8) = der_flags(5)
                  if (der_flags(8) .ne. 0) der_flags(5) = der_flags(8)
                  if (der_flags(6) .ne. 0) der_flags(9) = der_flags(6)
                  if (der_flags(9) .ne. 0) der_flags(6) = der_flags(9)
               endif

               if (atom(n) .eq. atom(s)) then
                  if (der_flags(4)  .ne. 0) der_flags(10)= der_flags(4)
                  if (der_flags(10) .ne. 0) der_flags(4) = der_flags(10)
                  if (der_flags(5)  .ne. 0) der_flags(11)= der_flags(5)
                  if (der_flags(11) .ne. 0) der_flags(5) = der_flags(11)
                  if (der_flags(6)  .ne. 0) der_flags(12)= der_flags(6)
                  if (der_flags(12) .ne. 0) der_flags(6) = der_flags(12)
               endif

               if (atom(r) .eq. atom(s)) then
                  if (der_flags(7)  .ne. 0) der_flags(10)= der_flags(7)
                  if (der_flags(10) .ne. 0) der_flags(7) = der_flags(10)
                  if (der_flags(8)  .ne. 0) der_flags(11)= der_flags(8)
                  if (der_flags(11) .ne. 0) der_flags(8) = der_flags(11)
                  if (der_flags(9)  .ne. 0) der_flags(12)= der_flags(9)
                  if (der_flags(12) .ne. 0) der_flags(9) = der_flags(12)
               endif 

c--------------------------------------------------------------------------
c   Set the derivative flag arguments.
c--------------------------------------------------------------------------
 
               call erd__memory_eri_derv_batch(nalpha, ncoeff, 
     *                 ncfps(m), ncfps(n), ncfps(r), ncfps(s),
     *                 npfps(m), npfps(n), npfps(r), npfps(s),
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
         enddo
 2000    continue

      enddo
      enddo
      enddo
      enddo

c------------------------------------------------------------------------
c   Now calculate the memory for the one-electron integrals.
c------------------------------------------------------------------------

      do m = m1, m2
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
