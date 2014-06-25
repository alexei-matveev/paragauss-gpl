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
      subroutine freq(kolat,kolvar,Bmat,Fmat,muext,G,Fstor,VR,n_00)
	use smlcom
	use zarsml
	use zzzsml
	use slspar_module, only: n_userdefined_atoms, userdefined_atom
      implicit real*8(a-h,o-z)
!	include 'comsim'
      real*8
     .   Bmat(3*kolat,kolvar),Fmat(kolvar,kolvar)
     .  ,G(kolvar,kolvar),amass(99),VR(kolvar,kolvar)
     >	,VIN(nlocv),F1(nlocv),F2(nlocv),Fstor(kolvar,kolvar)
     .  ,AMass1(nlocv),Freqr(nlocv),Fcoef
      integer   ifq1(nlocv), ifq2(nlocv),ipv(nlocv,3)

      integer   kolvar                 ! number of internal variables
      integer   kolat                  ! number of atoms
      integer   ierr                   ! error code returned by RS

c	amu=1.6605655*10.d-27
c	me= 0.9109534*10.d-30
c	2Ry=27.211652
c	1eV= 8065.465 cm-1
c	amu/me=1822.8874
c	5140.49=27.211652*8065.465/sqrt(1822.8874)

      DATA      Fcoef/5143.05D0/,ZERO/0.d0/ ! 7       8 
c      data      amass/1.007825, 4.002603, 7.016005, 9.012182,11.009305
      data      amass/1.007825, 2.000000, 7.016005, 9.012182,11.009305
     >              ,12.000000,14.003074,15.994915,18.998403,19.992439
     >              ,22.989770,23.985045,26.981541,27.976928,30.973763
     >              ,31.973909,34.968853,39.962383,38.963708,39.962591
     >              ,44.955914,47.947947,50.943963,51.940510,54.938046
     >              ,55.934939,58.933198,57.935347,62.929599,63.929145
     >              ,68.925581,73.921179,74.921596,79.916521,78.918336
     >              ,83.911506,84.911800,87.956250,88.905856,89.904708
     >             ,92.906378,97.905405,97.907110,101.904348,102.90550
     >          ,106.903480,106.905095,113.903361,114.90388,119.902199
     >          ,120.903824,129.906230,126.904477,131.90415,132.905770
     >          ,137.90524 ,138.90636 ,139.90544 ,140.90766,141.90773
     >          ,144.912691,151.91974 ,152.92124 ,157.92411,158.92535
     >          ,163.92918 ,164.93033 ,167.93031 ,168.93423,173.93887
     >          ,174.94079 ,177.94371 ,180.94801 ,183.95095,186.95577
     >          ,191.96149 ,192.96294 ,194.96479 ,196.96656,201.97063
     >          ,204.97441 ,207.97664 ,208.98039 ,208.98242,219.01130
     >          ,222.01757 ,223.019734,226.025406,227.02775,232.03805
     >          ,231.035881,238.050786,237.048169,244.06420,4*250,10000/
	amasj=1.d0
	do i=1,kolat
         ia=an(i)+0.01
	 if(ieq(i).lt.0) amasj=amasj+amass(ia)
	enddo
	print*,'total mass of molecule n_00',amasj,n_00
	
	
1     k=0
      do i=1,kolat
	if (i.eq.n_00) then
	 do jj=1,3 
	   k=k+1
	   AMass1(k)=1.d0/amasj
	 enddo
	else
         ia=an(i)+0.01
          do j=1,3
           k=k+1
           AMass1(k)=1.d0/amass(ia)
          enddo
	endif
      enddo
	if(n_userdefined_atoms.ne.0) then
	print*, 'correction with user defined atoms'
	 k=0
      do i=1,kolat
	do j=1,3
	k=k+1 
	do jj=1,n_userdefined_atoms
	if(i.eq.userdefined_atom(jj)%no) AMass1(k)=1.0_r8_kind/userdefined_atom(jj)%mass
	enddo
	enddo
	enddo
	endif
CAY get G=B*M^(-1)*Bt
      do i=1,muext   ! kolvar
        do j=1,muext
          G(i,j)=ZERO
          do k=1,3*kolat
            G(i,j)=G(i,j)+AMass1(k)*Bmat(k,i)*Bmat(k,j)
          enddo
        enddo ! j=1,kolvar
      enddo
	print*, 'mass_matt(:2,:2)'
	print*, g(1,1),g(1,2)
	print*, g(2,1),g(2,2)
	print*,'muext,kolvar', muext, kolvar
	call matinv(G,muext,det,kolvar,ipv)
        call matinv(G,kolvar,det,kolvar,ipv)
	print*, g(1,1),g(1,2)
	print*, g(2,1),g(2,2)
c calculation of square root of G-matrix
	MATZ=1
	call RS(kolvar,kolvar,G,Freqr,MATZ,VR,F1,F2,IERR)
	do ii = 1,kolvar
		do jj = 1,kolvar
			Fstor(ii,jj) = ZERO
			if (ii .eq. jj) Fstor(ii,jj) = sqrt(Freqr(ii))
		enddo
	enddo
	do i = 1,kolvar
		do j = 1,kolvar
			G(i,j) = ZERO
			do k = 1,kolvar
				do l = 1,kolvar
				G(i,j) = G(i,j) + VR(i,l)*Fstor(l,k)*VR(j,k)
				enddo
			enddo
		enddo
	enddo
c diagonalization of G**(1/2)FmatG**(1/2)-matrix
	do i1 = 1,kolvar
		do j1 = 1,kolvar
			Fstor(i1,j1) = ZERO
			do k1 = 1,kolvar
				do l1 = 1,kolvar
			         Fstor(i1,j1) = Fstor(i1,j1)   
     >	+		G(i1,l1)*Fmat(l1,k1)*G(k1,j1)
				enddo
			enddo
		enddo
	enddo
	call RS(kolvar,kolvar,Fstor,Freqr,MATZ,VR,F1,F2,IERR)
c standardization of of eigenvectors
	do i2 = 1,kolvar
		do j2 = 1,kolvar
			Fstor(i2,j2) = ZERO
			do k2 = 1,kolvar
			Fstor(i2,j2) = Fstor(i2,j2) + G(i2,k2)*VR(k2,j2)
			enddo
		enddo
	enddo
c calculation of intensity of vibrations
c
	call intens(VIN,Fstor,kolvar)
c
c identification of vibrations ( maximum value of eigenvectors - normal coordinates)
c
	do j=1,kolvar
		amaxi = ZERO
		ifq1(j) = 1
		do j1=1,kolvar
			if ( abs(Fstor(j1,j)) .gt. amaxi ) then
				amaxi = abs(Fstor(j1,j))
				ifq1(j) = j1
			endif
		enddo
	enddo
	write (*,*) 'eigenvectors'
	do k2=1,kolvar
		write (*,102) (Fstor(k1,k2),k1=1,kolvar)
	enddo
102	format (100f10.4)

c additional identification of vibrations (distribution of potential 
c energy on internal coordinates in %)
	do i=1,kolvar
		e = 0.0d0
		do i1=1,kolvar
			VR(i1,i)= fmat(i1,i1)*Fstor(i1,i)**2
			e = e+VR(i1,i)
		enddo
		do i2=1,kolvar
			VR(i2,i)= VR(i2,i)/e*100.d0
		enddo
	enddo
	write (*,*) 'distribution of potential energy in %'
	do k2=1,kolvar
		write (*,102) (VR(k1,k2),k1=1,kolvar)
	enddo
	do jj=1,kolvar
		amaxi = ZERO
		ifq2(jj) = 1
		do jj1=1,kolvar
			if ( abs(VR(jj1,jj)) .gt. amaxi ) then
				amaxi = abs(VR(jj1,jj))
				ifq2(jj) = jj1
			endif
		enddo
	enddo


	write(*,*)' Frequencies( cm**(-1) ), identification and'
	write(*,*)' intensities(m/mole)'
	do i=1,kolvar
		if(Freqr(i).gt.0.d0) then 
        		write(*,101)sqrt(Freqr(i))*Fcoef,ifq1(i),ifq2(i)
     *				,VIN(i), Freqr(i) 
		else
			write(*,101)Freqr(i)*Fcoef,ifq1(i),ifq2(i)
     *				,VIN(i)
		endif
	enddo ! i=1,kolvar
	if(n_00.ne.0) then
		n_00=0
		goto 1
	endif
101   FORMAT('  ',F15.4,1X,2I3,1X,F15.4,f15.6)
      e n d
c
	subroutine intens (VIN,VR1,kolvar)
	use smlcom, only: nlocv, sxo	
	implicit real*8 (A-H,O-Z)
c
c	It is the program of calculation of arbitrary IR intensities 
c	of fundamental vibrations of molecules
c
	real*8 VIN(kolvar),dpx(nlocv,2),dpy(nlocv,2),dpz(nlocv,2)
     &        ,VR1(kolvar,kolvar)
	real*8 dpqx(nlocv),dpqy(nlocv),dpqz(nlocv),dpsx(nlocv)
     &        ,dpsy(nlocv),dpsz(nlocv),Nav,
     *	PI,C
	logical ex
 	common /gen/igener(nlocv)
! 	common /sfreq/ sxo(300),/gen/igener(300)
	data PI/3.14159265358979324D+00/, C/2.997925D+08/,Nav/6.02205D+23/
c
	inquire (file='dipmom.sav',exist=ex)
        if (.not. ex) goto 1000
	open (40, file='dipmom.sav')
	do i=1,kolvar
		read (40,1) dpx(i,1),dpy(i,1),dpz(i,1)	! dipmom at x+delta
		read (40,1) dpx(i,2),dpy(i,2),dpz(i,2)	! dipmom at x-delta
	enddo
	close(40)
1	format (3f15.8)
c
c dPx,y,z/dS
	do j=1,kolvar
		dpsx(j)=(dpx(j,1)-dpx(j,2))/(2.d0*sxo(j)) ! /igener(j) 
		dpsy(j)=(dpy(j,1)-dpy(j,2))/(2.d0*sxo(j)) ! /igener(j) 
		dpsz(j)=(dpz(j,1)-dpz(j,2))/(2.d0*sxo(j)) ! /igener(j) 
	enddo
c
c dPx,y,z/dQ
	do k=1,kolvar
		dpqx(k)= 0.d0
		dpqy(k)= 0.d0
		dpqz(k)= 0.d0
c	write(*,*) k
		do k1=1,kolvar
			dpqx(k)=dpqx(k)+dpsx(k1)*VR1(k1,k)
			dpqy(k)=dpqy(k)+dpsy(k1)*VR1(k1,k)
			dpqz(k)=dpqz(k)+dpsz(k1)*VR1(k1,k)
c	write(*,*) dpsz(k1),VR1(k1,k),dpsz(k1)*VR1(k1,k)
		enddo
	enddo
c
c calculation of I~(dPx/dQ**2+dPy/dQ**2+dPz/dQ**2)
	write (*,*) '|dM/dQ|[(m**(3/2))/s or C/(kg**(1/2))]'
	do i=1,kolvar
		VIN(i)=(dpqx(i)**2 + dpqy(i)**2 + dpqz(i)**2)  
		write (*,*) sqrt(VIN(i))*3.7274d-1  
		VIN(i)= VIN(i)*3.7274d-1**2
		VIN(i)=(VIN(i)*PI*Nav/(3*C**2))
	enddo
1000	return
	end
