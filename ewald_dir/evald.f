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
	module ewald_main
	use ew_grid_module, only: r8_kind
        use sar_module, only:ndqt
	 real (r8_kind), dimension(ndqt)::an_epe_type =0.0_r8_kind
	 real (r8_kind), dimension(:),allocatable::an_epe
	end module ewald_main
	use scr_module		! 	the use statements go first
	use sar_module
	use LATT_module
	use ew_grid_module
	use ewald_main

        use mesh_module
        use proj_rest_module

      IMPLICIT REAL*8 (A-H,O-Z)    
	external system
      DIMENSION ABCENT(3,200)

	type(vertex)::latt_pr(1)
	real(kind=r8_kind),dimension(3):: latt_shift_pr=0.0_r8_kind
     &  ,lsxx=(/1.,0.,0./),lsyy=(/0.,1.,0.0/),lszz=(/0.,0.,1./)
	 real(r8_kind),dimension(3):: dipole	! dipole moment of the lattice
	real*8,  dimension(3)	::	dsh
	real*8,  dimension(3)	::	xx,yy,zz
	real*8,  pointer:: a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z
        real*8   vrest                   ! surf. and restored potential
        real*8, pointer                     :: vs(:)
        real*8, allocatable, dimension(:,:) :: c

          
      INTEGER*4,dimension(400):: CENTRS
	character*4 un,UNITS(5)
      DATA UNITS/'ANGS','A.U.','ATOM','LATT','AHDS'/
	DATA EPS/1.D0/,AHDS/1.d0/
	
	logical lgopt
        logical lgxcell
        logical lcellv 
        logical Surf_PM

	namelist /rrr/xx, yy, zz, dsh, lgopt	!hades lcgto interface shifts

      NAMELIST /PARMS/ 
     &                A,NATOMS,NPOINT,CENTRS,EPS,KSD,KSR,ZATOM,KLAT
     &		      ,xx, yy, zz, dsh
     &		      ,lhds,AHDS, lgxcell, lcellv, np_ions,npd
     &                ,pcmerge ,c1,c2, make_pcmerge   
     &                ,action_multipole
     &                ,action_screep,screep_xyz,r_screep
     &                ,action_screep_c3
     &                ,action_screep_c4
     &	              ,xmol_forms	
     &                ,latt_shift_pr,lsxx,lsyy,lszz
     &		      ,an_epe_type,NZ_screep
     &                ,lmax,ispher          ! multipole only
      NAMELIST /unit/ UN,NATOMS	! units in which lattice is defined
      NAMELIST /cluster/ klat 

 
	CENTRS=0
	lhds=.false.
	lgxcell=.true.
	make_pcmerge=.false.
	pcmerge=.false.
	lcellv = .false.

	a1x=>a(1,1)
	a1y=>a(2,1)
	A1Z=>A(3,1)
	A2X=>A(1,2)
	A2Y=>A(2,2)
	A2Z=>A(3,2)
	A3X=>A(1,3)
	A3Y=>A(2,3)
	A3Z=>A(3,3)

	ZATOM(1:20)=0.D0

      KSD=3
      KSR=2
      NPOINT=0
      READ 1000,UN                                         
	write(6,1000) UN
1000  FORMAT(A4)                                          
      DO 1 IU=1,5                                        
         IF (UN .EQ. UNITS(IU)) GOTO 2                  
   1  CONTINUE                                         
      PRINT 9000                                      
9000  FORMAT ('Wrong lenth units'  )  
      STOP 9000                                               
   2  READ (5,PARMS)
      WRITE(6,PARMS)
        
        if(action_screep .and. action_multipole) then
      write(*,*)  
     *'======================================================='
      write(*,*) 
     *' Both action_screep and action_multipole are choosen ! '
      write(*,*) 
     *'======================================================='
        stop
        end if

	if(action_screep_c3) action_screep=action_screep_c3
	if(make_pcmerge) pcmerge=.true.

        if(lgxcell .and. lcellv) then
           print*,'Either LGXCELL or LCELLVEC can be TRUE no both'
           stop
        endif

	if(lgxcell) inquire (file='gxcell',exist=lgxcell)
        if(lcellv) inquire (file='cellvec',exist=lcellv)
	if(lhds) then
	read(5,rrr)
	write(6,rrr) 
	endif ! lhds
	rr(:,1)=xx
	rr(:,2)=yy
	rr(:,3)=zz
	ds=dsh/angsau

      NPT=NPOINT
      IF(NATOMS.gt.400) then
      PRINT 117
117   FORMAT(10X,' NUMBER OF ATOMS MUST BE LESS THEN 400 ')
      STOP
	endif ! (NATOMS.gt.200
  17  PRINT 2001,UN
2001  FORMAT ('  TRANSLATION VECTORS OF DIRECT LATTACE IN ',A4//
     *        10X,'X',10X,'Y',10X,'Z')
      PRINT 2002,A
2002  FORMAT ((3(5X,F10.5)))

      IF (IU .eq. 1) then
	a=a/ANGSAU
	elseif(iu.eq.5) then
          A=A*AHDS
	endif    ! IU .eq. 1

	if(lhds) then	! use Hades input of cell
	READ(5,unit)
	write(6,unit)
	else
        READ 1000,UN                                     
	endif	! lhds/else

      DO 5 iiu=1,5                                     
         IF (UN .EQ. UNITS(iiu)) then
	iu=iiu
	GOTO 6               
	endif	! UN .EQ. UNITS(IU)
   5  CONTINUE                                      
      PRINT 9001                                   
9001  FORMAT ('   wrong specification of units ')
      STOP 9000                                             

    6 continue
	READ (5,*) (NA(JAT),IQT(JAT),SYMVL(JAT),
     *             RATOM(:,JAT),JAT=1,NATOMS)
1006  FORMAT (I3,I2,1X,A4,3F12.7)

	print*,'lattice shifted with latt_shift_pr and dipole moment'
	print*,        ' of the lattice'
	dipole=0.0_r8_kind
	   rr_x =lsxx
	   rr_y =lsyy
	   rr_z =lszz 
	do jat=1,NATOMS
	   latt_pr(1)%r=RATOM(:,JAT)+latt_shift_pr
           call scrot_vertex(latt_pr,1)	   
	 print 1006, NA(JAT),IQT(JAT),SYMVL(JAT),latt_pr(1)
	
	dipole=dipole+RATOM(:,JAT)*zatom(IQT(JAT))
	enddo ! i=1,NATOMS
	print*,'dipole moment of the translation cell'
	print*,dipole

	
	if(lgxcell .or. lcellv) then
           if (lgxcell) then 
            open(76, file = 'gxcell',  form='formatted', 
     *         status='unknown')
                 call read_cell_v
           endif
           if(lcellv) then
              open(76, file = 'cellvec',  form='formatted', 
     *        status='unknown')
              call read_cell_p
           endif
           close(76)
 	un='ANGS'
	if(lgxcell)write(6,*) ' gx.c coordinates in atomic units'
	if(lcellv)write(6,*) ' gx.cv coordinates in atomic units'
	do i=1,natoms
	   write(6,'(4d15.8,i3,4i2)') RATOM(:,i)/angsau,0.0,i,0,0,0,0
	enddo	! i=1,natoms

	dipole=0.0_r8_kind
	do jat=1,NATOMS
	   latt_pr(1)%r=RATOM(:,JAT)+latt_shift_pr
           call scrot_vertex(latt_pr,1)	   
	 print 1006, NA(JAT),IQT(JAT),SYMVL(JAT),latt_pr(1)
	
	dipole=dipole+RATOM(:,JAT)*zatom(IQT(JAT))
	enddo ! i=1,NATOMS
	endif ! lgxcell
	print*, 'Ortep input'
	do jat=1,NATOMS
	  print '(A6,3X,18x,3F9.4,i9)',SYMVL(JAT),RATOM(:,JAT),2
	  print '(1x,F8.2,45x,7X,i2)',0.20,7
	enddo
        
	print *,'dipole moments of the translation  cell'
	print*, dipole



      IF(NPOINT.NE.0)then	! points with explicitely defined coordinates
	do JP=1,NPOINT
	 READ(5,*) RP(:,JP)
	enddo	! JP=1,NPOINT
	endif	! NPOINT.NE.0


      NN=1
1500  NC=CENTRS(NN)
      IF(NC.EQ.0) GO TO 1515
	  RP(:,NPOINT+NN)=RATOM(:,NC)
      NN=NN+1
      GO TO 1500
1515  NPOINT=NPOINT+NN-1

        
        if(action_screep) then    
           if (NZ_screep.EQ.32) then
               write(6,*) 'screep sphere coordinates in input units'
	       write(6,*) screep_xyz%r(:)
	       call first_div ! make first division of sphere
	       db%sec_to_rp=reshape(sec,(/32/))
	       ind_scr=NPOINT+1
	       do i=1,32
	       rp(:,NPOINT+i)=db(i)%sec_to_rp%r
	       enddo ! i=1,32
	       NPOINT=NPOINT+32
           endif	
     
           if(NZ_screep.eq.128) then
               write(6,*) 'screep sphere coordinates in input units'
               write(6,*) screep_xyz%r(:)
               call first_second_div ! make first and second  division of sphere
               db_128%sec_to_rp=reshape(sec128,(/128/))
               ind_scr=NPOINT+1
               do i=1,128
               rp(:,NPOINT+i)=db_128(i)%sec_to_rp%r
               enddo ! i=1,128
               NPOINT=NPOINT+128
           endif 
       end if ! action_screep

       if(action_multipole .or. action_screep_c4) then         
        write(6,*) 'genmesh'
        write(6,*) 'screep_xyz%r',screep_xyz%r,
     *              'r_screep', r_screep, 'ispher', ispher
               if(action_multipole)  allocate( c(lmax+1,2*lmax+1) )
               call genmesh(screep_xyz%r,r_screep,ispher)
               ind_scr=NPOINT+1
               do i=1,mesh%npoints
               rp(1,NPOINT+i)=mesh%xp(i)
               rp(2,NPOINT+i)=mesh%yp(i)
               rp(3,NPOINT+i)=mesh%zp(i)
                  if(action_screep_c4) then
                    db_128(i)%sec_to_rp%r(1)=mesh%xp(i)
                    db_128(i)%sec_to_rp%r(2)=mesh%yp(i)
                    db_128(i)%sec_to_rp%r(3)=mesh%zp(i)
                  end if ! action_screep_c4
               end do
               NPOINT=NPOINT+mesh%npoints
               NZ_screep=mesh%npoints
       end if ! action_multipole .or. action_screep_c4
       

      PRINT 118,NPOINT
 118  FORMAT(10X,'NPOINT=',I5)
      IF(NPOINT.LE.400) GO TO 119
      PRINT 120
  120 FORMAT(10X,'NUMBER OF NPOINTS MUST BE LESS TNEN 400')
      STOP
  119 NQT=0
      DO 19 I=1,NATOMS
      IF(IQT(I).GT.NQT) NQT=IQT(I)
  19  CONTINUE
      IF(NQT.LE.40) GO TO 219
      PRINT 220,NQT
 220  FORMAT(10X,'NUMBER OF TYPES MUST WE LESS TNEN 40 NQT=',I5)
      STOP
1007  FORMAT (3G10.4)                                      
 219  PRINT 2061,UN                                       
2061  FORMAT (/'  TYPES  AND COORD. OF CELL ATOMS IN ',A4//
     *        2X,'NN',5X,'CHARGE',9X,'X(1)',13X,'Y(2)',13X,'Z(3)')
      PRINT 2062,(JAT,IQT(JAT),SYMVL(JAT)
     *,(RATOM(I,JAT),I=1,3),JAT=1,NATOMS)
2062  FORMAT ((2X,I3,2X,I2,2X,A4,3(2X,G10.4)))

      PRINT 2071,UN
2071  FORMAT (/'  COORD. OF POINTS IN ',A4//
     *        4X,'NN',11X,'X(1)',13X,'Y(2)',13X,'Z(3)')
      PRINT 2072,(JP,(RP(I,JP),I=1,3),JP=1,NPOINT)       
2072  FORMAT ((2X,I5,3X,3(2X,G15.7)))

	select case(iu)
	case(1)
!	angstrems
          RATOM(:,1:NATOMS)=RATOM(:,1:NATOMS)/ANGSAU       
          RP(:,1:NPOINT)=RP(:,1:NPOINT)/ANGSAU           
          screep_xyz%r=screep_xyz%r/ANGSAU       
            if(action_multipole) then 
            r_screep=r_screep/ANGSAU
            mesh%xp=mesh%xp/ANGSAU
            mesh%yp=mesh%yp/ANGSAU
            mesh%zp=mesh%zp/ANGSAU
            end if    
	case(5)
!	convert from AHDS units
          RATOM(:,1:NATOMS)=RATOM(:,1:NATOMS)*AHDS
          RP(:,1:NPOINT)=RP(:,1:NPOINT)*AHDS
          screep_xyz%r=screep_xyz%r*AHDS
	case(4)
	 write(6,*)' Lattice coordinates transformed in A.U.'
      DO  JAT=1,NATOMS                                 
         RL1=RATOM(1,JAT)                               
         RL2=RATOM(2,JAT)                              
         RL3=RATOM(3,JAT)                             
         RATOM(1,JAT)=RL1*A1X+RL2*A2X+RL3*A3X        
         RATOM(2,JAT)=RL1*A1Y+RL2*A2Y+RL3*A3Y       
         RATOM(3,JAT)=RL1*A1Z+RL2*A2Z+RL3*A3Z      
	write(6,*) symvl(jat),  ratom(:,jat)*angsau
	enddo ! JAT=1,NATOMS

      DO  JP=1,NPOINT                          
         RL1=RP(1,JP)                           
         RL2=RP(2,JP)                          
         RL3=RP(3,JP)                         
         RP(1,JP)=RL1*A1X+RL2*A2X+RL3*A3X    
         RP(2,JP)=RL1*A1Y+RL2*A2Y+RL3*A3Y   
         RP(3,JP)=RL1*A1Z+RL2*A2Z+RL3*A3Z  
	enddo ! JP=1,NPOINT
         RL1=screep_xyz%r(1)
         RL2=screep_xyz%r(2)
         RL3=screep_xyz%r(3)
	screep_xyz%r(1)=RL1*A1X+RL2*A2X+RL3*A3X
	screep_xyz%r(2)=RL1*A1Y+RL2*A2Y+RL3*A3Y
	screep_xyz%r(3)=RL1*A1Z+RL2*A2Z+RL3*A3Z
	end select


      if(action_screep) then   ! point charge screep sphere
         do i=1,NZ_screep
         if (NZ_screep.eq.32) db(i)%sec_to_rp%r=rp(:,ind_scr-1+i)
         if (NZ_screep.eq.128) db_128(i)%sec_to_rp%r=rp(:,ind_scr-1+i)
         enddo ! i=1,32 or 1,128
      end if   
       
      if(action_screep_c4) then
         do i=1,NZ_screep
         db_128(i)%sec_to_rp%r=rp(:,ind_scr-1+i)
         enddo ! i=1,32 or 1,128
      end if              


      PRINT 2066                                       
2066  FORMAT (/'  TYPES  AND COORD. OF CELL ATOMS IN A.U.'//
     *        2X,'NN',5X,'CHARGE',9X,'X(1)',13X,'Y(2)',13X,'Z(3)')
      PRINT 2067,(JAT,IQT(JAT),SYMVL(JAT)
     *,RATOM(1:3,JAT),JAT=1,NATOMS)
2067  FORMAT ((2X,I3,2X,I2,2X,A4,3(2X,G15.8)))

      CALL EVALD(EPS)

	write(6,*) '  specifiacation of cluster, type of charge	'
	write(6,*) '  coordinates and type of the periodic potential'
	write(6,*)
      PRINT110                           
  110 FORMAT(22X,'MQ',6X,'X',11X,'Y',11X,'Z',7X,'MT'/)

      DO I=1,KLAT                     
      READ (5,*)  MQ(I),ABCENT(:,I),MT(I)
!      READ (5,*)  MT(I)
	if( action_screep .or. action_multipole .or. action_screep_c4 ) then
	 ABCENT(:,I)= rp(:,MT(I))
	endif ! action_screep -- action_multipole -- action_screep_c4
  202 FORMAT(3X,I3,4X,3F10.5,4X,I3)
	  PRINT130,MQ(I),ABCENT(:,I),MT(I)
	enddo ! I=1,KLAT

  130 FORMAT(20X,I3,4X,3F10.5,4X,I3)
      INQT=NQT
      DO  I=1,KLAT
	  IF(MQ(I).GT.NQT) NQT=MQ(I)
	enddo ! I=1,KLAT

      DO 140 I=1,KLAT                  
      DO  IQ=1,NQT
	   VM(I,IQ)=0.D0
	enddo !	IQ=1,NQT
      DO 140 J=1,KLAT                 
      rw=0.D0
      DO  K=1,3                   
	  rw=rw+(ABCENT(K,I)-ABCENT(K,J))**2
	enddo ! K=1,3
      R=DSQRT(rw)
      IQ=MQ(J)
      IF(I.NE.J.AND.R.LT.1.D-4) PRINT 137,I,ABCENT(:,I),
     *J,ABCENT(:,J)
      GO TO 136
  137 FORMAT(10X,'ERROR IN ATOM',I3,2X,3F10.5/15X,'ATOM',I3,2X,3F10.5)
  136 IF(I.NE.J) VM(I,IQ)=VM(I,IQ)+1.D0/R
  140 CONTINUE

      DO I=1,KLAT
      II=MT(I)
      DO  J=1,NQT
	      IF(J.GT.INQT) PSUM(J,II)=0.D0
		VM(I,J)=PSUM(J,II)-VM(I,J)
	enddo ! J=1,NQT
	enddo ! I=1,KLAT
      ZSUM=0.D0
      DO I=1,NATOMS
	      IQ=IQT(I)
		ZSUM=ZSUM+ZATOM(IQ)
	enddo ! I=1,NATOMS
      PRINT 18,ZSUM
   18 FORMAT(10X,'totel charge of cell is equal to ',F7.3)
      IF(DABS(ZSUM).GT.0.001D0) STOP
      PRINT 2130
      DO JP=1,NPOINT            
	Psuma(JP)=0.d0
      PMADL(JP)=0.D0
	PMADL(JP)=dot_product(PSUM(1:INQT,JP),ZATOM(1:INQT))
 	 Psuma(JP)= dot_product(PSUMnp(1:INQT,JP),ZATOM(1:INQT))
      NT=0
      IF(JP.GT.NPT) NT=CENTRS(JP-NPT)
      PRINT 2131,JP,NT,PMADL(JP),Psuma(JP),(PSUM(I,JP),I=1,NQT),vpcm(jp)
	enddo ! JP=1,NPOINT
2130  FORMAT ( /'   madelung potential and its partial contributions'/
     &         '   on given points in a.u.'     /
     &         '   in the forth column for finite PC arrays potentials '/
     &         '   are presented'/
     *  4X,'NN', 2X,'NT',3X,'MAD.PT.    PC Arr    partial contributions'
     &     ' MPC ' /)
2131  FORMAT(2X,I3,1X,I3,10G12.5)
2132  PRINT 2137
2137  FORMAT ( /' MADELUNG POTENTIAL ON GIVEN ATOMS  IN A.U.'/)

!     point charge screep sphere    
        if(action_screep .and. NZ_screep.eq.32) then
	if(pcmerge) then
	db%vel=PMADL(ind_scr:ind_scr+31)-vpcm(ind_scr:ind_scr+31)
	print*,'db%vel'
	do i=ind_scr,ind_scr+31
	print*,db(i)%vel
	enddo
	! substract field from merged array
	else
	db%vel=PMADL(ind_scr:ind_scr+31)-Psuma(ind_scr:ind_scr+31)
	endif ! else

        do i=1,NZ_screep
	db_a(i,i)=0.560*sqrt(NZ_screep*1.0D0)/r_screep
	if(i.eq.NZ_screep) exit
	do j=i+1,NZ_screep
	xc%r=db(i)%sec_to_rp%r-db(j)%sec_to_rp%r	
	db_a(i,j)=1.d0/sqrt(dot_product(xc%r,xc%r))
	db_a(j,i)=db_a(i,j)
	enddo
	enddo

	call matinv(db_a,NZ_screep,det,NZ_screep,db_ipv)

	db%q=matmul(db_a,db%vel)
	db_qtot=sum(db%q)
	write(6,*) 'total screep charge', db_qtot
	write(6,*) 'screep surface potentials'
!	check solutions
	xm%r=0.d0
	do i=1,NZ_screep
	xm%r= xm%r+db(i)%sec_to_rp%r/NZ_screep
	db(i)%vch=0.d0
	do j=1,NZ_screep
	if(i.eq.j) then
        db(i)%vch=db(i)%vch+
     *   db(j)%q*0.560*sqrt(NZ_screep*1.0D0)/r_screep
	else
	xc%r=db(i)%sec_to_rp%r-db(j)%sec_to_rp%r
	db(i)%vch=db(i)%vch+db(j)%q/sqrt(dot_product(xc%r,xc%r))
	endif

	enddo ! j=1,32
	enddo ! i=1,32
	
	write(6,*) 'position of the screep sphere center in a.u'
	write(6,*) xm%r
	write(6,*)
        write(6,*) 'exect and fitted potentials, SCREEP charges'
	do i=1,NZ_screep
	 write(6,'(i5,3f16.6)') i,db(i)%vel,db(i)%vch,db(i)%q
	enddo ! i=1,32

	endif ! action_screep

	if(action_screep .and. NZ_screep.eq.128) then
	if(pcmerge) then
        db_128%vel=PMADL(ind_scr:ind_scr+127)-
     * vpcm(ind_scr:ind_scr+127)
	else
        db_128%vel=PMADL(ind_scr:ind_scr+127)-
     *   Psuma(ind_scr:ind_scr+127)
	endif ! else
	do i=1,NZ_screep
	db_a_128(i,i)=0.560*sqrt(NZ_screep*1.0D0)/r_screep
	if(i.eq.NZ_screep) exit
	do j=i+1,NZ_screep
        xc%r=db_128(i)%sec_to_rp%r-db_128(j)%sec_to_rp%r
	db_a_128(i,j)=1.d0/sqrt(dot_product(xc%r,xc%r))
	db_a_128(j,i)=db_a_128(i,j)
	enddo
	enddo
        call matinv(db_a_128,NZ_screep,det,NZ_screep,db_ipv_128)
	db_128%q=matmul(db_a_128,db_128%vel)
	db_qtot=sum(db_128%q)
	write(6,*) 'total screep charge', db_qtot
	write(6,*) 'screep surface potentials'
!	check solutions
	xm%r=0.d0
	do i=1,NZ_screep
	xm%r= xm%r+db_128(i)%sec_to_rp%r/NZ_screep
	db_128(i)%vch=0.d0
	do j=1,NZ_screep
	if(i.eq.j) then
       db_128(i)%vch=
     *db_128(i)%vch+db_128(j)%q*0.560*sqrt(NZ_screep*1.0D0)/r_screep
	else
	xc%r=db_128(i)%sec_to_rp%r-db_128(j)%sec_to_rp%r
       db_128(i)%vch=
     *  db_128(i)%vch+db_128(j)%q/sqrt(dot_product(xc%r,xc%r))
	endif

	enddo ! j=1,128
	enddo ! i=1,128
	
	write(6,*) 'position of the screep sphere center in a.u'
	write(6,*) xm%r
	write(6,*)
        write(6,*) 'exect and fitted potentials, SCREEP charges'
	do i=1,NZ_screep
	 write(6,'(i5,3f16.6)')i,db_128(i)%vel,db_128(i)%vch,db_128(i)%q
	enddo ! i=1,128 

        end if ! action_screep
!
! = = = = = = = = = = = = = = = = = = SCREEP-sphere of C4 symmetry = = = = = = = = = = =
! variant action_screep_c4 for NZ_screpp .le. 128 (actually 38, 50, 110 points of
! Lebedev's grid generated  by < genmesh > subroutine 
!
	if( action_screep_c4 ) then
	if(pcmerge) then
        db_128%vel=PMADL(ind_scr:ind_scr+NZ_screep-1)-
     * vpcm(ind_scr:ind_scr+NZ_screep-1)
	else
        db_128%vel=PMADL(ind_scr:ind_scr+NZ_screep-1)-
     *   Psuma(ind_scr:ind_scr+NZ_screep-1)
	endif ! else
	do i=1,NZ_screep
	db_a_128(i,i)=0.560*sqrt(NZ_screep*1.0D0)/r_screep
	if(i.eq.NZ_screep) exit
	do j=i+1,NZ_screep
        xc%r=db_128(i)%sec_to_rp%r-db_128(j)%sec_to_rp%r
	db_a_128(i,j)=1.d0/sqrt(dot_product(xc%r,xc%r))
	db_a_128(j,i)=db_a_128(i,j)
	enddo
	enddo
        print*, 'call matinv'
        call matinv(db_a_128,NZ_screep,det,128,db_ipv_128)
        print*, 'out matinv'
	db_128%q=matmul(db_a_128,db_128%vel)
	db_qtot=sum(db_128%q)
	write(6,*) 'total screep charge', db_qtot
	write(6,*) 'screep surface potentials'
!	check solutions
	xm%r=0.d0
	do i=1,NZ_screep
	xm%r= xm%r+db_128(i)%sec_to_rp%r/NZ_screep
	db_128(i)%vch=0.d0
	do j=1,NZ_screep
	if(i.eq.j) then
       db_128(i)%vch=
     *db_128(i)%vch+db_128(j)%q*0.560*sqrt(NZ_screep*1.0D0)/r_screep
	else
	xc%r=db_128(i)%sec_to_rp%r-db_128(j)%sec_to_rp%r
       db_128(i)%vch=
     *  db_128(i)%vch+db_128(j)%q/sqrt(dot_product(xc%r,xc%r))
	endif

	enddo ! j=1,128
	enddo ! i=1,128
	
	write(6,*) 'position of the screep sphere center in a.u'
	write(6,*) xm%r
	write(6,*)
        write(6,*)'exact  and fitted potentials, SCREEP charges'
	do i=1,NZ_screep
	 write(6,'(i5,3f16.6)') i,db_128(i)%vel,db_128(i)%vch,db_128(i)%q
	enddo ! i=1,128 

        end if ! action_screep
!
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


      if(action_multipole) then     !     multipole SCREEP sphere  
        vs => db_128(:mesh%npoints)%vel
        if(pcmerge) then
        vs=PMADL(ind_scr:ind_scr+mesh%npoints-1)  
     *   - vpcm(ind_scr:ind_scr+mesh%npoints-1)
        else
        vs=PMADL(ind_scr:ind_scr+mesh%npoints-1)  
     *   - Psuma(ind_scr:ind_scr+mesh%npoints-1)
        end if

        call project(vs,lmax,screep_xyz%r,c)
!            check solutions
        write(6,*) 'exect and fitted potentials (SCREEP multipoles)'
        do i=1,mesh%npoints
        call rest(mesh%xp(i),mesh%yp(i),mesh%zp(i),
     *            r_screep,screep_xyz%r,lmax,c,vrest )
        write(6,'(i3,3x,f10.7,2x,f10.7)') i, vs(i), vrest 
        end do 
      end if  ! action_multipole

!	now estimate ouside potentials of cluster

	if(action_screep) print *, 'Exact and SCREEP potentials'
      DO I=1,KLAT

      if (action_screep. or.action_screep_c4) then

	vel(i)=0.d0
	do j=1,NZ_screep

        if (NZ_screep.eq.32 .and. .not.action_screep_c4) then
	xc%r=abcent(:,i)-db(j)%sec_to_rp%r
	vel(i)= db(j)%q/sqrt(dot_product(xc%r,xc%r))+vel(i)
         endif	

        if (NZ_screep.eq.128 .and. .not.action_screep_c4) then
        xc%r=abcent(:,i)-db_128(j)%sec_to_rp%r
        vel(i)= db_128(j)%q/sqrt(dot_product(xc%r,xc%r))+vel(i)
        endif

        if (action_screep_c4) then
        xc%r=abcent(:,i)-db_128(j)%sec_to_rp%r
        vel(i)= db_128(j)%q/sqrt(dot_product(xc%r,xc%r))+vel(i)
        endif
	enddo ! j=1,32,128
      end if ! action_screep_c4

      if(action_multipole) then
           call rest(abcent(1,i),abcent(2,i),abcent(3,i),
     *     r_screep,screep_xyz%r,lmax,c,vel(i) )             
      end if    ! action_multipole

      VMS(I)=0.D0
	if(action_screep .or. action_multipole .or.action_screep_c4) then
	if(pcmerge) then
       PRINT 4020,I,pmadl(MT(I)),vpcm(MT(I))+vel(i)
	print *,'vpcm', vpcm(MT(I)),vel(i)
	else
       PRINT 4020,I,pmadl(MT(I)),Psuma(MT(I))+vel(i)
	endif

	else
	   VMS(I)=VMS(I)+dot_product(VM(I,1:NQT),ZATOM(1:NQT))
        PRINT 4020,I,VMS(I),(VM(I,IQ),IQ=1,NQT)
	endif
 4020 FORMAT(2X,I3,10G14.7) ! old G12.5
	enddo ! I=1,KLAT
  212 format(4f15.8,I3,4I2,f10.2) ! LCGTO inpt format

!     ******  PC array
	inquire (file='evald.mpcr',exist=other_pc_arr)
	if(.not.other_pc_arr.and.make_pcmerge) then
	print*,'file evald.mpcr is created'
	ierr=system('cp evald.pcr evald.mpcr')
	endif

      if(action_screep .or. action_screep_c4) then
	close(10)
	if(pcmerge) then
	open(10,file='evald.mpcr',
     &  form='formatted',status='old')

	else ! not pcmerge
	open(10,file='evald.pcr',
     &  form='formatted',status='old')
	allocate (rrpm(4,numion+NZ_screep),an_epe(numion),stat=istat)
	if(istat.ne.0) stop 'allocation failed'
	endif !pcmerge /else

	print*, 'write merged and SCREEP charges '
	read (10,*) numion
	if(istat.ne.0) stop 'an_epe allocation failed'
	print*,numion,'numion'
	do i=1,numion
	read(10,*) rrpm(1:4,i),i1,i2,i3,i4,i5, an_epe(i)
	enddo ! i=1,numion
	close(10)
	open(10,file='evald.pcr', form='formatted',status='old')
	write(10,*) numion+NZ_screep
	do i=1,numion
	write(10,212) rrpm(1:4,i),1,0,0,1,0,an_epe(i)
	enddo ! i=1,numion
	do j=1,NZ_screep	! screep PC array

        if (NZ_screep.eq.32) then
         print 212,db(j)%sec_to_rp%r,db(j)%q,1,0,0,1,0
         write(10,212) db(j)%sec_to_rp%r,db(j)%q,1,0,0,1,0	,0.0_r8_kind
        endif

        if (NZ_screep.eq.128) then
         print 212,db_128(j)%sec_to_rp%r,db_128(j)%q,1,0,0,1,0
         write(10,212) db_128(j)%sec_to_rp%r,
     &        db_128(j)%q,1,0,0,1,0,0.0_r8_kind	
        endif

        if (action_screep_c4) then
         print 212,db_128(j)%sec_to_rp%r,db_128(j)%q,1,0,0,1,0
         write(10,212) 
     &   db_128(j)%sec_to_rp%r,db_128(j)%q,1,0,0,1,0,0.0_r8_kind	
        endif        

	enddo ! j=1,32,128 or action_screep_c4
	close(10)

!       deallocate(symvl,IQT,MT,MQ,NA,RATOM,RP,PSUM,PSUMnp,VM)
!       deallocate(VMS,PMADL,PSUMA,vel,SUMA,SUMB,sump,rrpm)
      end if  ! action_screep

      if(action_multipole) then
	close(10)
	if(pcmerge) then
	open(10,file='evald.mpcr',
     &  form='formatted',status='old')

	else ! not pcmerge
	open(10,file='evald.pcr',
     &  form='formatted',status='old')
	allocate (rrpm(4,numion+NZ_screep),an_epe(numion),stat=istat)
	if(istat.ne.0) stop 'allocation failed'
	endif !pcmerge /else

	print*, 'write merged and SCREEP charges '
	read (10,*) numion
	if(istat.ne.0) stop 'an_epe allocation failed'
	print*,numion,'numion'
	do i=1,numion
	read(10,*) rrpm(1:4,i),i1,i2,i3,i4,i5, an_epe(i)
	enddo ! i=1,numion
	close(10)

	open(10,file='evald.pcr', form='formatted',status='old')
	write(10,*) numion+NZ_screep
	do i=1,numion
	write(10,212) rrpm(1:4,i),1,0,0,1,0,an_epe(i)
	enddo ! i=1,numion
	
!                 Print Multipole coefficients
        write(10,'(/t10, "Multipole coefficients" /)' )
        lm=0
        do l = 1, lmax+1
        do m = 1,2*l-1
        lm = lm + 1 
        write(10,'( t8,i2,t14,i2,t19,i2,t27,f19.15,3x,a18)')  
     &        lm,l,m,vlm(lm)
        end do  
        end do ! do_lm
        write(10,'(/t10," It is not for ParaGauss use as yet "/)' )
	close(10)
      end if   ! action multipole

      STOP

      contains
	subroutine read_gx_file
	real*8, dimension(ngxat):: an
	real*8, dimension(3,ngxat):: gxcent
	integer, dimension(8,ngxat)::indgx
	integer, dimension(ngxat)::impu
	ii=0
2000	ii=ii+1
	if(ii.gt.ngxat) then
	write(73,*) 'ngxat is too small'
	else	! ii.le.ngxat
	read(76,*)an(ii),gxcent(:,ii),indgx(:,ii),impu(ii)
	if(an(ii).gt.0.5d0) go to 2000
	write(73,*) 'read_gx_file : No of centers in gx.c file',ii-1
	jj=0
	do i=1,ii-1
	if(indgx(1,i).ne.0) then
	jj=jj+1
	RATOM(:,jj)=gxcent(:,i)	!!!/ ANGSAU           
	endif	! indgx(1,i).ne.0
	enddo	!i=1,ii-1 
	endif	! ii.gt.kamax/else
	return
	end subroutine read_gx_file

        subroutine read_cell_v
          real :: R_SHELL(ngxat,3), R_CORE(ngxat,3),an(ngxat)
          real :: ann, r_sh(3)
          logical :: type_number(ngxat,ngxat)
          integer :: i,j

          RATOM=0d0; R_SHELL=0d0; R_CORE=0d0
          do I=1,NATOMS
             read(76,'(f5.2,3f15.7)') an(i),r_core(i,:)
             iqt(i)=int(an(i))
             type_number(iqt(i),i)=.true.
          enddo
	  
          do
             read(76,*,end=22)ann,r_sh
             if(ann == zero) goto 22
             i=ann
             do j=1,ngxat
                if(type_number(i,j))then
                   R_SHELL(j,:)=r_sh
                   type_number(i,j)=.false.
                   goto 21
                endif
             enddo 
 21	     continue
          enddo 
 22	  continue
          write (6,*) '=============================================='
          write (6,*) 'Unit cell has been taken from GULPv'
          write (6,*) '     New coordinates (angstrom)'
          do i=1,natoms
             if(dot_product(R_SHELL(i,:),R_SHELL(i,:)) /= zero) then 
                ratom(:,i)=(r_core(i,:)+r_shell(i,:))/2.0_r8_kind
             else
                ratom(:,i)=r_core(i,:)
             endif
             write(6,1007) i,iqt(i),symvl(i),ratom(:,i)
          enddo
          write (6,*) '=============================================='
1007      FORMAT (1x,I3,6x,I2,4X,a5,3F10.5)
        end subroutine read_cell_v

        subroutine read_cell_p
          real :: R_SHELL(ngxat,3), R_CORE(ngxat,3),an(ngxat)
          real :: ann, r_sh(3)
          logical :: type_number(ngxat,ngxat)
          integer :: i,j

          write (6,*) '=============================================='
          write (6,*) 'Unit cell has been taken from GULPp'
          read(76,*) ((a(i,j),i=1,3),j=1,3)
          write (6,*) '     New cell vectors'
          write(6,'(3f15.7)') ((a(i,j),i=1,3),j=1,3)
          write(6,*) '----------------------------------------------'

	  IF (IU .eq. 1) then
	     a=a/ANGSAU
	  elseif(iu.eq.5) then
	     A=A*AHDS
	  endif			! IU .eq. 1

          RATOM=0d0; R_SHELL=0d0; R_CORE=0d0
          do I=1,NATOMS
             read(76,'(f5.2,3f15.7)') an(i),r_core(i,:)
             iqt(i)=int(an(i))
             type_number(iqt(i),i)=.true.
          enddo
          do
             read(76,*,end=22)ann,r_sh
             if(ann == zero) goto 22
             i=ann
             do j=1,ngxat
                if(type_number(i,j))then
                   R_SHELL(j,:)=r_sh
                   type_number(i,j)=.false.
                   goto 21
                endif
             enddo 
 21	     continue
          enddo 
 22	  continue
          write (6,*) '     New coordinates (angstrom)'
          do i=1,natoms
             if(dot_product(R_SHELL(i,:),R_SHELL(i,:)) /= zero) then 
                ratom(:,i)=(r_core(i,:)+r_shell(i,:))/2.0_r8_kind
             else
                ratom(:,i)=r_core(i,:)
             endif
             write(6,1008) i,iqt(i),symvl(i),ratom(:,i)
          enddo
          write (6,*) '=============================================='
1008      FORMAT (1x,I3,6x,I2,4X,a5,3F10.5)
        end subroutine read_cell_p
      end

      SUBROUTINE EVALD (EPS)

	use scr_module
	use sar_module, only:	IQT,RATOM,RP,PSUM,zatom,vpcm
     &                 ,suma,sumb,sump,psumnp,ndvt,SYMVL
	use latt_module ! numion=0
	use ew_grid_module, only: screep_xyz,r_screep,xc
	use ewald_main
      IMPLICIT REAL*8 (A-H,O-Z)                                 
	external system
	real*8, dimension(3) ::  rrp
	real*8 :: totzatom=0.d0	! 
	integer,dimension(3)::	npdir	! direct space counts
	real*8,dimension(ndvt):: sumd
	real*8, pointer:: A1X,A1Y,A1Z,A2X,A2Y,A2Z,A3X,A3Y,A3Z
	real*8, pointer:: b1X,b1Y,b1Z,b2X,b2Y,b2Z,b3X,b3Y,b3Z
	real*8, target:: B(3,3)
      DIMENSION BA(3),BB(3),BC(3),BL(3),AA(3,3),RRL(3)
      parameter ( PI=3.141592653D0 )                             
        if(np_ions.ne.0.or.lhds) then
        open(10, file='evald.pcr', form='formatted',status='unknown')
	endif

	A1X=>A(1,1)
	A1Y=>A(2,1)
	A1Z=>A(3,1)
	A2X=>A(1,2)
	A2Y=>A(2,2)
	A2Z=>A(3,2)
	A3X=>A(1,3)
	A3Y=>A(2,3)
	A3Z=>A(3,3)

	b1X=>b(1,1)
	b1Y=>b(2,1)
	b1Z=>b(3,1)
	b2X=>b(1,2)
	b2Y=>b(2,2)
	b2Z=>b(3,2)
	b3X=>b(1,3)
	b3Y=>b(2,3)
	b3Z=>b(3,3)

      B1X=A2Y*A3Z-A2Z*A3Y
      B1Y=A2Z*A3X-A2X*A3Z                            
      B1Z=A2X*A3Y-A2Y*A3X                           
      B2X=A3Y*A1Z-A3Z*A1Y                          
      B2Y=A3Z*A1X-A3X*A1Z                         
      B2Z=A3X*A1Y-A3Y*A1X                        
      B3X=A1Y*A2Z-A1Z*A2Y                              
      B3Y=A1Z*A2X-A1X*A2Z                             
      B3Z=A1X*A2Y-A1Y*A2X                            
      VCELL=ABS(A1X*B1X+A1Y*B1Y+A1Z*B1Z)           
	b=b/VCELL
	write(*,*) 'translation vectors of reciprocial cell'
	write(6,'((3(5X,F10.5)))') b 
      AA1=A1X*A1X+A1Y*A1Y+A1Z*A1Z                  
      SQAA1=DSQRT(AA1)                            
      AA2=A2X*A2X+A2Y*A2Y+A2Z*A2Z                
      SQAA2=DSQRT(AA2)                          
      AA3=A3X*A3X+A3Y*A3Y+A3Z*A3Z              
      SQAA3=DSQRT(AA3)                        
      BB1=B1X*B1X+B1Y*B1Y+B1Z*B1Z            
      SQBB1=DSQRT(BB1)                      
      BB2=B2X*B2X+B2Y*B2Y+B2Z*B2Z          
      SQBB2=DSQRT(BB2)                    
      BB3=B3X*B3X+B3Y*B3Y+B3Z*B3Z        
      SQBB3=DSQRT(BB3)                  
                                                                
      PI2=2.D0*PI                                                    
      SQPI=DSQRT(PI)                                                
      THR=1.D0/3.D0                                                
      TVCELL=VCELL**THR                                           
      EPSILO=SQPI/TVCELL*EPS                                     
	write(6,*) 'EPSILO:',EPSILO
      CBACK=2.D0/(PI*VCELL)                                     
      COEFF=-(PI/EPSILO)**2
      ACONST=PI/EPSILO**2/VCELL                                
	write(*,*) 'ACONST coeff cback ', ACONST, COEFF, CBACK
      APOL=2.D0*EPSILO/SQPI

      RS=IABS(KSD)*TVCELL                             
      BS=IABS(KSR)/TVCELL                            
      BBS=BS*BS                                     
      PRINT 2000,KSD,KSR,VCELL,rs
2000  FORMAT ('direct and back space tranlation paremeters '
     &        ,2I5,'DIR. AND REV.',
     *   /10X,'VCELL=',G12.4
     &  ,/10x,'DIR SPACE RAD',G12.4)

      IF(KSR.LT.0) GO TO 2222
      MAXB1=IDINT(BS*SQAA1+0.1D0)+1
      MAXB2=IDINT(BS*SQAA2+0.1D0)*2+1
      MAXB3=IDINT(BS*SQAA3+0.1D0)*2+1
      GO TO 3333
 2222 MAXB1=IABS(KSR)+1
      MAXB2=IABS(KSR)*2+1
      MAXB3=IABS(KSR)*2+1
 3333 IF(KSD.gT.0) then
      MAXA1=IDINT(RS*SQBB1+0.1D0)*2+1
      MAXA2=IDINT(RS*SQBB2+0.1D0)*2+1
      MAXA3=IDINT(RS*SQBB3+0.1D0)*2+1
      else
 4444 MAXA1=IABS(KSD)*2+1
      MAXA2=IABS(KSD)*2+1
      MAXA3=IABS(KSD)*2+1
	endif	! KSD.gT.0/else
C	if(npd(1,1).eq.0.and.npd(2,1).eq.0) then
C	npd(1,1)=MAXA1
C	npd(2,1)=-MAXA1
C	endif
C	if(npd(1,2).eq.0.and.npd(2,2).eq.0) then
C	npd(1,2)=MAXA2
C	npd(2,2)=-MAXA2
C	endif
C	if(npd(1,3).eq.0.and.npd(2,3).eq.0) then
C	npd(1,3)=MAXA3
C	npd(2,3)=-MAXA3
C	endif
 1414 PRINT 1415, MAXA1,MAXA2,MAXA3,MAXB1,MAXB2,MAXB3         
 1415 FORMAT(/10X,'MAXA(3),MAXB(3) '/10X,6I5)
      DO 110 JP=1,NPOINT                                     
	vpcm(jp)=0.d0	! potential due to merged arr
	sumd(1:NQT)=0.D0
	      SUMA(1:NQT)=0.D0
	      SUMB(1:NQT)=0.D0
		sump(1:NQT)=0.d0
      RX=RP(1,JP)
      RY=RP(2,JP)
      RZ=RP(3,JP)
      DO  K=1,3
	      BC(K)=0.D0
	      BA(K)=-B(K,1)
	enddo ! K=1,3
      KT=0
      IS=2
      ISS=2
!	sumbtt=0.d0
      DO 50 L1=1,MAXB1
      DO K=1,3
	      BB(K)=0.D0
	      BA(K)=BA(K)+B(K,1)
	enddo ! K=1,3
      IF(L1.EQ.2) IS=1
      IND=1
      DO L2=1,MAXB2,IS
      NID=IND*L2
      INDD=1
      IF(L2.EQ.3) ISS=1
      DO  L3=1,MAXB3,ISS
      KT=KT+1
      NIDD=INDD*L3
      BBL=0.D0
      DO K=1,3
	b3=Bc(3)/B(3,3)
	      BL(K)=BA(K)+BB(K)+BC(K)
	      BC(K)=BC(K)+NIDD*B(K,3)
	      BBL=BBL+BL(K)*BL(K)
	enddo	! K=1,3
      IF(KT.EQ.1) GO TO 56

C     IF (BBL.GT.BBS) GOTO 55
	b1=BA(1)/B(1,1)
	b2=Bb(2)/B(2,2)
!	sumbt=0.
      DO  I=1,NATOMS
	      RIX=RX-RATOM(1,I)
	      RIY=RY-RATOM(2,I)
	      RIZ=RZ-RATOM(3,I)
	      IQ=IQT(I)
	      BETA=BL(1)*RIX+BL(2)*RIY+BL(3)*RIZ
	      BETA=PI2*(BETA-IDINT(BETA))
	      SUMB(IQ)=SUMB(IQ)+EXP(COEFF*BBL)*COS(BETA)/BBL
!	sumbt=sumbt+EXP(COEFF*BBL)*COS(BETA)/BBL*zatom(iq)*cback
!	sumbtt=sumbtt+EXP(COEFF*BBL)*COS(BETA)/BBL*zatom(iq)*cback
!	if(jp.eq.1.and.b1.le.1.5) write(*,*) 'bet ',i,BETA
!     *               ,EXP(COEFF*BBL)*COS(BETA)/BBL*zatom(iq)*cback
	enddo ! I=1,NATOMS
!	if(jp.eq.1) then
!	write(*,*) b1,b2,b3,BBL/0.529177**2
!	write(*,*) 'sumbt ', sumbt,sumbtt
!	endif	! jp.eq.1
  56  INDD=-INDD
	enddo	!L3=1,MAXB3,ISS 
	BC=0.d0
	BB=BB+NID*B(:,2)
  55  IND=-IND
	enddo ! 	L2=1,MAXB2,IS
  50     CONTINUE                                       
             
         IND=1
	nn=0
	nnn=0
          AA=0.D0
	npdir=0
      DO  L1=1,MAXA1
      NID=IND*L1
      INDD=1
      DO L2=1,MAXA2
      NIDD=INDD*L2
      INDDD=1
      DO  L3=1,MAXA3
      NIDDD=INDDD*L3
	
      DO  K=1,3
         RRL(K)=AA(K,1)+AA(K,2)+AA(K,3)
         AA(K,3)=AA(K,3)+A(K,3)*NIDDD
	enddo !  K=1,3
	npdir(3)=npdir(3)+NIDDD
	nnn=nnn+1		! count number of translations
      DO 67 I=1,NATOMS
      IQ=IQT(I)
      RIX=RX-RATOM(1,I)
      RIY=RY-RATOM(2,I)
      RIZ=RZ-RATOM(3,I)
      RL=DSQRT((RRL(1)-RIX)**2+(RRL(2)-RIY)**2+(RRL(3)-RIZ)**2)

	if(lhds.or.np_ions.ne.0) then
	rrp=rrl+ratom(:,i)
!!!	if(numion.le.10) write(*,*) 'r ',rrp(1:3)
	call scrot_ev(1,rrp,1)	! convert to system of lcgto cluster
  64  format (i5,3x,3f10.5,3x,f10.5)

        totzatom=totzatom+zatom(iq)
      if(numion.le.np_ions) then
	if(npdir(3)-niddd.le.npd(1,3).and.npdir(3)-niddd.ge.npd(2,3)
     &         .and.npdir(2).le.npd(1,2).and.npdir(2).ge.npd(2,2)
     &      .and .npdir(1).le.npd(1,1).and .npdir(1).ge.npd(2,1))then
	if(rl.gt.1.D-4) SUMp(iq)=SUMp(iq)+1.d0/rl 
	if(jp.eq.1) then ! make output only while processing first point
	xc%r=screep_xyz%r-(rrl+ratom(:,i))
	if(dot_product(xc%r,xc%r).le.(r_screep/0.529177)**2) then
c$$$	write(6,*) SYMVL(i),(rrl+ratom(:,i))*0.529177d0!!!, zatom(iq),i
	endif ! dot_product(xc%r,xc%r).le.(screep_r/0.529177)**2
        numion=numion+1
!!! 	if(i.eq.1) write(6,*) npdir(1),npdir(2),npdir(3)-niddd
  202 format(4f15.8,I3,4I2,f10.2) ! LCGTO inpt format
	if(xmol_form) then
	write(10,*) symvl(abs(int(zatom(iq)))),RRp(1),RRp(2),RRp(3)
	else
	write(10,202) RRp(1:3),zatom(iq),1,0,0,1,0,an_epe_type(iq)
	endif ! xmol_form/else
	endif   ! jp.eq.1
	endif	!  abs(AA(K,3)...
	endif 	! (numion.le.112
	endif ! lhds
                            
C     IF (RL.GT.RS) GOTO 67
      EPSIRL=EPSILO*RL
       IF (RL .ge. 1.D-4) then
      SUMA(IQ)=SUMA(IQ)+ERFC(EPSIRL)/RL
      SUMd(IQ)=SUMd(IQ)+ERFC(EPSIRL)/RL
 	else
      SUMA(IQ)=SUMA(IQ)-APOL
!      SUMd(IQ)=SUMd(IQ)-APOL
 	endif
!	if(jp.eq.2.and.rl.lt.5.0) then
!	write(6,*)  ERFC(EPSIRL)/RL,rl
!	endif
  67  CONTINUE
	  INDDD=-INDDD
	enddo ! L3=1,MAXA3
      DO  K=1,3
	      AA(K,2)=AA(K,2)+A(K,2)*NIDD
	      AA(K,3)=0.D0
	enddo ! K=1,3
	npdir(3)=0
	npdir(2)=npdir(2)+NIDD
	  INDD=-INDD
	enddo !  L2=1,MAXA2
      DO  K=1,3
        AA(K,1)=AA(K,1)+A(K,1)*NID
        AA(K,2)=0.D0
      enddo ! K=1,3
	npdir(2)=0
	npdir(1)=npdir(1)+NID
  60  IND=-IND
      enddo	! L1=1,MAXA1
!      write(6,*) JP,(SUMd(IQ),IQ=1,NQT),epsilo
      DO  N=1,NATOMS
       SUMA(IQT(N))=SUMA(IQT(N))-ACONST
	enddo ! N=1,NATOMS
!      PRINT 2610,JP,(SUMA(IQ),IQ=1,NQT)
       DO  IQ=1,NQT
        PSUM(IQ,JP)=SUMB(IQ)*CBACK+SUMA(IQ)
	PSUMnp(IQ,JP)=sump(iq)
      enddo	!  IQ=1,NQT	 
!      write(*,*) JP,(cback*SUMb(IQ),IQ=1,NQT)
      write(6, 2610) JP,(PSUMnp(IQ,JP),IQ=1,NQT)
2610  FORMAT(5X,I5,10F12.7)
	if(jp.eq.1) then
	call atch_number
	if(pcmerge) call merge_pcarr()
	endif	! jp.eq.1
	if(pcmerge.and.other_pc_arr) then
	call calc_vpcm()
	endif ! pcmerge
 110  CONTINUE                                                    
	contains
	subroutine calc_vpcm()
!	write(*,*) 'pcm coordinates'
	do i=1,nopc_arr
	rrp(:)=rrpm(1:3,i)
	call scrot_ev(1,rrp,-1)	! convert to system of lcgto cluster
!	if(i.le.10) write(*,*) rrp(1:3)
	rl=sqrt((rx-rrp(1))**2+(ry-rrp(2))**2+(rz-rrp(3))**2)
	if(rl.gt.4.d-2) vpcm(jp)=vpcm(jp)+rrpm(4,i)/rl
	enddo ! i=1,nopc_arr
	end subroutine calc_vpcm

	subroutine atch_number()
	external system
	integer*4 system
	if(np_ions.ne.0.or.lhds) close (10)
	if(numion.ne.0) then	! add number of ions at start of file
	open(10, file='evald.npcr', form='formatted',status='unknown')
	write(10,*) numion 
	close(10)
	isy=system('cat evald.pcr>>evald.npcr;mv evald.npcr  evald.pcr')
	write(6,*)
	write(6,*) ' number of ions in the finite PC array', numion
	endif	! numion.ne.0
	end subroutine atch_number

	subroutine merge_pcarr()
!	makes a statistic mearge of several arrays
	real*8, dimension(3):: rdiff
	real *8 zat	! charge of current center in evald.pcr
	integer*4, dimension(5):: deck

	inquire (file='evald.mpcr',exist=other_pc_arr)
	write(*,*) 'other_pc_arr', other_pc_arr

	if(other_pc_arr) then
	write(6,*) 'file evald.mpcr is located'
	open(10,file='evald.mpcr',form='formatted',status='unknown')
	read (10,*) nopc_arr	! numb of ions in other_pc_arr
	allocate (rrpm(4,nopc_arr+numion)
     &                ,an_epe(nopc_arr+numion),stat=istat)
	if(istat.ne.0) stop 'allocation failed'
	do i=1,nopc_arr
	read(10,*) rrpm(1:4,i),deck,an_epe(i)
	enddo 	! nopc_arr
	close(10)
	if(.not.make_pcmerge)  return

	open(10,file='evald.pcr',form='formatted',status='unknown')
	read (10,*) numion
	rrpm(4,:nopc_arr)=rrpm(4,:nopc_arr)*c2
	outer:	do i=1,numion
	 read(10,*) RRp(1),RRp(2),RRp(3),zat,deck,an_epe(i)
	 inner: do j=1,nopc_arr
	  rdiff=rrp-rrpm(:3,j)
	  r2=dot_product(rdiff,rdiff)
	if(r2.lt.4.d-2) then
	rrpm(4,j)=rrpm(4,j)+zat*c1
	cycle outer
	endif ! r2.lt.1.d-6 /else
	 enddo inner
	nopc_arr=nopc_arr+1
	rrpm(:3,nopc_arr)=rrp
	rrpm(4,nopc_arr)=zat*c1
        an_epe(nopc_arr)=an_epe(i)
	enddo outer
	rrpm(4,:nopc_arr)=rrpm(4,:nopc_arr)/(c1+c2)

	close(10)
	open(10,file='evald.mpcr',form='formatted',status='unknown')
	write(10,*)nopc_arr
  202 format(4f15.8,I3,4I2,f10.2) ! LCGTO inpt format
	do i=1,nopc_arr
	write(10,202) rrpm(1:4,i),1,0,0,1,0,an_epe(i)
	enddo !i=1,nopc_arr
	close(10)
	endif	!  other_pc_arr
	end subroutine merge_pcarr

      E N D  
                                                     
