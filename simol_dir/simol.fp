c     keys  
c     :NMV_GX  - switch of projection mode(important in calc.
c		 with point charges
C     :TSLOC   - transition state locator

c simol.f:
c   MAIN simol:
c   converg:
c   seta:
c       entry    eta:
c   pdet:
c   prcyc:
c   linmin:
c   drv:
c   force:
c   dr:
c   r:
c   dt:
c   t:
c   df:
c   f:
c   defjdc:
c   trans:
c   update:
c   cnct1:
c   cnct2:
c   cnct3:
c   move:
c   rot:
c   overlp:
c   ddcopy:
c   newpnt:
	PROGRAM SIMOL
	use smlcom
	use optsml, only:	x=>xo,g
	use flpocm
	use simol_io, only:	eyead
	
 	implicit real*8 (a-h,o-z)
c       tapes :   9 - gross cycle points 
c		 10 - linmin shearch points 
c		 16 - previous hessian saved

	real*8, dimension(nlocv):: xlast,glast,gg
	dimension 
     .    h(nlocm)  ! hessian
     .   ,gd(nlocv)  ! work arrays to save current values of x-g coordinats
     .   ,xd(nlocv)  ! 
     .   ,dell(nlocv)!  initial steps
     .   ,dg(nlocv)  ! additional vector for BFGS update

	data 
     .  cncos/0.001d0/  ! restart parametor
     . ,pnorm/1.d0/      ! dummy initialization
	logical	ex
		write(*,*) 'you are using v180l version' 

!	open(9,file='ftn09',form='unformatted', status='unknown')
	open(9,file='fort.9',form='unformatted', status='unknown')
!	open(iu152,file='ftn92',form='formatted', status='old')
	open(iu152,file='fort.92',form='formatted', status='old')
        open(33, file='conv',status='unknown') !  this file will be deleted in
!                              		       !  the case of sucsesful update

        inquire (file='gsn',exist=ex)
	if (ex) then
		open (16,file='gsn',status='unknown',form='unformatted')   
	endif  



        call drv(iwork,n,dell) ! project Cartesian grad. on the intern. coord 

      if(iwork.ne.1) then
c             **   open file of linmin data   **
	open (10,file='fort.10', status='unknown',form='unformatted')
        alpha=1.d0
	fac=1.2d0
	end if ! iwork.ne.1

	gnorm=sqrt(OVERLP(1,g,g,g,n))

	go to(4,85,100,100),iwork

    4	do i=1,n   ! step aside to find zero approximation of hessian
	xd(i)=x(i)-sign(dell(i),g(i)) 
	enddo
	inquire (file='gsn',exist=ex)
	if(ex) call system('rm gsn')
	write(9) fiii,pnorm
        write(9) (x(i),i=1,n),(g(i),i=1,n)
	call update(2,n,xd) ! update geometry and exit the program

c____  based on the d_G/d_X values find H**(-1) & the search direction
   85	inquire (file='gsn',exist=ex)
	if (.not.ex) then
	write(*,*) ' file gsn is absent '
		open (16,file='gsn',status='unknown',form='unformatted')	
	     	read(9) fip1,pnorm
		read(9) (xd(i),i=1,n),(gd(i),i=1,n)
		do i=1,n   ! define actual dell values
			dell(i)=xd(i)-x(i)
		enddo ! i=1,n
		call seta(n,h,gd,fip1,dell,x,g)
	else
	        read(16) (h(i),i=1,n*(n+1)/2)
		open(26,file='gsn.txt',form='formatted',status='unknown')
		write(26,'(i5)') n*(n+1)/2
		do i=1,n*(n+1)/2
	 	 write(26,'(f15.6)') h(i)	
		enddo
	print *, 'file gsn.txt is generated'

		fip1=100000.0d0
	endif 

	if(fiii.ge.fip1) then ! the point with lesser energy is the starting one
	write(6,149) fip1,fiii
  149	format(1x,'function value ',f15.5, ' will not be replaced by'
     *  ,' value ',f15.5,/' calc. by restart procedure')
	fiii=fip1
	call ddcopy(x,xd,n)
        call ddcopy(g,gd,n)
	gnorm=SQRT(OVERLP(1,g,g,g,n))
	else
	write(6,141) fip1,fiii
  141	format(1x,'function value ',f15.5, ' is being  replaced by'
     *  ,' value ',f15.5,/' found in  restart procedure')
	end if ! fiii.ge.fip1/else

	call pdet(n,h,pnorm,gnorm,coss)
        call prcyc(n,x,g,gnorm,coss,fiii) ! print the cycl parametors
	pnlast=pnorm

    8   continue 
	if(abs(coss).le.cncos.and.iwork.ne.2)  then
	write(6,*) 'the search is restarted due too small coss value'
	go to 4 
	endif   ! coss.le.cncos

c---------- ** linmin initialization section **
	rewind 9
	write(9) fiii,pnorm 
	if(abs(pnorm).lt.1.0d-10) then
	alpha=0.0d0
	else
	alpha=alpha*pnlast/pnorm ! now the size of a new linmin step is defin.
	endif
	yead=dabs(alpha*dot)
	write(6,7136) yead
 7136	format(' ' ,' -alpha*p*g ', f18.6) 
	if(yead.lt.eyead) then
	write(6,7135)
 7135	format('  herberts test satisfied')
        write(6,*)'--------------------------------------------'
	stop 
	end if
c	save previous gross cycle generalized coord
	write (9)(x(i),i=1,n),(g(i),i=1,n)

	iwork=2
  100   call linmin(iwork,n)
c       **   analyze  linmin results **
	rewind  9 ! read energy & direction normalization at the cycle begin.
	read (9) smval,pnlast

	if(smval.lt.fiii.and.fiii-smval.gt.1.d-5) then
	print 480
  480	format(' no point lower in energy  than   the starting'
     *        /' point could be found in the line minimization')
	stop
	fiii=smval
	read(9) (x(i),i=1,n),(g(i),i=1,n)
	go to 4 ! linmin search is unsecsesfull, go to the restart section

	else ! chose new linmin direction 
        read(9) (xlast(i),i=1,n),(glast(i),i=1,n)
	end if ! smval.lt.fiii/else
	call converg(gnorm,pnorm,smval,n)   ! check convergence
        call eta(n,h,gg,xlast,glast,dg,x,g) ! update hessian matrix
        call pdet(n,h,pnorm,gnorm,coss) ! find new direcrionf of linear search
        call prcyc(n,x,g,gnorm,coss,fiii) ! print the cycl parametors
	go to 8 
	e n d

	subroutine converg(gnorm,pnorm,smval,n)
	use optsml, only: x=>xo,g
	use flpocm, only: alpha,fiii
 	 implicit real*8 (a-h,o-z)
!	common /flpocm/tol(305),alpha,fiii
!     ./opt/ x(300),g(600)/flpocm/tol(305),alpha,fiii
     	dimension tole(3)
	data tole/
     >  1.d-4   ! convergens by coordinate norm
     > ,5.d-7   ! energy convergence
     > ,5.d-5/  ! convergens by gradient value criteria
        xn=sqrt(OVERLP(1,x,x,x,n))
	tx=abs(alpha*pnorm)
	if(xn.ne.0.d0) tx=tx/xn
	tf=abs(fiii-smval) ! difference gained in the macrocycle
	if(tx.gt.tole(1)) go to 20
	print 511
  511   format( ' test on x sutisfied')
	go to 22
   20	if(tf.gt.tole(2)) go to 21
	print 512
  512	format( ' test on f satisfied')
	go to 22
   21	if(gnorm.gt.tole(3))  return
	print 513
  513	format( ' test on gradient satisfied')
   22	do i=1,n
	if (abs(g(i)).gt.tole(3)/2.d0) return
	end do
	print 515
  515	format( ' Peters test satisfied ')
	stop 
	e n d

	subroutine seta(n,h,gd
     > ,fiii     ! previous value of energy
     > ,dell,x,g)
c    this subroutine sets initial value of the hessian
	use smlcom, only:	nlocv,mu_type,pmste
	use flpocm, fip1=>fiii
 	implicit real*8 (a-h,o-z)
C!!!	common
C     &/steps/ pmste
	dimension h(n*(n+1)/2),gd(1),gg(n),xlast(n),glast(n)
     & ,dg(*),dell(*),x(n),g(n)
     & ,tdel(3)

	data zcons/5.d-5/,tdel/0.3d0,1.5d0,1.5d0/ 

	 do i=1,n*(n+1)/2
	h(i)=0.d0
	enddo  

	do i=1,n
	ii=i+n*(i-1)-((i*(i-1))/2)

	ggd=abs(gd(i)) ! ggd abs value of grad in the best point
	if(fiii.lt.fip1) ggd=abs(g(i))

	ggggg=gd(i)-g(i)
        if(abs(ggggg).lt.zcons) then
	h(ii)=abs(0.01/g(i))
	write(*,*) 'for variable ',i 
        write(*,*) 'h(ii)=0.01/g because diff of g is LT zcons'

	else 		! |del_g| > zcons
        h(ii)=dell(i)/ggggg

	if(h(ii).lt.0.d0) then
	write(*,*) 'for variable ',I, '  h-(ii) is lt 0'
		if(ggd.lt.zcons) then
		h(ii)=2.0
		else
		mut=mu_type(i)
		h(ii)=tdel(mut)*abs(dell(i))/ggd
		write(*,*) h(ii),abs(dell(i))*tdel(mut)
		endif ! ggd.lt.zcons/else
	endif ! h(ii).lt.0.d0
	endif ! dabs(ggggg).lt.zcons/else

     	if(ggd.lt.zcons) then
	pmstep = abs(pmste/zcons)
	else
	pmstep = abs(pmste/ggd)      
	endif
	if(h(ii).gt.pmstep)  then   ! decrease too big steps
	write(*,*) 'for variable ',i
        write(*,*) 'h(i,i)=pmstep= ',pmstep,' old h(i,i) ',h(ii)
        h(ii)=pmstep
	endif  ! h(ii).gt.pmstep

	enddo	
	write(16) (h(i),i=1,n*(n+1)/2)
	return

	entry eta(n,h,gg,xlast,glast,dg,x,g)
c **   current updating of revirsed hessian **
       	rewind 16
	read(16) (h(i),i=1,n*(n+1)/2)
	sy=0.d0
	yhy=0.d0
	do i=1,n
	gg(i)=0.d0
	     do k=1,n
	    ik=i+n*(k-1)-((k*(k-1))/2)
	   if(k.gt.i) ik=k+n*(i-1)-((i*(i-1))/2)
	  gg(i)=gg(i)+h(ik)*(g(k)-glast(k))
	end do ! do k=1,n
	   y=g(i)-glast(i)		! DG
	  yhy=yhy+gg(i)*y  		! FAE
	 sy=sy+(x(i)-xlast(i))*y 	! FAC
	end do ! do i=1,n


            do i=1,n ! dg makes BFGS different from DFP
	   dg(i)= (x(i)-xlast(i))/sy - gg(i)/yhy
	  enddo ! do i=1,n

	    do i=1,n
	   y=x(i)-xlast(i)
	  do k=i,n
	 ik=k+n*(i-1)-((i*(i-1))/2)
	h(ik)=h(ik)+y*(x(k)-xlast(k))/sy-gg(i)*gg(k)/yhy
     >                   +yhy*dg(i)*dg(k) ! BFGS additional term
	end do ! k=i,n
	end do ! i=1,n
 	rewind 16
 	write(16) (h(i),i=1,n*(n+1)/2)
	return

	entry hdeta(n,h,gg,xlast,glast,x,g)
c **   		current update of direct hessian **

       	rewind 16
	read(16) h
c          WRITE(6,'(10x,''hdeta: HESSIAN IN INTERNAL COORDINATES'')')
c          CALL VECPRT(H,N)

	write(*,*) xlast
        write(*,*) x
        write(*,*) glast
        write(*,*) g

	sy=0.d0
	yhy=0.d0
	do i=1,n
	gg(i)=0.d0
	ik=i*(i-1)/2
	     do k=1,i
 	    ik=ik+1
	  gg(i)=gg(i)+h(ik)*(x(k)-xlast(k))
	enddo ! do k=1,n
	     do k=i+1,n
	    ik=ik+k-1
	  gg(i)=gg(i)+h(ik)*(x(k)-xlast(k))
	enddo ! do k=1,n
	gg(i)=g(i)-glast(i)-gg(i)

	   y=x(i)-xlast(i)		! DG
	  yhy=yhy+gg(i)*y  		! FAE
	 sy=sy + y*y 			! FAC
	end do ! do i=1,n
	yhysy = yhy/sy

	    do i=1,n
	  do k=1,i
	ik=i*(i-1)/2+k
	h(ik)=h(ik)+(gg(i)*(x(k)-xlast(k)) + gg(k)*(x(i)-xlast(i))
     &		   - yhysy*(x(k)-xlast(k))*(x(i)-xlast(i)) )/sy
	end do ! k=i,n
	end do ! i=1,n
c	rewind 16
c	write(16) h
	return
	e n d

	subroutine pdet(n,h,pnorm,gnorm,coss)
c ** define direction of the new line search **
	use optsml, only: x=>xo,g
	use flpocm
	implicit real*8 (a-h,o-z)
!	common/flpocm/tolit(3),fac,p(300),dot,alpha,fiii
!     .  /opt/x(300),g(600)
	dimension h(*)
	pnorm=0.d0
	dot=0.d0
	do k=1,n
	s=0.d0
	do i=1,n
	ik=i+n*(k-1)-((k*(k-1))/2)
	if(k.gt.i) ik=k+n*(i-1)-((i*(i-1))/2)
	s=s-h(ik)*g(i)
c	if(k.eq.10) write(*,*) i,g(i),h(ik)
	end do ! i=1,n
	p(k)=s
	pnorm=pnorm+p(k)**2
	dot=dot+p(k)*g(k)
	enddo ! k=1,n
	pnorm=sqrt(pnorm)
	if(abs(pnorm*gnorm).lt.1.0d-17) then
		coss=1.0d0
	else
		coss=-dot/(pnorm*gnorm)
	endif
	return
	e n d

	subroutine prcyc(n,x,g,gnorm,coss,fip1)
	use flpocm, only: p
	implicit real*8 (a-h,o-z)
!	common/flpocm/tol(4),p(300),dummy(3)

	dimension g(n),x(n)
	write(6,500) fip1
 500	format(' at the begining of cycle the function value'
     *  ,' is',f20.6/' the current point is ...')
	nto=n/10
	nrem=n-nto*10
	iinc1=-9
	if(nto.ge.1) then      
	do  i=1,nto
	iinc1=iinc1+10
	iinc2=iinc1+9
	write(6,852) (j,j=iinc1,iinc2)
	write(6,853) (x(j),j=iinc1,iinc2)
	write(6,854) (g(j),j=iinc1,iinc2)
	write(6,855) (p(j),j=iinc1,iinc2)
	end do
	end if
	if (nrem.ge.1) then
	iinc1=iinc1+10
	iinc2=iinc1+(nrem-1)
 	write(6,852) (j,j=iinc1,iinc2)
	write(6,853) (x(j),j=iinc1,iinc2) ! values of internal coordinates
	write(6,854) (g(j),j=iinc1,iinc2) ! gradients
	write(6,855) (p(j),j=iinc1,iinc2) ! steps
	end if
	write(6,507) gnorm,coss
 507	format(' gradient norm ',f11.5, ' angle cosine ',f11.5)	
 852	format(4x,'i',6x,i3,9(9x,i3))  		
 853	format(2x,'x(i)',1x,f10.5,2x,9(f10.5,2x))  		
 854	format(2x,'g(i)',1x,f10.5,2x,9(f10.5,2x))  		
 855	format(2x,'p(i)',1x,f10.5,2x,9(f10.5,2x))  		
	return
	e n d

	subroutine linmin(iwork
     .  ,nn     	! number of variables
     .  )
	use optsml,only:	x=>xo,q=>g
	use flpocm, only: fac, p, g=>dot, t=>alpha, f=>fiii
	implicit real*8 (a-h,o-z)
	real *8 
     .  ita,lambda,tolit(3)
     . ,tlast ! best of the linmin previous results
     . ,flast ! the same for energy
	integer 
     .  ihk   ! =1 at initialization (iwork =1)
	      ! =2 if left bracket was set in previous inerpolation
	      ! =3 if right bracket was set in previous inerpolation

!	common/opt/  x(300)
!     &             ,q(600)
!	common
!     &     /flpocm/ tol(3),fac,p(300)
!     &  ,g
!     &  ,t ! at the start of the search t = alpha 
!     &  ,f

	data lambda/0.d0/, rlal/0.d0/,tolit/1.d-5,1.d-7,1.d-5/

	if(iwork.ge.3)	then
	rewind 10
	read(10) (p(i),i=1,nn)! direction vector
	endif ! iwork.ge.3

	xmaxm=0.d0
	  do i=1,nn 
	 if(dabs(p(i)).gt.xmaxm) xmaxm=dabs(p(i))
	end do
       xmaxm=0.08d0/xmaxm

        g=OVERLP(1,p,q,q,nn)
	
	go to (99,99,107,111), iwork

   99 	t=t*fac
	ihk=1
	tlast=0.d0
	flast=0.d0
	fb=f ! function value at the beginning of the interval
	gb=g ! the same for gradient value

	if(gb.gt.0.d0) then   ! change direction of search 
	do i=1,nn
	p(i)=-p(i)
	end do
	gb=-gb
	end if

	rewind 10
	write (10) (p(i),i=1,nn) ! direction vector
	ita=t
	t=0.d0

    7 	fa=fb ! here we are looking for the right bracket tb
	ga=gb
	ta=t
   71   if(iwork.gt.2.and.ita.gt.xmaxm) ita=xmaxm
	t=t+ita
	tb=t
	write(10) ihk,ta,fa,ga,tb,fb,gb,tlast,flast,t,ita,lambda,rlal
	backspace 10
        call NEWPNT (ita,p,x,x,nn)
	call update(3,nn,x) !iwork=3 with new data go to 107

  107  read(10) ihk,ta,fa,ga,tb,fb,gb,tlast,flast,t,ita,lambda,rlal
	backspace 10
	gb=g   
	fb=f
        print 86
   86   format(' ' ,2x,'i',5x,'xa',9x,'x',8x,'xb',13x,'fa',13x,'f',14x
     .  ,'fb',12x,'ga',8x,'g',9x,'gb')

	write(6,88)iwork,ta,tb,fa,fb,ga,gb
	if(gb.lt.0.d0.and.fb.gt.fa) return ! unsucsesful search due to
!					   ! inacuracy of gradient.
!					   ! try to change direction.
	if(.not.((gb.lt.0.d0).and.(fb.lt.fa))) go to 11
	ita=4.d0*ita
	go to 7 ! make on more but increased step

   11	z=3.d0*(fa-fb)/ita+ga+gb
	w=z*z-ga*gb
	if(w.lt.0.d0) w=0.d0
	w=sqrt(w)
	temp=ga+7.d0
	if(temp.ge.0.d0) lambda=ita*(1.d0-(ga+z+w)/(ga+gb+2.d0*z))	
	if(temp.lt.0.d0) lambda=ita*(1.d0-ga/(ga+z-w))	
	if((lambda.le.0.d0).or.(lambda.ge.ita)) lambda=ita/2.d0
	t=tb-lambda
	hkcon=abs((t-tlast)*g)
		if(ihk.ne.1.and.hkcon.lt.1.d-3) then
		t=tlast
		if(ihk.eq.2) call NEWPNT (-rlal,p,x,x,nn)
		return
		end if
	write(10) ihk,ta,fa,ga,tb,fb,gb,tlast,flast,t,ita,lambda,rlal
        call NEWPNT (-lambda,p,x,x,nn)
	call update(4,nn,x) ! iwork=4 with new data go to 111

  111	read(10)ihk,ta,fa,ga,tb,fb,gb,tlast,flast,t,ita,lambda,rlal 
	backspace 10

c---------------------------------------------------------------------
c	print results of cubic interpolation
c
	write(6,87)iwork,ta,t,tb,fa,f,fb,ga,g,gb
   87	format(2x,i2,3x,e10.3,2(2x,e10.3),3x,f15.6,2(2x,f15.6)
     * 	,3x,e10.3,2(2x,e10.3))
   88   format(1x,i2,3x,e10.3,2x,e10.3,10x,3x,f15.6,2x,f15.6,13x
     *  ,3x,e10.3,2x,e10.3)

c------------------------------------------------------------------
c	now interval for interpolation is just to big to have
c	good result from the first try
c
	if(fa.lt.f.or.fb.lt.f)  then
	if(fa.lt.f.and.g.gt.abs(ga)) then ! reset right bracket
	ita=ita-lambda
	fb=f
	gb=f
	tb=t
	goto 11
	else
	return     ! ##
	endif
	endif	! fa.lt.f.and.g.gt.abs(ga)
c-------------------------------------------------------------------

	ttest=abs(tlast-t)
	if(t.ne.0.d0) ttest=ttest/dabs(t)
		if(ttest.lt.tolit(1)) then 
	write(*,*) 'linmin: t-diff. criterion is satisfied'
	return
	endif ! ttest.lt.tolit(1)

	ftest=abs(flast-f)
	if(ftest.lt.tolit(2)) then
	write(*,*) 'linmin: f-diff. criterion is satisfied'
	return
	endif ! ftest.lt.tolit(2)

	if(abs(g).lt.tolit(3)) then
	write(*,*) 'linmin: g-diff. criterion is satisfied'
	return
	endif ! abs(g).lt.tolit(3)

	tlast=t
	flast=f
	rlal=lambda

	if(g.gt.0.d0) then ! set t as right bracket
	ita=ita-lambda
	ihk=3
	fb=f
	gb=g
	tb=t
	go to 11 ! make one more cubic interpolation
	end if

	if(f.lt.fa ) then ! set t as left bracket
	ita=lambda
	ihk=2
	fa=f
	ga=g
	ta=t
	t=tb
        call NEWPNT (lambda,p,x,x,nn)
	go to 11 ! make one more cubic interpolation
	end if

	ihk=1 ! fa is again left bracket but step is decreased
	z=ita-lambda
	ita=z/4.d0
c	fb=fa
c	gb=ga
c	t=ta
        call NEWPNT (-z,p,x,x,nn)
	go to 71
	e n d

	subroutine drv(iwork,kolvar,dell)
	use smlcom, only: sxo
	use xyzsml, only: x,y,z,kolat
	use optsml, only: xo,g
	use simol_io, only:	hessp,l_calc_h,jopx,icycle,mode,freqp
     &			       ,simol_inp
	implicit real*8 (a-h,o-z)
	logical ex    ! c:tsloc - true to calc. Hessian
	real*8, dimension(*):: dell
                                            
CC:IBM>>
C      character*26     str
C      external         fdate_
C      call fdate_(str)
C      write(6,*) '            ',str
CC:IBM<<

CC:HPUX>>
C      character*24     fdate
C      external         fdate
C      write(6,*) '              ',fdate()
CC:HPUX<<


CFPP:DEC!      character*24     fdate
CFPP:DEC!      external         fdate
CFPP:DEC!      write(6,*) '              ',fdate()


CC:CRAY>>
C      character*24     fdate
C      external         fdate
C      write(6,*) '            ',fdate()
CC:CRAY<<

CFPP:SGI!      character*24     fdate
CFPP:SGI!      external         fdate
CFPP:SGI!      write(6,*) '            ',fdate()

      call simol_inp(kolat,kolvar,iwork,dell,n_00)
	inquire (file='gsn',exist=ex)
        if(.not.ex.and.iwork.gt.2)
     &  open (16,file='gsn',status='unknown',form='unformatted')

	if(iwork.gt.kolvar+10.and.l_calc_h) then  !c:tsloc
	stop
C:FREQ>>
	elseif(iwork.gt.kolvar*2+10) then
        open(35,file='gx_grad')
        open(34,file='xo',form='unformatted')
        call freqp(kolat,kolvar,iwork,n_00)
        endif 		! iwork.gt.kolvar*2+10
C:FREQ<<

      call trans(kolat,iwork) ! project gradients on internal coordinates.

      if(iwork.eq.5) call  hessp(kolat,kolvar,iwork) ! c:tsloc

C:FREQ>>
        if(iwork.eq.10) then
        write(34) (xo(k),k=1,kolvar) 
        write(34) (g(k),k=1,kolvar) 
        write(34) (sxo(k),k=1,kolvar)
        rewind 34

	if(l_calc_h) then	       ! c:tsloc
        call hessp(kolat,kolvar,iwork) ! start freq. calc.; calc stops hear

	else
        call freqp(kolat,kolvar,iwork,n_00) ! start freq. calc.; calc stops hear
	endif ! l_calc_h
	endif ! iwork.eq.10

     	if(iwork.ge.11)  then
        open(34,file='xo',form='unformatted') 


CC:HPUX>>
C	if(iwork.eq.11) then
C	open(35,file='gx_grad.sav')
C	else
C        open(35,file='gx_grad.sav',access='append')
C	endif ! iwork.eq.11/else
CC:HPUX<<

C:DEC>>
	if(iwork.eq.11) then
	open(35,file='gx_grad.sav',status='unknown')
	else
!!!        open(35,file='gx_grad.sav',access='append',status='unknown')
        open(35,file='gx_grad.sav',position='append',status='unknown')
	endif ! iwork.eq.11/else
C:DEC<<

CC:SGI>>
C        if(iwork.eq.11) then
C       open(35,file='gx_grad.sav')
C       else
C            open(35,file='gx_grad.sav',access='append')
C       endif ! iwork.eq.11/else
CC:SGI<<

CC:CRAY>>
C        if(iwork.eq.11) then
C       open(35,file='gx_grad.sav')
C       else
C            open(35,file='gx_grad.sav',POSITION='APPEND')
C       endif ! iwork.eq.11/else
CC:CRAY<<

        write(35,3510) (g(i),i=1,kolvar)
 3510   format((5f15.8))
        close(35)

	if(l_calc_h) then	!c:tsloc
	call hessp(kolat,kolvar,iwork)
	else
	call freqp(kolat,kolvar,iwork,n_00)
	endif ! l_calc_h
	endif ! iwork.gt.10

C:FREQ<<

	if(iwork.eq.9) then
	write(6,*) ' iwork = -9: calc energy & gradients for a point'
	stop
	endif ! iwork.eq.9
      return
      e n d

	subroutine calc_type(incod_gx,icrel)
	if(icrel.ne.0) write(6,*) 'Relativistic '
	if(incod_gx.eq.90) 
     >  write(6,*) 'Perdew Wang and BLYP energies for the point are:'
	if(incod_gx.eq.88) 
     >  write(6,*) 'Becke Perdew  and BLYP energies for the point are:'
	if(incod_gx.eq.89) 
     >  write(6,*) 
     &  'RVWN Becke Perdew  and BLYP energies for the point are:'
        if(incod_gx.eq.81)
     >  write(6,*) 'VWN and BLYP energies for the point are:'
        if(incod_gx.eq.85)
     >  write(6,*) ' BLYP energy for the point is:'
	return
	e n d

	subroutine eq_update(xo,kolvar)
	real*8 xo(kolvar)
	logical xoexist

c ** with this option the geomentry can be ubdated to the equlibtium one **
        inquire (file='./xo',exist=xoexist)
        open(34,file='xo',form='unformatted',status='unknown')
        if(xoexist)  then
        read(34) (xo(k),k=1,kolvar)
        rewind 34
	call update(9,kolvar,xo) 		! -> stop
        else
        write(6,*) 
     &  'file with optimal values of internal coordinates (xo)'
        write(6,*) 'is not found'
        stop
        endif ! xoexist/else
	return
	e n d

C:XMOL>>
	subroutine xmol_o(kolat,an,x,y,z,sc_coeff)
	real*8 an(*),x(*),y(*),z(*),sc_coeff
	character*2  symv(100)
      DATA symv/'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE','NA',
     1 'MG','AL','SI',' P','S ','CL','AR',' K','CA','SC','TI','V ','CR',
     2 'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',
     3 'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN',
     4 'SN','SB','TE','I ','XE','CS','BA','LA','CE','PR','ND','PM','SM',
     5 'EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',' W','RE',
     6 'OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA',
     7 'AC','TH','PA',' U','NP','PU','AM','CM','BK','CF','X ','FM'/

	open(7,file='xmol.xyz',status='unknown')
	write(7,'(i3)') kolat
	write(7,*) 'comment line'
	do i=1,kolat
	write(7,777) symv(int(an(i)+0.0001))
     >  , x(i)*sc_coeff,y(i)*sc_coeff,z(i)*sc_coeff
 777    format(a2,3f15.6)
	enddo
	stop 'xmol'
	e n d
C:XMOL<<
	subroutine freq_open
	write(6,*) '------------------------------------------------'
	write(*,*) '               Frequency calculations '
        open(34,file='xo',form='unformatted',status='unknown')
	open(68,file='dipmom.sav',status='unknown')
	close(68,status='delete')
	return
	entry hess_open
        write(6,*) '------------------------------------------------'
	 write(*,*) '               Calculation of the Hessian      '
        open(34,file='xo',form='unformatted',status='unknown')
	return
	e n d

	subroutine scale_steps(kolvar,dell,xo_eqvl,xo,kolat
     &			       ,xx,yy,zz)
	use smlcom, only: nlocv, kamax
	use xyzsml, only: x,y,z,gra
	implicit real*8 (A-H,O-Z)
	dimension dell(kolvar),xo_eqvl(kolvar),xo(kolvar)
     &		 	,xx(kolat),yy(kolat),zz(kolat)

	xo_eqvl(1:kolvar)=xo(1:kolvar)
	do i=1,kolvar

	do k=1,kolat
	xx(k)=x(k)
	yy(k)=y(k)
	zz(k)=z(k)
	enddo ! k=1,kolat

	xo(i)=xo_eqvl(i)+dell(i)
	call update(0,kolvar,xo)

	dism=0.d0
	do k=1,kolat
	do l=k+1,kolat
	dis=sqrt(
     &           (xx(k)-xx(l))**2+(yy(k)-yy(l))**2+(zz(k)-zz(l))**2)
     &		                     -r(k,l)
	if(abs(dis)/r(k,l).gt.dism) then
	 dism=abs(dis)/r(k,l)
	 dd=abs(dis)
	 endif
	enddo
	enddo
	if(dd.gt.0.01d0) then
	dell(i)=dell(i)*(0.01/dd)
	endif ! dd.gt.0.01d0
	call update(0,kolvar,xo_eqvl)
	xo(i)=xo_eqvl(i)
	
        enddo ! i=1,kolvar
	return
	e n d

	subroutine def_i_tpr(iwork,i_tpr,kolvar)
        write(6,*)
	if(iwork.eq.0.or.iwork.eq.9) then
	write(6,*) 
     &   'the number of independent variables is equel to ',kolvar
	i_tpr= 1
	else
	i_tpr=0
	endif ! iwork.eq.0.or.iwork.eq.9
	return
	e n d

	subroutine pcgen(ii,kolat,numx,zmat,x,y,z,an,ieq,kamax)
	implicit real*8 (A-H,O-Z)
	dimension x(kamax),y(kamax),z(kamax),an(kamax)
	integer numx(kamax,3),zmat(kamax,3),ieq(kamax)
        open(7,file='pc.xyz',STATUS= 'UNKNOWN')
	do i=1,kolat
	id=sign(1,numx(i,3)-zmat(i,3))
	do iz=zmat(i,3),numx(i,3),id
	do ix=zmat(i,1),numx(i,1)
	do iy=zmat(i,2),numx(i,2)
	x_pc=x(i)+ix*x(ii)
	y_pc=y(i)+iy*y(ii)
	z_pc=z(i)+iz*z(ii)
	q_pc=(1-2*mod(abs(ix)+abs(iy)+abs(iz),2))*(an(i)*ieq(i))
	if( ix.eq.0.and.iy.eq.0) then
	k1=1
	elseif( ix.eq.0.and.iy.ne.0 .or. ix.ne.0.and.iy.eq.0) then
	k1=2
	else
	k1=4
	endif
	k2=0
	k3=0
	k4=1
	k5=0
	write(7,15000) x_pc,y_pc,z_pc,q_pc,k1,k2,k3,k4,k5
15000 FORMAT(4d15.8,I3,4I2) !cpps
	enddo
	enddo
	enddo
	enddo
	stop 'pc are generated'
	e n d


	subroutine hess(npoi,d2e,n,g,xo,EIGVAL,h,u
     &                 ,da,F,ENERGY,sxo
     &		       ,pg,g_line,h_line
     &		       ,h_save,SMAX,gud		    ) 	! -> stop here
	use simol_io, only:	hessp,l_calc_h,jopx,icycle,mode
	implicit real*8 (a-h,o-z)
	logical DEBUG,OK,over,CONVRG,over_step

	parameter (ZERO=0.d0
     & ,EIGMIN=0.001d0,EIGMAX=25.d0,IMAX=15
     & ,GNCRIT=0.0020D0, SNCRIT=0.0081 
     & ,DEBUG = .true.) 

	dimension d2e(n,n),gud(n,npoi),xo(n),da(n)
     &           ,EIGVAL(n),h((n*n+n)/2),u(n,n),G(N),F(N),sxo(n)
     &           ,xo_last(n),g_last(n),h_save((n*n+n)/2)

!	common
!     & /tsloc/ l_calc_h,jopx,ICYCLE,mode

        open(36,file='gx_grad.sav',STATUS= 'UNKNOWN')
c-------------------------------------------------------------------




	do i=1,n
	read(36,3510) (gud(k,1),k=1,n)  ! grad at x+delta
 3510 format((5f15.8))
        d2e(i,i)=(gud(i,1)-g(i))/sxo(i)
        do j=1,n
	  if(i.ne.j) then
          d2e(j,i)=(gud(j,1)-g(j))/sxo(i)
	  endif ! i.eq.j
	enddo ! j=1,n
        enddo ! i=1,n 
        close(36)

	write(6,*) ' H matrix'
        do i=1,n
        write(6,*) (d2e(j,i),j=1,n)
        enddo 
	write(6,*)
	h(1)=d2e(1,1)
        do i=2,n   		! make force matrix symmetric
	h((i*i+i)/2)=d2e(i,i)
        do j=1,i-1
	IJ=J+(I*I-I)/2
        h(ij)= (d2e(j,i)+d2e(i,j))/2.d0
        enddo ! j=1,n
        enddo ! i=1,n 
	write(16) h

      over = jopx.ge.2
      IF(JOPX.EQ.3)THEN
        DO I=1,N
          DA(I)=ZERO
        ENDDO
        IX=0
 1      READ(5,*) L,VAL
        IF(L.GT.0)THEN
        WRITE(6,'(5x,''L,VAL'',I5,F6.2)') L,VAL
         IX=IX+1
         DA(L)=VAL
         GOTO 1
       ENDIF  ! L.GT.0

       IF(IX.EQ.0)THEN
         OVER=.FALSE.
        ELSE
      WRITE(6,'(//10x,''SEARCH DIRECTION FROM INPUT'')')
      WRITE(6,'(/10x,''DA(I) : '',12f8.2)')(DA(I),I=1,N)
       ENDIF   ! IX.EQ.0

       CALL VECNRM(N,DA,SNORM)
       ENDIF    ! JOPX.EQ.5

	entry locmod(npoi,n,g,xo,EIGVAL,h,u
     &                 ,da,F,ENERGY,sxo
     &                 ,e_last,xo_last,g_last
     &		       ,pg,g_line,h_line
     &		       ,h_save,SMAX		    ) 	! -> stop here

      over = jopx.ge.2

      IF(DEBUG)THEN
          WRITE(6,'(10x,''HESSIAN IN INTERNAL COORDINATES'')')
          CALL VECPRT(H,N)
      ENDIF


 	do i=1,n*(n+1)/2
 	h_save(i)= h(i)
 	enddo

      CALL GIVENS(N,N,N,h,EIGVAL,U)
      CALL EVOUT( U,EIGVAL,N)

      NEGVAL=abs(mode)
      IF(NEGVAL.EQ.0)WRITE(6,'(10x,
     .   ''EIGENMODE FOLLOWING METHOD: '',A)')'MINIMUM SEARCHING '
      IF(NEGVAL.EQ.1)WRITE(6,'(10x,
     .   ''EIGENMODE FOLLOWING METHOD: '',A)')'TS  SEARCHING '

      NEG=0
      ITV=1
      DO  I=1,N
      IF(EIGVAL(I).LT.ZERO) THEN
          NEG=NEG+1
      ENDIF
      IF(DABS(EIGVAL(I)).LT.EIGMIN) THEN
           EIGVAL(I)=SIGN(EIGMIN,EIGVAL(I))
           IF(DEBUG) WRITE(6,1001) I,EIGVAL(I)
      ELSEIF(DABS(EIGVAL(I)).GT.EIGMAX) THEN
           EIGVAL(I)=SIGN(EIGMAX,EIGVAL(I))
           IF(DEBUG) WRITE(6,1002) I,EIGVAL(I)
      ENDIF
      ENDDO  ! I=1,N

C
        OK=NEGVAL.EQ.NEG

C     CHOICE OF EIGENMODE
      IF(OVER) THEN
	if(icycle.eq.1) then
        DO I=1,N
          DA(I)=U(I,ITV)
        ENDDO
	endif ! icycle.eq.1

        OVMAX=ZERO
        DO I=1,N
         OVERL=ZERO
         DO J=1,N
          OVERL=OVERL+ DA(J)*U(J,I)
         ENDDO
         OVERL=ABS(OVERL)
          IF(OVERL.GT.OVMAX)THEN
               OVMAX=OVERL
               ITV=I
          ENDIF  ! OVERL.GT.OVMAX
       ENDDO
      WRITE(6,'(10x,''OVERLAP WITH THE TV = '',F8.3)') OVMAX
      ENDIF   ! OVER
      if(icycle.eq.0) stop ' with icycle.eq.0'
	

      IF(JOPX.GE.2)THEN
        DO I=1,N
          DA(I)=U(I,ITV)
        ENDDO
      WRITE(6,'(10x,''COMPONENTS OF THE TRANSITIONS VECTOR: '')')
      WRITE(6,'(12x,''   U      VAR''/)')
      DO I=1,N
       UU=U(I,ITV)
       IF(ABS(UU).GT.0.3)THEN
          WRITE(6,'(12x,F6.2,(3x,i3),3x,F8.3)') UU,I
      ENDIF  ! ABS(UU).GT.0.3
      ENDDO
      ENDIF  ! JOPX.GE.4

      WRITE(6,'(20x,''CYCLE : '',I5)')ICYCLE

                   WRITE(6,1003) NEG
      WRITE(6,'(5x,''THE'',I3,'' th   EIGENMODE SELECTED '',
     .        ''AS THE TRANSITION VECTOR'')')ITV
      WRITE(6,'(12x,''   U     VAR''/)')
      DO I=1,N
       UU=U(I,ITV)
       IF(ABS(UU).GT.0.3)THEN
          WRITE(6,'(12x,F6.2,(3x,i3),3x,F8.3)')UU,I
      ENDIF
      ENDDO   ! I=1,N

 20   ICYCLE=ICYCLE+1

      CALL FUTG(F,U,G,N)
      CALL QNSTEP(ITV,MODE,U,EIGVAL,F,SXO,SNORM,G,N,SMAX)
       call gh_line(n,SXO,SNORM,g,h_save,g_line,h_line)
	pg=0.d0
	do i=1,n
	pg=pg+SXO(i)*g(i)
	enddo
       write(*,*) 'line g,h pg'
       write(*,*)  g_line,h_line,pg


      gnorm=sqrt(OVERLP(1,g,g,g,n))
      IF(ICYCLE.GT.IMAX) then
         WRITE(6,'(15x,''NO CONVERGENCE :AFTER '',I3,'' CYCLES''
     .   ,/5x,''GNORM = '',F9.5,3x,'' SNORM = '',F9.5,
     .   3x,''ENERGY = '',F16.6 )')ICYCLE,GNORM,SNORM,ENERGY
	stop 'NO CONVERGENCE'
	endif ! ICYCLE.GT.IMAX

 	over_step=.false.
 	if(ICYCLE.gt.2) then
 	do i=1,n
 	EIGVAL(i)=xo(i)-xo_last(i)
 	enddo
 	sno=sqrt(OVERLP(1,EIGVAL,EIGVAL,EIGVAL,n))
 	write(*,*)' pg,ENERGY-e_last,sno'
 	write(*,*) pg,ENERGY-e_last,sno
 	t=1.d0
   99	SR=t*t*abs(ENERGY-e_last-pg)/sno**3
 	SL=0.3d0*abs(g_line+t*h_line/2.d0)
 	if(sr.gt.sl.and.sno.gt.0.01d0) then
 	write(*,*) 'overstep ',SR,SL
c	over_step=.true.
c	t=t/2.d0
c	sno=sno/2.d0
c	goto 99
 	endif

c	if(over_step) then
c	write(*,*) 't ',t
c       DO I=1,N
c         XO(I)=xo_last(I)+EIGVAL(I)*t
c       ENDDO
c        call update(5,N,xo)
c	endif
 	endif ! ICYCLE.gt.2
	
c	if(debug) stop 'debug '
c***************************
       	rewind 16
	write(16) (h_save(i),i=1,n*(n+1)/2)
	rewind 9
	write(9) xo,g,da,ENERGY,pg,g_line,h_line
c***************************
       DO I=1,N
         XO(I)=XO(I)+SXO(I)
       ENDDO

      WRITE(6,'(10x,''g(I) =: '',10f9.5)')(G(I),I=1,N)
      WRITE(6,'(10x,''F(I) =: '',10f9.5)')(F(I),I=1,N)
      WRITE(6,'(10x,''S(I) =: '',10f9.5)')(SXO(I),I=1,N)

       CONVRG=SNORM.LT.SNCRIT.OR.GNORM.LT.GNCRIT
      IF(CONVRG)THEN
         WRITE(6,'(15x,''SEARCH FINISHED AFTER '',I3,'' CYCLES''
     .   ,/5x,''GNORM = '',F9.5,3x,'' SNORM = '',F9.5,
     .   3x,''ENERGY = '',F16.6 )')ICYCLE,GNORM,SNORM,ENERGY
        call update(9,N,xo)
	else
         WRITE(6,'('' AFTER '',I3,'' CYCLES    E= '',F16.6,
     .   3x,''GNORM = '',F9.5,''   SNORM = '',
     . F9.5 )')ICYCLE,ENERGY,GNORM,SNORM

         call update(5,N,xo)
      ENDIF   ! CONVRG/else

1001  FORMAT(' EIGENVALUE ',I2,' TOO SMALL. REPLACED BY ',F12.4)
C
1002  FORMAT(' EIGENVALUE ',I2,' TOO BIG. REPLACED BY ',F12.4)
1003  FORMAT('      NUMBER OF NEGATIVE EIGENVALUES  ',I5)
1011  FORMAT(' TRANSITION STATE SEARCH.')

1012  FORMAT(' MINIMUM SEARCH. ')

	return
      e n d
	subroutine gh_line(n,SXO,SNORM,g,h,g_line,h_line)
	implicit real*8 (A-H,O-Z)
	dimension SXO(n),g(n),h(*)

	g_line=0.d0
	h_line=0.d0
	do i=1,n
	aux=0.d0
	g_line=g_line+g(i)*SXO(i)
	ik=i*(i-1)/2
	do k=1,i
	    ik=ik+1
	aux=aux+SXO(k)*h(ik)
	enddo
	do k=i+1,n
	    ik=ik+K-1
	aux=aux+SXO(k)*h(ik)
	enddo
	h_line=h_line+SXO(i)*aux
	enddo
	g_line=g_line/SNORM
	h_line=h_line/SNORM**2
	return
	e n d
C:TSLOC<<

	subroutine force(d2e,kolvar,gr_eqvl,xo,sxo)
	use smlcom, only: nlocv,igener
	implicit real*8 (a-h,o-z)
	dimension d2e(kolvar,kolvar),gud(nlocv,2),gr_eqvl(kolvar)
     & ,xo(kolvar),sxo(kolvar)
c      common/gen/igener(nlocv)
c	common  /sfreq/ sxo(300)  ! steps on internal coord to calc frequencies
	parameter (deps = 1.d-11)
        open(36,file='gx_grad.sav',STATUS= 'UNKNOWN')
c-------------------------------------------------------------------
	do i=1,kolvar
	read(36,3510) (gud(k,1),k=1,kolvar)  ! grad at x+delta
        read(36,3510) (gud(k,2),k=1,kolvar)  ! grad at x-delta
 3510 format((5f15.8))
CAY Now we will calculate d2e/dx2 using interpolation to zero gradient
CAY point. It is needed for d2e/dx2 to be more presize in case when
CAY de/dx are substantially non-linear.
        dgrad=(gud(i,1)-gud(i,2))/2.d0 
        write(6,*) 'sxo;ig :' , sxo(i),igener(i),dgrad
        d2e(i,i)=dgrad/sxo(i)           !/igener(i)
	if (d2e(i,i).eq.0.d0) then
	do j=1,kolvar
		d2e(j,i)=0.d0
	enddo

	else	! d2e(i,i).ne.0
        alp=((gud(i,1)+gud(i,2))/2.d0-gr_eqvl(i))+deps ! non-linearity parameter

	if(1.-4.*gr_eqvl(i)*alp/(dgrad**2).ge.0.d0) then
        sqparm=sqrt(1.-4.*gr_eqvl(i)*alp/(dgrad**2))
        write(6,*)' Non-linearity factor(',i,',',i,'):',sqparm
        xtrue=xo(i) - sxo(i)*dgrad*(1.d0-sqparm)/(2.d0*alp) 
        write(6,*)' Zero-gradient coordinate''s value = ',xtrue
        d2e(i,i)=d2e(i,i)*sqparm 

        do j=1,kolvar

	  if(i.ne.j) then
          dgrad=(gud(j,1)-gud(j,2))/2.d0
          d2e(j,i)=dgrad/sxo(i)                   !/igener(j)
          alp=(gud(j,1)+gud(j,2))/2.d0-gr_eqvl(j)   ! non-linearity parameter
          addpar=2.d0*alp/(sxo(i)**2)*(xtrue-xo(i))
          d2e(j,i)=d2e(j,i)+addpar                !/igener(j)
	  endif ! i.eq.j
	enddo ! j=1,kolvar
	else 			! 1.-4.*gr_eqvl(i)*alp/(dgrad**2).lt.0.d0
	d2e(i,i) =0.d0
        do j=1,kolvar
	d2e(j,i) =0.d0
        enddo ! j=1,kolvar
	endif 			! 1.-4.*gr_eqvl(i)*alp/(dgrad**2).ge.0.d0
	endif	! i.ne.j
        enddo ! i=1,kolvar 

	write(6,*) ' F matrix'
        do i=1,kolvar
        write(6,*) (d2e(j,i),j=1,kolvar)
        enddo 
	write(6,*)
        do i=2,kolvar   ! make force matrix symmetric
        do j=1,i-1
        d2e(j,i)= (d2e(j,i)+d2e(i,j))/2.d0
	d2e(i,j)= d2e(j,i)
        enddo ! j=1,kolvar
        enddo ! i=1,kolvar 
        close(36)
	return
      e n d




      function dr(i,j,m,ndec)
	use xyzsml,only: x,y,z,kolat
      IMPLICIT real*8 (A-H,O-Z)
!      common/xyz/x(100),y(100),z(100),gr(300),kolat
      integer dim,djm
      dim=0
      djm=0
      if(i.eq.m) dim=1
      if(j.eq.m) djm=1
      go to (1,2,3),ndec
    1 dr=( djm-dim )*( x(j)-x(i) )/r(i,j)
      return
    2 dr=( djm-dim )*( y(j)-y(i) )/r(i,j)
      return
    3 dr=( djm-dim )*( z(j)-z(i) )/r(i,j)
      return
      e n d

      function r(i,j)
	use xyzsml,only: x,y,z,kolat
      IMPLICIT real*8 (A-H,O-Z)
c     bond lengts calculation
!      common/xyz/x(100),y(100),z(100),gr(300),kolat
      r=sqrt( (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2 )
      return
      e n d

      function dt(i,j,k,m,ndec)
	use xyzsml,only: x,y,z,kolat
      IMPLICIT real*8 (A-H,O-Z)
!      common/xyz/x,y,z,gr(300),kolat
      real*8
     .  r,dr,t,tet,dt,PI
      integer dim,djm,dkm
	data PI/3.14159265358979324D+00/
      tet=t(i,j,k)
      dim=0
      djm=0
      dkm=0
      if(i.eq.m) dim=1
      if(j.eq.m) djm=1
      if(k.eq.m) dkm=1
      go to (1,2,3),ndec
    1 if(abs(tet).lt.1.d-5.or.dabs(tet-pi).lt.1.d-5) then
        dt=0.d0 
	else 
	dt=(  cos(tet)*( r(k,j)*dr(i,j,m,1)  +  r(i,j)*dr(k,j,m,1) )
     .-( (dim-djm)*(x(k)-x(j)) +  (dkm-djm)*(x(i)-x(j)) ) ) /
     . (sin(tet)*r(j,i)*r(j,k) )
	endif ! dabs(tet).lt.1.d-5/else
      return
    2 if(abs(tet).lt.1.d-5.or.abs(tet-pi).lt.1.d-5) then
	dt=0.d0
        else 
	dt=(  cos(tet)*( r(k,j)*dr(i,j,m,2)  +  r(i,j)*dr(k,j,m,2) )
     .-( (dim-djm)*(y(k)-y(j)) +  (dkm-djm)*(y(i)-y(j))) ) /
     . ( dsin(tet)*r(j,i)*r(j,k) )
        endif ! dabs(tet).lt.1.d-5/else
      return
    3 if(abs(tet).lt.1.d-5.or.abs(tet-pi).lt.1.d-5) then
        dt=0.d0 
        else  
	dt=(  cos(tet)*( r(k,j)*dr(i,j,m,3)  +  r(i,j)*dr(k,j,m,3) )
     .-( (dim-djm)*(z(k)-z(j)) +  (dkm-djm)*(z(i)-z(j))) ) /
     . ( dsin(tet)*r(j,i)*r(j,k) )
        endif ! dabs(tet).lt.1.d-5/else
      return
      e n d

      function t(i,j,k)
	use xyzsml, only: x,y,z,kolat
      IMPLICIT real*8 (A-H,O-Z)
c     angles 
!      common/xyz/x,y,z,gr(300),kolat
      rji=sqrt( (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2 )
      rjk=sqrt( (x(k)-x(j))**2 + (y(k)-y(j))**2 + (z(k)-z(j))**2 )
      rjirjk=(x(i)-x(j))*(x(k)-x(j))
     .      +(y(i)-y(j))*(y(k)-y(j))
     .      +(z(i)-z(j))*(z(k)-z(j))
      t=rjirjk/(rji*rjk)
      t=acos(t)
      return
      e n d

      function df(i,j,k,l,m,ndec)
	use xyzsml, only: x,y,z,kolat
      IMPLICIT real*8 (A-H,O-Z)
c	** derivative of dihedral angle **
!      common/xyz/x(100),y(100),z(100),gr(300),kolat
C:TORS>>
	df=0.d0
      xij=-(x(i)-x(j))/r(i,j)
      yij=-(y(i)-y(j))/r(i,j)
      zij=-(z(i)-z(j))/r(i,j)
      xjk=(x(k)-x(j))/r(k,j)
      yjk=(y(k)-y(j))/r(k,j)
      zjk=(z(k)-z(j))/r(k,j)
      xkl=(x(l)-x(k))/r(k,l)
      ykl=(y(l)-y(k))/r(k,l)
      zkl=(z(l)-z(k))/r(k,l)

	xr1=yij*zjk-zij*yjk
	yr1=zij*xjk-xij*zjk
	zr1=xij*yjk-yij*xjk

	xr2=yjk*zkl-zjk*ykl
	yr2=zjk*xkl-xjk*zkl
	zr2=xjk*ykl-yjk*xkl
	if(m.eq.i) then
      go to (1,2,3),ndec
    1 df=-xr1/(r(i,j)*sin(t(i,j,k))**2)
      return
    2 df=-yr1/(r(i,j)*sin(t(i,j,k))**2)
      return
    3 df=-zr1/(r(i,j)*sin(t(i,j,k))**2)
      return
	endif	! m.eq.i
        if(m.eq.j) then
      go to (4,5,6),ndec
    4 s1=xr1*(r(j,k)-r(i,j)*cos(t(i,j,k)))
     &                 /(r(j,k)*r(i,j)*sin(t(i,j,k))**2)
      s2=(cos(t(j,k,l))*xr2)/(r(j,k)*sin(t(j,k,l))**2)
	df=s1-s2
      return
    5 s1=yr1*(r(j,k)-r(i,j)*cos(t(i,j,k)))
     &                /(r(j,k)*r(i,j)*sin(t(i,j,k))**2)
      s2=(cos(t(j,k,l))*yr2)/(r(j,k)*sin(t(j,k,l))**2)
        df=s1-s2
      return
    6 s1=zr1*(r(j,k)-r(i,j)*cos(t(i,j,k)))
     &                /(r(j,k)*r(i,j)*sin(t(i,j,k))**2)
      s2=(cos(t(j,k,l))*zr2)/(r(j,k)*sin(t(j,k,l))**2)
        df=s1-s2
      return
        endif   ! m.eq.j
        if(m.eq.l) then
      go to (7,8,9),ndec
    7 df=xr2/(r(k,l)*sin(t(j,k,l))**2)
      return
    8 df=yr2/(r(k,l)*sin(t(j,k,l))**2)
      return
    9 df=zr2/(r(k,l)*sin(t(j,k,l))**2)
      return
	endif    ! m.eq.l
        if(m.eq.k) then
      go to (10,11,12),ndec
   10 si = -xr1/(r(i,j)*sin(t(i,j,k))**2)
      s1=xr1*(r(j,k)-r(i,j)*cos(t(i,j,k)))
     &               /(r(j,k)*r(i,j)*sin(t(i,j,k))**2)
      s2=(cos(t(j,k,l))*xr2)/(r(j,k)*sin(t(j,k,l))**2)
      sj=s1-s2
      sl=xr2/(r(k,l)*sin(t(j,k,l))**2)
	df=-(si+sj+sl)
      return
   11 si = -yr1/(r(i,j)*sin(t(i,j,k))**2 )
      s1=yr1*(r(j,k)-r(i,j)*cos(t(i,j,k)))
     &               /(r(j,k)*r(i,j)*sin(t(i,j,k))**2)
      s2=(cos(t(j,k,l))*yr2)/(r(j,k)*sin(t(j,k,l))**2)
      sj=s1-s2
      sl=yr2/(r(k,l)*sin(t(j,k,l))**2)
        df=-(si+sj+sl)
      return
   12 si = -zr1/(r(i,j)*sin(t(i,j,k))**2 )
      s1=zr1*(r(j,k)-r(i,j)*cos(t(i,j,k)))
     &               /(r(j,k)*r(i,j)*sin(t(i,j,k))**2)
      s2=(cos(t(j,k,l))*zr2)/(r(j,k)*sin(t(j,k,l))**2)
      sj=s1-s2
      sl=zr2/(r(k,l)*sin(t(j,k,l))**2)
        df=-(si+sj+sl)
      return
	endif ! (m.eq.k)
C:TORS<<
      e n d
      function f(i,j,k,l )
	use xyzsml, only: x,y,z,kolat
      IMPLICIT real*8 (A-H,O-Z)
c 	** dihedral angle calc. **
!      common/xyz/x(100),y(100),z(100),gr(300),kolat
      real*8
     . a(3),b(3),c(3)

      a(1)=(y(i)-y(j))*(z(k)-z(j)) - (z(i)-z(j))*(y(k)-y(j))
      a(2)=(x(k)-x(j))*(z(i)-z(j)) - (x(i)-x(j))*(z(k)-z(j))
      a(3)=(x(i)-x(j))*(y(k)-y(j)) - (x(k)-x(j))*(y(i)-y(j))

      b(1)=(y(l)-y(j))*(z(k)-z(j)) - (z(l)-z(j))*(y(k)-y(j))
      b(2)=(x(k)-x(j))*(z(l)-z(j)) - (x(l)-x(j))*(z(k)-z(j))
      b(3)=(x(l)-x(j))*(y(k)-y(j)) - (x(k)-x(j))*(y(l)-y(j))
      am=sqrt( a(1)**2 + a(2)**2 + a(3)**2 )
      bm=sqrt( b(1)**2 + b(2)**2 + b(3)**2 )
	if(am.lt.1.d-5.or.bm.lt.1.d-5) then
        write(6,*) 'wrong definition of dihedral angle'
	write(6,*) i,j,k,l
	stop
	endif ! i,j,k,l
      ab=a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      f=ab/(am*bm)
      if ( abs(f) .gt. 1.d0 ) f=sign(1.d0,f)
      f=acos(f)
	c(1)=a(2)*b(3)-a(3)*b(2)
	c(2)=a(3)*b(1)-a(1)*b(3)
	c(3)=a(1)*b(2)-a(2)*b(1)
	cjk=c(1)*(x(k)-x(j))+c(2)*(y(k)-y(j))+c(3)*(z(k)-z(j))
 	if(cjk.lt.0.d0) f=-f
      return
      e n d

      subroutine defjdc(kolat,kolvar,jdc,mu)
	use smlcom, only:nlocv,igener
	use xyzsml, only: x,y,z
	use zzzsml
	use numsml
	use symsml
	implicit real*8 (a-h,o-z)
      real*8 jdc(3*kolat,kolvar)
	logical lmuext
      n3=3*kolat
      mumax=0
	igener(:n3)=0
	jdc(:n3,:n3)=0.d0
      do 100 n=1,3
      do i=1,kolat
      j=abs(zmat(i,1))
      k=abs(zmat(i,2))
      l=abs(zmat(i,3))
      if(numx(i,n).ne.0) then
      mu=numx(i,n)
      igener(mu)=igener(mu)+1
      mumax=max(mumax,mu)
      do 102 iat=1,kolat
      do  m=1,3
      irow=3*(iat-1)+m
      if(n.eq.1) jdc(irow,mu)=dr(i,j,iat,m)+jdc(irow,mu)
      if(n.eq.2) jdc(irow,mu)=dt(i,j,k,iat,m)+jdc(irow,mu)
      if(n.eq.3) jdc(irow,mu)=jdc(irow,mu)
     &                         +df(i,j,k,l,iat,m)*iszm(i) !cgx_
	enddo	! m=1,3
  102 continue
	endif ! numx(i,n).ne.0
	enddo ! i=1,kolat
  100 continue
      do i=1,mumax
        do j=1,n3
c         jdc(j,i)=jdc(j,i)/sqrt(DBLE(	igener(i)))
         jdc(j,i)=jdc(j,i)/igener(i)
        enddo
       enddo

      mu=mumax   ! the last variable directly considered
	write(6,*) 'mumax' , mumax
C:MUEXT>>
        lmuext= .false.
        if(lmuext) then

      do  n=1,3
      do 104 i=1,kolat
      j=abs(zmat(i,1))
      k=abs(zmat(i,2))
      l=abs(zmat(i,3))
      if(numx(i,n).eq.0) then
	if(n.eq.1.and.j.eq.0) goto 104
        if(n.eq.2.and.k.eq.0) goto 104
        if(n.eq.3.and.l.eq.0) goto 104
      mu=mu+1
      do 103 iat=1,kolat
      do  m=1,3
      irow=3*(iat-1)+m
      if(n.eq.1) jdc(irow,mu)=dr(i,j,iat,m)+jdc(irow,mu)
      if(n.eq.2) jdc(irow,mu)=dt(i,j,k,iat,m)+jdc(irow,mu)
      if(n.eq.3) jdc(irow,mu)=jdc(irow,mu)
     &                        +df(i,j,k,l,iat,m)*iszm(i) !cgx_
	enddo	! m=1,3
  103 continue
        endif ! numx(i,n).ne.0
  104 continue
        enddo ! n=1,3
	endif ! lmuext
C:MUEXT<<

	return
	e n d

      subroutine trans(kolat,iwork)
	use math_module
	use slspar_module
	use smlcom
	use xyzsml, only:x,y,z,gr=>gra
	use zzzsml
	use numsml
	use zarsml
	use symsml, only : isg
	use optsml, only: g
	implicit real*8 (a-h,o-z)
c**************************************************************c
c       transition matrix and inversion to internal coord      c 
c**************************************************************c

	    ! all things for pc
	integer*4 alloc_stat,start
      real(kind=r8_kind),pointer  :: grad(:) 
      real(kind=r8_kind),dimension(3):: con_vec,vec_1,vec_2,normal_vec
      real(kind=r8_kind),parameter,dimension(3) ::
     &  e1=(/1.0_r8_kind,0.0_r8_kind,0.0_r8_kind/),
     &  e2=(/0.0_r8_kind,1.0_r8_kind,0.0_r8_kind/),
     &  e3=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
        real(kind=r8_kind)             :: abs_e
    ! ----- executable code ------------------------------------ 

      real*8 jdc(nlocv,nlocv),s(nlocv,nlocv),b(nlocv,nlocv),ipv(nlocv,3)
	real*8 grs(nlocv)

	n3=0
	do i=1,kolat
	 if(impu(i).ne.0) n3=n3+3	
	enddo ! i=1,kolat

	jdc(:n3,:n3)=0.d0
c  calc. of matrix of partial derivatives dt/dr 
      mu=0
      do 100 n=1,3
      do 101 ip=1,kolat
      i=ip
      if(zmat(i,n).eq.0.or.impu(i).eq.0 ) goto 101 
      j=abs(zmat(i,1))
      k=abs(zmat(i,2))
      l=abs(zmat(i,3))
      mu=mu+1
C:TESTDEF>>
	if(i.eq.j) then
	write(6,*) 'bond variable is wrongly defined'
        write(6,*) i,j
	stop 'bond variable is wrongly defined'
	endif ! i.eq.j

        if(i.eq.k.or.j.eq.k) then
        write(6,*) 'ijk angle variable is wrongly defined'
        write(6,*) i,j,k
        stop 'ijk angle variable is wrongly defined'
        endif ! i.eq.k.or.j.eq.k

        select case(n)
        case(2)
        if(sin(t(i,j,k)).lt.1.d-5) then  ! delete variable block
        write(6,*) 'trans: ijk angle close to 180 '
        write(6,*) i,j,k
        stop 'angle variable is wrongly defined'
        endif ! dsin(t(i,j,k)).lt.1.d-5)
        case(3)
        if(sin(t(j,k,l)).lt.1.d-5) then  ! delete variable block
        write(6,*) 'angle variable is wrong defined'
        write(6,*) j,k,l
        stop 'angle variable is wrongly defined'
        endif ! dsin(t(i,j,k)).lt.1.d-5)
        end select

	if(mu.gt.nlocv) then
	write(6,*) 
     &       ' number of internal coordinates is limited by 300 now'
	stop
	endif ! mu.gt.300
C:TESTDEF<<
	iro=0
      do 102 iat=1,kolat
	if(impu(iat).ne.0) then
      do m=1,3
      irow=iro+m
      if(n.eq.1) jdc(irow,mu)=dr(i,j,iat,m)
      if(n.eq.2) jdc(irow,mu)=dt(i,j,k,iat,m)
      if(n.eq.3) jdc(irow,mu)=df(i,j,k,l,iat,m)
      enddo	! m=1,3
	iro=iro+3
	endif	! impu(iat).ne.0
  102 continue
  101 continue
  100 continue

        if(iwork.eq.9)  
     &   write(6,*) 'the total number of internal variables',mu

	if(mu.lt.n3-6) then
	write(6,*) 
     &     'WARNING: number of variables to define projection matrix'
        write(6,*) 
     &   '         is less then 3n-6 so you can obtain wrong results'
	endif ! mu.lt.n3-6

      kolvar=mu

      if(fixed_orientation) then
	   n_local=kolvar
	   call alloc_bmat(kolvar,kolat)
	   bmat=transpose(jdc)
	   constraint_mat(1:n_local,:) = bmat(1:n_local,:) 
       ! fixation of atom 'fixed_atom_1' to its original position --------
       j=(fixed_atom_1-1)*3+1
       do i = n_local+1,n_local+3
          ! erster Spaltenindex fuer atom i_atom ist j=(i_atom-1)*3+1
          constraint_mat(i,j)=one
          j=j+1
       enddo
	print *,'first atom done'
       ! -----------------------------------------------------------------  
	
       con_vec = atom(fixed_atom_2)%x - atom(fixed_atom_1)%x
       con_vec=con_vec/abs_value(con_vec)
       ! Using vectors e1=(1,0,0),e2=(0,1,0) and e3=(0,0,1)
       ! find two vectors orthogonal to con_vec with a little Schmidt-procedure---
       vec_1 = e1 - dot_product(e1,con_vec)*con_vec
       abs_e = abs_value(vec_1)
       if (abs_e<small) then 
       ! e1 was (anti)parallel to con_vec and thus cannont be used
          !                    which means that e2 and e3 should be ok
          ! test this hypothesis:
          abs_e = dot_product(e2,con_vec)
          if ( .not. (dot_product(e2,con_vec)<small .and. dot_product(e3,con_vec)<small) ) then
             call error_handler("invert_bmat: sth. fishy here")
          endif
          vec_1 = e2
          vec_2 = e3                                                            
        else ! this means e1 is not parallel (or antiparallel) to the connecvtion vector
          vec_1=vec_1/abs_e
          ! now try to construct vec_2 with e2
          vec_2 = e2 - dot_product(e2,con_vec)*con_vec - dot_product(e2,vec_1)*vec_1
          abs_e = abs_value(vec_2)
          if ( abs_e<small ) then ! this means that only e3 is left to make up vec_2
             vec_2 = e3 - dot_product(e3,con_vec)*con_vec - dot_product(e3,vec_1)*vec_1
             abs_e = abs_value(vec_2)
             if (abs_e < small ) then
                write(*,*)" invert_bmat: Vector VEC_2 orthogonal to the connection line"
                write(*,*)"              of atoms ",fixed_atom_1," and ",fixed_atom_2
                write(*,*)"              could not be constructed. Seek technical assistance"
                call error_handler(" ")
             endif
             vec_2 = vec_2/abs_e
          else !this means e2 was ok
             vec_2=vec_2/abs_e
          endif
          endif

       ! now find out column j of the b-matrix belonging to i_atom
       j=(fixed_atom_2-1)*3+1
       do i=1,3
          constraint_mat(n_local+4,j) = vec_1(i)
          constraint_mat(n_local+5,j) = vec_2(i)
          j=j+1
       enddo


       normal_vec = cross_product(atom(fixed_atom_2)%x
     &  -atom(fixed_atom_1)%x,atom(fixed_atom_3)%x-atom(fixed_atom_1)%x)
       abs_e = abs_value(normal_vec)
       if (abs_e<small) then
          write(*,*)" invert_bmat: the normal vector of the plane spanned by atoms"
          write(*,*)"             ",fixed_atom_1,", ",fixed_atom_2," and ",fixed_atom_3
          write(*,*)"              is rubbish. Seek technical assistance"
          call error_handler(" ")
       endif

       normal_vec=normal_vec/abs_e
       ! find out column for atom i_atom
       j=(fixed_atom_3-1)*3+1
       do i=1,3
          constraint_mat(n_local+6,j) = normal_vec(i)
          j=j+1
       enddo

       ! finally constraint_mat can be inverted using standard methods
       call invert_matrix(3*(kolat),constraint_mat)
       ! ... and the important stuff is copied to bmat_inv
       bmat_inv(:,:) = constraint_mat(:,1:n_local)
       deallocate(constraint_mat,STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler 
     &    ("invert_bmat: deallocation (4) failed")
	    allocate(grad(n_local),STAT=alloc_stat)
	    grad=0.d0
        do i=1,ubound(grad,1)
          do k=1,kolat
             start=(k-1)*3+1
             ende=k*3
             grad(i) = grad(i) + 
     &            bmat_inv(start,i)*gr(1,k) + 
     &            bmat_inv(start+1,i)*gr(2,k) + 
     &            bmat_inv(start+2,i)*gr(3,k)
          enddo
          if (abs(grad(i))<=1.0e-10_r8_kind ) grad(i)=0.d0
        enddo         
        write(*,*)" full internal gradient is:"
        do i=1,ubound(grad,1)
          write(*,'(I4,3x,f9.5)')i,grad(i)
        enddo 
	g(:kolvar)=grad
        
c$$$       endif                                                         
	else
                    
	   
c       construction of inverse matrix jdc
      do mu=1,kolvar
      do nu=1,kolvar
      s(mu,nu)=0.d0
      do i=1,n3
      s(mu,nu)=s(mu,nu) + jdc(i,mu)*jdc(i,nu)
	enddo
	enddo
	enddo

      call matinv(s,kolvar,det,nlocv,ipv)
      print*,'*** det=',det

      if(det.eq.0.d0) then
	print*,' internal variables are to be redefined ***'
      stop
	endif ! det.eq.0.d0

      do 400 is=1,kolvar
      do m=1,n3
      b(is,m)=0.d0
      do mu=1,kolvar
      b(is,m)=b(is,m) + s(is,mu)*jdc(m,mu)
	enddo	! mu=1,kolvar
	enddo	! m=1,n3
  400 continue

	do is=1,kolvar
	k=0
	do iat=1,kolat
	if(impu(iat).ne.0) then
	do kk=1,3
	k=k+1
      g(is)=g(is) + b(is,k)*gr(kk,iat) ! gradient in internal coord
	enddo ! kk=1,3
	endif ! impu(iat).ne.0
	enddo ! k=1,n3
        enddo ! is=1,kolvar

	endif ! fixed/else

      grs(:n3)=g(:n3)
      g(:n3)=0.0_r8_kind
 
c  		independent components of internal gradient
      numvar=0 ! variables counter
      mu=0     ! coord counter
      numxmx =0
      do 700 i=1,3
      do iat=1,kolat
	if(zmat(iat,i).eq.0.or.impu(iat).eq.0 ) cycle
      mu=mu+1  ! mu - total number of internal coordinates
      if(numx(iat,i).ne.0) then	 ! numx - matrix of changable vars
	if(numx(iat,i).gt.numxmx) numxmx=numx(iat,i)
      if(numx(iat,i).gt.numvar) then
      numvar=numvar+1  ! number of independent vars
      g(numvar)=grs(mu)
 	else	! numx(iat,i).le.numvar
 	g(numx(iat,i))=g(numx(iat,i))+grs(mu)*isg(mu)
        endif ! numx(iat,i).gt.numvar/ else
	endif ! numx(iat,i).ne.0
	enddo	! iat=1,kolat
  700 continue

C:TESTDEF>>
	if(numxmx.ne.numvar) then
	write(6,*) ' numx matrix is wrongly defined'
	write(6,*) ' max number is not eq to number of variables'
	stop ' numx matrix is wrongly defined'
	endif
C:TESTDEF<<

      kolvar=numvar

	write(6,*) ' gradients : '
      do i=1,kolvar
      write(6,13) i,g(i)
   13 format(i3,2x,f12.6,2x,f12.6)
	enddo

      return
      e n d

      subroutine update(iwork,kolvar,xo)
	use smlcom, only: kamax,nlocv
	use xyzsml, only: x,y,z,kolat
	use zzzsml
	use numsml
	use zarsml
	use symsml
	use simol_io, only: output,s_opt,xold
c	** update of cartesian coordinates **
      IMPLICIT real*8 (A-H,O-Z)

      integer ncon1(kamax),ncon2(kamax),ncon3(kamax)

      real*8
     . xo(nlocv)
     .,rij,dr,nx,ny,nz,dtet,dti,dtk,df,dfl,dfi,r,t,f,s
     .,rotmatr(3,3)


C:CHAINUPDT>>
	if(s_opt.ge.0.1d0.and.s_opt.lt.0.2d0) then	
        call chainupdt(iwork,kolvar,xo)
        write(6,*) ' update : iwork,kolvar',iwork,kolvar
        if(iwork-10.gt.kolvar*2) return
CAY--------- After update is finished restore orientation -----------
CAY     Now translate to old coordinates of first atom
        dxpr = x(1) - xold(1,1)
        dypr = y(1) - xold(2,1)
        dzpr = z(1) - xold(3,1)
        do ixnew=1,kolat
          x(ixnew) = x(ixnew) - dxpr
          y(ixnew) = y(ixnew) - dypr
          z(ixnew) = z(ixnew) - dzpr
        enddo
CAY     Calculate vector product of (1->2)new X (1->2)old
CAY     and then rotate
        if ( kolat .gt. 1 ) then
          x21n = x(2) - x(1)
          y21n = y(2) - y(1)
          z21n = z(2) - z(1)
          x21o = xold(1,2) - xold(1,1)
          y21o = xold(2,2) - xold(2,1)
          z21o = xold(3,2) - xold(3,1)
          xvp = y21n*z21o - z21n*y21o
          yvp = z21n*x21o - x21n*z21o
          zvp = x21n*y21o - y21n*x21o
          r21n=sqrt(x21n**2+y21n**2+z21n**2)
          r21o=sqrt(x21o**2+y21o**2+z21o**2)
          rvp=sqrt(xvp**2+yvp**2+zvp**2)
          if( rvp .gt. 1.d-7 ) then
            angle=asin(rvp/(r21o*r21n))
c            write(6,*)' Axes:',xvp/rvp,yvp/rvp,zvp/rvp
c            write(6,*)' Sin(a) , a =',rvp/(r21o*r21n),angle
            call rot(rotmatr,xvp/rvp,yvp/rvp,zvp/rvp,angle)
            call rotate(rotmatr,x,y,z,kolat)
          endif
c          write(6,*)' After 1st rotation'
c          do ixnew=1,kolat
c            write(6,'(3F15.8)')x(ixnew),y(ixnew),z(ixnew)
c          enddo
CAY       Find angle between planes (1-2-3)new & (1-2-3)old
          if ( kolat .gt. 2 ) then
            x32n = x(3) - x(2)
            y32n = y(3) - y(2)
            z32n = z(3) - z(2)
            x32o = xold(1,3) - xold(1,2)
            y32o = xold(2,3) - xold(2,2)
            z32o = xold(3,3) - xold(3,2)
            xvpn = y32n*z21o - z32n*y21o
            yvpn = z32n*x21o - x32n*z21o
            zvpn = x32n*y21o - y32n*x21o
            xvpo = y32o*z21o - z32o*y21o
            yvpo = z32o*x21o - x32o*z21o
            zvpo = x32o*y21o - y32o*x21o
            r21o = sqrt(x21o**2+y21o**2+z21o**2)
            r3n3o2=(xvpn-xvpo)**2+(yvpn-yvpo)**2+(zvpn-zvpo)**2
            rvpn2=xvpn**2+yvpn**2+zvpn**2
            rvpo2=xvpo**2+yvpo**2+zvpo**2
            if ( rvpn2 .le. 1.d-7 .or. rvpo2 .le. 1.d-7 ) goto 400
            coss=(rvpn2+rvpo2-r3n3o2)/(2*sqrt(rvpn2*rvpo2))
            if ( coss .ge. 1.d0 ) then
              angle=0.d0
            else if ( coss .le. -1.d0 ) then
              angle=3.14159265358979323846D0
            else
              angle = acos( coss )
            endif
CAY         And rotate now
c           write(6,*)' Axes:',x21o/r21o,y21o/r21o,z21o/r21o
c           write(6,*)' Cos(a) , a =',coss,angle
            call rot(rotmatr,x21o/r21o,y21o/r21o,z21o/r21o,angle)
            call rotate(rotmatr,x,y,z,kolat)
 400        continue
c            write(6,*)' After 2nd rotation'
c            do ixnew=1,kolat
c              write(6,'(3F15.8)')x(ixnew),y(ixnew),z(ixnew)
c            enddo
          endif
        endif
CAY--------------
	print*,'output(kolat,kolvar,iwork)'
	call output(kolat,kolvar,iwork) ! -> stop
	return
	endif ! s_opt.ge.0.1d0.and.and.s_opt.lt.0.2d0
C:CHAINUPDT<<


c
c *********************** bond updating  **********************
c
      do 1 iat=1,kolat ! cgx_ 2 -> 1
      if(numx(iat,1).ne.0.and.impu(iat).ne.0) then
      i=iat
      j=abs(zmat(i,1))
      rij=r(i,j)
      nx=(x(j)-x(i))/rij
      ny=(y(j)-y(i))/rij
      nz=(z(j)-z(i))/rij
      dr=( r(i,j) - xo(numx(i,1)) )/2
	dr1=dr
	dr2=dr

c define displacements of the atoms connected with the atom i
      call cnct1(i,j,ncon1,kolcon1,kolat)
      do  icon=1,kolcon1
      ncon1(kolcon1-icon+2)=ncon1(kolcon1-icon+1)
	enddo ! icon=1,kolcon1
      ncon1(1)=i
      kolcon1=kolcon1+1

c define displacements of the atoms conected with the atom j 
	kolcon2=0
    2   kolcon2=kolcon2+1
	ncon2(kolcon2)=j
	call cnct1(j,i,ncon2(kolcon2+1),kolcon,kolat)
      kolcon2=kolcon2+kolcon
        i=j
        j=abs(zmat(j,1))
        if(j.ne.0) go to 2

       if(zmat(iat,1).lt.0) then ! high priotity assymetric updating
         dr1=dr+dr
         dr2=0.d0
        endif ! zmat(iat,1).lt.0 

	if(i_tpr.eq.1) then
        write(6,*) 'bond updating:',dr1 ,dr2
        write(6,'(10i4)') (ncon1(ii),ii=1,kolcon1)
        write(6,'(10i4)') (ncon2(ii),ii=1,kolcon2)
	endif ! i_tpr

      do icon=1,kolcon1
      ii=ncon1(icon)
      x(ii)=dr1*nx + x(ii)
      y(ii)=dr1*ny + y(ii)
      z(ii)=dr1*nz + z(ii)
        enddo ! icon=1,kolcon1

      do jcon=1,kolcon2
      jj=ncon2(jcon)
      x(jj)=-dr2*nx + x(jj)
      y(jj)=-dr2*ny + y(jj)
      z(jj)=-dr2*nz + z(jj)
        enddo ! jcon=1,kolcon2

	endif ! numx(iat,1).ne.0.and.impu(iat).ne.0
    1 continue
c
c ****************** updating of angles ijk ******************
c
      do 10 iat=1,kolat ! cgx_ 3 -> 1
      if( numx(iat,2).eq.0.or.impu(iat).eq.0 ) go to 10
      i=iat
      j=abs(zmat(i,1))
      k=abs(zmat(i,2))
      dtet=t(i,j,k) - xo(numx(i,2))
      nx=(y(i)-y(j))*(z(k)-z(j)) - (y(k)-y(j))*(z(i)-z(j))
      ny=(x(k)-x(j))*(z(i)-z(j)) - (x(i)-x(j))*(z(k)-z(j))
      nz=(x(i)-x(j))*(y(k)-y(j)) - (x(k)-x(j))*(y(i)-y(j))
      s=sqrt( nx**2 + ny**2 + nz**2 )
      nx=nx/s
      ny=ny/s
      nz=nz/s
	if(i_tpr.eq.1) then
      write(6,7001)!nx,ny,nz,
     > i,j,k
 7001 format(! 'nx,ny,nz ',3(d12.6,1x)/
     > 'ijk angle ',3i2)
	endif ! i_tpr
c				    find atoms directly bond conected to i
      call cnct1(i,j,ncon3,kolcon3,kolat)
      call cnct2(i,j,k,ncon1,kolcon1,kolat)  !& angle connected to k
      do   icon=1,kolcon3
          do jcon=1,kolcon1
         if(ncon1(jcon).eq.ncon3(icon)) goto 11
        enddo ! jcon=1,kolcon1
	kolcon1 = kolcon1 + 1
       ncon1(kolcon1)=ncon3(icon)
   11 continue
      enddo ! icon=1,kolcon3

c       			find atoms directly bond connected to k
      call cnct1(k,j,ncon3,kolcon3,kolat)
      call cnct2(k,j,i,ncon2,kolcon2,kolat)  !& angle connected to k
      do   icon=1,kolcon3
	   do jcon=1,kolcon2
	  if(ncon2(jcon).eq.ncon3(icon)) goto 12
	 enddo ! jcon=1,kolcon2
        kolcon2 = kolcon2 +1
       ncon2(kolcon2)=ncon3(icon)
   12 continue
        enddo ! icon=1,kolcon3   

	
C:SUPDATE>>
	if(s_opt.ge.0.2d0) then
      dti=dtet/2.d0
      dtk=-dtet/2.d0

	if(zmat(iat,2).lt.0) then
	dti=dtet
        dtk=0.d0
	endif ! zmat(iat,2).lt.0

	else  ! s_opt.lt.0.2
        dti=dtet
        dtk=0.d0
        endif ! s_opt.ge.0.2/else
C:SUPDATE<<


      call move(i,j,k,dti,dtk,nx,ny,nz,ncon1,kolcon1,
     &            ncon2,kolcon2)
   10 continue
c
c ***************** updating of dihedral angles ijkl ************
c
      do 100 iat=1,kolat ! cgx_ 4 -> 1
      if( numx(iat,3).eq.0.or.impu(iat).eq.0 ) go to 100
      i=iat
      j=abs(zmat(i,1))
      k=abs(zmat(i,2))
      l=abs(zmat(i,3))
C:SIM>>
!	print*,i,j,k,l
!	print*,iszm(i),f(i,j,k,l),xo(numx(i,3))
      df=f(i,j,k,l) - xo(numx(i,3)) *iszm(i)
C:SIM<<
      nx=x(k)-x(j)
      ny=y(k)-y(j)
      nz=z(k)-z(j)
      s=sqrt( nx**2 + ny**2 + nz**2 )
      nx=nx/s
      ny=ny/s
      nz=nz/s
        if(i_tpr.eq.1) then 
	write(6,*) 'ijkl angle updating:',i,j,k,l
	endif ! i_tpr.eq.1

      call cnct1(i,j,ncon1,kolcon1,kolat) ! the down branch deff
	kolcon_p = kolcon1

       call cnct2(i,j,k,ncon1(kolcon1+1),kolcon1,kolat)  !& angle connected to k
 	kolcon1 = kolcon_p + kolcon1
 	kolcon_p = kolcon1

      call cnct3(i,k,j,i,ncon1(kolcon1+1),kolcon1,kolat)
	kolcon1 = kolcon_p + kolcon1
        kolcon_p = kolcon1

      call cnct3(i,j,k,l,ncon1(kolcon1+1),kolcon1,kolat) 
	kolcon1 = kolcon_p + kolcon1
       kolcon_p = kolcon1


      call cnct3(l,k,j,i,ncon2,kolcon2,kolat)  ! the up branch deff

C:SUPDATE>>
        if(s_opt.ge.0.2d0) then
      dfi=df/2.d0
      dfl=-df/2.d0
      if(zmat(iat,3).lt.0) then
	dfi=df
        dfl=0.d0
      endif ! zmat(iat,3).lt.0 

	else  ! s_opt.lt.0.2d0
        dfi=df
        dfl=0.d0
        endif ! s_opt.ge.0.2d0/else
C:SUPDATE<<

	if(i_tpr.eq.1) then 
      write(6,7009) dfi,dfl
 7009 format('dfi,dfl ',2f13.6)
	endif !   i_tpr.eq.1

      call move(i,j,l,dfi,dfl,nx,ny,nz,ncon1,kolcon1,
     .            ncon2,kolcon2)
  100 continue

C:CHECKUPDT>>
      do  iat=1,kolat ! cgx_ 2 -> 1
	if(impu(iat).ne.0) then
      i=iat
      j=abs(zmat(i,1))
      k=abs(zmat(i,2))
      l=abs(zmat(i,3))

      if(numx(iat,1).ne.0) then
      rij=r(i,j)
      dr= r(i,j) - xo(numx(i,1)) 
	if(abs(dr).gt.1.d-5) then
	Write(6,*) ' bond update is wrong ', numx(i,1),xo(numx(i,1)),r(i,j)
	stop
	endif !  (abs(dr).gt.1.d-5
        endif !  numx(iat,1).ne.0

      if( numx(iat,2).ne.0 ) then
      dtet=t(i,j,k) - xo(numx(i,2))
        if(abs(dtet).gt.1.d-5) then
        Write(6,*) 
     > ' ijk angle update is wrong', numx(i,2),xo(numx(i,2)),t(i,j,k)
	stop
        endif !  (abs(dtet).gt.1.d-5
        endif !  numx(iat,2).ne.0

      if( numx(iat,3).ne.0 ) then
      df=f(i,j,k,l)*iszm(i) - xo(numx(i,3))
        if(abs(df).gt.1.d-5) then
        Write(6,*)
     > ' ijkl angle update can be wrong'
     > , numx(i,3),xo(numx(i,3)), f(i,j,k,l)
        endif !  (abs(df).gt.1.d-5
        endif !  numx(iat,3).ne.0

	endif ! impu(iat).ne.0
	enddo ! iat=1,kolat

C:CHECKUPDT<<
	if(iwork.ne.0) write(6,*) ' update : iwork,kolvar',iwork,kolvar
C:FREQ>>
	if(iwork-10.gt.kolvar*2) return
C:FREQ<<
	if(iwork.ne.0) call output(kolat,kolvar,iwork) ! -> stop
	return
      e n d

      SUBROUTINE rotate(rmat,x,y,z,kolat)
      implicit real*8(a-h,o-z)
      dimension rmat(3,3),x(kolat),y(kolat),z(kolat)
      do i=1,kolat
        xt=rmat(1,1)*x(i)+rmat(1,2)*y(i)+rmat(1,3)*z(i)
        yt=rmat(2,1)*x(i)+rmat(2,2)*y(i)+rmat(2,3)*z(i)
        zt=rmat(3,1)*x(i)+rmat(3,2)*y(i)+rmat(3,3)*z(i)
        x(i)=xt
        y(i)=yt
        z(i)=zt
      enddo
      end
      SUBROUTINE rotateatom(rmat,x,y,z,i,j)
      implicit real*8(a-h,o-z)
      dimension rmat(3,3),x(1),y(1),z(1)
      xt=
     &rmat(1,1)*(x(i)-x(j))+rmat(1,2)*(y(i)-y(j))+rmat(1,3)*(z(i)-z(j))
      yt=
     &rmat(2,1)*(x(i)-x(j))+rmat(2,2)*(y(i)-y(j))+rmat(2,3)*(z(i)-z(j))
      zt=
     &rmat(3,1)*(x(i)-x(j))+rmat(3,2)*(y(i)-y(j))+rmat(3,3)*(z(i)-z(j))
      x(i)=xt+x(j)
      y(i)=yt+y(j)
      z(i)=zt+z(j)
      e n d

      SUBROUTINE chainupdt(iwork,kolvar,xo)
	use smlcom, only:nlocv
	use xyzsml, only: x,y,z,kolat
	use zzzsml
	use numsml
	use symsml
	use optsml, only: xoall
      IMPLICIT real*8 (A-H,O-Z)
!      common/xyz/x,y,z,gr,kolat
!      common/opt/dummy(600),xoall(300)
!      common/zzz/zmat,ieq(100)
!      common/num/numx,i_tpr
!      common/sym/isg(300),iszm(100)
      common/preserv/xold(3,3)
      real*8
     . xo(nlocv)
     .,rotmatr(3,3)
 
C     MAIN LOOP OVER ALL ATOMS
      do iat=1,kolat
        i = iat
        j = zmat( i, 1 )
        k = zmat( i, 2 )
        l = zmat( i, 3 )
        if ( i .eq. 1 ) then
          x(i) = xold(1,i)
          y(i) = xold(2,i)
          z(i) = xold(3,i)
        else
          rij = xoall(i-1)
C               ^^^^^^^^^^ this should be Rij
          theta = xoall( i + kolat - 3 )
C                 ^^^^^^^^^^^^^^^^^^^^^^ this should be Theta(i-j-k)
          phi = xoall( i + kolat*2 - 6 )
C               ^^^^^^^^^^^^^^^^^^^^^^^^ this should be Phi(i-j-k-l)
C In the case of optimized variable we use XO and NUMX arrays
          if (numx(i,1).ne.0) rij = xo(numx(i,1))
          if (numx(i,2).ne.0) theta =  xo(numx(i,2))
          if (numx(i,3).ne.0) phi = xo(numx(i,3))
 
C Now for i=2 we get atomic coords in the original direction
C and for i.ne.2 atom i is placed along axes j-k
          if ( i .eq. 2 ) then
            xijo = xold(1,i)-xold(1,j)
            yijo = xold(2,i)-xold(2,j)
            zijo = xold(3,i)-xold(3,j)
          else
            xijo = x(k) - x(j)
            yijo = y(k) - y(j)
            zijo = z(k) - z(j)
          endif
          rijo = sqrt( xijo**2 + yijo**2 + zijo**2 )
          xijo = xijo / rijo
          yijo = yijo / rijo
          zijo = zijo / rijo
          x(i) = x(j) + rij * xijo
          y(i) = y(j) + rij * yijo
          z(i) = z(j) + rij * zijo
c          write(6,'(A,I3,A,3I3,A,3F11.7)')
c     &    ' After 1st step: atom ',i,', zmat:',j,k,l
c     &   ,', coords: ',x(i),y(i),z(i)
          if ( i .eq. 2 ) goto 100  ! goto end of loop
 
C Now we have to rotate
C if i=3 we rotate bond i-j around atom j in the original plane i-j-k
C else we rotate this bond around atom j in the plane j-k-l
          if ( i .eq. 3 ) then
            x1 = -(xold(1,i) - xold(1,j))
            y1 = -(xold(2,i) - xold(2,j))
            z1 = -(xold(3,i) - xold(3,j))
            x2 = x(i) - x(j)
            y2 = y(i) - y(j)
            z2 = z(i) - z(j)
          else
            x1 = x(l) - x(k)
            y1 = y(l) - y(k)
            z1 = z(l) - z(k)
            x2 = x(j) - x(k)
            y2 = y(j) - y(k)
            z2 = z(j) - z(k)
          endif
          call vecmul(x1,y1,z1,x2,y2,z2,xp,yp,zp)
          rp = sqrt(xp*xp+yp*yp+zp*zp)
          if ( rp .lt. 1.e-7 ) then
            write(6,'(A/A)')
     &      ' !!! ERROR IN Z-MATRIX !!! WRONG ANGLE SPECIFIED.'
     &     ,' THETA EQUAL TO 0 OR 180 DEGREES IS NOT ALOWED HERE'
            stop 'Error'
          endif
          xp = xp / rp
          yp = yp / rp
          zp = zp / rp
          call rot(rotmatr,xp,yp,zp,theta)
          call rotateatom(rotmatr,x,y,z,i,j)
c          write(6,'(A,I3,A,3I3,A,3F11.7)')
c     &    ' After 2nd step: atom ',i,', zmat:',j,k,l
c     &   ,', coords: ',x(i),y(i),z(i)
          if ( i .eq. 2 ) goto 100  ! goto end of loop
 
C Now rotation around j-k axes on Phi-angle is left to do
          if ( i .gt. 3 ) then
            xp = x(k) - x(j)
            yp = y(k) - y(j)
            zp = z(k) - z(j)
            rp = sqrt(xp*xp+yp*yp+zp*zp)
            xp = xp / rp
            yp = yp / rp
            zp = zp / rp
            call rot(rotmatr,xp,yp,zp,-phi)
            call rotateatom(rotmatr,x,y,z,i,j)
c          write(6,'(A,I3,A,3I3,A,3F11.7)')
c     &    ' After 3rd step: atom ',i,', zmat:',j,k,l
c     &   ,', coords: ',x(i),y(i),z(i)
          endif
        endif
 100    continue
      enddo
      e n d

      subroutine vecmul(x1,y1,z1,x2,y2,z2,xp,yp,zp)
      IMPLICIT real*8 (A-H,O-Z)
      xp=y1*z2-z1*y2
      yp=z1*x2-x1*z2
      zp=x1*y2-y1*x2
      end

      subroutine cnct1(i1,j,ncon,kolcon,kolat)
	use zzzsml
	use zarsml
      IMPLICIT real*8 (A-H,O-Z)

c
c  determination of down bond connectivity set of i1 atom, exepting j-th
c
!      common/zzz/zmat,ieq(100)
!      common/zar/an(100),impu(100)
      integer ncon(*)
      i=i1
      kolcon=0
      kolall=0
      knew1=1
    5 knew2=0

      do k=1,knew1
      if(kolall.ne.0) i=ncon(kolall-knew1+k)

      do iat=2,kolat 
	if(impu(iat).ne.0) then
      if(abs(zmat(iat,1)).eq.i.and.iat.ne.j ) then
      knew2=knew2+1
      kolcon=kolcon+1
      ncon(kolcon)=iat
	endif ! zmat(iat,1).eq.i.iat.ne.j
	endif ! impu(iat).ne.0
	enddo ! iat=2,kolat

	enddo ! k=1,knew1

      if(knew2.ne.0) then
      knew1=knew2
      kolall=kolcon
      go to 5
	endif !knew2.ne.0
      return
      e n d

      subroutine cnct2(k1,j,i,ncon,kolcon,kolat)
	use zzzsml
	use numsml
	use zarsml
	implicit real*8 (a-h,o-z)
c**************************************************************c
c         determination of "connectivity set" of k-th atoms    c
c         through  j-th,  exepting i-th                        c 
c**************************************************************c

      integer ncon(*)
!      common/zzz/zmat,ieq(100)
!      common/num/numx,i_tpr
!      common/zar/an(100),impu(100)
      k=k1
      kolcon=0
      knew2=0

c       	** define  first neighbour jk down dependents **
      do iat=3,kolat 
	if(impu(iat).ne.0) then
      if(abs(zmat(iat,2)).eq.k.and.abs(zmat(iat,1)).eq.j.and.iat.ne.i)
     >  then
c	the current atom is conected to the cental atom and angle connected
c 	to the terminal atom k
      knew2=knew2+1  		! ** count number of branches **
      kolcon=kolcon+1
      ncon(kolcon)=iat
c       write(*,*) 'first neighbour kj down:',iat
	endif ! zmat(iat,2).eq.k .and. zmat(iat,1).eq.j.and.iat.ne.i
	endif ! impu(iat).ne.0
	enddo ! iat=3,kolat

c                   ** define k down bond dependents **
c                  (find all atoms on the branches)

    5 if(knew2.NE.0) THEN
      knew1=knew2
      kolall=kolcon
      knew2=0
 
      do  m=1,knew1
        k=ncon(kolall-knew1+m)
c	write(*,*) 'find atoms on the branches',k
         do  iat=1,kolat 
	 if(impu(iat).ne.0) then
c                                         __________________
        if(iat.ne.i.and.iat.ne.j.and.abs(zmat(iat,1)).eq.k) then
      knew2=knew2+1
      kolcon=kolcon+1
      ncon(kolcon)=iat
	endif ! iat.ne.i.and.zmat(iat,1).eq.k
	endif ! impu(iat).ne.0
	enddo ! iat=1,kolat
      enddo ! m=1,knew1
      go to 5
	ENDIF ! knew2.NE.0

c       ** define  first neighbour kj up dependent
c                  __________
        iat = abs(zmat(k1,2))
        if(iat.ne.0.and.iat.ne.i.and.iat.ne.j) then
      knew2=knew2+1
      kolcon=kolcon+1
      ncon(kolcon)=iat
c       write(*,*) 'first neighbour kj up:',iat

      call cnct1(iat,k1,ncon(kolcon+1),kolcn,kolat)  ! ???
	kolcon=kolcon+kolcn
	else
	goto  11
        endif ! iat.ne.0 .and.iat.ne.i

c       ** define  first neighbour k1-1 angle up->down dependents
        iat = abs(zmat(k1,2))
	jat = abs(zmat(iat,1))
   10   if(jat.ne.0.and.jat.ne.j.and.jat.ne.i) then
	do k=1,kolcon
	if(ncon(k).eq.jat) go to 11
	enddo ! k=1,kolcon
      kolcon=kolcon+1
      ncon(kolcon)=jat
      call cnct1(jat,iat,ncon(kolcon+1),kolcn,kolat)
	colcon=kolcon+kolcn
	iat=jat
	jat=abs(zmat(jat,1))
	go to 10
        endif ! jat.ne.0.and.jat.ne.j
 
c       ** define  first neighbour k1-1 bond up->down dependents
   11   jat = abs(zmat(k1,1))
        iat = k1
   20   if(jat.ne.0.and.jat.ne.j) then
        do k=1,kolcon 
        if(ncon(k).eq.jat) go to 21 
        enddo ! k=1,kolcon 

      kolcon=kolcon+1   
      ncon(kolcon)=jat
      call cnct1(jat,iat,ncon(kolcon+1),kolcn,kolat) 
        colcon=kolcon+kolcn
        iat=jat
        jat=abs(zmat(jat,1))
        go to 20 
        endif ! jat.ne.0 .jat.ne.j
   21 return
      e n d

      subroutine cnct3(l1,k1,j,i,ncon,kolcon,kolat)
	use zzzsml
	use numsml
	use zarsml
c**************************************************************c
c         detemination of "connectivity set" of l-th atom      c 
c          trough jk-axes, exepting i-th                       c 
c**************************************************************c

      IMPLICIT real*8 (A-H,O-Z)
!      common/zzz/zmat,ieq(100)
!      common/num/numx,i_tpr
!      common/zar/an(100),impu(100)
      integer ncon(*)
      logical rep
      l=l1
      k=k1
      kolcon=0
      knew2=0
      do  iat=1,kolat
	if(impu(iat).ne.0) then
	if(zmat(iat,3).ne.0) then
        rep = abs(zmat(iat,3)).eq.l 
     >			.and. abs(zmat(iat,2)).eq.k
     >				.and. abs(zmat(iat,1)).eq.j
        if( rep .and. iat.ne.i) then
        knew2=knew2+1
        kolcon=kolcon+1
        ncon(kolcon)=iat
       endif ! rep .and. zmat(iat,1).eq.j.and.iat.ne.i
      endif !zmat(iat,3).ne.0
	endif ! impu(iat).ne.0
      enddo ! iat=4,kolat

    5 if(knew2.ne.0) then
      knew1=knew2
      kolall=kolcon
      knew2=0
      do m=1,knew1
      k=ncon(kolall-knew1+m)
      do iat=1,kolat
	if(impu(iat).ne.0) then
      if(iat.ne.i.and.abs(zmat(iat,1)).eq.k) then
      knew2=knew2+1
      kolcon=kolcon+1
      ncon(kolcon)=iat
	endif ! iat.ne.i.and.zmat(iat,1).eq.k
	endif ! impu(iat).ne.0
	enddo ! iat=1,kolat
        enddo ! m=1,knew1
      go to 5
	endif ! knew2.ne.0
      return
      e n d
c**************************************************************c
c       rotation    i,k and "connectivity sets"                c
c**************************************************************c
      subroutine move(i,j,k,di,dk,nx,ny,nz,ncon1,kolcon1,
     .                ncon2,kolcon2)
	use xyzsml, only: x,y,z,kolat
	use numsml
	implicit real*8 (a-h,o-z)
!      common/num/numx,i_tpr
      integer  ncon1(*),ncon2(*)
!       common/xyz/x(100),y(100),z(100),gr(300),kolat
      real*8  nx,ny,nz,u(3,3)

c ********************** i ****************************
      call rot(u,nx,ny,nz,di)
      do 51 icon=1,kolcon1
      ncon1(kolcon1-icon+2)=ncon1(kolcon1-icon+1)
   51 continue
      ncon1(1)=i
      kolcon1=kolcon1+1
      do 1 ii=1,kolcon1
      iat=ncon1(ii)
      dx=x(iat)-x(j)
      dy=y(iat)-y(j)
      dz=z(iat)-z(j)
      x(iat)=u(1,1)*dx + u(1,2)*dy + u(1,3)*dz  + x(j)
      y(iat)=u(2,1)*dx + u(2,2)*dy + u(2,3)*dz  + y(j)
      z(iat)=u(3,1)*dx + u(3,2)*dy + u(3,3)*dz  + z(j)
    1 continue

c ********************** k ****************************
      call rot(u,nx,ny,nz,dk)
      do 52 icon=1,kolcon2
      ncon2(kolcon2-icon+2)=ncon2(kolcon2-icon+1)
   52 continue
      ncon2(1)=k
      kolcon2=kolcon2+1
      do 2 kk=1,kolcon2
      iat=ncon2(kk)
      dx=x(iat)-x(j)
      dy=y(iat)-y(j)
      dz=z(iat)-z(j)
      x(iat)=u(1,1)*dx + u(1,2)*dy + u(1,3)*dz  + x(j)
      y(iat)=u(2,1)*dx + u(2,2)*dy + u(2,3)*dz  + y(j)
      z(iat)=u(3,1)*dx + u(3,2)*dy + u(3,3)*dz  + z(j)
    2 continue

        if(i_tpr.eq.1) then  
        write(6,*) 'di,dk',di,dk 
      write(6,2001) (ncon1(m),m=1,kolcon1) 
 2001 format('ncon1 ',20i2) 
      write(6,2002) (ncon2(m),m=1,kolcon2)
 2002 format('ncon2 ',20i2)
        endif ! i_tpr.eq.1 
      return
      e n d

c********************************************c
c  rotation around axes (nx,ny,nz) on d      c
c********************************************c
      subroutine rot(u,nx,ny,nz,d)
      IMPLICIT real*8(A-H,O-Z)
      real*8  u(3,3),nx,ny,nz,d,cosd,sind

      cosd=cos(d)
      sind=sin(d)
      u(1,1)=cosd + (1.d0-cosd)*nx*nx
      u(2,1)=(1.d0-cosd)*nx*ny + sind*nz
      u(3,1)=(1.d0-cosd)*nx*nz - sind*ny
      u(1,2)=(1.d0-cosd)*nx*ny -sind*nz
      u(2,2)=cosd + (1.d0-cosd)*ny*ny
      u(3,2)=(1.d0-cosd)*ny*nz + sind*nx
      u(1,3)=(1.d0-cosd)*nx*nz + sind*ny
      u(2,3)=(1.d0-cosd)*ny*nz - sind*nx
      u(3,3)=cosd + (1.d0-cosd)*nz*nz
      return
      e n d
      FUNCTION OVERLP(IOV,VEC1,VEC2,SCALE,NVAR)
      IMPLICIT real*8(A-H,O-Z)
      DIMENSION VEC1(NVAR),VEC2(NVAR),SCALE(NVAR)
C     THE DOT PRODUCT OF VEC1 WITH VEC2 IF IOV = 1
C     THE DOT PRODUCT SCALED BY SCALE IF   IOV = 2
      OVERLP = 0.D0
      IF (IOV.EQ.1) THEN
        DO 10 I = 1,NVAR
        OVERLP = OVERLP + VEC1(I)*VEC2(I)
 10     CONTINUE
      ELSEIF (IOV.EQ.2) THEN
        DO 20 I = 1,NVAR
        OVERLP = OVERLP + VEC1(I)*VEC2(I)*SCALE(I)
 20     CONTINUE
      ELSE
        WRITE (6,'(A)') '**** ERROR IN OVERLP *****'
        stop
      ENDIF
      RETURN
      E N D
      SUBROUTINE DDCOPY(A,B,LENGTH)
      real*8
     . A(LENGTH),B(LENGTH)
      DO 100 I=1,LENGTH
100   A(I)=B(I)
      RETURN
      E N D
      SUBROUTINE NEWPNT (sc_fac,STEPP,XOLD,XNEW,NVAR)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION STEPP(NVAR),XOLD(NVAR),XNEW(NVAR)
      DO 10 I = 1,NVAR
      XNEW(I) = XOLD(I) + STEPP(I)*sc_fac
 10   CONTINUE
      RETURN
      E N D
