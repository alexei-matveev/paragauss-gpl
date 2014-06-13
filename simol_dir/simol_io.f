           module simol_io
!		tapes	*****

!	2 -  tape with regular positions of cluster centers
!	9 -line search tape
!	14 -  tape with coordinates of EPE shells
c	keys 
c	:tsloc>>
c------------------------------------------------------------
c 		** IWORK spesification  **
c       IWORK = 0 checking internal variable values
c       IWORK = 1 start of geometry optimization
c       IWORK = 5 TS location
c       IWORK = 6 point charges genaration
c       IWORK = 7 generation of xmol data
c       IWORK = 8 update of geometry to equilibrium after 
c				frequency calc (IWORK = 10)
c       IWORK = 9 one point only calculation
c	IWORK > 9 calculation of vibration frequencies
c       IWORK = 10 is used to invoke frequency calculation
c------------------------------------------------------------
	use smlcom
	use xyzsml, only:x,y,z,gra
	use zzzsml
	use numsml
	use zarsml
	use symsml
	use optsml
	use flpocm, only : etot=>fiii
	use slspar_module

	implicit real*8 (a-h,o-z)
	integer*4,dimension(2):: balance_atom_nbs 
	integer*4::n_balance_atoms
!	integer, parameter :: r8_kind = selected_real_kind(15)

	real(r8_kind),parameter:: gepe_const=evau*angsau**2
	real(r8_kind)::small_dist=0.05,mn,mx
	integer :: epe_lev
	type(gx_reg), dimension(:),allocatable:: gx
	logical::epe_pc_shells=.false.
        logical::gx_pos_regular=.false.	!regular pos of gx-file
        logical::epe_pos_regular=.false.!regular pos of epe lattice
	logical::sls_parameters=.false.

      dimension xx(kamax),yy(kamax),zz(kamax)
      dimension gud(nlocv*2)			! up & down gradients
!						! for hess and freq calc.
        real*8::s_opt
	real*8, dimension(3):: gtot
	real*8, dimension(3,3):: xold(3,3) ! x,y,z for atoms 1,2,3 on input
        real*8 eyead
	logical l_calc_h   ! true to calc hessian for the poin directly
	integer :: jopx,icycle,mode

C!!!      common       /steps/ pmste

	logical ch_xo,ny_updt

      real*8
     &  jdc(nlocv**2) ! projection matrix
     & ,d2e(nlocv**2) ! force constant matrix
     & ,EigVecR(nlocv**2) !  Vibrations Eigenvectors matrix
C>>:use of freq.f
     & ,GG(nlocv**2)
     & ,FF(nlocv**2)
C<<:use of freq.f 
cC>>:use of freq1.f
c     > ,fgm(nlocv**2) ! f*g matrix
cC<<:use of freq1.f
     & ,gr_eqvl(nlocv) ! stationary point gradients
     & ,DA(nlocv)      	! TS vector
     & ,F_TSLOC(nlocv)      ! TS vector

	dimension 
     >  xo_eqvl(nlocv)      ! eqvilibrium values of internal coordinates
C:PCFIELD>>
     >, vpch(kamax)	  ! potentials due to point charges
     >, qpch(30)  	  ! point charges
     >, fpch(3,kamax)       ! forces due to point charges
     >, fdip(3,kamax)       ! forces due to point charges
C:PCFIELD<<

        contains
         subroutine simol_inp(kolat,kolvar,iwork,dell,n_00)
         use slspar_module, only: gradients_3_body,atom,fixed_orientation 
	character*100 :: local_gsn
	real(r8_kind),dimension(3):: pos_temp
	 integer, intent(out):: kolat
	 integer, dimension(5)::ind_deck
	 integer*4, intent(out)::n_00
         integer(i4_kind):: status
         real*8 dell(*)
	ieq_mx = 0
	print *,'simol_inp entered'
	ny_updt =.false. ! to be true to change geometry by hands
	kolat=0
	ii=0
2000	ii=ii+1
c************  buf   ***************
         read(iu152,*)
     *   an(ii)
     *  ,x(ii),y(ii),z(ii)
     *  ,ieq(ii),impu(ii)
     *  ,zmat(ii,1),zmat(ii,2),zmat(ii,3)
     *  ,numx(ii,1),numx(ii,2),numx(ii,3),iepe(ii)
  203 format((f5.2,3(2X,f13.7),2i4,2x,3I3,2X,3I3,2x,2i3)) ! LCGTO format
	if(abs(ieq(ii)).gt.ieq_mx) ieq_mx = abs(ieq(ii))
2011    if(an(ii).gt.0.5d0) go to 2000
        kolat=ii-1
	if(kolat.eq.0) then
	write(6,*) 'geometry file is empty'
        stop 
	endif ! ii.eq.0
        iwork=-an(ii)+0.01d0


	l_calc_h = ieq(ii).gt.0  !c:tsloc
	jopx = ieq(ii)
	icycle = impu(ii)
	mode= zmat(ii,1)-1

	sls_forces=iepe(ii).ne.0
	sls_parameters=iepe(ii).gt.0
	call init_slsp()	!initialize known pair potentials
	if(sls_parameters) call make_epe_namelist
	if(iwork.eq.9.and.sls_parameters) call list_epe_par
        if(fixed_orientation) then
           allocate(atom(kolat),stat=status)
           if(status.ne.0) call error_handler 
     & (' allocation of atom failed')
           atom(:)%x(1)=x(:kolat)
           atom(:)%x(2)=y(:kolat)
           atom(:)%x(3)=z(:kolat)
           print *,'done allocate'
           endif


	s_opt = x(ii)
	eyead = y(ii) !c: cvparms
	if(eyead.lt.1.d-8) eyead = 1.d-4
	pmste = z(ii) !c: steps
	if(pmste.lt.1.d-8) pmste = 0.11d0

	selectcase(iwork)
	case(0)
!	save regular positions of ions for use in hades interface
!	to make proper relations between atoms
	 open(2, file='epe.r',status='new',form='formatted')
	 print*, 'file epe.r created'
	do i=1,kolat
	 if(iepe(i).ne.0) write(2,*) x(i),y(i),z(i) 
	enddo !i=1,kolat 
	close (2)
	case(7)

	if(zmat(ii,2).gt.0) then
	sc_coeff=1.d0
	elseif (zmat(ii,2).lt.0) then
	read (*,*) sc_coeff
	else
	sc_coeff=angsau
	endif
	call xmol_o(kolat,an,x,y,z,sc_coeff) ! -> xmol
	case(6)
	call pcgen(ii,kolat,numx,zmat,x,y,z,an,ieq,kamax)
        case(8)
        write(6,*) '** 8 - update to equilibrium geometry **'
	case(10)
	if(l_calc_h) then
	call hess_open  ! C:TSLOC
	else
	call freq_open   ! C:FREQ
	endif ! l_calc_h/else
	endselect

C:PCFIELD>>
	PCFIELD = s_opt*1.d5-int(s_opt*1.d4+.01)*10.d0
	call pcfield_cont()
C:PCFIELD<<
	


CAY----------- Orientation Preserving -----------------
      do ixold=1,min(3,kolat)
        xold(1,ixold) = x(ixold)
        xold(2,ixold) = y(ixold)
        xold(3,ixold) = z(ixold)
      enddo
CAY-----------


	print*,iwork,' iwork before entering read energy block'
	if(iwork.gt.0.and.iwork.lt.8) then

	selectcase(iwork)
	case(1)
!	call getenv("local_GSN_file",local_gsn)
!	open(55,file=trim(local_gsn),status='unknown')
!	close(55,status='delete')
!	call getenv("GSN_file", local_gsn)
!	open(55,file=trim(local_gsn),status='unknown')
!	close(55,status='delete')

	write(6,*) '----------------------------------------'
        write(6,*) '    ** Geometry optimization run **'
!        read(iu152,'(2f20.7,3i5)') etot,Et_BLYP,incod_gx,nygrd,icrel
        read(iu152,*) etot,Et_BLYP,incod_gx,nygrd,icrel
	call calc_type(incod_gx,icrel)

	casedefault
        read(iu152,*) etot,Et_BLYP
        write(6,*) ' The DFT energies for the point are:' 
	endselect

        write(6,*) etot,Et_BLYP

        endif !iwork.ne.0.and.iwork.lt.9



      write(6,*) ' the values of internal independent coordinates '
      kolvar=0
	mu=0 		     ! number of internel coordinates
      do  istr=1,3
      do 200 iat=1,kolat

	if(impu(iat).ne.0) then

        if(numx(iat,istr).lt.0) then ! change this coordinate
        numx(iat,istr)= -numx(iat,istr)
	ny_updt = .true.
	ch_xo = .true.
	else
	ch_xo = .false.
        endif ! numx(iat,istr).lt.0

      i=iat
      j=abs(zmat(i,1))
      k=abs(zmat(i,2))
      l=abs(zmat(i,3))

c ** independent variable are considered so numx control in action
      n=numx(iat,istr)
      go to (11,12,13),istr
c 		******************* bonds ****************
   11   if(j.ne.0) mu=mu+1
        isg(mu)=1
	mu_type(mu)=istr
	if(kolvar.ge.n.or.n.eq.0) go to 200
	xo(n)=r(i,j)
      kolvar=kolvar+1
	xp=xo(n)*angsau
	write (6,111) kolvar,xo(n),xp
        if(iwork.eq.10) then
	sxo(n)=0.05d0   ! c:freq
	write(6,*) ' step is ',sxo(n),sxo(n)*angsau
	else
	dell(n) = 0.005d0
	endif ! iwork.eq.10
	if(ch_xo) then
	read(iu152,*) xp
	if(iwork.eq.10) then
	sxo(n)=xp/angsau
        write(6,*) 'changed :' ,sxo(n),xp
	else
	xo(n) = xp/angsau
	write(6,*) 'changed :' ,xo(n),xp
	endif !iwork.eq.10/else
	endif ! ch_xo
      go to 200

   12  if(k.ne.0.and.j.ne.0) then 
c          ****************** ijk angles *****************
	mu=mu+1
        isg(mu)=1 
        mu_type(mu)=istr
	if(kolvar.ge.n.or.n.eq.0) go to 200
	xo(n)=t(i,j,k)

        if(sin(t(i,j,k)).lt.1.d-5) then
	write(6,*) 'input : ijk angle variable is wrong(too close to 180)'
	write(6,*) i,j,k
	stop 
	endif !sin(t(i,j,k)).lt.1.d-5)

      kolvar=kolvar+1
	xp=xo(n)*57.2957795130d0
 	write (6,111) kolvar,xo(n),xp

        if(iwork.eq.10)  then ! define step by ijk angle 
        sxo(n)=0.3d0/r(i,j)
	sxp=sxo(n)*57.29578
        if ((xp+sxp).ge.180.0d0) then
                sxp=179.9d0-xp
                sxo(n)=sxp/57.29578d0
        endif ! (xp+sxp).ge.180.0d0

        write(6,*) ' default step is ',sxo(n),sxo(n)*57.2957795130d0

	else	!iwork.ne.10  ** initialization step
	dell(n) = 0.03d0/sqrt(r(i,j)**3)
        endif ! iwork.eq.10

        if(ch_xo) then
	if(iwork.eq.10) then
	read(iu152,*) sxo(n)
        write(6,*) 'step is read :' ,sxo(n),sxo(n)*57.2957795130d0
        else ! iwork.ne.10 ! ch_xo
        read(iu152,*) xp
        xo(n) = xp/57.2957795130d0
        write(6,*) ' variable changed :' ,xo(n),xp
	        endif ! iwork.eq.10/else
        endif ! ch_xo

	endif ! k.ne.0.and.j.ne.0
      go to 200

c ************* dihedral angles  ***********
   13   if(l.ne.0.and.k.ne.0.and.j.ne.0) then !else the coordinate was deleted
        mu=mu+1  
        mu_type(mu)=istr
	if(n.eq.0) go to 200
	if(kolvar.ge.n) then
      if(abs(f(i,j,k,l)-xo(n)).lt.abs(xo(n))) then
        isg(mu)=1 
        iszm(i)=1
        else
	iszm(i)=-1
        isg(mu)=-1 
        endif ! abs(r(i,j)-xo(n)).lt.1.d-5 
        go to 200
	endif !kolvar.ge.n

	xo(n)=f(i,j,k,l)
	isg(n)=1
        iszm(i)=1
      kolvar=kolvar+1

	if(iwork.eq.10) then
	sxo(n)=0.010d0
	else
	dell(n)=0.005d0
	endif  ! iwork.eq.10/else

        if(iwork.eq.10.and.ch_xo) then
        read(iu152,*) xp
        sxo(n) = xp
        write(6,*) ' step is read :  ' ,sxo(n),xp
        endif ! ch_xo 

	xp=xo(n)*57.2957795130d0
 	write (6,111) kolvar,xo(n),xp
 111	format(2x,'xo(',i2,') ',2f12.6)
        if(ch_xo) then 
        read(iu152,*) xp 
        xo(n) = xp/57.2957795130d0
        write(6,*) 'changed :' ,xo(n),xp 
        endif ! ch_xo 
        endif ! l.ne.0.and.k.ne.0.and.j.ne.0
	endif ! impu(iat).ne.0
  200 continue
      enddo    ! istr=1,3
	print*, 'done input internal coordinates'

	if((iwork.eq.1.or.iwork.eq.9.or.iwork.eq.4)
     &                           .and..not.ny_updt) then
	call scale_steps(kolvar,dell,xo_eqvl,xo,kolat
     &                         ,xx,yy,zz)
	endif
	call def_i_tpr(iwork,i_tpr,kolvar)


	if(iwork.eq.0) stop
        if(iwork.eq.8) call eq_update(xo,kolvar)

        if(ny_updt.and.iwork.eq.9) call update(9,kolvar,xo)


        if(iwork.ge.9) then  ! **    frequency calc **

	if(iwork.gt.kolvar*2+10) return
c$$$!        read(iu152,'(2f20.7,3i5)') etot,Et_BLYP,incod_gx,nygrd,icrel
        read(iu152,*) etot,Et_BLYP,incod_gx,nygrd,icrel
	call calc_type(incod_gx,icrel)
        write(6,*) etot,Et_BLYP
	if(nygrd.eq.0) stop 

	call pc_fileld()
	endif !iwork.ge.9

!	now cartesian gradients are treated
	do n=1,kolat ! kolat number of centers inclusing dummy ones
	if(ieq(n).ne.0) read(iu152,*) iiii, gra(1:3,n)
	enddo ! n=1,kolat

	if(sls_forces) then
	inquire (file='epe.pcs',exist=epe_pc_shells)
	if(epe_pc_shells) then
	 print*,'file of EPE shells is found'
	open (14,file='epe.pcs',status='old',form='formatted')
        read (14,*) epe_nucen
	print*,'number of EPE shells in epe.pcs is equal to '
     &        ,epe_nucen
	allocate (epe(epe_nucen+kolat),stat=istat)
	if(istat.ne.0)  stop 'simol_inp: epe allocation failed'
	do i=1, epe_nucen
	read (14,*) epe(i+kolat)%s,z_dummy,ind_deck(1:4), 
     &              epe(i+kolat)%k,epe(i+kolat)%ant
	enddo ! i=1, nucen
	close(14)

	inquire (file='epe.r',exist=gx_pos_regular)
	if(gx_pos_regular) then
	print*, 'file with gx regular positions is found'
	open(2,file='epe.r',status='old',form='formatted')
	allocate(gx(kolat),stat=istat)
	if(istat.ne.0) stop 'simol_inp: gx_regular_pos alloc failed'
	gx(:)%ieq=ieq(1:kolat)
	do i=1,kolat
	 if(iepe(i).ne.0) read(2,*) gx(i)%r
	enddo ! i=1,kolat
	close(2)
	endif ! gx_pos_regular
        
	inquire (file='epe.pcr',exist=epe_pos_regular)
	if(epe_pos_regular) then
	print*, 'file with EPE regular positions epe.r is found'
	open(14,file='epe.pcr',status='old',form='formatted')
	read(14,*) nucen_pcr
	if(nucen_pcr.ne.epe_nucen) stop 'simol_inp: nucen_pcr.ne.epe_nucen'
	do i=1,nucen_pcr
	read(14,*) epe(i+kolat)%r
	enddo ! nucen_pcr
	else
	print*,' Warning: No file with regular positions'
	endif ! epe_pos_regular
           epe(1:kolat)%ant=an(1:kolat)
           epe(1:kolat)%s(1)=x(1:kolat)
           epe(1:kolat)%s(2)=y(1:kolat)
           epe(1:kolat)%s(3)=z(1:kolat)
	if(epe_pos_regular) then
	print*, 'check if regular positions of EPE and gx file coincide'
	n=1
	kl=0
	do while(n.le.epe_nucen)
	kl1=kl+1
	do k=kl1,kolat
	if(gx(k)%ieq.eq.0.or.iepe(k).eq.0) cycle
	if(dot_product(gx(k)%r-epe(n+kolat)%r,
     &                 gx(k)%r-epe(n+kolat)%r).lt.small_dist) then
	kl=kl+1
	pos_temp(:)=gx(kl)%r
	ieqt=gx(kl)%ieq
	gx(kl)%r=gx(k)%r
	gx(kl)%ieq=gx(k)%ieq
	gx(k)%r=pos_temp(:)
	gx(k)%ieq=ieqt

	print*,kl,' kl', gx(kl)%ieq
	pos_temp(:)=epe(kl+kolat)%r(:)
	ant=epe(kl)%ant
	epe(kl+kolat)%r(:)=epe(n+kolat)%r(:)
	epe(kl+kolat)%ant=epe(n+kolat)%ant
	epe(n+kolat)%r(:)=pos_temp(:)
	epe(n+kolat)%ant=ant
	pos_temp(:)=epe(kl+kolat)%s(:)
	epe(kl+kolat)%s(:)=epe(n+kolat)%s(:)
	epe(n+kolat)%s(:)=pos_temp(:)
	exit
	endif ! small_dist
	enddo !kolat
	n=n+1
	enddo ! i=1,nucen
	epe_kl=kl+1+kolat
	print*, 'done'
	if(kl.ne.0) then
	print*, 'number of EPE coinciding with centers in cluster ',kl
	else
	print*, 'no coinciding EPE with centers in cluster'
	endif
	endif

           if(n_types_central_atoms_3body.ne.0) then
           ks: do k_counter=1,kolat+epe_nucen
              do n_counter=1,max_type_ions
        if(abs(epe(k_counter)%ant-name_of_type(n_counter))
     &                      .lt.0.00001_r8_kind)  then
                    epe(k_counter)%k=n_counter

                    exit
                 end if
              end do

           end do ks
           endif


  
	call add_epe_f(kolat) ! i.e. due to outside epe links
        if(n_types_central_atoms_3body > 0) then
           call building_tet(kolat)
        etot=etot+energy_3_body(kolat)
        call gradients_3_body(kolat)           
           endif
	endif ! epe_pc_shells
	call add_sls_f(kolat) ! treate intercluster contributions
        print*, etot,' EPE correcred total energy'
	print*,'EPE correcred Cartesian gradients'
	grad_norm=0.0_r8_kind
	do i=1,kolat
	 print*,gra(:,i)
	grad_norm=grad_norm + dot_product(gra(:,i),gra(:,i))
	enddo ! i=1,kolat
	print*,'corrected gradient norm', grad_norm
	
	endif ! sls_forces

	n_00=0    ! c:nm_gx
	gtot(:)=0.d0
	n_balance_atoms=0
	do n=1,kolat ! kolat number of centers inclusing dummy ones
	if(impu(n).lt.0) then
	n_balance_atoms=n_balance_atoms+1
	balance_atom_nbs(n_balance_atoms)=n
	n_00=n  ! c:nm_gx
	endif ! impu(n).lt.0
	gtot(:)=gtot(:)+gra(:,n)
	if(ieq(n).le.0.or.impu(n).le.0) then ! nmv_gx
	gra(1:3,n) = 0.d0 !gradients are eq 0 for dummy atoms
	endif ! ieq(n).lt.0.or.an(n).gt.98.5
	enddo ! n=1,kolat
	if(iwork.eq.9) write(6,*) 'gtot ',gtot
	print*, 'number of balancing atoms' ,n_balance_atoms
	if(n_00.ne.0)  then
	if(n_balance_atoms.eq.1) then
	call sf_zero(gra,kolat,n_00)  ! c:NM_GX
	elseif(n_balance_atoms.eq.2) then
	call gett_balgrad(kolat)
	else
	stop 'wrong number of the balancing atoms'
	endif ! n_balance_atoms.eq.1/else/else
	endif! n_00.ne.0

	contains
	subroutine pcfield_cont
	if(PCFIELD.gt.0.1d0) then
        write(6,*) 'PCFIELD:', PCFIELD
	read(iu152,*) (qpch(k),k=1,ieq_mx)
	write(6,*) (qpch(k),k=1,ieq_mx)
	do i = 1, kolat

	ri = sqrt(x(i)**2+y(i)**2+z(i)**2)
	rix= x(i)/ri
	riy= y(i)/ri
	riz= z(i)/ri
        vpch(i) = 0.d0
	do k=1,3
	fpch(k,i) = 0.d0
        fdip(k,i) = 0.d0
	enddo ! k=1,3
	if(i.ne.1) then
	do j = 1, kolat
	if(i.ne.j.and.ieq(j).ne.0) then
	vpch(i) = vpch(i) - qpch(abs(ieq(j)))/r(i,j)
	fpch(1,i) = fpch(1,i) - (x(j)-x(i))*qpch(abs(ieq(j)))/r(i,j)**3
        fpch(2,i) = fpch(2,i) - (y(j)-y(i))*qpch(abs(ieq(j)))/r(i,j)**3
        fpch(3,i) = fpch(3,i) - (z(j)-z(i))*qpch(abs(ieq(j)))/r(i,j)**3
	cost=(rix*(x(j)-x(i))+riy*(y(j)-y(i))+riz*(z(j)-z(i)))/r(i,j)
        fdip(1,i) = fdip(1,i) 
     >             + 3.d0*cost*(x(j)-x(i))*qpch(abs(ieq(j)))/r(i,j)**4
        fdip(2,i) = fdip(2,i) 
     >             + 3.d0*cost*(y(j)-y(i))*qpch(abs(ieq(j)))/r(i,j)**4
        fdip(3,i) = fdip(3,i) 
     >             + 3.d0*cost*(z(j)-z(i))*qpch(abs(ieq(j)))/r(i,j)**4
	fdip(1,i) = fdip(1,i)- rix*qpch(abs(ieq(j)))/r(i,j)**3
	fdip(2,i) = fdip(2,i)- riy*qpch(abs(ieq(j)))/r(i,j)**3
	fdip(3,i) = fdip(3,i)- riz*qpch(abs(ieq(j)))/r(i,j)**3
	endif ! i.ne.j
	enddo ! j = 1, kolat
	if(PCFIELD.lt.1.1) write(6,*) i, vpch(i)
	if(i.eq.5.or.i.eq.9) write(6,*) (fpch(k,i),k=1,3)
 	if(i.eq.5.or.i.eq.9) write(6,*) (fdip(k,i),k=1,3)
	endif ! i.ne.1
        enddo ! 1, kolat
	if(PCFIELD.lt.1.1)stop ' PCFIELD'
	endif ! (s_opt*10.d5-int(s_opt*10.d4)*10.d0).gt.0.1d0)
	end subroutine pcfield_cont
	subroutine pc_fileld

        do n=1,kolat
	if(PCFIELD.lt.1.1) then
	elseif (PCFIELD.gt.1.1.and.PCFIELD.lt.2.1) then
	gra(:,n) = fpch(1:3,n)
        elseif (PCFIELD.gt.2.1.and.PCFIELD.lt.3.1) then
        do k=1,3
        gra(k,n) = fdip(k,n)
        enddo ! k=1,3
	endif ! PCFIELD.lt.1.1/else
        enddo ! n=1,kolat
	end subroutine pc_fileld
        end subroutine simol_inp

	subroutine  output(kolat,kolvar,iwork)
         use slspar_module, only: extended_format
	rewind iu152
	do ii=1,kolat
	if(extended_format) then
	write(iu152,1203) an(ii)
     * ,x(ii),y(ii),z(ii)
     * ,ieq(ii),impu(ii)
     * ,zmat(ii,1),zmat(ii,2),zmat(ii,3)
     * ,numx(ii,1),numx(ii,2),numx(ii,3),iepe(ii)
 1203 format((f5.2,3(f15.7),2i4,i5,2I4,i5,2I4,i5)) ! LCGTO format
	else
  203 format((f5.2,3(f15.7),2i4,i5,2I3,i5,2I3,i5)) ! LCGTO format
	write(iu152,203) an(ii)
     * ,x(ii),y(ii),z(ii)
     * ,ieq(ii),impu(ii)
     * ,zmat(ii,1),zmat(ii,2),zmat(ii,3)
     * ,numx(ii,1),numx(ii,2),numx(ii,3),iepe(ii)
	endif
	enddo
	if(extended_format) then
	write(iu152,1206) -DBLE(iwork),s_opt,eyead,pmste
     & ,ieq(kolat+1),icycle,zmat(kolat+1,1),0,0,0,0,0,iepe(kolat+1)
	else
	write(iu152,206) -DBLE(iwork),s_opt,eyead,pmste
     & ,ieq(kolat+1),icycle,zmat(kolat+1,1),0,0,0,0,0,iepe(kolat+1)   ! c:tsloc
	endif
  206 format((f5.1,3(f15.7),2i4,i5,2I3,i5,2I3,i5)) ! LCGTO format
 1206 format((f5.1,3(f15.7),2i4,i5,2I4,i5,2I4,i5)) ! LCGTO format

	if(iwork.ne.9)  close(33,status='delete') ! the run is OK so go on
        stop ' gx updated'
        end subroutine output

        subroutine  freqp(kolat,kolvar,iwork,n_00)
	iwork = iwork + 1
        read (34) (xo_eqvl(k),k=1,kolvar)
        read (34) (gr_eqvl(k),k=1,kolvar)
        read (34) (sxo(k),k=1,kolvar)
        rewind 34 
        iw = iwork - 10  

	if(iw.gt.kolvar*2) then
        call update(iwork,kolvar,xo_eqvl)
	call defjdc(kolat,kolvar,jdc,muext)
	call force(d2e,kolvar,gr_eqvl,xo_eqvl,sxo)
        call freq(kolat,kolvar,jdc,d2e,muext,GG,FF,EigVecR,n_00)
cC>>:use of freq1.f
c       call freq(kolat,kolvar,jdc,d2e,fgm,muext,EigVecR,n_00)
cC<<:use of freq1.f
	stop

	else
          do i=1,kolvar
         xo(i) = xo_eqvl(i)
        enddo
c	** define coordinates of new point and step there
        if(mod(iw,2).ne.0) then  
        xo(iw/2+1) = xo_eqvl(iw/2+1)+sxo(iw/2+1)
        else
        xo(iw/2) = xo_eqvl(iw/2)-sxo(iw/2)
        endif ! mod(iw,2).ne.0/else  
        call update(iwork,kolvar,xo)  ! -> stop
	endif !iw.gt.kolvar*2/else
        end subroutine  freqp
C:TSLOC>>
        subroutine hessp(kolat,kolvar,iwork)

	if(iwork.eq.5) then
        read(9) (xo_eqvl(i),i=1,kolvar) ! xlast
     &         ,(gr_eqvl(i),i=1,kolvar) ! glast
     &         ,(da(i)     ,i=1,kolvar) ! glast
     &         ,e_last,pg,g_line,h_line
	write(*,*)e_last,pg,g_line,h_line
	rewind 9
	call hdeta(kolvar,gg,EigVecR,xo_eqvl,gr_eqvl,xo,g)
	call locmod(1,kolvar,g,xo,EigVecR,GG,FF
     &             ,da,F_TSLOC,etot,sxo
     &             ,e_last,xo_eqvl,gr_eqvl
     &		   ,pg,g_line,h_line
     &		   ,d2e,pmste		     ) 	! -> stop here
	endif ! iwork.eq.5

	iwork = iwork + 1
        read (34) (xo_eqvl(k),k=1,kolvar)
        read (34) (gr_eqvl(k),k=1,kolvar)
        read (34) (sxo(k),k=1,kolvar)
        rewind 34 
        iw = iwork - 10  

	if(iw.gt.kolvar) then ! the dara to calc. H are ready
 	call hess(1	! one point to calc H
     &           ,d2e,kolvar,gr_eqvl,xo_eqvl,EigVecR,GG,FF
     &           ,da,F_TSLOC,etot,sxo
     &             ,pg,g_line,h_line
     &             ,d2e,pmste,gud             )  ! -> stop here

	stop

	else
	  do i=1,kolvar
	 xo(i) = xo_eqvl(i)
	enddo 
        xo(iw) = xo_eqvl(iw)+sxo(iw)
        call update(iwork,kolvar,xo)  ! ->stop
	endif !iw.gt.kolvar
        end subroutine hessp
 
	subroutine add_epe_f(kolat)
	type(ppdate)::rp
        real(r8_kind), dimension(3):: rr,rr_n1
	integer, intent(in):: kolat
	real(r8_kind)::sls_e,sls_g
	sls_e=0.0_r8_kind
	do n1=1,kolat
	rr_n1(1)=x(n1)
	rr_n1(2)=y(n1)
	rr_n1(3)=z(n1)
           
        if(an(n1)-aint(an(n1)).lt.0.00999.or.ieq(n1).eq.0) cycle        
	do n2=epe_kl,epe_nucen+kolat
        if(epe(n2)%ant-aint(epe(n2)%ant).lt.0.009999) cycle    
	mn=min(an(n1),epe(n2)%ant)
	mx=max(an(n1),epe(n2)%ant)
	call get_slsp(mn,mx,0,nr)
	
	if(nr.ne.0) then
           sls_g=0.0
           rr=rr_n1-epe(n2)%s(:)

           rss2=dot_product(rr,rr)
           rss1=sqrt(rss2)
           if(rss1*angsau.le.par%item%rqq) then
	print*, 'link ', n1, n2, rss1*angsau
              sls_e=sls_e- par%item%qq/rss1
              sls_g=par%item%qq/rss2/rss1
              
           end if
           gra(:,n1)=gra(:,n1)+rr(:)*sls_g
           rss1=rss1*angsau
           rss2=rss1**2


           sls_g=0.0
c$$$           if(rss2.le.rcuts/angsau**2) then
           if(rss2.le.par%item%cutoff**2) then
           rss4=rss2**2
           rss6=rss4*rss2
           rss8=rss4**2
           rss10=rss8*rss2
           sls_e=sls_e-evau*(par%item%c/rss6-par%item%d/rss8)
           sls_g=sls_g+gepe_const*(6.0_r8_kind*par%item%c/rss8
     &          +8.0_r8_kind*par%item%d/rss10)
              sls_e=sls_e+evau*par%item%b*exp(-rss1/par%item%r)
              sls_g=sls_g- 
     & gepe_const*par%item%b/(par%item%r*rss1)*exp(-rss1/par%item%r)
           endif! rss2.le.rcuts

           gra(:,n1)=gra(:,n1)+rr(:)*sls_g
           else
              print*,'simol_epe parameters are not found',an(n1),epe(n2)%ant
              stop
	endif                   ! nr.ne.0
	enddo ! n=2,kolat
	enddo ! n=1,kolat
	etot=etot+sls_e
	end subroutine add_epe_f

	subroutine add_sls_f(kolat)
! 	purpes: adds SLS forces for inter-cluster interactions
	type(ppdate)::rp
        real*8, dimension(3):: rr
	integer, intent(in):: kolat
	real(r8_kind)::sls_e,sls_g

	print*,'treate inter-cluster interactions'
	sls_e=0.0_r8_kind
	do n1=2,kolat
           
        if(an(n1)-aint(an(n1)).lt.0.0099.or.ieq(n1).eq.0) cycle        
	do n2=1,n1-1
        if(an(n2)-aint(an(n2)).lt.0.009999.or.ieq(n2).eq.0) cycle    
	mn=min(an(n1),an(n2))
	mx=max(an(n1),an(n2))
	call get_slsp(mn,mx,0,nr)
	if(nr.ne.0) then
	sls_g=0.0	
	rss1=r(n1,n2)*angsau
	rss2=rss1**2
	if(rss2.le.par%item%cutoff**2) then
	rss4=rss2**2
	rss6=rss4*rss2
	rss8=rss4**2
	rss10=rss8*rss2
	sls_e=sls_e-evau*(par%item%c/rss6-par%item%d/rss8)
	sls_g=sls_g+gepe_const*(6.0_r8_kind*par%item%c/rss8
     &                   +8.0_r8_kind*par%item%d/rss10)
	sls_e=sls_e+evau*par%item%b*exp(-rss1/par%item%r)
	sls_g=sls_g-
     &  gepe_const*par%item%b/(par%item%r*rss1)*exp(-rss1/par%item%r)
	endif ! rss2.le.rcuts
!!$        write(*,*) 'sls_',sls_e,sls_g,rss1/angsau,rcuts
        rr(1)=(x(n1)-x(n2))
        rr(2)=(y(n1)-y(n2))
        rr(3)=(z(n1)-z(n2))
        gra(:,n1)=gra(:,n1)+rr(:)*sls_g
        gra(:,n2)=gra(:,n2)-rr(:)*sls_g
	endif ! nr.ne.0
	enddo ! n=2,kolat
	enddo ! n=1,kolat
	etot=etot+sls_e
	end subroutine add_sls_f

	subroutine gett_balgrad(kolat)
	 type fbequ
!	  real(r8_kind),dimension(6,6)::a
!	  real(r8_kind),dimension(6)::c
!	  real(r8_kind),dimension(6)::x
	  real*8,dimension(6,6)::a,aa
	  real*8,dimension(6)::c
	  real*8,dimension(6)::x
	 end type fbequ
	type(fbequ) fb
	type(fbequ) fbs,fbr
	integer*4, intent(in)::kolat
	integer*4,dimension(6,3)::ipv
	real*8::det
	fb%a=0.0
	fb%a(:,1)=(/1.0,0.0,0.0,1.0,0.0,0.0/)
	fb%c(1)=-sum(gra(1,1:kolat))
	fb%a(:,2)=(/0.0,1.0,0.0,0.0,1.0,0.0/)
	fb%c(2)=-sum(gra(2,1:kolat))
	fb%a(:,3)=(/0.0,0.0,1.0,0.0,0.0,1.0/)
	fb%c(3)=-sum(gra(3,1:kolat))

	fb%a(2,4)=z(balance_atom_nbs(1))
	fb%a(3,4)=-y(balance_atom_nbs(1))
	fb%a(5,4)=z(balance_atom_nbs(2))
	fb%a(6,4)=-y(balance_atom_nbs(2))
	fb%c(4)=-sum(gra(2,1:kolat)*z(1:kolat))
     &          +sum(gra(3,1:kolat)*y(1:kolat))
	fb%a(1,5)=-z(balance_atom_nbs(1))
	fb%a(3,5)=x(balance_atom_nbs(1))
	fb%a(4,5)=-z(balance_atom_nbs(2))
	fb%a(6,5)=x(balance_atom_nbs(2))
	fb%c(5)=sum(gra(1,1:kolat)*z(1:kolat))
     &         -sum(gra(3,1:kolat)*x(1:kolat))
	fb%a(1,6)=y(balance_atom_nbs(1))
	fb%a(2,6)=-x(balance_atom_nbs(1))
	fb%a(4,6)=y(balance_atom_nbs(2))
	fb%a(5,6)=-x(balance_atom_nbs(2))
	fb%c(6)=-sum(gra(1,1:kolat)*y(1:kolat))
     &          +sum(gra(2,1:kolat)*x(1:kolat))
	print*,'fb%c'
	print*, fb%c
	do i=1,6
	fb%a(i,i)=fb%a(i,i)+1.d-8
	enddo
	fbs%a=fb%a
        call matinv(fb%a,6,det,6,ipv)
 	print*,'get_balgrad:',det
	fbr%a=matmul(fbs%a,fb%a)
	print*

	fb%x=matmul(fb%c,fb%a)
	print*,fb%x
	fbr%aa(:,1)=matmul(fb%x,fbs%a)
	gra(:,balance_atom_nbs(1))=fb%x(1:3)
	gra(:,balance_atom_nbs(2))=fb%x(4:6)
	print*,fbr%aa(:,1)
	
	end subroutine gett_balgrad

	subroutine get_balgrad(kolat)
	integer*4, intent(in)::kolat
	integer*4,dimension(2,3)::ipv
	 type fbequ
!!!	  real(r8_kind),dimension(6,6)::a
!!!	  real(r8_kind),dimension(6)::c
!!!	  real(r8_kind),dimension(6)::x
	  real*8,dimension(2,2)::a
	  real*8,dimension(2)::c
	  real*8,dimension(2)::x
	 end type fbequ
	type(fbequ) fb
	type(fbequ) fbs,fbr
	real*8::det
	fb%a=0.0
	fb%a(:,2)=(/0.0,1.0/)
	fb%c(2)=-sum(gra(3,1:kolat))

	fb%a(1,1)=-x(balance_atom_nbs(1))
	fb%c(1)=-sum(gra(1,1:kolat)*y(1:kolat))
     &          +sum(gra(2,1:kolat)*x(1:kolat))
	fbs%a=fb%a
	      call matinv(fb%a,2,det,2,ipv)
	fbr%a=matmul(fb%a,fbs%a)
	fb%x=matmul(fb%c,fb%a)
	fb%x=matmul(fb%a,fb%c)
	print*, fb%x,' fb%x'
	gra(2,balance_atom_nbs(1))=fb%x(1)
	gra(2,balance_atom_nbs(2))=-fb%x(1)
	gra(3,balance_atom_nbs(2))=fb%x(2)
	end subroutine get_balgrad


        subroutine sf_zero(gra,kolat,n_99)
c       ** make sum of forces equal to  zero
        real*8 gra(3,kolat)

        do n=1,kolat
        if(n.ne.n_99) then
        do k=1,3
         gra(k,n_99) =gra(k,n_99)-gra(k,n)
        enddo
        endif ! n.ne.n_99
        enddo
        return
        e n d subroutine sf_zero

      e n d module simol_io
