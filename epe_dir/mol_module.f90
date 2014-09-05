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
module mol_module
#include <def.h>
  use type_module
  use epecom_module,only: output_epe, epeit_count,N_VACANCIES,N_IMPURITIES

  implicit none
  save
  private

  integer(kind=i4_kind), public :: N_var,L_var
  integer(kind=i4_kind):: alloc_status
  real(kind=r8_kind), public :: dielectric_const

  public :: madelung, trans_vectors,  &
       mott_litleton, epe_generater, sorting_epeinout

!---------------------------------------------------------
contains


  subroutine madelung
    ! **calc madelung potentials and energies

    use epecom_module
    use culon_module, only: erfc

    real(kind=r8_kind) ::distance,dir_space,lim_contrib,back_space,gix2pi
    integer(kind=i4_kind) :: i_cell,j_epe,nmt,mad_index,n_epe_f_i_cell,k

    if(do_print) &
         write(output_epe, &
         "( ' No ion charge type first neighb. dist. energy(eV)',  &
         & ' potential(eV/a) madelung const.  potential(a.u) ')  ")

    do i_cell=1,n_ions_cell
       do j_epe=reg_2a_n_ions,1,-1
          if(epe(j_epe)%k.eq.type_of_ion(i_cell))nmt=j_epe
       enddo !
       mad_index=epe(nmt)%k
!!$ print*,'mad_index defined',i_cell,mad_index
! **nmt is the smallest number of ion with type type_of_ion(i)
       dir_space=zero
       madc%rmd(i_cell)=1.e20_r8_kind
       n_epe_f_i_cell=which_epe_ion(i_cell)%new

! **calculation of potential - summarizing over direct lattice
       do j_epe=reg_2a_n_ions,1,-1
          if(j_epe.eq.n_epe_f_i_cell) cycle
          distance=sqrt(dot_product(epe(n_epe_f_i_cell)%r-epe(j_epe)%r, &
               epe(n_epe_f_i_cell)%r-epe(j_epe)%r))
          if(distance.lt.0.0001) then
!!$             print*,'epe ions coincide' ,j_epe,n_epe_f_i_cell
!!$             print*,epe(j_epe)%r,epe(n_epe_f_i_cell)%r
          end if

          if(distance.lt.madc%rmd(i_cell))madc%rmd(i_cell)=distance
          if(distance.gt.epe(reg_2a_n_ions)%d) cycle
          dir_space=dir_space+q_zl(j_epe)*(erfc(error_function_parameter*distance))/distance
!       if(distance.le.3.75397976337778) print*,distance, &
!                             q_zl(j_epe)*(erfc(error_function_parameter*distance))/distance, &
!                             error_function_parameter,(erfc(error_function_parameter*distance))/distance
!
       enddo

       dir_space=dir_space*qau_qepe

       lim_contrib= erfo*q_zl(n_epe_f_i_cell)*qau_qepe

    ! **sum over reciprocial lattice
       back_space=0.d0
       do k=1,n_bs_points
          gix2pi=dot_product(gstr(k,:),epe(n_epe_f_i_cell)%r)*pi2
          back_space=back_space+rcos(k)*cos(gix2pi)+rsin(k)*sin(gix2pi)
       enddo ! k=1,n_bs_points

!      DPRINT 'back_space dir_space lim for i_cell ion',i_cell,back_space,dir_space,lim_contrib,qau_qepe
!      DPRINT 'q_zl(1)',q_zl(1)

!!$    write(output_epe,*)      back_space*qau_qepe/eau_ev  &
!!$         ,dir_space*qau_qepe/eau_ev&
!!$         ,lim_contrib*qau_qepe/eau_ev, &
!!$       error_function_parameter*auangs

       madc%fmd(i_cell)= (dir_space+lim_contrib+back_space)
       madc%emd(i_cell)=madc%fmd(i_cell)*(qau_qepe*q_zl(n_epe_f_i_cell))
       madc%amd(i_cell)=abs(madc%fmd(i_cell)* &
                       madc%rmd(i_cell)/(qau_qepe*q_zl(n_epe_f_i_cell)))

       if(nmt.eq.n_epe_f_i_cell) then
          mad%rmd(mad_index)=madc%rmd(i_cell)
          mad%fmd(mad_index)= (dir_space+lim_contrib+back_space)
          mad%emd(mad_index)=mad%fmd(mad_index)*(qau_qepe*q_zl(n_epe_f_i_cell))
          mad%amd(mad_index)=abs(mad%fmd(mad_index)*mad%rmd(mad_index) &
                                       /(qau_qepe*q_zl(n_epe_f_i_cell)))

       endif
       if(do_print) &
            write(output_epe, 100)i_cell,n_epe_f_i_cell,epe(n_epe_f_i_cell)%q, &
            !         q_zl(n_epe_f_i_cell), &
                             mad_index,madc%rmd(i_cell),madc%emd(i_cell)  &
                                      ,madc%fmd(i_cell),madc%amd(i_cell)  &
                                      ,madc%fmd(i_cell)*qau_qepe/eau_ev
100    format(1x,i5,i4,f9.4,i6,4x,f9.4,7x,f11.5,1x,f11.5,5x,f11.5,2x,f11.5)
    enddo

  end subroutine madelung
!-------------------------------------------------------------
!-------------------------------------------------------------

  subroutine trans_vectors
    ! **finding recirocal translation vectors and cell volume

    use epecom_module
    use culon_module, only: MINV

    real(kind=r8_kind) :: tol,sqaa1,sqaa2,sqaa3
    real(kind=r8_kind) :: gig,resm,qns,qnc,gcel
    real(kind=r8_kind) :: psit,res,psi,prsin,prcos,gix2pi
    integer(kind=i4_kind) :: L(3),M(3)
    integer(kind=i4_kind) :: maxb1,maxb2,maxb3,ii,jj,kk,i,j,mm,ng,io

    TOL=1.E-9
    sqaa1=sqrt(dot_product(VECTORS_TRANS(1,:),VECTORS_TRANS(1,:)))
    sqaa2=sqrt(dot_product(VECTORS_TRANS(2,:),VECTORS_TRANS(2,:)))
    sqaa3=sqrt(dot_product(VECTORS_TRANS(3,:),VECTORS_TRANS(3,:)))
    gmat=VECTORS_TRANS
    if(epeit_count == 0 .and. do_print) &
         write(output_epe, 50) GMAT(:, 1)
    if(epeit_count == 0 .and. do_print) &
         write(output_epe, 50) GMAT(:, 2)
    if(epeit_count == 0 .and. do_print) &
         write(output_epe, 50) GMAT(:, 3)
    call MINV(GMAT,3,VC,L,M)
    VC=abs(VC)
    MAXB1=int(20/VC**0.333333*SQAA1)*2+1
    MAXB2=int(20/VC**0.333333*SQAA2)*2+1
    MAXB3=int(20/VC**0.333333*SQAA3)*2+1
    if(epeit_count == 0 .and. do_print) &
         write(output_epe,*) 'No of back space translation ',MAXB1,MAXB2,MAXB3
    if(epeit_count == 0 .and. do_print) &
         write(output_epe, 30)VC
30  format(5x,'unit cell volume   ',e15.4,' A**3')
    if(epeit_count == 0 .and. do_print) &
         write(output_epe, 40)
40  format(5x,'translation vectors of reciprocial cell  ')
    if(epeit_count == 0 .and. do_print) then
       do J=1,3
          write(output_epe, 50)(GMAT(I,J),I=1,3)
50        format(5x,3E12.5)
       enddo !J=1,3
    end if
    QPIVC=Qau_qepe/(PI*VC)
    PTA2=(PI/ERROR_FUNCTION_PARAMETER)**2
    if(epeit_count == 0 .and. do_print) &
         write(output_epe,*) 'PTA2 QPIVC ', PTA2,QPIVC

! ** RESM - maximum value in sum over reciprocial lattice
    RESM=exp(-PTA2*dot_product(GMAT(1,:),GMAT(1,:)))/dot_product(GMAT(1,:),GMAT(1,:))

! **full set of reciprocial vectors having contribution in
! **function of gauss distribution of charge to GSTR.
! **calculated multiplicator of each vector to RSIN,RCOS

    n_bs_points=1
    do II=1,maxb1
       L(1)=II-maxb1/2+1
       do JJ=1,maxb2
          L(2)=JJ-maxb2/2+1
          do KK=1,maxb3
             L(3)=KK-maxb3/2+1
             GIG=0.
             do I=1,3
                GSTR(n_bs_points,I)=0.
                do J=1,3
                   GSTR(n_bs_points,I)=GSTR(n_bs_points,I)+L(J)*GMAT(j,i)
                enddo   ! J=1,3
                GIG=GIG+GSTR(n_bs_points,I)*GSTR(n_bs_points,I)
             enddo      ! I=1,3
             if(GIG.lt.1.E-10) cycle
             QNS=zero
             QNC=zero
             do MM=1,N_IONS_CELL
                GCEL=dot_product(GSTR(n_bs_points,:),R_ION_IN_CELL(MM,:))
                QNS=QNS+Q_Z(MM)*sin(PI2*GCEL)
                QNC=QNC+Q_Z(MM)*cos(PI2*GCEL)
             enddo ! MM=1,N_IONS_CELL

! **QNC, QNS - parts of structural factor
             if(abs(QNS).lt.1.d-9.and.abs(QNC).lt.1.d-9) cycle

! **GIG - square of recirocal vector does not equal zero

             RES=exp(-PTA2*GIG)/GIG
             if(RES/RESM.lt.TOL) cycle
! **RES/RESM  - precision of calculation of recirocal sum

             RSIN(n_bs_points)=RES*QNS*QPIVC
             RCOS(n_bs_points)=RES*QNC*QPIVC
             n_bs_points=n_bs_points+1
             if(n_bs_points.gt.ndngv) then
                print*, 'RES/RESM', RES,RESM
                write(*, "('  insufficient precision in STR on calculation of inver',  &
                     &  'sion vectors tol,res,n_bs_points =',/2e15.5,i10)" &
                     )tol,RES/RESM,n_bs_points
!!$10        format(1X,'insufficient precision in STR on calculation of inver',  &
!!$                    'sion vectors tol,res,n_bs_points =',/2e15.5,i10)
                ASSERT(n_bs_points.le.ndngv)
             endif ! n_bs_points.gt.ndngv
          enddo !KK=1,maxb3
       enddo    ! JJ=1,9
    enddo ! II=1,9

    n_bs_points=n_bs_points-1

    if(periodic_optimization)then
       do i=1,n_ions_cell
          gpsi(:,:n_ions_cell)=zero
          PSI=zero
          do ng=1,n_bs_points
             psit=0.
             do io=1,n_ions_cell
                GCEL=dot_product(GSTR(NG,:),R_ION_IN_CELL(io,:))
                gig=dot_product(GSTR(NG,:),GSTR(NG,:))
                QNS=Q_Z(io)*sin(PI2*GCEL)
                QNC=Q_Z(io)*cos(PI2*GCEL)
                RES=exp(-PTA2*GIG)/GIG  ! scalar
                PRSIN=RES*QNS*QPIVC
                PRCOS=RES*QNC*QPIVC
                GIX2PI=dot_product(GSTR(ng,:),R_ION_IN_CELL(i,:))*PI2
                PSI=PSI+pRSIN*sin(GIX2PI)+pRCOS*cos(GIX2PI)  ! sin & cos are functions of i only
                gPSI(:,i)=gPSI(:,i)+(pRSIN*cos(GIX2PI)-pRCOS*sin(GIX2PI))*GSTR(NG,:)*PI2
                gPSI(:,io)=gPSI(:,io)+GSTR(NG,:)*PRcos*PI2*sin(GIX2PI) &
                     -GSTR(NG,:)*PRsin*PI2*cos(GIX2PI)
             enddo ! 1,n_ions_cell
          enddo ! ng=1,n_bs_points
       enddo ! i=1,n_ions_cell
    endif       !periodic_optimization
    if(epeit_count == 0 .and. do_print) &
         write(output_epe, *) ' number of reciprocial lattice vectors', n_bs_points
  end subroutine trans_vectors
!---------------------------------------------------------------
!---------------------------------------------------------------

  subroutine mott_litleton
    ! **procedure of calculation of mott-litleton parameters for
    ! **finding nuclear and shall displacements in area 2a and
    ! **energy in area 2b

    use epecom_module,deltam=>delta
    use iounitadmin_module
    use culon_module, only: common_atom

    real(kind=r8_kind),dimension(ndt) :: DMOM
    real(kind=r8_kind) :: R2_first_coord_sphere,r_first_coord_sphere, &
         r2_curent,delta1,ml_dprime,sdmom,salf,reversed_fo
    real(kind=r8_kind) :: fo,alfo,alfh,pb,pbo,pbh,epsh,zeff
    integer(kind=i4_kind) :: ii_counter, i_counter, &
         j,i_ion_type,k_ion_type,k,kk,first_i_epe_ion,coord_number
    real(kind=r8_kind) :: ml_polarizability(max_type_ions) ! automatic arr
    real(kind=r8_kind), parameter::  one=1.0_r8_kind, three=3.0_r8_kind
    integer (kind=i4_kind):: ml_ten_unit,control_sum,ii,jj,jbuf

    integer (kind=i4_kind):: numat,ic
    namelist/ml_tensors_nml/ numat
    c1: do i_counter=1,N_IONS_CELL

       do j=1,reg_2a_n_ions
          first_i_epe_ion=j
          if(epe(J)%k.eq.TYPE_OF_ION(i_counter)) exit
       enddo
       i_ion_type=epe(first_i_epe_ion)%k
       ii=which_epe_ion(i_counter)%new

       do ii_counter=1,i_counter-1
          if(TYPE_OF_ION(ii_counter).eq.TYPE_OF_ION(i_counter)) cycle c1
       enddo

! **determining type of nearest neighbours having other type
! **than i_ion_type
       ml_dprime_factors(i_ion_type) = 0
       c2: do K=1,N_IONS_CELL

          k_ion_type=TYPE_OF_ION(K)
          if(TYPE_OF_ION(i_counter).eq.TYPE_OF_ION(K)) cycle
          if(Q_Z(i_counter).eq.Q_Z(K)) cycle
          do KK=1,K-1
             if(TYPE_OF_ION(kk).eq.TYPE_OF_ION(k)) cycle c2
!       check if ion of this type is already treated
          enddo
          jj=which_epe_ion(k)%new

!       define radius of first coordination sphre and its coordination number
          coord_number=0
          R2_first_coord_sphere=20.*scaling_factor
          do J=reg_2a_n_ions,1,-1
             if(i_ion_type.eq.epe(J)%k)cycle !no ions of the same type in first coord sphere
             if(TYPE_OF_ION(K).ne.epe(J)%k)cycle
             r2_curent=dot_product(epe(first_i_epe_ion)%r-epe(J)%r,epe(first_i_epe_ion)%r-epe(J)%r)
             jbuf=j
             delta1=r2_curent-R2_first_coord_sphere
             if(abs(delta1).gt.0.05_r8_kind) then
                if(delta1.lt.0.) then
                   R2_first_coord_sphere=r2_curent
                   coord_number=0
                else
                   cycle
                endif
             endif
             coord_number=coord_number+1
          enddo
          r_first_coord_sphere=sqrt(R2_first_coord_sphere)

! **second derivative of short-range potential
          if(R2_first_coord_sphere.GT.16.0_r8_kind) then
             ic=1
          else
             ic=common_atom(first_i_epe_ion,jbuf)+1                     !!!!!!!!!!!!!AS
             if(host%ro(i_ion_type,k_ion_type,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
          endif
          if(do_print) then
             if(abs(host%RO(i_ion_type,k_ion_type,ic)).lt.1.d-5.or.r_first_coord_sphere.lt.1.d-1) &
                  print*,'mott_litleton:', i_ion_type,k_ion_type,host%RO(i_ion_type,k_ion_type,ic)
          end if

          ml_dprime=coord_number/3.0_r8_kind* &
               (host%B(i_ion_type,k_ion_type,ic)/host%RO(i_ion_type,k_ion_type,ic)*  &
               (1.0_r8_kind/host%RO(i_ion_type,k_ion_type,ic)-2.0_r8_kind/r_first_coord_sphere)* &
               exp(-r_first_coord_sphere/host%RO(i_ion_type,k_ion_type,ic))- &
               30.0_r8_kind*host%C(i_ion_type,k_ion_type,ic)/R2_first_coord_sphere**4 &
               -56.0_r8_kind*host%D(i_ion_type,k_ion_type,ic)/R2_first_coord_sphere**5)
          write(output_epe, "(1X,'ion of type',I3,'  has  ',I3 &
               &,' neighbours of type',I3,  &
               &  '  ml_dprime=',F16.7)" &
               )i_ion_type,coord_number,k_ion_type,ml_dprime
          ml_dprime_factors(i_ion_type)=ml_dprime_factors(i_ion_type)+ml_dprime
!3     continue
       enddo  c2
    enddo  c1
    if(do_print) write(output_epe,*) 'ml_dprime_factors are defined'

! **calculation of polarization and dielectric permeability
    SDMOM=zero
    SALF=zero
    reversed_fo=zero
    do i_counter=1,N_IONS_CELL
       K=TYPE_OF_ION(i_counter)
       if( PK(TYPE_OF_ION(i_counter)).lt.0.00001_r8_kind) then
          print*, 'wrong pk for ion type', TYPE_OF_ION(i_counter)
          call error_handler('PK < 0.00001')
       endif
       DMOM(TYPE_OF_ION(i_counter))=Q_SHELL(TYPE_OF_ION(i_counter))/PK(TYPE_OF_ION(i_counter))
       ml_polarizability(TYPE_OF_ION(i_counter))=DMOM(TYPE_OF_ION(i_counter))*Q_SHELL(TYPE_OF_ION(i_counter))
       SDMOM=SDMOM+Q_Z(i_counter)*DMOM(TYPE_OF_ION(i_counter))
       SALF=SALF+ml_polarizability(TYPE_OF_ION(i_counter))
       reversed_fo=reversed_fo+Q_Z(i_counter)**2*(1./PK(K)+0.5_r8_kind/ml_dprime_factors(TYPE_OF_ION(i_counter)))
    enddo       ! I=1,N_IONS_CELL

    if(reversed_fo.lt.0.00000001_r8_kind) call error_handler('wrong reversed_fo')

    FO=one/reversed_fo
    ALFO=SALF+Qau_qepe**2/FO-two*qau_qepe*SDMOM
    ALFH=SALF-FO*SDMOM**2
    ZEFF=(qau_qepe-FO*SDMOM)/qau_qepe
    PB=4.*PI/(3.*VC)
    PBO=PB*ALFO
    PBH=PB*ALFH
    EPSH=(one+two*PBH)/(one-PBH)
    dielectric_const=(one+two*PBO)/(one-PBO)
    if(do_print) then
       write(output_epe,  &
            "(5X,'   static polarizability              ',F8.3,/,  &
            & 5X,'   static dielectric constant         ',F8.3,/,  &
            & 5X,'   high frequency polarizability      ',F8.3,/,  &
            & 5X,'   high frequency dielectric constant ',F8.3,/,  &
            & 5X,'   effective charge of Szigety        ',F8.3,/,  &
            & 5X,'   paramiter FO                       ',F8.3) ") &
            ALFO,dielectric_const,ALFH,EPSH,ZEFF,FO
    end if

    if(fixed_dielectric_const) then
       dielectric_const=fx_dielectric_const
       if(do_print) then
          print *,' dielectric conttant is fixed at value ', dielectric_const
          write(output_epe,*) ' dielectric conttant is fixed at value ', dielectric_const
       end if
    endif
! **calculation of mott-litleton parameters nuclears and shells
    do i_counter=1,N_IONS_CELL
       K=TYPE_OF_ION(i_counter)
       ml_fac(K)%s=(qau_qepe*Q_Z(i_counter))/(2.0*ml_dprime_factors(K))
       ml_fac(K)%c= (Q_NUCLEAR(K)/PK(K)+ml_fac(K)%s)*(dielectric_const+two)/(three*dielectric_const)
       ml_fac(K)%s=ml_fac(K)%s*(dielectric_const+two)/(three*dielectric_const)
       if(do_print) write(output_epe, 30)K,ml_fac(K)%c,ml_fac(K)%s,ml_polarizability(K)
30     format(1X,'ion type',I3,' mott-litleton paramiter for nuclear ',  &
            F8.4,' for shell ',F8.4,/1x,' polarizability',F8.3)
    enddo
        if(ml_tensors) then
        if(.not.allocated(ml_ten)) allocate (ml_ten(N_IONS_CELL))
        ml_ten_unit=openget_iounit(trim(epe_input_dir)//'/ml_tensors_in', &
        status='old')
        read (ml_ten_unit,nml=ml_tensors_nml)
        print*,'tensors from file ml_tensors_in'
        do i_counter=1,N_IONS_CELL
           read (ml_ten_unit,*)
           read (ml_ten_unit,*) ml_ten(i_counter)%c(:,1)
           read (ml_ten_unit,*) ml_ten(i_counter)%c(:,2)
           read (ml_ten_unit,*) ml_ten(i_counter)%c(:,3)
           print*,i_counter
           print*, ml_ten(i_counter)%c(:,1)
           print*, ml_ten(i_counter)%c(:,2)
           print*, ml_ten(i_counter)%c(:,3)
        end do
        control_sum=N_IONS_CELL
        do i_counter=1,N_IONS_CELL
           if(do_print) print*,pk(TYPE_OF_ION(i_counter)) ,i_counter
           if(pk(TYPE_OF_ION(i_counter)).lt.9999.0) then
            control_sum=  control_sum+1
            read (ml_ten_unit,*)
           read (ml_ten_unit,*) ml_ten(i_counter)%s(:,1)
           read (ml_ten_unit,*) ml_ten(i_counter)%s(:,2)
           read (ml_ten_unit,*) ml_ten(i_counter)%s(:,3)
           print*, control_sum
           print*, ml_ten(i_counter)%s(:,1)
           print*, ml_ten(i_counter)%s(:,2)
           print*, ml_ten(i_counter)%s(:,3)

           else
             ml_ten(i_counter)%s=ml_ten(i_counter)%c
           end if

        end do
        if(control_sum.ne.numat) then
           call error_handler('wrong structure of ml_tensors')
        end if
        call returnclose_iounit(ml_ten_unit)

        endif
  end subroutine mott_litleton
!------------------------------------------------------------
!------------------------------------------------------------
  subroutine sorting_epeinout(n_at_IIa)
    ! This subroutine resorts core and shell positions taken
    ! from epeinout file. They have to correspond to the order
    ! of regular atomic positions generated by EPE_GENERATER
    ! subroutine. Unfortunately on machines of different
    ! architecture the order of the epe array can be different.
    ! Such effect can prevent from transfering EPE tasks from one
    ! machine to another.
     use epecom_module

     integer(i4_kind) :: n_at_IIa
     integer(i4_kind) :: i,j,n_start,status,num,num1,n_checked
     real(r8_kind) :: shell_buf(3), core_buf(3),dist2,min_dist2
     logical, allocatable :: checked(:)

     allocate(checked(n_at_IIa),stat=status)
     ASSERT(status==0)
     checked=.false.

     n_checked=0
     do
        num1=0
        do j=1,n_at_IIa
           if(checked(j)) cycle
              num1=j
              exit
        enddo

        if(num1 == 0) exit

        min_dist2=1.0e1_r8_kind
        n_start=num1
        do i=n_start,n_at_IIa
           dist2=dot_product(epe(i)%o-epe(num1)%r,epe(i)%o-epe(num1)%r)
           if(dist2 .ge. min_dist2) cycle
              min_dist2=dist2
              num=i
        end do
        checked(num1)=.true.
        if(num.eq.num1) cycle
        n_checked=n_checked+1

        shell_buf=r_sh_ion(num1,:); core_buf=r_nuc_ion(num1,:)
        r_sh_ion(num1,:)=r_sh_ion(num,:); r_nuc_ion(num1,:)=r_nuc_ion(num,:)
        r_sh_ion(num,:)=shell_buf; r_nuc_ion(num,:)=core_buf
        core_buf=epe(num1)%o
        epe(num1)%o=epe(num)%o
        epe(num)%o=core_buf

!        print*, num1,num,min_dist2,n_checked
!        print*, epe(num1)%r
!        print*, r_sh_ion(num1,:)

     end do
     print*,'n_checked=',n_checked

     deallocate(checked,stat=status)
     ASSERT(status==0)

  end subroutine sorting_epeinout
!------------------------------------------------------------
!------------------------------------------------------------

  subroutine epe_generater(nions1,do_3b)
    ! **lattice generation

    use epecom_module

    logical :: do_3b
    integer(kind=i4_kind), intent(in) :: nions1
    integer(kind=i4_kind) :: limx,limy,limz
    integer(kind=i4_kind) :: lx,ly,lz,each_ion_counter,ii,jj,kk,j,i,io,ic,k,i1,j1,ll
    integer(kind=i4_kind) :: gen_center,nk,l0,l3,nmm
    integer(kind=i4_kind) :: n3,status
    real(kind=r8_kind) :: translation(3),u(3),min_product,product,counter(3), &
         position(3), position2(3),s(3),cr(3)
    real(kind=r8_kind), parameter:: one=1.0_r8_kind,three=3.0_r8_kind,coord_sph_eps=0.2_r8_kind
    integer(kind=i4_kind) :: dim
    real(kind=r8_kind) :: rr_max
    logical, allocatable:: check_ind(:,:,:)
    real(r8_kind), allocatable :: r_buf(:,:)

    type(ep_enviroment) :: epe_buf

! **coordinate centre
    limx=n_trans_primitive_cell(1)*2
    limy=n_trans_primitive_cell(2)*2
    limz=n_trans_primitive_cell(3)*2

    allocate(check_ind(limx,limy,limz),stat=status)
    ASSERT(status==0)
    check_ind=.false.

    dim=0
    do ll=1,max(limx,limy,limz)
       lx=ll; if(lx > limx) lx=limx
       ly=ll; if(ly > limy) ly=limy
       lz=ll; if(lz > limz) lz=limz
       do ii=1,lx
          if(mod(ii,2) /= zero) then
             counter(1)=real(ii-1,r8_kind)/2.0_r8_kind
          else
             counter(1)=-real(ii,r8_kind)/2.0_r8_kind
          end if
          do jj=1,ly
             if(mod(jj,2) /= zero) then
                counter(2)=real(jj-1,r8_kind)/2.0_r8_kind
             else
                counter(2)=-real(jj,r8_kind)/2.0_r8_kind
             end if
             m: do kk=1,lz
                if(mod(kk,2) /= zero) then
                   counter(3)=real(kk-1,r8_kind)/2.0_r8_kind
                else
                   counter(3)=-real(kk,r8_kind)/2.0_r8_kind
                end if

                if(check_ind(ii,jj,kk)) cycle m
                check_ind(ii,jj,kk) = .true.

                do j=1,3
                   translation(j)=dot_product(counter(:),vectors_trans(j,:))
                enddo ! j=1,3
                do io=1,n_ions_cell     !loop over all ions of the cell
                   u(:)=r_ion_in_cell(io,:)+translation(:)
                   if(core_shell) then
                      cr(:)=r_core_in_cell(io,:)+translation(:)
                      s(:)=r_shell_in_cell(io,:)+translation(:)
                   endif

                   min_product=1.e20_r8_kind
                   do ic=1,n_centres_of_generation
                      product=sqrt(dot_product(r_cent_gener(ic,:)-u,r_cent_gener(ic,:)-u))
                      if(product < min_product) then
                         min_product=product
                         gen_center=ic
                      end if
                   end do

                   if(dim == 0) then
                      epe(1)%r=u
                      epe(1)%c=cr
                      epe(1)%s=s
                      epe(1)%d=min_product
                      epe(1)%q=q_z(io)
                      epe(1)%k=type_of_ion(io)
                      epe(1)%m=io
                      epe(1)%fix=fixed(io)
                      epe(1)%gc=gen_center
                      rr_max=min_product
                      goto 1
                   endif

                   if(min_product < rr_max) then
                      do i1=1,dim
                         if (min_product < epe(i1)%d) then
                            if(i1 == nions1) then
                               epe(i1)%r=u
                               epe(i1)%c=cr
                               epe(i1)%s=s
                               epe(i1)%d=min_product
                               epe(i1)%q=q_z(io)
                               epe(i1)%k=type_of_ion(io)
                               epe(i1)%m=io
                               epe(i1)%fix=fixed(io)
                               epe(i1)%gc=gen_center
                               rr_max=min_product
                               exit
                            end if
                            if(dim < nions1) then
                               do j1=dim,i1,-1
                                  epe(j1+1)%r=epe(j1)%r
                                  epe(j1+1)%c=epe(j1)%c
                                  epe(j1+1)%s=epe(j1)%s
                                  epe(j1+1)%d=epe(j1)%d
                                  epe(j1+1)%q=epe(j1)%q
                                  epe(j1+1)%k=epe(j1)%k
                                  epe(j1+1)%m=epe(j1)%m
                                  epe(j1+1)%fix=epe(j1)%fix
                                  epe(j1+1)%gc=epe(j1)%gc
                               end do
                               rr_max=epe(dim+1)%d
                            else
                               do j1=dim-1,i1,-1
                                  epe(j1+1)%r=epe(j1)%r
                                  epe(j1+1)%c=epe(j1)%c
                                  epe(j1+1)%s=epe(j1)%s
                                  epe(j1+1)%d=epe(j1)%d
                                  epe(j1+1)%q=epe(j1)%q
                                  epe(j1+1)%k=epe(j1)%k
                                  epe(j1+1)%m=epe(j1)%m
                                  epe(j1+1)%fix=epe(j1)%fix
                                  epe(j1+1)%gc=epe(j1)%gc
                               end do
                               rr_max=epe(dim)%d
                            end if
                            epe(i1)%r=u
                            epe(i1)%c=cr
                            epe(i1)%s=s
                            epe(i1)%d=min_product
                            epe(i1)%q=q_z(io)
                            epe(i1)%k=type_of_ion(io)
                            epe(i1)%m=io
                            epe(i1)%fix=fixed(io)
                            epe(i1)%gc=gen_center
                            exit
                         end if
                      end do
                   else

                      if(dim < nions1)then
                         epe(dim+1)%r=u
                         epe(dim+1)%c=cr
                         epe(dim+1)%s=s
                         epe(dim+1)%d=min_product
                         epe(dim+1)%q=q_z(io)
                         epe(dim+1)%k=type_of_ion(io)
                         epe(dim+1)%m=io
                         epe(dim+1)%fix=fixed(io)
                         epe(dim+1)%gc=gen_center
                         rr_max=min_product
                      endif
                   endif

1                  dim=dim+1
                   if(dim > nions1) dim=nions1
                end do
             end do m
          end do
       end do
    end do

    deallocate(check_ind,stat=status)
    ASSERT(status==0)

    each_ion_counter=dim
    if(each_ion_counter == nions1) each_ion_counter=nions1-1

#ifdef NEW_EPE
    if(.not.allocated(which_epe_ion)) then
       allocate(which_epe_ion(n_ions_cell),stat=status)
       ASSERT(status == 0)
    endif
#endif
 !  print*,size(which_epe_ion),n_ions_cell,'size n_ions_cell'
    which_epe_ion(1:n_ions_cell)%new=nions1

    if(option_c3_symm) then
       if(.not.allocated(c3)) then
          allocate (c3(3),stat=status)
          ASSERT(status==0)
       endif
       call rot(c3(1)%rotmat,zero,zero,one,zero)
       call rot(c3(2)%rotmat,zero,zero,one,two*pi/three)
       call rot(c3(3)%rotmat,zero,zero,one,-two*pi/three)
!AS : I do not undestand a necessity of this part
!     Everything works correct without this part, but
!     with this part uncommented errors can take place
!!$       kl=1
!!$       do while(kl.le.each_ion_counter)
!!$          position=epe(kl)%r
!!$          nk=kl
!!$          do j=1,3
!!$             position2=matmul(c3(j)%rotmat,position)
!!$             if(j.eq.2.and.dot_product(position2-position,position2-position).lt.0.5_r8_kind) exit
!!$             do k=kl,each_ion_counter
!!$                if(dot_product(position2-epe(k)%r,position2-epe(k)%r).lt.0.5_r8_kind) then
!!$                   varepe=epe(kl)
!!$                   epe(kl)=epe(k)
!!$                   epe(k)=varepe
!!$                epe(kl)%sym_st_ind=nk
!!$                exit
!!$                endif
!!$             enddo ! each_ion_counter
!!$             kl=kl+1
!!$          enddo ! j=1,3
!!$       enddo ! while
    endif

! **finding number of spheres with ions same type (epe%k) and
! **number of last ions on this spheres
    l0=1
    do nk=1,each_ion_counter
       if(epe(nk)%k.eq.epe(nk+1)%k) then
          if(abs(epe(nk)%d-epe(nk+1)%d).le.0.01_r8_kind) cycle
       endif
       epe(l0)%u=nk
       l0=l0+1
    enddo
    l0=l0-1

! **define epe%l epe%u and print results of generation
    l_var=l0
    n_var=epe(l_var)%u
    if(epeit_count == 0 .and. do_print) &
         write(output_epe, 102)
102 format(11x,'x',11x,'y',11x,'z',  &
         11x,'epe%r',2x,'    q_zl    l_var  n center cell')
    do l3=1,l0
       if(l3.ne.1) then
          epe(l3)%l=epe(l3)%u-epe(l3-1)%u
          nmm=epe(l3-1)%u+1
       else
          epe(1)%l=epe(1)%u
          nmm=1
       endif ! l3.ne.1/else
       do n3=nmm,epe(l3)%u
          if(n3.lt.which_epe_ion(epe(n3)%m)%new) which_epe_ion(epe(n3)%m)%new=n3
       enddo
       if(epeit_count == 0 .and. do_print) then
          do n3=nmm,epe(l3)%u
             write(output_epe, '(5x,5f12.5,4i5,1x,l1)') epe(n3)%r/scaling_factor,  &
                  epe(n3)%d/scaling_factor ,epe(n3)%q,l3,n3,epe(n3)%gc,epe(n3)%m,epe(n3)%fix
          enddo
       end if
    enddo ! l3=1,l0

    q_zl(1:each_ion_counter)=epe(1:each_ion_counter)%q !!!
    nmvmax=0
    do i=1,n_ions_cell
       if(epeit_count == 0 .and. do_print) &
            write(output_epe,*) 'epe_generater:  which_epe_ion(i) ',i, which_epe_ion(i)%new
       nmvmax=max(nmvmax,which_epe_ion(i)%new)
       if(which_epe_ion(i)%new.gt.ndr1.and.periodic_optimization) then
          write(output_epe,*) ' epe_generater: size of region 1 should be incrised for '
          write(output_epe,*) ' periodic optimization'
          write(output_epe,*) which_epe_ion(i)%new,ndr1
          stop
       endif
    enddo       ! i=1,n_ions_cell
    if(epeit_count == 0 .and. do_print) &
         write(output_epe,*) 'epe_generater:  nmvmax ',nmvmax

    do  l3=1,l_var
       if(epe(epe(l3)%u)%d.lt.radius_first_sphere) reg_i_n_ions=epe(l3)%u
       if(epe(epe(l3)%u)%d.lt.radius_2a_sphere) reg_2a_n_ions=epe(l3)%u
    enddo

    !search fixed ions within I sphere and sent them to II sphere
    n_fixed=0
    do i=1,reg_i_n_ions
       if(epe(i)%fix) then
          j=reg_i_n_ions-n_fixed
          epe_buf=epe(i)
          epe(i)=epe(j)
          epe(j)=epe_buf
          n_fixed=n_fixed+1
       end if
       if(i==reg_i_n_ions-n_fixed)exit
    end do
    reg_i_n_ions=reg_i_n_ions-n_fixed

    if(epeit_count == 0 .and. do_print) &
         write(output_epe,26) reg_i_n_ions,reg_2a_n_ions
26  format(1x,'number of ions in area 1=',i4,', 1 & 2a=',i4)
    if (nmvmax.gt.reg_i_n_ions.and.periodic_optimization)  stop ' reg_i_n_ions have to be increased'
    if(reg_i_n_ions.gt.ndr1) call error_handler('*********stat: reg_i_n_ions is too big')

 !  print*,'call building_tet do_3b',do_3b
    if(n_types_central_atoms_3body > 0 .and. do_3b) call building_tet()
 !  print*,'done', option_c3_symm

#if 1
    if(option_c3_symm) then
       if(allocated(permut_ind)) then
          deallocate(permut_ind,stat=status)
          ASSERT(status==0)
       endif
       allocate(permut_ind(reg_2a_n_ions,3),stat=alloc_status)
       if(alloc_status.ne.0)  then
          call error_handler('failed allocate permut_ind')
       endif
       do i=1,reg_2a_n_ions
          position=epe(i)%r
!!$       if(epeit_count == 0)   print*,'i', i,epe(i)%sym_st_ind
          jj=0
          do j=1,3
             position2=matmul(c3(j)%rotmat,position)
!!$             do k=epe(i)%sym_st_ind,reg_2a_n_ions
             do k=1,reg_2a_n_ions
                if(dot_product(position2-epe(k)%r,position2-epe(k)%r).lt.0.9_r8_kind) then
                   permut_ind(i,j)=k
                   jj=jj+1
!!$print*,'permut_ind',i,j,permut_ind(i,j)
                   exit
                endif
             enddo
          enddo
          if (jj < 3) call error_handler("Mol_module: Your epe environment has no C3 symmetry")
       enddo
       allocate(r_buf(3,reg_2a_n_ions),stat=status)
       if(status /=0 ) call error_handler("Mol_module: allocation of r_buf failed")
       r_buf=zero

       do i=1,reg_2a_n_ions
          position=epe(i)%r
!!$print*,'permut_ind',i,permut_ind(i,:)
          do j=1,3
             position2=matmul(c3(j)%rotmat,position)
             r_buf(:,permut_ind(i,j))=r_buf(:,permut_ind(i,j))+position2
          enddo
       enddo
       epe(:reg_2a_n_ions)%r(1)=r_buf(1,:reg_2a_n_ions)/3.0_r8_kind
       epe(:reg_2a_n_ions)%r(2)=r_buf(2,:reg_2a_n_ions)/3.0_r8_kind
       epe(:reg_2a_n_ions)%r(3)=r_buf(3,:reg_2a_n_ions)/3.0_r8_kind

       r_buf=zero
       do i=1,reg_2a_n_ions
          position=epe(i)%s
          do j=1,3
             position2=matmul(c3(j)%rotmat,position)
             r_buf(:,permut_ind(i,j))=r_buf(:,permut_ind(i,j))+position2
          enddo
       enddo
       epe(:reg_2a_n_ions)%s(1)=r_buf(1,:reg_2a_n_ions)/3.0_r8_kind
       epe(:reg_2a_n_ions)%s(2)=r_buf(2,:reg_2a_n_ions)/3.0_r8_kind
       epe(:reg_2a_n_ions)%s(3)=r_buf(3,:reg_2a_n_ions)/3.0_r8_kind

       r_buf=zero
       do i=1,reg_2a_n_ions
          position=epe(i)%c
          do j=1,3
             position2=matmul(c3(j)%rotmat,position)
             r_buf(:,permut_ind(i,j))=r_buf(:,permut_ind(i,j))+position2
          enddo
       enddo
       epe(:reg_2a_n_ions)%c(1)=r_buf(1,:reg_2a_n_ions)/3.0_r8_kind
       epe(:reg_2a_n_ions)%c(2)=r_buf(2,:reg_2a_n_ions)/3.0_r8_kind
       epe(:reg_2a_n_ions)%c(3)=r_buf(3,:reg_2a_n_ions)/3.0_r8_kind

       deallocate(r_buf,stat=status)
       if(status /=0 ) call error_handler("Mol_module: deallocation of r_buf failed")
    endif
!!$    do i=1,reg_2a_n_ions
!!$       print*,epe(i)%r,i
!!$       print*,epe(i)%s,i
!!$       print*,epe(i)%c,i
!!$    end do
#endif

 !print*,' epe_generater return'
  contains

    subroutine building_tet
! **procedure looks for tetrahedrally coordinated atoms
! **and their neighbours

      integer(kind=i4_kind) :: i,i1,j,k,k1,l,m,n,ind
      integer(kind=i4_kind) :: astatus
      logical :: exit_cycle
      real(kind=r8_kind) :: dist
      real(kind=r8_kind), dimension(4) :: buf_dist
      integer(kind=i4_kind), dimension(4) :: buf_index
      integer(kind=i4_kind) :: max_ind(1)
      real(kind=r8_kind) :: max_val

      if (allocated(tetra_atoms)) then
          deallocate(tetra_atoms,stat=astatus)
          ASSERT(astatus==0)
      endif
      ! search of a number of the central atoms(ions) into 2a region
      ! and consiquently a number of the tetrahedrons
      n_tetrahedrons=0
      l1: do i=1,reg_2a_n_ions
!!$         if (periodic_optimization) then
!!$            do j=1,n_ions_cell
!!$               if(which_epe_ion(j)%new == i) goto 1
!!$            enddo
!!$            cycle l1
!!$         endif
1        k=epe(i)%k
         do j=1,n_types_central_atoms_3body
            if(k == types(j,1)) then
               n_tetrahedrons=n_tetrahedrons+1
               exit
            endif
         enddo
      enddo l1

      allocate(tetra_atoms(n_tetrahedrons,5),stat=status)
      if(status/=0) then
         call error_handler('building_tet:allocation TETRA_ATOMS is failed')
      endif
      tetra_atoms=0

      i1=0
      lab1: do i=1,reg_2a_n_ions
!!$         if (periodic_optimization) then
!!$            do j=1,n_ions_cell
!!$               if(which_epe_ion(j)%new == i) goto 2
!!$            enddo
!!$            cycle lab1
!!$         endif

2        k=epe(i)%k

         lab2: do j=1,n_types_central_atoms_3body
            if(k == types(j,1)) then
               exit_cycle=.false.
               ind=j
               exit lab2
            else
               exit_cycle=.true.
            endif
         enddo lab2

         if (exit_cycle) cycle lab1

         i1=i1+1
         tetra_atoms(i1,1)=i

         buf_dist=0.0_r8_kind
         buf_index=0
         lab3: do l=1,n_gen_ions
            if (l==i) cycle lab3

            k1=epe(l)%k
            lab4: do m=2,5
               if (k1 == types(ind,m)) then
                  exit_cycle=.false.
                  exit lab4
               else
                  exit_cycle=.true.
               endif
            enddo lab4

            if (exit_cycle) cycle lab3

            dist=sqrt(dot_product(epe(i)%r-epe(l)%r,epe(i)%r-epe(l)%r))
            if(dist > r3b(ind)) cycle lab3

            lab5: do n=1,4
               if(buf_dist(n) == 0.0_r8_kind) then
                  buf_dist(n)=dist
                  buf_index(n)=l
                  cycle lab3
               endif
            enddo lab5
            max_ind=maxloc(buf_dist)
            max_val=maxval(buf_dist)
            if(dist < max_val) then
               buf_dist(max_ind)=dist
               buf_index(max_ind)=l
            endif
         enddo lab3
         tetra_atoms(i1,2:5)=buf_index
!!$         do n=2,5
!!$            if(tetra_atoms(i1,n) == 0) then
!!$print*,tetra_atoms(i1,:),'!!!!!!!!'
!!$               call error_handler("The program cannot find complite list of atoms &
!!$                    & to calculate 3-body interaction. Please correct input variable &
!!$                    & R_3B into NAMELIST - three_body_interaction")
!!$            endif
!!$         enddo
      enddo lab1
      first_ind=1
      last_ind=n_tetrahedrons

#if 0
      do k=1,n_tetrahedrons
         i1=tetra_atoms(k,1)
         k1=tetra_atoms(k,2)
         l=tetra_atoms(k,3)
         m=tetra_atoms(k,4)
         n=tetra_atoms(k,5)
         print*,tetra_atoms(k,:),' tetra_atoms '
         print*,k,epe(i1)%k,epe(k1)%k,epe(l)%k,epe(m)%k,epe(n)%k
      enddo
#endif

    end subroutine building_tet

  end subroutine epe_generater
!---------------------------------------------
!---------------------------------------------
  subroutine rot(u,nx,ny,nz,d)
!********************************************c
!  rotation around axes (nx,ny,nz) on d      c
!********************************************c
    real(kind=r8_kind),intent(out)::  u(3,3)
    real(kind=r8_kind), intent(in):: nx,ny,nz,d
    real(kind=r8_kind) :: cosd,sind

    cosd=cos(d)
    sind=sin(d)
    u(1,1)=cosd + (1.0_r8_kind-cosd)*nx*nx
    u(2,1)=(1.0_r8_kind-cosd)*nx*ny + sind*nz
    u(3,1)=(1.0_r8_kind-cosd)*nx*nz - sind*ny
    u(1,2)=(1.0_r8_kind-cosd)*nx*ny -sind*nz
    u(2,2)=cosd + (1.0_r8_kind-cosd)*ny*ny
    u(3,2)=(1.0_r8_kind-cosd)*ny*nz + sind*nx
    u(1,3)=(1.0_r8_kind-cosd)*nx*nz + sind*ny
    u(2,3)=(1.0_r8_kind-cosd)*ny*nz - sind*nx
    u(3,3)=cosd + (1.0_r8_kind-cosd)*nz*nz
  end subroutine rot

end module mol_module
