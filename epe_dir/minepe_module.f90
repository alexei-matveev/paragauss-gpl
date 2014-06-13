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
module minepe_module

use type_module
use epecom_module, only:output_epe

implicit none
private
save

real(kind=r8_kind), public :: DELTA
public :: lattice_gradients_gopt_o
!---------------------------------------------------------
contains


SUBROUTINE lattice_gradients_gopt_o(XS,YS,ZS,XC,YC,ZC,DSX,DSY,DSZ, &
                 DCX,DCY,DCZ,DX,DY,DZ)
! **calculation of relaxation gradients of lattice ions 
! **(DS - shell, DC - core)

  use epecom_module
  use culon_module

  real(kind=r8_kind),intent(out), &
                    dimension(reg_I_n_ions) :: DSX,DSY,DSZ,DCX,DCY,DCZ 
  real(kind=r8_kind),intent(out),dimension(ndr1) :: DX,DY,DZ 
  real(kind=r8_kind),intent(in),dimension(n_gen_ions) :: XS,YS,ZS,XC,YC,ZC
  real(kind=r8_kind),dimension(3) :: gs,gc,gpsis,gpsic
  real(kind=r8_kind) :: DIST,F1,ERF1,RRR2,RSSX,RSSY,RSSZ,CSIX,CSIY,CSIZ
  real(kind=r8_kind) :: RSS2,RSSCSI,RSS1,RSSR2,RSSR8,RSSD,ESS1,ESS2,RSSRO
  real(kind=r8_kind) :: BRO,OSO,OSS,RSCX,RSCY,RSCZ,RCSX,RCSY,RCSZ,RCCX,RCCY
  real(kind=r8_kind) :: RCCZ,RSRX,RSRY,RSRZ,RCRX,RCRY,RCRZ,CCIX,CCIY,CCIZ
  real(kind=r8_kind) :: RSC2,RCS2,RCC2,RSR2,RCR2,RCSCSI,RSCCCI,RCCCCI,RSS3
  real(kind=r8_kind) :: RCS3,RSC3,RCC3,RSR1,RCR1,OSC,OSR,OCS,OCC,OCR,OSSCSI
  real(kind=r8_kind) :: OCSCSI,OSCCCI,OCCCCI,DISX,DISY,DISZ,DICX,DICY,DICZ
  real(kind=r8_kind) :: DIS3X,DIS3Y,DIS3Z,DIC3X,DIC3Y,DIC3Z,RRR(3)
  real(kind=r8_kind) :: RRSX,RRSY,RRSZ,RRCX,RRCY,RRCZ,RRS2,RRC2,RRS1,RRC1
  real(kind=r8_kind) :: OSX,OSY,OSZ,ORS,ORC,OSSX,OSSY,OSSZ,OSCX,OSCY,OSCZ
  real(kind=r8_kind) :: OCSX,OCSY,OCSZ,OCCX,OCCY,OCCZ,CMCSX,CMCSY,CMCSZ
  real(kind=r8_kind) :: CSIS2,CSIC2,ET2S,ET2C,CSIS,OEDS,CSIC,OEDC
  real(kind=r8_kind) :: GIS2,GIC2,PSIS,PSIC,PI2S,PI2C,OEPKX,OEPKY,OEPKZ
  real(kind=r8_kind) :: est,psi,gig,res,gcel,gcelpi,qns,qnc,prsin,prcos
  integer(kind=i4_kind) :: i,io,ig,jo,k,n,j,ic


  DIST=7.d0
  F1=0.026d0
  ERF1=4.d0/3.d0*ERROR_FUNCTION_PARAMETER**3/PIS
  DX=zero
  DY=zero
  DZ=zero
  DSX=zero
  DSY=zero
  DSZ=zero
  DCX=zero
  DCY=zero
  DCZ=zero

  est=zero
  do io=1,n_ions_cell 
    i=which_epe_ion(io)%new
    K=epe(I)%k
    DO J=reg_I_n_ions+1,reg_2a_n_ions   ! contrib due to reg. II discplacements
      if(epe(j)%m.eq.io) cycle
        RRR2=dot_product(epe(I)%r-epe(j)%r,epe(I)%r-epe(j)%r)
      IF(RRR2.GT.RADIUS_LONG_INTERACT) cycle
      N=epe(J)%k
      RSSX=XS(I)-XS(J)
      RSSY=YS(I)-YS(J)
      RSSZ=ZS(I)-ZS(J)
      CSIX=XS(J)-epe(J)%r(1)
      CSIY=YS(J)-epe(J)%r(2)
      CSIZ=ZS(J)-epe(J)%r(3)
      RSS2=RSSX**2+RSSY**2+RSSZ**2
      RSSCSI=RSSX*CSIX+RSSY*CSIY+RSSZ*CSIZ
      RSS1=SQRT(RSS2)
      if(RRR2.GT.16.0_r8_kind) then
         ic=1
      else
         ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
         if(host%ro(K,N,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
      endif
      IF(RRR2.le.host%sr1(K,N,ic)**2) then
        RSSR2=1.0_r8_kind/RSS2
        RSSR8=(RSSR2**2)**2
        RSSD=RSSR2*host%D(K,N,ic)
        ESS1=(6.0_r8_kind *host%C(K,N,ic) + 8.0_r8_kind*RSSD)*RSSR8
        ESS2=(48.0_r8_kind*host%C(K,N,ic) +80.0_r8_kind*RSSD)*(RSSR8*RSSR2)
        IF(RRR2.le.host%sr2(K,N,ic)**2) then
          RSSRO=1.0_r8_kind/(host%RO(K,N,ic)*RSS1)
          BRO=host%B(K,N,ic)*RSSRO*EXP(-RSS1/host%RO(K,N,ic))
          ESS1=-BRO+ESS1
          ESS2= BRO*(RSSR2+RSSRO)-ESS2
        endif 
        ESS2=0.5*ESS2*RSSCSI
        OSO=0.5*ESS1
        OSS=ESS1+ESS2

        DX(I)=DX(I)+OSS*RSSX+OSO*CSIX
        DY(I)=DY(I)+OSS*RSSY+OSO*CSIY
        DZ(I)=DZ(I)+OSS*RSSZ+OSO*CSIZ
      endif ! (RRR2.le.T

      IF(RRR2.GT.RADIUS_LONG_INTERACT) cycle
      RSCX=XS(I)-XC(J)
      RSCY=YS(I)-YC(J)
      RSCZ=ZS(I)-ZC(J)

      RCSX=XC(I)-XS(J)
      RCSY=YC(I)-YS(J)
      RCSZ=ZC(I)-ZS(J)

      RCCX=XC(I)-XC(J)
      RCCY=YC(I)-YC(J)
      RCCZ=ZC(I)-ZC(J)

      RSRX=XS(I)-epe(J)%r(1)
      RSRY=YS(I)-epe(J)%r(2)
      RSRZ=ZS(I)-epe(J)%r(3)

      RCRX=XC(I)-epe(J)%r(1)
      RCRY=YC(I)-epe(J)%r(2)
      RCRZ=ZC(I)-epe(J)%r(3)

      CCIX=XC(J)-epe(J)%r(1)
      CCIY=YC(J)-epe(J)%r(2)
      CCIZ=ZC(J)-epe(J)%r(3)

      RSS2=RSSX**2+RSSY**2+RSSZ**2
      RSC2=RSCX**2+RSCY**2+RSCZ**2
      RCS2=RCSX**2+RCSY**2+RCSZ**2
      RCC2=RCCX**2+RCCY**2+RCCZ**2
      RSR2=RSRX**2+RSRY**2+RSRZ**2
      RCR2=RCRX**2+RCRY**2+RCRZ**2

      RCSCSI=RCSX*CSIX+RCSY*CSIY+RCSZ*CSIZ
      RSCCCI=RSCX*CCIX+RSCY*CCIY+RSCZ*CCIZ
      RCCCCI=RCCX*CCIX+RCCY*CCIY+RCCZ*CCIZ

      RSS3=RSS2*RSS1
      RCS3=RCS2*SQRT(RCS2)
      RSC3=RSC2*SQRT(RSC2)
      RCC3=RCC2*SQRT(RCC2)
      RSR1=SQRT(RSR2)
      RCR1=SQRT(RCR2)

      OSS=-(Q_SHELL(K)*Q_SHELL(N))/RSS3
      OSC=-(Q_SHELL(K)*Q_NUCLEAR(N))/RSC3
      OSR=(Q_SHELL(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1+ERFO*EXP(-ET2*RSR2))/RSR2
      OCS=-(Q_NUCLEAR(K)*Q_SHELL(N))/RCS3
      OCC=-(Q_NUCLEAR(K)*Q_NUCLEAR(N))/RCC3
      OCR=(Q_NUCLEAR(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1+ERFO*EXP(-ET2*RCR2))/RCR2
      OSSCSI=-3.*OSS*RSSCSI/RSS2
      OCSCSI=-3.*OCS*RCSCSI/RCS2
      OSCCCI=-3.*OSC*RSCCCI/RSC2
      OCCCCI=-3.*OCC*RCCCCI/RCC2

      DISX=OSS*RSSX+OSC*RSCX+OSR*RSRX
      DISY=OSS*RSSY+OSC*RSCY+OSR*RSRY
      DISZ=OSS*RSSZ+OSC*RSCZ+OSR*RSRZ
      DICX=OCS*RCSX+OCC*RCCX+OCR*RCRX
      DICY=OCS*RCSY+OCC*RCCY+OCR*RCRY
      DICZ=OCS*RCSZ+OCC*RCCZ+OCR*RCRZ
      DIS3X=OSS*CSIX+OSC*CCIX+OSCCCI*RSCX+OSSCSI*RSSX
      DIS3Y=OSS*CSIY+OSC*CCIY+OSCCCI*RSCY+OSSCSI*RSSY
      DIS3Z=OSS*CSIZ+OSC*CCIZ+OSCCCI*RSCZ+OSSCSI*RSSZ
      DIC3X=OCS*CSIX+OCC*CCIX+OCSCSI*RCSX+OCCCCI*RCCX
      DIC3Y=OCS*CSIY+OCC*CCIY+OCSCSI*RCSY+OCCCCI*RCCY
      DIC3Z=OCS*CSIZ+OCC*CCIZ+OCSCSI*RCSZ+OCCCCI*RCCZ
      DSX(I)=DSX(I)+DISX+0.5*DIS3X
      DSY(I)=DSY(I)+DISY+0.5*DIS3Y
      DSZ(I)=DSZ(I)+DISZ+0.5*DIS3Z
      DCX(I)=DCX(I)+DICX+0.5*DIC3X
      DCY(I)=DCY(I)+DICY+0.5*DIC3Y
      DCZ(I)=DCZ(I)+DICZ+0.5*DIC3Z
    enddo !J=NN2A,NK2A  
! **done, contrib. due to region IIa 

! **treat contrib. from polarized region I
    DO J=1,reg_I_n_ions 
      if(epe(j)%m.eq.io.or.j.eq.i)  cycle
      N=epe(J)%k
        rrr=epe(i)%r-epe(j)%r

      RSSX=XS(I)-XS(J)
      RSSY=YS(I)-YS(J)
      RSSZ=ZS(I)-ZS(J)

      RSCX=XS(I)-XC(J)
      RSCY=YS(I)-YC(J)
      RSCZ=ZS(I)-ZC(J)

      RCSX=XC(I)-XS(J)
      RCSY=YC(I)-YS(J)
      RCSZ=ZC(I)-ZS(J)

      RCCX=XC(I)-XC(J)
      RCCY=YC(I)-YC(J)
      RCCZ=ZC(I)-ZC(J)

      RSRX=XS(I)-epe(J)%r(1)
      RSRY=YS(I)-epe(J)%r(2)
      RSRZ=ZS(I)-epe(J)%r(3)

      RCRX=XC(I)-epe(J)%r(1)
      RCRY=YC(I)-epe(J)%r(2)
      RCRZ=ZC(I)-epe(J)%r(3)

      RRSX=epe(I)%r(1)-XS(J)
      RRSY=epe(I)%r(2)-YS(J)
      RRSZ=epe(I)%r(3)-ZS(J)

      RRCX=epe(I)%r(1)-XC(J)
      RRCY=epe(I)%r(2)-YC(J)
      RRCZ=epe(I)%r(3)-ZC(J)
      RRR2=dot_product(rrr,rrr)
      RSS2=RSSX**2+RSSY**2+RSSZ**2
      RSC2=RSCX**2+RSCY**2+RSCZ**2
      RCS2=RCSX**2+RCSY**2+RCSZ**2
      RCC2=RCCX**2+RCCY**2+RCCZ**2
      RSR2=RSRX**2+RSRY**2+RSRZ**2
      RCR2=RCRX**2+RCRY**2+RCRZ**2
      RRS2=RRSX**2+RRSY**2+RRSZ**2
      RRC2=RRCX**2+RRCY**2+RRCZ**2

      RSS1=SQRT(RSS2)
      RSS3=RSS2*RSS1
      RCC3=RCC2*SQRT(RCC2)
      RSC3=RSC2*SQRT(RSC2)
      RCS3=RCS2*SQRT(RCS2)
      RSR1=SQRT(RSR2)
      RCR1=SQRT(RCR2)
      RRS1=SQRT(RRS2)
      RRC1=SQRT(RRC2)

! **short region interacting part of area 1 - start
      if(RRR2.GT.16.0_r8_kind) then
         ic=1
      else
         ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
         if(host%ro(K,N,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
      endif
      IF(RRR2.le.host%sr1(K,N,ic)**2) then
        RSSR2=1.0_r8_kind/RSS2
        ESS1=(6.0_r8_kind*host%C(K,N,ic)+8.0_r8_kind*host%D(K,N,ic)*RSSR2)*((RSSR2**2)**2)
        IF(RRR2.LE.host%sr2(K,N,ic)**2) ESS1=ESS1- &
        host%B(K,N,ic)/(host%RO(K,N,ic)*RSS1)*EXP(-RSS1/host%RO(K,N,ic))
        OSX=ESS1*RSSX
        OSY=ESS1*RSSY
        OSZ=ESS1*RSSZ
        DX(I)=DX(I)+OSX
        DY(I)=DY(I)+OSY
        DZ(I)=DZ(I)+OSZ
      endif 
! **short region interacting part of area 1 - end

      OSS=-(Q_SHELL(K)*Q_SHELL(N))/RSS3
      OSC=-(Q_SHELL(K)*Q_NUCLEAR(N))/RSC3
      OCS=-(Q_NUCLEAR(K)*Q_SHELL(N))/RCS3
      OCC=-(Q_NUCLEAR(K)*Q_NUCLEAR(N))/RCC3
      IF(ERROR_FUNCTION_PARAMETER*RSR1.LT.DIST) then  
        OSR=(Q_SHELL(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1+ERFO*EXP(-ET2*RSR2))/RSR2
        OCR=(Q_NUCLEAR(K)*Q_ion(epe(j)%k))*(ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1+ERFO*EXP(-ET2*RCR2))/RCR2
        ORS=(Q_ion(epe(i)%k)*Q_SHELL(N))*(ERF(ERROR_FUNCTION_PARAMETER*RRS1)/RRS1+ERFO*EXP(-ET2*RRS2))/RRS2
        ORC=(Q_ion(epe(i)%k)*Q_NUCLEAR(N))*(ERF(ERROR_FUNCTION_PARAMETER*RRC1)/RRC1+ERFO*EXP(-ET2*RRC2))/RRC2
      else
        OSR=(Q_SHELL(K)*Q_ion(epe(j)%k))/(RSR2*RSR1)
        OCR=(Q_NUCLEAR(K)*Q_ion(epe(j)%k))/(RCR2*RCR1)
        ORS=(Q_ion(epe(i)%k)*Q_SHELL(N))/(RRS2*RRS1)
        ORC=(Q_ion(epe(i)%k)*Q_NUCLEAR(N))/(RRC2*RRC1)
      endif !ERROR_FUNCTION_PARAMETER*RSR1.LT.DIST/else

      OSSX=OSS*RSSX
      OSSY=OSS*RSSY
      OSSZ=OSS*RSSZ
      OSCX=OSC*RSCX
      OSCY=OSC*RSCY
      OSCZ=OSC*RSCZ
      OCSX=OCS*RCSX
      OCSY=OCS*RCSY
      OCSZ=OCS*RCSZ
      OCCX=OCC*RCCX
      OCCY=OCC*RCCY
      OCCZ=OCC*RCCZ
      DSX(I)=DSX(I)+OSSX+OSCX+OSR*RSRX
      DSY(I)=DSY(I)+OSSY+OSCY+OSR*RSRY
      DSZ(I)=DSZ(I)+OSSZ+OSCZ+OSR*RSRZ
      DCX(I)=DCX(I)+OCSX+OCCX+OCR*RCRX
      DCY(I)=DCY(I)+OCSY+OCCY+OCR*RCRY
      DCZ(I)=DCZ(I)+OCSZ+OCCZ+OCR*RCRZ
    enddo !J=1,NNA
! **done polirized region I

    CSIX=XS(I)-epe(I)%r(1)
    CSIY=YS(I)-epe(I)%r(2)
    CSIZ=ZS(I)-epe(I)%r(3)
    CCIX=XC(I)-epe(I)%r(1)
    CCIY=YC(I)-epe(I)%r(2)
    CCIZ=ZC(I)-epe(I)%r(3)
    CMCSX=XC(I)-XS(I)
    CMCSY=YC(I)-YS(I)
    CMCSZ=ZC(I)-ZS(I)
    CSIS2=CSIX**2+CSIY**2+CSIZ**2
    CSIC2=CCIX**2+CCIY**2+CCIZ**2
    ET2S=ET2*CSIS2
    ET2C=ET2*CSIC2

    IF(ET2S.GT.F1) then
      CSIS=SQRT(CSIS2)
      OEDS=(ERF(ERROR_FUNCTION_PARAMETER*CSIS)/CSIS+ERFO*EXP(-ET2S))/CSIS2
    else
      OEDS=ERF1*(1.-ET2S*(0.6-0.2142856*ET2S))
    endif ! ET2S.GT.F1)/else
    
    IF(ET2C.GT.F1) then
      CSIC=SQRT(CSIC2)
      OEDC=(ERF(ERROR_FUNCTION_PARAMETER*CSIC)/CSIC+ERFO*EXP(-ET2C))/CSIC2
    else
      OEDC=ERF1*(1.-ET2C*(0.6-0.2142856*ET2C))
    endif       ! ET2C.GT.F1

    OEDS=Q_SHELL(K)*Q_ion(epe(i)%k)*OEDS
    OEDC=Q_NUCLEAR(K)*Q_ion(epe(i)%k)*OEDC

    gs(:)=zero
    gc(:)=zero
    psi=zero
    DO IG=1,n_bs_points
      gig=dot_product(gstr(ig,:),gstr(ig,:))
      res=exp(-pta2*gig)/gig*qpivc
      GIS2=PI2*(GSTR(IG,1)*XS(I)+GSTR(IG,2)*YS(I)+GSTR(IG,3)*ZS(I))
      GIC2=PI2*(GSTR(IG,1)*XC(I)+GSTR(IG,2)*YC(I)+GSTR(IG,3)*ZC(I))
      PSIS=RSIN(IG)*COS(GIS2)-RCOS(IG)*SIN(GIS2)
      PSIC=RSIN(IG)*COS(GIC2)-RCOS(IG)*SIN(GIC2)
      gs(:)=gs(:)+psis*gstr(ig,:)
      gc(:)=gc(:)+psic*gstr(ig,:)
      psi=psi+RSIN(IG)*SIN(GIC2)+RCOS(IG)*COS(GIC2)
      do jo=1,n_ions_cell       ! treat region one
        j=which_epe_ion(jo)%new
        gcel=dot_product(gstr(ig,:),R_ION_IN_CELL(jo,:))
        gcelpi=pi2*gcel
        qns=q_z(jo)*sin(gcelpi)
        qnc=q_z(jo)*cos(gcelpi)
        prsin=res*qns
        prcos=res*qnc
        gpsis(:)=gstr(ig,:)*(prcos*sin(gis2)-prsin*cos(gis2))
        gpsic(:)=gstr(ig,:)*(prcos*sin(gic2)-prsin*cos(gic2))
        DSX(j)=DSX(j)+gpsis(1)*PI*q_shell(k)
        DSY(j)=DSY(j)+gpsis(2)*PI*q_shell(k)
        DSZ(j)=DSZ(j)+gpsis(3)*PI*q_shell(k)
        DCX(j)=DCX(j)+gpsic(1)*PI*q_nuclear(k)
        DCY(j)=DCY(j)+gpsic(2)*PI*q_nuclear(k)
        DCZ(j)=DCZ(j)+gpsic(3)*PI*q_nuclear(k)
      enddo ! jo=1,n_ions_cell
    enddo !IG=1,n_bs_points

    PI2S=PI*Q_SHELL(K)  ! divided by two for gradients as for energy
    PI2C=PI*Q_NUCLEAR(K)
    est=est+psi*Q_NUCLEAR(K)/two
    OEPKX=PK(K)*CMCSX
    OEPKY=PK(K)*CMCSY
    OEPKZ=PK(K)*CMCSZ

    DSX(I)=DSX(I)+CSIX*OEDS+GS(1)*PI2S-OEPKX
    DSY(I)=DSY(I)+CSIY*OEDS+GS(2)*PI2S-OEPKY
    DSZ(I)=DSZ(I)+CSIZ*OEDS+GS(3)*PI2S-OEPKZ
    DCX(I)=DCX(I)+OEPKX+CCIX*OEDC+GC(1)*PI2C
    DCY(I)=DCY(I)+OEPKY+CCIY*OEDC+GC(2)*PI2C
    DCZ(I)=DCZ(I)+OEPKZ+CCIZ*OEDC+GC(3)*PI2C
  enddo !io=1,n_ions_cell


END SUBROUTINE lattice_gradients_gopt_o
!----------------------------------------------------
!----------------------------------------------------
end module minepe_module
