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
module culon_module

#include <def.h>
use type_module
use epecom_module,only:output_epe

implicit none
save
private

public :: energy_coulomb_po, energy_coulomb, energy_shortrange_po, &
          energy_shortrange, erf, erfc, minv, gauss_potential, &
          energy_3_body_po,energy_3_body,energy_vac_imp,common_atom
!----------------------------------------------------------
!----------------------------------------------------------
!
contains

  function common_atom(iatom1,iatom2)

    use epecom_module, only: n_tetrahedrons,epe,tetra_atoms

    integer(i4_kind) :: common_atom
    integer(i4_kind), intent(in) :: iatom1,iatom2
    integer(i4_kind) :: i,j,k
    logical :: yes1,yes2

    common_atom=0
    yes1=.false.; yes2=.false.
    l_i: do i=1,n_tetrahedrons
       l_j: do j=2,4
          if(tetra_atoms(i,j) /= iatom1 .and. tetra_atoms(i,j) /= iatom2) cycle l_j
          l_k: do k=j+1,5
             if(tetra_atoms(i,k) /= iatom1 .and. tetra_atoms(i,k) /= iatom2) cycle l_k
             common_atom = epe(tetra_atoms(i,1))%k
             exit l_i
          end do l_k
       end do l_j
    end do l_i

  end function common_atom
!----------------------------------------------------------
!----------------------------------------------------------

real(kind=r8_kind) function energy_coulomb_po(l_print)
! **program of the cell COULOMB part relaxation energy calc.
! **for periodic optimization

  use epecom_module,excoup=>ec

  logical, intent(in) :: l_print

  real(kind=r8_kind) :: E,EPK,e1,e2,e3,rrr2,rss2,rcs2,rsc2,rcc2
  real(kind=r8_kind) :: rsr2,rcr2,rrs2,rrc2,rsscsi,rcscsi,rrscsi
  real(kind=r8_kind) :: rsccci,rcccci,rrccci,rssa,rcs,rsc,rcca
  real(kind=r8_kind) :: cci,csi,rrr1,rss1,rsc1
  real(kind=r8_kind) :: rcs1,rcc1,rsr1,rcr1,rrs1,rrc1,rss3,rcs3
  real(kind=r8_kind) :: rsc3,rcc3,rrs3,rrc3,ea,eb,ec,csics,csic,csis
  real(kind=r8_kind) :: eprm,edc,eds,ed,eg,esum
  integer(kind=i4_kind) :: IO,i,k,j,j1,n,ind

  E=zero
  EPK=zero
  DO IO=1,n_ions_cell
    i=which_epe_ion(io)%new
    K=epe(I)%k

! **area 1 - start
! **area 2a - start

    E2=zero
    E3=zero
! **treat region 2a
    DO J1=reg_I_n_ions+1,reg_2a_n_ions
      J=reg_I_n_ions+1+reg_2a_n_ions-J1
      N=epe(J)%k
        RRR2=dot_product(epe(I)%r-epe(j)%r,epe(I)%r-epe(j)%r)
      IF(RRR2.GT.RADIUS_LONG_INTERACT) cycle

      RSS2=zero
      RCS2=zero
      RSC2=zero
      RCC2=zero
      RSR2=zero
      RCR2=zero
      RSSCSI=zero
      RCSCSI=zero
      RSCCCI=zero
      RCCCCI=zero

      DO IND=1,3
        RSSA=R_SH_ION(I,IND)-R_SH_ION(J,IND)
        RCS=R_NUC_ION(I,IND)-R_SH_ION(J,IND)
        RSC=R_SH_ION(I,IND)-R_NUC_ION(J,IND)
        RCCA=R_NUC_ION(I,IND)-R_NUC_ION(J,IND)
        CCI=R_NUC_ION(J,IND)-epe(J)%r(IND)
        CSI=R_SH_ION(J,IND)-epe(J)%r(ind)
        RSS2=RSS2+RSSA*RSSA
        RCS2=RCS2+RCS*RCS
        RSC2=RSC2+RSC*RSC
        RCC2=RCC2+RCCA*RCCA
        RSSCSI=RSSCSI+RSSA*CSI
        RCSCSI=RCSCSI+RCS*CSI
        RSCCCI=RSCCCI+RSC*CCI
        RCCCCI=RCCCCI+RCCA*CCI
      enddo !ind=1,3
        RRS2=dot_product(epe(I)%r-R_SH_ION(J,:),epe(I)%r-R_SH_ION(J,:))
        RRSCSI=dot_product(epe(I)%r-R_SH_ION(J,:),R_SH_ION(J,:)-epe(j)%r)
        RRC2=dot_product(epe(I)%r-R_NUC_ION(J,:),epe(I)%r-R_NUC_ION(J,:))
        RRCCCI=dot_product(epe(I)%r-R_SH_ION(J,:),R_NUC_ION(J,:)-epe(j)%r)
        RSR2=dot_product(R_SH_ION(I,:)-epe(j)%r,R_SH_ION(I,:)-epe(j)%r)
        RCR2=dot_product(R_NUC_ION(I,:)-epe(j)%r,R_NUC_ION(I,:)-epe(j)%r)


      RRR1=SQRT(RRR2)
      RSS1=SQRT(RSS2)
      RSC1=SQRT(RSC2)
      RCS1=SQRT(RCS2)
      RCC1=SQRT(RCC2)
      RSR1=SQRT(RSR2)
      RCR1=SQRT(RCR2)
      RRS1=SQRT(RRS2)
      RRC1=SQRT(RRC2)
      RSS3=RSS2*RSS1
      RCS3=RCS2*RCS1
      RSC3=RSC2*RSC1
      RCC3=RCC2*RCC1
      RRS3=RRS2*RRS1
      RRC3=RRC2*RRC1

      if(RRR1.le.epe(reg_2a_n_ions)%d) then
!!!        EA=Q_SHELL(K)*(Q_SHELL(N)/RSS1+Q_NUCLEAR(N)/RSC1-Q_ZL(J)*Q*ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1)
        EA=Q_SHELL(K)*(Q_SHELL(N)/RSS1+Q_NUCLEAR(N)/RSC1-Q_ion(n)*ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1)
        EB=Q_NUCLEAR(K)*(Q_SHELL(N)/RCS1+Q_NUCLEAR(N)/RCC1-Q_ion(n)*ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1)
        EC=-Q_ion(epe(I)%k)*(Q_SHELL(N)/RRS1+Q_NUCLEAR(N)/RRC1-Q_ion(n)/RRR1)
        EPRM=EA+EB
        E2=E2+EPRM
      endif

      EA=Q_SHELL(K)*(Q_SHELL(N)*RSSCSI/RSS3+Q_NUCLEAR(N)*RSCCCI/RSC3)
      EB=Q_NUCLEAR(K)*(Q_SHELL(N)*RCSCSI/RCS3+Q_NUCLEAR(N)*RCCCCI/RCC3)
      EC=-Q_ion(epe(I)%k)*(Q_SHELL(N)*RRSCSI/RRS3+Q_NUCLEAR(N)*RRCCCI/RRC3)
      EPRM=EA+EB
      E3=E3+EPRM
    enddo
! ** done treat region 2a

! **polarized region 1   - start
    E1=zero
    DO J=1,reg_I_n_ions
      if(j.eq.i) cycle
      N=epe(J)%k
      RSS2=zero
      RCC2=zero
      RSC2=zero
      RCS2=zero

      DO IND=1,3
        RSS2=RSS2+(R_SH_ION(I,IND)-R_SH_ION(J,IND))**2
        RCC2=RCC2+(R_NUC_ION(I,IND)-R_NUC_ION(J,IND))**2
        RCS2=RCS2+(R_NUC_ION(I,IND)-R_SH_ION(J,IND))**2
        RSC2=RSC2+(R_SH_ION(I,IND)-R_NUC_ION(J,IND))**2
      enddo !IND=1,3
        RRR2=dot_product(epe(I)%r-epe(J)%r,epe(I)%r-epe(J)%r)
        RRS2=dot_product(epe(I)%r-R_SH_ION(J,:),epe(I)%r-R_SH_ION(J,:))
        RCR2=dot_product(R_NUC_ION(I,:)-epe(J)%r,R_NUC_ION(I,:)-epe(J)%r)
        RRC2=dot_product(epe(I)%r-R_NUC_ION(J,:),epe(I)%r-R_NUC_ION(J,:))
        RSR2=dot_product(R_SH_ION(I,:)-epe(J)%r,R_SH_ION(I,:)-epe(J)%r)

      RRR1=SQRT(RRR2)
      RSS1=SQRT(RSS2)
      RRS1=SQRT(RRS2)
      RSR1=SQRT(RSR2)
      RCC1=SQRT(RCC2)
      RRC1=SQRT(RRC2)
      RCR1=SQRT(RCR2)
      RSC1=SQRT(RSC2)
      RCS1=SQRT(RCS2)

      if(RRR1.le.epe(reg_2a_n_ions)%d) then
        EA= Q_SHELL(epe(i)%k)*(Q_SHELL(N)/RSS1+Q_NUCLEAR(N)/RSC1-Q_ion(epe(J)%k)*ERF(ERROR_FUNCTION_PARAMETER*RSR1)/RSR1)
        EB=Q_NUCLEAR(K)*(Q_SHELL(N)/RCS1+Q_NUCLEAR(N)/RCC1-Q_ion(epe(J)%k)*ERF(ERROR_FUNCTION_PARAMETER*RCR1)/RCR1)
        EC=-Q_ion(epe(I)%k)*(Q_SHELL(N)*ERF(ERROR_FUNCTION_PARAMETER*RRS1)/RRS1+ &
            Q_NUCLEAR(N)*ERF(ERROR_FUNCTION_PARAMETER*RRC1)/RRC1-Q_ion(epe(J)%k)/RRR1)
        EPRM=EA+EB
        E1=E1+EPRM
      endif
    enddo

        CSICS=dot_product(R_NUC_ION(I,:)-R_SH_ION(I,:),R_NUC_ION(I,:)-R_SH_ION(I,:))
        CSIS=dot_product(R_SH_ION(I,:)-epe(i)%r,R_SH_ION(I,:)-epe(i)%r)
        CSIC=dot_product(R_NUC_ION(I,:)-epe(i)%r,R_NUC_ION(I,:)-epe(i)%r)

      IF(CSIC.GT.1.0E-9)  CSIC=SQRT(CSIC)
      IF(CSIS.GT.1.0E-9)  CSIS=SQRT(CSIS)
      EDC=ERFO
      EDS=ERFO
      IF(CSIC.GT.0.0001) EDC=-ERF(ERROR_FUNCTION_PARAMETER*CSIC)/CSIC
      IF(CSIS.GT.0.0001) EDS=-ERF(ERROR_FUNCTION_PARAMETER*CSIS)/CSIS
      ED=Q_ion(epe(I)%k)*(Q_NUCLEAR(epe(I)%k)*EDC+Q_SHELL(epe(I)%k)*EDS)
      EG=Q_SHELL(K)*gauss_potential(R_SH_ION,n_gen_ions,I)+Q_NUCLEAR(K)*gauss_potential(R_NUC_ION,n_gen_ions,I)
      EPK=EPK+0.5*PK(K)*CSICS
      esum=ed+eg

      EPRM=eg+ed+E2-0.5*E3+E1
      E=E+EPRM
    enddo !IO=1,n_ions_cell

! **polarized region 1   - end
    energy_coulomb_po=(e)/two
    if(l_print) then
       write(output_epe,106)E,E2BDIS,EPK,0.d0,E2BIND,energy_coulomb_po
106    FORMAT(1x,'DISPLACEMENT POLARIZATION FOR AREAS: 1 & 2a=',f13.4, &
            /,   1x,'                                         2b=',f13.4, &
            /,1x,      'INDUCED POLARIZATION OF AREAS     :  1=',f13.4, &
            /,1x,      '                                    2a=',f13.4, &
            /,1x,      '                                    2b=',f13.4, &
            /,1x,      'COULOMB ENERGY OF LATTICE         epol=',f13.4)
    endif
END function energy_coulomb_po
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

REAL(kind=r8_kind) FUNCTION energy_coulomb(l_print)
! **program of the cell COULOMB part relaxation energy calc.

  use epecom_module,excoup=>ec

  logical, intent(in) :: l_print

! real(kind=r8_kind) :: anion_dir,sum_mad ! to check direct contrib

  real(kind=r8_kind) :: reg1_coulomb
  real(kind=r8_kind) :: E,EPK
  real(kind=r8_kind) :: csics
  real(kind=r8_kind) :: epk2a
  integer(kind=i4_kind) :: i,j,n
  integer(kind=i4_kind) :: status

  type dot_products
     real(kind=r8_kind)::ss,sc,cs,cc,rs,rc,sr,cr,rr
  end type dot_products
  type (dot_products),allocatable, dimension(:):: sqdp
  real(kind=r8_kind),allocatable, dimension(:):: disp_core,disp_shell

  epe(1:reg_2a_n_ions)%s(1)=R_SH_ION(1:reg_2a_n_ions,1)
  epe(1:reg_2a_n_ions)%s(2)=R_SH_ION(1:reg_2a_n_ions,2)
  epe(1:reg_2a_n_ions)%s(3)=R_SH_ION(1:reg_2a_n_ions,3)
  epe(1:reg_2a_n_ions)%c(1)=R_nuc_ION(1:reg_2a_n_ions,1)
  epe(1:reg_2a_n_ions)%c(2)=R_nuc_ION(1:reg_2a_n_ions,2)
  epe(1:reg_2a_n_ions)%c(3)=R_nuc_ION(1:reg_2a_n_ions,3)

  do j=1,reg_I_n_ions
     epe(j)%qs=q_shell(epe(j)%k)
     epe(j)%qc=q_nuclear(epe(j)%k)
     epe(j)%q=q_ion(epe(j)%k)
  enddo

    allocate(sqdp(reg_I_n_ions),disp_core(reg_I_n_ions),disp_shell(reg_I_n_ions),&
         stat=status)
    if(status.ne.0) call error_handler("dp allocation failed")

    ! **polarized region 1   - start

    reg1_coulomb=0.0_r8_kind
    DO I=1,reg_I_n_ions
!      anion_dir=0.0_r8_kind
      DO J=I+1,reg_I_n_ions
        N=epe(J)%k

        sqdp(i)%ss=sqrt((epe(i)%s(1)-epe(j)%s(1))**2 + &
                        (epe(i)%s(2)-epe(j)%s(2))**2 + &
                        (epe(i)%s(3)-epe(j)%s(3))**2   )

        sqdp(i)%sc=sqrt((epe(i)%s(1)-epe(j)%c(1))**2 + &
                        (epe(i)%s(2)-epe(j)%c(2))**2 + &
                        (epe(i)%s(3)-epe(j)%c(3))**2   )

        sqdp(i)%cc=sqrt((epe(i)%c(1)-epe(j)%c(1))**2 + &
                        (epe(i)%c(2)-epe(j)%c(2))**2 + &
                        (epe(i)%c(3)-epe(j)%c(3))**2   )

        sqdp(i)%cs=sqrt((epe(i)%c(1)-epe(j)%s(1))**2 + &
                        (epe(i)%c(2)-epe(j)%s(2))**2 + &
                        (epe(i)%c(3)-epe(j)%s(3))**2   )

        sqdp(i)%cr=sqrt((epe(i)%c(1)-epe(j)%r(1))**2 + &
                        (epe(i)%c(2)-epe(j)%r(2))**2 + &
                        (epe(i)%c(3)-epe(j)%r(3))**2   )

        sqdp(i)%sr=sqrt((epe(i)%s(1)-epe(j)%r(1))**2 + &
                        (epe(i)%s(2)-epe(j)%r(2))**2 + &
                        (epe(i)%s(3)-epe(j)%r(3))**2   )

        sqdp(i)%rr=sqrt((epe(i)%r(1)-epe(j)%r(1))**2 + &
                        (epe(i)%r(2)-epe(j)%r(2))**2 + &
                        (epe(i)%r(3)-epe(j)%r(3))**2   )

        sqdp(i)%rs=sqrt((epe(i)%r(1)-epe(j)%s(1))**2 + &
                        (epe(i)%r(2)-epe(j)%s(2))**2 + &
                        (epe(i)%r(3)-epe(j)%s(3))**2   )

        sqdp(i)%rc=sqrt((epe(i)%r(1)-epe(j)%c(1))**2 + &
                        (epe(i)%r(2)-epe(j)%c(2))**2 + &
                        (epe(i)%r(3)-epe(j)%c(3))**2   )

! without dicplacements just sum of direct space contributions over all atoms in reg 1
  reg1_coulomb= reg1_coulomb &
       + epe(i)%qs *(Q_SHELL(N)/sqdp(i)%ss + Q_NUCLEAR(N)/sqdp(i)%sc &
                     -Q_ion(epe(j)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%sr)/sqdp(i)%sr &
                                                                                           ) &
       + epe(i)%qc*(Q_SHELL(N)/sqdp(i)%cs + Q_NUCLEAR(N)/sqdp(i)%cc &
                     -Q_ion(epe(J)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%cr)/sqdp(i)%cr &
                                                                                           ) &
       - epe(i)%q*( &
                   Q_SHELL(N)  *ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%rs)/sqdp(i)%rs &
                  +Q_NUCLEAR(N)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%rc)/sqdp(i)%rc &
                  -Q_ion(epe(J)%k)/sqdp(i)%rr )

!     if(i.eq.1) then
!      anion_dir=anion_dir &
!      + epe(i)%qs *(Q_SHELL(N)/sqdp(i)%ss + Q_NUCLEAR(N)/sqdp(i)%sc &
!                    -Q_ion(epe(j)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%sr)/sqdp(i)%sr &
!                                                                                          ) &
!      + epe(i)%qc*(Q_SHELL(N)/sqdp(i)%cs + Q_NUCLEAR(N)/sqdp(i)%cc &
!                    -Q_ion(epe(J)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%cr)/sqdp(i)%cr &
!                                                                                          )

!      if(j.eq.2) print*, &
!                -epe(i)%qs * Q_ion(epe(j)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%sr)/sqdp(i)%sr &
!                -epe(i)%qc * Q_ion(epe(J)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%cr)/sqdp(i)%cr, &
!                 (epe(i)%qs+epe(i)%qc),Q_ion(epe(J)%k),sqdp(i)%sr,'conrtib q q d', &
!                  ERROR_FUNCTION_PARAMETER,ERF(ERROR_FUNCTION_PARAMETER*sqdp(i)%sr)/sqdp(i)%sr
!      if(j.eq.2) print*,epe(i)%s,'i s'
!      if(j.eq.2) print*,epe(j)%r,'j r'
!
!     endif
      enddo
!     if(i.eq.1) then
!       DPRINT 'anion_dir for reg 1 charges', anion_dir,(epe(1)%qs+epe(1)%qc),Q_ion(epe(1)%k)
!     endif
   end DO

!DPRINT 'reg1_coulomb', reg1_coulomb



      EPK=zero
      DO I=1,reg_I_n_ions
         CSICS=dot_product(epe(i)%c-epe(i)%s,epe(i)%c-epe(i)%s)
         EPK=EPK+0.5_r8_kind*PK(epe(I)%k)*CSICS
      end DO

      E=zero

!     lim dicplacement contrins

      sqdp(1:reg_I_n_ions)%cc=sqrt((epe(1:reg_I_n_ions)%c(1)-epe(1:reg_I_n_ions)%r(1))**2 + &
                      (epe(1:reg_I_n_ions)%c(2)-epe(1:reg_I_n_ions)%r(2))**2 + &
                      (epe(1:reg_I_n_ions)%c(3)-epe(1:reg_I_n_ions)%r(3))**2   )
      where(sqdp(:)%cc.gt.0.0001_r8_kind)
            disp_core(:)=-arrerf(ERROR_FUNCTION_PARAMETER*sqdp(:)%cc)/sqdp(:)%cc
         elsewhere
            disp_core(:)= ERFO
      end where

      sqdp(1:reg_I_n_ions)%ss=sqrt((epe(1:reg_I_n_ions)%s(1)-epe(1:reg_I_n_ions)%r(1))**2 + &
                      (epe(1:reg_I_n_ions)%s(2)-epe(1:reg_I_n_ions)%r(2))**2 + &
                      (epe(1:reg_I_n_ions)%s(3)-epe(1:reg_I_n_ions)%r(3))**2   )

      where(sqdp(:)%ss.gt.0.0001_r8_kind)
            disp_shell(:)=-arrerf(ERROR_FUNCTION_PARAMETER*sqdp(:)%ss)/sqdp(:)%ss
         elsewhere
            disp_shell(:)= ERFO
      end where
      e=e+sum(epe(:reg_i_n_ions)%q*(epe(:reg_i_n_ions)%qs*disp_shell+epe(:reg_i_n_ions)%qc*disp_core))

      DPRINT 'where lim dicplacement contribs', E

!     done lim dicplacement contribs



! sum_mad=0.0_r8_kind
  DO I=1,reg_I_n_ions

!!$     CSIC=dot_product(epe(i)%c-epe(i)%r,epe(i)%c-epe(i)%r)
!!$     CSIS=dot_product(epe(i)%s-epe(i)%r,epe(i)%s-epe(i)%r)
!!$
!!$    IF(CSIC.GT.1.0E-9) CSIC=SQRT(CSIC)
!!$    IF(CSIS.GT.1.0E-9) CSIS=SQRT(CSIS)
!!$    EDC=ERFO
!!$    EDS=ERFO
!!$    EDC=disp_core(i)
!!$    EDS=disp_shell(i)
!!$    IF(CSIC.GT.0.0001) EDC=-ERF(ERROR_FUNCTION_PARAMETER*CSIC)/CSIC
!!$    IF(CSIS.GT.0.0001) EDS=-ERF(ERROR_FUNCTION_PARAMETER*CSIS)/CSIS
    E=E+ &
!!$         Q_ion(epe(I)%k)*(Q_NUCLEAR(epe(I)%k)*EDC+Q_SHELL(epe(I)%k)*EDS) +&
         Q_SHELL(epe(I)%k)  *gauss_potential( R_SH_ION,n_gen_ions,I)+&
         Q_NUCLEAR(epe(I)%k)*gauss_potential(R_NUC_ION,n_gen_ions,I) &
        -madc%emd(epe(i)%m)
!  sum_mad=sum_mad-madc%emd(epe(i)%m)
!  if(i.eq.1) then
!   print*,'bs-emd', i, Q_SHELL(  epe(I)%k)*gauss_potential( R_SH_ION,n_gen_ions,I),   &
!                     +Q_NUCLEAR(epe(I)%k)*gauss_potential(R_NUC_ION,n_gen_ions,I)    &
!        ,-madc%emd(epe(i)%m)
!   DPRINT epe(i)%q*(epe(i)%qs*disp_shell(i)+epe(i)%qc*disp_core(i)),'lim err func contrib'
!  endif
  enddo
! DPRINT sum_mad,'sum madc%emd over r'
! DPRINT 'energy local field contributions', E
      deallocate(sqdp,disp_core,disp_shell,stat=status)
      if(status.ne.0) then
         call error_handler("sqdp deallocation failed")
         end if

  E=E+reg2_coulomb()+reg1_coulomb
! DPRINT 'reg1_coulomb reg2_coulomb',reg1_coulomb, reg2_coulomb()

  EPK2A=zero
  DO I=reg_I_n_ions+1,reg_2a_n_ions
    EPK2A=EPK2A+0.5*PK(epe(I)%k)* &
         dot_product(R_NUC_ION(I,:)-R_SH_ION(I,:),R_NUC_ION(I,:)-R_SH_ION(I,:))
  enddo
  E=E -EPK2A+ ECRR

! **polarized region 1   - end

  energy_coulomb=E+EPK+EPK2A+E2BIND+E2BDIS
  DPRINT 'energy_coulomb', energy_coulomb
  if(l_print) then
     write(output_epe, 106)E,E2BDIS,EPK,EPK2A,E2BIND,energy_coulomb
106  FORMAT(1x,'DISPLACEMENT POLARIZATION FOR AREAS: 1 & 2a=',f20.11, &
          /,   1x,'                                         2b=',f13.4, &
          /,1x,      'INDUCED POLARIZATION OF AREAS     :  1=',f13.4, &
          /,1x,      '                                    2a=',f13.4, &
          /,1x,      '                                    2b=',f13.4, &
          /,1x,      'COULOMB ENERGY OF LATTICE         epol=',f13.4)
  endif
END function energy_coulomb
!----------------------------------------------------------------------
!----------------------------------------------------------------------

real(kind=r8_kind) function energy_shortrange_po(l_print)
! **calc. of short range part of the cell relaxation energy for
! **periodic optimization

  use epecom_module,excoup=>ec

  logical, intent(in) :: l_print

  real(kind=r8_kind) :: ea,eb,ec,ed,ex,ey,e1,e2,e3
  real(kind=r8_kind) :: rrr2,rss2,rsscsi,rrs2,rrscsi
  real(kind=r8_kind) :: rrr,rrs,csi,rss6,rrs6,rss8,rrs8
  real(kind=r8_kind) :: rrs1,rss1,expss,exprs,rssa,rrr1
  integer(kind=i4_kind) :: IO,i,k,j,j1,n,ind,ic

  EA=zero
  EB=zero
  EC=zero
  ED=zero
  EX=zero
  EY=zero
  DO IO=1,n_ions_cell
    i=which_epe_ion(io)%new
    K=epe(I)%k
!    write(output_epe,*) 'energy_shortrange_po', which_epe_ion(io)%new ,k

! **area 2a - start
    DO J1=reg_I_n_ions+1,reg_2a_n_ions
      J=reg_I_n_ions+1+reg_2a_n_ions-J1
      if(i.eq.j) cycle
      N=epe(J)%k
      RRR2=zero
      RSS2=zero
      RSSCSI=zero
      RRS2=zero
      RRSCSI=zero
      DO IND=1,3
        RRR=epe(I)%r(IND)-epe(J)%r(ind)
        RRS=epe(I)%r(IND)-R_SH_ION(J,IND)
        RSSA=R_SH_ION(I,IND)-R_SH_ION(J,IND)
        CSI=R_SH_ION(J,IND)-epe(J)%r(ind)
        RRR2=RRR2+RRR*RRR
        RSS2=RSS2+RSSA*RSSA
        RRS2=RRS2+RRS*RRS
        RSSCSI=RSSCSI+RSSA*CSI
        RRSCSI=RRSCSI+RRS*CSI
      enddo!IND=1,3

      if(RRR2.GT.16.0_r8_kind) then
         ic=1
      else
         ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
         if(host%ro(K,N,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
      endif

      IF(RRR2.GT.host%sr1(K,N,ic)**2) cycle
      RSS6=RSS2**3
      RRS6=RRS2**3
      RSS8=RSS6*RSS2
      RRS8=RRS6*RRS2
      EA=EA-host%C(K,N,ic)*(1.0_r8_kind/RSS6)-host%D(K,N,ic)*(1.0_r8_kind/RSS8)
      EC=EC-(6.0_r8_kind*host%C(K,N,ic)/RSS8+8.0_r8_kind*host%D(K,N,ic)/(RSS8*RSS2))*RSSCSI
      IF(RRR2.GT.host%sr2(K,N,ic)**2) cycle
      RRS1=SQRT(RRS2)
      RSS1=SQRT(RSS2)
      EXPSS=EXP(-RSS1/host%RO(K,N,ic))
      EXPRS=EXP(-RRS1/host%RO(K,N,ic))
      EB=EB+host%B(K,N,ic)*(EXPSS)!-EXPRS)
      ED=ED+host%B(K,N,ic)/host%RO(K,N,ic)*(EXPSS*RSSCSI/RSS1)
    enddo
!    write(output_epe,*) ea,eb,ec,ed
! **region  2a finished

! **Treat region I
    DO J=1,reg_I_n_ions
      if(j.eq.i) cycle
      RSS2=0.
      DO IND=1,3
        RSS2=RSS2+(R_SH_ION(I,IND)-R_SH_ION(J,IND))**2
      enddo! IND=1,3
      RRR2=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)

      if(RRR2.GT.16.0_r8_kind) then
         ic=1
      else
         ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
         if(host%ro(K,epe(J)%k,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
      endif

      IF(RRR2.GT.host%sr1(K,epe(J)%k,ic)**2) cycle
      EX=EX-host%C(K,epe(J)%k,ic)*(1.0_r8_kind/(RSS2**3)) &
           -host%D(K,epe(J)%k,ic)*(1.0_r8_kind/(RSS2**4))
      IF(RRR2.GT.host%sr2(K,epe(J)%k,ic)**2) cycle
      RRR1=SQRT(RRR2)
      RSS1=SQRT(RSS2)
      EY=EY+host%B(K,epe(J)%k,ic)*(EXP(-RSS1/host%RO(K,epe(J)%k,ic)))
    enddo
! **region 2a  - end
  enddo!IO=1,n_ions_cell

! **summarizing different parts of energy
   E1=(EX+EY)/two
   E2=(EA+EB)/two
   E3=-0.5*(EC+ED)/two
   energy_shortrange_po=(E1+E2+E3)
   if(l_print) then
      write(output_epe,104)energy_shortrange_po,E1,E2,E3
104   FORMAT(1x,'short range interaction energy     shr=',F15.6, &
           /,1x,'                 of         area    1 =',F15.6, &
           /,1x,'                 interaction of 1 & 2A=',F15.6, &
           /,1x,'                 of         area    2a=',F15.6)
   endif
END function energy_shortrange_po
!----------------------------------------------------------------------
!----------------------------------------------------------------------

REAL(kind=r8_kind) FUNCTION energy_shortrange(l_print)
! **calc. of short range part of the cell relaxation energy

  use epecom_module, excoup=>ec

  logical, intent(in) :: l_print

  real(kind=r8_kind) :: ea,eb,ec,ed,ex,ey,e1,e2,e3
  real(kind=r8_kind) :: rrr2,rss2,rsscsi,rrs2,rrscsi
  real(kind=r8_kind) :: rss6,rrs6,rss8,rrs8
  real(kind=r8_kind) :: rrs1,rss1,expss,exprs,rssa
  real(kind=r8_kind), parameter:: one=1.0_r8_kind
  integer(kind=i4_kind) :: i,k,j,j1,ind,ic



!!$call cpu_time(t1)
!!$print*,'<<<<<<<<<<< energy_shortrange ',t1
  EA=zero
  EB=zero
  EC=zero
  ED=zero
  EX=zero
  EY=zero
  DO I=1,reg_I_n_ions
    K=epe(I)%k
! **area 2a - start
    DO J1=reg_I_n_ions+1,reg_2a_n_ions
      J=reg_I_n_ions+1+reg_2a_n_ions-J1
      RSS2=zero
      DO IND=1,3
        RSSA=R_SH_ION(I,IND)-R_SH_ION(J,IND)
        RSS2=RSS2+RSSA*RSSA
      enddo! IND=1,3
      RSSCSI=dot_product(R_SH_ION(I,:)-R_SH_ION(J,:),R_SH_ION(J,:)-epe(j)%r)
      RRSCSI=dot_product(epe(i)%r-R_SH_ION(J,:),R_SH_ION(J,:)-epe(j)%r)
      RRR2=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)
      RRS2=dot_product(epe(i)%r-R_SH_ION(J,:),epe(i)%r-R_SH_ION(J,:))

      if(RRR2.GT.16.0_r8_kind) then
         ic=1
      else
         ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
         if(host%ro(K,epe(J)%k,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
      endif

      IF(RRR2.GT.host%sr1(epe(I)%k,epe(J)%k,ic)**2) cycle
      RSS6=RSS2**3
      RRS6=RRS2**3
      RSS8=RSS6*RSS2
      RRS8=RRS6*RRS2
      EA=EA-host%C(epe(I)%k,epe(J)%k,ic)*(1.0_r8_kind/RSS6- &
           1.0_r8_kind/RRS6)-host%D(epe(I)%k,epe(J)%k,ic)*(1.0_r8_kind/RSS8 &
           -1.0_r8_kind/RRS8)
      EC=EC-(6.0_r8_kind*host%C(epe(I)%k,epe(J)%k,ic)/RSS8 &
           +8.0_r8_kind*host%D(epe(I)%k,epe(J)%k,ic)/(RSS8*RSS2))*RSSCSI+ &
           (6.0_r8_kind*host%C(epe(I)%k,epe(J)%k,ic)/RRS8 &
           +8.0_r8_kind*host%D(epe(I)%k,epe(J)%k,ic)/(RRS8*RRS2))*RRSCSI
      IF(RRR2.GT.host%sr2(epe(I)%k,epe(J)%k,ic)**2) cycle
      RRS1=SQRT(RRS2)
      RSS1=SQRT(RSS2)
      EXPSS=EXP(-RSS1/host%RO(epe(I)%k,epe(J)%k,ic))
      EXPRS=EXP(-RRS1/host%RO(epe(I)%k,epe(J)%k,ic))
      EB=EB+host%B(epe(I)%k,epe(J)%k,ic)*(EXPSS-EXPRS)
      ED=ED+host%B(epe(I)%k,epe(J)%k,ic)/host%RO(epe(I)%k,epe(J)%k,ic)*(EXPSS*RSSCSI/RSS1-EXPRS*RRSCSI/RRS1)
    enddo
! **region  2a finished

! **Treat region I
      DO J=I+1,reg_I_n_ions
        RRR2=dot_product(epe(I)%r-epe(J)%r,epe(I)%r-epe(J)%r)
        RSS2=dot_product(R_SH_ION(I,:)-R_SH_ION(J,:),&
             R_SH_ION(I,:)-R_SH_ION(J,:))

        if(RRR2.GT.16.0_r8_kind) then
           ic=1
        else
           ic=common_atom(i,j)+1                     !!!!!!!!!!!!!AS
           if(host%ro(K,epe(J)%k,ic)==0.0_r8_kind) ic=1      !!!!!!!!!!!!AS
        endif

        IF(RRR2.GT.host%sr1(K,epe(J)%k,ic)**2) cycle
        EX=EX-host%C(K,epe(J)%k,ic)*(one/(RSS2**3)-one/(RRR2**3))-  &
              host%D(K,epe(J)%k,ic)*(one/(RSS2**4)-one/(RRR2**4))
        IF(RRR2.GT.host%sr2(K,epe(J)%k,ic)**2) cycle
        EY=EY+host%B(K,epe(J)%k,ic)*(EXP(-SQRT(RSS2)/host%RO(K,epe(J)%k,ic)) &
        -EXP(-SQRT(RRR2)/host%RO(K,epe(J)%k,ic)))
      enddo ! **region 1  - end
  enddo

! **summarizing different parts of energy
  E1=EX+EY
  E2=EA+EB
  E3=-0.5*(EC+ED)
  energy_shortrange=E1+E2+E3
  DPRINT 'energy_shortrange', energy_shortrange
  if(l_print) then
     write(output_epe,104)energy_shortrange,E1,E2,E3
104  FORMAT(1x,'short range interaction energy     shr=',F13.4, &
          /,1x,'                 of         area    1 =',F13.4, &
          /,1x,'                 interaction of 1 & 2A=',F13.4, &
          /,1x,'                 of         area    2a=',F13.4)
  endif
END function energy_shortrange
!-------------------------------------------------------------------

!--------------------------------------------------------
real(kind=r8_kind) function energy_3_body_po()
! **calc. of the 3-body interaction part of the relaxation energy of
! **lattice - optimization of a unit cell

  use epecom_module

  real(kind=r8_kind) :: e3b
  real(kind=r8_kind) :: deg2rad
  real(kind=r8_kind) :: rss1,rss2,scal1
  real(kind=r8_kind) :: theta
  integer(kind=i4_kind) :: i,j,k,i1,j1,k1
  integer(kind=i4_kind) :: ia_1,ia_2,ia_3,io,icl

  deg2rad=pi/180.0_r8_kind
  e3b=zero
  fst: do io=1,n_ions_cell
     icl=which_epe_ion(io)%new
     do i=1,n_tetrahedrons
        if(icl==tetra_atoms(i,1)) goto 1
     enddo
     cycle fst
1    ia_1=tetra_atoms(i,1)
     i1=epe(ia_1)%k
     scnd: do j=2,4
        ia_2=tetra_atoms(i,j)
        if(ia_2==0_i4_kind) cycle scnd
        j1=epe(ia_2)%k
        thrd: do k=j+1,5
           ia_3=tetra_atoms(i,k)
           if(ia_3==0_i4_kind) cycle thrd
           k1=epe(ia_3)%k
           rss1=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:))
           rss2=dot_product(r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:), &
                r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
           scal1=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
           rss1=sqrt(rss1)
           rss2=sqrt(rss2)
           theta=acos(scal1/(rss1*rss2))
           e3b=e3b+0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2
        enddo thrd!k=j+1,5
     enddo scnd!j=2,4
  enddo fst!io=1,n_ions_cell

  energy_3_body_po=e3b
end function energy_3_body_po
!-------------------------------------------------------

!------------------------------------------------------
real(kind=r8_kind) function energy_3_body()
! **calc. of the 3-body interaction part of the relaxation energy of
! **lattice

  use epecom_module

  real(kind=r8_kind) :: e3b
  real(kind=r8_kind) :: deg2rad
  real(kind=r8_kind) :: rss1,rss2,rsr1,rsr2,scal1,scal2
  real(kind=r8_kind) :: theta_a,theta_b
  integer(kind=i4_kind) :: i,j,k,i1,j1,k1,l
  integer(kind=i4_kind) :: ia_1,ia_2,ia_3

  deg2rad=pi/180.0_r8_kind
  e3b=zero
  fst:  do i=1,n_tetrahedrons
     do l=1,5
        if (tetra_atoms(i,l)<=reg_I_n_ions) goto 1
     enddo
     cycle fst
1    ia_1=tetra_atoms(i,1)
     i1=epe(ia_1)%k
     scnd: do j=2,4
        ia_2=tetra_atoms(i,j)
        if(ia_2==0) cycle scnd
        j1=epe(ia_2)%k
        thrd: do k=j+1,5
           ia_3=tetra_atoms(i,k)
           if(ia_3==0) cycle thrd
           if(ia_1 > reg_I_n_ions.and.ia_2 > reg_I_n_ions.and.ia_3 > reg_I_n_ions) cycle thrd
           k1=epe(ia_3)%k
           rss1=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:))
           rss2=dot_product(r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:), &
                r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
           rsr1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
                epe(ia_2)%r(:)-epe(ia_1)%r(:))
           rsr2=dot_product(epe(ia_3)%r(:)-epe(ia_1)%r(:), &
                epe(ia_3)%r(:)-epe(ia_1)%r(:))
           scal1=dot_product(r_sh_ion(ia_2,:)-r_sh_ion(ia_1,:), &
                r_sh_ion(ia_3,:)-r_sh_ion(ia_1,:))
           scal2=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
                epe(ia_3)%r(:)-epe(ia_1)%r(:))
           rss1=sqrt(rss1)
           rss2=sqrt(rss2)
           rsr1=sqrt(rsr1)
           rsr2=sqrt(rsr2)
           theta_a=acos(scal1/(rss1*rss2))
           theta_b=acos(scal2/(rsr1*rsr2))
           e3b=e3b+0.5*ki(j1,i1,k1)*((theta_a-theta_0(j1,i1,k1)*deg2rad)**2- &
                (theta_b-theta_0(j1,i1,k1)*deg2rad)**2)
        enddo thrd!k=j+1,5
     enddo scnd!j=2,4
  enddo fst!i=1,n_tetrahedrons

  energy_3_body=e3b
  DPRINT 'energy_3_body', e3b
end function energy_3_body
!------------------------------------------------------------

!------------------------------------------------------------
subroutine energy_vac_imp(use_ref_data,N_VACANCIES,N_IMPURITIES)
  ! **  energies due to defect - defect and defect - lattice interactions

  use epecom_module, n_vac=>N_VACANCIES,N_IMP=>N_IMPURITIES

  integer(kind=i4_kind), intent(in) ::N_VACANCIES,N_IMPURITIES
  logical, intent(in):: use_ref_data

  real(kind=r8_kind), allocatable, dimension(:)   :: reg_2a_ec
  real(kind=r8_kind) :: vs_ewaback,vc_ewaback
  real(kind=r8_kind) :: ecoul_vaccluster_epe,ec_cluster_epe
  real(kind=r8_kind) :: ec_reg2a,ec_reg2a_ref,ec_current,eshort_epecluster_epe
  real(kind=r8_kind) :: angle_fac
  real(kind=r8_kind) :: ec_ewadir_cluster, eshort_vacepe, &
        e_vac_3b=0.0_r8_kind,e_imp_3b=0.0_r8_kind
  integer(kind=i4_kind) :: i,ind,j,k,ii1,jj1
  integer(kind=i4_kind) :: status
  type epe_dot_products
     real(kind=r8_kind)::ss,sq_ss,sq_sc,sq_cs,sq_cc,ss3,ss6,ss8,ss10 &
          ,cs,cc,sc,sr,cr,rr,rc,rs_sr,rc_cr,ss_sr,sc_cr,cc_cr,cs_sr &
          ,rs
  end type epe_dot_products
  type factors
     real(kind=r8_kind)::ss,sc,cs,cc,sewa,cewa,gsewa,gcewa, &
          rs,rc
  end type factors
  real(kind=r8_kind),parameter::third=0.33333333_r8_kind,tenth=0.1_r8_kind,one=1.0_r8_kind
  real(kind=r8_kind)::serfunc_arg(3),cerfunc_arg(3),gewa_const
  type(epe_dot_products),target::reg,ref,var
  type(epe_dot_products),pointer::l
  type each_center_contributions
     real(kind=r8_kind) coul,short,regI,reg2a,eshort
     type(factors):: fac
  end type each_center_contributions
  type(each_center_contributions),allocatable,dimension(:)::clus,vaca
  real(kind=r8_kind), parameter::tail_exp_dist=7.0_r8_kind, &
       erfarg_th=0.035_r8_kind,gerfarg_th=0.026_r8_kind
  type epe_imp_tet
     real(kind=r8_kind) :: coord(5,3)
     character(len=3) :: center(5)
     integer(kind=i4_kind) :: number(5)
  end type epe_imp_tet

  epe(:)%s(1)=R_sh_ION(:,1)
  epe(:)%s(2)=R_sh_ION(:,2)
  epe(:)%s(3)=R_sh_ION(:,3)
  epe(:)%c(1)=R_nuc_ION(:,1)
  epe(:)%c(2)=R_nuc_ION(:,2)
  epe(:)%c(3)=R_nuc_ION(:,3)

  defect_energy_short = zero
  defect_energy_coul  = zero

  eshort_vaccluster=zero
  eshort_vacepe=zero
  eshort_epecluster=zero
  eshort_epecluster_epe=zero

  ecoul_vaccluster=zero
  ecoul_vaccluster_epe=zero
  ecoul_epecluster=zero
  ec_cluster_epe=zero

  gewa_const=4.0_r8_kind/3.0_r8_kind*ERROR_FUNCTION_PARAMETER**3/PIS

! **vacancies start monitoring
  allocate(vaca(N_VACANCIES), stat=status)
  if(status.ne.0 ) call error_handler("allocate vaca failed")
  vaca(1:n_vacancies)%coul  = zero ! ?
  vaca(1:n_vacancies)%short = zero

  vac: do I=1,N_VACANCIES
     k=epe(I)%k
     vaca(i)%reg2a=zero

!!$   vaca(i)%regi=-madc%emd(epe(i)%m)
     ! for vacancy in the shifted position one need recalculate
     ! Madelung energy contribution


     vaca(i)%regi=zero
     !*** direct space ewald contributions
     do j=1,reg_2a_n_ions

        var%sr=dot_product(epe(i)%s-epe(j)%r,epe(i)%s-epe(j)%r)
        var%cr=dot_product(epe(i)%c-epe(j)%r,epe(i)%c-epe(j)%r)

        serfunc_arg(1)=error_function_parameter*sqrt(var%sr)
        serfunc_arg(2)=serfunc_arg(1)**2
        serfunc_arg(3)=serfunc_arg(2)**2
        cerfunc_arg(1)=error_function_parameter*sqrt(var%cr)
        cerfunc_arg(2)=cerfunc_arg(1)**2
        cerfunc_arg(3)=cerfunc_arg(2)**2

        vaca(i)%fac%sewa=erfo*(one-third*serfunc_arg(2)+tenth*serfunc_arg(3))
        vaca(i)%fac%cewa=erfo*(one-third*cerfunc_arg(2)+tenth*cerfunc_arg(3))
        if(serfunc_arg(2).gt.erfarg_th) vaca(i)%fac%sewa=-erf(serfunc_arg(1))/sqrt(var%sr)
        if(cerfunc_arg(2).gt.erfarg_th) vaca(i)%fac%cewa=-erf(cerfunc_arg(1))/sqrt(var%cr)

        vaca(i)%regI=vaca(i)%regI - q_ion(epe(j)%k)* &
             (q_shell(epe(i)%k)*vaca(i)%fac%sewa+ &
             q_nuclear(epe(i)%k)*vaca(i)%fac%cewa)
     end do
     !*** back space ewald contributions
     vs_ewaback=zero
     do ind=1,n_bs_points
        angle_fac=dot_product(gstr(ind,:),epe(i)%s)*pi2
        vs_ewaback=vs_ewaback+rsin(ind)*sin(angle_fac)+rcos(ind)*cos(angle_fac)
     enddo!i=1,n_bs_points
     vaca(i)%regi=vaca(i)%regi-q_shell(epe(i)%k)*vs_ewaback

     vc_ewaback=zero
     do ind=1,n_bs_points
        angle_fac=dot_product(gstr(ind,:),epe(i)%c)*pi2
        vc_ewaback=vc_ewaback+rsin(ind)*sin(angle_fac)+rcos(ind)*cos(angle_fac)
     enddo!i=1,n_bs_points
     vaca(i)%regi=vaca(i)%regi-q_nuclear(epe(i)%k)*vc_ewaback

     ! shell model intra ionic energy contrib (4)
     vaca(i)%regi=vaca(i)%regi - pk(epe(i)%k)* &
          dot_product(epe(i)%c-epe(i)%s,epe(i)%c-epe(i)%s)/2.0_r8_kind

     call vac_to_vac_contibs

     ! **vacancy - lattice
     ! **treat contibutions from region I

     call vac_to_latt_contribs
     vaca(i)%coul=vaca(i)%reg2a+vaca(i)%regI
     eshort_vacepe=eshort_vacepe+vaca(i)%short
     ecoul_vaccluster_epe=ecoul_vaccluster_epe+vaca(i)%coul
  enddo vac

  if(n_types_central_atoms_3body > 0) then
     call vac_3_body(e_vac_3b)
     eshort_vacepe=eshort_vacepe+e_vac_3b  !!  ??
  end if

  deallocate(vaca,stat=status)
  if(status.ne.0)  call error_handler(" deallocate vaca failed")

  ! **start treatment of impurities
  diffpg_ec_ecref=zero
  diff_reg2a_ec_ecref=zero

  if(n_impurities /= 0) then
     if(.not.allocated(clus)) then
        allocate(clus(n_impurities),stat=status)
        if(status.ne.0)   call error_handler(" coul allocate clus failed")
     endif
     if(.not.allocated(reg_2a_ec)) then
        allocate(reg_2a_ec(n_impurities),stat=status)
        if(status.ne.0)   call error_handler("coul allocate reg_2a_ec failed")
     endif
     clus(1:n_impurities)%coul  = zero
     clus(1:n_impurities)%short = zero
     reg_2a_ec(:)               = zero
  end if

  impi: do I=1,N_IMPURITIES
     K=TYPE_IMPURITY(I)
     impj: do J=I+1,N_IMPURITIES
        ! **coulomb interaction of imurities : energies and gradients
        if(use_ref_data) then !only regular contribution (a constant)
           ii1=imp2center(i); jj1=imp2center(j)
           reg%ss=dot_product(reg_reference(i)%rs-reg_reference(j)%rs, &
                reg_reference(i)%rs-reg_reference(j)%rs)
           reg%cs=dot_product(reg_reference(i)%rc-reg_reference(j)%rs, &
                reg_reference(i)%rc-reg_reference(j)%rs)
           reg%sc=dot_product(reg_reference(i)%rs-reg_reference(j)%rc, &
                reg_reference(i)%rs-reg_reference(j)%rc)
           reg%cc=dot_product(reg_reference(i)%rc-reg_reference(j)%rc, &
                reg_reference(i)%rc-reg_reference(j)%rc)
           reg%rr=dot_product(R_IMP(I,:)-R_IMP(J,:),R_IMP(I,:)-R_IMP(J,:))
           clus(i)%coul=Q_SH_IMPURITY(I)*Q_SH_IMPURITY(J)/sqrt(reg%ss) &
                +Q_SH_IMPURITY(I)*Q_NUC_IMPURITY(J)/sqrt(reg%sc) &
                +Q_NUC_IMPURITY(I)*Q_SH_IMPURITY(J)/sqrt(reg%cs) &
                +Q_NUC_IMPURITY(I)*Q_NUC_IMPURITY(J)/sqrt(reg%cc)
           l=>reg
           CALL calc_shortrange &
                   (TYPE_IMPURITY(I),TYPE_IMPURITY(J),ii1,jj1, &
                   clus(i)%short)

        else    ! calculate energy of incluster interactions
           ii1=imp2center(i); jj1=imp2center(j)
           var%ss=dot_product(R_SH_IMP(I,:) &
                -R_SH_IMP(J,:),R_SH_IMP(I,:)-R_SH_IMP(J,:))
           var%cs=dot_product(R_NUC_IMP(I,:) &
                -R_SH_IMP(J,:),R_NUC_IMP(I,:)-R_SH_IMP(J,:))
           var%sc=dot_product(R_SH_IMP(I,:) &
                -R_NUC_IMP(J,:),R_SH_IMP(I,:)-R_NUC_IMP(J,:))
           var%cc=dot_product(R_NUC_IMP(I,:) &
                -R_NUC_IMP(J,:),R_NUC_IMP(I,:)-R_NUC_IMP(J,:))
           var%rr=dot_product(R_IMP(I,:)-R_IMP(J,:),R_IMP(I,:)-R_IMP(J,:))
           ! **short range interaction : energies and gradients

           l=>var
           CALL calc_shortrange &
                (TYPE_IMPURITY(I),TYPE_IMPURITY(J),ii1,jj1,clus(i)%short)

           clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SH_IMPURITY(J)/SQRT(var%ss)
           clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SH_IMPURITY(J)/SQRT(var%cs)
           clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUC_IMPURITY(J)/SQRT(var%sc)
           clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUC_IMPURITY(J)/SQRT(var%cc)

           ! calc  energy contrib (1)
           clus(i)%coul=clus(i)%fac%ss+clus(i)%fac%sc+clus(i)%fac%cs+clus(i)%fac%cc

        endif! use_ref_data/else

        eshort_epecluster=eshort_epecluster+clus(i)%short
        ecoul_epecluster=ecoul_epecluster+clus(i)%coul

     end do impj

     ! **done impurity - impurity iteraction

     ! **impurity - lattice interaction
     clus(i)%coul=zero
     clus(i)%short=zero

! **coulomb energy
! **a) calculate contribution to the coulomb energy
! **   due to the evald sumation in the coordinate space
     ec_ewadir_cluster=zero
     vacj: DO J=1,N_VACANCIES
        if(use_ref_data) then ! only regular data
           reg%sr=dot_product(reg_reference(I)%rs-epe(J)%r,reg_reference(I)%rs-epe(J)%r)
           reg%cr=dot_product(reg_reference(I)%rc-epe(J)%r,reg_reference(I)%rc-epe(J)%r)
        else
           reg%sr=dot_product(R_SH_IMP(I,:)-epe(J)%r,R_SH_IMP(I,:)-epe(J)%r )
           reg%cr=dot_product(R_NUC_IMP(I,:)-epe(J)%r,R_NUC_IMP(I,:)-epe(J)%r)
        endif

        serfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(reg%sr)
        serfunc_arg(2)=serfunc_arg(1)**2
        serfunc_arg(3)=serfunc_arg(2)**2
        cerfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(reg%cr)
        cerfunc_arg(2)=cerfunc_arg(1)**2
        cerfunc_arg(3)=cerfunc_arg(2)**2
        clus(i)%fac%sewa=ERFO*(one-third*serfunc_arg(2)+tenth*serfunc_arg(3))
        clus(i)%fac%cewa=ERFO*(one-third*cerfunc_arg(2)+tenth*cerfunc_arg(3))
        IF(serfunc_arg(2).GT.erfarg_th) clus(i)%fac%sewa=-ERF(serfunc_arg(1))/SQRT(reg%sr)
        if(cerfunc_arg(2).gt.erfarg_th) clus(i)%fac%cewa=-erf(cerfunc_arg(1))/sqrt(reg%cr)

        !calc ewald dir summation energy contrib (2)
        ec_ewadir_cluster=ec_ewadir_cluster+ Q_ion(epe(J)%k)*( &
             Q_SH_IMPURITY(I)*clus(i)%fac%sewa+ &
             Q_NUC_IMPURITY(I)*clus(i)%fac%cewa)
     enddo vacj!J=1,N_VACANCIES

     clus(i)%coul=clus(i)%coul+ec_ewadir_cluster

     ! **Evald method reciprocial space contribution to energy (clus(i)%coul)
     ! **acctually  clus(i)%coul contribution for use_ref_data.eq.true have to be calculated
     ! **quantumchemically
     !** treat energy contributions
     if(use_ref_data) then
        vs_ewaback=zero
        DO ind=1,n_bs_points
           angle_fac=dot_product(GSTR(ind,:),reg_reference(i)%rs)*PI2
           vs_ewaback=vs_ewaback+RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
        enddo
        clus(i)%coul=clus(i)%coul+Q_SH_IMPURITY(I)*vs_ewaback
        vc_ewaback=zero
        DO ind=1,n_bs_points
           angle_fac=dot_product(GSTR(ind,:),reg_reference(i)%rc)*PI2
           vc_ewaback=vc_ewaback+RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
        enddo!I=1,n_bs_points
        clus(i)%coul=clus(i)%coul+Q_NUC_IMPURITY(I)*vc_ewaback
        clus(i)%coul=clus(i)%coul + PK_IMPURITY(I)* &
             dot_product(reg_reference(i)%rc-reg_reference(i)%rs, &
             reg_reference(i)%rc-reg_reference(i)%rs)/2.0_r8_kind

     else ! regular run

        ! ewald backspace energy contrib (3)
        vs_ewaback=zero
        DO ind=1,n_bs_points
           angle_fac=dot_product(GSTR(ind,:),R_SH_IMP(i,:))*PI2
           vs_ewaback=vs_ewaback+RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
        enddo!I=1,n_bs_points
        clus(i)%coul=clus(i)%coul+Q_SH_IMPURITY(I)*vs_ewaback

        vc_ewaback=zero
        DO ind=1,n_bs_points
           angle_fac=dot_product(GSTR(ind,:),R_nuc_IMP(i,:))*PI2
           vc_ewaback=vc_ewaback+RSIN(ind)*SIN(angle_fac)+RCOS(ind)*COS(angle_fac)
        enddo!I=1,n_bs_points
        clus(i)%coul=clus(i)%coul+Q_NUC_IMPURITY(I)*vc_ewaback

        ! shell model intra ionic energy contrib (4)
        clus(i)%coul=clus(i)%coul + pk_impurity(i)* & !!!!!
             dot_product(r_nuc_imp(i,:)-r_sh_imp(i,:), &
             r_nuc_imp(i,:)-r_sh_imp(i,:))/2.0_r8_kind
     endif
     ! **done Ewald method reciprocial space contributions

     ! **a)calculate contributions to energy (clus(i)%coul) due to interaction
     ! **  of impurities with shell model ions of region A and
     ! **  corresponding contributions of Evald summation in the
     ! **  coordinate space
     ! for use_ref_data mode add regular and substract reference contribs
     reg_Ii: DO J=n_vacancies+1,reg_I_n_ions
        !calc dot products
        var%sr=dot_product(R_SH_IMP(I,:)-epe(J)%r,R_SH_IMP(I,:)-epe(J)%r)
        var%cr=dot_product(R_NUC_IMP(I,:)-epe(J)%r,R_NUC_IMP(I,:)-epe(J)%r)
        var%ss=dot_product(R_SH_IMP(I,:)-epe(J)%s,R_SH_IMP(I,:)-epe(J)%s)
        var%cc=dot_product(R_NUC_IMP(I,:)-epe(J)%c,R_NUC_IMP(I,:)-epe(J)%c)
        var%sc=dot_product(R_SH_IMP(I,:)-epe(J)%c,R_SH_IMP(I,:)-epe(J)%c)
        var%cs=dot_product(R_NUC_IMP(I,:)-epe(J)%s,R_NUC_IMP(I,:)-epe(J)%s)
        var%rr=dot_product(R_IMP(I,:)-epe(J)%r,R_IMP(I,:)-epe(J)%r)

        ! error_function contribs
        serfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(var%sr)
        serfunc_arg(2)=serfunc_arg(1)**2
        cerfunc_arg(1)=ERROR_FUNCTION_PARAMETER*SQRT(var%cr)
        cerfunc_arg(2)=cerfunc_arg(1)**2
        clus(i)%fac%sewa=ERF(serfunc_arg(1))/SQRT(var%sr)
        clus(i)%fac%cewa=ERF(cerfunc_arg(1))/SQRT(var%cr)

        ! energy prefactors
        clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%cc)
        clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%cs)
        clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%sc)
        clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%ss)

        if(use_ref_data) then
           reg%ss=dot_product(reg_reference(i)%rs-epe(j)%s, &
                reg_reference(i)%rs- epe(j)%s)
           reg%sc=sqrt(dot_product(reg_reference(i)%rs-epe(j)%c, &
                reg_reference(i)%rs-epe(j)%c ))
           reg%sr=sqrt(dot_product(reg_reference(i)%rs-epe(J)%r,  &
                reg_reference(i)%rs-epe(J)%r  ))
           reg%cr=sqrt(dot_product(reg_reference(i)%rc-epe(J)%r,  &
                reg_reference(i)%rc-epe(J)%r  ))
           reg%cs=sqrt(dot_product(reg_reference(i)%rc-epe(j)%s, &
                reg_reference(i)%rc-epe(j)%s ))
           reg%cc=sqrt(dot_product(reg_reference(i)%rc-epe(j)%c, &
                reg_reference(i)%rc-epe(j)%c ))
           reg%rr=var%rr
           clus(i)%fac%sewa= ERF(ERROR_FUNCTION_PARAMETER*reg%sr)/reg%sr
           clus(i)%fac%cewa= ERF(ERROR_FUNCTION_PARAMETER*reg%cr)/reg%cr
           !*** coulon contrib regular data
           clus(i)%coul=clus(i)%coul+Q_SH_IMPURITY(I)*( &
                Q_SHELL(epe(j)%k)/sqrt(reg%ss)+ Q_NUCLEAR(epe(j)%k)/reg%sc &
                -Q_ion(epe(J)%k)*clus(i)%fac%sewa)&
                +Q_NUC_IMPURITY(I)*( &
                Q_SHELL(epe(j)%k)/reg%cs+ Q_NUCLEAR(epe(j)%k)/reg%cc- &
                Q_ion(epe(J)%k)*clus(i)%fac%cewa)

           if(i.eq.n_impurities.and.use_epe_pgdata) diffpg_ec_ecref=diffpg_ec_ecref &
                -(reg_I_pg(j)%vs*Q_SHELL(epe(j)%k)+ &
                reg_I_pg(j)%vc*Q_NUCLEAR(epe(j)%k))*7.17115929581458_r8_kind
           ! add (PG_var-PG_ref) contrib

        else ! regular run not use_ref_data mode

           ! qq/r + ewald err function contrib
           clus(i)%coul=clus(i)%coul+ &
                clus(i)%fac%ss+clus(i)%fac%sc &
                -Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa &
                + clus(i)%fac%cc+clus(i)%fac%cs &
                -Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa
        endif!  use_ref_data/else

        ! **calculate short-range interaction contributions to
        ! **energy (clus(i)%short)
        ! **short-range interaction impurity-lattice
        l=>var
        if(use_ref_data.and.explicit_coupling) then
           ii1=imp2center(i)
           CALL calc_shortrange(explicit_coupling_type(i),epe(j)%k,ii1,j, &
                clus(i)%eshort,harmonic_bondfix=.true.)
        else
           ii1=imp2center(i)
           CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,ii1,j,clus(i)%eshort) !var
        endif
        clus(i)%short=clus(i)%short+clus(i)%eshort ! No 1

        if(use_ref_data) then
           ii1=imp2center(i)
           l=>reg
           CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,ii1,j,clus(i)%eshort)
           clus(i)%short=clus(i)%short+clus(i)%eshort

           ref%ss=dot_product(epe_reference(i)%rs-R_SH_ION(J,:), &
                epe_reference(i)%rs-R_SH_ION(J,:) )
           ref%rr=var%rr
           l=>ref
           if(explicit_coupling) then
              ii1=imp2center(i)
              CALL calc_shortrange(explicit_coupling_type(i),epe(j)%k,ii1,j, &
                   clus(i)%eshort)
           else
              ii1=imp2center(i)
              CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,ii1,j,clus(i)%eshort)
           end if
           clus(i)%short=clus(i)%short-clus(i)%eshort
        else
        endif
     enddo reg_Ii

     ! ** 2a-region start
     ec_reg2a=zero
     reg_2aj: DO J=reg_I_n_ions+1,reg_2a_n_ions

        ! calc dot products
        var%rr=dot_product(r_imp(i,:)-epe(j)%r,r_imp(i,:)-epe(j)%r)
        var%sr=sqrt(dot_product(R_SH_IMP(I,:)-epe(j)%r,R_SH_IMP(I,:)-epe(j)%r))
        var%cr=sqrt(dot_product(R_NUC_IMP(I,:)-epe(j)%r,R_NUC_IMP(I,:)-epe(j)%r))
        var%ss=dot_product(R_SH_IMP(I,:)-epe(j)%s,R_SH_IMP(I,:)-epe(j)%s)
        var%cc=dot_product(R_NUC_IMP(I,:)-epe(j)%c,R_NUC_IMP(I,:)-epe(j)%c)
        var%sc=dot_product(R_SH_IMP(I,:)-epe(j)%c,R_SH_IMP(I,:)-epe(j)%c)
        var%cs=dot_product(R_NUC_IMP(I,:)-epe(j)%s,R_NUC_IMP(I,:)-epe(j)%s)
        var%ss_sr=dot_product(R_SH_IMP(I,:)-epe(j)%s,epe(j)%s-epe(j)%r)
        var%sc_cr=dot_product(R_SH_IMP(I,:)-epe(j)%c,epe(j)%c-epe(j)%r)
        var%cc_cr=dot_product(R_NUC_IMP(I,:)-epe(j)%c,epe(j)%c-epe(j)%r)
        var%cs_sr=dot_product(R_NUC_IMP(I,:)-epe(j)%s,epe(j)%s-epe(j)%r)

        ! var shortrange contribs
        l=>var
        if(use_ref_data.and.explicit_coupling) then
           ii1=imp2center(i)
           CALL calc_shortrange(explicit_coupling_type(I),epe(j)%k,ii1,j, &
                clus(i)%eshort)
        else
           ii1=imp2center(i)
           CALL calc_shortrange(TYPE_IMPURITY(I),epe(j)%k,ii1,j,clus(i)%eshort)
        end if
        clus(i)%short=clus(i)%short+clus(i)%eshort  ! No 2

        ! var longrange prefactors
        clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%cc)
        clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%cs)
        clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(var%sc)
        clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(var%ss)
        clus(i)%fac%sewa=ERF(ERROR_FUNCTION_PARAMETER*var%sr)/var%sr
        clus(i)%fac%cewa=ERF(ERROR_FUNCTION_PARAMETER*var%cr)/var%cr

        !var  qq/r + ewald error function contrib
        ! the same contribution is calculated with PG
        ! when cluster is embedded in epe PC
        ec_current= &
             clus(i)%fac%ss+clus(i)%fac%sc + &
             clus(i)%fac%cc+clus(i)%fac%cs  &
             -Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa &
             -Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa

        ec_reg2a=ec_reg2a+ec_current

        if(use_ref_data) then ! modify code for energy only
           reg%ss=dot_product(reg_reference(i)%rs-epe(J)%s, &
                reg_reference(i)%rs-epe(J)%s)
           reg%sr=sqrt(dot_product(reg_reference(i)%rs-epe(j)%r, &
                reg_reference(i)%rs-epe(j)%r))
           reg%cr=sqrt(dot_product(reg_reference(i)%rc-epe(j)%r, &
                reg_reference(i)%rc-epe(j)%r))
           reg%cc=dot_product(reg_reference(i)%rc-epe(J)%c, &
                reg_reference(i)%rc-epe(J)%c)
           reg%sc=dot_product(reg_reference(i)%rs-epe(J)%c, &
                reg_reference(i)%rs-epe(J)%c)
           reg%cs=dot_product(reg_reference(i)%rc-epe(J)%s, &
                reg_reference(i)%rc-epe(J)%s)
           reg%rr=dot_product(r_imp(i,:)-epe(j)%r,r_imp(i,:)-epe(j)%r)

           clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(reg%cc)
           clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(reg%cs)
           clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(reg%sc)
           clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(reg%ss)
           clus(i)%fac%sewa=ERF(ERROR_FUNCTION_PARAMETER*reg%sr)/reg%sr
           clus(i)%fac%cewa=ERF(ERROR_FUNCTION_PARAMETER*reg%cr)/reg%cr

           ! var longrange reg contribs for use_ref_data mode
           ec_reg2a=ec_reg2a+clus(i)%fac%ss+clus(i)%fac%sc- &
                Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa &
                +clus(i)%fac%cc+clus(i)%fac%cs- &
                Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa

           ! reg shortrange contribs for use_ref_data mode
           reg%ss_sr=dot_product(reg_reference(i)%rs-epe(j)%s,epe(j)%s-epe(j)%r)
           l=>reg
           ii1=imp2center(i)
           CALL calc_shortrange(TYPE_IMPURITY(I),epe(J)%k,ii1,j,clus(i)%eshort)
           clus(i)%short=clus(i)%short+clus(i)%eshort

           ref%ss=dot_product(epe_reference(i)%rs-epe(j)%s,epe_reference(i)%rs-epe(j)%s)
           ref%sr=sqrt(dot_product(epe_reference(i)%rs-epe(j)%r,epe_reference(i)%rs-epe(j)%r))
           ref%cr=sqrt(dot_product(epe_reference(i)%rc-epe(j)%r,epe_reference(i)%rc-epe(j)%r))
           ref%cc=dot_product(epe_reference(i)%rc-epe(j)%c,epe_reference(i)%rc-epe(j)%c)
           ref%sc=dot_product(epe_reference(i)%rs-epe(j)%c,epe_reference(i)%rs-epe(j)%c)
           ref%cs=dot_product(epe_reference(i)%rc-epe(j)%s,epe_reference(i)%rc-epe(j)%s)
           ref%rr=var%rr

           clus(i)%fac%cc=Q_NUC_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(ref%cc)
           clus(i)%fac%cs=Q_NUC_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(ref%cs)
           clus(i)%fac%sc=Q_SH_IMPURITY(I)*Q_NUCLEAR(epe(j)%k)/sqrt(ref%sc)
           clus(i)%fac%ss=Q_SH_IMPURITY(I)*Q_SHELL(epe(j)%k)/sqrt(ref%ss)
           clus(i)%fac%sewa=ERF(ERROR_FUNCTION_PARAMETER*ref%sr)/ref%sr
           clus(i)%fac%cewa=ERF(ERROR_FUNCTION_PARAMETER*ref%cr)/ref%cr

           ! ref lonrange contrib for use_ref_data mode
           ec_reg2a_ref=clus(i)%fac%ss+clus(i)%fac%sc- &
                Q_SH_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%sewa  &
                +clus(i)%fac%cs+clus(i)%fac%cc- &
                Q_NUC_IMPURITY(I)*Q_ion(epe(J)%k)*clus(i)%fac%cewa
           ec_reg2a=ec_reg2a-ec_reg2a_ref

           ! ref shortrange contribs for use_ref_data mode
           ref%ss_sr=dot_product(epe_reference(i)%rs-R_SH_ION(J,:),R_SH_ION(J,:)-epe(j)%r)
           l=>ref
           if(explicit_coupling) then
              ii1=imp2center(i)
              CALL calc_shortrange(explicit_coupling_type(I),epe(J)%k,ii1,j, &
                   clus(i)%eshort)
           else
              ii1=imp2center(i)
              CALL calc_shortrange(TYPE_IMPURITY(I),epe(J)%k,ii1,j,clus(i)%eshort)
           end if
           clus(i)%short=clus(i)%short-clus(i)%eshort
        endif
     enddo reg_2aj

!     clus(i)%tot=clus(i)%coul+ec_reg2a+clus(i)%short
     eshort_epecluster_epe=eshort_epecluster_epe+clus(i)%short
     ec_cluster_epe=ec_cluster_epe+ec_reg2a+clus(i)%coul
  enddo impi

  if(qm_interfaced_mode) then
     if(n_types_central_atoms_3body > 0) then
        call imp_3_body_1cul(e_imp_3b)
     endif
  else
     if(n_types_central_atoms_3body_im > 0) then
        call imp_3_body_cul(e_imp_3b)
     endif
  endif

  DPRINT 'e_imp_3b in culon', e_imp_3b

  if(n_impurities /= 0) then
     deallocate(clus,stat=status)
     if(status.ne.0) call error_handler("deallocate clus failed")
  end if
  ec_cluster_epe=ec_cluster_epe+diffpg_ec_ecref
  ! **region 2a finished
  ! **impurities done

  defect_energy_short=eshort_vaccluster+eshort_vacepe+ &
       eshort_epecluster+eshort_epecluster_epe+e_imp_3b
  defect_energy_coul=ec_cluster_epe &
       +ecoul_vaccluster_epe+ecoul_vaccluster+ecoul_epecluster+E2BDEF

  DPRINT 'defect_energies vac imp', defect_energy_short,defect_energy_coul

  !..............................................
contains

  SUBROUTINE calc_shortrange(K,N,ik,jn,EP,harmonic_bondfix)
    ! **potential of interaction of ions with a lattice and
    ! **each other.

    use epecom_module, only: host

    logical, optional, intent(in):: harmonic_bondfix
    real(kind=r8_kind) :: EP
    integer(kind=i4_kind) :: K,N,M,ik,jn

    l%sq_ss=sqrt(l%ss)
    l%ss3=l%sq_ss*l%ss
    l%ss6=l%ss3**2
    l%ss8=l%ss6*l%ss
    l%ss10=l%ss8*l%ss

    ! **calculation of energy and gradients for remote neighbours
    if(l%ss .gt. 16.0_r8_kind) then
       M=1
    else
       M=common_atom(ik,jn)+1                 !!!!!!!!!!!!!!AS
       if(host%ro(K,N,M)==0.0_r8_kind) M=1   !!!!!!!!!!!!!!AS
    end if
    if(l%ss.le.host%sr1(K,N,M)**2) then
       EP =-host%C(K,N,M)/(l%ss6)-host%D(K,N,M)/(l%ss8)
       IF(l%ss.le.host%sr2(K,N,M)**2) then
          EP=EP+host%B(K,N,M)*EXP(-l%sq_ss/host%RO(K,N,M))
       endif
    else
       EP = 0.0_r8_kind
    end if

          if(present(harmonic_bondfix)) then
           ! r in a.u , energy in eV
           if(host%k(K,N,m).ne.0.0.and.abs(sqrt(l%rr)-host%r0(K,N,m)).lt.0.1) then
           DPRINT 'culon harmonic_bond_fix',l%sq_ss
           EP=EP+host%k(K,N,m)*(eau_ev/auangs**2)*(l%sq_ss-host%r0(K,N,m))**2
           endif
          endif

  END SUBROUTINE calc_shortrange
  !...........................................

  !..........................................
  subroutine vac_to_vac_contibs

    vaca(i)%coul=zero
    do j=1,n_vacancies
       if(i.eq.j) cycle
       var%rr=dot_product(epe(I)%r-epe(J)%r,epe(I)%r-epe(J)%r)
       var%ss=dot_product(epe(I)%s-epe(J)%s,epe(I)%s-epe(J)%s)
       var%cs=dot_product(epe(I)%c-epe(J)%s,epe(I)%c-epe(J)%s)
       var%sc=dot_product(epe(I)%s-epe(J)%c,epe(I)%s-epe(J)%c)
       var%cc=dot_product(epe(I)%c-epe(J)%c,epe(I)%c-epe(J)%c)

       vaca(i)%fac%cc=q_nuclear(epe(i)%k)*q_nuclear(epe(j)%k)/sqrt(var%cc)
       vaca(i)%fac%cs=q_nuclear(epe(i)%k)*q_shell(epe(j)%k)/sqrt(var%cs)
       vaca(i)%fac%sc=q_shell(epe(i)%k)*q_nuclear(epe(j)%k)/sqrt(var%sc)
       vaca(i)%fac%ss=q_shell(epe(i)%k)*q_shell(epe(j)%k)/sqrt(var%ss)

       if(i.ge.j) cycle
       vaca(i)%coul=-(vaca(i)%fac%ss+vaca(i)%fac%sc+vaca(i)%fac%cs+vaca(i)%fac%cc)
       l=>var
       call calc_shortrange(epe(i)%k,epe(j)%k,i,j,vaca(i)%short)
       eshort_vaccluster=eshort_vaccluster-vaca(i)%short
       ecoul_vaccluster=ecoul_vaccluster+vaca(i)%coul
    enddo
  end subroutine vac_to_vac_contibs
  !..............................................

  !.............................................
  subroutine vac_to_latt_contribs

    vaca(i)%short=zero

    do j=n_vacancies+1,reg_i_n_ions

       ! dot products
       var%rr=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)
       var%ss=dot_product(epe(i)%s-epe(j)%s,epe(i)%s-epe(j)%s)
       var%cs=dot_product(epe(i)%c-epe(j)%s,epe(i)%c-epe(j)%s)
       var%sc=dot_product(epe(i)%s-epe(j)%c,epe(i)%s-epe(j)%c)
       var%cc=dot_product(epe(i)%c-epe(j)%c,epe(i)%c-epe(j)%c)

       l=>var
       call calc_shortrange(epe(i)%k,epe(j)%k,i,j,vaca(i)%eshort)
       vaca(i)%short=vaca(i)%short-vaca(i)%eshort
       var%rc=dot_product(epe(i)%r-epe(j)%c, epe(i)%r-epe(j)%c)

       var%sq_sc=sqrt(var%sc)
       var%sq_cs=sqrt(var%cs)
       var%sq_cc=sqrt(var%cc)

       vaca(i)%regI=vaca(i)%regI &
            - q_shell(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_ss+q_nuclear(epe(j)%k)/var%sq_sc) &
            - q_nuclear(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_cs+q_nuclear(epe(j)%k)/var%sq_cc)

    enddo
    ! ** area    i  finished

    ! **treat region 2a contributions
    DO J=reg_I_n_ions+1,reg_2a_n_ions

       var%rr=dot_product(epe(i)%r-epe(j)%r,epe(i)%r-epe(j)%r)

       var%ss=dot_product(epe(i)%s-epe(j)%s,epe(i)%s-epe(j)%s)
       var%cs=dot_product(epe(i)%c-epe(j)%s,epe(i)%c-epe(j)%s)
       var%sc=dot_product(epe(i)%s-epe(j)%c,epe(i)%s-epe(j)%c)
       var%cc=dot_product(epe(i)%c-epe(j)%c,epe(i)%c-epe(j)%c)
       var%ss_sr=dot_product(epe(i)%s-epe(j)%s,epe(j)%s-epe(j)%r)

       l=>var
       call calc_shortrange(epe(I)%k,epe(J)%k,i,j,vaca(i)%eshort)
       vaca(i)%short=vaca(i)%short-vaca(i)%eshort
       ! contrib No 1 vacancy

       var%rc=dot_product(epe(i)%r-epe(J)%c,epe(i)%r-epe(J)%c)
       var%rs_sr=dot_product(epe(i)%r-epe(J)%s,epe(j)%s-epe(J)%r)
       var%rc_cr=dot_product(epe(i)%r-epe(j)%c,epe(j)%c-epe(j)%r)

       ! **monopol contribution        core and shell j with the vacancy
       var%sq_cc=sqrt(var%cc)
       var%sq_sc=sqrt(var%sc)
       var%sq_cs=sqrt(var%cs)
       vaca(i)%reg2a=vaca(i)%reg2a &
!!$            + q_ion(epe(i)%k)*q_ion(epe(j)%k)/sqrt(var%rr) &
            - q_shell(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_ss+q_nuclear(epe(j)%k)/var%sq_sc) &
            - q_nuclear(epe(i)%k)*(q_shell(epe(j)%k)/var%sq_cs+q_nuclear(epe(j)%k)/var%sq_cc)

       ! **(5) dipol contribution, j-dipol - i-vacanvy
       ecoul_lat_relaxation=ecoul_lat_relaxation + 0.5*Q_ion(epe(I)%k) * &
            (Q_SHELL(epe(J)%k)*var%rs_sr/(var%ss3) &
            +Q_NUCLEAR(epe(J)%k)*var%rc_cr/(var%rc*SQRT(var%rc)))
    enddo
    ! **area 2a finished
  end subroutine vac_to_latt_contribs
  !.......................................

  !.......................................
  subroutine imp_3_body_1cul(e3bi)

    real(kind=r8_kind), intent(out) :: e3bi
    real(kind=r8_kind) :: theta,cos_th,rss1,rss2,scal1
    real(kind=r8_kind) :: r1(3),r2(3),r3(3)
    real(kind=r8_kind) :: deg2rad
    integer(kind=i4_kind) :: i,j,k,ia_1,ia_2,ia_3,i1,j1,k1,l,ii

    deg2rad=pi/180.0_r8_kind
    e3bi=zero
    if( n_vacancies == 0 ) return
    fst:  do i=1,n_tetrahedrons
       do l=1,5
          if (tetra_atoms(i,l) <= n_vacancies) goto 1
       enddo
       cycle fst
1      ia_1=tetra_atoms(i,1)
       i1=epe(ia_1)%k
       scnd: do j=2,4
          ia_2=tetra_atoms(i,j)
          if(ia_2 == 0) cycle scnd
          j1=epe(ia_2)%k
          thrd: do k=j+1,5
             ia_3=tetra_atoms(i,k)
             if(ia_3 == 0) cycle thrd
             k1=epe(ia_3)%k
             if(ia_1 > n_vacancies.and.ia_2 >n_vacancies .and.ia_3 > n_vacancies) cycle thrd
             if(ia_1 > n_vacancies) then
                r1=epe(ia_1)%s(:)
             else
                ii=epe(ia_1)%gc
                if(use_ref_data) then
                   r1=reg_reference(ii)%rs
                else
                   r1=r_sh_imp(ii,:)
                endif
             endif
             if(ia_2 > n_vacancies) then
                r2=epe(ia_2)%s(:)
             else
                ii=epe(ia_2)%gc
                if(use_ref_data) then
                   r2=reg_reference(ii)%rs
                else
                   r2=r_sh_imp(ii,:)
                endif
             endif
             if(ia_3 > n_vacancies) then
                r3=epe(ia_3)%s(:)
             else
                ii=epe(ia_3)%gc
                if(use_ref_data) then
                   r3=reg_reference(ii)%rs
                else
                   r3=r_sh_imp(ii,:)
                endif
             endif

             rss1=dot_product(r2-r1,r2-r1)
             rss2=dot_product(r3-r1,r3-r1)
             scal1=dot_product(r2-r1,r3-r1)

             rss1=sqrt(rss1)
             rss2=sqrt(rss2)
             cos_th=scal1/(rss1*rss2)
             theta=acos(cos_th)

             e3bi=e3bi+0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2

          enddo thrd!k=j+1,4
       enddo scnd!j=1,3
    enddo fst!i=1,n_tetra_atoms

  end subroutine imp_3_body_1cul
  !.......................................

  !.......................................
  subroutine imp_3_body_cul(e3bi)

    real(kind=r8_kind), intent(out) :: e3bi
    type(epe_imp_tet) :: tetra_imp

    integer(kind=i4_kind) :: i,j,k,i1,jj

    e3bi=zero
    if( n_impurities == 0 ) return
    fst: do i=1,n_types_central_atoms_3body_im
       i1=types_im(i,1)

       ! impurity-impurity + impurity-lattice
       jj=0
       scnd: do j=1,n_impurities

          if(i1 == type_impurity(j)) then
             call select_centers(tetra_imp,j,'imp',i)
             call energy(e3bi,tetra_imp)
          else
             jj=jj+1
          endif
       enddo scnd

       ! lattice-impurity
       if(jj == n_impurities) then
          thrd: do j=n_vacancies+1,reg_I_n_ions
             if(i1 ==  epe(j)%k) then
                call select_centers(tetra_imp,j,'epe',i)
                do k=2,5
                   if(tetra_imp%center(k) == 'imp') goto 1
                enddo
                cycle thrd
1               call energy(e3bi,tetra_imp)
             endif
          enddo thrd
       endif
    enddo fst

  end subroutine imp_3_body_cul
  !.......................................

  !.......................................
  subroutine energy(e3bi,tetra_imp)

    real(kind=r8_kind), intent(inout) :: e3bi
    type(epe_imp_tet),intent(in) :: tetra_imp

    real(kind=r8_kind) :: theta,cos_th,rss1,rss2,scal1
    real(kind=r8_kind) :: deg2rad
    integer(kind=i4_kind) :: ia_1,ia_2,ia_3,i1,j1,k1,j,k

    deg2rad=pi/180.0_r8_kind
    ia_1=tetra_imp%number(1)
    if(tetra_imp%center(1) == 'imp') then
       i1=type_impurity(ia_1)
    else
       i1=epe(ia_1)%k
    endif

    do j=2,4
       ia_2=tetra_imp%number(j)
       if(ia_2 == 0) cycle
       if(tetra_imp%center(j) == 'imp') then
          j1=type_impurity(ia_2)
       else
          j1=epe(ia_2)%k
       endif

       do k=j+1,5
          ia_3=tetra_imp%number(k)
          if(ia_3 == 0) cycle
          if(tetra_imp%center(k) == 'imp') then
             k1=type_impurity(ia_3)
          else
             k1=epe(ia_3)%k
          endif
          if(tetra_imp%center(1) /='imp'.and. tetra_imp%center(j) /='imp'.and. &
               tetra_imp%center(k) /= 'imp') cycle

          rss1=dot_product(tetra_imp%coord(j,:)-tetra_imp%coord(1,:), &
               tetra_imp%coord(j,:)-tetra_imp%coord(1,:))
          rss2=dot_product(tetra_imp%coord(k,:)-tetra_imp%coord(1,:), &
               tetra_imp%coord(k,:)-tetra_imp%coord(1,:))
          scal1=dot_product(tetra_imp%coord(j,:)-tetra_imp%coord(1,:), &
               tetra_imp%coord(k,:)-tetra_imp%coord(1,:))

          rss1=sqrt(rss1)
          rss2=sqrt(rss2)
          cos_th=scal1/(rss1*rss2)
          theta=acos(cos_th)

          e3bi=e3bi+0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2
       enddo
    enddo

  end subroutine energy
  !.......................................

  !.......................................
  subroutine select_centers(tetra_imp,num,who_am_i,n_type_3b)

    type(epe_imp_tet),intent(out) :: tetra_imp
    integer(kind=i4_kind),intent(in) :: num,n_type_3b
    character(len=3),intent(in) :: who_am_i

    real(kind=r8_kind) :: dist,rrr(3)
    real(kind=r8_kind), dimension(4) :: buf_dist
    integer(kind=i4_kind), dimension(4) :: buf_index
    character(len=3), dimension(4) :: buf_type
    integer(kind=i4_kind) :: max_ind(1)
    real(kind=r8_kind) :: max_val
    integer(kind=i4_kind) :: i,j,k

    tetra_imp%coord=0.0_r8_kind
    tetra_imp%number(1)=num
    if(who_am_i == 'imp') then
       tetra_imp%center(1)='imp'
       if(use_ref_data) then
          tetra_imp%coord(1,:)=reg_reference(num)%rs
          rrr=reg_reference(num)%rs
       else
          tetra_imp%coord(1,:)=r_sh_imp(num,:)
          rrr=r_sh_imp(num,:)
       endif
    else
       tetra_imp%center(1)='epe'
!       tetra_imp%coord(1,:)=epe(num)%s
!       rrr=epe(num)%s
       tetra_imp%coord(1,:)=epe(num)%r
       rrr=epe(num)%r
    endif

    buf_dist=0.0_r8_kind
    buf_index=0
    fst: do i=1,n_impurities
       do j=2,5
          if(type_impurity(i)==types_im(n_type_3b,j)) goto 1
       enddo
       cycle fst

1      if(use_ref_data) then
          dist=dot_product(reg_reference(i)%rs-rrr, &
               reg_reference(i)%rs-rrr)
!               reg_reference(j)%rs-rrr) ! bug fixed
       else
          dist=dot_product(r_sh_imp(i,:)-rrr,r_sh_imp(i,:)-rrr)
       endif
       dist=sqrt(dist)
       if(dist > r3b_im(n_type_3b)) cycle fst

       do k=1,4
          if(buf_dist(k) == 0.0_r8_kind) then
             buf_dist(k)=dist
             buf_index(k)=i
             buf_type(k)='imp'
             cycle fst
          endif
       enddo
       max_ind=maxloc(buf_dist)
       max_val=maxval(buf_dist)
       if(dist < max_val) then
          buf_dist(max_ind)=dist
          buf_index(max_ind)=i
          buf_type(max_ind)='imp'
       endif
    enddo fst

    scnd: do i=n_vacancies+1,reg_I_n_ions
       if(i == num) cycle scnd
       do j=2,5
          if(epe(i)%k == types_im(n_type_3b,j)) goto 2
       enddo
       cycle scnd

2      dist=dot_product(epe(i)%s-rrr,epe(i)%s-rrr)
       dist=sqrt(dist)
       if(dist > r3b_im(n_type_3b)) cycle scnd

       do k=1,4
          if(buf_dist(k) == 0.0_r8_kind) then
             buf_dist(k)=dist
             buf_index(k)=i
             buf_type(k)='epe'
             cycle scnd
          endif
       enddo
       max_ind=maxloc(buf_dist)
       max_val=maxval(buf_dist)
       if(dist < max_val) then
          buf_dist(max_ind)=dist
          buf_index(max_ind)=i
          buf_type(max_ind)='epe'
       endif
    enddo scnd

    tetra_imp%number(2:5)=buf_index
    tetra_imp%center(2:5)=buf_type
    do i=2,5
       j=tetra_imp%number(i)
       if(tetra_imp%center(i) == 'imp') then
          if(use_ref_data) then
             tetra_imp%coord(i,:)=reg_reference(j)%rs
          else
             tetra_imp%coord(i,:)=r_sh_imp(j,:)
          endif
       else
          tetra_imp%coord(i,:)=epe(j)%r
       endif
    enddo

  end subroutine select_centers
  !.......................................

  !.......................................
  subroutine vac_3_body(e3bv)
    real(kind=r8_kind), intent(out) :: e3bv
    real(kind=r8_kind) :: theta,cos_th,rss1,rss2,scal1
    real(kind=r8_kind) :: deg2rad
    integer(kind=i4_kind) :: i,j,k,ia_1,ia_2,ia_3,i1,j1,k1,l

    deg2rad=pi/180.0_r8_kind
    e3bv=zero
    if( n_vacancies == 0 ) return
    fst:  do i=1,n_tetrahedrons
       do l=1,5
          if (tetra_atoms(i,l) <= n_vacancies) goto 1
       enddo
       cycle fst
1      ia_1=tetra_atoms(i,1)
       i1=epe(ia_1)%k
       scnd: do j=2,4
          ia_2=tetra_atoms(i,j)
          if(ia_2 == 0) cycle scnd
          j1=epe(ia_2)%k
          thrd: do k=j+1,5
             ia_3=tetra_atoms(i,k)
             if(ia_3 == 0) cycle thrd
             if(ia_1 > n_vacancies.and.ia_2 >n_vacancies .and.ia_3 > n_vacancies) cycle thrd
             k1=epe(ia_3)%k

             rss1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
                  epe(ia_2)%s(:)-epe(ia_1)%s(:))
!!$             if(ia_1 <= n_vacancies.and.ia_2 <= n_vacancies) &
!!$                  rss1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_2)%r(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 <= n_vacancies.and.ia_2 > n_vacancies) &
!!$                  rss1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_2)%s(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_2 <= n_vacancies) &
!!$                  rss1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_2)%r(:)-epe(ia_1)%s(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_2 > n_vacancies) &
!!$                  rss1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_2)%s(:)-epe(ia_1)%s(:))

             rss2=dot_product(epe(ia_3)%s(:)-epe(ia_1)%s(:), &
                  epe(ia_3)%s(:)-epe(ia_1)%s(:))
!!$             if(ia_1 <= n_vacancies.and.ia_3 <= n_vacancies) &
!!$                  rss2=dot_product(epe(ia_3)%r(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_3)%r(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 <= n_vacancies.and.ia_3 > n_vacancies) &
!!$                  rss2=dot_product(epe(ia_3)%s(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_3)%s(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_3 <= n_vacancies) &
!!$                  rss2=dot_product(epe(ia_3)%r(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_3)%r(:)-epe(ia_1)%s(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_3 > n_vacancies) &
!!$                  rss2=dot_product(epe(ia_3)%s(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_3)%s(:)-epe(ia_1)%s(:))

             scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
                  epe(ia_3)%s(:)-epe(ia_1)%s(:))
!!$             if(ia_1 <= n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 <= n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_3)%r(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 <= n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 > n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_3)%s(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 <= n_vacancies.and.ia_2 > n_vacancies.and.ia_3 <= n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_3)%r(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 <= n_vacancies.and.ia_2 > n_vacancies.and.ia_3 > n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%r(:), &
!!$                  epe(ia_3)%s(:)-epe(ia_1)%r(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 <= n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_3)%r(:)-epe(ia_1)%s(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_2 <= n_vacancies.and.ia_3 > n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%r(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_3)%s(:)-epe(ia_1)%s(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_2 > n_vacancies.and.ia_3 <= n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_3)%r(:)-epe(ia_1)%s(:))
!!$
!!$             if(ia_1 > n_vacancies.and.ia_2 > n_vacancies.and.ia_3 > n_vacancies) &
!!$                  scal1=dot_product(epe(ia_2)%s(:)-epe(ia_1)%s(:), &
!!$                  epe(ia_3)%s(:)-epe(ia_1)%s(:))

             rss1=sqrt(rss1)
             rss2=sqrt(rss2)
             cos_th=scal1/(rss1*rss2)
             theta=acos(cos_th)
             e3bv=e3bv-0.5*ki(j1,i1,k1)*(theta-theta_0(j1,i1,k1)*deg2rad)**2
          enddo thrd!k=j+1,4
       enddo scnd!j=1,3
    enddo fst!i=1,n_tetra_atoms

  end subroutine vac_3_body
END SUBROUTINE energy_vac_imp
!-------------------------------------------------------------------

!-------------------------------------------------------------------
REAL(kind=r8_kind) FUNCTION erf(X)
! **error-function
! **MARK 4 RELEASE NAG COPYRIGHT 1974.
! **MARK 4.5 REVISED
! *MARK 8 REVISED. IER-221 (MAR 1980)
! **SINGLE PRECISION VERSION, FOR IBM

  integer(kind=i4_kind),parameter :: NCFC=18, &
                                     NCFD=17
  integer(kind=i4_kind) :: J
  real(kind=r8_kind),parameter ::  XUP=6.25, &
                                   SQRTPI=1.7724538509055160, &
                                   ZERO=0.0, &
                                   ONE=1.0, &
                                   TWO=2.0, &
                                   THREE=3.0, &
                                   TWENTY=20.0, &
                                   HALF=0.5
  real(kind=r8_kind) :: XV, X, X2, BJP2, BJP1, BJ
  REAL(kind=r8_kind) :: C(18)= &
      (/1.9449071068178803,4.20186582324414E-2,-1.86866103976769E-2  &
      ,5.1281061839107E-3,-1.0683107461726E-3,1.744737872522E-4  &
      ,-2.15642065714E-5,1.7282657974E-6,-2.00479241E-8  &
      ,-1.64782105E-8,2.0008475E-9,2.57716E-11,-3.06343E-11  &
      ,1.9158E-12,3.703E-13,-5.43E-14,-4.0E-15,1.2E-15/), &
                        D(17)= &
      (/1.4831105640848036,-3.010710733865950E-1,6.89948306898316E-2  &
      ,-1.39162712647222E-2,2.4207995224335E-3,-3.658639685849E-4  &
      ,4.86209844323E-5,-5.7492565580E-6,6.113243578E-7  &
      ,-5.89910153E-8,5.2070091E-9,-4.232976E-10,3.18811E-11  &
      ,-2.2361E-12,1.467E-13,-9.0E-15,5.0E-16/)
! **NO FAILURE EXITS
  XV = ABS(X)
  IF (XV.lt.XUP) then
    IF (XV.gt.TWO) then
      X2 = TWO-TWENTY/(XV+THREE)
! **SUMMATION
      BJP2 = ZERO
      BJP1 = C(NCFC)
      J = NCFC - 1
      do
        BJ = X2*BJP1 - BJP2 + C(J)
        IF (J.EQ.1) exit
        BJP2 = BJP1
        BJP1 = BJ
        J = J - 1
      enddo
      X2 = HALF*(BJ-BJP2)/XV*EXP(-X*X)/SQRTPI
      erf = (ONE-X2)*SIGN(ONE,X)
      return
    endif !XV.lt.XUP
    X2 = X*X - TWO
! **SUMMATION
    BJP2 = ZERO
    BJP1 = D(NCFD)
    J = NCFD - 1
    do
      BJ = X2*BJP1 - BJP2 + D(J)
      IF (J.EQ.1) exit
      BJP2 = BJP1
      BJP1 = BJ
      J = J - 1
    enddo
    erf = HALF*(BJ-BJP2)*X
    return
  endif !XV.lt.XUP
  erf = SIGN(ONE,X)
END function erf
!--------------------------------------------------------------------
!--------------------------------------------------------------------

real(kind=r8_kind) FUNCTION ERFC(X)
  !
  ! **COMPLEMETARY ERRORFUNCTION USING CHEBYSHEV APPROXIMATION (YL LUKE
  ! **1975 PP123-4) D,DD,SV AS IN PRESS, NUM. REC. 1988 "CHEBEV"
  !
  ! See also modules/fermi_module.f90, utilities/ewald.f90,
  ! utilities/ewald_new.f90
  !
  integer(kind=i4_kind),parameter :: NA=25, &
                                     NC=22
  real(kind=r8_kind),parameter ::  ZERO=0.0, &
                                   ONE=1.0, &
                                   TWO=2.0, &
                                   THREE=3.0, &
                                   FOUR=4.0, &
                                   HALF=0.5, &
                                   SQPI2=1.128379167095513

   real(kind=r8_kind) :: A(0:NA)=  &
       (/.109547129977762D+1, -.289175401126989D+0,  .110456398633795D+0,  &
        -.412531882278565D-1,  .140828380706516D-1, -.432929544743143D-2,  &
         .119827190159228D-2, -.299972962353249D-3,  .683258603788747D-4,  &
        -.142469884548677D-4,  .273540877283989D-5, -.048619128719754D-5,  &
         .008038727621172D-5, -.001241841831213D-5,  .000179953258879D-5,  &
        -.000024547948775D-5,  .000003162508603D-5, -.000000385902200D-5,  &
         .000000044720291D-5, -.000000004933613D-5,  .000000000519303D-5,  &
        -.000000000052258D-5,  .000000000005037D-5, -.000000000000466D-5,  &
         .000000000000041D-5, -.000000000000004D-5/), &
                         C(0:NC)=  &
       (/.975083423708556D+0, -.240493938504146D-1,  .820452240880432D-3,  &
        -.434293081303427D-4,  .301844703403493D-5, -.025447331925082D-5,  &
         .002485835302051D-5, -.000273172013238D-5,  .000033084722797D-5,  &
        -.000004350549080D-5,  .000000614121457D-5, -.000000092236928D-5,  &
         .000000014635665D-5, -.000000002439278D-5,  .000000000424976D-5,  &
        -.000000000077084D-5,  .000000000014507D-5, -.000000000002824D-5,  &
         .000000000000567D-5, -.000000000000117D-5,  .000000000000025D-5,  &
        -.000000000000005D-5,  .000000000000001D-5/)
    real(kind=r8_kind) :: X,D,DD,Z,ALPHA,SV
    integer(kind=i4_kind) :: J
    D=ZERO
    DD=ZERO
    IF(ABS(X).LT.THREE) THEN
! **CALCULATE VIA ERF
      Z=X/THREE
      ALPHA=TWO-FOUR*Z*Z
      DO J=NA,0,-1
         SV=D
         D=-ALPHA*D-DD+A(J)
         DD=SV
      ENDDO
      ERFC=ONE-SQPI2*Z*(D-DD)
    ELSE
! **CALCULATE DIRECTLY
      Z=ABS(THREE/X)
      ALPHA=TWO-FOUR*Z*Z
      DO J=NC,0,-1
        SV=D
        D=-ALPHA*D-DD+C(J)
        DD=SV
      ENDDO
      IF(X.GT.ZERO) THEN
        ERFC=HALF*EXP(-X*X)/X*(D+HALF*ALPHA*DD)*SQPI2
      ELSE
        ERFC=TWO-HALF*EXP(-X*X)/(-X)*(D+HALF*ALPHA*DD)*SQPI2
      ENDIF
    ENDIF
END function erfc
!----------------------------------------------------------------------

SUBROUTINE minv(A,N,D,L,M)

  real(kind=r8_kind),intent(out):: d
  real(kind=r8_kind) :: a(*),biga,hold
  integer(kind=i4_kind) :: n,l(3),m(3),nk,k,kk,iz,i,ij
  integer(kind=i4_kind) :: ki,ji,jp,jk,ik,kj,jq,jr,j
  D=1.0
  NK=-N
  DO K=1,N
    NK=NK+N
    L(K)=K
    M(K)=K
    KK=NK+K
    BIGA=A(KK)
    DO J=K,N
      IZ=N*(J-1)
      DO I=K,N
        IJ=IZ+I
        if(ABS(BIGA).ge.ABS(A(IJ))) cycle
        BIGA=A(IJ)
        L(K)=I
        M(K)=J
      enddo !I=K,N
    enddo !J=K,N
    J=L(K)
    if(j.gt.k) then
      KI=K-N
      DO I=1,N
        KI=KI+N
        HOLD=-A(KI)
        JI=KI-K+J
        A(KI)=A(JI)
        A(JI) =HOLD
      enddo !I=1,N
    endif
    I=M(K)
    IF(I.gt.K) then
      JP=N*(I-1)
      DO J=1,N
        JK=NK+J
        JI=JP+J
        HOLD=-A(JK)
        A(JK)=A(JI)
        A(JI) =HOLD
      enddo !J=1,N
    endif
    if(BIGA.eq.0.0) then
      D=0.0
      RETURN
    endif
    DO I=1,N
        if(i.eq.k) cycle
        IK=NK+I
        A(IK)=A(IK)/(-BIGA)
    enddo !I=1,N
    DO I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO J=1,N
        IJ=IJ+N
        IF(I.eq.K) cycle
        IF(J.eq.K) cycle
        KJ=IJ-I+K
        A(IJ)=HOLD*A(KJ)+A(IJ)
      enddo !J=1,N
    enddo !I=1,N
    KJ=K-N
    DO J=1,N
      KJ=KJ+N
      IF(J.eq.K) cycle
      A(KJ)=A(KJ)/BIGA
    enddo !J=1,N
    D=D*BIGA
    A(KK)=1.0/BIGA
  enddo !K=1,N
  K=N
  do
    K=(K-1)
    IF(K.le.0) exit
    I=L(K)
    IF(I.gt.K) then
      JQ=N*(K-1)
      JR=N*(I-1)
      DO J=1,N
        JK=JQ+J
        HOLD=A(JK)
        JI=JR+J
        A(JK)=-A(JI)
        A(JI) =HOLD
      enddo !J=1,N
    endif
    J=M(K)
    IF(J.gt.K) then
      KI=K-N
      DO I=1,N
        KI=KI+N
        HOLD=A(KI)
        JI=KI-K+J
        A(KI)=-A(JI)
        A(JI) =HOLD
      enddo !I=1,N
    endif
  enddo
END subroutine minv
!------------------------------------------
!------------------------------------------

REAL(kind=r8_kind) FUNCTION gauss_potential(RG,NDIM,NI)
! **procedure of calculation of gauss charge potential

  use epecom_module

  real(kind=r8_kind), intent(in),dimension(*) :: rg
  real(kind=r8_kind) :: X(3),psi,gix,GIX2
  integer(kind=i4_kind) :: i,i1,i2,i3,NDIM,NI
  PSI=zero
  I1=NI
  I2=NI+NDIM
  I3=NI+2*NDIM
  X(1)=RG(I1)
  X(2)=RG(I2)
  X(3)=RG(I3)

  DO I=1,n_bs_points
    gix=dot_product(GSTR(I,:),x(:))
    GIX2=GIX*PI2
    PSI=PSI+RSIN(I)*SIN(GIX2)+RCOS(I)*COS(GIX2)
  enddo !I=1,n_bs_points
  gauss_potential=PSI
END FUNCTION gauss_potential
!-----------------------------------------------------
!-----------------------------------------------------
real(KIND=r8_kind) FUNCTION reg2_coulomb()
! **program of the cell COULOMB part relaxation energy calc.

use epecom_module
use type_module

logical, allocatable::cutoff(:)
! real(kind=r8_kind):: anion_dir
real(kind=r8_kind) :: e2,e3
integer(kind=i4_kind) :: j,n,j1,n1,status
real(kind=r8_kind), allocatable, dimension(:)::dp_reg1
type dot_products
real(kind=r8_kind)::ss,sc,cs,cc,rs,rc,sr,cr,rr,sssr,cssr,rssr,rccr,sccr,cccr
end type dot_products
type (dot_products),allocatable, dimension(:):: dp
type sq_dot_products
real(kind=r8_kind)::rr,ss,cs,sc,cc,sr,cr,rs,rc
end type sq_dot_products
type (sq_dot_products),allocatable, dimension(:)::sqdp
type spec_ion
real(kind=r8_kind)::r(3),s(3),c(3),q,qs,qc
end type spec_ion
type (spec_ion),allocatable, dimension(:):: reg1

real(kind=r8_kind),pointer::temp(:)
allocate(dp_reg1(reg_I_n_ions),cutoff(reg_I_n_ions),stat=status)
if(status.ne.0) call error_handler("dp_reg1 allocate failed")

E2=zero
E3=zero


! anion_dir=0.0_r8_kind

DO J1=reg_I_n_ions+1,reg_2a_n_ions
   J=reg_I_n_ions+1+reg_2a_n_ions-J1

   dp_reg1(:)=(epe(:reg_I_n_ions)%r(1)-epe(J)%r(1))**2+ &
              (epe(:reg_I_n_ions)%r(2)-epe(J)%r(2))**2+ &
              (epe(:reg_I_n_ions)%r(3)-epe(J)%r(3))**2

   N=epe(J)%k

   where(dp_reg1(:).gt.RADIUS_LONG_INTERACT)
      cutoff=.false.
      elsewhere
      cutoff=.true.
   end where

   n1=count(cutoff)
   if(n1.ne.0) then
      allocate(dp(n1),sqdp(n1),reg1(n1))
      temp=>epe(:reg_I_n_ions)%r(1)
      reg1(:)%r(1)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%r(2)
      reg1(:)%r(2)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%r(3)
      reg1(:)%r(3)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%s(1)
      reg1(:)%s(1)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%s(2)
      reg1(:)%s(2)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%s(3)
      reg1(:)%s(3)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%c(1)
      reg1(:)%c(1)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%c(2)
      reg1(:)%c(2)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%c(3)
      reg1(:)%c(3)=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%q
      reg1(:)%q=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%qs
      reg1(:)%qs=pack(temp,cutoff)
      temp=>epe(:reg_I_n_ions)%qc
      reg1(:)%qc=pack(temp,cutoff)

      dp(:)%rr=pack(dp_reg1(:),cutoff)


      dp(:)%ss=(reg1(:n1)%s(1)-epe(J)%s(1))**2+ &
               (reg1(:n1)%s(2)-epe(J)%s(2))**2+ &
               (reg1(:n1)%s(3)-epe(J)%s(3))**2
      dp(:)%cs=(reg1(:n1)%c(1)-epe(J)%s(1))**2+ &
               (reg1(:n1)%c(2)-epe(J)%s(2))**2+ &
               (reg1(:n1)%c(3)-epe(J)%s(3))**2
      dp(:)%sc=(reg1(:n1)%s(1)-epe(J)%c(1))**2+ &
               (reg1(:n1)%s(2)-epe(J)%c(2))**2+ &
               (reg1(:n1)%s(3)-epe(J)%c(3))**2
      dp(:)%cc=(reg1(:n1)%c(1)-epe(J)%c(1))**2+ &
               (reg1(:n1)%c(2)-epe(J)%c(2))**2+ &
               (reg1(:n1)%c(3)-epe(J)%c(3))**2
      dp(:)%rs=(reg1(:n1)%r(1)-epe(J)%s(1))**2+ &
               (reg1(:n1)%r(2)-epe(J)%s(2))**2+ &
               (reg1(:n1)%r(3)-epe(J)%s(3))**2
      dp(:)%rc=(reg1(:n1)%r(1)-epe(J)%c(1))**2+ &
               (reg1(:n1)%r(2)-epe(J)%c(2))**2+ &
               (reg1(:n1)%r(3)-epe(J)%c(3))**2
      dp(:)%sr=(reg1(:n1)%s(1)-epe(J)%r(1))**2+ &
               (reg1(:n1)%s(2)-epe(J)%r(2))**2+ &
               (reg1(:n1)%s(3)-epe(J)%r(3))**2
      dp(:)%cr=(reg1(:n1)%c(1)-epe(J)%r(1))**2+ &
               (reg1(:n1)%c(2)-epe(J)%r(2))**2+ &
               (reg1(:n1)%c(3)-epe(J)%r(3))**2
      dp(:)%sssr=(reg1(:n1)%s(1)-epe(J)%s(1))*(epe(J)%s(1)-epe(J)%r(1))+ &
                 (reg1(:n1)%s(2)-epe(J)%s(2))*(epe(J)%s(2)-epe(J)%r(2))+ &
                 (reg1(:n1)%s(3)-epe(J)%s(3))*(epe(J)%s(3)-epe(J)%r(3))

      dp(:)%cssr=(reg1(:n1)%c(1)-epe(J)%s(1))*(epe(J)%s(1)-epe(J)%r(1))+ &
                 (reg1(:n1)%c(2)-epe(J)%s(2))*(epe(J)%s(2)-epe(J)%r(2))+ &
                 (reg1(:n1)%c(3)-epe(J)%s(3))*(epe(J)%s(3)-epe(J)%r(3))
      dp(:)%rssr=(reg1(:n1)%r(1)-epe(J)%s(1))*(epe(J)%s(1)-epe(J)%r(1))+ &
                 (reg1(:n1)%r(2)-epe(J)%s(2))*(epe(J)%s(2)-epe(J)%r(2))+ &
                 (reg1(:n1)%r(3)-epe(J)%s(3))*(epe(J)%s(3)-epe(J)%r(3))
      dp(:)%sccr=(reg1(:n1)%s(1)-epe(J)%c(1))*(epe(J)%c(1)-epe(J)%r(1))+ &
                 (reg1(:n1)%s(2)-epe(J)%c(2))*(epe(J)%c(2)-epe(J)%r(2))+ &
                 (reg1(:n1)%s(3)-epe(J)%c(3))*(epe(J)%c(3)-epe(J)%r(3))
      dp(:)%cccr=(reg1(:n1)%c(1)-epe(J)%c(1))*(epe(J)%c(1)-epe(J)%r(1))+ &
                 (reg1(:n1)%c(2)-epe(J)%c(2))*(epe(J)%c(2)-epe(J)%r(2))+ &
                 (reg1(:n1)%c(3)-epe(J)%c(3))*(epe(J)%c(3)-epe(J)%r(3))
      dp(:)%rccr=(reg1(:n1)%r(1)-epe(J)%c(1))*(epe(J)%c(1)-epe(J)%r(1))+ &
                 (reg1(:n1)%r(2)-epe(J)%c(2))*(epe(J)%c(2)-epe(J)%r(2))+ &
                 (reg1(:n1)%r(3)-epe(J)%c(3))*(epe(J)%c(3)-epe(J)%r(3))

      sqdp(:)%rr=SQRT(dp(:)%rr)
      sqdp(:)%ss=SQRT(dp(:)%ss)
      sqdp(:)%sc=SQRT(dp(:)%sc)
      sqdp(:)%cs=SQRT(dp(:)%cs)
      sqdp(:)%cc=SQRT(dp(:)%cc)
      sqdp(:)%sr=SQRT(dp(:)%sr)
      sqdp(:)%cr=SQRT(dp(:)%cr)
      sqdp(:)%rs=SQRT(dp(:)%rs)
      sqdp(:)%rc=SQRT(dp(:)%rc)


      E2=E2+sum( &
      reg1(:)%qs*(Q_SHELL(epe(J)%k)/sqdp(:)%ss+Q_NUCLEAR(epe(J)%k)/sqdp(:)%sc &
                   -Q_ion(epe(J)%k)*arrERF(ERROR_FUNCTION_PARAMETER*sqdp(:)%sr)/sqdp(:)%sr &
                                                                                          )  &
     +reg1(:)%qc*(Q_SHELL(epe(J)%k)/sqdp(:)%cs+Q_NUCLEAR(N)/sqdp(:)%cc  &
                  -Q_ion(epe(J)%k)*arrERF(ERROR_FUNCTION_PARAMETER*sqdp(:)%cr)/sqdp(:)%cr &
                                                                                          ) &
     -reg1(:)%q*(  Q_SHELL(epe(J)%k)/sqdp(:)%rs+Q_NUCLEAR(epe(J)%k)/sqdp(:)%rc &
                                                              -Q_ion(epe(J)%k)/sqdp(:)%rr) )

!    anion_dir=anion_dir &
!             +reg1(1)%qs*(Q_SHELL(epe(J)%k)/sqdp(1)%ss+Q_NUCLEAR(epe(J)%k)/sqdp(1)%sc &
!                         -Q_ion(epe(J)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(1)%sr)/sqdp(1)%sr &
!                                                                                         )  &
!             +reg1(1)%qc*(Q_SHELL(epe(J)%k)/sqdp(1)%cs+Q_NUCLEAR(N)/sqdp(1)%cc  &
!                         -Q_ion(epe(J)%k)*ERF(ERROR_FUNCTION_PARAMETER*sqdp(1)%cr)/sqdp(1)%cr &
!                                                                                         )


      E3=E3+sum(&
           reg1(:)%qs*(  Q_SHELL(N)*dp(:)%sssr/dp(:)%ss/sqdp(:)%ss +   &
                         Q_NUCLEAR(N)*dp(:)%sccr/dp(:)%sc/sqdp(:)%sc ) &
          +reg1(:)%qc*(  Q_SHELL(N)*dp(:)%cssr/dp(:)%cs/sqdp(:)%cs +   &
                         Q_NUCLEAR(N)*dp(:)%cccr/dp(:)%cc/sqdp(:)%cc)  &
          -reg1(:)%q *(  Q_SHELL(N)*dp(:)%rssr/dp(:)%rs/sqdp(:)%rs +   &
                         Q_NUCLEAR(N)*dp(:)%rccr/dp(:)%rc/sqdp(:)%rc ) )


      deallocate(dp,reg1,sqdp,stat=status)
      if(status.ne.0) call error_handler("deallocate reg1 failed")
   end if

enddo
!     DPRINT 'anion_dir contrib',anion_dir
!     DPRINT 'E2 first contrib',E2
!     DPRINT 'E2 second contrib',E3


deallocate(dp_reg1,cutoff)
reg2_coulomb=e2-e3*0.5_r8_kind
end FUNCTION reg2_coulomb

  FUNCTION arrerf(X) result(a)

  integer(kind=i4_kind),parameter :: NCFC=18, &
                                     NCFD=17
  integer(kind=i4_kind) :: J,k
  real(kind=r8_kind),parameter ::  XUP=6.25, &
                                   SQRTPI=1.7724538509055160, &
                                   ZERO=0.0, &
                                   ONE=1.0, &
                                   TWO=2.0, &
                                   THREE=3.0, &
                                   TWENTY=20.0, &
                                   HALF=0.5
  real(kind=r8_kind) :: XV,  X2, BJP2, BJP1, BJ
  real(kind=r8_kind) ,intent(in)::x(:)
  real(kind=r8_kind),dimension(size(x))::a
  REAL(kind=r8_kind) :: C(18)= &
      (/1.9449071068178803,4.20186582324414E-2,-1.86866103976769E-2  &
      ,5.1281061839107E-3,-1.0683107461726E-3,1.744737872522E-4  &
      ,-2.15642065714E-5,1.7282657974E-6,-2.00479241E-8  &
      ,-1.64782105E-8,2.0008475E-9,2.57716E-11,-3.06343E-11  &
      ,1.9158E-12,3.703E-13,-5.43E-14,-4.0E-15,1.2E-15/), &
                        D(17)= &
      (/1.4831105640848036,-3.010710733865950E-1,6.89948306898316E-2  &
      ,-1.39162712647222E-2,2.4207995224335E-3,-3.658639685849E-4  &
      ,4.86209844323E-5,-5.7492565580E-6,6.113243578E-7  &
      ,-5.89910153E-8,5.2070091E-9,-4.232976E-10,3.18811E-11  &
      ,-2.2361E-12,1.467E-13,-9.0E-15,5.0E-16/)
! **NO FAILURE EXITS
  do k=1,size(x)
  XV = ABS(X(k))
  IF (XV.lt.XUP) then
    IF (XV.gt.TWO) then
      X2 = TWO-TWENTY/(XV+THREE)
! **SUMMATION
      BJP2 = ZERO
      BJP1 = C(NCFC)
      J = NCFC - 1
      do
        BJ = X2*BJP1 - BJP2 + C(J)
        IF (J.EQ.1) exit
        BJP2 = BJP1
        BJP1 = BJ
        J = J - 1
      enddo
      X2 = HALF*(BJ-BJP2)/XV*EXP(-X(k)*X(k))/SQRTPI
      a(k) = (ONE-X2)*SIGN(ONE,X(k))
    else !XV.lt.XUP
    X2 = X(k)*X(k)- TWO
! **SUMMATION
    BJP2 = ZERO
    BJP1 = D(NCFD)
    J = NCFD - 1
    do
      BJ = X2*BJP1 - BJP2 + D(J)
      IF (J.EQ.1) exit
      BJP2 = BJP1
      BJP1 = BJ
      J = J - 1
    enddo
    a(k) = HALF*(BJ-BJP2)*X(k)
 end IF

  else !XV.lt.XUP
  a(k) = SIGN(ONE,X(k))
end IF
end do
END function arrerf
end module culon_module


