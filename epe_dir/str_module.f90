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
module str_module
use type_module
use epecom_module, only: output_epe

implicit none
private
save

real(kind=r8_kind), public :: shift_cent_gener(3,3)

public :: reg2a_ml_displacements
!------------------------------------------------------
contains




SUBROUTINE reg2a_ml_displacements

! **the program makes initial mott-litleton displacements of all lattice
! **ions on the first step and then define displacments of ions in region
! **2A. The displacments are propotional to gradient of electrostatic field.

  use epecom_module
  use culon_module

  real(kind=r8_kind) :: shell_displacement(3),core_displacement(3),QEFF, &
       shell_fac,core_fac
  integer(kind=i4_kind) :: i_epe,iimp,ivac
  integer(kind=i4_kind) :: counter
! **displacment of ions in region 2a

  IF(basic_action.EQ.2.or.(N_IMPURITIES.EQ.0.and.n_VACANCIES.eq.0 &
        .and..not.ml_cluster_simulated)) return

  DO I_epe=reg_I_n_ions+1,reg_2a_n_ions
     shell_displacement(:)=zero
     core_displacement(:)=zero
     if(ml_tensors) then
        if(ml_cluster_simulated) then
           do counter =1,n_ml_cluster_simulators
              if(i_epe.eq.reg_2a_n_ions) &
                   write(output_epe,*)'  ml_cluster_simulators(counter)%q', &
                   ml_cluster_simulators(counter)%q

              QEFF=ml_cluster_simulators(counter)%q &
                   /sqrt(dot_product(ml_cluster_simulators(counter)%r-epe(I_epe)%r, &
                   ml_cluster_simulators(counter)%r-epe(I_epe)%r))**3

              shell_displacement(:)=shell_displacement(:)+ &
                   QEFF* matmul(ml_ten(epe(I_epe)%m)%s, &
                   epe(I_epe)%r- ml_cluster_simulators(counter)%r)

              core_displacement(:)=shell_displacement(:)+ &
                   QEFF*matmul(ml_ten(epe(I_epe)%m)%c-ml_ten(epe(I_epe)%m)%s, &
                   epe(I_epe)%r-ml_cluster_simulators(counter)%r)
           end do

        elseif(use_epe_pgdata) then
           shell_displacement(:)=shell_displacement(:)+13.5515309044_r8_kind* &
                matmul(ml_ten(epe(I_epe)%m)%s, reg_I_pg(I_epe)%gs(:))
           core_displacement(:)=shell_displacement(:)+13.5515309044_r8_kind* &
                matmul(ml_ten(epe(I_epe)%m)%c-ml_ten(epe(I_epe)%m)%s, &
                reg_I_pg(I_epe)%gc(:))
        end if
           if(i_epe.le.300) then
                write(output_epe,*) 'mld s c ', i_epe,&
                sqrt(dot_product(shell_displacement,shell_displacement)), &
                sqrt(dot_product(core_displacement,core_displacement))
               write(output_epe,*) ml_ten(epe(I_epe)%m)%s(1,1), ml_ten(epe(I_epe)%m)%c(1,1)
            endif
            if(.not.fixed_reg2a_relcontribs) then
              ! **cycle by vacancies and impurities
              ! **displacements of the lattice ions will be proportional to
              ! **gradients of induced  electristatic fild
              DO IVAC=1,N_VACANCIES

                 shell_fac=-(q_shell(epe(ivac)%k)/sqrt(dot_product( &
                      epe(IVAC)%s-epe(I_epe)%r, &
                      epe(IVAC)%s-epe(I_epe)%r))**3)
                 if (i_epe.eq.reg_2a_n_ions) then
                    write(output_epe,*)'ivac',shell_fac,dot_product( &
                         epe(IVAC)%s-epe(I_epe)%r, &
                         epe(IVAC)%s-epe(I_epe)%r), epe(IVAC)%s
                 end if
                 core_fac=-(q_nuclear(epe(ivac)%k)/sqrt(dot_product( &
                      epe(IVAC)%c-epe(I_epe)%r, &
                      epe(IVAC)%c-epe(I_epe)%r))**3)

                 shell_displacement(:)=shell_displacement(:)+ shell_fac*&
                      matmul(ml_ten(epe(I_epe)%k)%s , epe(I_epe)%r-epe(IVAC)%s) + &
                      core_fac*matmul(ml_ten(epe(I_epe)%k)%s , epe(I_epe)%r-epe(IVAC)%c)

                 core_displacement(:)=core_displacement(:) &
                      +shell_fac*matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                      epe(I_epe)%r-epe(IVAC)%s) &
                      +core_fac*matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                      epe(I_epe)%r-epe(IVAC)%c)

              enddo!IVAC=1,N_VACANCIES

              ! **cycle by impurities
              ! **displacments will be proportional to
              ! **electrostatic field gradients in positions R()
              DO iimp=1,N_IMPURITIES
                 shell_fac=Q_SH_IMPURITY(iimp)/sqrt(dot_product(R_SH_IMP(iimp,:)-epe(I_epe)%r, &
                      R_SH_IMP(iimp,:)-epe(I_epe)%r))**3
                 core_fac=Q_NUC_IMPURITY(iimp)/sqrt(dot_product(R_NUC_IMP(iimp,:)-epe(I_epe)%r, &
                      R_NUC_IMP(iimp,:)-epe(I_epe)%r))**3

                 shell_displacement(:)=shell_displacement(:)+ &
                      shell_fac* matmul(ml_ten(epe(I_epe)%k)%s,epe(I_epe)%r-R_SH_IMP(iimp,:)) + &
                      core_fac*  matmul(ml_ten(epe(I_epe)%k)%s,epe(I_epe)%r-R_NUC_IMP(iimp,:))
                 core_displacement(:)=core_displacement(:)+ &
                      shell_fac* matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                      epe(I_epe)%r-R_SH_IMP(iimp,:))+ &
                      core_fac *  matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                      epe(I_epe)%r-R_NUC_IMP(iimp,:))

                 if(use_epe_reference) then
                    shell_fac=Q_SH_IMPURITY(iimp)/sqrt(dot_product( &
                         reg_reference(iimp)%rs-epe(I_epe)%r, &
                         reg_reference(iimp)%rs-epe(I_epe)%r ))**3
                    core_fac=Q_NUC_IMPURITY(iimp)/sqrt(dot_product( &
                         reg_reference(iimp)%rc-epe(I_epe)%r, &
                         reg_reference(iimp)%rc-epe(I_epe)%r ))**3

                 end if
                    shell_displacement(:)=shell_displacement(:)+ &
                         shell_fac* matmul(ml_ten(epe(I_epe)%k)%s,epe(I_epe)%r-reg_reference(iimp)%rs) + &
                         core_fac*  matmul(ml_ten(epe(I_epe)%k)%s,epe(I_epe)%r-reg_reference(iimp)%rc)

                    core_displacement(:)=core_displacement(:)+ &
                         shell_fac* matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                         epe(I_epe)%r-reg_reference(iimp)%rs)+ &
                         core_fac *  matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                         epe(I_epe)%r-reg_reference(iimp)%rc)

                    shell_fac=-Q_SH_IMPURITY(iimp)/sqrt(dot_product(epe_reference(iimp)%rs-epe(I_epe)%r, &
                         epe_reference(iimp)%rs-epe(I_epe)%r ))**3
                    core_fac=-Q_NUC_IMPURITY(iimp)/sqrt(dot_product(epe_reference(iimp)%rc-epe(I_epe)%r, &
                         epe_reference(iimp)%rc-epe(I_epe)%r ))**3

                    shell_displacement(:)=shell_displacement(:)+ &
                         shell_fac* matmul(ml_ten(epe(I_epe)%k)%s,epe(I_epe)%r-epe_reference(iimp)%rs) + &
                         core_fac*  matmul(ml_ten(epe(I_epe)%k)%s,epe(I_epe)%r-epe_reference(iimp)%rc)
                    core_displacement(:)=core_displacement(:)+ &
                         shell_fac* matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                         epe(I_epe)%r-epe_reference(iimp)%rs)+ &
                         core_fac *  matmul(ml_ten(epe(I_epe)%k)%c-ml_ten(epe(I_epe)%k)%s, &
                         epe(I_epe)%r-epe_reference(iimp)%rc)

                 enddo

              end if

           else !calc with ml_fac
              if(.not.fixed_reg2a_relcontribs) then
              ! **cycle by vacancies and impurities
              ! **displacements of the lattice ions will be proportional to
              ! **gradients of induced  electristatic fild
              DO IVAC=1,N_VACANCIES

                 shell_fac=-(q_shell(epe(ivac)%k)/sqrt(dot_product( &
                      epe(IVAC)%s-epe(I_epe)%r, &
                      epe(IVAC)%s-epe(I_epe)%r))**3)
                 if (i_epe.eq.reg_2a_n_ions) then
                    write(output_epe,*)'ivac',shell_fac,dot_product( &
                         epe(IVAC)%s-epe(I_epe)%r, &
                         epe(IVAC)%s-epe(I_epe)%r), epe(IVAC)%s
                 end if
                 core_fac=-(q_nuclear(epe(ivac)%k)/sqrt(dot_product( &
                      epe(IVAC)%c-epe(I_epe)%r, &
                      epe(IVAC)%c-epe(I_epe)%r))**3)

                 shell_displacement(:)=shell_displacement(:)+ &
                      ml_fac(epe(I_epe)%k)%s * shell_fac*(epe(I_epe)%r-epe(IVAC)%s) + &
                      ml_fac(epe(I_epe)%k)%s * core_fac*(epe(I_epe)%r-epe(IVAC)%c)

                 core_displacement(:)=core_displacement(:) &
                      +(ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s) &
                      * shell_fac*(epe(I_epe)%r-epe(IVAC)%s) &
                      +(ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s) &
                      * core_fac*(epe(I_epe)%r-epe(IVAC)%c)

              enddo !IVAC=1,N_VACANCIES

              ! **cycle by impurities
              ! **displacments will be proportional to
              ! **electrostatic field gradients in positions R()
              DO iimp=1,N_IMPURITIES
                 shell_fac=Q_SH_IMPURITY(iimp)/sqrt(dot_product(R_SH_IMP(iimp,:)-epe(I_epe)%r, &
                      R_SH_IMP(iimp,:)-epe(I_epe)%r))**3
                 core_fac=Q_NUC_IMPURITY(iimp)/sqrt(dot_product(R_NUC_IMP(iimp,:)-epe(I_epe)%r, &
                      R_NUC_IMP(iimp,:)-epe(I_epe)%r))**3
                 shell_displacement(:)=shell_displacement(:)+ &
                      shell_fac* ml_fac(epe(I_epe)%k)%s*(epe(I_epe)%r-R_SH_IMP(iimp,:)) + &
                      core_fac*  ml_fac(epe(I_epe)%k)%s*(epe(I_epe)%r-R_NUC_IMP(iimp,:))
                 core_displacement(:)=core_displacement(:)+ &
                      shell_fac* (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)* &
                      (epe(I_epe)%r-R_SH_IMP(iimp,:))+ &
                      core_fac *  (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)* &
                      (epe(I_epe)%r-R_NUC_IMP(iimp,:))

                 if(use_epe_reference) then
                    shell_fac=Q_SH_IMPURITY(iimp)/sqrt(dot_product( &
                         reg_reference(iimp)%rs-epe(I_epe)%r, &
                         reg_reference(iimp)%rs-epe(I_epe)%r ))**3
                    core_fac=Q_NUC_IMPURITY(iimp)/sqrt(dot_product( &
                         reg_reference(iimp)%rc-epe(I_epe)%r, &
                         reg_reference(iimp)%rc-epe(I_epe)%r ))**3

                    if (i_epe.eq.reg_2a_n_ions) then
                       write(output_epe,*)'iimp',shell_fac,dot_product( &
                            reg_reference(iimp)%rs-epe(I_epe)%r, &
                            reg_reference(iimp)%rs-epe(I_epe)%r )
                    end if
                    shell_displacement(:)=shell_displacement(:)+ &
                         shell_fac* ml_fac(epe(I_epe)%k)%s*(epe(I_epe)%r-reg_reference(iimp)%rs) + &
                         core_fac*  ml_fac(epe(I_epe)%k)%s*(epe(I_epe)%r-reg_reference(iimp)%rc)

                    core_displacement(:)=core_displacement(:)+ &
                         shell_fac* (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)* &
                         (epe(I_epe)%r-reg_reference(iimp)%rs)+ &
                         core_fac *  (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)* &
                         (epe(I_epe)%r-reg_reference(iimp)%rc)

                    shell_fac=-Q_SH_IMPURITY(iimp)/sqrt(dot_product(epe_reference(iimp)%rs-epe(I_epe)%r, &
                         epe_reference(iimp)%rs-epe(I_epe)%r ))**3
                    core_fac=-Q_NUC_IMPURITY(iimp)/sqrt(dot_product(epe_reference(iimp)%rc-epe(I_epe)%r, &
                         epe_reference(iimp)%rc-epe(I_epe)%r ))**3
                    shell_displacement(:)=shell_displacement(:)+ &
                         shell_fac* ml_fac(epe(I_epe)%k)%s*(epe(I_epe)%r-epe_reference(iimp)%rs) + &
                         core_fac*  ml_fac(epe(I_epe)%k)%s*(epe(I_epe)%r-epe_reference(iimp)%rc)
                    core_displacement(:)=core_displacement(:)+ &
                         shell_fac* (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)* &
                         (epe(I_epe)%r-epe_reference(iimp)%rs)+ &
                         core_fac *  (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)* &
                         (epe(I_epe)%r-epe_reference(iimp)%rc)

                 end if
              enddo

           end if

  if(ml_cluster_simulated) then
     do counter =1,n_ml_cluster_simulators
        if(i_epe.eq.reg_2a_n_ions) &
             write(output_epe,*)'  ml_cluster_simulators(counter)%q', &
             ml_cluster_simulators(counter)%q

        QEFF=ml_cluster_simulators(counter)%q &
             /sqrt(dot_product(ml_cluster_simulators(counter)%r-epe(I_epe)%r, &
             ml_cluster_simulators(counter)%r-epe(I_epe)%r))**3
        shell_displacement(:)=shell_displacement(:)+ &
             ml_fac(epe(I_epe)%k)%s*QEFF*(epe(I_epe)%r- ml_cluster_simulators(counter)%r)
        core_displacement(:)=core_displacement(:)+ &
             (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)*QEFF* &
             (epe(I_epe)%r-ml_cluster_simulators(counter)%r)
     end do
     elseif(use_epe_pgdata) then
!!$        print*,I_epe
!!$        print*,reg_I_pg(I_epe)%gs(:)
!!$        print*,reg_I_pg(I_epe)%gc(:)
     shell_displacement(:)=shell_displacement(:)+13.5515309044_r8_kind* &
                ml_fac(epe(I_epe)%k)%s* reg_I_pg(I_epe)%gs(:)
!!$     print*,shell_displacement(:)
     core_displacement(:)=core_displacement(:)+13.5515309044_r8_kind* &
                (ml_fac(epe(I_epe)%k)%c-ml_fac(epe(I_epe)%k)%s)* reg_I_pg(I_epe)%gc(:)
     if(i_epe.le.300) &
          write(output_epe,*) 'mld', i_epe,&
          sqrt(dot_product(shell_displacement,shell_displacement)), &
          sqrt(dot_product(core_displacement,core_displacement))

  end if

end if

     R_SH_ION(I_epe,:)=epe(I_epe)%r+shell_displacement(:)
     R_NUC_ION(I_epe,:)=R_SH_ION(I_epe,:)+core_displacement(:)
     epe(i_epe)%s=epe(I_epe)%r+shell_displacement(:)
     epe(i_epe)%c=epe(I_epe)%s+core_displacement(:)

  enddo

END SUBROUTINE reg2a_ml_displacements

end module str_module








