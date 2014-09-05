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
  subroutine calc_cpksQai(th_mo_coul,ilower,iupper)
    !
    !   Purpose:
    !   calculates G_ai_kl*S_kl of Qai using  cpks_gvec intermediate
    !
    !  Changed by subroutine:
    !  pert_vec         righthand-side for density fit
    !  pert_mat            coefficients
    !  pert_vec_xc      exchange fit projection for the extended MDA
    !
    !  Subroutine called by: pcorr_calc
    !
    !
#include "def.h"
    use comm_module
    use msgtag_module
    use symmetry_data_module
    use iounitadmin_module ! provides IO-units
    use readwriteblocked_module    ! contains the routines for buffered I/O
    use filename_module      ! set I/O-Filenames
    use prepare_integralfiles_module
    use options_module, only: options_integrals_on_file
    use cpksdervs_matrices
    use fit_coeff_module, only: get_fit,fit
    implicit none
    !------------ Declaration of formal parameters -----------------
    !  ilower/iupper : Range of charge fit functions to work on
    !  jlower/jupper : Range of exchange fit functions to work on
    integer(kind=i4_kind),           intent(in) :: ilower,iupper
    !** End of interface *****************************************


    !  Modifications:
    !  1. The double sum over orbitals has been changed considerably.
    !     the sum now runs ONLY over the lower triangular matrix.
    !     The reason for this change is, that the 3-center integrals
    !     only exist for unique pairs of orbitals (because the integrals
    !     are symmetric!). Thus if the indices run over the full square
    !     matrix, the integrals are not read in correctly.
    !     Unfortunately this means, that one if-statement has to be
    !     included. I have looked this up from the OLD lcgto - obviously
    !     its authors also didn`t find a workaround for the if-statement.
    !     If anybody has a better idea, please let me know.
    !     FN 29.8.95
    !

    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind)  :: i_i, i_s, &
         i_step, n_spin, n_irrep

    integer(kind=i4_kind)  :: i_grad ! required to make dummy allocation
    integer(kind=i4_kind)  :: occ_dim,vir_dim

    logical :: integrals_on_file
    integer(kind=i4_kind)  :: stat,k,mo
    real(kind=r8_kind),allocatable:: ai(:,:)

    type(readwriteblocked_tapehandle),intent(inout):: th_mo_coul
    integer(kind=i4_kind) :: klower, kupper, k_step, n_split

    n_split=n_cpks_split
    n_spin =symmetry_data_n_spin()
    n_irrep = symmetry_data_n_irreps()
    integrals_on_file = options_integrals_on_file()


#define mo_coul_onfile
#ifdef mo_coul_onfile
  call readwriteblocked_startread(trim(tmpfile('mo_coul.dat')), th_mo_coul)
#endif

#define split_co_ai
#ifndef split_co_ai
   do i_s=1,n_spin
   do i_i=1,n_irrep
             occ_dim=size(cpks(1,i_i,i_s)%Qai,1)
             vir_dim=size(cpks(1,i_i,i_s)%Qai,2)
#ifdef mo_coul_onfile
     allocate(ai(occ_dim,vir_dim),stat=stat)
     ASSERT(stat.eq.0)
     do k=1,iupper-ilower+1
     do mo=1,size(ai,1)
     call readwriteblocked_read(ai(mo,:),th_mo_coul)
     enddo
     do i_grad=1,size(cpks,1)
      cpks(i_grad,i_i,i_s)%Qai=cpks(i_grad,i_i,i_s)%Qai &
       + cpks_fitcoeff_grads(ilower+k-1,i_grad)*ai
     enddo
     enddo
    deallocate(ai,stat=stat)
    ASSERT(stat.eq.0)
#else

                do n1=1,occ_dim
                do n2=1,vir_dim
                 do i_grad=1,size(cpks,1)
                 cpks(i_grad,i_i,i_s)%Qai(n1,n2)=cpks(i_grad,i_i,i_s)%Qai(n1,n2) &
                 +dot_product(cpks_fitcoeff_grads(ilower:iupper,i_grad),cpks3c(i_i,i_s)%co_ai(:,n1,n2))
                enddo ! i_grad
                enddo
                enddo
#endif
   enddo
   enddo

#else
!i.e. split_co_ai

    ! number of fit functions assigned to this worker:
    i_step = iupper - ilower + 1

    !
    ! These iterations over fit function indices define format
    ! of the on-disk file mo_coul.dat and therefore have to be
    ! done identically at the place where it is written, see
    ! pert_coeff_module.f90p.
    !
    ! k_step == 1 means processing one fit-function at a time:
    k_step = i_step / n_split + 1
    kupper = 0
    do while ( kupper < i_step )
        ! process in pieces the range of k from klower to kupper:
        klower = kupper + 1
        kupper = kupper + k_step

        ! last piece may be smaller:
        kupper = min(kupper, i_step)

        ASSERT(kupper-klower>=0)
!       print *, "READ3: klower=", klower, "kupper=", kupper

        do i_s = 1, n_spin
            do i_i = 1, n_irrep
                occ_dim = size(cpks(1, i_i, i_s)%Qai, 1)
                vir_dim = size(cpks(1, i_i, i_s)%Qai, 2)

                allocate(ai(occ_dim, vir_dim), stat=stat)
                ASSERT(stat.eq.0)

                do k = klower, kupper
                    do mo = 1, size(ai, 1)
                        call readwriteblocked_read(ai(mo, :), th_mo_coul)
                    enddo

                    do i_grad = 1, size(cpks, 1)
                        cpks(i_grad, i_i, i_s)%Qai = cpks(i_grad, i_i, i_s)%Qai &
                                                   + cpks_fitcoeff_grads(ilower + k - 1, i_grad) * ai
                    enddo
                enddo

                deallocate(ai, stat=stat)
                ASSERT(stat.eq.0)
            enddo
        enddo
    enddo ! while
#endif

#ifdef mo_coul_onfile
     call readwriteblocked_stopread(th_mo_coul,'keep')
#endif

  end subroutine calc_cpksQai
