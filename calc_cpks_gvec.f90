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
subroutine calc_cpks_gvec(ilower, iupper)
  !
  !   Purpose:
  !   calculates cpks_gvec as intermediate in calculation G_ai_kl*S_kl
  !   cpks_gvec is part of fit coeff 1st derivatives
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
  use eigen_data_module,only: eigvec
  use prepare_integralfiles_module
  use integralstore_module, only: integralstore_3c_co
  use options_module, only: options_integrals_on_file
  use cpksdervs_matrices
  use fit_coeff_module, only: get_fit,fit
  use error_module
  USE_MEMLOG
  implicit none
  !------------ Declaration of formal parameters -----------------
  !  ilower/iupper : Range of charge fit functions to work on
  integer(i4_kind), intent(in) :: ilower, iupper
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
  type(fit)                   :: n_fit

  integer(kind=i4_kind)  :: i_i, i_s, i_a, i_b, n1, n2, &
       i_meta, i_last, i_step, n_spin, n_irrep, n_proj

  integer(kind=i4_kind)  :: i_grad ! required to make dummy allocation
  integer(kind=i4_kind)  :: n_partners

  real(kind=r8_kind), pointer :: coul_int(:)
  type(readwriteblocked_tapehandle) :: th_ch
  logical :: integrals_on_file

  real(kind=r8_kind) :: fact
  real(kind=r8_kind),allocatable :: cpks_gvec_temp(:,:),co_ii(:)
  real(kind=r8_kind), allocatable:: co_i(:,:,:)
  integer(kind=i4_kind)  :: eig_dim,occ_dim
  integer(kind=i4_kind)  :: k_step, klower, kupper, max_co_i
  integer(i4_kind), parameter :: n_split = 4
  character(len=256) :: trace_message

  n_spin =symmetry_data_n_spin()
  n_irrep = symmetry_data_n_irreps()
  integrals_on_file = options_integrals_on_file()
  n_proj = 1

  call get_fit(n_fit)

  i_step = iupper - ilower + 1
  allocate (cpks_gvec_temp(i_step,size(cpks,1)),co_ii(i_step), &
       stat=cpksalloc(13))
  ASSERT(cpksalloc(13).eq.0)
  MEMLOG(size(cpks_gvec_temp)+size(co_ii))

  allocate (cpks_gvec(n_fit%n_ch,size(cpks,1)),stat=cpksalloc(3))
  ASSERT(cpksalloc(3).eq.0)
  MEMLOG(size(cpks_gvec))
  cpks_gvec=0.0_r8_kind

  !
  ! FIXME: integer  arithmetics can be  tricky for zero  sized arrays,
  ! i_step == 0.   Here k_step == 1 means  processing one fit-function
  ! at a time. This number should be positive to ensure progress:
  !
  k_step = i_step / n_split + 1

  !
  ! This appears to be for debug printing only:
  !
  if(comm_i_am_master()) then
     max_co_i=0
     do i_i=1,n_irrep
        eig_dim=size(eigvec(i_i)%m,1)
        occ_dim=cpks(1,i_i,1)%occ_dim
        max_co_i=max(max_co_i,(k_step*occ_dim*eig_dim*8)/1000)
     enddo
     write(trace_message,*) max_co_i
     call write_to_trace_unit(" split gvec co_i"//trim(trace_message))
  endif

  !
  ! Process in pieces klower:kupper the range of k from 0 to i_step:
  !
  kupper = 0
  do while (kupper < i_step)
     klower = kupper + 1
     kupper = kupper + k_step

     ! last piece may be smaller:
     kupper = min(kupper, i_step)

     ASSERT(kupper-klower>=0)
!!$     print *, "READ3: klower=", klower, "kupper=", kupper

     spins: do i_s=1,n_spin
        if (integrals_on_file) then
           allocate(coul_int(i_step),STAT=cpksalloc(6))
           ASSERT(cpksalloc(6).eq.0)
           MEMLOG(size(coul_int))
           call prepare_integral_open('coul', th_ch)
        else
           i_meta = 1
        endif

        if( .not. integrals_on_file )then
           ASSERT(allocated(integralstore_3c_co))
        endif

        irrs: do i_i=1,n_irrep

           ASSERT(allocated(eigvec(i_i)%m))
           eig_dim=size(eigvec(i_i)%m,1)
           occ_dim=cpks(1,i_i,i_s)%occ_dim
           !
           ! FIXME: this is the largest data structure here, O(N^3), formally:
           !
           allocate(co_i(kupper - klower + 1, occ_dim, eig_dim), stat=cpksalloc(168))
           ASSERT(cpksalloc(168).eq.0)
           co_i=0.0_r8_kind


           n_partners=symmetry_data_n_partners(i_i)

           do i_a=1,symmetry_data_dimension(i_i)
              do i_b=1,i_a
                 if (i_a.eq.i_b) then
                    fact = 1.0_r8_kind
                 else
                    fact = 2.0_r8_kind
                 endif

                 if (integrals_on_file) then
                    call readwriteblocked_read(coul_int,th_ch)
                 else
                    i_last = i_meta + i_step - 1
                    coul_int => integralstore_3c_co(i_meta:i_last)
                    i_meta = i_last + 1
                 endif

                 do n1=1,occ_dim
                    co_i(:, n1, i_b) = co_i(:, n1, i_b) &
                         + fact * coul_int(klower:kupper) * eigvec(i_i)%m(i_a, n1, i_s)
                    co_i(:, n1, i_a) = co_i(:, n1, i_a) &
                         + fact * coul_int(klower:kupper) * eigvec(i_i)%m(i_b, n1, i_s)
                 enddo

              enddo! i_b loop
           enddo!i_a loop

           do i_grad=1,size(cpks,1)
              cpks_gvec_temp(klower:kupper,i_grad)=0.0_r8_kind
           enddo

           do n1=1,occ_dim
              do n2=1,occ_dim
                 co_ii(klower:kupper)=0.0_r8_kind
                 do i_a=1,size(co_i,3)
                    !co_ii(n1,n2)
                    co_ii(klower:kupper) = co_ii(klower:kupper) &
                         + co_i(:,n1,i_a) * n_partners * eigvec(i_i)%m(i_a, n2, i_s)
                 enddo
                 do i_grad=1,size(cpks,1)
                    cpks_gvec_temp(klower:kupper, i_grad) = cpks_gvec_temp(klower:kupper, i_grad) &
                         + co_ii(klower:kupper) * cpks(i_grad, i_i, i_s)%s1(n1, n2)
                 enddo
              enddo
           enddo

           do i_grad=1,size(cpks,1)
              cpks_gvec(:,i_grad)=cpks_gvec(:,i_grad)+&
                   matmul(cpks_gvec_temp(klower:kupper,i_grad), &
                   cpks_gmat(ilower+klower-1:ilower+kupper-1,:))
           enddo

           deallocate(co_i,stat=cpksalloc(168))
           ASSERT(cpksalloc(168).eq.0)
           cpksalloc(168)=1
        enddo irrs

        if (integrals_on_file) then
           call prepare_integral_close(th_ch)
           MEMLOG(-size(coul_int))
           deallocate(coul_int,STAT=cpksalloc(6))
           ASSERT(cpksalloc(6).eq.0)
           cpksalloc(6)=1
        else
           nullify(coul_int)
        endif
     enddo spins
  enddo

  MEMLOG(-size(cpks_gvec_temp)-size(co_ii))
  deallocate (cpks_gvec_temp,co_ii, stat=cpksalloc(13))
  ASSERT(cpksalloc(13).eq.0)
  cpksalloc(13)=1

end subroutine calc_cpks_gvec
