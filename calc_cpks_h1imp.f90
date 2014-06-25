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
subroutine calc_cpks_h1imp(ilower, iupper)
  !
  !   Purpose:
  !   calculates G_ai_kl*S_kl of Qai using  cpks_gvec intermediate
  !
  !  Changed by subroutine:
  !
  !  Subroutine called by: cpks_g4constructs()
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
  use integralstore_module, only: integralstore_3c_co
  use options_module, only: options_integrals_on_file
  use cpksdervs_matrices
  use fit_coeff_module, only: get_fit, fit !, coeff_charge
  use eigen_data_module, only: eigvec
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

  integer(kind=i4_kind)  :: i_i,i_s,i_a,i_b,n1,n2, &
       i_meta, i_last, i_step,  n_spin, n_irrep

  integer(kind=i4_kind)  :: i_grad ! required to make dummy allocation
  integer(kind=i4_kind)  :: occ_dim,eig_dim
  integer(i4_kind) :: k_step, klower, kupper
  integer(i4_kind), parameter :: n_split = 4

  real(kind=r8_kind),pointer :: coul_int(:), eigvec1(:)!!!, eigvec2(:)
  type(readwriteblocked_tapehandle) :: th_ch
  logical ::  integrals_on_file

  real(kind=r8_kind) :: fact
  real(kind=r8_kind),allocatable :: co_i(:,:,:),co_n2i(:,:)

  n_spin =symmetry_data_n_spin()
  n_irrep = symmetry_data_n_irreps()
  integrals_on_file = options_integrals_on_file()
  i_step = iupper - ilower + 1

  !
  ! FIXME: integer  arithmetics can be  tricky for zero  sized arrays,
  ! i_step == 0.   Here k_step == 1 means  processing one fit-function
  ! at a time. This number should be positive to ensure progress:
  !
  k_step = i_step / n_split + 1

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
           allocate(coul_int(i_step),STAT=cpksalloc(109))
           ASSERT(cpksalloc(109).eq.0)
           MEMLOG(size(coul_int))
           call prepare_integral_open('coul', th_ch)
        else
           i_meta = 1
        endif

        !  print*,'h1 imp: cpks_fitcoeff_grads(:,1)',sum(cpks_fitcoeff_grads(:,1)*coeff_charge)
        do i_i=1,n_irrep

           occ_dim=size(cpks(1,i_i,i_s)%h1,1)
           eig_dim=size(eigvec(i_i)%m,1)
           !
           ! FIXME: this is the largest data structure here, O(N^3), formally:
           !
           allocate(co_i(kupper - klower + 1, occ_dim, eig_dim), &
                co_n2i(kupper - klower + 1, occ_dim), stat=cpksalloc(161))
           ASSERT(cpksalloc(161).eq.0)
           co_i=0.0_r8_kind

           do i_a=1,symmetry_data_dimension(i_i)
              do i_b=1,i_a
                 if (i_a.eq.i_b) then
                    fact = 0.5_r8_kind
                 else
                    fact = 1.0_r8_kind
                 endif

                 if (integrals_on_file) then
                    call readwriteblocked_read(coul_int,th_ch)
                 else
                    i_last = i_meta + i_step - 1
                    coul_int => integralstore_3c_co(i_meta:i_last)
                    i_meta = i_last + 1
                 endif


                 do n1=1,occ_dim
                    eigvec1 => eigvec(i_i)%m(:,n1,i_s)
                    co_i(:,n1,i_b)=co_i(:,n1,i_b)+fact*coul_int(klower:kupper) * eigvec1(i_a)
                    co_i(:,n1,i_a)=co_i(:,n1,i_a)+fact*coul_int(klower:kupper) * eigvec1(i_b)
                 enddo

              enddo! i_b loop
           enddo!i_a loop

           do n2=1,occ_dim
              co_n2i=0.0_r8_kind
              do i_a=1,size(co_i,3)
                 co_n2i(:,:)=co_n2i(:,:)+co_i(:,:occ_dim,i_a)*eigvec(i_i)%m(i_a,n2,i_s)
              enddo
              do i_grad=1,size(cpks,1)
                 do n1=1,occ_dim
                    cpks(i_grad,i_i,i_s)%h1(n1,n2)=cpks(i_grad,i_i,i_s)%h1(n1,n2)+ &
                         dot_product(cpks_fitcoeff_grads(ilower+klower-1:ilower+kupper-1,i_grad),co_n2i(:,n1))
                 enddo
              enddo
           enddo

           deallocate(co_i,co_n2i,stat=cpksalloc(161))
           ASSERT(cpksalloc(161).eq.0)
           cpksalloc(161)=1
        enddo! i_ir

        if (integrals_on_file) then
           call prepare_integral_close(th_ch)
           MEMLOG(-size(coul_int))
           deallocate(coul_int,STAT=cpksalloc(109))
           ASSERT(cpksalloc(109).eq.0)
           cpksalloc(109)=1
        else
           nullify(coul_int)
        endif
     enddo  spins
  enddo

end subroutine calc_cpks_h1imp
