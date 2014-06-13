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
!===============================================================
! Public interface of module
!===============================================================
subroutine cpks_g4constructs()
  ! output :   cpks%b

# define b_on_file
  ! define FPP_TIMERS 2
# include "def.h"
  use type_module    ! standard data types
  use interfaces, only: calc_cpks_gvec, calc_cpks_h1imp
  use symmetry_data_module
  use datatype       ! user defined data types
  use comm_module
  use comm, only: comm_allreduce
  use msgtag_module
  use iounitadmin_module  ! provides io units
  use print_module   !interface for outroutines
  use filename_module      ! set I/O-Filenames
  use occupation_module, only: get_n_elec
  use eigen_data_module, ONLY : eigval,eigvec
  use fit_coeff_module
  use bounds_module, only: bounds_ch, bounds_calc
  use pairs_module
  use pert_coeff_module
  use linsys_module
  use init_module, only: init
  use mat_charge_module
  use virtual_levels_module, only: eigvec_vir_dealloc
  use occupied_levels_module, ONLY: eigval_occ
  use calc3c_switches
  use cpksdervs_matrices
  use xc_cntrl, only: xc_is_on=>is_on, xc_ANY
  use post_scf_module
  use gradient_data_module, only: cpks_gradient_totalsym, gradient_index
  use comm, only: comm_bcast
  use grid_module,only:grid_close
  use operations_module, only: operations_solvation_effect !!!!!!!!!!!!AS
  use solv_2nd_deriv_module, only: cpks_solv_main,cpks_solv_Qai,cpks_solv_H1,cpks_solv_H1_1 !!!!!!!!!!!!!!!AS
  use readwriteblocked_module
  USE_MEMLOG

  implicit none

  type(fit) :: n_fit

  integer(kind=i4_kind) :: n_ai, n_occ, n_vir, j_cpks, n_holes

  ! ilower/iupper : charge fit functions
  integer(kind=i4_kind) :: info, iupper, ilower,&
       i_i,i_g,i_spn, &
       n_spin, n_irrep, n_ch, n_xc, n_mat, n_proj
  integer(kind=i4_kind) :: i_grad
  real(kind=r8_kind) :: target_charge
  real(kind=r8_kind),dimension(:,:),allocatable :: d1_temp
  real(kind=r8_kind),allocatable :: QxCOai(:),GxQxCOai(:,:),cpks_gamma(:,:),cpks_betta(:,:,:)
  ! the core fit functions are distributed precisely as the charge fit functions

  external error_handler

  real(kind=r8_kind), allocatable:: help_cpks_gvec(:,:)
  integer(kind=i4_kind):: k
  real(kind=r8_kind):: gamma_eps=10.0e-10
  integer(kind=i4_kind)::msgtag_start_cpksgvec=186
  integer(kind=i4_kind)::msgtag_done_cpksgvec=188
  type(readwriteblocked_tapehandle):: cpks_b_th(0:n_cpks+1)
  type(readwriteblocked_tapehandle):: th_mo_coul,th_s1,th_h1ai,th_h1
  integer(kind=i4_kind)::n_split
  character(len=2) :: b_fn
  logical:: cpks_ser_limited
  ! ----------------------------------------------------------------

  ! context: all processors
  FPP_TIMER_START(t_cpks_const)
  FPP_TIMER_START(t_cpks_init)
  call print_cpksalloc()
  MEMSET(0)

  n_split=n_cpks_split
  call get_cpks_dims() ! get sizes before h1ai deallocate

  if(.not.comm_i_am_master()) then
     call dealocate_h1ai_s1ai()
  else
     call tofile_h1es1_ai() ! deallocates h1ai and s1ai
  endif
  call say(" cpks h1ai and s1ai deallocated ")

  !  s1 and h1 here is stil in memory
  call tofile_h1(th_h1)
  call print_cpksalloc()

  comm_dervs_xc=.false.
  n_spin = symmetry_data_n_spin()
  n_irrep = symmetry_data_n_irreps()
  n_proj = 1

  cpks(:,:,:)%n_fin_cyc=-1  ! n_fin_cyc nedds to be defined
  ! this is a mask, if ge.0 then these
  ! data are not treated futher on
  cpks_gradient_index=>gradient_index

  call get_fit(n_fit)
  n_ch = n_fit%n_ch
  n_xc = n_fit%n_xc
  n_mat = (n_ch*(n_ch+1))/2

  if(comm_i_am_master()) then
     ! FIXME: cpks_gamma(0:n_cpks+1,...) and at the end of sub
     ! there is s access to 1+i_cpks element of it
     allocate(cpks_gamma(0:n_cpks+1,size(cpks,1)), &
          cpks_betta(0:n_cpks,0:n_cpks,size(cpks,1)), &
          stat=cpksalloc(118))
     ASSERT(cpksalloc(118).eq.0)
     MEMLOG(size(cpks_gamma)+size(cpks_betta))
     cpks_gamma=0.0_r8_kind
     cpks_betta=0.0_r8_kind
  endif



  if(comm_i_am_master()) call hole_corrected_s1() ! s1 is first required here

  ! load eigvec_vir and eigval_vir on each slave
  ASSERT(allocated(eigval))
  ASSERT(allocated(eigvec))

  do i_g = 1, n_irrep
     call comm_bcast(eigvec(i_g)%m)
     do i_grad = 1, size(cpks,1)
        do i_spn = 1, n_spin
           call comm_bcast(cpks(i_grad, i_g, i_spn)%s1)
        enddo
     enddo
  enddo

  call tofile_s1(th_s1)
  call print_cpksalloc

  if(.not.allocated(cpks_gmat)) then
     allocate(cpks_gmat(n_ch,n_ch),stat=cpksalloc(4))
     MEMLOG(size(cpks_gmat))
     ASSERT(cpksalloc(4).eq.0)
  endif

  if(comm_i_am_master()) call calc_cpks_gmat()
  call comm_bcast(cpks_gmat)

  call fromfile_s1(th_s1)

#if 1 /* calculate cpks_gvec */

  master_distribute_work: if(comm_i_am_master()) then
     call say(" calculate cpks_gvec ")

     ! FIXME: shouldnt it have been already done?
     call bounds_calc(n_fit)

     ! First start the slaves
     to_slaves: if (comm_parallel() ) then
        !
        ! Skip the share of master:
        !
        ilower = bounds_ch(1) + 1

        slaves :do i_i = 2, comm_get_n_processors()
           iupper = ilower + bounds_ch(i_i) - 1

           call comm_init_send(i_i,msgtag_start_cpksgvec) !msgtag_cpks_gvec) start calc cpks_gvec
           call commpack(ilower,info)
           ASSERT(info.eq.0)
           call commpack(iupper,info)
           ASSERT(info.eq.0)
           call comm_send()

           !
           ! Move forward:
           !
           ilower = ilower + bounds_ch(i_i)
        enddo slaves
     endif to_slaves

     ! Then start the master
     ! pert_vec actually holds the charge density projections [rho_spin|f_k]

     if( comm_parallel() ) then
        ilower = 1
        iupper = bounds_ch(1)
     else
        ilower=1
        iupper=n_ch
     endif

     call calc_cpks_gvec(ilower, iupper) ! now on master,  depends on cpks%s1

#ifndef memprof
     if(comm_parallel()) then ! receive cpks_gvec  contribs from slave
#endif
        allocate(help_cpks_gvec(size(cpks_gvec,1),size(cpks_gvec,2)),stat=cpksalloc(124))
        ASSERT(cpksalloc(124).eq.0)
        MEMLOG(size(help_cpks_gvec))
#ifndef memprof
        from_slaves :do i_i=2,comm_get_n_processors()
           call comm_save_recv(i_i,msgtag_packed_message)
           call communpack(help_cpks_gvec(1,1),size(cpks_gvec),1,info)
           ASSERT(info.eq.0)
           cpks_gvec=cpks_gvec+help_cpks_gvec
        enddo from_slaves
#endif

        MEMLOG(-size(help_cpks_gvec))
        deallocate(help_cpks_gvec,stat=cpksalloc(124))
        ASSERT(cpksalloc(124).eq.0)
        cpksalloc(124)=1
#ifndef memprof
        call comm_init_send(comm_all_other_hosts,msgtag_done_cpksgvec) !all slave contribs in
        call comm_send
     endif
#endif

  else master_distribute_work !i.e. now on slave

     call comm_save_recv(comm_master_host,msgtag_start_cpksgvec)
     call communpack(ilower,info)
     ASSERT(info.eq.0)
     call communpack(iupper,info)
     ASSERT(info.eq.0)

     call calc_cpks_gvec(ilower, iupper) ! s1 imp now  on slave

     !
     !      result: cpks_gvec to be used in Q calcs
     !
     !  occupied orbital coefficiens direvatives contibs calculated with
     !  use of 3 center inegrals and overlap derivative s1 matrix
     !  array of (n_ch,n_grad) shape where all occ MO contribs summed up
     !
     !  result to be communicated to all hosts to be used for Qai calcs
     call comm_init_send(comm_master_host,msgtag_packed_message)
     call commpack(cpks_gvec(1,1),size(cpks_gvec),1,info)
     call comm_send()
     call comm_save_recv(comm_master_host,msgtag_done_cpksgvec) !master done with cpks_gvec

  endif master_distribute_work !/else
#endif /* calculate cpks_gvec */



  FPP_TIMER_STOP(t_cpks_init)
  DPRINT 't_cpks_init calc_cpks_gvec=', FPP_TIMER_VALUE(t_cpks_init)

  FPP_TIMER_START(t_cpks_init)

  call say(" calculate co_ai ")
#define mo_coul_onfile
#define no_co_ai

#ifndef no_co_ai
  call allocate_co_ai()
#endif

  call calc_cpks3c_ai(ilower,iupper, th_mo_coul)
#define split_co_ai

#ifdef mo_coul_onfile
#ifndef no_co_ai
  call write_mo_coul()
  call deallocate_co_ai() !!!!!!!
#endif
  call say(" cpks mo_coul stored on file ")
#endif
  FPP_TIMER_STOP(t_cpks_init)
  DPRINT 't_cpks_init 3c_ai calculated', FPP_TIMER_VALUE(t_cpks_init)




  call allocate_Q_B() ! in reg run B is on file
#if 0
  call allocate_AB()

  DPRINT 'allocate_Q_AB_B_HBH done'
  print*,'max_cpks_coul_grad_size  co_ai_size cpks_size cpks_fitmat_size ', &
       max_cpks_coul_grad_size, co_ai_size,cpks_size,cpks_fitmat_size
#endif

  FPP_TIMER_START(t_cpks_init)

!!$        DPRINT '    '
!!$        DPRINT  'sum s1(i_grad) and sum abs s1(i_grad) to calc xc part of GxSkl'
!!$        do i_spn=1,size(cpks,3)
!!$           do i_grad=1,size(cpks,1)
!!$!!!          cpks(i_grad,1,i_spn)%h1=0.0_r8_kind
!!$              DPRINT  sum(cpks(i_grad,1,i_spn)%s1),sum(abs(cpks(i_grad,1,i_spn)%s1)),cpks(i_grad,1,i_spn)%s1(1,1)
!!$              !            if(size(cpks(i_grad,1,i_spn)%s1).gt.1)  &
!!$              !                 print*,cpks(i_grad,1,i_spn)%s1(2,2),cpks(i_grad,1,i_spn)%s1(1,2),cpks(i_grad,1,i_spn)%s1(2,1)
!!$           enddo
!!$        enddo
!!$
!!$
!!$        DPRINT '      '

  call fromfile_h1(th_h1)

  ![[== XC call 1 ===================================================
  if( xc_is_on(xc_ANY) )then
     FPP_TIMER_START(t_dervs_xc)
     call say(" cpks Qxc factors ")
     call print_cpksalloc()
     cpks_Qxc_calculated=.false.
     call cpks_xc_main(21)  ! output: Qai xc contrib, func h1 contrib
     cpks_Qxc_calculated=.true.
     FPP_TIMER_STOP(t_dervs_xc)

     DPRINT 't_dervs_xc Qxc',FPP_TIMER_VALUE(t_dervs_xc)

#if 0
     call print_Qai()
#endif

!!$          DPRINT ' h1 xc contrib',hai
!!$          do i_grad=1,size(cpks,1)
!!$             if(n_spin.eq.2) then
!!$                DPRINT sum(cpks(i_grad,1,1)%h1),sum(cpks(i_grad,1,2)%h1) ,i_grad
!!$                DPRINT cpks(i_grad,1,1)%h1(1,1),cpks(i_grad,1,2)%h1(1,1)
!!$             else
!!$                DPRINT sum(cpks(i_grad,1,1)%h1)
!!$             endif
!!$          enddo
  endif
  !]]== eof XC call 1 ===============================================

  FPP_TIMER_STOP(t_cpks_init)

  FPP_TIMER_START(t_cpks_init)

#if 1
  call allocate_AB()

  DPRINT 'allocate_Q_AB_B_HBH done'
  DPRINT 'max_cpks_coul_grad_size=', max_cpks_coul_grad_size
  DPRINT 'co_ai_size=', co_ai_size
  DPRINT 'cpks_size=', cpks_size
  DPRINT 'cpks_fitmat_size=', cpks_fitmat_size
#endif
  allocate (cpks_fitcoeff_grads(n_fit%n_ch,size(cpks,1)),stat=cpksalloc(108))
  !            MEMLOG(size(cpks_fitcoeff_grads)) ! deallocated in main_integral
  ASSERT(cpksalloc(108).eq.0)

  master_cpks_grad_totalsym: if(comm_i_am_master()) then

     ! explicit dependence fit coeff derivatives from 3 center and 2 center
     ! coulombs  integral variation


!!$     DPRINT '[cpks_gradient_totalsym control sum] and [derivatives of a] - (cpks_gradient_totalsym*cpks_Gmat) '
!!$     do i_grad=1,size(cpks,1)
!!$        DPRINT sum(cpks_gradient_totalsym(:,i_grad)*coeff_charge)
!!$        DPRINT sum(matmul(cpks_Gmat,cpks_gradient_totalsym(:,i_grad)))
!!$     enddo
     !--------------------------------------------------------------


     allocate(d1_temp(size(cpks_gradient_totalsym,1), &
          size(cpks_gradient_totalsym,2)),stat=cpksalloc(24))
     MEMLOG(size(d1_temp))
     ASSERT(cpksalloc(24).eq.0)

     d1_temp=cpks_gradient_totalsym  ! both 2 and 3 center contribs
     call dgemm('n','n', n_ch,size(cpks_gradient_totalsym,2),n_ch, &
          1.0_r8_kind,cpks_Gmat,n_ch,d1_temp,n_ch, &
          0.0_r8_kind,cpks_gradient_totalsym,n_ch)

!!$             DPRINT ' d1(cpks_gradient_totalsym) contribs to a_k  derivatives'
!!$             do i_grad=1,size(cpks,1)
!!$                DPRINT sum(cpks_gradient_totalsym(:,i_grad)),i_grad,sum(d1_temp(:,i_grad))
!!$             enddo

     MEMLOG(-size(d1_temp))
     deallocate(d1_temp,stat=cpksalloc(24))
     ASSERT(cpksalloc(24).eq.0)
     cpksalloc(24)=1
  endif master_cpks_grad_totalsym

  cpks_fitcoeff_grads(:,:)=cpks_gradient_totalsym-cpks_gvec/n_spin
  ! cpks_gvec is not split by spins as a contrib to cpks_fitcoeff_grads
  ! which should be the same for both spins ????

  call comm_bcast(cpks_fitcoeff_grads)
  ! broadcasted to be used in calc_cpksQai

  call calc_cpksQai(th_mo_coul,ilower,iupper) ! coul contrib
  if(operations_solvation_effect) call cpks_solv_Qai() !!!!!!!!!!!!!!!AS depends on s1

  call tofile_s1(th_s1)
  !           s1 is not to be used to the end of subroutine and thus can be
  !           presently stored on tape and temporary deallocated
  !           next use to calculate complite first derivatives of hamiltonian

  call reduce_Qai()

  if(comm_i_am_master()) then
     DPRINT 'call fromfile_h1ai'
#if 0
     call fromfile_h1es1_ai()
#endif
     call assembl_Qai() ! Qai+h1ai+s1ai, h1ai and s1ai are not used any more
  endif
#if 0
  DPRINT 'Qai assembled'
  call print_Qai()
  if(size(cpks(1,1,1)%h1ai).gt.0) then
     DPRINT 'h1ai(1,1)', cpks(1,1,1)%h1ai(1,1),cpks(1,1,1)%h1ai(size(cpks(1,1,1)%h1ai,1),1)
  endif
#endif
#if 0
  if(comm_i_am_master()) call dealocate_h1ai_s1ai()
#endif

  MEMLOG(-size(cpks_gvec))
  deallocate (cpks_gvec,stat=cpksalloc(3)) ! allocated in calc_cpks_gvec
  ASSERT(cpksalloc(3).eq.0)
  cpksalloc(3)=1
  FPP_TIMER_STOP(t_cpks_init)
  DPRINT 't_cpks_init Q coul calculated', FPP_TIMER_VALUE(t_cpks_init)

  allocate(QxCOai(ilower:iupper),GxQxCOai(n_ch,size(cpks,1)),stat=cpksalloc(34))
  ASSERT(cpksalloc(34).eq.0)
  MEMLOG(size(QxCOai)+size(GxQxCOai))


#if 0 /* debug tool */
  call tofix_mo()
#endif

  if(comm_i_am_master()) call b0() ! cpks_gamma(0)

#ifdef b_on_file
  call deallocate_Q()
#endif
  call allocate_HBH()

#ifndef b_on_file
  call hole_corrected_b(0)
#else
  call hole_corrected_b()
#endif

  FPP_TIMER_START(t_cpks_init)
  call coul_imp_grads(GxQxCOai,th_mo_coul,cpks_gmat)   ! (1)

  call comm_allreduce(GxQxCOai)
#if 0
  call allocate_AB()

  DPRINT 'allocate_Q_AB_B_HBH done'
  print*,'max_cpks_coul_grad_size  co_ai_size cpks_size cpks_fitmat_size ', &
       max_cpks_coul_grad_size, co_ai_size,cpks_size,cpks_fitmat_size
#endif
  call calc_ab(th_mo_coul) !(1)
  FPP_TIMER_STOP(t_cpks_init)


  i_cpks=0
  j_cpks=0

  ![[== XC call 2 ===================================================
  if( xc_is_on(xc_ANY) )then
     FPP_TIMER_START(t_dervs_xc)
     call cpks_xc_main(22) ! contrib to A matrix, i_cpks=0
     FPP_TIMER_STOP(t_dervs_xc)
     DPRINT 't_dervs_xc i_cpks=0 contrib to A matrix',FPP_TIMER_VALUE(t_dervs_xc)

  endif
  !]]== eof XC call 2 ===============================================

  if(operations_solvation_effect) call cpks_solv_main() !!!!!!!!!!!!!!!AS

  call reduce_ABai() !co_ai

  if(comm_i_am_master()) then
#ifdef b_on_file

     ! print warning
     do i_grad=1,size(cpks,1)
        if(abs(cpks_gamma(0,i_grad)).ge.gamma_eps)  cycle
        DPRINT 'cpks calculations for projection are skeeped'
     enddo

     call eigval_weight_ab()

     call readwriteblocked_startread(trim(tmpfile('cpks_b_0')), cpks_b_th(0))

     ! calculate betta(0,0)

     betta_spin_loop: do i_spn=1,n_spin
        betta_irr:  do i_g=1,n_irrep
           n_ai=cpks3c(i_g,i_spn)%n_ai
           if(n_ai <= 0) cycle
           betta_grads: do i_grad=1,size(cpks,1)
              do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                 call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%HBH(:,k),cpks_b_th(0))
              enddo

              if(abs(cpks_gamma(0,i_grad)).lt.gamma_eps)  then
                 cpks(i_grad,i_g,i_spn)%n_fin_cyc=0
              else ! exist for cpks_gamma(0) .ne. 0
                 cpks_betta(0,0,i_grad)=cpks_betta(0,0,i_grad)+&
                      real(symmetry_data_n_partners(i_g),r8_kind)/cpks_gamma(0,i_grad)* &
                      sum(cpks(i_grad,i_g,i_spn)%ABi*cpks(i_grad,i_g,i_spn)%HBH)
              endif

           enddo betta_grads
        enddo betta_irr
     enddo betta_spin_loop
     call readwriteblocked_stopread(cpks_b_th(0),'keep')

     ! calculate cpks_gamma(1)

     call readwriteblocked_startwrite(trim(tmpfile('cpks_b_1')), cpks_b_th(1))
     bplus1_spin_loop: do i_spn=1,n_spin
        bplus1_irr:  do i_g=1,n_irrep
           n_ai=cpks3c(i_g,i_spn)%n_ai
           if(n_ai.le.0) cycle
           bplus1_grads: do i_grad=1,size(cpks,1)

              cpks(i_grad,i_g,i_spn)%m=0.0_r8_kind
              cpks(i_grad,i_g,i_spn)%m(0,0)=1.0_r8_kind-cpks_betta(0,0,i_grad)

              cpks(i_grad,i_g,i_spn)%ABi=cpks(i_grad,i_g,i_spn)%ABi &
                   -cpks_betta(0,0,i_grad)*cpks(i_grad,i_g,i_spn)%HBH
              cpks_gamma(1,i_grad)=cpks_gamma(1,i_grad)+symmetry_data_n_partners(i_g)*&
                   sum(cpks(i_grad,i_g,i_spn)%ABi**2)

              do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                 call readwriteblocked_write(cpks(i_grad,i_g,i_spn)%ABi(:,k),cpks_b_th(1)) !{2)
              enddo

           enddo bplus1_grads
        enddo bplus1_irr
     enddo bplus1_spin_loop
     call readwriteblocked_stopwrite(cpks_b_th(1))

#else
     betta_spin_loop: do i_spn=1,n_spin
        betta_irr:  do i_g=1,n_irrep
           n_ai=cpks3c(i_g,i_spn)%n_ai
           if(n_ai <= 0) cycle
           betta_grads: do i_grad=1,size(cpks,1)

              if(abs(cpks_gamma(0,i_grad)).lt.gamma_eps)  then
                 cpks(i_grad,i_g,i_spn)%n_fin_cyc=0
              else ! exist for cpks_gamma(0) .ne. 0
                 cpks_betta(0,0,i_grad)=cpks_betta(0,0,i_grad)+&
                      real(symmetry_data_n_partners(i_g),r8_kind)/cpks_gamma(0,i_grad)* &
                      sum(cpks(i_grad,i_g,i_spn)%ABi(:,:)*cpks(i_grad,i_g,i_spn)%B(:,:,0))
              endif
           enddo betta_grads
        enddo betta_irr
     enddo betta_spin_loop

     bplus1_spin_loop: do i_spn=1,n_spin
        bplus1_irr:  do i_g=1,n_irrep
           n_ai=cpks3c(i_g,i_spn)%n_ai
           if(n_ai.le.0) cycle
           bplus1_grads: do i_grad=1,size(cpks,1)

              cpks(i_grad,i_g,i_spn)%m=0.0_r8_kind
              cpks(i_grad,i_g,i_spn)%m(0,0)=1.0_r8_kind-cpks_betta(0,0,i_grad)

              cpks(i_grad,i_g,i_spn)%B(:,:,1)=cpks(i_grad,i_g,i_spn)%ABi(:,:) &
                   -cpks_betta(0,0,i_grad)*cpks(i_grad,i_g,i_spn)%B(:,:,0)

              cpks_gamma(1,i_grad)=cpks_gamma(1,i_grad)+symmetry_data_n_partners(i_g)*&
                   sum(cpks(i_grad,i_g,i_spn)%B(:,:,1)**2)

           enddo bplus1_grads
        enddo bplus1_irr
     enddo bplus1_spin_loop
#endif

  endif ! i_am_master






  fill_cpks_i: do i_cpks=1,n_cpks


     DPRINT i_cpks,'i_cpks_________________________________________________________'
     DPRINT 'HBH sum'
     ! ** main result B-matrix for (i_cpks+1)
     FPP_TIMER_START(t_cpks_init)
#ifdef b_on_file
     call hole_corrected_b()
#else
     call hole_corrected_b(i_cpks)
#endif
     call coul_imp_grads(GxQxCOai,th_mo_coul,cpks_gmat) !(2)

     call comm_allreduce(GxQxCOai)

     call calc_ab(th_mo_coul)  !(2)
     FPP_TIMER_STOP(t_cpks_init)

     ![[== XC call 3 ===================================================
     if( xc_is_on(xc_ANY) )then
        FPP_TIMER_START(t_dervs_xc)
        call cpks_xc_main(23) ! contrib to A, i_cpks
        FPP_TIMER_STOP(t_dervs_xc)
        DPRINT 't_dervs_xc i_cpks contrib to A matrix',FPP_TIMER_VALUE(t_dervs_xc),i_cpks

     endif
     !]]== eof XC call 3 ===============================================

     if(operations_solvation_effect) call cpks_solv_main() !!!!!!!!!!!!!!!AS

     call reduce_ABai()

     fill_cpks_master: if(comm_i_am_master()) then


        fill_gamma_spin: do i_spn=1,n_spin
           fill_gamma_irr: do i_g=1,n_irrep
              n_ai=cpks3c(i_g,i_spn)%n_ai
              if(n_ai.le.0) cycle

              fill_gamma_grads: do i_grad=1,size(cpks,1)
                 do n_occ=1,size(cpks(1,i_g,i_spn)%ABi,1)
                    do n_vir=1,size(cpks(1,i_g,i_spn)%ABi,2)
                       if(abs(eigval_occ(i_g)%m(n_occ,i_spn) &
                            -cpks3c(i_g,i_spn)%eigval_vir_and_holes(n_vir)).gt.0.000000001_r8_kind)  then
                          cpks(i_grad,i_g,i_spn)%ABi(n_occ,n_vir)=cpks(i_grad,i_g,i_spn)%ABi(n_occ,n_vir) &
                               /(eigval_occ(i_g)%m(n_occ,i_spn)-cpks3c(i_g,i_spn)%eigval_vir_and_holes(n_vir))
                       else
                          cpks(i_grad,i_g,i_spn)%ABi(n_occ,n_vir)=0.0_r8_kind
                       endif
                    enddo
                 enddo

              enddo fill_gamma_grads
           enddo fill_gamma_irr
        enddo fill_gamma_spin


        ! fill cpks_betta

#ifndef b_on_file
        fill_betta_spin: do i_spn=1,n_spin
           fill_betta_irr: do i_g=1,n_irrep
              n_ai=cpks3c(i_g,i_spn)%n_ai
              if(cpks3c(i_g,i_spn)%n_ai.le.0) cycle
              fill_betta_grads: do i_grad=1,size(cpks,1)

                 cpks(i_grad,i_g,i_spn)%B(:,:,i_cpks+1)=cpks(i_grad,i_g,i_spn)%ABi(:,:)
                 ! ABi will not be used below exept this loop thus it can be hole occupation
                 ! modified to calc betta factors

                 if(cpks(i_grad,i_g,i_spn)%n_fin_cyc.ge.0) cycle
                 betta_cpks_j: do j_cpks=0,i_cpks
                    call dgemm('t','n', 1,1,n_ai, &
                         real(symmetry_data_n_partners(i_g),r8_kind)/cpks_gamma(j_cpks,i_grad), &
                         cpks(i_grad,i_g,i_spn)%ABi(:,:),n_ai, &
                         cpks(i_grad,i_g,i_spn)%B(:,:,j_cpks),n_ai, &
                         1.0_r8_kind,cpks_betta(j_cpks,i_cpks,i_grad),1)
                 enddo betta_cpks_j

              enddo fill_betta_grads
           enddo fill_betta_irr
        enddo fill_betta_spin

#else


        betta_cpks_j: do j_cpks=0,i_cpks
           write(b_fn,'(i2)') j_cpks
           call readwriteblocked_startread(trim(tmpfile('cpks_b_'//adjustl(b_fn))), cpks_b_th(j_cpks))
           fill_betta_spin: do i_spn=1,n_spin
              fill_betta_irr: do i_g=1,n_irrep
                 n_ai=cpks3c(i_g,i_spn)%n_ai
                 if(cpks3c(i_g,i_spn)%n_ai.le.0) cycle
                 fill_betta_grads: do i_grad=1,size(cpks,1)
                    do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                       call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%HBH(:,k),cpks_b_th(j_cpks))
                    enddo
                    if(cpks(i_grad,i_g,i_spn)%n_fin_cyc.ge.0) cycle !no beta calc if converged by i_cpks
                    if(abs(cpks_gamma(i_cpks,i_grad)).le.gamma_eps) cycle

                    call dgemm('t','n', 1,1,n_ai, &
                         real(symmetry_data_n_partners(i_g),r8_kind)/cpks_gamma(j_cpks,i_grad), &
                         cpks(i_grad,i_g,i_spn)%ABi(:,:),n_ai, &
                         cpks(i_grad,i_g,i_spn)%HBH(:,:),n_ai, &
                         1.0_r8_kind,cpks_betta(j_cpks,i_cpks,i_grad),1)

                 enddo fill_betta_grads
              enddo fill_betta_irr
           enddo fill_betta_spin
           call readwriteblocked_stopread(cpks_b_th(j_cpks),'keep')
        enddo betta_cpks_j

#endif /* b_on_file */

        do i_grad=1,size(cpks,1)
           do i_g=1,n_irrep
              do i_spn=1,n_spin
                 if(cpks(i_grad,i_g,i_spn)%n_fin_cyc.ge.0) cycle
                 cpks(i_grad,i_g,i_spn)%m(i_cpks,i_cpks)=1.0_r8_kind-cpks_betta(i_cpks,i_cpks,i_grad)
                 cpks(i_grad,i_g,i_spn)%m(i_cpks,i_cpks-1)=-1.0_r8_kind
                 do j_cpks=0,i_cpks-1
                    cpks(i_grad,i_g,i_spn)%m(j_cpks,i_cpks)=-cpks_betta(j_cpks,i_cpks,i_grad)
                 enddo
              enddo
           enddo
        enddo

        call b_plus(i_cpks)

        ! set dependent on gamma(i_cpks) mask

        cpks_ser_limited=.false.
        do i_grad=1,size(cpks,1)
           if(abs(cpks_gamma(i_cpks,i_grad)).ge.gamma_eps ) cycle
           do i_spn=1,n_spin
              do i_g=1,n_irrep
                 if(cpks(i_grad,i_g,i_spn)%n_fin_cyc.ge.0) cycle
                 cpks(i_grad,i_g,i_spn)%n_fin_cyc=i_cpks
                 cpks_ser_limited=.true.
              enddo
           enddo
        enddo

        if(cpks_ser_limited) then
           DPRINT 'cpks_ser_limited i_cpks', i_cpks
        endif

     endif fill_cpks_master
  enddo fill_cpks_i
  !-------------------------------------------------------------------------------------------


  if(comm_i_am_master()) then
     call solve_cpks() !!! result: B(NCPKS+1)
  endif


  FPP_TIMER_START(t_cpks_init)
#ifndef b_on_file
  call hole_corrected_b(n_cpks+1)
#else
  call hole_corrected_b()
#endif
  call coul_imp_grads(GxQxCOai,th_mo_coul,cpks_gmat)  !(3)

  call comm_allreduce(GxQxCOai)

  if(comm_i_am_master()) then
     MEMLOG(-size(cpks_gamma)-size(cpks_betta))
     deallocate(cpks_gamma, cpks_betta, stat=cpksalloc(118))
     ASSERT(cpksalloc(118).eq.0)
     cpksalloc(118)=1
  endif


  if(comm_i_am_master())  then
     cpks_fitcoeff_grads=cpks_fitcoeff_grads+GxQxCOai
  endif

  call comm_bcast(cpks_fitcoeff_grads)
  !   ! brodcated to be used in  calc_cpks_h1imp

  call calc_ab(th_mo_coul,'delete')   !(3)  here co_ai is last time used
#ifndef mo_coul_onfile
  call deallocate_co_ai()
#endif

  FPP_TIMER_STOP(t_cpks_init)

  i_cpks=n_cpks+1 ! now check obtained result
  ! and calculate density matrix contrib
  ! to fitcoeff derivatives


!!$             DPRINT '%h1 func before xc imp calcs------',sum(abs(cpks(1,1,1)%h1(:,:)))


  ![[== XC call 4 ===================================================
  if( xc_is_on(xc_ANY) )then
     call say("implicite xc dervs ")
     FPP_TIMER_START(t_dervs_xc)
     comm_dervs_xc=.true. ! set on to communicate dervs_xc
     call cpks_xc_main(24) ! contrib to A matrix, i_cpks=n_cpks+1 (last of 3 calls)
     ! result: dependent on B h1 contrib
     comm_dervs_xc=.false. ! set off for next postscf call
     FPP_TIMER_STOP(t_dervs_xc)
     DPRINT 't_dervs_xc n_cpks+1 contrib to A matrix',FPP_TIMER_VALUE(t_dervs_xc),i_cpks
  endif
  !]]== eof XC call 4 ===============================================

  if(operations_solvation_effect) call cpks_solv_main() !!!!!!!!!!!!!!!AS

  call reduce_ABai()

  call cpks_convergence()

#ifndef b_on_file
  if(comm_i_am_master())  call deallocate_B()
  call deallocate_Q()
#endif
  call deallocate_AB()

  DPRINT 't_cpks_init coul_imp_grads', FPP_TIMER_VALUE(t_cpks_init)
  FPP_TIMER_START(t_cpks_init)

  !     s1 is used here for h1 check  thus is to be swapt to memory
  !     next use is wor calculation of w1
  call fromfile_s1(th_s1)

  if(operations_solvation_effect) call cpks_solv_H1() !!!!!!!!!!!!!!!AS
  if(operations_solvation_effect) call cpks_solv_H1_1() !!!!!!!!!!!!!!!AS
  call calc_cpks_h1imp(ilower, iupper) !!! s1 undependent

  FPP_TIMER_STOP(t_cpks_init)
  DPRINT 't_cpks_init calc_cpks_h1imp', FPP_TIMER_VALUE(t_cpks_init)

  if(comm_parallel()) call collect_h1()


#if 0  /* convenient check on h1 */
  h1_irr: do i_g=1,n_irrep
     n_occ=size(cpks(1,i_g,1)%h1,1)
     if(n_occ.gt.0) then
        print*
        print*,'%h1 func+imp',cpks(1,i_g,1)%h1(1,1),cpks(1,i_g,1)%h1(n_occ,n_occ)
        do n_occ=1,size(cpks(1,i_g,1)%h1,1)
           print*,cpks(1,i_g,1)%h1(n_occ,:)
           print*,'%h1-%s1*eig_p check',n_occ,cpks(1,i_g,1)%h1(n_occ,n_occ) &
                -eigval(i_g)%m(n_occ,1)*cpks(1,i_g,1)%s1(n_occ,n_occ), &
                eigval(i_g)%m(n_occ,1)
        enddo
     endif !n_occ.gt.0
  enddo h1_irr
#endif
  call hole_corrected_h1()

  !         if(operations_solvation_effect) then
  !          call cpks_solv_ham()
  !          print*,'solv h0',sum(cpks(1,1,1)%Qai),sum(cpks(2,1,1)%Qai)
  !          call calc_solv_ham_mo('ai')
  !         endif


  if(allocated(cpks_gmat))  then
     MEMLOG(-size(cpks_gmat))
     deallocate(cpks_gmat,stat=cpksalloc(4))
     ASSERT(cpksalloc(4).eq.0)
     cpksalloc(4)=1
  endif

  MEMLOG(-size(QxCOai)-size(GxQxCOai))
  deallocate(QxCOai,GxQxCOai,stat=cpksalloc(34))
  ASSERT(cpksalloc(34).eq.0)
  cpksalloc(34)=1



  FPP_TIMER_STOP(t_cpks_const)

#ifdef FPP_TIMERS
  print*,'t_s1_imp_grarho  =',FPP_TIMER_VALUE(t_s1_imp_grarho)
  print*,'t_s1_imp_dervsrho=',FPP_TIMER_VALUE(t_s1_imp_dervsrho)
  print*,'t_xc_qh_expl     =',FPP_TIMER_VALUE(t_xc_qh_expl)
  print*,'t_h1q_dvdrho     =',FPP_TIMER_VALUE(t_h1q_dvdrho)
  print*,'t_dervs_helps     =',FPP_TIMER_VALUE(t_dervs_helps)
  print*,'t_cpks_init      =',FPP_TIMER_VALUE(t_cpks_init)
  print*,'t_xc_bimp_grarho gragamma =',FPP_TIMER_VALUE(t_xc_bimp_grarho)
  print*,'t_xc_ab          =',FPP_TIMER_VALUE(t_xc_ab)
  print*,'t_dervs_imp      =',FPP_TIMER_VALUE(t_dervs_imp)
  print*,'t_calc_imp_dervs=', FPP_TIMER_VALUE(t_calc_imp_dervs)
  print*,'t_grid_unique tot   =',FPP_TIMER_VALUE(t_grid_uni)
  print*,'t_orbital_calculate = ',FPP_TIMER_VALUE(t_orbital_calc)
  print*,'t_density_calculate=', FPP_TIMER_VALUE(t_density_calc)
  print*,'t_dervs_xc tot   =',FPP_TIMER_VALUE(t_dervs_xc)
  print*,'t_cpks_const     =',FPP_TIMER_VALUE(t_cpks_const)
  print*, 'cpks_constructs mem controls'
#endif
  MEMSET(0)
  call grid_close(.false.)

contains
  subroutine calc_cpks_gmat()
    real(kind=r8_kind),allocatable:: cpks_inv_chmat(:,:)
    real(kind=r8_kind),dimension(:  ),allocatable :: inv_linsys_mat ! for cpks equations
    real(kind=r8_kind),dimension(:  ),allocatable :: linsys_mat,linsys_dual
    call get_n_elec(target_charge)

    ! Linear solver
    ! note that size of linsys_dual is n_ch + 1 !
    allocate(linsys_mat((n_ch*(n_ch+1))/2),linsys_dual(n_ch+1),STAT=cpksalloc(5))
    MEMLOG(size(linsys_mat)+size(linsys_dual))
    ASSERT(cpksalloc(5).eq.0)
    linsys_mat = mat_charge

    call decompose_linsys(n_ch,linsys_mat,charge_norm,linsys_dual)

    allocate(inv_linsys_mat(n_ch*(n_ch+1)/2),stat=cpksalloc(116))
    ASSERT(cpksalloc(116).eq.0)
    MEMLOG(size(inv_linsys_mat))


    allocate(cpks_inv_chmat(n_ch,n_ch),stat=cpksalloc(19))
    ASSERT(cpksalloc(19).eq.0)
    MEMLOG(size(cpks_inv_chmat))

    inv_linsys_mat=linsys_mat
    call gmat_linsys(n_ch,inv_linsys_mat,charge_norm, &
         linsys_dual, cpks_gmat,cpks_inv_chmat)

    MEMLOG(-size(inv_linsys_mat)-size(cpks_inv_chmat))
    deallocate(inv_linsys_mat,cpks_inv_chmat,stat=cpksalloc(19))
    ASSERT(cpksalloc(19).eq.0)
    cpksalloc(19)=1
    cpksalloc(116)=1

    DPRINT ' cpks_g4constructs: cpks_gmat is calculated' ,sum(cpks_gmat)

    MEMLOG(-size(linsys_dual)-size(linsys_mat))
    deallocate(linsys_dual,linsys_mat,STAT=cpksalloc(5))
    ASSERT(cpksalloc(5).eq.0)
    cpksalloc(5)=1
  end subroutine calc_cpks_gmat

#ifdef mo_coul_onfile
#ifndef no_co_ai
  subroutine write_mo_coul()
    integer(kind=i4_kind):: i_spn,i_g,k,i_occ
    call readwriteblocked_startwrite(trim(tmpfile('mo_coul.dat')), th_mo_coul)
    do i_spn=1,size(cpks3c,2)
       do i_g=1,size(cpks3c,1)
          do k=1,size(cpks3c(i_g,i_spn)%co_ai,1)
             do i_occ=1,size(cpks3c(i_g,i_spn)%co_ai,2)
                call readwriteblocked_write(cpks3c(i_g,i_spn)%co_ai(k,i_occ,:),th_mo_coul)
             enddo
          enddo
       enddo
    enddo
    call readwriteblocked_stopwrite(th_mo_coul)
  end subroutine write_mo_coul
#endif
#endif

  subroutine tofile_h1es1_ai()
    integer(kind=i4_kind):: i_spn,i_g,i_occ,i_grad,stat,j_vir
    call readwriteblocked_startwrite(trim(tmpfile('h1ai.dat')), th_h1ai)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)

          do i_grad=1,size(cpks,1)

             if(size(cpks(i_grad,i_g,i_spn)%h1ai).ne.0) then
                do i_occ=1,size(cpks(i_grad,i_g,i_spn)%h1ai,1)
                   do j_vir=1,size(cpks(i_grad,i_g,i_spn)%h1ai,2)
                      cpks(i_grad,i_g,i_spn)%h1ai(i_occ,j_vir) = cpks(i_grad,i_g,i_spn)%h1ai(i_occ,j_vir) &
                           -cpks(i_grad,i_g,i_spn)%s1ai(i_occ,j_vir)*eigval(i_g)%m(i_occ,i_spn)
                   enddo
                   call readwriteblocked_write(cpks(i_grad,i_g,i_spn)%h1ai(i_occ,:),th_h1ai)
                enddo
             endif

             deallocate(cpks(i_grad,i_g,i_spn)%h1ai,cpks(i_grad,i_g,i_spn)%s1ai,stat=stat)
             ASSERT(stat.eq.0)
             cpksalloc(167)=1 ! h1ai
             cpksalloc(162)=1
          enddo

       enddo
    enddo
    call readwriteblocked_stopwrite(th_h1ai)
  end subroutine tofile_h1es1_ai

  subroutine fromfile_h1es1_ai()
    integer(kind=i4_kind):: i_spn, i_g, i_occ, i_grad, stat, occ_dim, vir_dim
    call readwriteblocked_startread(trim(tmpfile('h1ai.dat')), th_h1ai)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          occ_dim= size(cpks(1,i_g,i_spn)%Qai,1)
          vir_dim= size(cpks(1,i_g,i_spn)%Qai,2)
          do i_grad=1,size(cpks,1)
             allocate(cpks(i_grad,i_g,i_spn)%h1ai(occ_dim,vir_dim),cpks(i_grad,i_g,i_spn)%s1ai(0,0),stat=stat)
             ASSERT(stat.eq.0)
             do i_occ=1,occ_dim
                call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%h1ai(i_occ,:),th_h1ai)
             enddo
          enddo
       enddo
    enddo
    call readwriteblocked_stopread(th_h1ai,'delete')
  end subroutine fromfile_h1es1_ai

  subroutine tofile_h1(th_h1)
    type(readwriteblocked_tapehandle),intent(inout) :: th_h1
    integer(kind=i4_kind):: i_spn,i_g,i_occ,i_grad,stat
    call readwriteblocked_startwrite(trim(tmpfile('h1.dat')), th_h1)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             do i_occ=1,size(cpks(i_grad,i_g,i_spn)%h1,1)
                call readwriteblocked_write(cpks(i_grad,i_g,i_spn)%h1(i_occ,:),th_h1)
             enddo
          enddo
       enddo
    enddo
    call readwriteblocked_stopwrite(th_h1)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             deallocate(cpks(i_grad,i_g,i_spn)%h1,stat=stat)
             ASSERT(stat.eq.0)
          enddo
       enddo
    enddo
  end subroutine tofile_h1

  subroutine tofile_s1(th_s1)
    type(readwriteblocked_tapehandle),intent(inout) :: th_s1
    integer(kind=i4_kind):: i_spn,i_g,i_occ,i_grad,stat
    call readwriteblocked_startwrite(trim(tmpfile('s1.dat')), th_s1)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             do i_occ=1,size(cpks(i_grad,i_g,i_spn)%s1,1)
                call readwriteblocked_write(cpks(i_grad,i_g,i_spn)%s1(i_occ,:),th_s1)
             enddo
          enddo
       enddo
    enddo
    call readwriteblocked_stopwrite(th_s1)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             deallocate(cpks(i_grad,i_g,i_spn)%s1,stat=stat)
             ASSERT(stat.eq.0)
          enddo
       enddo
    enddo
  end subroutine tofile_s1

  subroutine fromfile_s1(th_s1)
    type(readwriteblocked_tapehandle),intent(inout) :: th_s1
    integer(kind=i4_kind):: i_spn,i_g,i_occ,i_grad,n_occ,stat
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          n_occ=cpks(1,i_g,i_spn)%occ_dim
          do i_grad=1,size(cpks,1)
             allocate(cpks(i_grad,i_g,i_spn)%s1(n_occ,n_occ),stat=stat)
             ASSERT(stat.eq.0)
          enddo
       enddo
    enddo
    call readwriteblocked_startread(trim(tmpfile('s1.dat')), th_s1)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             do i_occ=1,size(cpks(i_grad,i_g,i_spn)%s1,1)
                call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%s1(i_occ,:),th_s1)
             enddo
          enddo
       enddo
    enddo
    call readwriteblocked_stopread(th_s1,'delete')
  end subroutine fromfile_s1

  subroutine fromfile_h1(th_h1)
    type(readwriteblocked_tapehandle),intent(inout) :: th_h1
    integer(kind=i4_kind):: i_spn,i_g,i_occ,i_grad,n_occ,stat
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          n_occ=cpks(1,i_g,i_spn)%occ_dim
          do i_grad=1,size(cpks,1)
             allocate(cpks(i_grad,i_g,i_spn)%h1(n_occ,n_occ),stat=stat)
             ASSERT(stat.eq.0)
          enddo
       enddo
    enddo
    call readwriteblocked_startread(trim(tmpfile('h1.dat')), th_h1)
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             do i_occ=1,size(cpks(i_grad,i_g,i_spn)%h1,1)
                call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%h1(i_occ,:),th_h1)
             enddo
          enddo
       enddo
    enddo
    call readwriteblocked_stopread(th_h1,'delete')
  end subroutine fromfile_h1

#ifndef no_co_ai
  subroutine deallocate_co_ai()
    integer(i4_kind):: i_g,i_spn
    do i_g=1,size(cpks,2)
       do i_spn=1,size(cpks,3)
          MEMLOG(-size(cpks3c(i_g,i_spn)%co_ai))
          deallocate(cpks3c(i_g,i_spn)%co_ai,stat=cpksalloc(8))
          ASSERT(cpksalloc(8).eq.0)
          cpksalloc(8)=1
       enddo
    enddo
  end subroutine deallocate_co_ai
#endif

  subroutine cpks_convergence()
    integer(i4_kind):: stat
    real(r8_kind) ,allocatable:: Qai(:,:)
    real(r8_kind):: U_AB_B0
#ifdef b_on_file
    if(comm_i_am_master()) then
       call readwriteblocked_startread(trim(tmpfile('cpks_b_0')), cpks_b_th(0))
       write(b_fn,'(i2)') n_cpks+1
       call readwriteblocked_startread(trim(tmpfile('cpks_b_'//adjustl(b_fn))), cpks_b_th(n_cpks+1))
    endif
#endif
    spin_use_cpks: do i_spn=1,n_spin
       irr_use_cpks: do i_g=1,n_irrep

!!$          if(comm_i_am_master()) then
!!$             DPRINT 'AB + xc before energy weighting fin'
!!$             DPRINT sum(abs(cpks(1,i_g,i_spn)%ABi)),i_g,'i_g'
!!$             DPRINT '     '
!!$             do n_occ=1,size(cpks(1,i_g,i_spn)%ABi,1)
!!$                DPRINT cpks(1,i_g,i_spn)%ABi(n_occ,:)
!!$             enddo
!!$          endif

          n_ai=cpks3c(i_g,i_spn)%n_ai
          if(n_ai.gt.0) then
             gr6: do i_grad=1,size(cpks,1)
                do n_occ=1,size(cpks(i_grad,i_g,i_spn)%ABi,1)
                   do n_vir=1,size(cpks(i_grad,i_g,i_spn)%ABi,2)

                      if(comm_i_am_master().and. &
                           abs(eigval_occ(i_g)%m(n_occ,i_spn)- &
                           cpks3c(i_g,i_spn)%eigval_vir_and_holes(n_vir)).gt.0.000000001_r8_kind ) then
                         cpks(i_grad,i_g,i_spn)%ABi(n_occ,n_vir)=cpks(i_grad,i_g,i_spn)%ABi(n_occ,n_vir) &
                              /(eigval_occ(i_g)%m(n_occ,i_spn)-cpks3c(i_g,i_spn)%eigval_vir_and_holes(n_vir))
                      else
                         cpks(i_grad,i_g,i_spn)%ABi(n_occ,n_vir)=0.0_r8_kind
                      endif

                   enddo
                enddo

                if(comm_i_am_master()) then
#ifndef b_on_file
!!$             DPRINT                                   sum(abs(cpks(i_grad,i_g,i_spn)%B(:,:,n_cpks+1)))
!!$             DPRINT                                   sum(abs(cpks(i_grad,i_g,i_spn)%B(:,:,0)))

                   ! here Qai holds U-AB-BO to be equal to zero
                   cpks(i_grad,i_g,i_spn)%Qai= &
                        cpks(i_grad,i_g,i_spn)%B(:,:,n_cpks+1)-cpks(i_grad,i_g,i_spn)%ABi &
                        -cpks(i_grad,i_g,i_spn)%B(:,:,0)

                   U_AB_B0=sum(cpks(i_grad,i_g,i_spn)%Qai**2)
!!$    print*,  'B',  sum(cpks(i_grad,i_g,i_spn)%B(:,:,n_cpks+1)),cpks(i_grad,i_g,i_spn)%B(1,1,n_cpks+1)
!!$    DPRINT  'abs B', sum(abs(cpks(i_grad,i_g,i_spn)%B(:,:,n_cpks+1)))
#else
                   allocate(Qai(size(cpks(i_grad,i_g,i_spn)%ABi,1),size(cpks(i_grad,i_g,i_spn)%ABi,2)),stat=stat)
                   ASSERT(stat.eq.0)

                   ! cpks(i_grad,i_g,i_spn)%Qai=-cpks(i_grad,i_g,i_spn)%ABi
                   Qai=-cpks(i_grad,i_g,i_spn)%ABi

                   do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                      call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%ABi(:,k),cpks_b_th(n_cpks+1))
                   enddo
                   ! cpks(i_grad,i_g,i_spn)%Qai=cpks(i_grad,i_g,i_spn)%Qai+cpks(i_grad,i_g,i_spn)%ABi
                   Qai=Qai+cpks(i_grad,i_g,i_spn)%ABi
                   do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                      call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%ABi(:,k),cpks_b_th(0))
                   enddo
                   ! cpks(i_grad,i_g,i_spn)%Qai=cpks(i_grad,i_g,i_spn)%Qai-cpks(i_grad,i_g,i_spn)%ABi
                   ! U_AB_B0=sum(cpks(i_grad,i_g,i_spn)%Qai**2)
                   Qai=Qai-cpks(i_grad,i_g,i_spn)%ABi
                   U_AB_B0=sum(Qai**2)
                   deallocate(Qai,stat=stat)
                   ASSERT(stat.eq.0)
#endif
                   print*, 'U-AB-BO and fin cpks sol',i_g,i_grad,U_AB_B0
                endif ! i_am_master
             enddo gr6
          endif


       enddo irr_use_cpks
    enddo spin_use_cpks
#ifdef b_on_file
    if(comm_i_am_master()) then
       call readwriteblocked_stopread(cpks_b_th(0),'delete')
       call readwriteblocked_stopread(cpks_b_th(n_cpks+1),'delete')
    endif
#endif
  end subroutine cpks_convergence

  subroutine b_plus(i_cpks)

    ! calculates ABi,  m and gamma(i_cpks+1)
    ! sets n_fin_cyc to particular value if converged by gamma

    integer(i4_kind), intent(in):: i_cpks
#ifdef b_on_file
    cpks_j: do j_cpks=0,i_cpks

       write(b_fn,'(i2)') j_cpks
       call readwriteblocked_startread(trim(tmpfile('cpks_b_'//adjustl(b_fn))), cpks_b_th(j_cpks))

       fill_bplus1_spin: do i_spn=1,n_spin
          fill_bplus1_irr: do i_g=1,n_irrep
             n_ai=cpks3c(i_g,i_spn)%n_ai
             if(n_ai.le.0) cycle
             fill_bplus1_grads: do i_grad=1,size(cpks,1)

                do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                   call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%HBH(:,k),cpks_b_th(j_cpks))
                enddo

                if(cpks(i_grad,i_g,i_spn)%n_fin_cyc.ge.0) cycle
                cpks(i_grad,i_g,i_spn)%ABi=cpks(i_grad,i_g,i_spn)%ABi &
                     -cpks_betta(j_cpks,i_cpks,i_grad)*cpks(i_grad,i_g,i_spn)%HBH

             enddo fill_bplus1_grads

          enddo fill_bplus1_irr
       enddo fill_bplus1_spin
       call readwriteblocked_stopread(cpks_b_th(j_cpks),'keep')
    enddo cpks_j


    !print*,'write 1+i_cpks',1+i_cpks
    write(b_fn,'(i2)') 1+i_cpks
    call readwriteblocked_startwrite(trim(tmpfile('cpks_b_'//adjustl(b_fn))), cpks_b_th(1+i_cpks))
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                call readwriteblocked_write(cpks(i_grad,i_g,i_spn)%ABi(:,k),cpks_b_th(1+i_cpks))
             enddo
             if(cpks(i_grad,i_g,i_spn)%n_fin_cyc.ge.0) cycle
             cpks_gamma(1+i_cpks,i_grad)=cpks_gamma(1+i_cpks,i_grad) &
                  +symmetry_data_n_partners(i_g)*sum(cpks(i_grad,i_g,i_spn)%ABi**2)
          enddo
       enddo
    enddo
    call readwriteblocked_stopwrite(cpks_b_th(1+i_cpks))
#else
    cpks_j: do j_cpks=0,i_cpks


       fill_bplus1_spin: do i_spn=1,n_spin
          fill_bplus1_irr: do i_g=1,n_irrep
             n_ai=cpks3c(i_g,i_spn)%n_ai
             fill_bplus1_n_ai: if(n_ai.le.0) cycle
             fill_bplus1_grads: do i_grad=1,size(cpks,1)

                if(cpks(i_grad,i_g,i_spn)%n_fin_cyc.ge.0) cycle
                cpks(i_grad,i_g,i_spn)%B(:,:,i_cpks+1)=cpks(i_grad,i_g,i_spn)%B(:,:,i_cpks+1) &
                     -cpks_betta(j_cpks,i_cpks,i_grad)*cpks(i_grad,i_g,i_spn)%B(:,:,j_cpks)
             enddo fill_bplus1_grads

          enddo fill_bplus1_irr
       enddo fill_bplus1_spin
    enddo cpks_j

    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          do i_grad=1,size(cpks,1)
             cpks_gamma(1+i_cpks,i_grad)=cpks_gamma(1+i_cpks,i_grad) &
                  +symmetry_data_n_partners(i_g)*sum(cpks(i_grad,i_g,i_spn)%B(:,:,1+i_cpks)**2)
          enddo
       enddo
    enddo
#endif

  end subroutine b_plus

  subroutine allocate_HBH()
    integer(i4_kind):: i_grad,i_ir,i_spn,occ_dim,vir_dim
    do i_grad=1,size(cpks,1)
       do i_ir=1,size(cpks,2)
          do i_spn=1,size(cpks,3)
             occ_dim=size(cpks(i_grad,i_ir,i_spn)%h1ai,1)
             vir_dim=size(cpks(i_grad,i_ir,i_spn)%h1ai,2)
             allocate( cpks(i_grad,i_ir,i_spn)%HBH(occ_dim,vir_dim),stat=cpksalloc(166))
             ASSERT(cpksalloc(166).eq.0) ! HBH
             MEMLOG(size(cpks(i_grad,i_ir,i_spn)%HBH))
             cpks_size=cpks_size+size(cpks(i_grad,i_ir,i_spn)%HBH)
          enddo
       enddo
    enddo
  end subroutine allocate_HBH

  subroutine allocate_AB()
    integer(i4_kind):: i_grad,i_ir,i_spn,occ_dim,vir_dim
    do i_grad=1,size(cpks,1)
       do i_ir=1,size(cpks,2)
          do i_spn=1,size(cpks,3)
             occ_dim=cpks3c(i_ir,i_spn)%occ_dim
             vir_dim=cpks3c(i_ir,i_spn)%vir_dim
             allocate( cpks(i_grad,i_ir,i_spn)%ABi(occ_dim,vir_dim), &
                  stat=cpksalloc(165))
             ASSERT(cpksalloc(165).eq.0) ! ABi
             MEMLOG(size(cpks(i_grad,i_ir,i_spn)%ABi))
             cpks_size=cpks_size+size(cpks(i_grad,i_ir,i_spn)%ABi)
          enddo
       enddo
    enddo
  end subroutine allocate_AB

  subroutine allocate_Q_B()
    integer(i4_kind):: i_grad,i_ir,i_spn,occ_dim,vir_dim
    do i_grad=1,size(cpks,1)
       do i_ir=1,size(cpks,2)
          do i_spn=1,size(cpks,3)
             occ_dim=cpks3c(i_ir,i_spn)%occ_dim
             vir_dim=cpks3c(i_ir,i_spn)%vir_dim
             allocate( cpks(i_grad,i_ir,i_spn)%Qai(occ_dim,vir_dim), &
                  stat=cpksalloc(169))
             ASSERT(cpksalloc(169).eq.0)
             MEMLOG(size(cpks(i_grad,i_ir,i_spn)%Qai))
             cpks_size=cpks_size+size(cpks(i_grad,i_ir,i_spn)%Qai)
             cpks(i_grad,i_ir,i_spn)%Qai=0.0_r8_kind
#ifndef b_on_file
             if(comm_i_am_master()) then
                allocate( cpks(i_grad,i_ir,i_spn)%B(occ_dim,vir_dim,0:n_cpks+1),   &
                     stat=cpksalloc(163))    ! B
                ASSERT(cpksalloc(163).eq.0)
                cpks(i_grad,i_ir,i_spn)%B=0.0_r8_kind
                MEMLOG(size(cpks(i_grad,i_ir,i_spn)%B))
                cpks_size=cpks_size+size(cpks(i_grad,i_ir,i_spn)%B)
                call say(" cpks B matrix allocated ")
             endif
#else
             if(i_grad.eq.1.and.i_ir.eq.1.and.i_spn.eq.1) &
                  call say(" cpks B matrix stored to file")
#endif
          enddo
       enddo
    enddo
  end subroutine allocate_Q_B

#ifndef no_co_ai
  subroutine allocate_co_ai()
    integer(i4_kind):: i_g,i_spn,n_occ,n_vir
    co_ai_size=0
    do i_g=1,n_irrep
       do i_spn=1,n_spin
#if 0
          n_occ=size(cpks(1,i_g,i_spn)%h1ai,1)
          n_vir=size(cpks(1,i_g,i_spn)%h1ai,2)
#else
          n_occ=cpks3c(i_g,i_spn)%occ_dim
          n_vir=cpks3c(i_g,i_spn)%vir_dim
#endif
          allocate(cpks3c(i_g,i_spn)%co_ai(iupper-ilower+1,n_occ,n_vir),stat=cpksalloc(8))
          ASSERT(cpksalloc(8).eq.0)
          MEMLOG(size(cpks3c(i_g,i_spn)%co_ai)) ! to be deallocated in this program
          co_ai_size=co_ai_size+size(cpks3c(i_g,i_spn)%co_ai)
          cpks3c(i_g,i_spn)%co_ai=0.0_r8_kind
       enddo
    enddo
  end subroutine allocate_co_ai
#endif

  subroutine print_Qai()
    integer(kind=i4_kind):: i
    do i_g=1,n_irrep
       do i_grad=1,size(cpks,1)
          if(n_spin.eq.2) then
             DPRINT 'sum Q',sum(cpks(i_grad,i_g,1)%Qai),sum(cpks(i_grad,i_g,2)%Qai) ,i_grad,i_g
             if(i_grad.eq.1.and.size(cpks(i_grad,i_g,1)%Qai).gt.0) then
                do i=1,size(cpks(i_grad,i_g,1)%Qai,2)
                   DPRINT cpks(i_grad,i_g,1)%Qai(:,i),1
                enddo
             endif
             if(i_grad.eq.1.and.size(cpks(i_grad,i_g,2)%Qai).gt.0) then
                do i=1,size(cpks(i_grad,i_g,2)%Qai,2)
                   DPRINT cpks(i_grad,i_g,2)%Qai(:,i),2
                enddo
             endif
          else
             DPRINT 'sum Q',sum(cpks(i_grad,i_g,1)%Qai),i_grad,i_g
             if(i_grad.eq.1.and.size(cpks(i_grad,i_g,1)%Qai).gt.0) then
                do i=1,size(cpks(i_grad,i_g,1)%Qai,2)
                   DPRINT cpks(i_grad,i_g,1)%Qai(:,i),1
                enddo
             endif
          endif
       enddo
    enddo
  end subroutine print_Qai

  subroutine deallocate_Q()
    do i_grad=1,size(cpks,1)
       do i_g=1,size(cpks,2)
          do i_spn=1,size(cpks,3)
             MEMLOG(-size(cpks(i_grad,i_g,i_spn)%Qai))
             deallocate(cpks(i_grad,i_g,i_spn)%Qai,stat=cpksalloc(169))
             ASSERT(cpksalloc(169).eq.0)
             cpksalloc(169)=1
          enddo
       enddo
    enddo
  end subroutine deallocate_Q

  subroutine deallocate_AB()
    do i_grad=1,size(cpks,1)
       do i_g=1,size(cpks,2)
          do i_spn=1,size(cpks,3)
             MEMLOG(-size(cpks(i_grad,i_g,i_spn)%Abi))
             deallocate(cpks(i_grad,i_g,i_spn)%Abi,stat=cpksalloc(165))
             ASSERT(cpksalloc(165).eq.0)
             cpksalloc(165)=1
          enddo
       enddo
    enddo
  end subroutine deallocate_AB

  subroutine deallocate_B()
    do i_grad=1,size(cpks,1)
       do i_g=1,size(cpks,2)
          do i_spn=1,size(cpks,3)
             MEMLOG(-size(cpks(i_grad,i_g,i_spn)%B))
             deallocate(cpks(i_grad,i_g,i_spn)%B,stat=cpksalloc(163))
             ASSERT(cpksalloc(163).eq.0)
             cpksalloc(163)=1
          enddo
       enddo
    enddo
  end subroutine deallocate_B

  subroutine get_cpks_dims()
    do i_g=1,size(cpks,2)
       do i_spn=1,size(cpks,3)
          cpks3c(i_g,i_spn)%occ_dim=size(cpks(1,i_g,i_spn)%h1ai,1)
          cpks3c(i_g,i_spn)%vir_dim=size(cpks(1,i_g,i_spn)%h1ai,2)
       enddo
    enddo
  end subroutine get_cpks_dims

  subroutine dealocate_h1ai_s1ai()
    do i_grad=1,size(cpks,1)
       do i_g=1,size(cpks,2)
          do i_spn=1,size(cpks,3)
             MEMLOG(-size(cpks(i_grad,i_g,i_spn)%s1ai))
#if 0 /* immediate h1ai deallocate */
             deallocate(cpks(i_grad,i_g,i_spn)%s1ai,stat=cpksalloc(162))
#else
             MEMLOG(-size(cpks(i_grad,i_g,i_spn)%h1ai))
             deallocate(cpks(i_grad,i_g,i_spn)%h1ai,cpks(i_grad,i_g,i_spn)%s1ai,stat=cpksalloc(162))
#endif
             ASSERT(cpksalloc(162).eq.0)
             cpksalloc(162)=1
          enddo
       enddo
    enddo
  end subroutine dealocate_h1ai_s1ai

#if 0
  subroutine tofix_mo()
    integer(kind=i4_kind):: i_g,n_occ,fxmo_un
    fxmo_un=openget_iounit(file=trim(input_dir)//'/ocmo.dat',&
         form='FORMATTED',status='unknown')

    do i_g=1,n_irrep
       print*,' MO coefficients'
       do n_occ=1,size(eigval_occ(i_g)%m,1)
          write(fxmo_un,*) eigvec_occ(i_g)%m(:,n_occ,1)
          print*,eigvec_occ(i_g)%m(:,n_occ,1)
       enddo
    enddo
    do i_g=1,n_irrep
       do n_occ=1,size(eigval(i_g)%m,1)
          write(fxmo_un,*) eigvec(i_g)%m(:,n_occ,1)
          print*,eigvec(i_g)%m(:,n_occ,1)
       enddo
    enddo
    call returnclose_iounit(fxmo_un,status='keep')
  end subroutine tofix_mo
#endif

  subroutine coul_imp_grads(GxQxCOai,th_mo_coul,cpks_gmat)
    type(readwriteblocked_tapehandle),intent(inout) :: th_mo_coul
    real(kind=r8_kind),intent(out):: GxQxCOai(:,:)
    real(kind=r8_kind),intent(in):: cpks_gmat(:,:)
    ! *** end of interface ***

    integer(kind=i4_kind):: i_step
    integer(kind=i4_kind):: k, mo, stat
    integer(kind=i4_kind) :: klower, kupper, k_step, k_dim

    real(kind=r8_kind),allocatable :: QxCOai(:,:),ai(:,:)

    GxQxCOai=0.0_r8_kind

#ifdef mo_coul_onfile
    call readwriteblocked_startread(trim(tmpfile('mo_coul.dat')), th_mo_coul)
#endif

#ifndef split_co_ai
    spin: do i_spn=1,n_spin
       irr: do i_g=1,n_irrep
#ifdef mo_coul_onfile
          allocate(QxCOai(iupper-ilower+1,size(cpks,1)), &
               ai(size(cpks(1,i_g,i_spn)%HBH,1),size(cpks(1,i_g,i_spn)%HBH,2)), &
               stat=stat)
#else
          allocate(QxCOai(iupper-ilower+1,size(cpks,1)), stat=stat)
#endif
          ASSERT(stat.eq.0)

          do k=1,size(QxCOai,1)
#ifdef mo_coul_onfile
             do mo=1,size(ai,1)
                call readwriteblocked_read(ai(mo,:),th_mo_coul)
             enddo
#endif
             grads: do i_grad=1,size(cpks,1)
                if(size(cpks(i_grad,i_g,i_spn)%HBH).gt.0) then

                   QxCOai(k,i_grad)=0.0_r8_kind

#ifdef mo_coul_onfile
                   QxCOai(k,i_grad)=QxCOai(k,i_grad)+real(symmetry_data_n_partners(i_g),r8_kind)* &
                        sum(cpks(i_grad,i_g,i_spn)%HBH*ai)
#else
                   QxCOai(k,i_grad)=QxCOai(k,i_grad)+real(symmetry_data_n_partners(i_g),r8_kind)* &
                        sum(cpks(i_grad,i_g,i_spn)%HBH*cpks3c(i_g,i_spn)%co_ai(k,:,:))
#endif
                endif
             enddo grads
          enddo !k


          if(size(cpks(1,i_g,i_spn)%HBH).gt.0) then
             do i_grad=1,size(cpks,1)
                call dgemm('n','n', n_ch,1,iupper-ilower+1, &
                     2.0_r8_kind*(3-n_spin),cpks_gmat(:,ilower:iupper),n_ch, &
                     QxCOai(:,i_grad),iupper-ilower+1, &
                     1.0_r8_kind,GxQxCOai(:,i_grad),n_ch)
                ! this GxQxCOai is used to calculate cpks_fitcoeff_grads

             enddo
          endif

          deallocate(QxCOai,stat=stat)
          ASSERT(stat.eq.0)
#ifdef mo_coul_onfile
          deallocate(ai,stat=stat)
          ASSERT(stat.eq.0)
#endif

       enddo irr
    enddo spin
#else
    !i.e. split_co_ai

    ! number of fit functions assigned to this worker:
    i_step = iupper - ilower + 1

    !
    ! These iterations over fit function indices define format
    ! of the on-disk file mo_coul.dat and therefore have to be
    ! done identically at the place where it is written, see
    ! pert_coeff_module.f90.
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
       !           print *, "READ1: klower=", klower, "kupper=", kupper

       k_dim  = kupper - klower + 1 ! <= k_step

       split_spin: do i_spn=1,n_spin
          do i_g=1,n_irrep
             allocate(QxCOai(k_dim,size(cpks,1)), &
                  ai(size(cpks(1,i_g,i_spn)%HBH,1),size(cpks(1,i_g,i_spn)%HBH,2)), &
                  stat=stat)
             ASSERT(stat.eq.0)

             do k=1,k_dim
                do mo=1,size(ai,1)
                   call readwriteblocked_read(ai(mo,:),th_mo_coul)
                enddo
                do i_grad=1,size(cpks,1)
                   if(size(cpks(i_grad,i_g,i_spn)%HBH).gt.0) then

                      QxCOai(k,i_grad)=0.0_r8_kind

                      QxCOai(k,i_grad)=QxCOai(k,i_grad)+real(symmetry_data_n_partners(i_g),r8_kind)* &
                           sum(cpks(i_grad,i_g,i_spn)%HBH*ai)
                   endif
                enddo
             enddo !k


             if(size(cpks(1,i_g,i_spn)%HBH).gt.0) then
                do i_grad=1,size(cpks,1)
                   ASSERT(k_step>0)
                   ASSERT(k_dim>0)
!                       print*,'dbg: from',klower,'to',kupper,'dim',kupper-klower+1
                   call dgemm('n','n', n_ch,1,k_dim, &
                        2.0_r8_kind*(3-n_spin),cpks_Gmat(:,ilower+klower-1:ilower+kupper-1),n_ch, &
                        QxCOai(:,i_grad),size(QxCOai,1), &
                        1.0_r8_kind,GxQxCOai(:,i_grad),n_ch)
                   ! this GxQxCOai is used to calculate cpks_fitcoeff_grads

                enddo
             endif

             deallocate(QxCOai,ai,stat=stat)
             ASSERT(stat.eq.0)

          enddo
       enddo split_spin
    enddo ! while
#endif


#ifdef mo_coul_onfile
    call readwriteblocked_stopread(th_mo_coul,'keep')
#endif

  end subroutine coul_imp_grads

  function solve1(A) result(x)
    !
    ! Solves linear equation system
    !
    !     A * x = b
    !
    ! for fixed rhs vector b:
    !
    !     b = ( 1, 0, 0, ..., 0)^T
    !
    use matrix_module, only: linsolve
    implicit none
    real(r8_kind), intent(in) :: A(:,:)
    real(r8_kind)             :: x(size(A,2)) ! result
    ! *** end of interface ***

    real(r8_kind)    :: rhs(size(A,1))

    ! fake rhs with one non-zero element:
    rhs(:) = 0.0_r8_kind
    rhs(1) = 1.0_r8_kind

    ! function linsolve provides the result

    x = linsolve(A,rhs)
  end function solve1

  subroutine solve_cpks()
    !
    ! Solves linear equation systems each prametrized by cpks(:,:,:)%m
    ! and using the resulting  coefficients forms a linear combination
    ! over CPKS iterations  in either %B (core version)  or %ABi (file
    ! version)
    !
    ! INPUT: ...  %m, %B(:,:,0:n_fin) or files
    ! OUTPUT: ... %B(:,:,n_cpks+1) or a file
    !
    implicit none
    ! *** end of interface **

    character(len=2) :: b_fn
    ! eclipse global vars:
    integer(i4_kind) :: i_spn, i_g, i_grad
    integer(i4_kind) :: n_fin
#ifndef b_on_file
    integer(i4_kind) :: n_fin1 ! debug check only
#endif
    real(r8_kind)    :: coeffs(0:n_cpks) ! replacement for %rhs

    ! Some of the code below seems  to assume that "n_fin" is the same
    ! for all combinations of  (i_grad,i_g,i_spn) Some code was always
    ! using  "n_fin"  from the  last  iteration  over  that triple  of
    ! indices in "column major" order.
    !
    ! We here will take the first instead. I cannot do better now as I
    ! dont understand what %n_fin_cylce is, AM
    !
    ! NV:  %n_fin_cylce  is  a  mask,  when  ge  0  it  indicate  that
    ! calculations for these data are not required any more
    !

    n_fin = cpks(1,1,1)%n_fin_cyc
    ! FIXME: n_fin is what?  NV:  maximum i_cpks for which all chanels
    ! contribute

    if(n_fin.lt.0) n_fin = n_cpks
    ! FIXME: in which  case is it less than  zero?  NV: if convergence
    ! by gamma is not reached

    ! so that local coeffs(:) is large enough:
    ASSERT(n_fin<=n_cpks)

    ! see unconditional accesses to cpks(1,...) below:
    ASSERT(size(cpks,1)>0)

#ifdef b_on_file
    ! first open all tapes, FIXME:  which n_fin is supposed to be used
    ! here?:
    ! do i_cpks=0,n_fin
    do i_cpks=0,n_cpks
       write(b_fn,'(i2)') i_cpks
       call readwriteblocked_startread(trim(tmpfile('cpks_b_'//adjustl(b_fn))), cpks_b_th(i_cpks))
    enddo

    do i_spn=1,n_spin
       do i_g=1,n_irrep
          ! FIXME:  dimensions  dont depend  on  i_grad,  but what  if
          ! size(cpks,1)==0?  cpks is allocated with first index equal
          ! to  number of  gradient projections  NV: the  case without
          ! gradient projections is out of consideration

          if(size(cpks(1,i_g,i_spn)%HBH).eq.0) cycle
          do i_grad=1,size(cpks,1)

             ![[==== solve CPKS linear equations for (i_grad,i_g,i_spn) ====
#if 1 /* debug checks */
             !
             ! Check  that the assumption  n_fin is  the same  for all
             ! (i_grad,i_g,i_spn)
             !
             n_fin = cpks(i_grad,i_g,i_spn)%n_fin_cyc
             if(n_fin.lt.0) n_fin = n_cpks
             !               ASSERT(n_fin1==n_fin)
#endif

#if 1 /* debug prints */
             print*, '=================================================='
             print*, 'CPKS matrix, irrep=',i_g,'grad=',i_grad
             print*, '--------------------------------------------------'
             do i_cpks=0,n_fin
                write(*,'(10E10.2)') cpks(i_grad,i_g,i_spn)%m(i_cpks,0:n_fin)
             enddo
#endif

             ! solve a special linear  equation system with a provided
             ! matrix and fixed rhs of the form b = (1, 0, ..., 0)
             coeffs(0:n_fin) = solve1( cpks(i_grad,i_g,i_spn)%m(0:n_fin,0:n_fin) )

#if 1 /* debug prints */
             print*, 'Solution:'
             write(*,'(10E10.2)') coeffs(0:n_fin)
#endif
             !]]============================================================

             ![[======== sum up using the resulting coefficients ===========
             !
             ! Now read and sum up data from tapes 0:n_fin.  That data
             ! is stored  in 'column  major' order in  indices i_grad,
             ! i_g  and i_spn  of cpks(i_grad,i_g,i_spn)  so  that the
             ! loops above may not be reordered!
             !
             ! First  zero accumulator,  use %ABi  not %B  as  for "in
             ! core" code:
             cpks(i_grad,i_g,i_spn)%ABi = 0.0_r8_kind
             !               do i_cpks=0,n_fin
             do i_cpks=0,n_cpks
                do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                   call readwriteblocked_read(cpks(i_grad,i_g,i_spn)%HBH(:,k),cpks_b_th(i_cpks))
                enddo

                if(i_cpks.gt.n_fin) cycle
                cpks(i_grad,i_g,i_spn)%ABi = cpks(i_grad,i_g,i_spn)%ABi &
                     + cpks(i_grad,i_g,i_spn)%HBH * coeffs(i_cpks)
             enddo
             !]]============================================================

          enddo
       enddo
    enddo

    ! last close all tapes, delete all but tape i_cpks==0:
    ! do i_cpks=0,n_fin
    do i_cpks=0,n_cpks
       if(i_cpks.ne.0) then
          call readwriteblocked_stopread(cpks_b_th(i_cpks),'delete')
       else
          call readwriteblocked_stopread(cpks_b_th(i_cpks),'keep')
       endif
    enddo

    ! write summed-up result to a file corresponding to i_cpks==n_cpks+1
    write(b_fn,'(i2)') n_cpks+1
    call readwriteblocked_startwrite(trim(tmpfile('cpks_b_'//adjustl(b_fn))), cpks_b_th(1+n_cpks))
    do i_spn=1,n_spin
       do i_g=1,n_irrep
          ! FIXME: dimensions dont depend on i_grad, but what if size(cpks,1)==0?
          ! some gradient projections should exist
          if(size(cpks(1,i_g,i_spn)%HBH).eq.0) cycle
          do i_grad=1,size(cpks,1)
             do k=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                call readwriteblocked_write(cpks(i_grad,i_g,i_spn)%ABi(:,k),cpks_b_th(n_cpks+1))
             enddo
          enddo
       enddo
    enddo
    call readwriteblocked_stopwrite(cpks_b_th(n_cpks+1))
#else

    do i_spn=1,n_spin
       do i_g=1,n_irrep
          ! FIXME:  dimensions  dont depend  on  i_grad,  but what  if
          ! size(cpks,1)==0?  NV:  cpks is allocated  with first index
          ! equal to number of gradient projections

          if(size(cpks(1,i_g,i_spn)%HBH).eq.0) cycle
          do i_grad=1,size(cpks,1)

             ![[==== solve CPKS linear equations for (i_grad,i_g,i_spn) ====
#if 1 /* debug checks */
             !
             ! Check  that the assumption  n_fin is  the same  for all
             ! (i_grad,i_g,i_spn)
             !
             n_fin1 = cpks(i_grad,i_g,i_spn)%n_fin_cyc
             if(n_fin1.lt.0) n_fin1 = n_cpks
             ASSERT(n_fin1==n_fin)
#endif

#if 1 /* debug prints */
             print*, '=================================================='
             print*, 'CPKS matrix, irrep=',i_g,'grad=',i_grad
             print*, '--------------------------------------------------'
             do i_cpks=0,n_fin
                write(*,'(10E10.2)') cpks(i_grad,i_g,i_spn)%m(i_cpks,0:n_fin)
             enddo
#endif

             ! solve a special linear  equation system with a provided
             ! matrix and fixed rhs of the form b = (1, 0, ..., 0)
             coeffs(0:n_fin) = solve1( cpks(i_grad,i_g,i_spn)%m(0:n_fin,0:n_fin) )

#if 1 /* debug prints */
             print*, 'Solution:'
             write(*,'(10E10.2)') coeffs(0:n_fin)
#endif
             !]]============================================================

             ![[======== sum up using the resulting coefficients ===========
             WARN('branch untested after rewrite!')
             ! First zero accumulator, use %B not %ABi as for "on file" code:
             cpks(i_grad,i_g,i_spn)%B(:,:,n_cpks+1) = 0.0_r8_kind
             do i_cpks=0,n_fin
                cpks(i_grad,i_g,i_spn)%B(:,:,n_cpks+1) = cpks(i_grad,i_g,i_spn)%B(:,:,n_cpks+1) &
                     + cpks(i_grad,i_g,i_spn)%B(:,:,i_cpks) * coeffs(i_cpks)
             enddo
             !]]============================================================

          enddo
       enddo
    enddo

#endif
  end subroutine solve_cpks

  subroutine b0()

#ifdef b_on_file
    call readwriteblocked_startwrite(trim(tmpfile('cpks_b_0')), cpks_b_th(0))
    spin_loop: do i_spn=1,n_spin
       irr_Q: do i_g=1,n_irrep

          n_ai=size(cpks(1,i_g,i_spn)%Qai)
          n_holes=cpks3c(i_g,i_spn)%n_holes

          if(n_ai.le.0) cycle
          gr1: do i_grad=1,size(cpks,1)
             cpks(i_grad,i_g,i_spn)%ABi=cpks(i_grad,i_g,i_spn)%Qai
             do k=1,size(cpks(i_grad,i_g,i_spn)%Qai,2)
                call readwriteblocked_write(cpks(i_grad,i_g,i_spn)%ABi(:,k),cpks_b_th(0)) !(1)
             enddo
             cpks_gamma(0,i_grad)=cpks_gamma(0,i_grad)+symmetry_data_n_partners(i_g)* &
                  sum(cpks(i_grad,i_g,i_spn)%Qai**2)
          enddo gr1
       enddo irr_Q
    enddo spin_loop
    call readwriteblocked_stopwrite(cpks_b_th(0))
#else
    do i_spn=1,n_spin
       do i_g=1,n_irrep

          n_ai=size(cpks(1,i_g,i_spn)%Qai)
          n_holes=cpks3c(i_g,i_spn)%n_holes

          if(n_ai.le.0) cycle
          do i_grad=1,size(cpks,1)
             cpks(i_grad,i_g,i_spn)%B(:,:,0)=cpks(i_grad,i_g,i_spn)%Qai
             cpks_gamma(0,i_grad)=cpks_gamma(0,i_grad)+symmetry_data_n_partners(i_g)* &
                  sum(cpks(i_grad,i_g,i_spn)%Qai**2)
          enddo
       enddo
    enddo
#endif
  end subroutine b0

  subroutine hole_corrected_b(i_cpks)
    use occupation_module, only: occ_num
    integer(kind=i4_kind),optional,intent(in):: i_cpks
    integer(kind=i4_kind):: i_occ,i_spn,i_g,i_grad,n_ai,n_holes
    integer(kind=i4_kind):: j_vir,eig_dim,occ_dim,vir_dim
    real(r8_kind):: dgen,pocc,pvir
    spin_loop: do i_spn=1,n_spin
       irr_Q: do i_g=1,n_irrep

          n_ai=size(cpks(1,i_g,i_spn)%HBH)
          n_holes=cpks3c(i_g,i_spn)%n_holes
          eig_dim=size(occ_num(i_g)%m,1)
          occ_dim=size(cpks(1,i_g,i_spn)%HBH,1)
          vir_dim=size(cpks(1,i_g,i_spn)%HBH,2)

          abi: if(n_ai.gt.0) then

             gr1: do i_grad=1,size(cpks,1)

                if(comm_i_am_master())  then
                   dgen=(3-n_spin)*symmetry_data_n_partners(i_g)

                   if(present(i_cpks)) then
                      do i_occ=1,occ_dim
                         pocc=occ_num(i_g)%m(i_occ,i_spn)/dgen
                         do j_vir=1,vir_dim
                            pvir=1.0_r8_kind-occ_num(i_g)%m(j_vir+eig_dim-vir_dim,i_spn)/dgen
                            if(i_occ.eq.eig_dim-vir_dim+j_vir) then
                               cpks(i_grad,i_g,i_spn)%HBH(i_occ,j_vir)=0.0_r8_kind
                            else
                               cpks(i_grad,i_g,i_spn)%HBH(i_occ,j_vir)=cpks(i_grad,i_g,i_spn)%B(i_occ,j_vir,i_cpks)*pocc*pvir
                            endif
                         enddo
                      enddo
                   else
                      do i_occ=1,occ_dim
                         pocc=occ_num(i_g)%m(i_occ,i_spn)/dgen
                         do j_vir=1,vir_dim
                            pvir=1.0_r8_kind-occ_num(i_g)%m(j_vir+eig_dim-vir_dim,i_spn)/dgen
                            if(i_occ.eq.eig_dim-vir_dim+j_vir) then
                               cpks(i_grad,i_g,i_spn)%HBH(i_occ,j_vir)=0.0_r8_kind
                            else
                               cpks(i_grad,i_g,i_spn)%HBH(i_occ,j_vir)=cpks(i_grad,i_g,i_spn)%ABi(i_occ,j_vir)*pocc*pvir
                            endif
                         enddo
                      enddo
                   endif

                endif

                call comm_bcast(cpks(i_grad, i_g, i_spn)%HBH)
#if 0
                if(i_grad.eq.1) then
                   print*,'HBH spin i_g i_grad',i_spn,i_g,1
                   do j_vir=1,size(cpks(i_grad,i_g,i_spn)%HBH,2)
                      print*,j_vir,cpks(i_grad,i_g,i_spn)%HBH(:,j_vir)
                   enddo
                endif
#endif
                !  call comm_bcast(cpks(i_grad, i_g, i_spn)%B(:,:,i_cpks))
             enddo gr1
          endif abi
       enddo irr_Q
    enddo spin_loop
  end subroutine hole_corrected_b

  subroutine calc_ab(th_mo_coul,status)
    implicit none
    character(len=6), optional, intent(in):: status
    type(readwriteblocked_tapehandle),intent(inout) :: th_mo_coul
    ! *** end of interface ***

    integer(kind=i4_kind):: i_spn,i_g,i_grad,mo,k,n_ai,stat
    integer(kind=i4_kind):: klower, kupper, k_step, i_step
    real(kind=r8_kind), allocatable:: ai(:,:)

#ifdef mo_coul_onfile
    call readwriteblocked_startread(trim(tmpfile('mo_coul.dat')), th_mo_coul)
#endif

#ifndef split_co_ai
    AB_spin_loop: do i_spn=1,n_spin
       AB_irr: do i_g=1, n_irrep
          n_ai=cpks3c(i_g,i_spn)%n_ai
          if(n_ai.gt.0) then
#ifdef mo_coul_onfile
             allocate(ai(size(cpks(1,i_g,i_spn)%ABi,1),size(cpks(1,i_g,i_spn)%ABi,2)),stat=stat)
             ASSERT(stat.eq.0)

             do i_grad=1,size(cpks,1)
                cpks(i_grad,i_g,i_spn)%ABi(:,:)=0.0_r8_kind
             enddo

             do k=1,iupper-ilower+1
                do mo=1,size(cpks(1,i_g,i_spn)%ABi,1)
                   call readwriteblocked_read(ai(mo,:),th_mo_coul)
                enddo

#endif

                AB_grads: do i_grad=1,size(cpks,1)
#ifdef mo_coul_onfile
                   cpks(i_grad,i_g,i_spn)%ABi(:,:)=cpks(i_grad,i_g,i_spn)%ABi(:,:)+ &
                        ai*GxQxCOai(ilower+k-1,i_grad)
#else
                   call dgemm('t','n', n_ai,1,iupper-ilower+1,   &
                        1.0_r8_kind, cpks3c(i_g,i_spn)%co_ai,iupper-ilower+1, &
                        GxQxCOai(ilower:iupper,i_grad),iupper-ilower+1, &
                        0.0_r8_kind,cpks(i_grad,i_g,i_spn)%ABi(:,:),n_ai)
#endif


#if 0 /* test print */
                   if(i_grad.eq.1.and.size(cpks(i_grad,i_g,i_spn)%ABi).ne.0) then
                      print*,i_grad,i_g,i_spn,'ABi'
                      do i=1,size(cpks(i_grad,i_g,i_spn)%ABi,2)
                         print*,cpks(i_grad,i_g,i_spn)%ABi(:,i),i
                      enddo
                   endif
#endif
                enddo AB_grads
#ifdef mo_coul_onfile
             enddo
             deallocate(ai,stat=stat)
             ASSERT(stat.eq.0)
#endif

          endif ! n_ai
       enddo AB_irr
    enddo AB_spin_loop

#else
    !i.e. split_co_ai
    do i_spn=1,n_spin
       do i_g=1, n_irrep
          do i_grad=1,size(cpks,1)
             cpks(i_grad,i_g,i_spn)%ABi(:,:)=0.0_r8_kind
          enddo
       enddo
    enddo

    ! number of fit functions assigned to this worker:
    i_step = iupper - ilower + 1

    !
    ! These iterations over fit  function indices define format of the
    ! on-disk  file   mo_coul.dat  and  therefore  have   to  be  done
    ! identically   at   the   place   where  it   is   written,   see
    ! pert_coeff_module.f90.
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
       !       print *, "READ2: klower=", klower, "kupper=", kupper

       do i_spn = 1, n_spin
          do i_g = 1, n_irrep
             n_ai = cpks3c(i_g, i_spn)%n_ai
             if ( n_ai <= 0 ) cycle
             allocate(ai(size(cpks(1, i_g, i_spn)%ABi, 1), size(cpks(1, i_g, i_spn)%ABi, 2)), stat=stat)
             ASSERT(stat.eq.0)


             do k = klower, kupper
                do mo = 1, size(cpks(1, i_g, i_spn)%ABi, 1)
                   call readwriteblocked_read(ai(mo, :), th_mo_coul)
                enddo

                do i_grad = 1, size(cpks, 1)
                   cpks(i_grad, i_g, i_spn)%ABi(:, :) = cpks(i_grad, i_g, i_spn)%ABi(:, :) &
                        + ai * GxQxCOai(ilower + k - 1, i_grad)
                enddo
             enddo

             deallocate(ai, stat=stat)
             ASSERT(stat.eq.0)
          enddo
       enddo
    enddo ! while
#endif

#ifdef mo_coul_onfile
    if(present(status)) then
       call readwriteblocked_stopread(th_mo_coul,status)
    else
       call readwriteblocked_stopread(th_mo_coul,'keep')
    endif
#endif
#if 0 /* debug print */
    do i_spn=1,n_spin
       do i_g=1, n_irrep
          do i_grad=1,size(cpks,1)
             print*,'abi',sum( cpks(i_grad,i_g,i_spn)%ABi(:,:))
          enddo
       enddo
    enddo
#endif

  end subroutine calc_ab

  subroutine hole_corrected_s1()
    use occupation_module, only: occ_num
    integer(kind=i4_kind):: occ_dim,n_holes,i_occ,i_spn,i_g,i_grad
    integer(kind=i4_kind):: j_occ
    real(r8_kind):: dgen

    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          dgen=(3-n_spin)*symmetry_data_n_partners(i_g)
          occ_dim=size(cpks(1,i_g,i_spn)%s1,1)
          n_holes=cpks3c(i_g,i_spn)%n_holes
          do i_grad=1,size(cpks,1)
#if 0
             do i_occ=1,n_holes

                cpks(i_grad,i_g,i_spn)%s1(occ_dim-n_holes+i_occ,:)= &
                     cpks(i_grad,i_g,i_spn)%s1(occ_dim-n_holes+i_occ,:)* &
                     cpks3c(i_g,i_spn)%parocc_occ(occ_dim-n_holes+i_occ)
             enddo
             do i_occ=1,occ_dim-n_holes
                cpks(i_grad,i_g,i_spn)%s1(i_occ,occ_dim)= &
                     cpks(i_grad,i_g,i_spn)%s1(i_occ,occ_dim)*cpks3c(i_g,i_spn)%parocc_occ(occ_dim)
             enddo
#else
             do i_occ=1,occ_dim
                do j_occ=1,i_occ
                   cpks(i_grad,i_g,i_spn)%s1(i_occ,j_occ)= &
                        cpks(i_grad,i_g,i_spn)%s1(i_occ,j_occ)*occ_num(i_g)%m(i_occ,i_spn)/dgen
                   cpks(i_grad,i_g,i_spn)%s1(j_occ,i_occ)=cpks(i_grad,i_g,i_spn)%s1(i_occ,j_occ)
                enddo
             enddo
#endif
          enddo
       enddo
    enddo
  end subroutine hole_corrected_s1

  subroutine hole_corrected_h1()
    use occupation_module, only: occ_num
    integer(kind=i4_kind):: occ_dim,n_holes,i_occ,i_spn,i_g,i_grad
    integer(kind=i4_kind):: j_occ
    real(r8_kind):: dgen

    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          dgen=(3-n_spin)*symmetry_data_n_partners(i_g)
          occ_dim=size(cpks(1,i_g,i_spn)%h1,1)
          if(occ_dim.gt.0) then
             n_holes=cpks3c(i_g,i_spn)%n_holes
             do i_grad=1,size(cpks,1)
                if(comm_i_am_master()) then
#if 0 /* this is code for one hole in spin-irrep */
                   do i_occ=1,n_holes
                      cpks(i_grad,i_g,i_spn)%h1(occ_dim-n_holes+i_occ,:)= &
                           cpks(i_grad,i_g,i_spn)%h1(occ_dim-n_holes+i_occ,:)* &
                           cpks3c(i_g,i_spn)%parocc_occ(occ_dim-n_holes+i_occ)
                   enddo
                   do i_occ=1,occ_dim-n_holes
                      cpks(i_grad,i_g,i_spn)%h1(i_occ,occ_dim)= &
                           cpks(i_grad,i_g,i_spn)%h1(i_occ,occ_dim)*cpks3c(i_g,i_spn)%parocc_occ(occ_dim)
                   enddo
#else
                   do i_occ=1,occ_dim
                      do j_occ=1,i_occ
                         cpks(i_grad,i_g,i_spn)%h1(i_occ,j_occ)= &
                              cpks(i_grad,i_g,i_spn)%h1(i_occ,j_occ)*occ_num(i_g)%m(i_occ,i_spn)/dgen
                         cpks(i_grad,i_g,i_spn)%h1(j_occ,i_occ)=cpks(i_grad,i_g,i_spn)%h1(i_occ,j_occ)
                      enddo
                   enddo
#endif
                endif
                call comm_bcast(cpks(i_grad, i_g, i_spn)%h1)
             enddo
          endif
       enddo
    enddo
  end subroutine hole_corrected_h1

  subroutine eigval_weight_ab()
    integer(kind=i4_kind):: spin,occ_dim,n_holes

    spin_loop: do spin=1,n_spin
       weightAB_irr:  do i_g=1,n_irrep
          n_ai=cpks3c(i_g,spin)%n_ai
          if(n_ai.gt.0) then
             n_holes=cpks3c(i_g,spin)%n_holes
             occ_dim=size(cpks(1,i_g,spin)%ABi,1)
             do n_occ=1,occ_dim
                do n_vir=1,size(cpks(1,i_g,spin)%ABi,2)
                   if(n_occ-occ_dim+n_holes.eq.n_vir) then
                      do i_grad=1,size(cpks,1)
                         cpks(i_grad,i_g,spin)%ABi(n_occ,n_vir)=0.0_r8_kind
                      enddo
                   else
                      do i_grad=1,size(cpks,1)
                         cpks(i_grad,i_g,spin)%ABi(n_occ,n_vir)=cpks(i_grad,i_g,spin)%ABi(n_occ,n_vir) &
                              /(eigval_occ(i_g)%m(n_occ,spin)-cpks3c(i_g,spin)%eigval_vir_and_holes(n_vir))
                      enddo
                   endif

                enddo
             enddo
          endif
       enddo weightAB_irr
    enddo spin_loop
  end subroutine eigval_weight_ab

  subroutine assembl_Qai()
    use eigen_data_module, ONLY : eigval
    integer(kind=i4_kind):: eig_dim,occ_dim,vir_dim,i_occ,stat
#if 1
    real(kind=r8_kind), allocatable:: h1ai(:,:)
    call readwriteblocked_startread(trim(tmpfile('h1ai.dat')), th_h1ai)
#endif
    do i_spn=1,size(cpks,3)
       do i_g=1,size(cpks,2)
          eig_dim=size(eigval(i_g)%m,1)
          occ_dim=size(cpks(1,i_g,i_spn)%Qai,1)
          vir_dim=size(cpks(1,i_g,i_spn)%Qai,2)
          allocate(h1ai(occ_dim,vir_dim),stat=stat)
          ASSERT(stat.eq.0)
          if(occ_dim*vir_dim.ne.0) then
             do i_grad=1,size(cpks,1)
                do i_occ=1,occ_dim
                   call readwriteblocked_read(h1ai(i_occ,:),th_h1ai)
                enddo


                do n_occ=1,occ_dim
                   do n_vir=1,vir_dim
#if 0
                      cpks(i_grad,i_g,i_spn)%h1ai(n_occ,n_vir)=cpks(i_grad,i_g,i_spn)%h1ai(n_occ,n_vir) &
                           -cpks(i_grad,i_g,i_spn)%s1ai(n_occ,n_vir)*eigval(i_g)%m(n_occ,i_spn)
#endif
                      if(abs(eigval(i_g)%m(n_occ,i_spn)- &
                           eigval(i_g)%m(eig_dim-vir_dim+n_vir,i_spn)).gt.0.000000001_r8_kind) then
                         cpks(i_grad,i_g,i_spn)%Qai(n_occ,n_vir)= ( cpks(i_grad,i_g,i_spn)%Qai(n_occ,n_vir) &
#if 1
                         +h1ai(n_occ,n_vir) &
#else
                         +cpks(i_grad,i_g,i_spn)%h1ai(n_occ,n_vir) &
#endif
                         !              -cpks(i_grad,i_g,i_spn)%s1ai(n_occ,n_vir)*eigval(i_g)%m(n_occ,i_spn) &
                         ) /(eigval(i_g)%m(n_occ,i_spn)-eigval(i_g)%m(eig_dim-vir_dim+n_vir,i_spn))
                      else
                         cpks(i_grad,i_g,i_spn)%Qai(n_occ,n_vir)=0.0_r8_kind
                      endif
                   enddo
                enddo

             enddo
          endif

#if 1
          deallocate(h1ai,stat=stat)
          ASSERT(stat.eq.0)
#endif
       enddo
    enddo
#if 1
    call readwriteblocked_stopread(th_h1ai,'delete')
#endif
  end subroutine assembl_Qai

  subroutine collect_h1()
    use msgtag_module, only: msgtag_packed_message
    real(kind=r8_kind),allocatable:: help(:,:)
    integer(kind=i4_kind):: spin,i
    if(comm_i_am_master()) then
       do spin=1,size(cpks,3)
          do i_g=1,symmetry_data_n_irreps()
             n_occ=size(cpks(1,i_g,spin)%h1,1)
             if(n_occ.gt.0) then
                allocate(help(n_occ,n_occ),stat=cpksalloc(129))
                MEMLOG(size(help))
                ASSERT(cpksalloc(129).eq.0)
                cpksalloc(129)=1
                do i_grad=1,size(cpks,1)
                   do i=2,comm_get_n_processors()
                      call comm_save_recv(i,msgtag_packed_message)
                      call communpack(help(1,1),n_occ**2,1,info)
                      cpks(i_grad,i_g,spin)%h1=cpks(i_grad,i_g,spin)%h1+help
                   end do
                enddo
                MEMLOG(-size(help))
                deallocate(help,stat=cpksalloc(129))
                ASSERT(cpksalloc(129).eq.0)
                cpksalloc(129)=1
             endif
          enddo
       enddo
    else
       do spin=1,size(cpks,3)
          do i_g=1,symmetry_data_n_irreps()
             n_occ=size(cpks(1,i_g,spin)%h1,1)
             if(n_occ.gt.0) then
                do i_grad=1,size(cpks,1)
                   call comm_init_send(comm_master_host,msgtag_packed_message)
                   call commpack(cpks(i_grad,i_g,spin)%h1(1,1),n_occ**2,1,info)
                   call comm_send()
                enddo
             endif
          enddo
       enddo
    endif
  end subroutine collect_h1

  subroutine reduce_ABai()
    !
    ! Looks  like   this  is   the  reduciton  of   cpks(i_grad,  i_g,
    ! i_spn)%ABi(:, :) on the rank-0
    !
    use comm, only: comm_reduce
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: i, j, k

    !
    ! Arrays %ABi may be 0-sized, indexing into them may be illegal:
    !
    do k = 1, size(cpks, 3) ! spin
       do j = 1, size(cpks, 2) ! irrep
          do i = 1, size(cpks, 1) ! gradient component
             call comm_reduce(cpks(i, j, k)%ABi)
          enddo
       enddo
    enddo
  end subroutine reduce_ABai

  subroutine reduce_Qai()
    !
    ! Looks like this is the reduciton of cpks(i_grad, i_g,
    ! i_spn)%Qai(:, :) on the rank-0
    !
    use comm, only: comm_reduce
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: i, j, k

    !
    ! Arrays %Qai may be 0-sized, indexing into them may be illegal:
    !
    do k = 1, size(cpks, 3) ! spin
       do j = 1, size(cpks, 2) ! irrep
          do i = 1, size(cpks, 1) ! gradient component
             call comm_reduce(cpks(i, j, k)%Qai)
          enddo
       enddo
    enddo
  end subroutine reduce_Qai

  subroutine say(phrase)
    use comm, only: comm_rank
    use iounitadmin_module, only: write_to_trace_unit
    implicit none
    character(len=*), intent(in) :: phrase
    ! *** end of interface ***

    if( comm_rank() == 0 ) then
       call write_to_trace_unit("cpks_g4constructs: "//phrase)
    endif
  end subroutine say

end subroutine cpks_g4constructs
