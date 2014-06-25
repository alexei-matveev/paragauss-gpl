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
!===============================================================
! Public interface of module
!===============================================================
#define FPP_DEBUG
#include <def.h>
Subroutine ll_calculate_hfc (na, nb, la, lb)
  !
  !  Purpose: calculation of all primitive 2 center orbital
  !           for hfc
  !
  !
  !  Author: DG
  !  Date:   1/2003
  !
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
      Use unique_atom_module
      Use gamma_module
      Use type_module
      Use datatype
      Use solid_harmonics_module, Only: solid_harmonics_calc, &
     & solid_harmonics_scalar
      Use solhrules_module
      Use integralpar_module
      Use options_module, Only: options_integral_expmax, options_finite_nucleus
      Use int_data_dipole_module, Only: prim_int_2cob_hfc, center1,center2
      Implicit None
!
  !== Interrupt end of public interface of module =================
!
      Integer (Kind=i4_kind), Intent (In) :: na ! number of unique atom a
      Integer (Kind=i4_kind), Intent (In) :: nb ! number of unique atom b
      Integer (Kind=i4_kind), Intent (In) :: la ! angular momentum of unique atom a
      Integer (Kind=i4_kind), Intent (In) :: lb ! angular momentum of unique atom b
!================================================================
! End of public interface of module
!================================================================
!
      Integer (Kind=i4_kind) :: naexps, nbexps
      Real (Kind=r8_kind), Pointer :: aexps (:), bexps (:)
      Integer (Kind=i4_kind) :: max_order
!
  ! constants

      Real (Kind=r8_kind), Parameter :: pi = &
     & 3.14159265358979324_r8_kind, very_small = 1.0e-100_r8_kind, &
     & very_big = 1.0e100_r8_kind, zero = 0.0_r8_kind, one = &
     & 1.0_r8_kind, two = 2.0_r8_kind, four = 4.0_r8_kind, six = &
     & 6.0_r8_kind

      Real (Kind=r8_kind), Dimension (0:8) :: dfac = (/ 1.0_r8_kind, &
     & 1.0_r8_kind, 3.0_r8_kind, 15.0_r8_kind, 105.0_r8_kind, &
     & 945.0_r8_kind, 10395.0_r8_kind, 135135.0_r8_kind, &
     & 2027025.0_r8_kind /)
      Integer (Kind=i4_kind) :: one_i, zero_i
      Integer (Kind=i4_kind) :: num, counter, m, ma, mb, alloc_stat
!
      Logical, Allocatable :: cutoff (:, :)
!
  ! help factors
      Real (Kind=r8_kind), Allocatable, Dimension (:, :) :: fact0_arr, &
     & fact1_arr, fact2_arr, fact10
      Real (Kind=r8_kind), Allocatable, Dimension (:) :: fact0, fact1, &
     & fact2, fact4, fact5, fact6, fact7, fact8,fact81, rcsabc, tau
!
  ! help arrays for gamma-function
      Real (Kind=r8_kind), Allocatable, Dimension (:, :) :: gamma_arg, &
     &  gamma_help, gamma_arg_fin_a,gamma_arg_fin_b
      Real (Kind=r8_kind), Allocatable, Dimension (:) :: gamma_arg2
  !arrays of some derivatives
      Real (Kind=r8_kind), Allocatable, Dimension (:, :, :) :: &
     & gamma_first_der
      Real (Kind=r8_kind), Allocatable, Dimension (:, :, :) :: &
     & gamma_second_der
  ! help arrays for solid harmincs
      Real (Kind=r8_kind), Allocatable :: yl_arr (:, :), yl_arr2 (:, &
     & :),yl_arr_fin_a (:, :),yl_arr_fin_b (:, :), clmamb (:, :), clmamb2 (:, :), clmamb_scalar (:),&
     clmamb_scalar_xyz (:,:), clmamb_xyz(:,:,:)
!
  ! help arrays for product_rule and diff_rule
      Real (Kind=r8_kind), Allocatable :: prod_arr (:, :, :, :, :), prod_arr1 (:, :, :, :, :), &
     & diff_arr (:, :, :),diff_arr_a (:, :, :), diff_arr0 (:, :), diff_arr0_xyz(:,:,:), diff_arr_fin(:,:,:)
!
  ! help arrays for product rule
      Real (Kind=r8_kind), Allocatable :: help_arr0 (:, :)

      Real (Kind=r8_kind) :: arg

  ! cartesian coordinates
      Real (Kind=r8_kind), Dimension (3) :: xa, xb, xc, xd
!
      Integer (Kind=i4_kind) :: i_ua, j_ua, i_l, k, l
      Integer (Kind=i4_kind) :: ly_max
      Integer (Kind=i4_kind) :: n_equals

      Real (Kind=r8_kind), Allocatable :: aexp_arr (:), bexp_arr (:)
      Type (unique_atom_type), Pointer :: ua_pointer
!
      Real (Kind=r8_kind), Allocatable :: diff_arr_xyz (:, :, :, :),diff_arr_axyz(:,:,:,:), diff_arr_xyz_a(:,:,:,:), &
     & yl_arr_xyz (:, :, :),yl_arr_axyz(:,:,:), yl_arr_2xyz (:, :, :),yl_arr_2axyz(:,:,:), diff_arr_2xyz (:, &
     & :, :, :), diff_arr_2axyz(:,:,:,:)!second derivatives
!
      Real (Kind=r8_kind), Allocatable :: overlap (:, :, :), overlap_xyz (:, :, :, :),&
           help_mat(:,:,:)
!
!
!
!
      Real (Kind=r8_kind), Allocatable :: nuc_second_der (:, :, :, :), &
     & aqapb (:, :), nuc_2grad_ac(:,:,:,:)
      Real (Kind=r8_kind), Allocatable :: finite_nuc_iso (:, :, :)
      Integer (Kind=i4_kind) :: k_gr, nm_la, nm_lb, lasq, lbsq, laposq, &
     & lbposq, l2pm
      Integer (Kind=i4_kind), Dimension (6, 2) :: meta2xyz, meta2xyzL
      Integer (Kind=i4_kind), Dimension (3) :: xyz_map
       Real (Kind=r8_kind):: cexp
      Intrinsic Max
!
!




!--------------------- Initialisation-----------------------
! using mapping x=1 , y=2, z=3 for two indecses use metaindex
! Integrals type1                          type2
!  x y z                                    x y z
!x 1 4 5                                  x   4 5
!y   2 6                                  y 1   6
!z     3                                  z 2 3
! meta2xyz(k_gr,{1,2})
      xyz_map (1) = 3_i4_kind
      xyz_map (2) = 4_i4_kind
      xyz_map (3) = 2_i4_kind

      meta2xyz (1, 1) = 1_i4_kind ; meta2xyz (1, 2) = 1_i4_kind
      meta2xyz (2, 1) = 2_i4_kind ; meta2xyz (2, 2) = 2_i4_kind
      meta2xyz (3, 1) = 3_i4_kind ; meta2xyz (3, 2) = 3_i4_kind
      meta2xyz (4, 1) = 1_i4_kind ; meta2xyz (4, 2) = 2_i4_kind
      meta2xyz (5, 1) = 1_i4_kind ; meta2xyz (5, 2) = 3_i4_kind
      meta2xyz (6, 1) = 2_i4_kind ; meta2xyz (6, 2) = 3_i4_kind
!---------------------------------
      meta2xyzL (1, 1) = 2_i4_kind ; meta2xyzL (1, 2) = 1_i4_kind
      meta2xyzL (2, 1) = 3_i4_kind ; meta2xyzL (2, 2) = 1_i4_kind
      meta2xyzL (3, 1) = 3_i4_kind ; meta2xyzL (3, 2) = 2_i4_kind
      meta2xyzL (4, 1) = 1_i4_kind ; meta2xyzL (4, 2) = 2_i4_kind
      meta2xyzL (5, 1) = 1_i4_kind ; meta2xyzL (5, 2) = 3_i4_kind
      meta2xyzL (6, 1) = 2_i4_kind ; meta2xyzL (6, 2) = 3_i4_kind

      nm_la = 2 * la + 1
      nm_lb = 2 * lb + 1
      lasq = la ** 2
      lbsq = lb ** 2
      laposq = (la+1) ** 2
      lbposq = (lb+1) ** 2
      one_i = 1_i4_kind
      zero_i = 0_i4_kind
      ly_max = Max (la, lb)
      max_order = 3 + la + lb
      naexps = unique_atoms(na)%l_ob(la)%n_exponents
      nbexps = unique_atoms(nb)%l_ob(lb)%n_exponents
!----------------------------------------------------------
      Allocate (fact0_arr(nbexps, naexps), fact1_arr(nbexps, naexps), &
     & fact2_arr(nbexps, naexps), cutoff(nbexps, naexps), &
     & Stat=alloc_stat)
      If (alloc_stat .Ne. 0) Call error_handler ("LL_CALCULATE_HFC: allocat&
     &ion (1) failed")
!
!
      xa = center1
      xb = center2

      xd = xa - xb
      aexps => unique_atoms(na)%l_ob(la)%exponents(:)
      bexps => unique_atoms(nb)%l_ob(lb)%exponents(:)
!
      arg = sum (xd**2)
!
      fact0_arr = (spread(aexps, 1, nbexps)+spread(bexps, 2, naexps))
      fact1_arr = (spread(aexps, 1, nbexps)*spread(bexps, 2, naexps))
!
      Where (fact0_arr >= very_small)! prevent division by zero
         fact2_arr = fact1_arr / fact0_arr
      Elsewhere
         fact2_arr = very_big
      End Where
!
      Where (fact2_arr*arg > options_integral_expmax())! cutoff: where almost no overlap
         cutoff = .False. ! is present calculation is not necessary
      Elsewhere
         cutoff = .True.
      End Where
!
      num = count (cutoff)
      If (num == 0) Then ! all integrals are equal zero
         prim_int_2cob_hfc = 0.0_r8_kind
         Deallocate (fact0_arr, fact1_arr, fact2_arr, cutoff, &
        & Stat=alloc_stat)
         If (alloc_stat .Ne. 0) Call error_handler ("LL_CALCULATE_HFC: deal&
        &location (1) failed")
         Return
      End If
!
      Allocate (fact0(num),rcsabc(num),&
     & fact1(num), fact2(num), fact4(num), &
     & fact5(num), fact6(num), fact7(num), &
     & fact8(num),fact81(num), &
     & tau(num), gamma_arg_fin_a(num, 3),gamma_arg_fin_b(num, 3),gamma_arg(num, 3),&
     & aexp_arr(num), bexp_arr(num), &
     & overlap(num, (la+1)**2, (lb+1)**2),&
     & overlap_xyz(num,laposq,lbposq,3),&
     & help_mat(num,laposq,lbposq),&
     & clmamb_scalar((Max(la,lb)+1)**2),&
     & clmamb_scalar_xyz((Max(la,lb)+1)**2,3),&
     & clmamb(num, (la+1)**2),&
     & clmamb_xyz(num,laposq,3),&
     & clmamb2(num, (la+1)**2), &
     & diff_arr_a(num,laposq,lbposq),&
     & diff_arr_xyz_a(num,laposq,lbposq,3),&
     & diff_arr0((la+1)**2,(lb+1)**2),&
     & diff_arr0_xyz((la+1)**2,(lb+1)**2,3),&
     & aqapb(num, 0:ly_max+2),fact10(num,(ly_max+1)**2),&
     &  Stat=alloc_stat)
      If (alloc_stat .Ne. 0) Call error_handler&
           ("LL_CALCULATE_HFC: allocation (2) failed")
!
!-----------array allocation---------------------

!------------------------------------------------
  ! List of *facts* at the beginning
  ! fact0 = a + b
  ! fact1 = a * b
  ! fact2 = a*b/(a+b)
  ! fact7= 1/sqrt(a**l*(2l-1)!!)
      fact0 = pack (fact0_arr, cutoff)
      fact1 = pack (fact1_arr, cutoff)
      fact2 = pack (fact2_arr, cutoff)
!
      aexp_arr = pack (spread(aexps, 1, nbexps), cutoff)
      bexp_arr = pack (spread(bexps, 2, naexps), cutoff)
!
!
      Deallocate (fact0_arr, fact1_arr, fact2_arr, Stat=alloc_stat)
      If (alloc_stat /= 0) Call error_handler ("LL_CALCULATE_HFC: deallocat&
     &ion (2) failed")
!
  ! gamma_arg = (a*vec_a + b*vec_b)/(a + b)
      gamma_arg (:, 1) = (pack(spread(aexps*xa(1), 1, &
     & nbexps)+spread(bexps*xb(1), 2, naexps), cutoff)) / fact0
!
      gamma_arg (:, 2) = (pack(spread(aexps*xa(2), 1, &
     & nbexps)+spread(bexps*xb(2), 2, naexps), cutoff)) / fact0
!
      gamma_arg (:, 3) = (pack(spread(aexps*xa(3), 1, &
     & nbexps)+spread(bexps*xb(3), 2, naexps), cutoff)) / fact0
!
!
  ! precalculation of solid harmonics
      clmamb_scalar = solid_harmonics_scalar (Max(la, lb), xd)

  clmamb_scalar_xyz=0.0_r8_kind

  do l=1,(ly_max+1)**2
     do k_gr=1,3
        do k=1,solhrules_differential(xyz_map(k_gr),l)%n_summands
           clmamb_scalar_xyz(l,k_gr)=clmamb_scalar_xyz(l,k_gr)+&
                solhrules_differential(xyz_map(k_gr),l)%coef(k)*&
                clmamb_scalar(solhrules_differential(xyz_map(k_gr),l)%lm_sh(k))
        end do
     end do
  end do!loop over l

      fact4 = 1.0_r8_kind
      counter = 1
      tau = fact2 * arg ! a*b/(a+b)*(A-B)**2
      Do l = 0, la
         Do m = 1, 2 * l + 1
            clmamb (:, counter) = clmamb_scalar (counter) * fact4
            clmamb2 (:, counter) = clmamb_scalar (counter) * fact4 * &
           & (tau-real(l, kind=r8_kind))
            do k_gr=1,3
               clmamb_xyz(:,counter,k_gr)=clmamb_scalar_xyz(counter,k_gr)*fact4

            end do
            counter = counter + 1
         End Do
         fact4 = - fact4 * fact2 * 2.0_r8_kind
      End Do
  ! calculate derivatives of solid harmonics

  ! first calculating 2-center integrals----------------
  ! fact5=fact2*(3.0_r8_kind-2.0_r8_kind*tau+2.0_r8_kind*la)
  ! a*b/(a+b)(3-2*tau+2*l)
      fact6 = 1.0_r8_kind / Sqrt (aexp_arr**la*dfac(la)) / Sqrt &
     & (bexp_arr**lb*dfac(lb)) * Exp (-tau) * (4.0_r8_kind*fact2/fact0) &
     & ** 0.75_r8_kind
      fact5 = fact2 * fact6
      fact7 = (fact2*2.0_r8_kind) ** lb
      counter = 1
      Do i_l = 0, lb
         Do mb = 1, 2 * i_l + 1
       ! calculate C1m(grad a(i))Clm(a-b)
            diff_arr0 (:, counter) = reshape &
           & (diff_rule(spread(clmamb_scalar, 1, 1), 1, (la+1)**2, counter), (/ (la+1)**2 /))
             ! now calculate derivatives of diff_arr0
            do k_gr=1,3
               diff_arr0_xyz(:,counter,k_gr)=reshape(diff_rule(&
                    spread(clmamb_scalar_xyz(:,k_gr),1,1),1,laposq,counter),(/laposq/))
            end do
            counter = counter + 1
         End Do
      End Do
    counter=1
    help_mat(:,:,1)= 0.0_r8_kind
    diff_arr_xyz_a=spread(diff_arr0_xyz(:,:,:),1,num)
    diff_arr_a=spread(diff_arr0(:,:),1,num)
    do i_l=0,lb
       help_mat(:,:,1)=spread(fact6*(2.0_r8_kind*fact2)**i_l,&
            2,laposq)
       magnetic_number_b: do mb=1,2*i_l+1
          ! overlap
          overlap(:,1:laposq,counter)=help_mat(:,:,1)*&
               prod_rule(spread(diff_arr0(:,counter),1,num),clmamb(:,:),1,laposq)
          ! now calculate derivatives of overlap
          do k_gr=1,3
             do l=1,laposq
                overlap_xyz(:,l,counter,k_gr)=overlap(:,l,counter)*&
                     (-2.0_r8_kind*fact2)*xd(k_gr)
             end do
             overlap_xyz(:,1:laposq,counter,k_gr)=overlap_xyz(:,1:laposq,counter,k_gr)+&
                  help_mat(:,:,1)*&
                  (prod_rule(diff_arr_xyz_a(:,:,counter,k_gr),clmamb(:,:),1,laposq)+&
                  & prod_rule(diff_arr_a(:,:,counter), clmamb_xyz(:,:,k_gr),1,laposq))
          end do
          counter=counter+1
       end do magnetic_number_b
    enddo
  !    Call integral_interrupt_2cob3c ()


      fact8 = 2.0_r8_kind * Sqrt (fact0/pi)
      aqapb (:, 0) = 1.0_r8_kind
      aqapb (:, 1) = aexp_arr / fact0
      Do i_l = 2, ly_max + 2
         aqapb (:, i_l) = aqapb (:, i_l-1) * aqapb (:, 1)
      End Do

!
         Allocate (gamma_help(num, max_order), Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CACLULATE_HFC : alloc&
              &ation (3.1) failed")
         Allocate ( gamma_second_der(num,1+la+lb, 6) , Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CACLULATE_HFC : alloc&
              &ation (3.2) failed")
         Allocate ( gamma_first_der(num, 1+la+lb, 3) , Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CACLULATE_HFC : alloc&
              &ation (3.3) failed")
         Allocate ( gamma_arg2(num) , Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CACLULATE_HFC : alloc&
              &ation (3.4) failed")
         Allocate (  help_arr0(num, laposq), Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CACLULATE_HFC : alloc&
              &ation (3.5) failed")
         Allocate (nuc_second_der(num, nm_la, nm_lb, 9),nuc_2grad_ac(num, nm_la, nm_lb, 6), Stat=alloc_stat)
         If (alloc_stat .Ne. 0) Call error_handler ("LL_CALCULATE_HFC: allocat&
              &ion (4) failed")
       !  Allocate (finite_nuc_iso(num, nm_la, nm_lb), Stat=alloc_stat)
          Allocate (finite_nuc_iso(num, (la+1)**2, (lb+1)**2), Stat=alloc_stat)
         If (alloc_stat .Ne. 0) Call error_handler ("LL_CALCULATE_HFC: allocat&
              &ion (4) failed")
        !7,8,9 - include isotropic part
        ! nuc_second_der = 0.0_r8_kind
        ! nuc_2grad_ac = 0.0_r8_kind
!

         Allocate (yl_arr(num, (ly_max+1)**2), yl_arr2(num, (ly_max+1)**2),&
              & yl_arr_fin_a(num, (ly_max+1)**2),yl_arr_fin_b(num, (ly_max+1)**2),&
              & prod_arr(num,(la+1)**2, (lb+1)**2, 0:la+lb, 6),&
              & prod_arr1(num,(la+1)**2, (lb+1)**2, 0:la+lb, 6),&
              & yl_arr_xyz(num, (ly_max+1)**2, 3),&
              & yl_arr_axyz(num, (ly_max+1)**2, 3),&
              & yl_arr_2axyz(num,(ly_max+1)**2, 6),&
              & yl_arr_2xyz(num, (ly_max+1)**2, 6), &
              & diff_arr(num, (la+1)**2, (lb+1)**2),&
              & diff_arr_fin(num, (la+1)**2, (lb+1)**2),&
              & diff_arr_xyz(num, laposq, lbposq, 3),&
              & diff_arr_axyz(num, laposq, lbposq, 3),&
              & diff_arr_2xyz(num, laposq, lbposq, 6),&
              & diff_arr_2axyz(num, laposq, lbposq, 6),&
         &Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CACLULATE_HFC : alloc&
              &ation (5) failed")
         ! Initialisation

  ! do a precalculation of a factor needed for the
  ! product rule
         counter = 1
         fact4 = 1.0_r8_kind
         Do i_l = 0, ly_max
            Do ma = 1, 2 * i_l + 1
               fact10 (:, counter) = fact4
               counter = counter + 1
            End Do
            fact4 = fact4 * aexp_arr / (fact0)
         End Do
!
      !   main_ua_loop: Do i = 1, n_unique_atoms

         ualoop: Do i_ua = 1,  n_unique_atoms ! loop over third center

            gamma_help = 0.0_r8_kind
            gamma_second_der = 0.0_r8_kind
            gamma_first_der = 0.0_r8_kind
            gamma_arg2 = 0.0_r8_kind
            help_arr0 = 0.0_r8_kind
            yl_arr = 0.0_r8_kind
            yl_arr_fin_a = 0.0_r8_kind
            yl_arr_fin_b = 0.0_r8_kind
            yl_arr2 = 0.0_r8_kind
            prod_arr = 0.0_r8_kind
            prod_arr1 = 0.0_r8_kind
            yl_arr_xyz = 0.0_r8_kind
            yl_arr_axyz = 0.0_r8_kind
            yl_arr_2axyz = 0.0_r8_kind
            diff_arr = 0.0_r8_kind
            diff_arr_fin = 0.0_r8_kind
            diff_arr_xyz = 0.0_r8_kind
            diff_arr_axyz = 0.0_r8_kind
            diff_arr_2xyz = 0.0_r8_kind
            diff_arr_2axyz = 0.0_r8_kind
            nuc_second_der = 0.0_r8_kind
            finite_nuc_iso =  0.0_r8_kind
            nuc_2grad_ac = 0.0_r8_kind
            fact8 = 0.0_r8_kind
            fact81 = 0.0_r8_kind
         ua_pointer => unique_atoms (i_ua)

         n_equals = ua_pointer%n_equal_atoms
  ! precalculate prod_arr and calculate nuclear attraction
    ! precalculate prod_arr and calculate nuclear attraction
         equal_atoms: Do j_ua = 1, n_equals
       !  Call integral_interrupt_2cob3c ()
!
         xc = unique_atoms(i_ua)%position(:, j_ua)
        ! xc = unique_atoms(i_ua)%position(:, 1)
        ! xc = 0.0_r8_kind


  !yl_arr([alfabeta],[lm]) = Clm(d)
            yl_arr (:, :) = solid_harmonics_calc &
                 (ly_max, gamma_arg(:,:)-spread(xc, 1, num))
  !gamma_arg2([alfabeta]) = (alfa+beta)*d*d
            gamma_arg2 (:) = ((gamma_arg(:, 1)-xc(1))**2+&
                 (gamma_arg(:, 2)-xc(2))**2+(gamma_arg(:,3)-xc(3))**2) * fact0
  !gamma_help([]) = I_l((alfa+beta)*d*d)
  !old version          gamma_help (:, 1:1+la+lb) = gamma (1+la+lb, gamma_arg2(:, j))
  ! a new version for the term 1
  ! we need I_l+2((alfa+beta)*d*d)*(2(alfa+beta))^2 * d_i * d_j
            gamma_help (:, 1:1+la+lb+2) = gamma (1+la+lb+2, gamma_arg2(:))
  ! for the double derivatives the following mapping is valued
  !
  ! Actual mapping is : C11 = z; C12 = x; C13 = y

        opt_finite:    if (options_finite_nucleus) then

               cexp =  (3.0_r8_kind/2.0_r8_kind/ua_pointer%nuclear_radius/ua_pointer%nuclear_radius)
              !  cexp = cexp + 0.000001_r8_kind
             !   xc(1) = xc(1) + 1.0_r8_kind

               !((b+c)*vec_a-c*vec_c)/b
               gamma_arg_fin_a (:, 1) = ((bexp_arr+cexp)*xa(1) - cexp * xc(1)) / bexp_arr
               !
               gamma_arg_fin_a (:, 2) = ((bexp_arr+cexp)*xa(2) - cexp * xc(2)) / bexp_arr
               !
               gamma_arg_fin_a (:, 3) = ((bexp_arr+cexp)*xa(3) - cexp * xc(3)) / bexp_arr
               !((a+c)*vec_b-c*vec_c)/a
               gamma_arg_fin_b (:, 1) = ((aexp_arr+cexp)*xb(1) - cexp * xc(1)) / aexp_arr
               !
               gamma_arg_fin_b (:, 2) = ((aexp_arr+cexp)*xb(2) - cexp * xc(2)) / aexp_arr
               !
               gamma_arg_fin_b (:, 3) = ((aexp_arr+cexp)*xb(3) - cexp * xc(3)) / aexp_arr

               ! Clm(vec_a-(b*vec_b+c*vec_c)/(b+c)) ! [ab],[lm]
                yl_arr_fin_b (:, :) = solid_harmonics_calc(ly_max, -gamma_arg_fin_b(:,:) + spread(xa, 1, num))
               ! Clm(vec_b-(a*vec_a+c*vec_c)/(a+c))
                yl_arr_fin_a (:, :) = solid_harmonics_calc(ly_max, -gamma_arg_fin_a(:,:) + spread(xb, 1, num))

                fact4 = 1.0_r8_kind
                counter = 1
                Do l = 0, la
                   Do m = 1, 2 * l + 1
                      yl_arr_fin_a (:, counter) = yl_arr_fin_a (:, counter) * fact4
                      counter = counter + 1
                   end Do
                  ! fact4 = -fact4 * 2.0_r8_kind * fact1/(fact0 + cexp)
                   fact4 = fact4 * 2.0_r8_kind * fact1/(fact0 + cexp)
                end Do
                !calculate diff array

                counter = 1
                Do i_l = 0, lb
                   Do mb = 1, 2 * i_l + 1
                      ! calculate Clm(grad a)Clm(vec_a -((a + c)*vec_b-c*vec_c)/a )
                      diff_arr_fin (:, :, counter) = diff_rule &
                           & (yl_arr_fin_b(:, :), 1, laposq, counter)
                      counter = counter + 1
                   end Do
                end Do
               ! precalculation of two factors

               fact81 =  (4.0_r8_kind*fact1/PI/PI) ** 0.75_r8_kind &
                    / Sqrt (aexp_arr**la*dfac(la))/Sqrt(bexp_arr**lb*dfac(lb))*&
                    Sqrt (cexp/(fact0+cexp)) * &
                    & (cexp/(fact0+cexp)) * Exp (-fact1/(fact0+cexp)*sum ((xa-xb)**2)-&
                    & bexp_arr*cexp/(fact0+cexp)*sum ((xb-xc)**2)-aexp_arr*cexp/(fact0+cexp)*sum ((xa-xc)**2) )


               ! now the acutal calculation starts
               counter = 1
               help_mat(:,:,1)= 0.0_r8_kind
               Do i_l = 0, lb
                  help_mat(:,:,1)=spread(fact81 *(2.0_r8_kind * fact1/(fact0 + cexp))**i_l,2,laposq)
                  Do mb = 1, 2 * i_l + 1
                     finite_nuc_iso(:,:, counter) = finite_nuc_iso(:,:,counter)&
                          + help_mat(:,:,1) * prod_rule(diff_arr_fin(:,:,counter), yl_arr_fin_a(:,:),1,laposq)
                     counter = counter + 1
                  End Do! loop over l
               End Do
            end if opt_finite


            Do i_l = 1, 1 + la + lb
               Do k_gr = 1, 3
                  gamma_first_der (:, i_l, k_gr) =&
                       & gamma_help (:, i_l+1) * 2 * fact0 * (gamma_arg(:, k_gr)-xc(k_gr))
                  gamma_second_der (:, i_l, k_gr) =&
                       & gamma_help (:, i_l+2) * 4 * fact0 * fact0 *&
                       & (gamma_arg(:,k_gr)-xc(k_gr)) * (gamma_arg(:, k_gr)-xc(k_gr))
               End Do
               Do k_gr = 4, 6
                  gamma_second_der (:, i_l, k_gr) =&
                       gamma_help (:,i_l+2) * 4 * fact0 * fact0 *&
                       & (gamma_arg(:, meta2xyz(k_gr, 1))-xc(meta2xyz(k_gr, 1))) * &
                       & (gamma_arg(:, meta2xyz(k_gr, 2))-xc(meta2xyz(k_gr, 2)))
               End Do
            End Do
!
            fact8 = 2.0_r8_kind * Sqrt (fact0/pi)
  !yl_arr2([alfabeta],[lm]) =Clm(d)/(alfa/alfa+beta)^2
            yl_arr2 (:, :) = yl_arr (:, :) / fact10 (:, 1:((ly_max+1)**2))
!
!

 ! now calculation of derivatives of yl_arr with respect to -c
            yl_arr_xyz = 0.0_r8_kind
            yl_arr_axyz = 0.0_r8_kind
        !   fact6 = spread(-aexp_arr(:)/fact0(:),1,num)
               Do l = 0, ly_max
                  Do m = 1, 2 * l + 1
                     l2pm = l * l + m
                     Do k_gr = 1, 3
                        Do k = 1, solhrules_differential(xyz_map(k_gr), l2pm)%n_summands
                           yl_arr_xyz (:, l2pm, k_gr) = yl_arr_xyz (:, l2pm, k_gr) - &    ! + & Now with respect to c
                                & solhrules_differential(xyz_map(k_gr), l2pm)%coef(k) * &
                                & yl_arr (:, solhrules_differential(xyz_map(k_gr), l2pm)%lm_sh(k))
                           yl_arr_axyz(:,l2pm, k_gr) =-aexp_arr/fact0 * yl_arr_xyz (:, l2pm, k_gr)

                        End Do
                     End Do
                  End Do
               End Do!loop over l

 !now calculate the second derivative of yl_arr with reepect to -c
            yl_arr_2xyz = 0.0_r8_kind
            yl_arr_2axyz = 0.0_r8_kind
            Do l = 0, ly_max
               Do m = 1, 2 * l + 1
                  l2pm = l * l + m

                  Do k_gr = 1, 6
                     Do k = 1, solhrules_differential(xyz_map(meta2xyz(k_gr, 2)), l2pm)%n_summands
                        yl_arr_2xyz (:, l2pm, k_gr) =yl_arr_2xyz (:,l2pm, k_gr)  - &! Now with respect to c
                       & solhrules_differential(xyz_map(meta2xyz(k_gr,2)), l2pm)%coef(k) *&
                       & yl_arr_xyz (:,solhrules_differential(xyz_map(meta2xyz(k_gr,2)), l2pm)%lm_sh(k), meta2xyz(k_gr, 1))
                        ! Debug Wrong!!!!!
                       !  yl_arr_2axyz (:, l2pm, k_gr) = -aexp_arr/fact0 *  yl_arr_2xyz (:, l2pm, k_gr)
                        ! yl_arr_2axyz (:, l2pm, k_gr) =  yl_arr_2xyz (:, l2pm, k_gr)
                     End Do
                  End Do
               End Do
            End Do !loop over l

            ! FIXME: this cannot  be correct.  See how l2pm  = l^2 + m
            ! is  used here  as a  fixed  constant, not  as a  running
            ! index.   Also cf.   the warning  issued by  the compiler
            ! that  l2pm may  be  used uninitialed  here. The  comment
            ! cannot be right:
            ABORT('needs some work')
            ! tipa pravilno ("appears correct", rus.)
             Do k_gr = 4, 6
                 yl_arr_2axyz (:, l2pm, k_gr) = -aexp_arr/fact0 *  yl_arr_2xyz (:, l2pm, k_gr)
                 yl_arr_2axyz (:, l2pm, k_gr-3) = -aexp_arr/fact0 *  yl_arr_2xyz (:, l2pm, k_gr)
             End Do
!

            !calculation of diff_arr and diff_arr_xyz
            counter = 1
            Do i_l = 0, lb
               help_arr0 = spread (aqapb(:, i_l), 2, laposq)
               Do mb = 1, 2 * i_l + 1
                  diff_arr (:, :, counter) = help_arr0 * diff_rule &
                 & (yl_arr2(:, :), 1, laposq, counter)
                   ! now let us differentiate diff_arr and second der of diff_arr
                  Do k_gr = 1, 3
                     diff_arr_xyz (:, :, counter,k_gr ) = help_arr0 * &
                    & diff_rule (yl_arr_xyz(:, :, k_gr)/fact10, 1, &
                    & laposq, counter)

                     diff_arr_axyz (:, :, counter,k_gr ) = help_arr0 * &
                    & diff_rule (yl_arr_axyz(:, :, k_gr)/fact10, 1, &
                    & laposq, counter)


                      diff_arr_2xyz (:, :, counter, k_gr) = help_arr0 * & !Coefficient
                    & diff_rule(yl_arr_2xyz(:, :, k_gr)/fact10, &
                    & 1, laposq, counter)

                     diff_arr_2axyz (:, :, counter, k_gr) = help_arr0 * & !Coefficient
                    & diff_rule(yl_arr_2axyz(:, :, k_gr)/fact10, &
                    & 1, laposq, counter)
                  End Do!
                  Do k_gr = 4, 6

                     diff_arr_2xyz (:, :, counter, k_gr) = help_arr0 * & !Coefficient
                          & diff_rule(yl_arr_2xyz(:, :, k_gr)/fact10, 1, laposq, counter)
                     diff_arr_2axyz (:, :, counter, k_gr) = help_arr0 * & !Coefficient
                          & diff_rule(yl_arr_2axyz(:, :, k_gr)/fact10, 1, laposq, counter)
                  End Do
                  counter = counter + 1
!
               End Do
            End Do
             !calculation of second derivative of diff_arr_xyz
!

!
!------------------ The first term----------------------------
!d2I((a+b)d^2)/dcidcj * X
            prod_arr = 0.0_r8_kind
            Do i_l = 0, lb
               Do mb = 1, 2 * i_l + 1
                  prod_arr (:, :, i_l**2+mb, 0:la+i_l, 1) = &
                 & prod_rule_nested2 (yl_arr(:, :), diff_arr(:, :, :), &
                 & overlap(:, :, :), la, i_l, mb, &
                 & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr))
               End Do
            End Do

               Do mb = 1, 2 * lb + 1
                  Do i_l = 0, lb + la
                     Do k_gr = 1,3
                        nuc_second_der (:, :, mb, k_gr) = nuc_second_der &
                             & (:, :, mb, k_gr) +&
                             & prod_arr (:, la**2+1:(la+1)**2, lb**2+mb, i_l, 1) *&
                             & spread (fact8*gamma_second_der(:, i_l+1, k_gr), 2, 2*la+1)
                        !isotropic part
                        nuc_second_der (:, :, mb, k_gr+6) =&
                             & nuc_second_der (:, :, mb, k_gr+6) +&
                             & prod_arr (:, la**2+1:(la+1)**2, lb**2+mb, i_l, 1) *&
                             & spread  ((-2*fact0*fact8*gamma_help(:, i_l+2)), 2,  2*la+1)
                        !------------------------------!
                     End Do
                     Do k_gr = 4,6
                        nuc_second_der (:, :, mb, k_gr) =&
                             nuc_second_der (:, :, mb, k_gr) +&
                             & prod_arr (:, la**2+1:(la+1)**2, lb**2+mb, i_l, 1) *&
                             & spread  (fact8*gamma_second_der(:, i_l+1, k_gr), 2, 2*la+1)
                        !
                     End Do
                  End Do
               End Do
! spread (fact8*gamma_help(:, i_l+1), 2, 2*la+1)
!Extra term

!---------------- The second and the third term --------------------
!dI((a+b)d^2)/dci * dX/dcj +  dI((a+b)d^2)/dcj * dX/dci
!dI((a+b)d^2)/dci * dX/daj +  dI((a+b)d^2)/daj * dX/dci
            prod_arr = 0.0_r8_kind
            prod_arr1 = 0.0_r8_kind
            Do k_gr = 1, 3
               Do i_l = 0, lb
                  Do mb = 1, 2 * i_l + 1
!dX/dci
                     prod_arr (:, :, i_l**2+mb, 0:la+i_l, k_gr) = &
                    & prod_rule_nested2&
                    & (yl_arr_xyz(:, :, k_gr),diff_arr(:, :, :), overlap(:, :, :), la, i_l, mb, &
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2&
                    & (yl_arr(:, :), diff_arr_xyz(:, :, :, k_gr), overlap(:, :, :), la, i_l, mb,&
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr))
!dX/daj
                     prod_arr1 (:, :, i_l**2+mb, 0:la+i_l, k_gr) = &
                          & prod_rule_nested2&
                          & (yl_arr(:, :),diff_arr(:, :, :), overlap_xyz(:, :, :,k_gr), la, i_l, mb, &
                          & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr)) &
                          & + prod_rule_nested2&
                          & (yl_arr_axyz(:, :,k_gr), diff_arr(:, :, :), overlap(:, :, :), la, i_l, mb,&
                          & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr))  &
                          & + prod_rule_nested2&
                          & (yl_arr(:, :), diff_arr_axyz(:, :, :,k_gr), overlap(:, :, :), la, i_l, mb,&
                          & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr))

                  End Do
               End Do
            End Do
!
!

               Do mb = 1, 2 * lb + 1
                  Do i_l = 0, lb + la
                     Do k_gr = 1, 6
                        nuc_second_der (:, :, mb, k_gr) = nuc_second_der (:, :, mb, k_gr) +&
                             & prod_arr (:, la**2+1:(la+1)**2, lb**2+mb, i_l, meta2xyz(k_gr, 1)) *&
                             & spread (fact8*gamma_first_der(:, i_l+1,meta2xyz(k_gr, 2)), 2, 2*la+1) +&
                             & prod_arr (:,la**2+1:(la+1)**2, lb**2+mb, i_l, meta2xyz(k_gr, 2)) *&
                             & spread (fact8*gamma_first_der(:, i_l+1, meta2xyz(k_gr, 1)), 2, 2*la+1)

                        nuc_2grad_ac (:, :, mb, k_gr) = nuc_2grad_ac (:, :, mb, k_gr) +&
                             & prod_arr1 (:, la**2+1:(la+1)**2, lb**2+mb, i_l, meta2xyzL(k_gr, 1)) *&
                             & spread (fact8*gamma_first_der(:, i_l+1,meta2xyzL(k_gr, 2)), 2, 2*la+1) +&
                             & prod_arr (:,la**2+1:(la+1)**2, lb**2+mb, i_l, meta2xyzL(k_gr, 2)) *&
                             & spread (fact8*(-aexp_arr/fact0)*gamma_first_der(:, i_l+1, meta2xyzL(k_gr, 1)), 2, 2*la+1)
                     End Do
               End Do
            End Do
!

!------------------the last term-----------------
!I((a+b)d^2) * d2X/dcidcj

            prod_arr = 0.0_r8_kind
            prod_arr1 = 0.0_r8_kind
            Do k_gr = 1, 6
               Do i_l = 0, lb
                  Do mb = 1, 2 * i_l + 1
                     prod_arr (:, :, i_l**2+mb, 0:la+i_l, k_gr) = &
                    & prod_rule_nested2&
                    & (yl_arr_2xyz(:, :, k_gr),diff_arr(:, :, :), overlap(:, :, :), la, i_l, mb, &
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2&
                    & (yl_arr_xyz(:, :, meta2xyz(k_gr, 1)), diff_arr_xyz(:, :, :, meta2xyz(k_gr, 2)),overlap(:, :, :), la, i_l,mb,&
                    & (-2.0_r8_kind*aexp_arr),(-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2 &
                    & ( yl_arr_xyz(:, :, meta2xyz(k_gr, 2)) , diff_arr_xyz(:,:,:,meta2xyz(k_gr, 1)),overlap(:, :, :), la, i_l, mb,&
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2&
                    & (yl_arr(:, :),diff_arr_2xyz(:, :, :, k_gr),overlap(:, :, :),la, i_l, mb,&
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr))
!
!d2X/daidcj

                      prod_arr1 (:, :, i_l**2+mb, 0:la+i_l, k_gr) = &
                    & prod_rule_nested2&
                    & (yl_arr_xyz(:, :,meta2xyzL(k_gr,2)),diff_arr(:, :, :),&
                    & overlap_xyz(:, :, :,meta2xyzL(k_gr,1)), la, i_l, mb, &
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2&
                    & (yl_arr(:, :), diff_arr_xyz(:, :, :, meta2xyzL(k_gr,2)),&
                    & overlap_xyz(:, :, :,meta2xyzL(k_gr,1)), la, i_l,mb,&
                    & (-2.0_r8_kind*aexp_arr),(-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2 &
                    & ( yl_arr_2axyz(:, :, k_gr) , diff_arr(:,:,:),overlap(:, :, :), la, i_l, mb,&
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2&
                    & (yl_arr_axyz(:, :,meta2xyzL(k_gr,1)),diff_arr_xyz(:, :, :, meta2xyzL(k_gr,2)),&
                    & overlap(:, :, :),la, i_l, mb,&
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr))&
                         & + prod_rule_nested2 &
                    & ( yl_arr_xyz(:, :, meta2xyzL(k_gr,2)) , diff_arr_axyz(:,:,:,meta2xyzL(k_gr,1)),&
                    & overlap(:, :, :), la, i_l, mb,&
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr)) &
                    & + prod_rule_nested2&
                    & (yl_arr(:, :),diff_arr_2axyz(:, :, :, k_gr),overlap(:, :, :),la, i_l, mb,&
                    & (-2.0_r8_kind*aexp_arr), (-2.0_r8_kind*bexp_arr))
                  End Do
               End Do
            End Do



            Do k_gr = 1, 6
               Do mb = 1, 2 * lb + 1
                  Do i_l = 0, lb + la
                     nuc_second_der (:, :, mb, k_gr) =&
                          nuc_second_der (:, :, mb, k_gr) +&
                          prod_arr (:, la**2+1:(la+1)**2, lb**2+mb, i_l, k_gr) *&
                          spread (fact8*gamma_help(:, i_l+1), 2, 2*la+1)
                     nuc_2grad_ac (:, :, mb, k_gr) =&
                          nuc_2grad_ac (:, :, mb, k_gr) +&
                          prod_arr1 (:, la**2+1:(la+1)**2, lb**2+mb, i_l, k_gr) *&
                          spread (fact8*gamma_help(:, i_l+1), 2, 2*la+1)
                  End Do
               End Do
            End Do
!


end Do equal_atoms


    !  End Do ualoop
!--------------------- Initialisation-----------------------
! using mapping x=1 , y=2, z=3 for two indecses use metaindex
! Integrals type1                          type2
!  x y z                                    x y z
!x 1 4 5                                  x   4 5
!y   2 6                                  y 1   6
!z     3                                  z 2 3
! meta2xyz(k_gr,{1,2})
         ! Mapping to contracted primitive integrals
         ! d2/dx2 means d2(1/R)/dx2 = (3xx -rr)/r^5
         ! 1 = d2/dx2   4 = d2/dxdy
         ! 2 = d2/dy2   5 = d2/dxdz
         ! 3 = d2/dz2   6 = d2/dydz
         ! 7 =  Lx/r^3
         ! 8 =  Ly/r^3
         ! 9 =  Lz/r^3
         do k_gr = 1,6
            do mb=1,2*lb+1
               do ma=1,2*la+1
                  prim_int_2cob_hfc(:,:,mb,ma,k_gr,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,k_gr,i_ua)+&
                       unpack(nuc_second_der(:,ma,mb,k_gr)/n_equals,cutoff,zero)!

               enddo
            end do
         end do

         if (options_finite_nucleus) then
         ! coef_finite_nucleus = (3/2/ua_pointer%nuclear_radius/ua_pointer%nuclear_radius/PI)**(1.5_r8_kind)
          do mb=1,2*lb+1
               do ma=1,2*la+1
                  prim_int_2cob_hfc(:,:,mb,ma,7,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,7,i_ua)+& !Lx/r^3
                       unpack((nuc_2grad_ac (:, ma, mb, 6)-nuc_2grad_ac (:, ma, mb, 3))/n_equals , cutoff,zero)! / n_equals
                  prim_int_2cob_hfc(:,:,mb,ma,8,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,8,i_ua)+& !Ly/r^3
                       unpack((nuc_2grad_ac (:, ma, mb, 2)-nuc_2grad_ac (:, ma, mb, 5))/n_equals , cutoff,zero)! / n_equals
                  prim_int_2cob_hfc(:,:,mb,ma,9,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,9,i_ua)+& !Lz/r^3
                       unpack((nuc_2grad_ac (:, ma, mb, 4)-nuc_2grad_ac (:, ma, mb, 1))/n_equals ,cutoff,zero)! / n_equals

                 prim_int_2cob_hfc(:,:,mb,ma,10,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,10,i_ua)+&
                     unpack(finite_nuc_iso(:,(la)**2+ma,(lb)**2+mb)/n_equals,cutoff,zero)

               enddo
            end do
         else
             do mb=1,2*lb+1
               do ma=1,2*la+1
                  prim_int_2cob_hfc(:,:,mb,ma,7,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,7,i_ua)+& !Lx/r^3
                       unpack((nuc_2grad_ac (:, ma, mb, 6)-nuc_2grad_ac (:, ma, mb, 3))/n_equals , cutoff,zero)! / n_equals
                  prim_int_2cob_hfc(:,:,mb,ma,8,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,8,i_ua)+& !Ly/r^3
                       unpack((nuc_2grad_ac (:, ma, mb, 2)-nuc_2grad_ac (:, ma, mb, 5))/n_equals , cutoff,zero)! / n_equals
                  prim_int_2cob_hfc(:,:,mb,ma,9,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,9,i_ua)+& !Lz/r^3
                       unpack((nuc_2grad_ac (:, ma, mb, 4)-nuc_2grad_ac (:, ma, mb, 1))/n_equals ,cutoff,zero)! / n_equals
                  prim_int_2cob_hfc(:,:,mb,ma,10,i_ua)=prim_int_2cob_hfc(:,:,mb,ma,10,i_ua)+& ! d2/dx2 + d2/dy2+d2/dz2 (1/r)
                       unpack((nuc_second_der(:,ma,mb,1)+nuc_second_der(:,ma,mb,2)+&          ! avoiding delta function operator
                       nuc_second_der(:,ma,mb,3)+nuc_second_der(:,ma,mb,7)+nuc_second_der(:,ma,mb,8)+&
                       nuc_second_der(:,ma,mb,9))/n_equals,cutoff,zero)!

               enddo
            end do
         end if
         end Do ualoop
!
 Deallocate (yl_arr,yl_arr_fin_a,yl_arr_fin_b, yl_arr2,&
              & yl_arr_xyz,yl_arr_axyz,&
              & yl_arr_2xyz,yl_arr_2axyz, &
              & gamma_help,&
              & gamma_second_der,&
              & gamma_first_der,&
              & gamma_arg2, &
              & prod_arr,&
              & prod_arr1,&
              & help_arr0, Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CACLULATE_HFC : deall&
        &ocation (4) failed")
         Deallocate (nuc_second_der,&
              & finite_nuc_iso,&
              & nuc_2grad_ac,&
              & diff_arr, &
              & diff_arr_fin,&
              & diff_arr_xyz,&
              & diff_arr_axyz,&
              & diff_arr_2xyz,&
              & diff_arr_2axyz, Stat=alloc_stat)
         If (alloc_stat /= 0) Call error_handler ("LL_CALCULATE_HFC: deallocat&
              &ion (5) failed")

      Deallocate (clmamb, clmamb2, tau, clmamb_scalar, clmamb_scalar_xyz,clmamb_xyz, fact10, fact6, &
     & fact7, fact2, fact1, Stat=alloc_stat)
      If (alloc_stat /= 0) Call error_handler ("LL_CALCULATE_HFC: deallocat&
     &ion (6) failed")
      Deallocate (fact8,fact81, fact5, fact0,rcsabc, fact4, overlap,overlap_xyz,help_mat,&
           gamma_arg,gamma_arg_fin_a, gamma_arg_fin_b, &
     & Stat=alloc_stat)
      If (alloc_stat /= 0) Call error_handler ("LL_CALCULATE_HFC: deallocat&
     &ion (7) failed")
      Deallocate (aexp_arr, bexp_arr,  diff_arr0, diff_arr0_xyz, &
     & diff_arr_a, diff_arr_xyz_a, Stat=alloc_stat)
      If (alloc_stat /= 0) Call error_handler ("LL_CALCULATE_HFC: deallocat&
     &ion (8) failed")
      Deallocate (cutoff, Stat=alloc_stat)
      If (alloc_stat /= 0) Call error_handler ("LL_CALCULATE_HFC: deallocat&
     &ion (9) failed")
      Deallocate (aqapb, Stat=alloc_stat)
      If (alloc_stat /= 0) Call error_handler ("LL_CALCULATE_HFC: deallocat&
     &ion (10) failed")
End Subroutine ll_calculate_hfc
