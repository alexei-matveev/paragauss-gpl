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
module calc3c_switches

! define FPP_TIMERS 2
# include "def.h"
FPP_TIMER_DECL(t_p1_sumup)
FPP_TIMER_DECL(t_w1_sumup)
FPP_TIMER_DECL(t_2c_total)
FPP_TIMER_DECL(t_cpks_grad_fitmat)
FPP_TIMER_DECL(t_contract_2c_dervs)
FPP_TIMER_DECL(t_calc_2c_dervs)
FPP_TIMER_DECL(t_calc_3center)
FPP_TIMER_DECL(t_graphi_sec_der)
FPP_TIMER_DECL(t_cross_nuc_dervsrho)
FPP_TIMER_DECL(t_cross_dervs_grarho)
FPP_TIMER_DECL(t_help_nuc_arr2)
FPP_TIMER_DECL(t_calc_imp_dervs)
FPP_TIMER_DECL(t_init_mocalc)
FPP_TIMER_DECL(t_grid_uni)
FPP_TIMER_DECL(t_dervs_xc)
FPP_TIMER_DECL(t_cpks_init)
FPP_TIMER_DECL(t_s1_imp_grarho)
FPP_TIMER_DECL(t_s1_imp_dervsrho)
FPP_TIMER_DECL(t_xc_qh_expl)
FPP_TIMER_DECL(t_xc_bimp_grarho)
FPP_TIMER_DECL(t_xc_ab)
FPP_TIMER_DECL(t_cpks_const)
FPP_TIMER_DECL(t_dervs_quad_sums)
FPP_TIMER_DECL(t_ph_sums)
FPP_TIMER_DECL(t_nuc3rd)
FPP_TIMER_DECL(t_dervs_imp)
FPP_TIMER_DECL(t_orbital_calc)
FPP_TIMER_DECL(t_density_calc)
FPP_TIMER_DECL(t_nuc_dervsrho)
FPP_TIMER_DECL(t_dervs_helps)
FPP_TIMER_DECL(t_h1q_dvdrho)
FPP_TIMER_DECL(t_dervs_weight)
! switches common to NEW and OLD integrals
!logical:: integralpar_cpksdervs=.false. ! off to check 3co 2nd dervs
!logical:: integralpar_2dervs=.false.
logical:: ewpcints=.false.
logical:: ewpcgrads=.false.
logical:: ewpcdervs=.false.
!logical:: integralpar_dervs=.false.
logical:: print_epe

#ifdef NEW_INTEGRALS
! switches only for the NEW integrals
logical:: overlap_var_fixed=.false.
logical:: gamma_var_fixed=.false.
logical:: harmonic_var_fixed=.false.

logical:: shift_xa=.false. ! now on to check derivatives
logical:: shift_xb=.false.

logical:: print_quad=.false.

logical:: old_3c_co=.false.
logical:: old_3c_fc=.false.
logical:: new_3c_co=.true.
logical:: new_3c_fc=.true. ! default=t, (f if cpks_noco)
logical:: new_3c_co_grad=.true.
logical:: radang_lco_grad=.true.
logical:: radang_sco_grad=.true.

logical:: new_ewpc=.false.
logical:: old_ewpc=.false.

logical:: old_solv_grad=.false.
logical:: new_solv_grad=.true.

logical:: old_potential=.false.
logical:: old_elfield=.false.
logical:: old_nuc_grad=.true.
#else
! switches only for the OLD integrals
logical:: overlap_var_fixed=.false.
logical:: gamma_var_fixed=.false.
logical:: harmonic_var_fixed=.false.
logical:: shift_xa=.false.
logical:: old_3c_co=.true.
logical:: old_3c_fc=.true.
logical:: new_3c_co=.false.
logical:: new_3c_fc=.false.
logical:: new_3c_co_grad=.false.
logical:: new_ewpc=.false.
logical:: old_ewpc=.true.
logical:: old_solv_grad=.true.
logical:: new_solv_grad=.false.
logical:: old_potential=.true.
logical:: old_elfield=.true.
logical:: old_nuc_grad=.true.
#endif
end module calc3c_switches
!================================================================
! End of public interface of module
!================================================================
