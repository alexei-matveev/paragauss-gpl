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
module msgtag_module
!---------------------------------------------------------------
!msgtag_module.f90
!  Purpose: This  module contains all the necessary message tags
!           for parallel work
!
!  Author: AS
!  Date: 7/98
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
!------------ Modules used --------------------------------------
#include "def.h"
use type_module
!---------------------------------------------------------------
implicit none
save
!------------ Declaration of constants -------------------------
!  These are message tags, for instance:
!     msgtag_pert   "Start perturbation theory, son!"
!     msgtag_end    "Sabbath, my son!"
#define INTEGER_T       integer(kind=i4_kind), parameter, public
INTEGER_T :: msgtag_pert=                   2
INTEGER_T :: msgtag_ham_calc_main =         3
INTEGER_T :: msgtag_dens=                   4
INTEGER_T :: msgtag_finito=                 5
INTEGER_T :: msgtag_eigensolve=             6
!NTEGER_T :: msgtag_charge_fit=             7
INTEGER_T :: msgtag_borders=                8
!NTEGER_T :: msgtag_finalize_geometry=      9
! do not use INTEGER_T :: msgtag_filenames=             10
INTEGER_T :: msgtag_commdata=              11
!NTEGER_T :: msgtag_init=                  12
INTEGER_T :: msgtag_debugorbitals=         13
!NTEGER_T :: msgtag_pcalc_pert_coeff=      14
INTEGER_T :: msgtag_SendBackEig=           15
INTEGER_T :: msgtag_occ_levels=            16
INTEGER_T :: msgtag_eigvec=                17
!NTEGER_T :: msgtag_gr_send=               18
INTEGER_T :: msgtag_xcham_send=            19
INTEGER_T :: msgtag_build_xc=              20
!NTEGER_T :: msgtag_xc_setup=              21
INTEGER_T :: msgtag_SendBackHam=           22
INTEGER_T :: msgtag_xcfitham_send=         23
INTEGER_T :: msgtag_build_xcfit=           24
!NTEGER_T :: msgtag_xcfit_setup=           25
!NTEGER_T :: msgtag_send_slavetiming=      26
INTEGER_T :: msgtag_slavetiming=           27
!NTEGER_T :: msgtag_grid_close=            28
! do not use 29 and 30 for msgtags, - because of compiler bugs:
! INTEGER_T :: msgtag_int_setup=             29
! INTEGER_T :: msgtag_int_setup=             30
INTEGER_T :: msgtag_int_shutdown=          31
!INTEGER_T :: msgtag_int_2cff_setup=        32
!INTEGER_T :: msgtag_int_2cff_shutdown=     33
!INTEGER_T :: msgtag_int_2cff_do=           34
!INTEGER_T :: msgtag_int_2cff_done=         35
! WAS     :: msgtag_int_2cff_norm=         36
INTEGER_T :: msgtag_int_2cff_result=       37
!INTEGER_T :: msgtag_int_2cff_rep_back=     38
!INTEGER_T :: msgtag_int_2cob3c_setup=      39
!INTEGER_T :: msgtag_int_2cob3c_shutdown=   40
!INTEGER_T :: msgtag_int_2cob3c_do=         41
!INTEGER_T :: msgtag_int_2cob3c_done=       42
!INTEGER_T :: msgtag_int_2cob_norm=         43
INTEGER_T :: msgtag_int_2cob3c_result=     44
!INTEGER_T :: msgtag_int_2cob3c_rep_back=   45
INTEGER_T :: msgtag_int_rel_setup=         46
INTEGER_T :: msgtag_int_rel_shutdown=      47
INTEGER_T :: msgtag_int_2cob3c_filesize=   48
!NTEGER_T :: msgtag_post_scf=              50
INTEGER_T :: msgtag_ph_send=               51
!NTEGER_T :: msgtag_gridph_close=          53
!NTEGER_T :: msgtag_main_gradient=         54
!NTEGER_T :: msgtag_gradient_send_3c=      55
INTEGER_T :: msgtag_grad_3c=               56
!NTEGER_T :: msgtag_gradient_send_fit_ch=  57
INTEGER_T :: msgtag_grad_ch=               58
INTEGER_T :: msgtag_occ_levels_eigval=     59
INTEGER_T :: msgtag_fit_coeff_send=        60
!NTEGER_T :: msgtag_fit_coeff_shutdown=    61
!NTEGER_T :: msgtag_xc_close=              62
!NTEGER_T :: msgtag_xcfit_close=           63
!NTEGER_T :: msgtag_openspecialunits=      64
!NTEGER_T :: msgtag_gradinfo_dealloc=      65
!NTEGER_T :: msgtag_grad_datashutdown=     66
!NTEGER_T :: msgtag_eigvec_occ_dealloc=    67
!NTEGER_T :: msgtag_eigvec_occ_dealloc_e=  68
!NTEGER_T :: msgtag_ua_close=              69
!NTEGER_T :: msgtag_op_dealloc=            70
!NTEGER_T :: msgtag_intstore_dealloc=      71
INTEGER_T :: msgtag_rel_Ttrans =           72
INTEGER_T :: msgtag_rel_Ttrans_back =      73
INTEGER_T :: msgtag_rel_momspace =         74
INTEGER_T :: msgtag_rel_momspace_back =    75
INTEGER_T :: msgtag_rel_send_overlap=      76
INTEGER_T :: msgtag_rel_trafomat=          77
INTEGER_T :: msgtag_SendBackTrafomat=      78
!NTEGER_T :: msgtag_rel_gradient_trafo=    79
INTEGER_T :: msgtag_rel_gradient_send=     80
INTEGER_T :: msgtag_xcmdaham_send=         81
INTEGER_T :: msgtag_build_xcmda=           82
!NTEGER_T :: msgtag_xcmda_setup=           83
!NTEGER_T :: msgtag_pre_dens_slave=        84
INTEGER_T :: msgtag_vir_levels=            85
INTEGER_T :: msgtag_eigvec_vir_dealloc=    86
INTEGER_T :: msgtag_rot_levels=            87
!NTEGER_T :: msgtag_xcmda_close=           88
!INTEGER_T :: msgtag_int_dipole_setup=      89
!INTEGER_T :: msgtag_int_dipole_shutdown=   90
!INTEGER_T :: msgtag_int_dipole_do=         91
!INTEGER_T :: msgtag_int_dipole_done=       92
INTEGER_T :: msgtag_int_dipole_result=     93
!INTEGER_T :: msgtag_int_dipole_rep_back=   94
INTEGER_T :: msgtag_response_setup=        95
INTEGER_T :: msgtag_response_level_send=   96
INTEGER_T :: msgtag_response_chmat=        97
INTEGER_T :: msgtag_response_2index=       98
INTEGER_T :: msgtag_resp_2in_sl_start=     99
INTEGER_T :: msgtag_response_2in_send=    100
INTEGER_T :: msgtag_response_2in_send2=   101
INTEGER_T :: msgtag_response_3index=      102
INTEGER_T :: msgtag_resp_sl_start_3ind=   103
INTEGER_T :: msgtag_resp_3ind_fitdim=     104
INTEGER_T :: msgtag_response_3ind_send=   105
INTEGER_T :: msgtag_response_3Clb_start=  106
INTEGER_T :: msgtag_resp_sl_start_3Clb=   107
INTEGER_T :: msgtag_response_3Clb_send=   108
INTEGER_T :: msgtag_send_all_levels=      111
INTEGER_T :: msgtag_resp_3Clb_tmptape=    112
INTEGER_T :: msgtag_resp_3Clb_tmptape_2=  113
!NTEGER_T :: msgtag_orb_plot_send_start=  114
!NTEGER_T :: msgtag_orbital_plot_send=    115
!NTEGER_T :: msgtag_properties_main =     116
!NTEGER_T :: msgtag_prescf_finalize =     117
!NTEGER_T :: msgtag_eigen_data_send=      118
!NTEGER_T :: msgtag_free_eigen=           119
!NTEGER_T :: msgtag_alloc_eigen=          120
INTEGER_T :: msgtag_density_data_free =   121
INTEGER_T :: msgtag_grid_send_n_points =  122
INTEGER_T :: msgtag_grid_send_points=     123
INTEGER_T :: msgtag_int_setup=            124
INTEGER_T :: msgtag_space_point=          125
!NTEGER_T :: msgtag_bounds_poten=         126
INTEGER_T :: msgtag_start_poten=          127
INTEGER_T :: msgtag_finish_poten=         128
!NTEGER_T :: msgtag_prescf_init =         129
INTEGER_T :: msgtag_geom_grad=            130
!NTEGER_T :: msgtag_Q_elec=               131
INTEGER_T :: msgtag_solv_ham=             132
!NTEGER_T :: msgtag_back_ham_solv=        133
INTEGER_T :: msgtag_surf_point=           134
INTEGER_T :: msgtag_surf_point_sa=        135
INTEGER_T :: msgtag_sp_grinfo_dealloc=    136
INTEGER_T :: msgtag_start_field=          137
INTEGER_T :: msgtag_finish_field=         138
INTEGER_T :: msgtag_start_rho=            139
INTEGER_T :: msgtag_finish_rho=           140
!NTEGER_T :: msgtag_bounds_field=         141
INTEGER_T :: msgtag_free_bnds_ptn=        142
INTEGER_T :: msgtag_free_bnds_fld=        143
!NTEGER_T :: msgtag_del_poten=            144
!NTEGER_T :: msgtag_del_field=            145
INTEGER_T :: msgtag_epe_data_sent=        146
INTEGER_T :: msgtag_epe_do_gradients=     147
INTEGER_T :: msgtag_epe_grad_done=        148
INTEGER_T :: msgtag_epe_grads_init=       149
INTEGER_T :: msgtag_epe_grads_cons=       150
INTEGER_T :: msgtag_epe_forces_done=      151
INTEGER_T :: msgtag_epe_init_slave=       152
INTEGER_T :: msgtag_epe_finish_slave=     153
INTEGER_T :: msgtag_epe_defects=          154
INTEGER_T :: msgtag_epe_send_only=        155
INTEGER_T :: msgtag_vac_intermediate=     156
INTEGER_T :: msgtag_imp_intermediate=     157
INTEGER_T :: msgtag_epe_send_def=         158
INTEGER_T :: msgtag_epe_consdef=          159
INTEGER_T :: msgtag_epe_def_fin=          160
!
INTEGER_T :: msgtag_stop_daemon    =      161
! <<< slave deamon should exit
INTEGER_T :: msgtag_packed_message =      163
! <<< any simple message, do unpacking
!NTEGER_T :: msgtag_filenames=            164
INTEGER_T :: msgtag_open_densmat =        165
!NTEGER_T :: msgtag_orb_plot_fit_sent =   166
! added by HH in pg_V2.1_hh (was 124-139, changed to ):
INTEGER_T :: msgtag_response_rhs_h=       167
INTEGER_T :: msgtag_resp_sl_start_rhs=    168
INTEGER_T :: msgtag_resp_rhs_fitdim=      169
INTEGER_T :: msgtag_response_rhs_send=    170
INTEGER_T :: msgtag_resp_rho1dat=         171
INTEGER_T :: msgtag_resp_rho1dat_1b=      172
INTEGER_T :: msgtag_resp_all3index=       173
INTEGER_T :: msgtag_resp_all3Clb_start=   174
INTEGER_T :: msgtag_resp_addFcoul=        175
INTEGER_T :: msgtag_response_FG=          176
INTEGER_T :: msgtag_resp_FGstart=         177
INTEGER_T :: msgtag_resp_FGslave=         178
INTEGER_T :: msgtag_resp_FGmaster=        179
INTEGER_T :: msgtag_resp_Fcoulstart=      180
INTEGER_T :: msgtag_resp_Fcoulslave=      181
INTEGER_T :: msgtag_resp_Fcoulmaster=     182
!NTEGER_T :: msgtag_eigs_entry      =     183
!NTEGER_T :: msgtag_main_integral   =     184
!NTEGER_T :: msgtag_cpks_constructs   =   185
!NTEGER_T :: msgtag_main_integral_2   =   186
!NTEGER_T :: msgtag_dealloc_gradxc    =   187
INTEGER_T :: msgtag_cpks_xc           =   188
INTEGER_T :: msgtag_cpks_solv1        =   189
INTEGER_T :: msgtag_cpks_solv2        =   190
INTEGER_T :: msgtag_cpks_solv3        =   191
INTEGER_T :: msgtag_cpks_solv4        =   192

INTEGER_T :: msgtag_resp_2center=         189
INTEGER_T :: msgtag_resp_3center_co=      190

INTEGER_T :: msgtag_tddft_eps_eta=        191
INTEGER_T :: msgtag_tddft_clshell=        192
INTEGER_T :: msgtag_tddft_opshell=        193
INTEGER_T :: msgtag_tddft_diagonl=        194
INTEGER_T :: msgtag_tddft_diagnup=        195
INTEGER_T :: msgtag_tddft_diagndn=        196
!NTEGER_T :: msgtag_tddft_parmlt1=        197
INTEGER_T :: msgtag_tddft_parmlt2=        198
INTEGER_T :: msgtag_tddft_parmlt3=        199
INTEGER_T :: msgtag_tddft_parmlt4=        200
INTEGER_T :: msgtag_tddft_dipolem=        201
INTEGER_T :: msgtag_tddft_diagsnd=        202
INTEGER_T :: msgtag_tddft_sendeig=        203
INTEGER_T :: msgtag_tddft_dealloM=        204

INTEGER_T :: msgtag_send_2c_Colmb=        205
INTEGER_T :: msgtag_send_3c_Colmb=        206
INTEGER_T :: msgtag_send_2c_start=        207
INTEGER_T :: msgtag_tmp_3co_send =        208
INTEGER_T :: msgtag_nori_2c_send =        209

INTEGER_T :: msgtag_ind_dipmom   =        210
INTEGER_T :: msgtag_Pol_ham      =        211
INTEGER_T :: msgtag_back_Pol_ham =        212

! entry point for legacy parallel eigensolver, slaves should
! call eigen_data_solve()
INTEGER_T :: msgtag_eigen_data_solve =    213
!NTEGER_T :: msgtag_eigs_entry      =     214

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
INTEGER_T :: msgtag_end=                  254
INTEGER_T :: msgtag_error=                255
!..............................................
!message tags >= 1000 are reserved for molecular
!mechanics part, see molmech_slave_module

public :: msgtag_name

!================================================================
! End of public interface of module
!================================================================

contains

  function msgtag_name(tag) result(name)
    implicit none
    integer(i4_kind), intent(in) :: tag
    character(len=32)            :: name
    ! *** end of interface ***

#define CMP(x) if(tag==x) name=__STRING(x)
    name = 'undef'
    CMP(msgtag_pert)
    CMP(msgtag_ham_calc_main)
    CMP(msgtag_dens)
    CMP(msgtag_finito)
    CMP(msgtag_eigensolve)
    !CMP(msgtag_coeff)
    CMP(msgtag_borders)
    !CMP(msgtag_corr)
    CMP(msgtag_commdata)
    !MP(msgtag_init)
    CMP(msgtag_debugorbitals)
    !MP(msgtag_bounds)
    CMP(msgtag_SendBackEig)
    CMP(msgtag_occ_levels)
    CMP(msgtag_eigvec)
    !MP(msgtag_gr_send)
    CMP(msgtag_xcham_send)
    CMP(msgtag_build_xc)
    !MP(msgtag_xc_setup)
    CMP(msgtag_SendBackHam)
    CMP(msgtag_xcfitham_send)
    CMP(msgtag_build_xcfit)
    !MP(msgtag_xcfit_setup)
    !CMP(msgtag_send_slavetiming)
    CMP(msgtag_slavetiming)
    !CMP(msgtag_grid_close)
    CMP(msgtag_int_shutdown)
    !CMP(msgtag_int_2cff_setup)
    !CMP(msgtag_int_2cff_shutdown)
    !CMP(msgtag_int_2cff_do)
    !CMP(msgtag_int_2cff_done)
    !CMP(msgtag_int_2cff_norm)
    CMP(msgtag_int_2cff_result)
    !CMP(msgtag_int_2cff_rep_back)
    !CMP(msgtag_int_2cob3c_setup)
    !CMP(msgtag_int_2cob3c_shutdown)
    !CMP(msgtag_int_2cob3c_do)
    !CMP(msgtag_int_2cob3c_done)
    !CMP(msgtag_int_2cob_norm)
    CMP(msgtag_int_2cob3c_result)
    !CMP(msgtag_int_2cob3c_rep_back)
    CMP(msgtag_int_rel_setup)
    CMP(msgtag_int_rel_shutdown)
    CMP(msgtag_int_2cob3c_filesize)
    !CMP(msgtag_post_scf)
    CMP(msgtag_ph_send)
    !CMP(msgtag_gridph_close)
    !CMP(msgtag_main_gradient)
    CMP(msgtag_grad_3c)
    CMP(msgtag_grad_ch)
    CMP(msgtag_occ_levels_eigval)
    CMP(msgtag_fit_coeff_send)
!   CMP(msgtag_fit_coeff_shutdown)
!   CMP(msgtag_xc_close)
!   CMP(msgtag_xcfit_close)
!   CMP(msgtag_openspecialunits)
!   CMP(msgtag_gradinfo_dealloc)
!   CMP(msgtag_eigvec_occ_dealloc)
!   CMP(msgtag_ua_close)
!   CMP(msgtag_op_dealloc)
!   CMP(msgtag_intstore_dealloc)
    CMP(msgtag_rel_Ttrans)
    CMP(msgtag_rel_Ttrans_back)
    CMP(msgtag_rel_momspace)
    CMP(msgtag_rel_momspace_back)
    CMP(msgtag_rel_send_overlap)
    CMP(msgtag_rel_trafomat)
    CMP(msgtag_SendBackTrafomat)
    CMP(msgtag_rel_gradient_send)
    CMP(msgtag_xcmdaham_send)
    CMP(msgtag_build_xcmda)
    !MP(msgtag_xcmda_setup)
    !MP(msgtag_pre_dens_slave)
    CMP(msgtag_vir_levels)
    CMP(msgtag_eigvec_vir_dealloc)
    CMP(msgtag_rot_levels)
    !MP(msgtag_xcmda_close)
    !CMP(msgtag_int_dipole_setup)
    !CMP(msgtag_int_dipole_shutdown)
    !CMP(msgtag_int_dipole_do)
    !CMP(msgtag_int_dipole_done)
    CMP(msgtag_int_dipole_result)
    !CMP(msgtag_int_dipole_rep_back)
    CMP(msgtag_response_setup)
    CMP(msgtag_response_level_send)
    CMP(msgtag_response_chmat)
    CMP(msgtag_response_2index)
    CMP(msgtag_resp_2in_sl_start)
    CMP(msgtag_response_2in_send)
    CMP(msgtag_response_2in_send2)
    CMP(msgtag_response_3index)
    CMP(msgtag_resp_sl_start_3ind)
    CMP(msgtag_resp_3ind_fitdim)
    CMP(msgtag_response_3ind_send)
    CMP(msgtag_response_3Clb_start)
    CMP(msgtag_resp_sl_start_3Clb)
    CMP(msgtag_response_3Clb_send)
    CMP(msgtag_resp_2center)
    CMP(msgtag_send_all_levels)
    CMP(msgtag_resp_3Clb_tmptape)
    CMP(msgtag_resp_3Clb_tmptape_2)
    CMP(msgtag_resp_3center_co)

    CMP(msgtag_tddft_eps_eta)
    CMP(msgtag_tddft_clshell)
    CMP(msgtag_tddft_opshell)
    CMP(msgtag_tddft_diagonl)
    CMP(msgtag_tddft_diagnup)
    CMP(msgtag_tddft_diagndn)
    !CMP(msgtag_tddft_parmlt1)
    CMP(msgtag_tddft_parmlt2)
    CMP(msgtag_tddft_parmlt3)
    CMP(msgtag_tddft_parmlt4)
    CMP(msgtag_tddft_dipolem)
    CMP(msgtag_tddft_diagsnd)
    CMP(msgtag_tddft_sendeig)
    CMP(msgtag_tddft_dealloM)

    CMP(msgtag_send_2c_Colmb)
    CMP(msgtag_send_3c_Colmb)
    CMP(msgtag_send_2c_start)

    CMP(msgtag_tmp_3co_send)
    CMP(msgtag_nori_2c_send)

!   CMP(msgtag_properties_main)
!   CMP(msgtag_eigen_data_send)
!   CMP(msgtag_free_eigen)
!   CMP(msgtag_alloc_eigen)
    CMP(msgtag_density_data_free)
!   CMP(msgtag_prescf_finalize)
    CMP(msgtag_grid_send_n_points)
    CMP(msgtag_grid_send_points)
    CMP(msgtag_int_setup)
    CMP(msgtag_space_point)
!   CMP(msgtag_bounds_poten)
    CMP(msgtag_start_poten)
    CMP(msgtag_finish_poten)
!   CMP(msgtag_prescf_init)
    CMP(msgtag_geom_grad)
!   CMP(msgtag_Q_elec)
    CMP(msgtag_solv_ham)
!   CMP(msgtag_back_ham_solv)
    CMP(msgtag_surf_point)
    CMP(msgtag_surf_point_sa)
    CMP(msgtag_sp_grinfo_dealloc)
    CMP(msgtag_start_field)
    CMP(msgtag_finish_field)
    CMP(msgtag_start_rho)
    CMP(msgtag_finish_rho)
!   CMP(msgtag_bounds_field)
    CMP(msgtag_free_bnds_ptn)
    CMP(msgtag_free_bnds_fld)
!   CMP(msgtag_del_poten)
!   CMP(msgtag_del_field)
    CMP(msgtag_epe_data_sent)
    CMP(msgtag_epe_do_gradients)
    CMP(msgtag_epe_grad_done)
    CMP(msgtag_epe_grads_init)
    CMP(msgtag_epe_grads_cons)
    CMP(msgtag_epe_forces_done)
    CMP(msgtag_epe_init_slave)
    CMP(msgtag_epe_finish_slave)
    CMP(msgtag_epe_defects)
    CMP(msgtag_epe_send_only)
    CMP(msgtag_vac_intermediate)
    CMP(msgtag_imp_intermediate)
    CMP(msgtag_epe_send_def)
    CMP(msgtag_epe_consdef)
    CMP(msgtag_epe_def_fin)
    CMP(msgtag_stop_daemon)
    CMP(msgtag_packed_message)
!   CMP(msgtag_filenames)
    CMP(msgtag_open_densmat)
    CMP(msgtag_response_rhs_h)
    CMP(msgtag_resp_sl_start_rhs)
    CMP(msgtag_resp_rhs_fitdim)
    CMP(msgtag_response_rhs_send)
    CMP(msgtag_resp_rho1dat)
    CMP(msgtag_resp_rho1dat_1b)
    CMP(msgtag_resp_all3index)
    CMP(msgtag_resp_all3Clb_start)
    CMP(msgtag_resp_addFcoul)
    CMP(msgtag_response_FG)
    CMP(msgtag_resp_FGstart)
    CMP(msgtag_resp_FGslave)
    CMP(msgtag_resp_FGmaster)
    CMP(msgtag_resp_Fcoulstart)
    CMP(msgtag_resp_Fcoulslave)
    CMP(msgtag_resp_Fcoulmaster)
!   CMP(msgtag_main_integral)
    CMP(msgtag_end)
    CMP(msgtag_error)
    CMP(msgtag_eigen_data_solve)

    if(name=='undef')then
       write(name,'(i5)') tag
    endif
  end function msgtag_name
!--------------- End of module ----------------------------------
end module msgtag_module
