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
!=====================================================================
! Public interface of module
!=====================================================================
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
!NTEGER_T :: msgtag_pert=                   2
INTEGER_T :: msgtag_finito=                 5
INTEGER_T :: msgtag_eigensolve=             6
!NTEGER_T :: msgtag_borders=                8
INTEGER_T :: msgtag_commdata=              11
!NTEGER_T :: msgtag_debugorbitals=         13
INTEGER_T :: msgtag_SendBackEig=           15
!NTEGER_T :: msgtag_eigvec=                17
!NTEGER_T :: msgtag_xcham_send=            19
!NTEGER_T :: msgtag_SendBackHam=           22
INTEGER_T :: msgtag_xcfitham_send=         23
INTEGER_T :: msgtag_slavetiming=           27
!NTEGER_T :: msgtag_int_shutdown=          31
INTEGER_T :: msgtag_int_2cff_result=       37
INTEGER_T :: msgtag_int_2cob3c_result=     44
!NTEGER_T :: msgtag_int_rel_setup=         46
!NTEGER_T :: msgtag_int_rel_shutdown=      47
INTEGER_T :: msgtag_int_2cob3c_filesize=   48
!NTEGER_T :: msgtag_ph_send=               51
INTEGER_T :: msgtag_grad_3c=               56
INTEGER_T :: msgtag_grad_ch=               58
!NTEGER_T :: msgtag_occ_levels_eigval=     59
INTEGER_T :: msgtag_fit_coeff_send=        60
!NTEGER_T :: msgtag_rel_Ttrans =           72
!NTEGER_T :: msgtag_rel_Ttrans_back =      73
!NTEGER_T :: msgtag_rel_momspace =         74
!NTEGER_T :: msgtag_rel_momspace_back =    75
!NTEGER_T :: msgtag_rel_send_overlap=      76
!NTEGER_T :: msgtag_rel_trafomat=          77
!NTEGER_T :: msgtag_SendBackTrafomat=      78
!NTEGER_T :: msgtag_rel_gradient_send=     80
INTEGER_T :: msgtag_xcmdaham_send=         81
INTEGER_T :: msgtag_build_xcmda=           82
INTEGER_T :: msgtag_vir_levels=            85
INTEGER_T :: msgtag_rot_levels=            87
INTEGER_T :: msgtag_int_dipole_result=     93
!NTEGER_T :: msgtag_response_setup=        95
INTEGER_T :: msgtag_response_level_send=   96
!NTEGER_T :: msgtag_response_chmat=        97
!NTEGER_T :: msgtag_response_2index=       98
!NTEGER_T :: msgtag_resp_2in_sl_start=     99
INTEGER_T :: msgtag_response_2in_send=    100
!NTEGER_T :: msgtag_response_2in_send2=   101
!NTEGER_T :: msgtag_response_3index=      102
!NTEGER_T :: msgtag_resp_sl_start_3ind=   103
!NTEGER_T :: msgtag_resp_3ind_fitdim=     104
!NTEGER_T :: msgtag_response_3ind_send=   105
!NTEGER_T :: msgtag_response_3Clb_start=  106
!NTEGER_T :: msgtag_resp_sl_start_3Clb=   107
INTEGER_T :: msgtag_response_3Clb_send=   108
!NTEGER_T :: msgtag_send_all_levels=      111
!NTEGER_T :: msgtag_resp_3Clb_tmptape=    112
!NTEGER_T :: msgtag_resp_3Clb_tmptape_2=  113
!NTEGER_T :: msgtag_grid_send_n_points =  122
!NTEGER_T :: msgtag_grid_send_points=     123
!NTEGER_T :: msgtag_int_setup=            124
INTEGER_T :: msgtag_space_point=          125
INTEGER_T :: msgtag_finish_poten=         128
INTEGER_T :: msgtag_geom_grad=            130
INTEGER_T :: msgtag_surf_point=           134
INTEGER_T :: msgtag_surf_point_sa=        135
INTEGER_T :: msgtag_sp_grinfo_dealloc=    136
INTEGER_T :: msgtag_start_field=          137
INTEGER_T :: msgtag_finish_field=         138
!NTEGER_T :: msgtag_start_rho=            139
!NTEGER_T :: msgtag_finish_rho=           140
INTEGER_T :: msgtag_free_bnds_fld=        143
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
INTEGER_T :: msgtag_stop_daemon    =      161
INTEGER_T :: msgtag_packed_message =      163
! added by HH in pg_V2.1_hh (was 124-139, changed to ):
!NTEGER_T :: msgtag_response_rhs_h=       167
!NTEGER_T :: msgtag_resp_sl_start_rhs=    168
!NTEGER_T :: msgtag_resp_rhs_fitdim=      169
!NTEGER_T :: msgtag_response_rhs_send=    170
!NTEGER_T :: msgtag_resp_rho1dat=         171
!NTEGER_T :: msgtag_resp_rho1dat_1b=      172
!NTEGER_T :: msgtag_resp_all3index=       173
!NTEGER_T :: msgtag_resp_all3Clb_start=   174
!NTEGER_T :: msgtag_resp_addFcoul=        175
!NTEGER_T :: msgtag_response_FG=          176
!NTEGER_T :: msgtag_resp_FGstart=         177
!NTEGER_T :: msgtag_resp_FGslave=         178
!NTEGER_T :: msgtag_resp_FGmaster=        179
!NTEGER_T :: msgtag_resp_Fcoulstart=      180
!NTEGER_T :: msgtag_resp_Fcoulslave=      181
!NTEGER_T :: msgtag_resp_Fcoulmaster=     182
!NTEGER_T :: msgtag_cpks_xc           =   188
INTEGER_T :: msgtag_cpks_solv1        =   189
INTEGER_T :: msgtag_cpks_solv2        =   190
INTEGER_T :: msgtag_cpks_solv3        =   191
INTEGER_T :: msgtag_cpks_solv4        =   192
INTEGER_T :: msgtag_resp_2center=         189
!NTEGER_T :: msgtag_resp_3center_co=      190
!NTEGER_T :: msgtag_tddft_eps_eta=        191
!NTEGER_T :: msgtag_tddft_clshell=        192
!NTEGER_T :: msgtag_tddft_opshell=        193
INTEGER_T :: msgtag_tddft_diagonl=        194
!NTEGER_T :: msgtag_tddft_diagnup=        195
!NTEGER_T :: msgtag_tddft_diagndn=        196
!NTEGER_T :: msgtag_tddft_parmlt2=        198
INTEGER_T :: msgtag_tddft_parmlt3=        199
!NTEGER_T :: msgtag_tddft_parmlt4=        200
INTEGER_T :: msgtag_tddft_dipolem=        201
!NTEGER_T :: msgtag_tddft_diagsnd=        202
INTEGER_T :: msgtag_tddft_sendeig=        203
!NTEGER_T :: msgtag_tddft_dealloM=        204
INTEGER_T :: msgtag_send_2c_Colmb=        205
!NTEGER_T :: msgtag_send_3c_Colmb=        206
!NTEGER_T :: msgtag_send_2c_start=        207
INTEGER_T :: msgtag_tmp_3co_send =        208
!NTEGER_T :: msgtag_nori_2c_send =        209
INTEGER_T :: msgtag_ind_dipmom   =        210
INTEGER_T :: msgtag_Pol_ham      =        211
INTEGER_T :: msgtag_back_Pol_ham =        212
!NTEGER_T :: msgtag_end=                  254
INTEGER_T :: msgtag_error=                255
!..............................................
!message tags >= 1000 are reserved for molecular
!mechanics part, see molmech_slave_module

public :: msgtag_name

  !===================================================================
  ! End of public interface of module
  !===================================================================

contains

  function msgtag_name(tag) result(name)
    implicit none
    integer(i4_kind), intent(in) :: tag
    character(len=32)            :: name
    ! *** end of interface ***

#define CMP(x) if(tag==x) name=__STRING(x)
    name = 'undef'
    !MP(msgtag_pert)
    CMP(msgtag_finito)
    CMP(msgtag_eigensolve)
    !MP(msgtag_borders)
    CMP(msgtag_commdata)
    !MP(msgtag_debugorbitals)
    CMP(msgtag_SendBackEig)
    !MP(msgtag_eigvec)
    !MP(msgtag_xcham_send)
    !MP(msgtag_SendBackHam)
    CMP(msgtag_xcfitham_send)
    CMP(msgtag_slavetiming)
    !MP(msgtag_int_shutdown)
    CMP(msgtag_int_2cff_result)
    CMP(msgtag_int_2cob3c_result)
    !MP(msgtag_int_rel_setup)
    !MP(msgtag_int_rel_shutdown)
    CMP(msgtag_int_2cob3c_filesize)
    !MP(msgtag_ph_send)
    CMP(msgtag_grad_3c)
    CMP(msgtag_grad_ch)
    !MP(msgtag_occ_levels_eigval)
    CMP(msgtag_fit_coeff_send)
    !MP(msgtag_rel_Ttrans)
    !MP(msgtag_rel_Ttrans_back)
    !MP(msgtag_rel_momspace)
    !MP(msgtag_rel_momspace_back)
    !MP(msgtag_rel_send_overlap)
    !MP(msgtag_rel_trafomat)
    !MP(msgtag_SendBackTrafomat)
    !MP(msgtag_rel_gradient_send)
    CMP(msgtag_xcmdaham_send)
    CMP(msgtag_build_xcmda)
    CMP(msgtag_vir_levels)
    CMP(msgtag_rot_levels)
    CMP(msgtag_int_dipole_result)
    !MP(msgtag_response_setup)
    CMP(msgtag_response_level_send)
    !MP(msgtag_response_chmat)
    !MP(msgtag_response_2index)
    !MP(msgtag_resp_2in_sl_start)
    CMP(msgtag_response_2in_send)
    !MP(msgtag_response_2in_send2)
    !MP(msgtag_response_3index)
    !MP(msgtag_resp_sl_start_3ind)
    !MP(msgtag_resp_3ind_fitdim)
    !MP(msgtag_response_3ind_send)
    !MP(msgtag_response_3Clb_start)
    !MP(msgtag_resp_sl_start_3Clb)
    CMP(msgtag_response_3Clb_send)
    !MP(msgtag_send_all_levels)
    !MP(msgtag_resp_3Clb_tmptape)
    !MP(msgtag_resp_3Clb_tmptape_2)
    !MP(msgtag_grid_send_n_points)
    !MP(msgtag_grid_send_points)
    !MP(msgtag_int_setup)
    CMP(msgtag_space_point)
    CMP(msgtag_finish_poten)
    CMP(msgtag_geom_grad)
    CMP(msgtag_surf_point)
    CMP(msgtag_surf_point_sa)
    CMP(msgtag_sp_grinfo_dealloc)
    CMP(msgtag_start_field)
    CMP(msgtag_finish_field)
    !MP(msgtag_start_rho)
    !MP(msgtag_finish_rho)
    CMP(msgtag_free_bnds_fld)
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
    !MP(msgtag_response_rhs_h)
    !MP(msgtag_resp_sl_start_rhs)
    !MP(msgtag_resp_rhs_fitdim)
    !MP(msgtag_response_rhs_send)
    !MP(msgtag_resp_rho1dat)
    !MP(msgtag_resp_rho1dat_1b)
    !MP(msgtag_resp_all3index)
    !MP(msgtag_resp_all3Clb_start)
    !MP(msgtag_resp_addFcoul)
    !MP(msgtag_response_FG)
    !MP(msgtag_resp_FGstart)
    !MP(msgtag_resp_FGslave)
    !MP(msgtag_resp_FGmaster)
    !MP(msgtag_resp_Fcoulstart)
    !MP(msgtag_resp_Fcoulslave)
    !MP(msgtag_resp_Fcoulmaster)
    !MP(msgtag_cpks_xc)
    CMP(msgtag_cpks_solv1)
    CMP(msgtag_cpks_solv2)
    CMP(msgtag_cpks_solv3)
    CMP(msgtag_cpks_solv4)
    CMP(msgtag_resp_2center)
    !MP(msgtag_resp_3center_co)
    !MP(msgtag_tddft_eps_eta)
    !MP(msgtag_tddft_clshell)
    !MP(msgtag_tddft_opshell)
    CMP(msgtag_tddft_diagonl)
    !MP(msgtag_tddft_diagnup)
    !MP(msgtag_tddft_diagndn)
    !MP(msgtag_tddft_parmlt2)
    CMP(msgtag_tddft_parmlt3)
    !MP(msgtag_tddft_parmlt4)
    CMP(msgtag_tddft_dipolem)
    !MP(msgtag_tddft_diagsnd)
    CMP(msgtag_tddft_sendeig)
    !MP(msgtag_tddft_dealloM)
    CMP(msgtag_send_2c_Colmb)
    !MP(msgtag_send_3c_Colmb)
    !MP(msgtag_send_2c_start)
    CMP(msgtag_tmp_3co_send)
    !MP(msgtag_nori_2c_send)
    CMP(msgtag_ind_dipmom)
    CMP(msgtag_Pol_ham)
    CMP(msgtag_back_Pol_ham)
    !MP(msgtag_end)
    CMP(msgtag_error)

    if(name=='undef')then
       write(name,'(i5)') tag
    endif
  end function msgtag_name
!--------------- End of module ----------------------------------
end module msgtag_module
