#-*- makefile -*-
#
# ParaGauss, a program package for high-performance computations
# of molecular systems
# Copyright (C) 2014
# T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
# M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
# A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
# T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
# M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
# M. Roderus, N. Rösch
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 as published
# by the Free Software Foundation [1].
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# [1] http://www.gnu.org/licenses/gpl-2.0.html
#
# Please see the accompanying LICENSE file for further information.
#
# Makefile for new ttfs (lcgto), Thomas Belling 2/95
#
# AM, 31.01.2004:
#
# 1) Before build, find suitable settings in
#
#	./make/*.mk
#
#    directory. Either copy, or create a symbolic link:
#
# 	machine.mk -> make/HOST.mk
#
# 2) Build the rest:
#
# 	make -r
#
#    Parallel make also works:
#
#       make -srkjN
#
#    where N is  the maximum job number (any  if omitted).  Hint: some
#    compilers require licences to run every instance (regatta).
#
# The dependencies  are build automatically, but they  are not checked
# for loops, you may issue
#
#   make chkdepend
#
# to do that from time to time.
#

# Set default  to ``build'' in  order not to overwrite  executables in
# the  installation  directory.   Then  you'll need  to  issue  ``make
# install'' manually
default: build

# Version as  printed in  the legal  header (and as  part of  the file
# name, occasionally):
MAJOR =  3
MINOR = .2
# Maintenance number OR branch abbreviation (e.g. your initials):
EXTRA = .1
# Change or increment this:
BUILD =

VERS = V$(MAJOR)$(MINOR)$(EXTRA)$(BUILD)

# version for ttfs, read_operations, ( and get_mpi_config ):
# scriptvers := $(VERS)
scriptvers := V$(MAJOR)$(MINOR)

#
# The *printed* version of ParaGauss.  Use no spaces, or you will deal
# with shell-escape sequences:
#
paragauss_vers = $(VERS)-$(shell date +%Y%m%d)


# compile optimizer into PG:
WITH_OPTIMIZER = 1

#
# Guile is a Scheme implementation intended as a GNU extension language:
#
WITH_GUILE = 1

#
# Use memlog facility:
#
WITH_MEMLOG = 0

# yet another integral package
WITH_SHGI = 1

# use new integral implementation
NEW_INTEGRALS = 1

# compile 2nd derivatives (req NEW_INTEGRALS):
WITH_SECDER = 1

#
# Use    parallel   eigensolver    with   scheduling    (called   from
# eigen_data_module.f90).  If  you want  to use this  function, please
# make sure that a correct  cost function for your machine is provied.
# To   see   how   to   generate   a  cost   function,   please   read
# schedeig/se_runrec/README.    When   set,  WITH_SCHEDEIG   Overrides
# WITH_MATRIX_PARALLEL    for     eigensolver    (see    #ifdefs    in
# eigen_data_module.f90):
#
WITH_SCHEDEIG = 0

#
# Use  PBLAS,   ScaLAPACK  for  relativistic   transformation  of  the
# hamiltonian  in  pre-SCF  (relgrads.f90)  and  the  solver  for  the
# generalized   eigenvalue   problem   (eigen_data_module.f90).    The
# eigensolver option is overriden by WITH_SCHEDEIG = 1:
#
WITH_MATRIX_PARALLEL = 0

#
# Compile EPE into PG:
#
WITH_EPE = 1
NEW_EPE = 1
ifeq ($(WITH_EPE),0)
	NEW_EPE = 0
endif

#
# Compile solvation into PG:
#
WITH_SOLV = 1

#
# If you want  to use 2nd derivatives within  solvation part.  That is
# just    temporary,    as    now    solvation   cannot    run    with
# "no_cpks_coul_grads" variant. The bug will be fixed a bit later.  On
# the other hand,  without this switch set allocation  of the 3-center
# gradients will become a bottleneck:
#
no_cpks_coul_grads=1

#
# Compile Molecular Mechanics into PG:
#
WITH_MOLMECH = 1

#
# Compile Effective Fragment Potential method into PG:
#
WITH_EFP = 0
ifeq ($(WITH_MOLMECH),0)
	WITH_EFP = 0
endif

#
# Compile response module into PG:
#
WITH_RESPONSE = 1

#
# Compile GTEN and HFCC module into PG:
#
WITH_GTENSOR = 1

#
# Core density intended for use with pseudopotentials (incomplete):
#
WITH_CORE_DENS = 1

#
# For debug use CCP DFT repository:
#
WITH_LIBDFTAUTO = 0

#
# Use eri4c library for 4-center integrals:
#
WITH_ERI4C = 0

#
# Use DFT+U and DFT+Umol:
#
WITH_DFTPU = 0

#
# Compile the old (pre-hitachi) input pre-processing:
#
WITH_OLD_INPUT = 1

#
# If checks for NaNs and Infs are available:
#
WITH_ISNAN = 1

#
# Read/write  scfcontrol file for  interactive control  of convergence
# parameters:
#
WITH_SCFCONTROL = 0

#
# Compile experimental code:
#
WITH_EXPERIMENTAL = 0

#
#
# BGY3D  is linked  in  just as  any  other library.   It  is not  the
# responsibility of (this) Make to  assemble the library. So beware of
# version incompatibilities:
#
WITH_BGY3D = 0

#
# Any of these require ScaLAPACK and BLACS libraries:
#
ifeq ($(WITH_SCHEDEIG),1)
	WITH_SCALAPACK = 1
endif
ifeq ($(WITH_MATRIX_PARALLEL),1)
	WITH_SCALAPACK = 1
endif

#
# Include path for FORTRAN files/module search path will contain these
# directoris unconditionally, first  those relative to the $(BASEDIR).
# $(COMMDIR) is  set in  machine.mk depending on  the build  type. The
# list will be augmented, according to the options set above:
#
dirs =  modules \
        $(COMMDIR) \
        epe_dir \
        solvation \
        effect_frag_poten \
        molmech \
        molmech/solvation \
        external_centers \
        elec_pot_and_field \
        optimizer \
        lapack \
        gaussq \
        cpks \
        autodiff \
        libdftauto \
        punchfile

# misc options (debug_targets) here:
-include misc.mk

#### Machine-dependent stuff here:
include machine.mk
machine.mk:
	@echo -e \
	 "\n=============================================" \
	 "\n=== YOU SEEM TO HAVE FORGOTTEN TO CREATE  ===" \
	 "\n=== A LINK: ln -s make/HOST.mk machine.mk ===" \
	 "\n=============================================" \
	 "\n"
	@false

# for makedependent: OBJS MODOBJS EISOBJS COMMOBJS

#### OBJECTS ####

libttfs_callcomm.a += \
		assert_failed.o \
		trace_hook.o \
		chargefit.o \
		error_handler.o \
		elec_pot_and_field/field_calculate.o \
		integral_calc_quad_2cff.o \
		integral_calc_quad_2cob3c.o \
		integral_calc_quad_dipole.o \
		modules/integral_calc_quad_module.o \
		integral_interrupt_2cob3c.o \
		integral_main_2cff.o \
		integral_main_2cob3c.o \
		integral_main_dipole.o \
		integral_setup_2cff.o \
		integral_setup_2cob3c.o \
		integral_setup_dipole.o \
		integral_setup.o \
		integral_shutdown_2cob3c.o \
		integral_shutdown.o \
		main_gradient.o \
		main_integral.o \
		main_master.o \
		modules/paragauss.o \
		main.o \
		main_scf.o \
		main_slave.o \
		elec_pot_and_field/potential_calculate.o \
		write_input.o \
		modules/initialization.o \
		modules/back_trafo_module.o \
		modules/bounds_module.o \
		modules/casc_logic_module.o \
		modules/clebsch_gordan.o \
		modules/density_data_module.o \
		modules/eigen_data_module.o \
		elec_pot_and_field/elec_static_field_module.o \
		modules/filename_module.o \
		modules/fit_coeff_module.o \
		modules/gradient_data_module.o \
		modules/grid_module.o \
		modules/ham_calc_module.o \
		modules/int_distribute_module.o \
		modules/integralpar_module.o \
		modules/int_send_2cff_module.o \
		modules/int_send_2cob3c_module.o \
		modules/int_send_2cob3c_spor_module.o \
		modules/int_send_2cob3c_spor.o \
		modules/int_send_aux_module.o \
		modules/int_send_dipole_module.o \
		modules/machineparameters_module.o \
		modules/occupied_levels_module.o \
		modules/operations_module.o \
		modules/options_module.o \
		modules/orbital_plot_module.o \
		tddft/nto_plot_module.o \
		modules/orbitalprojection_module.o \
		modules/output_module.o \
		modules/pert_coeff_module.o \
		modules/ph_cntrl.o \
		external_centers/pointcharge_module.o \
		external_centers/point_dqo_module.o \
		external_centers/induced_dipoles_module.o \
		external_centers/calc_id_module.o \
		modules/post_scf_module.o \
		elec_pot_and_field/potential_calc_module.o \
		elec_pot_and_field/potential_module.o \
		modules/prepare_integralfiles_module.o \
		modules/prescf_module.o \
		modules/properties_module.o \
		modules/quadrupel_fname.o \
		modules/quadrupel_module.o \
		modules/reltrafo.o \
		modules/relgrads.o \
		modules/relgrads_store.o \
		solvation/solv_cavity_module.o \
		solvation/solv_electrostat_module.o \
		modules/symm_adapt_module.o \
		modules/symm_adapt_xpack.o \
		modules/symmetry_data_module.o \
		modules/time_module.o \
		modules/timer_module.o \
		modules/unique_atom_methods.o \
		modules/virtual_levels_module.o \
		modules/xcfit_hamiltonian.o \
		modules/xc_hamiltonian.o \
		modules/xcmda_hamiltonian.o \
		modules/xpack.o \
		modules/back_trafo_tapes.o \
		modules/pseudo_potential_module.o

libttfs_dontcallcomm.a += \
		gengrp.o \
		solvation/grad_solv_calculate.o \
		integral_trafo.o \
		lapack/f77_blas.o \
		lapack/f77_lapack.o \
		lapack/f77_scalapack.o \
		main_dipole.o \
		modules/symmetry.o \
		precision_check.o \
		read_input.o \
		autodiff/becke_step_func.o \
		modules/atoms_data_module.o \
		modules/baerends_module.o \
		modules/becke_perdew_module.o \
		modules/bessel_module.o \
		modules/calc3c_switches.o \
		modules/contraction_module.o \
		modules/convergence_module.o \
		modules/datatype.o \
		modules/density_calc_module.o \
		modules/density_calc_cpks.o \
		modules/dimensions.o \
		modules/dipole_module.o \
		modules/dip_prim_module.o \
		solvation/build_mol_surfaces.o \
		solvation/energy_and_grad_of_cavity.o \
		solvation/disp_rep_wrap.o \
		solvation/polyhedron_module.o \
		solvation/solv_charge_mixing_module.o \
		solvation/disp_rep_module.o \
		solvation/help_cavity_module.o \
		solvation/solv_2nd_deriv_module.o \
		solvation/cavity_image_module.o \
		modules/echo_input_module.o \
		modules/efield_module.o \
		modules/efm_decl.o \
		modules/efm_module.o \
		modules/energy_calc_module.o \
		modules/error_module.o \
		modules/fermi_module.o \
		modules/fitcontract_2c_grad_module.o \
		modules/fitcontract_2c_module.o \
		modules/fitcontract_module.o \
		modules/fit_trafo_module.o \
		modules/fit_trafo_tapes.o \
		modules/frag_orb_analysis_module.o \
		modules/gamma_module.o \
		modules/gradient_2c_fit_ch_module.o \
		modules/group_module.o \
		modules/hamiltonian_module.o \
		modules/init_module.o \
		modules/input_module.o \
		modules/int_data_2cff_module.o \
		modules/int_data_2cob3c_module.o \
		modules/int_data_dipole_module.o \
		modules/integerstack_module.o \
		modules/integral_2c_fit_ch_module.o \
		modules/integral_2c_fit_xc_module.o \
		modules/istore.o \
		modules/integralstore_module.o \
		modules/interfaces.o \
		modules/io.o \
		modules/iounitadmin_module.o \
		modules/linsys_module.o \
		modules/ll_calculate_grads_module.o \
		modules/mat_charge_module.o \
		modules/math_module.o \
		modules/matrix_check.o \
		modules/matrix_eigenval.o \
		modules/matrix_linsolve.o \
		modules/matrix_matmult.o \
		modules/matrix_methods.o \
		modules/matrix_functions.o \
		modules/matrix_module.o \
		modules/matrix_sparse.o \
		modules/matrix_types.o \
		modules/matrix_xpack.o \
		modules/mixing_module.o \
		modules/msgtag_module.o \
		modules/occupation_module.o \
		modules/orbital_module.o \
		modules/orbitalstore_module.o \
		modules/overlap_module.o \
		modules/pairs_module.o \
		modules/pbe_ggcxc_module.o \
		modules/perdew_wang_module.o \
		modules/population_module.o \
		modules/print_module.o \
		modules/pw_ldac_module.o \
		modules/quadrupelstore_module.o \
		modules/readwriteblocked_module.o \
		modules/rspace_plot.o \
		modules/solhrules_module.o \
		modules/solid_harmonics_module.o \
		modules/spectrum_module.o \
		modules/spin_orbit_module.o \
		modules/strings.o \
		modules/prim_int_store.o \
		modules/symm_adapt_int.o \
		modules/symm_adapt_struct.o \
		modules/symmetry_element.o \
		modules/symm_module.o \
		modules/symm_positions.o \
		modules/type_module.o \
		modules/uatom_symmadapt.o \
		modules/unique_atom_module.o \
		modules/vwnc.o \
		modules/xc_cntrl.o \
		modules/xc_ham_trafo.o \
		optimizer/eigensolver.o \
		optimizer/gxfile.o \
		modules/s2_expect.o \
		modules/relxc.o \
		modules/hcth.o \
		modules/debug.o \
		modules/xc_func.o \
		modules/exchange.o \
		modules/vdw_dft.o \
		modules/cpksdervs_matrices.o \
		modules/ch_response_module.o \
		modules/diis_fock_module.o \
		modules/gga_response_module.o \
		modules/m06.o \
		modules/lyp.o \
		modules/vsxc_mgga_module.o \
		modules/tpss_mgga_module.o \
		modules/tpss.o \
		modules/pbe_gga_module.o \
		modules/pw_lda.o

libttfs_data.a = \
		optimizer/atom_data_module.o \
		modules/constants.o \
		gaussq/gaussq_data.o

libttfs_pun.a = punchfile/punchfile.o

libttfs_epe.a = \
		main_epe_block.o \
		symm_epe.o \
		symm_ewpc_gen.o \
		epe_dir/lin_search_epe.o \
		epe_dir/epe_driver.o \
		epe_dir/epe_lattice_optimization.o \
		epe_dir/main_epe_module.o \
		epe_dir/atoms_parameters_module.o \
		epe_dir/culon_module.o \
		epe_dir/epecom_module.o \
		epe_dir/epepar_module.o \
		epe_dir/epe_pg_module.o \
		epe_dir/epe_pot_module.o \
		epe_dir/minepe_module.o \
		epe_dir/mol_module.o \
		epe_dir/str_module.o \
		modules/epe_module.o \
		modules/ewaldpc_module.o

ifeq ($(NEW_EPE),1)
	libttfs_epe.a += \
		epe_dir/read_epeinput.o \
		epe_dir/qm_epe_interface_module.o

endif

ifeq ($(NEW_EPE),0)
	libttfs_epe.a += \
		epe_dir/read_epe_input.o
endif

libttfs_efp.a = \
		effect_frag_poten/efp_data_module.o \
		effect_frag_poten/efp_module.o \
		effect_frag_poten/efp_efp_module.o \
		effect_frag_poten/efp_polar_module.o \
		effect_frag_poten/efp_rep_module.o \
		effect_frag_poten/efp_only_opt_module.o \
		effect_frag_poten/efp_solv_module.o \
		effect_frag_poten/efp_solv_grad_module.o

libttfs_molmech.a = \
		molmech/qmmm_interface_module.o \
		molmech/string_qmmm_module.o \
		molmech/qmmm1_interface_module.o \
		molmech/inp_out_module.o \
		molmech/tasks_main_options_module.o \
		molmech/common_data_module.o \
		molmech/species_module.o \
		molmech/element_data_module.o \
		molmech/n_body_lists_module.o \
		molmech/potentials_module.o \
		molmech/energy_and_forces_module.o \
		molmech/covalent_module.o \
		molmech/van_der_waals_module.o \
		molmech/coulomb_module.o \
		molmech/calc_energy_module.o \
		molmech/hess_and_opt_module.o \
		molmech/read_molmech_input.o \
		molmech/slab_module.o \
		molmech/ewald2d_module.o \
		molmech/ewald_module.o \
		molmech/molmech_msgtag_module.o \
		molmech/molmech_slave_module.o \
		molmech/mm_timer_module.o \
		molmech/external_field_module.o \
		molmech/pc_array_module.o \
		molmech/solvation/vdwcm_module.o \
		molmech/solvation/cavity_module.o \
                molmech/solvation/solv_elec_stat_module.o \
                molmech/solvation/ewald_solv_module.o \
                molmech/solvation/energy_and_grad_of_cavity_mm.o \
                molmech/solvation/disp_rep_wrap_mm.o \
                molmech/solvation/calc_solv_eff.o \
		molmech/main_molmech.o 

libttfs_opt.a = optimizer/opt_data_module.o \
		optimizer/hesse_module.o \
		optimizer/coordinates_module.o \
		optimizer/coortype_module.o \
		optimizer/frequency_module.o \
		optimizer/geo_operations_module.o \
		optimizer/ts_module.o \
		optimizer/valence_coord_module.o \
		optimizer/step_module.o \
		optimizer/line_search_module.o \
		optimizer/gradient_module.o \
		optimizer/slspar_module.o \
		optimizer/allocopt_module.o \
		optimizer/optimizer.o \
                optimizer/vff_hessian.o \
                optimizer/thermodyn_prop_module.o \

libttfs_shgi.a = shgi/shgi.o \
		shgi/shgi_cntrl.o \
		shgi/shgi_common.o \
		shgi/shgi_utils.o \
		shgi/shgi_dnf.o \
		shgi/shgi_rad.o \
		shgi/shgi_shr.o \
		shgi/shgi_ang.o \
		shgi/shgi_ab.o \
		shgi/shgi_sym.o \
		shgi/shgi_pseudo.o \
		shgi/shgi_relfit.o \
		shgi/shgi_relnuc.o \
		shgi/shgi_pcm.o \
		shgi/shgi_slv.o \
		shgi/shgi_ep_ef.o \
		shgi/shgi_adkh.o \
		shgi/shgi_ext_c.o \
		shgi/shgi_dip.o \

libttfs_resp.a = \
		modules/response_module.o \
		modules/int_send_2c_resp.o \
		modules/int_send_3c_resp.o \
		modules/noRI_module.o \
		tddft/dvdson_module.o \
		tddft/eigensolve_module.o \
		tddft/global_module.o \
		tddft/read_module.o \
		tddft/eigenblock_module.o \
		tddft/linalg_module.o \
		tddft/phys_param_module.o \
		tddft/result_module.o \
		tddft/init_tddft_module.o \
		tddft/tddft_diag.o \
		tddft/resp_util_module.o \
		tddft/resp_dipole_module.o \
		modules/int_resp_module.o \
		tddft/lan_solve_module.o \
		tddft/laso/dilaso.o\
		tddft/laso/dippla.o\
		tddft/laso/diwla.o\
		tddft/laso/dlabax.o\
		tddft/laso/dlabcm.o\
		tddft/laso/dlabfc.o\
		tddft/laso/dlaeig.o\
		tddft/laso/dlager.o\
		tddft/laso/dlaran.o\
		tddft/laso/dmvpc.o\
		tddft/laso/dnlaso.o\
		tddft/laso/dnppla.o\
		tddft/laso/dortqr.o\
		tddft/laso/dvsort.o\
		tddft/laso/urand.o\
		tddft/nto_module.o\
		cpks/cpks_grid_utils.o

libttfs_gten.a = \
		modules/gtensor_module.o \
		modules/hfc_module.o \
		ll_calculate_dipoleg.o \
		ll_calculate_hfc.o

libttfs_secder.a = \
		cpks_g4constructs.o \
		cpks/cpks_common.o \
		cpks/cpks_utils.o \
		calc_cpks_gvec.o \
		calc_cpksQai.o \
		calc_cpks_h1imp.o

ifeq ($(WITH_SCHEDEIG),1)
	dirs += schedeig
	FPPOPTIONS += -DWITH_SCHEDEIG
	libttfs_callcomm.a += \
		schedeig/se_eigen_module.o \
		schedeig/se_scheduling_module.o \
		schedeig/se_timefunction_module.o
endif

ifeq ($(WITH_MATRIX_PARALLEL),1)
	FPPOPTIONS += -DWITH_MATRIX_PARALLEL
	libttfs_callcomm.a += \
		modules/matrix_parallel.o
endif

ifeq ($(NEW_INTEGRALS),1)
	libttfs_dontcallcomm.a += \
		calc_3center.o \
		calc_radial_r2.o \
		calc_radial_3cGa.o \
                calc_radial_3cSA.o \
		calc_radial_nucp.o \
		calc_radial_nucGa.o \
		calc_radial_potG.o \
		modules/calc_3center_module.o \
		calc_3c_colc_setup.o \
		calc_3c_codervs_setup.o \
		calc_3c_fitcontract.o \

	FPPOPTIONS += -DNEW_INTEGRALS
else
	libttfs_dontcallcomm.a += \
		ss_calculate_grads.o \
		ls_calculate_grads.o \
		ll_calculate_grads.o \
		modules/ls_calculate_grads_module.o
endif

ifeq ($(WITH_ERI4C),1)
	libttfs_callcomm.a += \
		eri/eri_main.o \
		eri/eri_block.o \
		eri/eri_schedule.o \
		eri/eri_conversion.o \
		modules/direct_scf.o
	libttfs_dontcallcomm.a += \
		eri/eri_types.o \
		eri/eri_auxiliary.o \
	FPPOPTIONS += -DWITH_ERI4C
endif

ifeq ($(WITH_MEMLOG),1)
	libttfs_dontcallcomm.a += modules/memlog_module.o
	FPPOPTIONS += -DWITH_MEMLOG
endif

ifeq ($(WITH_OLD_INPUT),1)
	FPPOPTIONS += -DWITH_OLD_INPUT
endif

ifeq ($(WITH_ISNAN),1)
	FPPOPTIONS += -DWITH_ISNAN
	FPPOPTIONS += -DWITH_ISINF
endif

ifeq ($(WITH_NANINFCHK),1)
	libttfs_dontcallcomm.a += naninfchk.o
	FPPOPTIONS += -DFPP_HAVE_NANINFCHK

# for absoft only:
naninfchk.o: FFLAGS += -YBOZTYPE=INT
endif

#
# COMMOBJS/COMMCOBJS for communication (message passing):
#
# less if serial:
COMMOBJS =	$(COMMDIR)/comm.o \
		$(COMMDIR)/comm_module.o \
		$(COMMDIR)/commpack_module.o \
		$(if $(serial),,\
		$(COMMDIR)/commparameter_module.o \
		)

# nothing if serial:
COMMCOBJS =	$(if $(serial),,\
		$(COMMCDIR)/mpi_comm.o \
		$(COMMCDIR)/mpipack.o \
		$(COMMCDIR)/comm_variables.o \
		$(COMMCDIR)/errout.o \
		)

EISOBJS =	eis/ch.o \
		eis/htribk.o \
		eis/htridi.o \
		eis/rebak.o \
		eis/reduc.o \
		eis/rs.o \
		eis/rsg.o \
		eis/scopy.o \
		eis/sscal.o \
		eis/tql2.o \
		eis/tqlrat.o \
		eis/tred1.o \
		eis/tred2.o


DOCFILES = $(OBJS:.o=.doc) $(MODOBJS:.o=.doc)

# f90objs  = $(OBJS) $(MODOBJS) $(EISOBJS) $(COMMOBJS)
f90objs = \
	$(libttfs_callcomm.a) \
	$(libttfs_dontcallcomm.a) \
	$(libttfs_pun.a) \
	$(libttfs_data.a) \
	$(libttfs_pun.a) \
	$(COMMOBJS) \
	$(EISOBJS) 

ttfs_libs = \
	libttfs_comm.a \
	libttfs_commc.a \
	libttfs_eis.a \
	libttfs_callcomm.a \
	libttfs_dontcallcomm.a \
	libttfs_data.a \
	libttfs_pun.a

ifeq ($(WITH_OPTIMIZER),1)
	ttfs_libs += libttfs_opt.a
	extra_ttfs_libflags += -lttfs_opt
	f90objs += $(libttfs_opt.a)
	FPPOPTIONS += -DWITH_OPTIMIZER
endif

ifeq ($(WITH_EPE),1)
	ttfs_libs += libttfs_epe.a
	extra_ttfs_libflags += -lttfs_epe
	f90objs += $(libttfs_epe.a)
	FPPOPTIONS += -DWITH_EPE
	ifeq ($(WITH_MOLMECH),0)
		# FIXME: molmech/inp_out_module is used also in epe,
		#        so append it in case MOLMECH happens to be
		#        disabled
		f90objs += molmech/inp_out_module.o
		libttfs_epe.a += molmech/inp_out_module.o
	endif
endif

ifeq ($(NEW_EPE),1)
	FPPOPTIONS += -DNEW_EPE
endif

ifeq ($(no_cpks_coul_grads),1)
	FPPOPTIONS += -Dno_cpks_coul_grads
endif

ifeq ($(WITH_EFP),1)
	ttfs_libs += libttfs_efp.a
	extra_ttfs_libflags += -lttfs_efp
	f90objs += $(libttfs_efp.a)
	FPPOPTIONS += -DWITH_EFP
endif

ifeq ($(WITH_MOLMECH),1)
	ttfs_libs += libttfs_molmech.a
	extra_ttfs_libflags += -lttfs_molmech -lttfs_callcomm
	f90objs += $(libttfs_molmech.a)
	FPPOPTIONS += -DWITH_MOLMECH
endif

ifeq ($(WITH_SHGI),1)
        dirs += shgi
	ttfs_libs += libttfs_shgi.a
	extra_ttfs_libflags += -lttfs_shgi
	f90objs += $(libttfs_shgi.a)
	FPPOPTIONS += -DWITH_SHGI
endif

ifeq ($(WITH_RESPONSE),1)
        dirs += tddft
	ttfs_libs += libttfs_resp.a
	extra_ttfs_libflags += -lttfs_resp
	f90objs += $(libttfs_resp.a)
	FPPOPTIONS += -DWITH_RESPONSE
endif

ifeq ($(WITH_GTENSOR),1)
	ttfs_libs += libttfs_gten.a
	extra_ttfs_libflags += -lttfs_gten
	f90objs += $(libttfs_gten.a)
	FPPOPTIONS += -DWITH_GTENSOR
endif

ifeq ($(WITH_CORE_DENS),1)
	FPPOPTIONS += -DWITH_CORE_DENS
endif

ifeq ($(WITH_SECDER),1)
ifeq ($(NEW_INTEGRALS),0)
$(error WITH_SECDER requires NEW_INTEGRALS )
endif
	ttfs_libs += libttfs_secder.a
	extra_ttfs_libflags += -lttfs_secder
	f90objs += $(libttfs_secder.a)
	FPPOPTIONS += -DWITH_SECDER
ifeq ($(WITH_EXPERIMENTAL),1)
	experimental_objs = \
			cpks/cpks_grid_utils.o \
			cpks/cpks_xc_resp.o \
			cpks/cpks_dens_calc.o
	libttfs_secder.a += $(experimental_objs)
	f90objs += $(experimental_objs)
	FPPOPTIONS += -DWITH_EXPERIMENTAL
endif

ifeq ($(WITH_SCFCONTROL),1)
	FPPOPTIONS += -DWITH_SCFCONTROL
endif
endif

# # # # # # # # # # # # # LIBDFTAUTO SECTION # # # # # # # # # # # # #
ifeq ($(WITH_LIBDFTAUTO),1)
	ttfs_libs += libdftauto/libdftauto.a
	extra_ttfs_libflags += -L./libdftauto -ldftauto
	FPPOPTIONS += -DWITH_LIBDFTAUTO

FORCE:
libdftauto/libdftauto.a: FORCE
	$(MAKE) -C libdftauto

libdftauto_clean:
	$(MAKE) -C libdftauto clean

clean += libdftauto_clean
endif

# # # # # # # # # # # # # AUTODIFF SECTION # # # # # # # # # # # # # #

#
# Makefile  for the  AD lib,  set $(AD)  before include,  This defines
# $(libad.a) target among ather things:
#
AD = ./autodiff
include $(AD)/module.mk

#
# These  are supposed  to  be  the list  of  dependencies, and  linker
# options, respectively:
#
ttfs_libs += $(libad.a)
extra_ttfs_libflags += -L$(AD) -lad

# # # # # # # # # # Dynamic Load Balancing (DLB) # # # # # # # # # # #

#
# Makefile for  the DLB lib,  set $(DLB) before include,  This defines
# $(libdlb.a) target among ather things:
#
DLB = ./dlb
include $(DLB)/module.mk

#
# These  are supposed  to  be  the list  of  dependencies, and  linker
# options, respectively:
#
dirs += $(DLB)
ttfs_libs += $(libdlb.a)
extra_ttfs_libflags += -L$(DLB) -ldlb

# # # # # Electron Repulsion Integrals on 4 Center (ERI4C) # # # # # #

ifeq ($(WITH_ERI4C),1)
ERI4C = ./eri4c
include $(ERI4C)/module.mk
f90objs += $(ERI4C-iface)
dirs += $(ERI4C)
ttfs_libs += $(liberi4c.a)
ttfs_liberi4c += -L$(ERI4C) -leri4c
endif

# # # # # Electron Repulsion Integrals on 4 Center (ERI4C) # # # # # #

ifeq ($(WITH_DFTPU),1)
f90objs += modules/dft_plus_u_module.o
FPPOPTIONS += -DWITH_DFTPU
endif

# # # # # # # # # # # # GUILE SECTION # # # # # # # # # # # # # # # #

#
# Set WITH_GUILE = 1 to build a parallel Scheme interpreter "guile-qm"
# capable of running QM calculations with PG:
#
ifeq ($(WITH_GUILE),0)
main_objs = ./main.o
else

#
# Makefile for the Guile  interface, set $(GUILE) before include, This
# defines  $(libguile-comm.a) and  $(guile-qm.o)  targets among  ather
# things:
#
GUILE = ./guile
include $(GUILE)/module.mk

#
# $(guile-qm.o) contais  main() function.  FIXME: because  of a single
# fortran object  we do not  introduce another group/lib,  just append
# scheme.o here:
#
main_objs = $(guile-qm.o) modules/scheme.o
exe = guile-qm

f90objs += modules/scheme.o

dirs += $(GUILE)
ttfs_libs += $(libguile-comm.a)
extra_ttfs_libflags += -L$(GUILE) -lguile-comm
FPPOPTIONS += -DWITH_GUILE

#
# Two more  target to build.  These  scripts need to know  the path to
# "bindir" and/or "pkgdatadir". Inline relevant paths:
#
scripts = runqm #baslib/bin/qm-find-basis
build: $(scripts)

runqm: guile/runqm.scm
#baslib/bin/qm-find-basis: baslib/bin/qm-find-basis.scm

$(scripts):
	sed -e 's|@bindir[@]|$(BASEDIR)|g' -e 's|@pkgdatadir[@]|$(BASEDIR)|g' $(^) > $(@)
	chmod +x $(@)
endif

# # # # # # # # # # # # BGY3D SECTION # # # # # # # # # # # # # # # #
ifeq ($(WITH_BGY3D),1)
  FPPOPTIONS += -DWITH_BGY3D
  libttfs_callcomm.a += modules/bgy3d.o
endif
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#
# A large chunk of dependancies is historically included from there:
#
include Make.rules

##### PRIVATE LIBS #####

libttfs_callcomm.a: $(libttfs_callcomm.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_dontcallcomm.a: $(libttfs_dontcallcomm.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_opt.a: $(libttfs_opt.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_data.a: $(libttfs_data.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_pun.a: $(libttfs_pun.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_epe.a: $(libttfs_epe.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_efp.a: $(libttfs_efp.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_molmech.a: $(libttfs_molmech.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_shgi.a: $(libttfs_shgi.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_resp.a: $(libttfs_resp.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_gten.a: $(libttfs_gten.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_secder.a: $(libttfs_secder.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_comm.a: $(COMMOBJS)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_commc.a: $(COMMCOBJS)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libttfs_eis.a: $(EISOBJS)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

##### EXECUTABLE #####

# $(exe) points to an executable in local dir
# $(EXE) points to an executable in the installation dir
#
# Note that
#     $(FC) $(LINKFLAGS) -o $(EXE) ...
# will most probably overwrite inode (as in "cat NEWCODE > $(EXE)")
# but
#     mv $(exe) $(EXE)
# will change the directry entry to point to a new inode
# Maybe that will fix SIGBUS in running jobs while re-compiling.
#
define do_link
	$(FC) $(LINKFLAGS) -o $(exe) \
		-L. \
		$(main_objs) \
		-lttfs_callcomm \
		-lttfs_comm -lttfs_commc \
		-lttfs_dontcallcomm \
		-lttfs_callcomm \
		$(extra_ttfs_libflags) \
		-lttfs_dontcallcomm \
		-lttfs_data \
		-lttfs_pun \
		-lttfs_eis \
		$(ttfs_liberi4c) \
                $(LIBS)
endef

#INTERMEDIATE: $(exe)
#PRECIOUS: $(exe)

build: $(exe)
	@echo DONT FORGET TO make install $(exe)
install:
	@if [ ! -f $(exe) ]; then echo COMPILE FIRST by make $(exe); false; fi
	rm -f $(EXE)
	mv $(exe) $(EXE)

$(exe): $(main_objs) $(ttfs_libs)
	$(do_link)

# 'make dontlink' compiles but dosnt link
dontlink: $(main_objs) $(ttfs_libs)

# 'make link' should re-link the executable
link:
	$(do_link)

symlink-default      = mainscf_V$(MAJOR)$(MINOR)
symlink-default-arch = mainscf_V$(MAJOR)$(MINOR)-$(ARCH)
symlink-branch       = mainscf_V$(MAJOR)$(MINOR)$(EXTRA)
symlink-branch-arch  = mainscf_V$(MAJOR)$(MINOR)$(EXTRA)-$(ARCH)
symlink-default: symlink-default-arch
	if [ -L $(EXEDIR)/$(symlink-default) ]; then rm $(EXEDIR)/$(symlink-default); fi
	( cd $(EXEDIR); ln -s $(symlink-default-arch) $(symlink-default) )
symlink-default-arch: $(EXE)
	if [ -L $(EXEDIR)/$(symlink-default-arch) ]; then rm $(EXEDIR)/$(symlink-default-arch); fi
	( cd $(EXEDIR); ln -s $(exe) $(symlink-default-arch) )
symlink-branch: symlink-branch-arch
	if [ -L $(EXEDIR)/$(symlink-branch) ]; then rm $(EXEDIR)/$(symlink-branch); fi
	( cd $(EXEDIR); ln -s $(symlink-branch-arch) $(symlink-branch) )
symlink-branch-arch: $(EXE)
	if [ -L $(EXEDIR)/$(symlink-branch-arch) ]; then rm $(EXEDIR)/$(symlink-branch-arch); fi
	( cd $(EXEDIR); ln -s $(exe) $(symlink-branch-arch) )
symlinks: symlink-branch
# By "make symlinks" one creates  links to default executables for the
# current ?development? branch.   Create default !production! links by
# manual "make symlink-default"

#### INSTALLATION OF LIBRARIES
libraries = \
	epe_dir/AlumoSilicate-FF.epe \
	epe_dir/Alumina-FF.epe \
	epe_dir/MgO-FF.epe \
	epe_dir/bases.epe

libinstall: $(libraries)
	cp -r $(libraries) $(LIBDIR)

##### PRAEPROZESSOR BATCH FILES #####

##### SPECIAL COMMANDS #####


doku:	documentation/interfaces

documentation/interfaces:	$(DOCFILES)
	perl $(BINDIR)/collect_interfaces.perl $(DOCFILES)

reference:	documentation/reference.html

documentation/reference.html:	$(DOCFILES)
	perl $(BINDIR)/make_html_reference.perl -v $(VERS) $(DOCFILES)

realclean: clean libclean backupclean backupclean-cvs depclean

#
# Some of the include files (Petsc,  i am looking at you) may define a
# target named  clean already.  To avoid the  warnings we  augment its
# prerequisites here:
#
clean: clean-all

clean-all: libclean depclean $(clean)
	rm -f *.o */*.o */*/*.o
	rm -f $(cobjs) $(cobjs:.o=.d)
	rm -f *.mod */*.mod */*/*.mod
	rm -f *.F90 */*.F90 */*/*.F90
	rm -f *__genmod.f90 */*__genmod.f90 */*/*__genmod.f90

prepclean:
	rm -f $(f90objs:.o=.$(src_ext))

libclean:	#rm all library objects
	rm -f *.a

depclean:	#rm all .dep files that contain dependencies of files
	rm -f *.dep */*.dep */*/*.dep
	rm -f $(depfile) $(objlist)

backupclean:	#rm backup files automatically created by emacs, merges, makedependent
	rm -f *~ */*~ 
	rm -f findmerge.* *.oldmod */*.oldmod
	rm -f *.contrib* */*.contrib* *.merge* */*.merge*
	rm -f *.old */*.old
	rm -f *.orig */*.orig */*/*.orig
	rm -f *.bak */*.bak */*/*.bak

backupclean-cvs:
	rm -f .\#* */.\#* */*/.\#*

docclean:	# rm *.doc
	rm -f $(DOCFILES)

tags:
	ctags -l fortran *.f90 */*.f90 comm/mpi_dir/*.f90

TAGS:
	etags *.f90 */*.f90 comm/mpi_dir/*.f90

.PHONY: TAGS tags
# DO NOT DELETE

