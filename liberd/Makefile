
all: liberd.a

liberd.a = \
        erd.o \
        erd__1111_csgto.o \
        erd__1111_def_blocks.o \
        erd__2d_atom_coefficients.o \
        erd__2d_atom_pq_integrals.o \
        erd__2d_coefficients.o \
        erd__2d_pq_integrals.o \
        erd__atom_int2d_to_e000.o \
        erd__atom_int2d_to_e0f0.o \
        erd__cartesian_norms.o \
        erd__csgto.o \
        erd__ctr_1st_half.o \
        erd__ctr_2nd_half_new.o \
        erd__ctr_2nd_half_update.o \
        erd__ctr_4index_block.o \
        erd__ctr_4index_reorder.o \
        erd__ctr_rs_expand.o \
        erd__ctr_tu_expand.o \
        erd__dsqmin_line_segments.o \
        erd__e0f0_def_blocks.o \
        erd__e0f0_pcgto_block.o \
        erd__gener_eri_batch.o \
        erd__hrr_matrix.o \
        erd__hrr_step.o \
        erd__hrr_transform.o \
        erd__int2d_to_e000.o \
        erd__int2d_to_e0f0.o \
        erd__map_ijkl_to_ikjl.o \
        erd__memory_1111_csgto.o \
        erd__memory_csgto.o \
        erd__memory_eri_batch.o \
        erd__move_ry.o \
        erd__normalize_cartesian.o \
        erd__pppp_pcgto_block.o \
        erd__prepare_ctr.o \
        erd__rys_1_roots_weights.o \
        erd__rys_2_roots_weights.o \
        erd__rys_3_roots_weights.o \
        erd__rys_4_roots_weights.o \
        erd__rys_5_roots_weights.o \
        erd__rys_roots_weights.o \
        erd__rys_x_roots_weights.o \
        erd__set_abcd.o \
        erd__set_ij_kl_pairs.o \
        erd__spherical_transform.o \
        erd__sppp_pcgto_block.o \
        erd__sspp_pcgto_block.o \
        erd__sssp_pcgto_block.o \
        erd__ssss_pcgto_block.o \
        erd__transpose_batch.o \
        erd__xyz_to_ry_abcd.o \
        erd__xyz_to_ry_matrix.o \
#       erd__2d_pcd_integrals.o \
#       erd__2d_abcd_integrals.o \
#       erd__2d_derv_abcd_integrals.o \
#       erd__2d_derv_pcd_integrals.o \
#       erd__gener_eri_derv_batch.o \
#       erd__memory_derv_csgto.o \
#       erd__memory_eri_derv_batch.o \
#       erd__set_derv_abcd.o \
#       erd__set_derv_sequence.o \
#       erd__derv_csgto.o \
#       erd__derv_def_blocks.o \
#       erd__derv_int2d_to_0000.o \
#       erd__derv_int2d_to_a000.o \
#       erd__derv_int2d_to_ab00.o \
#       erd__derv_int2d_to_abc0.o \
#       erd__derv_int2d_to_abcd.o \
#       erd__derv_int2d_to_e000.o \
#       erd__derv_int2d_to_e0c0.o \
#       erd__derv_int2d_to_e0cd.o \
#       erd__derv_pcgto_block.o \
#       erd__print_batch.o \

# defines $(FC) and other:
include ../machine.inc

liberd.a : $(liberd.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

FLAGS = $(FBASEFLAGS) $(F90FLAGS) # -fbounds-check
%.o : %.f
	$(FC) $(FLAGS) -c $(<F)

%.o : %.f90
	$(FC) $(FLAGS) -c $(<F)

clean :
	rm -f $(liberd.a) liberd.a

