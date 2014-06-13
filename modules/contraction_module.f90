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
module contraction_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
  !
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

  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind; ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: contract

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function n_indep_func(UA,L,IRR)
    use unique_atom_module, only: unique_atoms
    use options_module, only: options_spin_orbit
    implicit none
    integer(IK), intent(in) :: UA,L,IRR
    integer(IK)             :: n_indep_func ! result
    ! *** end of interface ***

    if(options_spin_orbit)then
       n_indep_func = unique_atoms(UA) &
            %symadapt_spor_partner(IRR,L) &
            %n_independent_fcts
    else
       n_indep_func = unique_atoms(UA) &
            %symadapt_partner(IRR,L) &
            %n_independent_fcts
    endif
  end function n_indep_func

  subroutine contract(i_irrep, matrix_uc, matrix_c)
    !
    ! written by. M.Mayer
    ! purpose: contract one irrep of a matrix necessary for SCF-part
    ! (overlap,kin,nuc)
    !
    use unique_atom_module
    implicit none
    integer(IK), intent(in) :: i_irrep ! dimension of irrep
    real(RK), intent(in), target :: matrix_uc(:, :) ! uncontracted matrix
    real(RK), intent(out), target :: matrix_c(:, :) ! contracted matrix
    ! *** end of interface ***

    integer(IK) :: i_unique_1,i_unique_2,i_l_1,i_l_2,i_ind_1,i_ind_2,i_n_1,i_n_2,&
         n_exponents1,n_exponents2,n_uncontracted1,n_uncontracted2,n_contracted1,n_contracted2,&
         i_uc_1,i_uc_2,i_c_1,i_c_2,&
         pointer_uc_2,pointer_c_2,pointer_uc_1,pointer_c_1
    real(RK),pointer,dimension(:,:) :: contraction1,contraction2,mat_p_uc,mat_p_c

    !------------ Executable code ------------------------------
    matrix_c=0.0_rk
    pointer_uc_1=1
    pointer_c_1=1
    do i_unique_1=1,n_unique_atoms
       do i_l_1=0,unique_atoms(i_unique_1)%lmax_ob
          contraction1=>unique_atoms(i_unique_1)%l_ob(i_l_1)%contractions
          n_contracted1=unique_atoms(i_unique_1)%l_ob(i_l_1)%n_contracted_fcts
          n_uncontracted1=unique_atoms(i_unique_1)%l_ob(i_l_1)%n_uncontracted_fcts
          n_exponents1=unique_atoms(i_unique_1)%l_ob(i_l_1)%n_exponents
             do i_ind_1=1,n_indep_func(i_unique_1,i_l_1,i_irrep)
!!$                unique_atoms(i_unique_1)%&
!!$                  symadapt_spor_partner(i_irrep,i_l_1)%n_independent_fcts
                pointer_uc_2=1
                pointer_c_2=1

                do i_unique_2=1,n_unique_atoms
                   do i_l_2=0,unique_atoms(i_unique_2)%lmax_ob
                         contraction2=>unique_atoms(i_unique_2)%l_ob(i_l_2)%contractions
                         n_contracted2=unique_atoms(i_unique_2)%l_ob(i_l_2)%n_contracted_fcts
                         n_uncontracted2=unique_atoms(i_unique_2)%l_ob(i_l_2)%n_uncontracted_fcts
                         n_exponents2=unique_atoms(i_unique_2)%l_ob(i_l_2)%n_exponents
                         do i_ind_2=1,n_indep_func(i_unique_2,i_l_2,i_irrep)
!!$                            unique_atoms(i_unique_2)%&
!!$                              symadapt_spor_partner(i_irrep,i_l_2)%n_independent_fcts

                            mat_p_uc=>matrix_uc(pointer_uc_2:pointer_uc_2+&
                                 n_exponents2-1,pointer_uc_1:pointer_uc_1+n_exponents1-1)
                            mat_p_c=>matrix_c(pointer_c_2:pointer_c_2+n_contracted2+&
                                 n_uncontracted2-1,pointer_c_1:&
                                 pointer_c_1+n_contracted1+n_uncontracted1-1)
                            ! now actual contraction starts
                            do i_uc_1=1,n_uncontracted1
                               do  i_uc_2=1,n_uncontracted2
                                  mat_p_c(i_uc_2,i_uc_1)=mat_p_uc(i_uc_2,i_uc_1)

                               end do
                               do i_c_2=1,n_contracted2
                                  do i_n_2=1,n_exponents2
                                     mat_p_c(i_c_2+n_uncontracted2,i_uc_1)=&
                                          mat_p_c(i_c_2+n_uncontracted2,i_uc_1)+mat_p_uc(i_n_2,i_uc_1)*&
                                          contraction2(i_n_2,i_c_2)

                                  end do
                               end do
                            end do
                            do i_c_1=1,n_contracted1
                               do i_n_1=1,n_exponents1
                                  do  i_uc_2=1,n_uncontracted2
                                     mat_p_c(i_uc_2,i_c_1+n_uncontracted1)=&
                                          mat_p_c(i_uc_2,i_c_1+n_uncontracted1)+mat_p_uc(i_uc_2,i_n_1)&
                                          *contraction1(i_n_1,i_c_1)
                                  end do
                                  do i_c_2=1,n_contracted2
                                     do i_n_2=1,n_exponents2
                                        mat_p_c(i_c_2+n_uncontracted2,i_c_1+n_uncontracted1)=&
                                             mat_p_c(i_c_2+n_uncontracted2,i_c_1+n_uncontracted1)+&
                                             mat_p_uc(i_n_2,i_n_1)*contraction2(i_n_2,i_c_2)*&
                                             contraction1(i_n_1,i_c_1)
                                     end do
                                  end do
                               end do
                            end do
                            pointer_uc_2=pointer_uc_2+n_exponents2
                            pointer_c_2=pointer_c_2+n_contracted2+n_uncontracted2
                         end do! loop over i_ind_2
                   enddo! loop over i_l_2
                end do! loop over i_unique_2
                pointer_uc_1=pointer_uc_1+n_exponents1
                pointer_c_1=pointer_c_1+n_contracted1+n_uncontracted1
             end do! loop over i_ind_1
       end do! loop over i_l_1
    end do! loop over i_unique_1
  end subroutine contract
  !*************************************************************

  !--------------- End of module ----------------------------------
end module contraction_module
