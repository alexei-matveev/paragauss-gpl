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
module symm_adapt_xpack
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

# include "def.h"
  use type_module, only:&
       & I2K=>i2_kind,&
       & IK=>i4_kind,&
       & RK=>r8_kind  ! type specification parameters
  !use error_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  interface cpack
     module procedure pack_symm_adapt_ul
     module procedure pack_symm_adapt_irrbas
     module procedure pack_symm_adapt_spinor
     module procedure pack_symm_adapt_ao
     module procedure pack_arr_symm_adapt_ul
  end interface

  interface cunpack
     module procedure unpack_symm_adapt_ul
     module procedure unpack_symm_adapt_irrbas
     module procedure unpack_symm_adapt_spinor
     module procedure unpack_symm_adapt_ao
     module procedure unpack_arr_symm_adapt_ul
  end interface

  !------------ public functions and subroutines ------------------

  public symm_adapt_send,&
       & symm_adapt_receive

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----



  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine symm_adapt_send()
    use comm_module
    use msgtag_module
    use spin_orbit_module, only: is_on, op_SpinOrbit
    implicit none
    ! *** end of interface ***

    if(.not.is_on(op_SpinOrbit)) return

    call comm_init_send(comm_all_other_hosts,msgtag_packed_message)
    call pack_everything()
    call comm_send()
  end subroutine symm_adapt_send

  subroutine symm_adapt_receive()
    use comm_module
    use msgtag_module
    use spin_orbit_module, only: is_on, op_SpinOrbit
    implicit none
    ! *** end of interface ***

    integer(IK) :: tag

    if(.not.is_on(op_SpinOrbit)) return

    call comm_save_recv(comm_master_host,msgtag_packed_message)

    tag = comm_msgtag()

    call unpack_everything()
  end subroutine symm_adapt_receive

  subroutine pack_everything()
    use error_module
    use xpack
    use symm_adapt_struct
    use spin_orbit_module, only: is_on, op_FitTrafo
    implicit none
    ! *** end of interface ***

    integer(IK) :: test=777

    DPRINT MyID,'saxp/pack_everything: entered'

    call cpack(LSymAdp )
    call cpack(LSymAdpL)
    if(is_on(op_FitTrafo))then
       call cpack(LSymAdpS)
    endif

    call pck(test)

    DPRINT MyID,'saxp/pack_everything: exited'
  end subroutine pack_everything

  subroutine unpack_everything()
    use error_module
    use xpack
    use symm_adapt_struct
    use spin_orbit_module, only: is_on, op_FitTrafo
    implicit none
    ! *** end of interface ***

    integer(IK) :: test

    DPRINT MyID,'saxp/unpack_everything: entered'

    call cunpack(LSymAdp )
    call cunpack(LSymAdpL)
    if(is_on(op_FitTrafo))then
       call cunpack(LSymAdpS)
    endif

    call upck(test)
    call error(test/=777,"'saxp/unpack_everything: test failed")

    DPRINT MyID,'saxp/unpack_everything: exited'
  end subroutine unpack_everything

  subroutine pack_arr_symm_adapt_ul(sa)
    use symm_adapt_struct
    implicit none
    type(symm_adapt_ul),intent(in) :: sa(:)
    ! *** end of interface ***

    integer(IK) :: i

    do i=1,size(sa)
       call cpack(sa(i))
    enddo
  end subroutine pack_arr_symm_adapt_ul

  subroutine unpack_arr_symm_adapt_ul(sa)
    use symm_adapt_struct
    implicit none
    type(symm_adapt_ul),intent(inout) :: sa(:)
    ! *** end of interface ***

    integer(IK) :: i

    do i=1,size(sa)
       call cunpack(sa(i))
    enddo
  end subroutine unpack_arr_symm_adapt_ul

  subroutine pack_symm_adapt_ul(sa)
    use xpack
    use symm_adapt_struct
    implicit none
    type(symm_adapt_ul),intent(in) :: sa
    ! *** end of interface ***

    integer(IK) :: l, irr, ea

    call pck(sa%n_ea)
    call pck(sa%lmax)
    call pck(sa%n_irr)

    call pck(sa%exists)

    do l = 0, sa%lmax
       do irr = 1, sa%n_irr
          if (sa%exists(l, irr))then
             do ea = 1, sa%n_ea
                call cpack (sa%aos(ea, l, irr))
             enddo
          endif
       enddo
    enddo
  end subroutine pack_symm_adapt_ul

  subroutine unpack_symm_adapt_ul(sa)
    use xpack
    use symm_adapt_struct
    implicit none
    type(symm_adapt_ul),intent(inout) :: sa
    ! *** end of interface ***

    integer(IK) :: l, irr, ea
    integer(IK) :: n_ea,lmax,n_irr

    !DPRINT MyID,'saxp/unpack_symm_adapt_ul: entered'

    call upck(n_ea)
    call upck(lmax)
    call upck(n_irr)

    call alloc(n_ea,lmax,n_irr,sa)

    call upck(sa%exists)

    do l = 0, lmax
       do irr = 1, n_irr
          if (sa%exists(l, irr))then
             do ea = 1, n_ea
                call cunpack (sa%aos(ea, l, irr))
             enddo
          endif
       enddo
    enddo
  end subroutine unpack_symm_adapt_ul
  !
  !--------------------------------------------

  subroutine pack_symm_adapt_irrbas(bs)
    use xpack
    use symm_adapt_struct
    implicit none
    type(symm_adapt_irrbas),intent(in) :: bs
    ! *** end of interface ***

    integer(IK) :: p,f

    call pck(bs%dim_irr)
    call pck(bs%n_indep)

    do p=1,bs%dim_irr
       do f=1,bs%n_indep
          call cpack(bs%bas(p,f))
       enddo
    enddo
  end subroutine pack_symm_adapt_irrbas

  subroutine unpack_symm_adapt_irrbas(bs)
    use xpack
    use symm_adapt_struct
    implicit none
    type(symm_adapt_irrbas),intent(inout) :: bs
    ! *** end of interface ***

    integer(IK) :: p,f
    integer(IK) :: dim_irr,n_indep

    call upck(dim_irr)
    call upck(n_indep)

    call alloc(dim_irr,n_indep,bs)

    do p=1,dim_irr
       do f=1,n_indep
          call cunpack(bs%bas(p,f))
       enddo
    enddo
  end subroutine unpack_symm_adapt_irrbas

  subroutine pack_symm_adapt_spinor(sp)
    use symm_adapt_struct
    implicit none
    type(symm_adapt_spinor),intent(in) :: sp
    ! *** end of interface ***

    integer(IK) :: c

    do c=1,2
       call cpack(sp%psi(c))
    enddo
  end subroutine pack_symm_adapt_spinor

  subroutine unpack_symm_adapt_spinor(sp)
    use symm_adapt_struct
    implicit none
    type(symm_adapt_spinor),intent(inout) :: sp
    ! *** end of interface ***

    integer(IK) :: c

    do c=1,2
       call cunpack(sp%psi(c))
    enddo
  end subroutine unpack_symm_adapt_spinor

  subroutine pack_symm_adapt_ao(ao)
    use xpack
    use symm_adapt_struct
    implicit none
    type(symm_adapt_ao),intent(in) :: ao
    ! *** end of interface ***

    call pck(ao%n)
    call pck(ao%m)
    call pck(ao%c)  !<<< both re and im
  end subroutine pack_symm_adapt_ao

  subroutine unpack_symm_adapt_ao(ao)
    use xpack
    use symm_adapt_struct
    implicit none
    type(symm_adapt_ao),intent(inout) :: ao
    ! *** end of interface ***

    integer(IK) :: n

    call upck(n)

    call alloc(n,ao)

    call upck(ao%m)
    call upck(ao%c)  !<<< both re and im
  end subroutine unpack_symm_adapt_ao

  !--------------- End of module ----------------------------------
end module symm_adapt_xpack
