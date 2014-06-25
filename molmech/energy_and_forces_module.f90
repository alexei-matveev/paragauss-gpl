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
module energy_and_forces_module
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use tasks_main_options_module
  use slab_module
  use species_module
  use potentials_module
  use qmmm_interface_module
  use molmech_msgtag_module
  use comm_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  real(kind=r8_kind), public :: E_total,E(100),E_coulomb,E_ew_d,E_ew_r
  real(kind=r8_kind), public :: E_solv_tot,E_solv_el,E_solv_dr,E_solv_cav
  real(kind=r8_kind), allocatable, public :: Grad(:,:)
  real(kind=r8_kind), allocatable, public :: Grad_s(:) ! Energy gradients in respect to strain
  real(kind=r8_kind), allocatable, public :: H(:,:)
  !------------ public functions and subroutines ------------------
  public init_energy_and_forces,write_energy_and_gradients,write_gxfile
  public shutdown_gradients, grads_to_qmmm, send_receive_e_g_h, send_receive_e_g_h_cascad
  public grad2frac, print_lat_param_and_grad
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine init_energy_and_forces()

    integer(kind=i4_kind) :: status

    E_total=zero
    E=zero
    E_coulomb=zero
    E_ew_d=zero
    E_ew_r=zero
    E_solv_tot=zero
    E_solv_el=zero
    E_solv_dr=zero
    E_solv_cav=zero
    
    if(calc_gradients) then
       if(calc_strain)  then 
          if(lattice_calc) then
             allocate(Grad_s(6),stat=status)
             if(status /= 0) call error_handler("MolMech: failed GRAD_S allocation(l)")
             Grad_s=zero
          else if(slab_calc) then
             allocate(Grad_s(3),stat=status)
             if(status /= 0) call error_handler("MolMech: failed GRAD_S allocation(l)")
             Grad_s=zero
          end if
       end if

       allocate(Grad(3,n_species),stat=status)
       if(status /= 0) call error_handler("MolMech: failed GRAD allocation")
       Grad=zero
    end if

  end subroutine init_energy_and_forces
  !****************************************************************

  !****************************************************************
  subroutine grad2frac()

    real(kind=r8_kind) :: bmat(3,3),bmat_s(2,2)
    integer(kind=i4_kind) :: i

    if(lattice_calc) then
       bmat(1,:)=vect%v1
       bmat(2,:)=vect%v2
       bmat(3,:)=vect%v3

       do i=1,n_species
          Grad(:,i)=matmul(bmat,Grad(:,i))
       end do
    else if(slab_calc) then
       bmat_s(1,:)=vect_s%v1
       bmat_s(2,:)=vect_s%v2

       do i=1,n_species
          Grad(1:2,i)=matmul(bmat_s,Grad(1:2,i))
       end do
    end if

  end subroutine grad2frac
  !****************************************************************

  !****************************************************************
  subroutine send_receive_e_g_h()
    
    integer(i4_kind) :: info,status,n_pr,i,ns
    real(r8_kind) :: E_t,E_c,E_r
    real(r8_kind), allocatable :: Gr(:,:),Gs(:),HH(:,:)

    ns=0 
    if(calc_strain) then
       if(lattice_calc) ns=6
       if(slab_calc) ns=3
    end if

    if(comm_i_am_master()) then
        n_pr=comm_get_n_processors()
        if(calc_gradients) then 
           allocate(Gr(3,n_species),stat=status)
           if(status /= 0) call error_handler("MolMech: failed GR allocation")
        end if
        if(calc_hessian) then
           allocate(HH(3*n_species+ns,3*n_species+ns),stat=status)
           if(status /= 0) call error_handler("MolMech: failed HH allocation")
        end if

        do i=2,n_pr
           call comm_save_recv(i,msgtag_mm_collect_gh)
           call communpack(E_t,info)
           E_total=E_total+E_t
           call communpack(E_c,info)
           E_coulomb=E_coulomb+E_c
           call communpack(E_r,info)
           E_ew_r=E_ew_r+E_r

           if(calc_gradients) then
              call communpack(Gr(1,1),3*n_species,1,info)
              Grad=Grad+Gr
              if(calc_strain) then
                 allocate(Gs(ns),stat=status)
                 if(status /= 0) call error_handler("MolMech: failed Gs allocation")
                 call communpack(Gs(1),ns,1,info)
                 Grad_s=Grad_s+Gs
              end if
           end if
           if(calc_hessian) then
              call communpack(HH(1,1),(3*n_species+ns)*(3*n_species+ns),1,info)
              H=H+HH
           end if
        end do

        if(calc_gradients) then 
           deallocate(Gr,stat=status)
           if(status /= 0) call error_handler("MolMech: failed GR deallocation")
           if(calc_strain) then
              deallocate(Gs,stat=status)
              if(status /= 0) call error_handler("MolMech: failed Gs dealloation")
           end if
        end if
        if(calc_hessian) then
           deallocate(HH,stat=status)
           if(status /= 0) call error_handler("MolMech: failed HH deallocation")
        end if
    else
       call comm_init_send(comm_master_host,msgtag_mm_collect_gh)
       call commpack(E_total,info)
       call commpack(E_coulomb,info)
       call commpack(E_ew_r,info)
       if(calc_gradients) then
          call commpack(Grad(1,1),3*n_species,1,info)
          if(calc_strain) then
             call commpack(Grad_s(1),ns,1,info)
          end if
       end if
       if(calc_hessian) then
          call commpack(H(1,1),(3*n_species+ns)*(3*n_species+ns),1,info)
       end if
       call comm_send()
    end if

  end subroutine send_receive_e_g_h
  !****************************************************************

  !****************************************************************
  subroutine send_receive_e_g_h_cascad()

    use casc_logic_module

    integer(i4_kind) :: np,status,ip,info,ns
    integer(i4_kind), allocatable :: vpn(:),rpn(:)
    integer(i4_kind) :: vpn1,vpn2,lower_node,upper_node
    real(r8_kind) :: E_t,E_c,E_r
    real(r8_kind), allocatable :: Gr(:,:),Gs(:),HH(:,:)

    ns=0
    if(calc_strain) then
       if(lattice_calc) ns=6
       if(slab_calc) ns=3
    end if

    ip=comm_myindex()
    np=comm_get_n_processors()
    allocate(vpn(np),rpn(2*np-1),stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech: allocation of VPN is failed")

    call comp_vpn(vpn)
    vpn1=vpn(ip); vpn2=vpn1
    call comp_rpn(rpn)

    if(2*(vpn1/2) == vpn1) then
       if(calc_gradients) then 
          allocate(Gr(3,n_species),stat=status)
          if(status /= 0) call error_handler("MolMech: failed GR allocation(Cascad)")
          if(calc_strain) then
             allocate(Gs(ns),stat=status)
             if(status /= 0) call error_handler("MolMech: failed Gs allocation(Cascad)")
          end if
       end if
       if(calc_hessian) then
          allocate(HH(3*n_species+ns,3*n_species+ns),stat=status)
          if(status /= 0) call error_handler("MolMech: failed HH allocation(Cascad)")
       end if
    end if

    vpn2_lab: do while(vpn2 > 1)
       if(mod(vpn2,2) == 0) then
          vpn2=vpn2/2
          lower_node=rpn(2*vpn2+1)

          call comm_save_recv(lower_node,msgtag_mm_collect_gh)
          call communpack(E_t,info)
          E_total=E_total+E_t
          call communpack(E_c,info)
          E_coulomb=E_coulomb+E_c
          call communpack(E_r,info)
          E_ew_r=E_ew_r+E_r
          if(calc_gradients) then
             call communpack(Gr(1,1),3*n_species,1,info)
             Grad=Grad+Gr
             if(calc_strain) then
                call communpack(Gs(1),ns,1,info)
                Grad_s=Grad_s+Gs
             end if
          end if
          if(calc_hessian) then
             call communpack(HH(1,1),(3*n_species+ns)*(3*n_species+ns),1,info)
             H=H+HH
          end if
       else
          upper_node=rpn((vpn2-1)/2)

          call comm_init_send(upper_node,msgtag_mm_collect_gh)
          call commpack(E_total,info)
          call commpack(E_coulomb,info)
          call commpack(E_ew_r,info)
          if(calc_gradients) then
             call commpack(Grad(1,1),3*n_species,1,info)
             if(calc_strain) then
                call commpack(Grad_s(1),ns,1,info)
             end if
          end if
          if(calc_hessian) then
             call commpack(H(1,1),(3*n_species+ns)*(3*n_species+ns),1,info)
          end if
          call comm_send()
          
          vpn2=0
       endif
    enddo vpn2_lab

    if(mod(vpn1,2)==0) then
       if(calc_gradients) then 
          deallocate(Gr,stat=status)
          if(status /= 0) call error_handler("MolMech: failed GR deallocation(Cascad)")
          if(calc_strain) then
             deallocate(Gs,stat=status)
             if(status /= 0) call error_handler("MolMech: failed Gs deallocation(Cascad)")
          end if
       end if
       if(calc_hessian) then
          deallocate(HH,stat=status)
          if(status /= 0) call error_handler("MolMech: failed HH deallocation(Cascad)")
       end if
    endif

    deallocate(vpn,rpn,stat=status)
    if ( status /= 0) call error_handler( &
         "MolMech: allocation of VPN is failed")

  end subroutine send_receive_e_g_h_cascad
  !****************************************************************

  !****************************************************************
  subroutine shutdown_gradients()

    integer(kind=i4_kind) :: status

    if(calc_gradients) then
       if(allocated(Grad)) then 
          deallocate(Grad,stat=status)
          if(status /= 0) call error_handler("MolMech: failed GRAD deallocation")
       end if
       if(calc_strain)  then 
          if(allocated(Grad_s)) then 
             deallocate(Grad_s,stat=status)
             if(status /= 0) call error_handler("MolMech: failed GRAD_S deallocation")
          end if
       end if
    end if

  end subroutine shutdown_gradients
  !****************************************************************

  !****************************************************************
  subroutine print_lat_param_and_grad()

    real(r8_kind) :: econv

    if(trim(energy_unit) == "KJ/MOL") then
       econv=one
    elseif(trim(energy_unit) == "KCAL/MOL") then 
       econv=j2c
    elseif(trim(energy_unit) == "EV") then
       econv=kjm2ev
    end if

    write(output_device,'(80("*"))')
    write(output_device,'(a65)') " Final unit cell parameters and gradients with respect to strain:"
    if(lattice_calc) then
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " a     = ",cel%a,    " angstrom","   Grad(1) = ",Grad_s(1)*econv
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " b     = ",cel%b,    " angstrom","   Grad(2) = ",Grad_s(2)*econv
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " c     = ",cel%c,    " angstrom","   Grad(3) = ",Grad_s(3)*econv
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " alpha = ",cel%alpha," degrees ","   Grad(4) = ",Grad_s(4)*econv
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " beta  = ",cel%beta, " degrees ","   Grad(5) = ",Grad_s(5)*econv
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " gamma = ",cel%gamma," degrees ","   Grad(6) = ",Grad_s(6)*econv
    else if(slab_calc) then
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " a     = ",cel_s%a,    " angstrom","   Grad(1) = ",Grad_s(1)*econv
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " b     = ",cel_s%b,    " angstrom","   Grad(2) = ",Grad_s(2)*econv
       write(output_device,'(a9,f13.7,a9,a13,f13.7)') &
            " alpha = ",cel_s%alpha," degrees ","   Grad(3) = ",Grad_s(3)*econv
    end if
    
  end subroutine print_lat_param_and_grad
  !****************************************************************

  !****************************************************************
  subroutine write_energy_and_gradients()

    integer(kind=i4_kind) :: i,ty
    real(r8_kind) :: econv
    character(len=26) :: format1
    character(len=9) :: en_u
    character(len=14) :: g_u
    character(len=1) :: c_s

    if(trim(energy_unit) == "KJ/MOL") then
       econv=one
       format1='(a12,1x,f24.8)'
       en_u=' '
       g_u='kJ/(mol*ang)  '
    elseif(trim(energy_unit) == "KCAL/MOL") then 
       econv=j2c
       format1='(a12,1x,f24.8,a8,f24.8,a9)'
       en_u=' kcal/mol'
       g_u='kcal/(mol*ang)'
    elseif(trim(energy_unit) == "EV") then
       econv=kjm2ev !0.01036433615d0 !kjm2ev
       format1='(a12,1x,f24.8,a8,f24.8,a9)'
       en_u=' eV      '
       g_u='eV/ang        '
    end if

    write(output_device,'(80("*"))')
    write(output_device,'(a18)') 'Energy of system:'
    write(output_device,'(80("-"))')
    write(output_device,'(a17)') 'Covalent energy:'

    write(output_device,'(80("."))')
    write(output_device,'(a17)') 'Two body energy:'
    do i=1,n_poten+n_user_pot
       if(i <= n_poten .and. E(i) /= zero) then
          if(poten_list%type(i)==2 .and. poten_list%bonded_type(i) ==0) &
               write(output_device,format1) &
               poten_list%name(i),E(i),' kJ/mol,',econv*E(i),en_u
       end if
    end do

    write(output_device,'(80("."))')
    write(output_device,'(a19)') 'Three body energy:'
    do i=1,n_poten+n_user_pot
       if(i <= n_poten .and. E(i) /= zero) then
          if(poten_list%type(i)==3) &
               write(output_device,format1) &
               poten_list%name(i),E(i),' kJ/mol,',econv*E(i),en_u
       end if
    end do

    write(output_device,'(80("."))')
    write(output_device,'(a18)') 'Four body energy:'
    do i=1,n_poten+n_user_pot
       if(i <= n_poten .and. E(i) /= zero) then
          if(poten_list%type(i)==4) &
               write(output_device,format1) &
               poten_list%name(i),E(i),' kJ/mol,',econv*E(i),en_u
       end if
    end do

    write(output_device,'(80("."))')
    write(output_device,'(a19)') 'Of diagonal terms:'
    do i=1,n_poten+n_user_pot
       if(i <= n_poten .and. E(i) /= zero) then
          if(poten_list%type(i)==5) &
               write(output_device,format1) &
               poten_list%name(i),E(i),' kJ/mol,',econv*E(i),en_u
       end if
    end do

    write(output_device,'(80("-"))')
    write(output_device,'(a22)') 'Van der Waals energy:'
    do i=1,n_poten+n_user_pot
       if(i <= n_poten .and. E(i) /= zero) then
          if(poten_list%type(i)==2 .and. poten_list%bonded_type(i) ==1) &
               write(output_device,format1) &
               poten_list%name(i),E(i),' kJ/mol,',econv*E(i),en_u
       else if(i > n_poten .and. E(i) /= zero) then
          write(output_device,format1) &
               'User defined',E(i),' kJ/mol,',econv*E(i),en_u
       end if
    end do

    write(output_device,'(80("-"))')
    write(output_device,'(a19)') ' Core-shell energy:'
    do i=1,n_poten+n_user_pot
       if(i <= n_poten .and. E(i) /= zero) then
          if(poten_list%type(i)==6 .and. poten_list%bonded_type(i) ==0) &
               write(output_device,format1) &
               poten_list%name(i),E(i),' kJ/mol,',econv*E(i),en_u
       end if
    end do

    if(coulomb) then
       write(output_device,'(80("-"))')
       write(output_device,'(a16)') 'Coulomb energy:'
       write(output_device,format1) &
            "",E_coulomb,' kJ/mol,',econv*E_coulomb,en_u
       if(lattice_calc .or. slab_calc) then
          write(output_device,'(a13)') 'Direct part:'
          write(output_device,format1) &
               "",E_ew_d,' kJ/mol,',econv*E_ew_d,en_u
          write(output_device,'(a17)') 'Reciprocal part:'
          write(output_device,format1) &
               "",E_ew_r,' kJ/mol,',econv*E_ew_r,en_u
       end if
    end if

    if(solvent) then
       write(output_device,'(80("-"))')
       write(output_device,'(a18)') 'Solvation energy:'
       write(output_device,format1) &
            "",E_solv_tot,' kJ/mol,',econv*E_solv_tot,en_u
       write(output_device,'(a20)') 'Electrostatic part:'
       write(output_device,format1) &
            "",E_solv_el,' kJ/mol,',econv*E_solv_el,en_u
       write(output_device,'(a19)') 'Cavitation energy:'
       write(output_device,format1) &
            "",E_solv_cav,' kJ/mol,',econv*E_solv_cav,en_u
       write(output_device,'(a17)') 'Disp-rep energy:'
       write(output_device,format1) &
            "",E_solv_dr,' kJ/mol,',econv*E_solv_dr,en_u
    end if

    write(output_device,'(80("="))')
    write(output_device,'(a14)') 'Total energy:'
    write(output_device,'(a12,1x,f24.8,a8,f24.8,a9)') &
         "",E_total,' kJ/mol,',econv*E_total,en_u
    write(output_device,'(80("*"))')

    if(calc_gradients) then
       write(output_device,'(a)') 'Gradients ('//trim(g_u)//'):'
       write(output_device,'(a76)') &
            '  Num  Init_num Name c/s        dE/dx              dE/dy              dE/dz '
       write(output_device,'(80("-"))')
       do i=1,n_species
          ty=atoms_cart(i)%type
          c_s="c"
          if(atoms(ty)%c_s == 1) c_s="s"
          write(output_device,'(i7,1x,i7,1x,a4,1x,a1,1x,1x,3f18.8)') &
               i,atoms_cart(i)%initial_number,atoms(ty)%name,c_s,Grad(:,i)*econv
       end do
       write(output_device,'(80("*"))')
    end if

  end subroutine write_energy_and_gradients
  !****************************************************************

  !****************************************************************
  subroutine write_gxfile()

    integer :: gx_file
    integer :: i,j,n_gx_atoms,mm_a
    real(r8_kind) :: atnum(max_gxatoms),coor(max_gxatoms,3),grad_gx(3)
    integer :: unum(max_gxatoms),num(max_gxatoms),znum(max_gxatoms,3),cnum(max_gxatoms,3)  


    atnum=zero; coor=zero; unum=0; num=0; znum=0; cnum=0

    if(with_optimizer) then
       call get_file_device(gx_file,'gxfile','inp')
    else
       call get_file_device(gx_file,'gx.mm','inp')
    end if

    n_gx_atoms=0; mm_a=0
    loop_1: do i=1,max_gxatoms
       read(gx_file,*,end=200) atnum(i),coor(i,1:3), &
            unum(i),num(i),znum(i,1:3),cnum(i,1:3)
       if(atnum(i) /= gx_dummy_atom .and. atnum(i) > zero) mm_a=mm_a+1
       n_gx_atoms=n_gx_atoms+1
       if(mm_a == n_species) exit loop_1
    end do loop_1
    read(gx_file,*,end=200) atnum(n_gx_atoms+1)
    

    rewind gx_file

    mm_a=0
    loop_2: do i=1,n_gx_atoms
       if(atnum(i) == gx_dummy_atom) then
          write(gx_file,100) atnum(i),coor(i,1:3), &
               unum(i),num(i),znum(i,1:3),cnum(i,1:3)
       else
          mm_a=mm_a+1
          loop_3: do j=1,n_species
             if(atoms_cart(j)%initial_number == mm_a) then
!!$                coor(i,1:3)=atoms_cart(j)%r*a2b
                write(gx_file,100) atnum(i),coor(i,1:3), &
                     unum(i),num(i),znum(i,1:3),cnum(i,1:3)
                exit loop_3
             end if
          end do loop_3
       end if
    end do loop_2
    
    write(gx_file,110) atnum(n_gx_atoms+1),coor(n_gx_atoms+1,1:3), &
         unum(n_gx_atoms+1),num(n_gx_atoms+1), &
         znum(n_gx_atoms+1,1:3),cnum(n_gx_atoms+1,1:3)
    
    write(gx_file,'(2F24.12,2x,3I5)') E_total/h2kJm,E_total/h2kJm, &
         n_gx_atoms,n_gx_atoms,n_gx_atoms

    loop_4: do i=1,n_species
       loop_5: do j=1,n_species
          if(atoms_cart(j)%initial_number == i) then
             grad_gx=Grad(:,j)/(h2kJm*a2b)
             if(i <= n_fixed) grad_gx =zero
             write(gx_file,'(I5,5x,3F17.12)') i,grad_gx
             exit loop_5
          end if
       end do loop_5
    end do loop_4

    call close_file_device(gx_file)

    return

100 format(f5.2,3(2x,f21.12),2i3,2x,3I4,2X,3I4)
110 format(f5.0,3(2X,f21.12),2i3,2x,3I4,2X,3I4)

200 call error_handler("MolMech (write_gxfile): Error during reading GX file")

  end subroutine write_gxfile
  !****************************************************************

  !****************************************************************
  subroutine grads_to_qmmm()
    !for IMOMM

    integer(i4_kind) :: i

    if(qmmm == 2) then
       energy_mm2=E_total/h2kJm
       do i=1,n_species
          grad_mm2(i)%x=Grad(1,i)/(h2kJm*a2b)
          grad_mm2(i)%y=Grad(2,i)/(h2kJm*a2b)
          grad_mm2(i)%z=Grad(3,i)/(h2kJm*a2b)
       end do
    elseif(qmmm == 3) then
       energy_mm3=E_total/h2kJm
       do i=1,n_species
          grad_mm3(i)%x=Grad(1,i)/(h2kJm*a2b)
          grad_mm3(i)%y=Grad(2,i)/(h2kJm*a2b)
          grad_mm3(i)%z=Grad(3,i)/(h2kJm*a2b)
       end do
    end if

  end subroutine grads_to_qmmm
  !****************************************************************
end module energy_and_forces_module
