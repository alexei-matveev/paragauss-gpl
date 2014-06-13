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

module coulomb_module
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use tasks_main_options_module
  use slab_module, only: slab_calc
  use species_module
  use external_field_module
  use energy_and_forces_module
  use n_body_lists_module
  use comm_module, only: comm_get_n_processors

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  character(len=13), public :: Coulomb_type

  type, public :: dipole_pair
     integer(kind=i4_kind) :: at_num(2)
     real(kind=r8_kind) :: D
  end type dipole_pair
  type(dipole_pair), public, allocatable :: dip_pair(:)
  integer(kind=i4_kind), public :: N_dm_pair
  integer(kind=i4_kind), public, allocatable :: N_interact(:)
  !------------ public functions and subroutines ------------------
  public read_coulomb_options, write_coulomb_to_output, read_dipoles, dip_atom_pairs, & 
       direct_coulomb_E_and_F, dipol_dipol_E_and_F, write_dipoles_to_output, &
       shutdown_electrostatic
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  type dipole
     integer(kind=i4_kind) :: atom_type(2)
     real(kind=r8_kind) :: D
  end type dipole
  type(dipole), allocatable :: dipole_bond(:)
  integer(kind=i4_kind) :: n_dip

  character(len=13), parameter :: df_coulomb_type="DIRECT" !or "DIPOLE_DIPOLE", "EWALD" 
  ! or "ewald"
  namelist /electrostatic/ coulomb_type

  character(len=len_name) :: atom_name(2)
  real(kind=r8_kind) :: dip_mom
  character(len=len_name) :: df_atom_name(2)=(/'      ','      '/)
  real(kind=r8_kind) :: df_dip_mom=zero
  namelist /bond_dipole/ atom_name, dip_mom

#ifdef FPP_GFORTRAN_BUGS
  ! Gfortran 4.3 doesnt like this typedef inside a subroutine:
  type inp_data
     type (dipole) :: data
     type (inp_data), pointer :: next_data
  end type inp_data
#endif
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  function read_coulomb_options()

    logical :: read_coulomb_options
    integer(i4_kind) :: i

    coulomb_type=df_coulomb_type

    call go_to_first_input_line()
    read_coulomb_options=find_namelist("&ELECTROSTATIC",i)
    if(.not.read_coulomb_options) goto 100
    read(input_device,nml=electrostatic, end=100, err=200)
    call upcase(coulomb_type)

    if(lattice_calc .and. trim(coulomb_type) /= "EWALD") coulomb_type="EWALD"
    if(slab_calc .and. trim(coulomb_type) /= "EWALD") coulomb_type="EWALD"

    read_coulomb_options=.true.
    return

100 read_coulomb_options=.false.
    if(lattice_calc .and. trim(coulomb_type) /= "EWALD") coulomb_type="EWALD"
    if(slab_calc .and. trim(coulomb_type) /= "EWALD") coulomb_type="EWALD"
    return

200 call input_nm_error(0,"ELECTROSTATIC")

  end function read_coulomb_options
  !****************************************************************

  !****************************************************************
  subroutine write_coulomb_to_output()

    write(output_device,'(a43,a13)') " Description of Electrostatic interaction: ",coulomb_type
    if(comm_get_n_processors() > 1 .and. trim(coulomb_type) == "EWALD") then
       write(output_device,'(a37,i3,a11)') " Reciprocal part runs in parallel on ", &
            comm_get_n_processors()," processors"
    end if
    write(output_device,'(80("-"),/)')

  end subroutine write_coulomb_to_output
  !****************************************************************

  !****************************************************************
  function read_dipoles()

    logical :: read_dipoles
    
    integer(kind=i4_kind) :: count,type_a(2),i,status,ii
    character(len=len_name) :: name
#ifndef FPP_GFORTRAN_BUGS
    ! Gfortran 4.3 doesnt like this typedef inside a subroutine:
    type inp_data
       type (dipole) :: data
       type (inp_data), pointer :: next_data
    end type inp_data
#endif
    type (inp_data), target :: first_data
    type (inp_data), pointer :: current_data, tmp_data

    current_data=>first_data
    nullify(current_data%next_data)

    count=0_i4_kind
    call go_to_first_input_line()
    do
       atom_name=df_atom_name
       dip_mom=df_dip_mom

       read_dipoles=find_namelist("&BOND_DIPOLE",ii)
       if(.not.read_dipoles) goto 100
       read(input_device,nml=bond_dipole, end=100, err=200)
       if(dip_mom == zero) cycle
       if(atom_name(1)=="      " .or. atom_name(2)=="      ") goto 200

       call name_without_cs(atom_name(1),name)
       type_a(1)=name2type(name,"C")
       if(type_a(1)==0) cycle
       call name_without_cs(atom_name(2),name)
       type_a(2)=name2type(name,"C")
       if(type_a(2)==0) cycle

       allocate(tmp_data)
       tmp_data%data=dipole(type_a,dip_mom)
       nullify(tmp_data%next_data)
       current_data%next_data =>tmp_data
       current_data => tmp_data

       count=count+1
    end do
100 n_dip=count

    if(n_dip==0) then
       read_dipoles=.false.
       return
    else
       allocate(dipole_bond(n_dip),stat=status)
       if(status /= 0) call error_handler("MolMech: failed DIPOLE_BOND allocation")

       current_data => first_data

       do i=1,n_dip
          current_data => current_data%next_data
          dipole_bond(i)=current_data%data
       end do
       read_dipoles=.true.
       return
    end if

200 call input_nm_error(count+1,"BOND_DIPOLE")
  end function read_dipoles
  !****************************************************************

  !****************************************************************
  subroutine write_dipoles_to_output()

    integer(kind=i4_kind) :: i,j,k
    character(len=len_name) :: nm1,nm2

    write(output_device,'(a24)') "Dipole moments of bonds:"

    do i=1,n_dip
       j=dipole_bond(i)%atom_type(1)
       k=dipole_bond(i)%atom_type(2)
       nm1=atoms(j)%name
       nm2=atoms(k)%name
       write(output_device,'(a6,1x,a6,f15.7,a6)') nm1, nm2, dipole_bond(i)%D, " Debay"
    end do
    write(output_device,'(80("-"),/)')

  end subroutine write_dipoles_to_output
  !****************************************************************

  !****************************************************************
  subroutine direct_coulomb_E_and_F

    integer(kind=i4_kind) :: i,j,k,length,nn,n,k1,k2,l,m,m1,l1,l2
    integer(kind=i4_kind) :: a_type1,a_type2,ns
    real(kind=r8_kind) :: q1,q2,r(3),r1(3),dr,E_buf
    real(kind=r8_kind) :: f1(3),f2(3),dE_dr
    real(kind=r8_kind) :: d2E_dr2,ff(6,6),fstore

    ns=3*n_species

    do i=1,n_species
       a_type1=atoms_cart(i)%type
       q1=atoms(a_type1)%charge
       length=size(bonded_spcs(i)%list,1)

       nn=1
       j_n: do j=1,N_total(i)+n_ext_pc
          if(j <= N_total(i)) then !regular centers
             k=i+j
             if(k > n_species) k=k-n_species
             a_type2=atoms_cart(k)%type
             q2=atoms(a_type2)%charge
             r1=atoms_cart(k)%r
          else      !external centers
             k=j-N_total(i)+n_species
             a_type2=pc_cart(k-n_species)%type
             q2=p_charges(a_type2-n_species_types)%charge
             r1=pc_cart(k-n_species)%r
          end if

          n_l: do n=nn,length
             if(k == bonded_spcs(i)%list(n,1)) then
                nn=n+1
                cycle j_n
             end if
          end do n_l

          r=r1-atoms_cart(i)%r
          dr=sqrt(dot_product(r,r))

          E_buf=coulomb_factor*q1*q2/(dr*eps)
          E_coulomb=E_coulomb+E_buf
          if(calc_gradients) dE_dr=-coulomb_factor*q1*q2/(dr**2*eps)
          if(calc_hessian) d2E_dr2=two*coulomb_factor*q1*q2/(dr**3*eps)

          E_total=E_total+E_buf

          if(calc_gradients) then
             f1=r/dr
             f2=-f1
             if(i <= n_species) Grad(:,i)=Grad(:,i)+dE_dr*f2
             if(k <= n_species) Grad(:,k)=Grad(:,k)+dE_dr*f1
          end if

          if(calc_hessian) then
             ff(1,1)=one/dr-r(1)*r(1)/dr**3
             ff(1,2)=-r(1)*r(2)/dr**3
             ff(1,3)=-r(1)*r(3)/dr**3
             ff(1,4)=-ff(1,1)
             ff(1,5)=-ff(1,2)
             ff(1,6)=-ff(1,3)
             ff(2,1)=-r(2)*r(1)/dr**3
             ff(2,2)=one/dr-r(2)*r(2)/dr**3
             ff(2,3)=-r(2)*r(3)/dr**3
             ff(2,4)=-ff(2,1)
             ff(2,5)=-ff(2,2)
             ff(2,6)=-ff(2,3)
             ff(3,1)=-r(3)*r(1)/dr**3
             ff(3,2)=-r(3)*r(2)/dr**3
             ff(3,3)=one/dr-r(3)*r(3)/dr**3
             ff(3,4)=-ff(3,1)
             ff(3,5)=-ff(3,2)
             ff(3,6)=-ff(3,3)
             ff(4,1:3)=ff(1:3,4)
             ff(5,1:3)=ff(1:3,5)
             ff(6,1:3)=ff(1:3,6)
             ff(4,4:6)=ff(1,1:3)
             ff(5,4:6)=ff(2,1:3)
             ff(6,4:6)=ff(3,1:3)
             
             k1=3*(i-1)
             k2=3*(k-1)
             do l=1,3
                do m=1,6
                   if(m <= 3) then
                      fstore=f1(m)
                      m1=k1+m
                   elseif(m > 3) then 
                      m1=k2+(m-3)
                      fstore=f2(m-3)
                   end if
                   l1=k1+l; l2=k2+l
                   if(l1<=ns .and. m1<=ns) H(l1,m1)=H(l1,m1)-d2E_dr2*fstore*f2(l)-dE_dr*ff(l+3,m)
                   if(l2<=ns .and. m1<=ns) H(l2,m1)=H(l2,m1)-d2E_dr2*fstore*f1(l)-dE_dr*ff(l,m)
                end do
             end do
          end if
 
       end do j_n
    end do

  end subroutine direct_coulomb_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine dip_atom_pairs()

    integer(kind=i4_kind) :: i,j,k,status
    integer(kind=i4_kind) :: ty(2),num(2),count
    real(kind=r8_kind) :: r(3),dr,r_cov

    type inp_data
       type (dipole_pair) :: data
       type (inp_data), pointer :: next_data
    end type inp_data
    type (inp_data), target :: first_data
    type (inp_data), pointer :: current_data, tmp_data

    !definition of all possible atom pairs forming chemical bond and
    !having specified dipole moment
    current_data=>first_data
    nullify(current_data%next_data)

    count=0
    l_i: do i=1,n_species-1
       ty(1)=atoms_cart(i)%type

       l_j: do j=i+1,n_species
          ty(2)=atoms_cart(j)%type
          r=atoms_cart(j)%r-atoms_cart(i)%r
          dr=sqrt(dot_product(r,r))
          r_cov=atoms(ty(1))%r_coval+atoms(ty(2))%r_coval
          if(dr > rcf*r_cov) cycle l_j

          l_k: do k=1,n_dip
             if(ty(1) == dipole_bond(k)%atom_type(1) .and. &
                  ty(2) == dipole_bond(k)%atom_type(2)) then
                num(1)=i; num(2)=j
             else if(ty(1) == dipole_bond(k)%atom_type(2) .and. &
                  ty(2) == dipole_bond(k)%atom_type(1)) then
                num(2)=i; num(1)=j
             else
                cycle l_k
             end if

             allocate(tmp_data)
             tmp_data%data=dipole_pair(num,dipole_bond(k)%D)
             nullify(tmp_data%next_data)
             current_data%next_data =>tmp_data
             current_data => tmp_data
             
             count=count+1
             exit l_k
          end do l_k
       end do l_j
    end do l_i
    N_dm_pair=count

    allocate(dip_pair(N_dm_pair),stat=status)
    if(status /= 0) call error_handler("MolMech: failed DIP_PAIR allocation")

    current_data => first_data
    do i=1,N_dm_pair
       current_data => current_data%next_data
       dip_pair(i)=current_data%data
    end do

    deallocate(dipole_bond,stat=status)
    if(status /= 0) call error_handler("MolMech: failed DIPOLE_BOND deallocation")

    allocate(N_interact(N_dm_pair),stat=status)
    if(status /= 0) call error_handler("MolMech: failed N_interact allocation")

    if(mod(N_dm_pair,2)==0) then
       do i=1,N_dm_pair
          if(i<=N_dm_pair/2 .and. mod(i,2)==1) N_interact(i)=N_dm_pair/2
          if(i<=N_dm_pair/2 .and. mod(i,2)==0) N_interact(i)=N_dm_pair/2-1
          if(i >N_dm_pair/2 .and. mod(i-N_dm_pair/2,2)==1) N_interact(i)=N_dm_pair/2-1
          if(i >N_dm_pair/2 .and. mod(i-N_dm_pair/2,2)==0) N_interact(i)=N_dm_pair/2
       end do
    else if(mod(N_dm_pair-1,2)==0) then
       N_interact=(N_dm_pair-1)/2
    end if

  end subroutine dip_atom_pairs
  !****************************************************************

  !****************************************************************
  subroutine dipol_dipol_E_and_F()

    integer(kind=i4_kind) :: i,j,jj,k1,k2,k3,k4,l,m,m1
    integer(kind=i4_kind) :: ia1,ia2,ia3,ia4
    real(kind=r8_kind) :: r1(3),r2(3),r3(3),r4(3) 
    real(kind=r8_kind) :: Ra(3),Rb(3),ra2(3),rb2(3),Rab(3),dra,drb,drab
    real(kind=r8_kind) :: xci,cos_xci,alpha,cos_alpha,beta,cos_beta
    real(kind=r8_kind) :: E_buf0,E_buf,da,db
    real(kind=r8_kind) :: dE_dr,dE_dx,dE_da,dE_db
    real(kind=r8_kind) :: fr(4,3),fx(4,3),fa(4,3),fb(4,3)
    real(r8_kind) :: d2E_dr2,d2E_drdx,d2E_drda,d2E_drdb
    real(r8_kind) :: d2E_dxdr,d2E_dx2,d2E_dxda,d2E_dxdb
    real(r8_kind) :: d2E_dadr,d2E_dadx,d2E_da2,d2E_dadb
    real(r8_kind) :: d2E_dbdr,d2E_dbdx,d2E_dbda,d2E_db2
    real(r8_kind) :: ffr(12,12),ffx(12,12),ffa(12,12),ffb(12,12),fst_r,fst_x,fst_a,fst_b
    real(r8_kind) :: d_sin(12),drr(12),ddra(12),ddrb(12),rr,rra(3),rrb(3)

    E_buf0=dipole_factor*c2J/dip_eps
!!$print*,N_dm_pair,'N_dm_pair'
    l_i: do i=1,N_dm_pair
       ia1=dip_pair(i)%at_num(1)
       ia2=dip_pair(i)%at_num(2)
       da=dip_pair(i)%d
       r2=atoms_cart(ia2)%r
       r1=atoms_cart(ia1)%r
       Ra=r2-r1
       dra=sqrt(dot_product(Ra,Ra))
       ra2=(atoms_cart(ia2)%r+atoms_cart(ia1)%r)/two

       l_j: do j=1,N_interact(i)
          jj=i+j
          if(jj > N_dm_pair) jj=jj-N_dm_pair
          ia3=dip_pair(jj)%at_num(1)
          ia4=dip_pair(jj)%at_num(2)
          db=dip_pair(jj)%d
          if(ia1==ia3 .or. ia1==ia4) cycle l_j
          if(ia2==ia3 .or. ia2==ia4) cycle l_j

          r3=atoms_cart(ia3)%r
          r4=atoms_cart(ia4)%r
          Rb=r4-r3
          drb=sqrt(dot_product(Rb,Rb))
          rb2=(atoms_cart(ia4)%r+atoms_cart(ia3)%r)/two
          Rab=rb2-ra2
          drab=sqrt(dot_product(Rab,Rab))

          cos_xci=dot_product(Ra,Rb)/(dra*drb)
          xci=acos(cos_xci)
          cos_alpha=dot_product(Ra,Rab)/(dra*drab)
          alpha=acos(cos_alpha)
          cos_beta=dot_product(Rb,Rab)/(drb*drab)
          beta=acos(cos_beta)

          E_buf=E_buf0*(da*db/drab**3)*(cos_xci-three*cos_alpha*cos_beta)
          E_coulomb=E_coulomb+E_buf
          E_total=E_total+E_buf

          if(calc_gradients) then
             dE_dr=-E_buf0*three*(da*db/drab**4)*(cos_xci-three*cos_alpha*cos_beta)
             dE_dx=-E_buf0*(da*db/drab**3)*sin(xci)
             dE_da=E_buf0*(da*db/drab**3)*three*sin(alpha)*cos_beta
             dE_db=E_buf0*(da*db/drab**3)*three*cos_alpha*sin(beta)

             fr(1,:)=-Rab/(two*drab)
             fr(2,:)=-Rab/(two*drab)
             fr(3,:)=Rab/(two*drab)
             fr(4,:)=Rab/(two*drab)

             fx(1,:)= (Rb/(dra*drb)-Ra*dot_product(Ra,Rb)/(dra**3*drb))/sin(xci)
             fx(2,:)=(-Rb/(dra*drb)+Ra*dot_product(Ra,Rb)/(dra**3*drb))/sin(xci)
             fx(3,:)= (Ra/(dra*drb)-Rb*dot_product(Ra,Rb)/(dra*drb**3))/sin(xci)
             fx(4,:)=(-Ra/(dra*drb)+Rb*dot_product(Ra,Rb)/(dra*drb**3))/sin(xci)

             fa(1,:)=-(-(r4+r3-two*r1)/(two*dra*drab)+Ra*dot_product(Ra,Rab)/(dra**3*drab)+ &
                  Rab*dot_product(Ra,Rab)/(two*dra*drab**3))/sin(alpha)
             fa(2,:)= -((r4+r3-two*r2)/(two*dra*drab)-Ra*dot_product(Ra,Rab)/(dra**3*drab)+ &
                  Rab*dot_product(Ra,Rab)/(two*dra*drab**3))/sin(alpha)
             fa(3,:)= -(Ra/(two*dra*drab)-Rab*dot_product(Ra,Rab)/(two*dra*drab**3))/sin(alpha)
             fa(4,:)= -(Ra/(two*dra*drab)-Rab*dot_product(Ra,Rab)/(two*dra*drab**3))/sin(alpha)

             fb(1,:)=-(-Rb/(two*drb*drab)+Rab*dot_product(Rb,Rab)/(two*drb*drab**3))/sin(beta)
             fb(2,:)=-(-Rb/(two*drb*drab)+Rab*dot_product(Rb,Rab)/(two*drb*drab**3))/sin(beta)
             fb(3,:)=-(-(two*r3-r2-r1)/(two*drb*drab)+Rb*dot_product(Rb,Rab)/(drb**3*drab)- &
                  Rab*dot_product(Rb,Rab)/(two*drb*drab**3))/sin(beta)
             fb(4,:)= -((two*r4-r2-r1)/(two*drb*drab)-Rb*dot_product(Rb,Rab)/(drb**3*drab)- &
                  Rab*dot_product(Rb,Rab)/(two*drb*drab**3))/sin(beta)

             Grad(:,ia1)=Grad(:,ia1)+dE_dr*fr(1,:)+dE_dx*fx(1,:)+dE_da*fa(1,:)+dE_db*fb(1,:)
             Grad(:,ia2)=Grad(:,ia2)+dE_dr*fr(2,:)+dE_dx*fx(2,:)+dE_da*fa(2,:)+dE_db*fb(2,:)
             Grad(:,ia3)=Grad(:,ia3)+dE_dr*fr(3,:)+dE_dx*fx(3,:)+dE_da*fa(3,:)+dE_db*fb(3,:)
             Grad(:,ia4)=Grad(:,ia4)+dE_dr*fr(4,:)+dE_dx*fx(4,:)+dE_da*fa(4,:)+dE_db*fb(4,:)
          end if
          
          if(calc_hessian) then
             d2E_dr2=E_buf0*twelve*(da*db/drab**5)*(cos_xci-three*cos_alpha*cos_beta)
             d2E_drdx=E_buf0*three*(da*db/drab**4)*sin(xci)
             d2E_drda=-E_buf0*three*(da*db/drab**4)*three*sin(alpha)*cos_beta
             d2E_drdb=-E_buf0*three*(da*db/drab**4)*three*cos_alpha*sin(beta)
             !...............................
             d2E_dxdr=d2E_drdx
             d2E_dx2=-E_buf0*(da*db/drab**3)*cos_xci
             d2E_dxda=zero
             d2E_dxdb=zero
             !...............................
             d2E_dadr=d2E_drda
             d2E_dadx=d2E_dxda
             d2E_da2=E_buf0*(da*db/drab**3)*three*cos_alpha*cos_beta
             d2E_dadb=-E_buf0*(da*db/drab**3)*three*sin(alpha)*sin(beta)
             !...............................
             d2E_dbdr=d2E_drdb
             d2E_dbdx=d2E_dxdb
             d2E_dbda=d2E_dadb
             d2E_db2=E_buf0*(da*db/drab**3)*three*cos_alpha*cos_beta
             !...............................
             !...............................
             ffr(1,1)=one/(four*drab)+Rab(1)*fr(1,1)/(two*drab**2)
             ffr(1,2)=Rab(1)*fr(1,2)/(two*drab**2)
             ffr(1,3)=Rab(1)*fr(1,3)/(two*drab**2)
             ffr(1,4)=one/(four*drab)+Rab(1)*fr(2,1)/(two*drab**2)
             ffr(1,5)=Rab(1)*fr(2,2)/(two*drab**2)
             ffr(1,6)=Rab(1)*fr(2,3)/(two*drab**2)
             ffr(1,7)=-one/(four*drab)+Rab(1)*fr(3,1)/(two*drab**2)
             ffr(1,8)=Rab(1)*fr(3,2)/(two*drab**2)
             ffr(1,9)=Rab(1)*fr(3,3)/(two*drab**2)
             ffr(1,10)=-one/(four*drab)+Rab(1)*fr(4,1)/(two*drab**2)
             ffr(1,11)=Rab(1)*fr(4,2)/(two*drab**2)
             ffr(1,12)=Rab(1)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(2,1)=Rab(2)*fr(1,1)/(two*drab**2)
             ffr(2,2)=one/(four*drab)+Rab(2)*fr(1,2)/(two*drab**2)
             ffr(2,3)=Rab(2)*fr(1,3)/(two*drab**2)
             ffr(2,4)=Rab(2)*fr(2,1)/(two*drab**2)
             ffr(2,5)=one/(four*drab)+Rab(2)*fr(2,2)/(two*drab**2)
             ffr(2,6)=Rab(2)*fr(2,3)/(two*drab**2)
             ffr(2,7)=Rab(2)*fr(3,1)/(two*drab**2)
             ffr(2,8)=-one/(four*drab)+Rab(2)*fr(3,2)/(two*drab**2)
             ffr(2,9)=Rab(2)*fr(3,3)/(two*drab**2)
             ffr(2,10)=Rab(2)*fr(4,1)/(two*drab**2)
             ffr(2,11)=-one/(four*drab)+Rab(2)*fr(4,2)/(two*drab**2)
             ffr(2,12)=Rab(2)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(3,1)=Rab(3)*fr(1,1)/(two*drab**2)
             ffr(3,2)=Rab(3)*fr(1,2)/(two*drab**2)
             ffr(3,3)=one/(four*drab)+Rab(3)*fr(1,3)/(two*drab**2)
             ffr(3,4)=Rab(3)*fr(2,1)/(two*drab**2)
             ffr(3,5)=Rab(3)*fr(2,2)/(two*drab**2)
             ffr(3,6)=one/(four*drab)+Rab(3)*fr(2,3)/(two*drab**2)
             ffr(3,7)=Rab(3)*fr(3,1)/(two*drab**2)
             ffr(3,8)=Rab(3)*fr(3,2)/(two*drab**2)
             ffr(3,9)=-one/(four*drab)+Rab(3)*fr(3,3)/(two*drab**2)
             ffr(3,10)=Rab(3)*fr(4,1)/(two*drab**2)
             ffr(3,11)=Rab(3)*fr(4,2)/(two*drab**2)
             ffr(3,12)=-one/(four*drab)+Rab(3)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(4,1:3)=ffr(1:3,4)
             ffr(4,4)=one/(four*drab)+Rab(1)*fr(2,1)/(two*drab**2)
             ffr(4,5)=Rab(1)*fr(2,2)/(two*drab**2)
             ffr(4,6)=Rab(1)*fr(2,3)/(two*drab**2)
             ffr(4,7)=-one/(four*drab)+Rab(1)*fr(3,1)/(two*drab**2)
             ffr(4,8)=Rab(1)*fr(3,2)/(two*drab**2)
             ffr(4,9)=Rab(1)*fr(3,3)/(two*drab**2)
             ffr(4,10)=-one/(four*drab)+Rab(1)*fr(4,1)/(two*drab**2)
             ffr(4,11)=Rab(1)*fr(4,2)/(two*drab**2)
             ffr(4,12)=Rab(1)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(5,1:3)=ffr(1:3,5)
             ffr(5,4)=Rab(2)*fr(2,1)/(two*drab**2)
             ffr(5,5)=one/(four*drab)+Rab(2)*fr(2,2)/(two*drab**2)
             ffr(5,6)=Rab(2)*fr(2,3)/(two*drab**2)
             ffr(5,7)=Rab(2)*fr(3,1)/(two*drab**2)
             ffr(5,8)=-one/(four*drab)+Rab(2)*fr(3,2)/(two*drab**2)
             ffr(5,9)=Rab(2)*fr(3,3)/(two*drab**2)
             ffr(5,10)=Rab(2)*fr(4,1)/(two*drab**2)
             ffr(5,11)=-one/(four*drab)+Rab(2)*fr(4,2)/(two*drab**2)
             ffr(5,12)=Rab(2)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(6,1:3)=ffr(1:3,6)
             ffr(6,4)=Rab(3)*fr(2,1)/(two*drab**2)
             ffr(6,5)=Rab(3)*fr(2,2)/(two*drab**2)
             ffr(6,6)=one/(four*drab)+Rab(3)*fr(2,3)/(two*drab**2)
             ffr(6,7)=Rab(3)*fr(3,1)/(two*drab**2)
             ffr(6,8)=Rab(3)*fr(3,2)/(two*drab**2)
             ffr(6,9)=-one/(four*drab)+Rab(3)*fr(3,3)/(two*drab**2)
             ffr(6,10)=Rab(3)*fr(4,1)/(two*drab**2)
             ffr(6,11)=Rab(3)*fr(4,2)/(two*drab**2)
             ffr(6,12)=-one/(four*drab)+Rab(3)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(7,1:6)=ffr(1:6,7)
             ffr(7,7)=one/(four*drab)-Rab(1)*fr(3,1)/(two*drab**2)
             ffr(7,8)=-Rab(1)*fr(3,2)/(two*drab**2)
             ffr(7,9)=-Rab(1)*fr(3,3)/(two*drab**2)
             ffr(7,10)=one/(four*drab)-Rab(1)*fr(4,1)/(two*drab**2)
             ffr(7,11)=-Rab(1)*fr(4,2)/(two*drab**2)
             ffr(7,12)=-Rab(1)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(8,1:6)=ffr(1:6,8)
             ffr(8,7)=-Rab(2)*fr(3,1)/(two*drab**2)
             ffr(8,8)=one/(four*drab)-Rab(2)*fr(3,2)/(two*drab**2)
             ffr(8,9)=-Rab(2)*fr(3,3)/(two*drab**2)
             ffr(8,10)=-Rab(2)*fr(4,1)/(two*drab**2)
             ffr(8,11)=one/(four*drab)-Rab(2)*fr(4,2)/(two*drab**2)
             ffr(8,12)=-Rab(2)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(9,1:6)=ffr(1:6,9)
             ffr(9,7)=-Rab(3)*fr(3,1)/(two*drab**2)
             ffr(9,8)=-Rab(3)*fr(3,2)/(two*drab**2)
             ffr(9,9)=one/(four*drab)-Rab(3)*fr(3,3)/(two*drab**2)
             ffr(9,10)=-Rab(3)*fr(4,1)/(two*drab**2)
             ffr(9,11)=-Rab(3)*fr(4,2)/(two*drab**2)
             ffr(9,12)=one/(four*drab)-Rab(3)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(10,1:9)=ffr(1:9,10)
             ffr(10,10)=one/(four*drab)-Rab(1)*fr(4,1)/(two*drab**2)
             ffr(10,11)=-Rab(1)*fr(4,2)/(two*drab**2)
             ffr(10,12)=-Rab(1)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(11,1:9)=ffr(1:9,11)
             ffr(11,10)=-Rab(2)*fr(4,1)/(two*drab**2)
             ffr(11,11)=one/(four*drab)-Rab(2)*fr(4,2)/(two*drab**2)
             ffr(11,12)=-Rab(2)*fr(4,3)/(two*drab**2)
             !...............................
             ffr(12,1:9)=ffr(1:9,12)
             ffr(12,10)=-Rab(3)*fr(4,1)/(two*drab**2)
             ffr(12,11)=-Rab(3)*fr(4,2)/(two*drab**2)
             ffr(12,12)=one/(four*drab)-Rab(3)*fr(4,3)/(two*drab**2)
             !...............................
             !...............................
             d_sin(1:3)=-fx(1,:)*cos_xci/(sin(xci)**2)
             d_sin(4:6)=-fx(2,:)*cos_xci/(sin(xci)**2)
             d_sin(7:9)=-fx(3,:)*cos_xci/(sin(xci)**2)
             d_sin(10:12)=-fx(4,:)*cos_xci/(sin(xci)**2)

             rr=dot_product(Ra,Rb)

             drr(1:3)=-Rb
             drr(4:6)=Rb
             drr(7:9)=-Ra
             drr(10:12)=Ra

             ddra(1:3)=-Ra/dra
             ddra(4:6)=Ra/dra
             ddra(7:9)=zero
             ddra(10:12)=zero

             ddrb(1:3)=zero
             ddrb(4:6)=zero
             ddrb(7:9)=-Rb/drb
             ddrb(10:12)=Rb/drb
             !...............................
             ffx(1,1)=fx(1,1)*sin(xci)*d_sin(1)+ &
                  (Rb(1)*(-ddra(1)/(dra**2*drb)-ddrb(1)/(dra*drb**2))+ &
                  rr/(dra**3*drb)-Ra(1)*drr(1)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(1)/(dra**4*drb)-ddrb(1)/(dra**3*drb**2)))/sin(xci)
             ffx(1,2)=fx(1,1)*sin(xci)*d_sin(2)+ &
                  (Rb(1)*(-ddra(2)/(dra**2*drb)-ddrb(2)/(dra*drb**2))- &
                  Ra(1)*drr(2)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(2)/(dra**4*drb)-ddrb(2)/(dra**3*drb**2)))/sin(xci)
             ffx(1,3)=fx(1,1)*sin(xci)*d_sin(3)+ &
                  (Rb(1)*(-ddra(3)/(dra**2*drb)-ddrb(3)/(dra*drb**2))- &
                  Ra(1)*drr(3)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(3)/(dra**4*drb)-ddrb(3)/(dra**3*drb**2)))/sin(xci)
             ffx(1,4)=fx(1,1)*sin(xci)*d_sin(4)+ &
                  (Rb(1)*(-ddra(4)/(dra**2*drb)-ddrb(4)/(dra*drb**2))- &
                  rr/(dra**3*drb)-Ra(1)*drr(4)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(4)/(dra**4*drb)-ddrb(4)/(dra**3*drb**2)))/sin(xci)
             ffx(1,5)=fx(1,1)*sin(xci)*d_sin(5)+ &
                  (Rb(1)*(-ddra(5)/(dra**2*drb)-ddrb(5)/(dra*drb**2))- &
                  Ra(1)*drr(5)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(5)/(dra**4*drb)-ddrb(5)/(dra**3*drb**2)))/sin(xci)
             ffx(1,6)=fx(1,1)*sin(xci)*d_sin(6)+ &
                  (Rb(1)*(-ddra(6)/(dra**2*drb)-ddrb(6)/(dra*drb**2))- &
                  Ra(1)*drr(6)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(6)/(dra**4*drb)-ddrb(6)/(dra**3*drb**2)))/sin(xci)
             ffx(1,7)=fx(1,1)*sin(xci)*d_sin(7)+ &
                  (-one/(dra*drb)+Rb(1)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))- &
                  Ra(1)*drr(7)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(7)/(dra**4*drb)-ddrb(7)/(dra**3*drb**2)))/sin(xci)
             ffx(1,8)=fx(1,1)*sin(xci)*d_sin(8)+ &
                  (Rb(1)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))- &
                  Ra(1)*drr(8)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(8)/(dra**4*drb)-ddrb(8)/(dra**3*drb**2)))/sin(xci)
             ffx(1,9)=fx(1,1)*sin(xci)*d_sin(9)+ &
                  (Rb(1)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))- &
                  Ra(1)*drr(9)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(9)/(dra**4*drb)-ddrb(9)/(dra**3*drb**2)))/sin(xci)
             ffx(1,10)=fx(1,1)*sin(xci)*d_sin(10)+ &
                  (one/(dra*drb)+Rb(1)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))- &
                  Ra(1)*drr(10)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(10)/(dra**4*drb)-ddrb(10)/(dra**3*drb**2)))/sin(xci)
             ffx(1,11)=fx(1,1)*sin(xci)*d_sin(11)+ &
                  (Rb(1)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))- &
                  Ra(1)*drr(11)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(11)/(dra**4*drb)-ddrb(11)/(dra**3*drb**2)))/sin(xci)
             ffx(1,12)=fx(1,1)*sin(xci)*d_sin(12)+ &
                  (Rb(1)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))- &
                  Ra(1)*drr(12)/(dra**3*drb)- &
                  Ra(1)*rr*(-three*ddra(12)/(dra**4*drb)-ddrb(12)/(dra**3*drb**2)))/sin(xci)
             !...............................
             ffx(2,1)=fx(1,2)*sin(xci)*d_sin(1)+ &
                  (Rb(2)*(-ddra(1)/(dra**2*drb)-ddrb(1)/(dra*drb**2))- &
                  Ra(2)*drr(1)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(1)/(dra**4*drb)-ddrb(1)/(dra**3*drb**2)))/sin(xci)
             ffx(2,2)=fx(1,2)*sin(xci)*d_sin(2)+ &
                  (Rb(2)*(-ddra(2)/(dra**2*drb)-ddrb(2)/(dra*drb**2))+ &
                  rr/(dra**3*drb)-Ra(2)*drr(2)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(2)/(dra**4*drb)-ddrb(2)/(dra**3*drb**2)))/sin(xci)
             ffx(2,3)=fx(1,2)*sin(xci)*d_sin(3)+ &
                  (Rb(2)*(-ddra(3)/(dra**2*drb)-ddrb(3)/(dra*drb**2))- &
                  Ra(2)*drr(3)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(3)/(dra**4*drb)-ddrb(3)/(dra**3*drb**2)))/sin(xci)
             ffx(2,4)=fx(1,2)*sin(xci)*d_sin(4)+ &
                  (Rb(2)*(-ddra(4)/(dra**2*drb)-ddrb(4)/(dra*drb**2))- &
                  Ra(2)*drr(4)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(4)/(dra**4*drb)-ddrb(4)/(dra**3*drb**2)))/sin(xci)
             ffx(2,5)=fx(1,2)*sin(xci)*d_sin(5)+ &
                  (Rb(2)*(-ddra(5)/(dra**2*drb)-ddrb(5)/(dra*drb**2))- &
                  rr/(dra**3*drb)-Ra(2)*drr(5)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(5)/(dra**4*drb)-ddrb(5)/(dra**3*drb**2)))/sin(xci)
             ffx(2,6)=fx(1,2)*sin(xci)*d_sin(6)+ &
                  (Rb(2)*(-ddra(6)/(dra**2*drb)-ddrb(6)/(dra*drb**2))- &
                  Ra(2)*drr(6)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(6)/(dra**4*drb)-ddrb(6)/(dra**3*drb**2)))/sin(xci)
             ffx(2,7)=fx(1,2)*sin(xci)*d_sin(7)+ &
                  (Rb(2)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))- &
                  Ra(2)*drr(7)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(7)/(dra**4*drb)-ddrb(7)/(dra**3*drb**2)))/sin(xci)
             ffx(2,8)=fx(1,2)*sin(xci)*d_sin(8)+ &
                  (-one/(dra*drb)+Rb(2)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))- &
                  Ra(2)*drr(8)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(8)/(dra**4*drb)-ddrb(8)/(dra**3*drb**2)))/sin(xci)
             ffx(2,9)=fx(1,2)*sin(xci)*d_sin(9)+ &
                  (Rb(2)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))- &
                  Ra(2)*drr(9)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(9)/(dra**4*drb)-ddrb(9)/(dra**3*drb**2)))/sin(xci)
             ffx(2,10)=fx(1,2)*sin(xci)*d_sin(10)+ &
                  (Rb(2)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))- &
                  Ra(2)*drr(10)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(10)/(dra**4*drb)-ddrb(10)/(dra**3*drb**2)))/sin(xci)
             ffx(2,11)=fx(1,2)*sin(xci)*d_sin(11)+ &
                  (one/(dra*drb)+Rb(2)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))- &
                  Ra(2)*drr(11)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(11)/(dra**4*drb)-ddrb(11)/(dra**3*drb**2)))/sin(xci)
             ffx(2,12)=fx(1,2)*sin(xci)*d_sin(12)+ &
                  (Rb(2)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))- &
                  Ra(2)*drr(12)/(dra**3*drb)- &
                  Ra(2)*rr*(-three*ddra(12)/(dra**4*drb)-ddrb(12)/(dra**3*drb**2)))/sin(xci)
             !...............................
             ffx(3,1)=fx(1,3)*sin(xci)*d_sin(1)+ &
                  (Rb(3)*(-ddra(1)/(dra**2*drb)-ddrb(1)/(dra*drb**2))- &
                  Ra(3)*drr(1)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(1)/(dra**4*drb)-ddrb(1)/(dra**3*drb**2)))/sin(xci)
             ffx(3,2)=fx(1,3)*sin(xci)*d_sin(2)+ &
                  (Rb(3)*(-ddra(2)/(dra**2*drb)-ddrb(2)/(dra*drb**2))- &
                  Ra(3)*drr(2)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(2)/(dra**4*drb)-ddrb(2)/(dra**3*drb**2)))/sin(xci)
             ffx(3,3)=fx(1,3)*sin(xci)*d_sin(3)+ &
                  (Rb(3)*(-ddra(3)/(dra**2*drb)-ddrb(3)/(dra*drb**2))+ &
                  rr/(dra**3*drb)-Ra(3)*drr(3)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(3)/(dra**4*drb)-ddrb(3)/(dra**3*drb**2)))/sin(xci)
             ffx(3,4)=fx(1,3)*sin(xci)*d_sin(4)+ &
                  (Rb(3)*(-ddra(4)/(dra**2*drb)-ddrb(4)/(dra*drb**2))- &
                  Ra(3)*drr(4)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(4)/(dra**4*drb)-ddrb(4)/(dra**3*drb**2)))/sin(xci)
             ffx(3,5)=fx(1,3)*sin(xci)*d_sin(5)+ &
                  (Rb(3)*(-ddra(5)/(dra**2*drb)-ddrb(5)/(dra*drb**2))- &
                  Ra(3)*drr(5)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(5)/(dra**4*drb)-ddrb(5)/(dra**3*drb**2)))/sin(xci)
             ffx(3,6)=fx(1,3)*sin(xci)*d_sin(6)+ &
                  (Rb(3)*(-ddra(6)/(dra**2*drb)-ddrb(6)/(dra*drb**2))- &
                  rr/(dra**3*drb)-Ra(3)*drr(6)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(6)/(dra**4*drb)-ddrb(6)/(dra**3*drb**2)))/sin(xci)
             ffx(3,7)=fx(1,3)*sin(xci)*d_sin(7)+ &
                  (Rb(3)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))- &
                  Ra(3)*drr(7)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(7)/(dra**4*drb)-ddrb(7)/(dra**3*drb**2)))/sin(xci)
             ffx(3,8)=fx(1,3)*sin(xci)*d_sin(8)+ &
                  (Rb(3)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))- &
                  Ra(3)*drr(8)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(8)/(dra**4*drb)-ddrb(8)/(dra**3*drb**2)))/sin(xci)
             ffx(3,9)=fx(1,3)*sin(xci)*d_sin(9)+ &
                  (-one/(dra*drb)+Rb(3)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))- &
                  Ra(3)*drr(9)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(9)/(dra**4*drb)-ddrb(9)/(dra**3*drb**2)))/sin(xci)
             ffx(3,10)=fx(1,3)*sin(xci)*d_sin(10)+ &
                  (Rb(3)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))- &
                  Ra(3)*drr(10)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(10)/(dra**4*drb)-ddrb(10)/(dra**3*drb**2)))/sin(xci)
             ffx(3,11)=fx(1,3)*sin(xci)*d_sin(11)+ &
                  (Rb(3)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))- &
                  Ra(3)*drr(11)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(11)/(dra**4*drb)-ddrb(11)/(dra**3*drb**2)))/sin(xci)
             ffx(3,12)=fx(1,3)*sin(xci)*d_sin(12)+ &
                  (one/(dra*drb)+Rb(3)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))- &
                  Ra(3)*drr(12)/(dra**3*drb)- &
                  Ra(3)*rr*(-three*ddra(12)/(dra**4*drb)-ddrb(12)/(dra**3*drb**2)))/sin(xci)
             !...............................
             ffx(4,1:3)=ffx(1:3,4)
             ffx(4,4)=fx(2,1)*sin(xci)*d_sin(4)+ &
                  (-Rb(1)*(-ddra(4)/(dra**2*drb)-ddrb(4)/(dra*drb**2))+ &
                  rr/(dra**3*drb)+Ra(1)*drr(4)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(4)/(dra**4*drb)-ddrb(4)/(dra**3*drb**2)))/sin(xci)
             ffx(4,5)=fx(2,1)*sin(xci)*d_sin(5)+ &
                  (-Rb(1)*(-ddra(5)/(dra**2*drb)-ddrb(5)/(dra*drb**2))+ &
                  Ra(1)*drr(5)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(5)/(dra**4*drb)-ddrb(5)/(dra**3*drb**2)))/sin(xci)
             ffx(4,6)=fx(2,1)*sin(xci)*d_sin(6)+ &
                  (-Rb(1)*(-ddra(6)/(dra**2*drb)-ddrb(6)/(dra*drb**2))+ &
                  Ra(1)*drr(6)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(6)/(dra**4*drb)-ddrb(6)/(dra**3*drb**2)))/sin(xci)
             ffx(4,7)=fx(2,1)*sin(xci)*d_sin(7)+ &
                  (one/(dra*drb)-Rb(1)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))+ &
                  Ra(1)*drr(7)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(7)/(dra**4*drb)-ddrb(7)/(dra**3*drb**2)))/sin(xci)
             ffx(4,8)=fx(2,1)*sin(xci)*d_sin(8)+ &
                  (-Rb(1)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))+ &
                  Ra(1)*drr(8)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(8)/(dra**4*drb)-ddrb(8)/(dra**3*drb**2)))/sin(xci)
             ffx(4,9)=fx(2,1)*sin(xci)*d_sin(9)+ &
                  (-Rb(1)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))+ &
                  Ra(1)*drr(9)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(9)/(dra**4*drb)-ddrb(9)/(dra**3*drb**2)))/sin(xci)
             ffx(4,10)=fx(2,1)*sin(xci)*d_sin(10)+ &
                  (-one/(dra*drb)-Rb(1)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))+ &
                  Ra(1)*drr(10)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(10)/(dra**4*drb)-ddrb(10)/(dra**3*drb**2)))/sin(xci)
             ffx(4,11)=fx(2,1)*sin(xci)*d_sin(11)+ &
                  (-Rb(1)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))+ &
                  Ra(1)*drr(11)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(11)/(dra**4*drb)-ddrb(11)/(dra**3*drb**2)))/sin(xci)
             ffx(4,12)=fx(2,1)*sin(xci)*d_sin(12)+ &
                  (-Rb(1)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))+ &
                  Ra(1)*drr(12)/(dra**3*drb)+ &
                  Ra(1)*rr*(-three*ddra(12)/(dra**4*drb)-ddrb(12)/(dra**3*drb**2)))/sin(xci)
             !...............................
             ffx(5,1:3)=ffx(1:3,5)
             ffx(5,4)=fx(2,2)*sin(xci)*d_sin(4)+ &
                  (-Rb(2)*(-ddra(4)/(dra**2*drb)-ddrb(4)/(dra*drb**2))+ &
                  Ra(2)*drr(4)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(4)/(dra**4*drb)-ddrb(4)/(dra**3*drb**2)))/sin(xci)
             ffx(5,5)=fx(2,2)*sin(xci)*d_sin(5)+ &
                  (-Rb(2)*(-ddra(5)/(dra**2*drb)-ddrb(5)/(dra*drb**2))+ &
                  rr/(dra**3*drb)+Ra(2)*drr(5)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(5)/(dra**4*drb)-ddrb(5)/(dra**3*drb**2)))/sin(xci)
             ffx(5,6)=fx(2,2)*sin(xci)*d_sin(6)+ &
                  (-Rb(2)*(-ddra(6)/(dra**2*drb)-ddrb(6)/(dra*drb**2))+ &
                  Ra(2)*drr(6)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(6)/(dra**4*drb)-ddrb(6)/(dra**3*drb**2)))/sin(xci)
             ffx(5,7)=fx(2,2)*sin(xci)*d_sin(7)+ &
                  (-Rb(2)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))+ &
                  Ra(2)*drr(7)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(7)/(dra**4*drb)-ddrb(7)/(dra**3*drb**2)))/sin(xci)
             ffx(5,8)=fx(2,2)*sin(xci)*d_sin(8)+ &
                  (one/(dra*drb)-Rb(2)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))+ &
                  Ra(2)*drr(8)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(8)/(dra**4*drb)-ddrb(8)/(dra**3*drb**2)))/sin(xci)
             ffx(5,9)=fx(2,2)*sin(xci)*d_sin(9)+ &
                  (-Rb(2)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))+ &
                  Ra(2)*drr(9)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(9)/(dra**4*drb)-ddrb(9)/(dra**3*drb**2)))/sin(xci)
             ffx(5,10)=fx(2,2)*sin(xci)*d_sin(10)+ &
                  (-Rb(2)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))+ &
                  Ra(2)*drr(10)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(10)/(dra**4*drb)-ddrb(10)/(dra**3*drb**2)))/sin(xci)
             ffx(5,11)=fx(2,2)*sin(xci)*d_sin(11)+ &
                  (-one/(dra*drb)-Rb(2)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))+ &
                  Ra(2)*drr(11)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(11)/(dra**4*drb)-ddrb(11)/(dra**3*drb**2)))/sin(xci)
             ffx(5,12)=fx(2,2)*sin(xci)*d_sin(12)+ &
                  (-Rb(2)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))+ &
                  Ra(2)*drr(12)/(dra**3*drb)+ &
                  Ra(2)*rr*(-three*ddra(12)/(dra**4*drb)-ddrb(12)/(dra**3*drb**2)))/sin(xci)
             !...............................
             ffx(6,1:3)=ffx(1:3,6)
             ffx(6,4)=fx(2,3)*sin(xci)*d_sin(4)+ &
                  (-Rb(3)*(-ddra(4)/(dra**2*drb)-ddrb(4)/(dra*drb**2))+ &
                  Ra(3)*drr(4)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(4)/(dra**4*drb)-ddrb(4)/(dra**3*drb**2)))/sin(xci)
             ffx(6,5)=fx(2,3)*sin(xci)*d_sin(5)+ &
                  (-Rb(3)*(-ddra(5)/(dra**2*drb)-ddrb(5)/(dra*drb**2))+ &
                  Ra(3)*drr(5)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(5)/(dra**4*drb)-ddrb(5)/(dra**3*drb**2)))/sin(xci)
             ffx(6,6)=fx(2,3)*sin(xci)*d_sin(6)+ &
                  (-Rb(3)*(-ddra(6)/(dra**2*drb)-ddrb(6)/(dra*drb**2))+ &
                  rr/(dra**3*drb)+Ra(3)*drr(6)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(6)/(dra**4*drb)-ddrb(6)/(dra**3*drb**2)))/sin(xci)
             ffx(6,7)=fx(2,3)*sin(xci)*d_sin(7)+ &
                  (-Rb(3)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))+ &
                  Ra(3)*drr(7)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(7)/(dra**4*drb)-ddrb(7)/(dra**3*drb**2)))/sin(xci)
             ffx(6,8)=fx(2,3)*sin(xci)*d_sin(8)+ &
                  (-Rb(3)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))+ &
                  Ra(3)*drr(8)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(8)/(dra**4*drb)-ddrb(8)/(dra**3*drb**2)))/sin(xci)
             ffx(6,9)=fx(2,3)*sin(xci)*d_sin(9)+ &
                  (one/(dra*drb)-Rb(3)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))+ &
                  Ra(3)*drr(9)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(9)/(dra**4*drb)-ddrb(9)/(dra**3*drb**2)))/sin(xci)
             ffx(6,10)=fx(2,3)*sin(xci)*d_sin(10)+ &
                  (-Rb(3)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))+ &
                  Ra(3)*drr(10)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(10)/(dra**4*drb)-ddrb(10)/(dra**3*drb**2)))/sin(xci)
             ffx(6,11)=fx(2,3)*sin(xci)*d_sin(11)+ &
                  (-Rb(3)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))+ &
                  Ra(3)*drr(11)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(11)/(dra**4*drb)-ddrb(11)/(dra**3*drb**2)))/sin(xci)
             ffx(6,12)=fx(2,3)*sin(xci)*d_sin(12)+ &
                  (-one/(dra*drb)-Rb(3)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))+ &
                  Ra(3)*drr(12)/(dra**3*drb)+ &
                  Ra(3)*rr*(-three*ddra(12)/(dra**4*drb)-ddrb(12)/(dra**3*drb**2)))/sin(xci)
             !...............................
             ffx(7,1:6)=ffx(1:6,7)
             ffx(7,7)=fx(3,1)*sin(xci)*d_sin(7)+ &
                  (Ra(1)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))+ &
                  rr/(drb**3*dra)-Rb(1)*drr(7)/(drb**3*dra)- &
                  Rb(1)*rr*(-three*ddrb(7)/(drb**4*dra)-ddra(7)/(drb**3*dra**2)))/sin(xci)
             ffx(7,8)=fx(3,1)*sin(xci)*d_sin(8)+ &
                  (Ra(1)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))- &
                  Rb(1)*drr(8)/(drb**3*dra)- &
                  Rb(1)*rr*(-three*ddrb(8)/(drb**4*dra)-ddra(8)/(drb**3*dra**2)))/sin(xci)
             ffx(7,9)=fx(3,1)*sin(xci)*d_sin(9)+ &
                  (Ra(1)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))- &
                  Rb(1)*drr(9)/(drb**3*dra)- &
                  Rb(1)*rr*(-three*ddrb(9)/(drb**4*dra)-ddra(9)/(drb**3*dra**2)))/sin(xci)
             ffx(7,10)=fx(3,1)*sin(xci)*d_sin(10)+ &
                  (Ra(1)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))- &
                  rr/(drb**3*dra)-Rb(1)*drr(10)/(drb**3*dra)- &
                  Rb(1)*rr*(-three*ddrb(10)/(drb**4*dra)-ddra(10)/(drb**3*dra**2)))/sin(xci)
             ffx(7,11)=fx(3,1)*sin(xci)*d_sin(11)+ &
                  (Ra(1)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))- &
                  Rb(1)*drr(11)/(drb**3*dra)- &
                  Rb(1)*rr*(-three*ddrb(11)/(drb**4*dra)-ddra(11)/(drb**3*dra**2)))/sin(xci)
             ffx(7,12)=fx(3,1)*sin(xci)*d_sin(12)+ &
                  (Ra(1)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))- &
                  Rb(1)*drr(12)/(drb**3*dra)- &
                  Rb(1)*rr*(-three*ddrb(12)/(drb**4*dra)-ddra(12)/(drb**3*dra**2)))/sin(xci)
             !...............................
             ffx(8,1:6)=ffx(1:6,8)
             ffx(8,7)=fx(3,2)*sin(xci)*d_sin(7)+ &
                  (Ra(2)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))- &
                  Rb(2)*drr(7)/(drb**3*dra)- &
                  Rb(2)*rr*(-three*ddrb(7)/(drb**4*dra)-ddra(7)/(drb**3*dra**2)))/sin(xci)
             ffx(8,8)=fx(3,2)*sin(xci)*d_sin(8)+ &
                  (Ra(2)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))+ &
                  rr/(drb**3*dra)-Rb(2)*drr(8)/(drb**3*dra)- &
                  Rb(2)*rr*(-three*ddrb(8)/(drb**4*dra)-ddra(8)/(drb**3*dra**2)))/sin(xci)
             ffx(8,9)=fx(3,2)*sin(xci)*d_sin(9)+ &
                  (Ra(2)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))- &
                  Rb(2)*drr(9)/(drb**3*dra)- &
                  Rb(2)*rr*(-three*ddrb(9)/(drb**4*dra)-ddra(9)/(drb**3*dra**2)))/sin(xci)
             ffx(8,10)=fx(3,2)*sin(xci)*d_sin(10)+ &
                  (Ra(2)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))- &
                  Rb(2)*drr(10)/(drb**3*dra)- &
                  Rb(2)*rr*(-three*ddrb(10)/(drb**4*dra)-ddra(10)/(drb**3*dra**2)))/sin(xci)
             ffx(8,11)=fx(3,2)*sin(xci)*d_sin(11)+ &
                  (Ra(2)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))- &
                  rr/(drb**3*dra)-Rb(2)*drr(11)/(drb**3*dra)- &
                  Rb(2)*rr*(-three*ddrb(11)/(drb**4*dra)-ddra(11)/(drb**3*dra**2)))/sin(xci)
             ffx(8,12)=fx(3,2)*sin(xci)*d_sin(12)+ &
                  (Ra(2)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))- &
                  Rb(2)*drr(12)/(drb**3*dra)- &
                  Rb(2)*rr*(-three*ddrb(12)/(drb**4*dra)-ddra(12)/(drb**3*dra**2)))/sin(xci)
             !...............................
             ffx(9,1:6)=ffx(1:6,9)
             ffx(9,7)=fx(3,3)*sin(xci)*d_sin(7)+ &
                  (Ra(3)*(-ddra(7)/(dra**2*drb)-ddrb(7)/(dra*drb**2))- &
                  Rb(3)*drr(7)/(drb**3*dra)- &
                  Rb(3)*rr*(-three*ddrb(7)/(drb**4*dra)-ddra(7)/(drb**3*dra**2)))/sin(xci)
             ffx(9,8)=fx(3,3)*sin(xci)*d_sin(8)+ &
                  (Ra(3)*(-ddra(8)/(dra**2*drb)-ddrb(8)/(dra*drb**2))- &
                  Rb(3)*drr(8)/(drb**3*dra)- &
                  Rb(3)*rr*(-three*ddrb(8)/(drb**4*dra)-ddra(8)/(drb**3*dra**2)))/sin(xci)
             ffx(9,9)=fx(3,3)*sin(xci)*d_sin(9)+ &
                  (Ra(3)*(-ddra(9)/(dra**2*drb)-ddrb(9)/(dra*drb**2))+ &
                  rr/(drb**3*dra)-Rb(3)*drr(9)/(drb**3*dra)- &
                  Rb(3)*rr*(-three*ddrb(9)/(drb**4*dra)-ddra(9)/(drb**3*dra**2)))/sin(xci)
             ffx(9,10)=fx(3,3)*sin(xci)*d_sin(10)+ &
                  (Ra(3)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))- &
                  Rb(3)*drr(10)/(drb**3*dra)- &
                  Rb(3)*rr*(-three*ddrb(10)/(drb**4*dra)-ddra(10)/(drb**3*dra**2)))/sin(xci)
             ffx(9,11)=fx(3,3)*sin(xci)*d_sin(11)+ &
                  (Ra(3)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))- &
                  Rb(3)*drr(11)/(drb**3*dra)- &
                  Rb(3)*rr*(-three*ddrb(11)/(drb**4*dra)-ddra(11)/(drb**3*dra**2)))/sin(xci)
             ffx(9,12)=fx(3,3)*sin(xci)*d_sin(12)+ &
                  (Ra(3)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))- &
                  rr/(drb**3*dra)-Rb(3)*drr(12)/(drb**3*dra)- &
                  Rb(3)*rr*(-three*ddrb(12)/(drb**4*dra)-ddra(12)/(drb**3*dra**2)))/sin(xci)
             !...............................
             ffx(10,1:9)=ffx(1:9,10)
             ffx(10,10)=fx(4,1)*sin(xci)*d_sin(10)+ &
                  (-Ra(1)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))+ &
                  rr/(drb**3*dra)+Rb(1)*drr(10)/(drb**3*dra)+ &
                  Rb(1)*rr*(-three*ddrb(10)/(drb**4*dra)-ddra(10)/(drb**3*dra**2)))/sin(xci)
             ffx(10,11)=fx(4,1)*sin(xci)*d_sin(11)+ &
                  (-Ra(1)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))+ &
                  Rb(1)*drr(11)/(drb**3*dra)+ &
                  Rb(1)*rr*(-three*ddrb(11)/(drb**4*dra)-ddra(11)/(drb**3*dra**2)))/sin(xci)
             ffx(10,12)=fx(4,1)*sin(xci)*d_sin(12)+ &
                  (-Ra(1)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))+ &
                  Rb(1)*drr(12)/(drb**3*dra)+ &
                  Rb(1)*rr*(-three*ddrb(12)/(drb**4*dra)-ddra(12)/(drb**3*dra**2)))/sin(xci)
             !...............................
             ffx(11,1:9)=ffx(1:9,11)
             ffx(11,10)=fx(4,2)*sin(xci)*d_sin(10)+ &
                  (-Ra(2)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))+ &
                  Rb(2)*drr(10)/(drb**3*dra)+ &
                  Rb(2)*rr*(-three*ddrb(10)/(drb**4*dra)-ddra(10)/(drb**3*dra**2)))/sin(xci)
             ffx(11,11)=fx(4,2)*sin(xci)*d_sin(11)+ &
                  (-Ra(2)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))+ &
                  rr/(drb**3*dra)+Rb(2)*drr(11)/(drb**3*dra)+ &
                  Rb(2)*rr*(-three*ddrb(11)/(drb**4*dra)-ddra(11)/(drb**3*dra**2)))/sin(xci)
             ffx(11,12)=fx(4,2)*sin(xci)*d_sin(12)+ &
                  (-Ra(2)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))+ &
                  Rb(2)*drr(12)/(drb**3*dra)+ &
                  Rb(2)*rr*(-three*ddrb(12)/(drb**4*dra)-ddra(12)/(drb**3*dra**2)))/sin(xci)
             !...............................
             ffx(12,1:9)=ffx(1:9,12)
             ffx(12,10)=fx(4,3)*sin(xci)*d_sin(10)+ &
                  (-Ra(3)*(-ddra(10)/(dra**2*drb)-ddrb(10)/(dra*drb**2))+ &
                  Rb(3)*drr(10)/(drb**3*dra)+ &
                  Rb(3)*rr*(-three*ddrb(10)/(drb**4*dra)-ddra(10)/(drb**3*dra**2)))/sin(xci)
             ffx(12,11)=fx(4,3)*sin(xci)*d_sin(11)+ &
                  (-Ra(3)*(-ddra(11)/(dra**2*drb)-ddrb(11)/(dra*drb**2))+ &
                  Rb(3)*drr(11)/(drb**3*dra)+ &
                  Rb(3)*rr*(-three*ddrb(11)/(drb**4*dra)-ddra(11)/(drb**3*dra**2)))/sin(xci)
             ffx(12,12)=fx(4,3)*sin(xci)*d_sin(12)+ &
                  (-Ra(3)*(-ddra(12)/(dra**2*drb)-ddrb(12)/(dra*drb**2))+ &
                  rr/(drb**3*dra)+Rb(3)*drr(12)/(drb**3*dra)+ &
                  Rb(3)*rr*(-three*ddrb(12)/(drb**4*dra)-ddra(12)/(drb**3*dra**2)))/sin(xci)
             !...............................
             !...............................
             d_sin(1:3)=-fa(1,:)*cos_alpha/(sin(alpha)**2)
             d_sin(4:6)=-fa(2,:)*cos_alpha/(sin(alpha)**2)
             d_sin(7:9)=-fa(3,:)*cos_alpha/(sin(alpha)**2)
             d_sin(10:12)=-fa(4,:)*cos_alpha/(sin(alpha)**2)

             rr=dot_product(Ra,Rab)

             drr(1:3)=-Rab-Ra*half
             drr(4:6)=Rab-Ra*half
             drr(7:9)=Ra*half
             drr(10:12)=Ra*half

             rra=r4+r3-two*r1
             rrb=r4+r3-two*r2

             ffa(1,1)=fa(1,1)*sin(alpha)*d_sin(1)- &
                  (one/(dra*drab)-rra(1)*(-ddra(1)/(two*dra**2*drab)-fr(1,1)/(two*dra*drab**2))- &
                  rr/(dra**3*drab)+Ra(1)*drr(1)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(1)/(dra**4*drab)-fr(1,1)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(1)*drr(1)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(1)/(two*dra**2*drab**3)-three*fr(1,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,2)=fa(1,1)*sin(alpha)*d_sin(2)- &
                  (-rra(1)*(-ddra(2)/(two*dra**2*drab)-fr(1,2)/(two*dra*drab**2))+ &
                  Ra(1)*drr(2)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(2)/(dra**4*drab)-fr(1,2)/(dra**3*drab**2))+ &
                  Rab(1)*drr(2)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(2)/(two*dra**2*drab**3)-three*fr(1,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,3)=fa(1,1)*sin(alpha)*d_sin(3)- &
                  (-rra(1)*(-ddra(3)/(two*dra**2*drab)-fr(1,3)/(two*dra*drab**2))+ &
                  Ra(1)*drr(3)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(3)/(dra**4*drab)-fr(1,3)/(dra**3*drab**2))+ &
                  Rab(1)*drr(3)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(3)/(two*dra**2*drab**3)-three*fr(1,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,4)=fa(1,1)*sin(alpha)*d_sin(4)- &
                  (-rra(1)*(-ddra(4)/(two*dra**2*drab)-fr(2,1)/(two*dra*drab**2))+ &
                  rr/(dra**3*drab)+Ra(1)*drr(4)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(4)/(dra**4*drab)-fr(2,1)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(1)*drr(4)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(4)/(two*dra**2*drab**3)-three*fr(2,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,5)=fa(1,1)*sin(alpha)*d_sin(5)- &
                  (-rra(1)*(-ddra(5)/(two*dra**2*drab)-fr(2,2)/(two*dra*drab**2))+ &
                  Ra(1)*drr(5)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(5)/(dra**4*drab)-fr(2,2)/(dra**3*drab**2))+ &
                  Rab(1)*drr(5)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(5)/(two*dra**2*drab**3)-three*fr(2,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,6)=fa(1,1)*sin(alpha)*d_sin(6)- &
                  (-rra(1)*(-ddra(6)/(two*dra**2*drab)-fr(2,3)/(two*dra*drab**2))+ &
                  Ra(1)*drr(6)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(6)/(dra**4*drab)-fr(2,3)/(dra**3*drab**2))+ &
                  Rab(1)*drr(6)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(6)/(two*dra**2*drab**3)-three*fr(2,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,7)=fa(1,1)*sin(alpha)*d_sin(7)- &
                  (-one/(two*dra*drab)-rra(1)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))+ &
                  Ra(1)*drr(7)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(7)/(dra**4*drab)-fr(3,1)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(1)*drr(7)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,8)=fa(1,1)*sin(alpha)*d_sin(8)- &
                  (-rra(1)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))+ &
                  Ra(1)*drr(8)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(8)/(dra**4*drab)-fr(3,2)/(dra**3*drab**2))+ &
                  Rab(1)*drr(8)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,9)=fa(1,1)*sin(alpha)*d_sin(9)- &
                  (-rra(1)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))+ &
                  Ra(1)*drr(9)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(9)/(dra**4*drab)-fr(3,3)/(dra**3*drab**2))+ &
                  Rab(1)*drr(9)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,10)=fa(1,1)*sin(alpha)*d_sin(10)- &
                  (-one/(two*dra*drab)-rra(1)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))+ &
                  Ra(1)*drr(10)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(10)/(dra**4*drab)-fr(4,1)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(1)*drr(10)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,11)=fa(1,1)*sin(alpha)*d_sin(11)- &
                  (-rra(1)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))+ &
                  Ra(1)*drr(11)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(11)/(dra**4*drab)-fr(4,2)/(dra**3*drab**2))+ &
                  Rab(1)*drr(11)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(1,12)=fa(1,1)*sin(alpha)*d_sin(12)- &
                  (-rra(1)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))+ &
                  Ra(1)*drr(12)/(dra**3*drab)+ &
                  Ra(1)*rr*(-three*ddra(12)/(dra**4*drab)-fr(4,3)/(dra**3*drab**2))+ &
                  Rab(1)*drr(12)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(2,1)=fa(1,2)*sin(alpha)*d_sin(1)- &
                  (-rra(2)*(-ddra(1)/(two*dra**2*drab)-fr(1,1)/(two*dra*drab**2))+ &
                  Ra(2)*drr(1)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(1)/(dra**4*drab)-fr(1,1)/(dra**3*drab**2))+ &
                  Rab(2)*drr(1)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(1)/(two*dra**2*drab**3)-three*fr(1,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,2)=fa(1,2)*sin(alpha)*d_sin(2)- &
                  (one/(dra*drab)-rra(2)*(-ddra(2)/(two*dra**2*drab)-fr(1,2)/(two*dra*drab**2))- &
                  rr/(dra**3*drab)+Ra(2)*drr(2)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(2)/(dra**4*drab)-fr(1,2)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(2)*drr(2)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(2)/(two*dra**2*drab**3)-three*fr(1,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,3)=fa(1,2)*sin(alpha)*d_sin(3)- &
                  (-rra(2)*(-ddra(3)/(two*dra**2*drab)-fr(1,3)/(two*dra*drab**2))+ &
                  Ra(2)*drr(3)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(3)/(dra**4*drab)-fr(1,3)/(dra**3*drab**2))+ &
                  Rab(2)*drr(3)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(3)/(two*dra**2*drab**3)-three*fr(1,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,4)=fa(1,2)*sin(alpha)*d_sin(4)- &
                  (-rra(2)*(-ddra(4)/(two*dra**2*drab)-fr(2,1)/(two*dra*drab**2))+ &
                  Ra(2)*drr(4)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(4)/(dra**4*drab)-fr(2,1)/(dra**3*drab**2))+ &
                  Rab(2)*drr(4)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(4)/(two*dra**2*drab**3)-three*fr(2,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,5)=fa(1,2)*sin(alpha)*d_sin(5)- &
                  (-rra(2)*(-ddra(5)/(two*dra**2*drab)-fr(2,2)/(two*dra*drab**2))+ &
                  rr/(dra**3*drab)+Ra(2)*drr(5)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(5)/(dra**4*drab)-fr(2,2)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(2)*drr(5)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(5)/(two*dra**2*drab**3)-three*fr(2,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,6)=fa(1,2)*sin(alpha)*d_sin(6)- &
                  (-rra(2)*(-ddra(6)/(two*dra**2*drab)-fr(2,3)/(two*dra*drab**2))+ &
                  Ra(2)*drr(6)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(6)/(dra**4*drab)-fr(2,3)/(dra**3*drab**2))+ &
                  Rab(2)*drr(6)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(6)/(two*dra**2*drab**3)-three*fr(2,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,7)=fa(1,2)*sin(alpha)*d_sin(7)- &
                  (-rra(2)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))+ &
                  Ra(2)*drr(7)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(7)/(dra**4*drab)-fr(3,1)/(dra**3*drab**2))+ &
                  Rab(2)*drr(7)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,8)=fa(1,2)*sin(alpha)*d_sin(8)- &
                  (-one/(two*dra*drab)-rra(2)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))+ &
                  Ra(2)*drr(8)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(8)/(dra**4*drab)-fr(3,2)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(2)*drr(8)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,9)=fa(1,2)*sin(alpha)*d_sin(9)- &
                  (-rra(2)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))+ &
                  Ra(2)*drr(9)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(9)/(dra**4*drab)-fr(3,3)/(dra**3*drab**2))+ &
                  Rab(2)*drr(9)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,10)=fa(1,2)*sin(alpha)*d_sin(10)- &
                  (-rra(2)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))+ &
                  Ra(2)*drr(10)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(10)/(dra**4*drab)-fr(4,1)/(dra**3*drab**2))+ &
                  Rab(2)*drr(10)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,11)=fa(1,2)*sin(alpha)*d_sin(11)- &
                  (-one/(two*dra*drab)-rra(2)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))+ &
                  Ra(2)*drr(11)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(11)/(dra**4*drab)-fr(4,2)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(2)*drr(11)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(2,12)=fa(1,2)*sin(alpha)*d_sin(12)- &
                  (-rra(2)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))+ &
                  Ra(2)*drr(12)/(dra**3*drab)+ &
                  Ra(2)*rr*(-three*ddra(12)/(dra**4*drab)-fr(4,3)/(dra**3*drab**2))+ &
                  Rab(2)*drr(12)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(3,1)=fa(1,3)*sin(alpha)*d_sin(1)- &
                  (-rra(3)*(-ddra(1)/(two*dra**2*drab)-fr(1,1)/(two*dra*drab**2))+ &
                  Ra(3)*drr(1)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(1)/(dra**4*drab)-fr(1,1)/(dra**3*drab**2))+ &
                  Rab(3)*drr(1)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(1)/(two*dra**2*drab**3)-three*fr(1,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,2)=fa(1,3)*sin(alpha)*d_sin(2)- &
                  (-rra(3)*(-ddra(2)/(two*dra**2*drab)-fr(1,2)/(two*dra*drab**2))+ &
                  Ra(3)*drr(2)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(2)/(dra**4*drab)-fr(1,2)/(dra**3*drab**2))+ &
                  Rab(3)*drr(2)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(2)/(two*dra**2*drab**3)-three*fr(1,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,3)=fa(1,3)*sin(alpha)*d_sin(3)- &
                  (one/(dra*drab)-rra(3)*(-ddra(3)/(two*dra**2*drab)-fr(1,3)/(two*dra*drab**2))- &
                  rr/(dra**3*drab)+Ra(3)*drr(3)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(3)/(dra**4*drab)-fr(1,3)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(3)*drr(3)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(3)/(two*dra**2*drab**3)-three*fr(1,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,4)=fa(1,3)*sin(alpha)*d_sin(4)- &
                  (-rra(3)*(-ddra(4)/(two*dra**2*drab)-fr(2,1)/(two*dra*drab**2))+ &
                  Ra(3)*drr(4)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(4)/(dra**4*drab)-fr(2,1)/(dra**3*drab**2))+ &
                  Rab(3)*drr(4)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(4)/(two*dra**2*drab**3)-three*fr(2,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,5)=fa(1,3)*sin(alpha)*d_sin(5)- &
                  (-rra(3)*(-ddra(5)/(two*dra**2*drab)-fr(2,2)/(two*dra*drab**2))+ &
                  Ra(3)*drr(5)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(5)/(dra**4*drab)-fr(2,2)/(dra**3*drab**2))+ &
                  Rab(3)*drr(5)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(5)/(two*dra**2*drab**3)-three*fr(2,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,6)=fa(1,3)*sin(alpha)*d_sin(6)- &
                  (-rra(3)*(-ddra(6)/(two*dra**2*drab)-fr(2,3)/(two*dra*drab**2))+ &
                  rr/(dra**3*drab)+Ra(3)*drr(6)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(6)/(dra**4*drab)-fr(2,3)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(3)*drr(6)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(6)/(two*dra**2*drab**3)-three*fr(2,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,7)=fa(1,3)*sin(alpha)*d_sin(7)- &
                  (-rra(3)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))+ &
                  Ra(3)*drr(7)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(7)/(dra**4*drab)-fr(3,1)/(dra**3*drab**2))+ &
                  Rab(3)*drr(7)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,8)=fa(1,3)*sin(alpha)*d_sin(8)- &
                  (-rra(3)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))+ &
                  Ra(3)*drr(8)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(8)/(dra**4*drab)-fr(3,2)/(dra**3*drab**2))+ &
                  Rab(3)*drr(8)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,9)=fa(1,3)*sin(alpha)*d_sin(9)- &
                  (-one/(two*dra*drab)-rra(3)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))+ &
                  Ra(3)*drr(9)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(9)/(dra**4*drab)-fr(3,3)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(3)*drr(9)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,10)=fa(1,3)*sin(alpha)*d_sin(10)- &
                  (-rra(3)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))+ &
                  Ra(3)*drr(10)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(10)/(dra**4*drab)-fr(4,1)/(dra**3*drab**2))+ &
                  Rab(3)*drr(10)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,11)=fa(1,3)*sin(alpha)*d_sin(11)- &
                  (-rra(3)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))+ &
                  Ra(3)*drr(11)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(11)/(dra**4*drab)-fr(4,2)/(dra**3*drab**2))+ &
                  Rab(3)*drr(11)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(3,12)=fa(1,3)*sin(alpha)*d_sin(12)- &
                  (-one/(two*dra*drab)-rra(3)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))+ &
                  Ra(3)*drr(12)/(dra**3*drab)+ &
                  Ra(3)*rr*(-three*ddra(12)/(dra**4*drab)-fr(4,3)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(3)*drr(12)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(4,1:3)=ffa(1:3,4)
             ffa(4,4)=fa(2,1)*sin(alpha)*d_sin(4)- &
                  (-one/(dra*drab)+rrb(1)*(-ddra(4)/(two*dra**2*drab)-fr(2,1)/(two*dra*drab**2))- &
                  rr/(dra**3*drab)-Ra(1)*drr(4)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(4)/(dra**4*drab)-fr(2,1)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(1)*drr(4)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(4)/(two*dra**2*drab**3)-three*fr(2,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,5)=fa(2,1)*sin(alpha)*d_sin(5)- &
                  (rrb(1)*(-ddra(5)/(two*dra**2*drab)-fr(2,2)/(two*dra*drab**2))- &
                  Ra(1)*drr(5)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(5)/(dra**4*drab)-fr(2,2)/(dra**3*drab**2))+ &
                  Rab(1)*drr(5)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(5)/(two*dra**2*drab**3)-three*fr(2,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,6)=fa(2,1)*sin(alpha)*d_sin(6)- &
                  (rrb(1)*(-ddra(6)/(two*dra**2*drab)-fr(2,3)/(two*dra*drab**2))- &
                  Ra(1)*drr(6)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(6)/(dra**4*drab)-fr(2,3)/(dra**3*drab**2))+ &
                  Rab(1)*drr(6)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(6)/(two*dra**2*drab**3)-three*fr(2,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,7)=fa(2,1)*sin(alpha)*d_sin(7)- &
                  (one/(two*dra*drab)+rrb(1)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))- &
                  Ra(1)*drr(7)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(7)/(dra**4*drab)-fr(3,1)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(1)*drr(7)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,8)=fa(2,1)*sin(alpha)*d_sin(8)- &
                  (rrb(1)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))- &
                  Ra(1)*drr(8)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(8)/(dra**4*drab)-fr(3,2)/(dra**3*drab**2))+ &
                  Rab(1)*drr(8)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,9)=fa(2,1)*sin(alpha)*d_sin(9)- &
                  (rrb(1)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))- &
                  Ra(1)*drr(9)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(9)/(dra**4*drab)-fr(3,3)/(dra**3*drab**2))+ &
                  Rab(1)*drr(9)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,10)=fa(2,1)*sin(alpha)*d_sin(10)- &
                  (one/(two*dra*drab)+rrb(1)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))- &
                  Ra(1)*drr(10)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(10)/(dra**4*drab)-fr(4,1)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(1)*drr(10)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,11)=fa(2,1)*sin(alpha)*d_sin(11)- &
                  (rrb(1)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))- &
                  Ra(1)*drr(11)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(11)/(dra**4*drab)-fr(4,2)/(dra**3*drab**2))+ &
                  Rab(1)*drr(11)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(4,12)=fa(2,1)*sin(alpha)*d_sin(12)- &
                  (rrb(1)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))- &
                  Ra(1)*drr(12)/(dra**3*drab)- &
                  Ra(1)*rr*(-three*ddra(12)/(dra**4*drab)-fr(4,3)/(dra**3*drab**2))+ &
                  Rab(1)*drr(12)/(two*dra*drab**3)+ &
                  Rab(1)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(5,1:3)=ffa(1:3,5)
             ffa(5,4)=fa(2,2)*sin(alpha)*d_sin(4)- &
                  (rrb(2)*(-ddra(4)/(two*dra**2*drab)-fr(2,1)/(two*dra*drab**2))- &
                  Ra(2)*drr(4)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(4)/(dra**4*drab)-fr(2,1)/(dra**3*drab**2))+ &
                  Rab(2)*drr(4)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(4)/(two*dra**2*drab**3)-three*fr(2,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,5)=fa(2,2)*sin(alpha)*d_sin(5)- &
                  (-one/(dra*drab)+rrb(2)*(-ddra(5)/(two*dra**2*drab)-fr(2,2)/(two*dra*drab**2))- &
                  rr/(dra**3*drab)-Ra(2)*drr(5)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(5)/(dra**4*drab)-fr(2,2)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(2)*drr(5)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(5)/(two*dra**2*drab**3)-three*fr(2,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,6)=fa(2,2)*sin(alpha)*d_sin(6)- &
                  (rrb(2)*(-ddra(6)/(two*dra**2*drab)-fr(2,3)/(two*dra*drab**2))- &
                  Ra(2)*drr(6)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(6)/(dra**4*drab)-fr(2,3)/(dra**3*drab**2))+ &
                  Rab(2)*drr(6)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(6)/(two*dra**2*drab**3)-three*fr(2,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,7)=fa(2,2)*sin(alpha)*d_sin(7)- &
                  (rrb(2)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))- &
                  Ra(2)*drr(7)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(7)/(dra**4*drab)-fr(3,1)/(dra**3*drab**2))+ &
                  Rab(2)*drr(7)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,8)=fa(2,2)*sin(alpha)*d_sin(8)- &
                  (one/(two*dra*drab)+rrb(2)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))- &
                  Ra(2)*drr(8)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(8)/(dra**4*drab)-fr(3,2)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(2)*drr(8)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,9)=fa(2,2)*sin(alpha)*d_sin(9)- &
                  (rrb(2)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))- &
                  Ra(2)*drr(9)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(9)/(dra**4*drab)-fr(3,3)/(dra**3*drab**2))+ &
                  Rab(2)*drr(9)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,10)=fa(2,2)*sin(alpha)*d_sin(10)- &
                  (rrb(2)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))- &
                  Ra(2)*drr(10)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(10)/(dra**4*drab)-fr(4,1)/(dra**3*drab**2))+ &
                  Rab(2)*drr(10)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,11)=fa(2,2)*sin(alpha)*d_sin(11)- &
                  (one/(two*dra*drab)+rrb(2)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))- &
                  Ra(2)*drr(11)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(11)/(dra**4*drab)-fr(4,2)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(2)*drr(11)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(5,12)=fa(2,2)*sin(alpha)*d_sin(12)- &
                  (rrb(2)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))- &
                  Ra(2)*drr(12)/(dra**3*drab)- &
                  Ra(2)*rr*(-three*ddra(12)/(dra**4*drab)-fr(4,3)/(dra**3*drab**2))+ &
                  Rab(2)*drr(12)/(two*dra*drab**3)+ &
                  Rab(2)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(6,1:3)=ffa(1:3,6)
             ffa(6,4)=fa(2,3)*sin(alpha)*d_sin(4)- &
                  (rrb(3)*(-ddra(4)/(two*dra**2*drab)-fr(2,1)/(two*dra*drab**2))- &
                  Ra(3)*drr(4)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(4)/(dra**4*drab)-fr(2,1)/(dra**3*drab**2))+ &
                  Rab(3)*drr(4)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(4)/(two*dra**2*drab**3)-three*fr(2,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,5)=fa(2,3)*sin(alpha)*d_sin(5)- &
                  (rrb(3)*(-ddra(5)/(two*dra**2*drab)-fr(2,2)/(two*dra*drab**2))- &
                  Ra(3)*drr(5)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(5)/(dra**4*drab)-fr(2,2)/(dra**3*drab**2))+ &
                  Rab(3)*drr(5)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(5)/(two*dra**2*drab**3)-three*fr(2,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,6)=fa(2,3)*sin(alpha)*d_sin(6)- &
                  (-one/(dra*drab)+rrb(3)*(-ddra(6)/(two*dra**2*drab)-fr(2,3)/(two*dra*drab**2))- &
                  rr/(dra**3*drab)-Ra(3)*drr(6)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(6)/(dra**4*drab)-fr(2,3)/(dra**3*drab**2))- &
                  rr/(four*dra*drab**3)+Rab(3)*drr(6)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(6)/(two*dra**2*drab**3)-three*fr(2,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,7)=fa(2,3)*sin(alpha)*d_sin(7)- &
                  (rrb(3)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))- &
                  Ra(3)*drr(7)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(7)/(dra**4*drab)-fr(3,1)/(dra**3*drab**2))+ &
                  Rab(3)*drr(7)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,8)=fa(2,3)*sin(alpha)*d_sin(8)- &
                  (rrb(3)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))- &
                  Ra(3)*drr(8)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(8)/(dra**4*drab)-fr(3,2)/(dra**3*drab**2))+ &
                  Rab(3)*drr(8)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,9)=fa(2,3)*sin(alpha)*d_sin(9)- &
                  (one/(two*dra*drab)+rrb(3)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))- &
                  Ra(3)*drr(9)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(9)/(dra**4*drab)-fr(3,3)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(3)*drr(9)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,10)=fa(2,3)*sin(alpha)*d_sin(10)- &
                  (rrb(3)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))- &
                  Ra(3)*drr(10)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(10)/(dra**4*drab)-fr(4,1)/(dra**3*drab**2))+ &
                  Rab(3)*drr(10)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,11)=fa(2,3)*sin(alpha)*d_sin(11)- &
                  (rrb(3)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))- &
                  Ra(3)*drr(11)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(11)/(dra**4*drab)-fr(4,2)/(dra**3*drab**2))+ &
                  Rab(3)*drr(11)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(6,12)=fa(2,3)*sin(alpha)*d_sin(12)- &
                  (one/(two*dra*drab)+rrb(3)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))- &
                  Ra(3)*drr(12)/(dra**3*drab)- &
                  Ra(3)*rr*(-three*ddra(12)/(dra**4*drab)-fr(4,3)/(dra**3*drab**2))+ &
                  rr/(four*dra*drab**3)+Rab(3)*drr(12)/(two*dra*drab**3)+ &
                  Rab(3)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(7,1:6)=ffa(1:6,7)
             ffa(7,7)=fa(3,1)*sin(alpha)*d_sin(7)- &
                  (Ra(1)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))- &
                  rr/(four*dra*drab**3)-Rab(1)*drr(7)/(two*dra*drab**3)- &
                  Rab(1)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(7,8)=fa(3,1)*sin(alpha)*d_sin(8)- &
                  (Ra(1)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))- &
                  Rab(1)*drr(8)/(two*dra*drab**3)- &
                  Rab(1)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(7,9)=fa(3,1)*sin(alpha)*d_sin(9)- &
                  (Ra(1)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))- &
                  Rab(1)*drr(9)/(two*dra*drab**3)- &
                  Rab(1)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(7,10)=fa(3,1)*sin(alpha)*d_sin(10)- &
                  (Ra(1)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))- &
                  rr/(four*dra*drab**3)-Rab(1)*drr(10)/(two*dra*drab**3)- &
                  Rab(1)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(7,11)=fa(3,1)*sin(alpha)*d_sin(11)- &
                  (Ra(1)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))- &
                  Rab(1)*drr(11)/(two*dra*drab**3)- &
                  Rab(1)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(7,12)=fa(3,1)*sin(alpha)*d_sin(12)- &
                  (Ra(1)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))- &
                  Rab(1)*drr(12)/(two*dra*drab**3)- &
                  Rab(1)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(8,1:6)=ffa(1:6,8)
             ffa(8,7)=fa(3,2)*sin(alpha)*d_sin(7)- &
                  (Ra(2)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))- &
                  Rab(2)*drr(7)/(two*dra*drab**3)- &
                  Rab(2)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(8,8)=fa(3,2)*sin(alpha)*d_sin(8)- &
                  (Ra(2)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))- &
                  rr/(four*dra*drab**3)-Rab(2)*drr(8)/(two*dra*drab**3)- &
                  Rab(2)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(8,9)=fa(3,2)*sin(alpha)*d_sin(9)- &
                  (Ra(2)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))- &
                  Rab(2)*drr(9)/(two*dra*drab**3)- &
                  Rab(2)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(8,10)=fa(3,2)*sin(alpha)*d_sin(10)- &
                  (Ra(2)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))- &
                  Rab(2)*drr(10)/(two*dra*drab**3)- &
                  Rab(2)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(8,11)=fa(3,2)*sin(alpha)*d_sin(11)- &
                  (Ra(2)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))- &
                  rr/(four*dra*drab**3)-Rab(2)*drr(11)/(two*dra*drab**3)- &
                  Rab(2)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(8,12)=fa(3,2)*sin(alpha)*d_sin(12)- &
                  (Ra(2)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))- &
                  Rab(2)*drr(12)/(two*dra*drab**3)- &
                  Rab(2)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(9,1:6)=ffa(1:6,9)
             ffa(9,7)=fa(3,3)*sin(alpha)*d_sin(7)- &
                  (Ra(3)*(-ddra(7)/(two*dra**2*drab)-fr(3,1)/(two*dra*drab**2))- &
                  Rab(3)*drr(7)/(two*dra*drab**3)- &
                  Rab(3)*rr*(-ddra(7)/(two*dra**2*drab**3)-three*fr(3,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(9,8)=fa(3,3)*sin(alpha)*d_sin(8)- &
                  (Ra(3)*(-ddra(8)/(two*dra**2*drab)-fr(3,2)/(two*dra*drab**2))- &
                  Rab(3)*drr(8)/(two*dra*drab**3)- &
                  Rab(3)*rr*(-ddra(8)/(two*dra**2*drab**3)-three*fr(3,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(9,9)=fa(3,3)*sin(alpha)*d_sin(9)- &
                  (Ra(3)*(-ddra(9)/(two*dra**2*drab)-fr(3,3)/(two*dra*drab**2))- &
                  rr/(four*dra*drab**3)-Rab(3)*drr(9)/(two*dra*drab**3)- &
                  Rab(3)*rr*(-ddra(9)/(two*dra**2*drab**3)-three*fr(3,3)/(two*dra*drab**4)))/sin(alpha)
             ffa(9,10)=fa(3,3)*sin(alpha)*d_sin(10)- &
                  (Ra(3)*(-ddra(10)/(two*dra**2*drab)-fr(4,1)/(two*dra*drab**2))- &
                  Rab(3)*drr(10)/(two*dra*drab**3)- &
                  Rab(3)*rr*(-ddra(10)/(two*dra**2*drab**3)-three*fr(4,1)/(two*dra*drab**4)))/sin(alpha)
             ffa(9,11)=fa(3,3)*sin(alpha)*d_sin(11)- &
                  (Ra(3)*(-ddra(11)/(two*dra**2*drab)-fr(4,2)/(two*dra*drab**2))- &
                  Rab(3)*drr(11)/(two*dra*drab**3)- &
                  Rab(3)*rr*(-ddra(11)/(two*dra**2*drab**3)-three*fr(4,2)/(two*dra*drab**4)))/sin(alpha)
             ffa(9,12)=fa(3,3)*sin(alpha)*d_sin(12)- &
                  (Ra(3)*(-ddra(12)/(two*dra**2*drab)-fr(4,3)/(two*dra*drab**2))- &
                  rr/(four*dra*drab**3)-Rab(3)*drr(12)/(two*dra*drab**3)- &
                  Rab(3)*rr*(-ddra(12)/(two*dra**2*drab**3)-three*fr(4,3)/(two*dra*drab**4)))/sin(alpha)
             !...............................
             ffa(10:12,1:12)=ffa(7:9,1:12)
             !...............................
             !...............................
             d_sin(1:3)=-fb(1,:)*cos_beta/(sin(beta)**2)
             d_sin(4:6)=-fb(2,:)*cos_beta/(sin(beta)**2)
             d_sin(7:9)=-fb(3,:)*cos_beta/(sin(beta)**2)
             d_sin(10:12)=-fb(4,:)*cos_beta/(sin(beta)**2)

             rr=dot_product(Rb,Rab)

             drr(1:3)=-Rb*half
             drr(4:6)=-Rb*half
             drr(7:9)=-Rab+Rb*half
             drr(10:12)=Rab+Rb*half

             rra=two*r3-r2-r1
             rrb=two*r4-r2-r1

             ffb(1,1)=fb(1,1)*sin(beta)*d_sin(1)- &
                  (-Rb(1)*(-ddrb(1)/(two*drb**2*drab)-fr(1,1)/(two*drb*drab**2))- &
                  rr/(four*drb*drab**3)+Rab(1)*drr(1)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(1)/(two*drb**2*drab**3)-three*fr(1,1)/(two*drb*drab**4)))/sin(beta)
             ffb(1,2)=fb(1,1)*sin(beta)*d_sin(2)- &
                  (-Rb(1)*(-ddrb(2)/(two*drb**2*drab)-fr(1,2)/(two*drb*drab**2))+ &
                  Rab(1)*drr(2)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(2)/(two*drb**2*drab**3)-three*fr(1,2)/(two*drb*drab**4)))/sin(beta)
             ffb(1,3)=fb(1,1)*sin(beta)*d_sin(3)- &
                  (-Rb(1)*(-ddrb(3)/(two*drb**2*drab)-fr(1,3)/(two*drb*drab**2))+ &
                  Rab(1)*drr(3)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(3)/(two*drb**2*drab**3)-three*fr(1,3)/(two*drb*drab**4)))/sin(beta)
             ffb(1,4)=fb(1,1)*sin(beta)*d_sin(4)- &
                  (-Rb(1)*(-ddrb(4)/(two*drb**2*drab)-fr(2,1)/(two*drb*drab**2))- &
                  rr/(four*drb*drab**3)+Rab(1)*drr(4)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(4)/(two*drb**2*drab**3)-three*fr(2,1)/(two*drb*drab**4)))/sin(beta)
             ffb(1,5)=fb(1,1)*sin(beta)*d_sin(5)- &
                  (-Rb(1)*(-ddrb(5)/(two*drb**2*drab)-fr(2,2)/(two*drb*drab**2))+ &
                  Rab(1)*drr(5)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(5)/(two*drb**2*drab**3)-three*fr(2,2)/(two*drb*drab**4)))/sin(beta)
             ffb(1,6)=fb(1,1)*sin(beta)*d_sin(6)- &
                  (-Rb(1)*(-ddrb(6)/(two*drb**2*drab)-fr(2,3)/(two*drb*drab**2))+ &
                  Rab(1)*drr(6)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(6)/(two*drb**2*drab**3)-three*fr(2,3)/(two*drb*drab**4)))/sin(beta)
             ffb(1,7)=fb(1,1)*sin(beta)*d_sin(7)- &
                  (one/(two*drb*drab)-Rb(1)*(-ddrb(7)/(two*drb**2*drab)-fr(3,1)/(two*drb*drab**2))+ &
                  rr/(four*drb*drab**3)+Rab(1)*drr(7)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(7)/(two*drb**2*drab**3)-three*fr(3,1)/(two*drb*drab**4)))/sin(beta)
             ffb(1,8)=fb(1,1)*sin(beta)*d_sin(8)- &
                  (-Rb(1)*(-ddrb(8)/(two*drb**2*drab)-fr(3,2)/(two*drb*drab**2))+ &
                  Rab(1)*drr(8)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(8)/(two*drb**2*drab**3)-three*fr(3,2)/(two*drb*drab**4)))/sin(beta)
             ffb(1,9)=fb(1,1)*sin(beta)*d_sin(9)- &
                  (-Rb(1)*(-ddrb(9)/(two*drb**2*drab)-fr(3,3)/(two*drb*drab**2))+ &
                  Rab(1)*drr(9)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(9)/(two*drb**2*drab**3)-three*fr(3,3)/(two*drb*drab**4)))/sin(beta)
             ffb(1,10)=fb(1,1)*sin(beta)*d_sin(10)- &
                  (-one/(two*drb*drab)-Rb(1)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))+ &
                  rr/(four*drb*drab**3)+Rab(1)*drr(10)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(1,11)=fb(1,1)*sin(beta)*d_sin(11)- &
                  (-Rb(1)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))+ &
                  Rab(1)*drr(11)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(1,12)=fb(1,1)*sin(beta)*d_sin(12)- &
                  (-Rb(1)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))+ &
                  Rab(1)*drr(12)/(two*drb*drab**3)+ &
                  Rab(1)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(2,1)=fb(1,2)*sin(beta)*d_sin(1)- &
                  (-Rb(2)*(-ddrb(1)/(two*drb**2*drab)-fr(1,1)/(two*drb*drab**2))+ &
                  Rab(2)*drr(1)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(1)/(two*drb**2*drab**3)-three*fr(1,1)/(two*drb*drab**4)))/sin(beta)
             ffb(2,2)=fb(1,2)*sin(beta)*d_sin(2)- &
                  (-Rb(2)*(-ddrb(2)/(two*drb**2*drab)-fr(1,2)/(two*drb*drab**2))- &
                  rr/(four*drb*drab**3)+Rab(2)*drr(2)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(2)/(two*drb**2*drab**3)-three*fr(1,2)/(two*drb*drab**4)))/sin(beta)
             ffb(2,3)=fb(1,2)*sin(beta)*d_sin(3)- &
                  (-Rb(2)*(-ddrb(3)/(two*drb**2*drab)-fr(1,3)/(two*drb*drab**2))+ &
                  Rab(2)*drr(3)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(3)/(two*drb**2*drab**3)-three*fr(1,3)/(two*drb*drab**4)))/sin(beta)
             ffb(2,4)=fb(1,2)*sin(beta)*d_sin(4)- &
                  (-Rb(2)*(-ddrb(4)/(two*drb**2*drab)-fr(2,1)/(two*drb*drab**2))+ &
                  Rab(2)*drr(4)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(4)/(two*drb**2*drab**3)-three*fr(2,1)/(two*drb*drab**4)))/sin(beta)
             ffb(2,5)=fb(1,2)*sin(beta)*d_sin(5)- &
                  (-Rb(2)*(-ddrb(5)/(two*drb**2*drab)-fr(2,2)/(two*drb*drab**2))- &
                  rr/(four*drb*drab**3)+Rab(2)*drr(5)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(5)/(two*drb**2*drab**3)-three*fr(2,2)/(two*drb*drab**4)))/sin(beta)
             ffb(2,6)=fb(1,2)*sin(beta)*d_sin(6)- &
                  (-Rb(2)*(-ddrb(6)/(two*drb**2*drab)-fr(2,3)/(two*drb*drab**2))+ &
                  Rab(2)*drr(6)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(6)/(two*drb**2*drab**3)-three*fr(2,3)/(two*drb*drab**4)))/sin(beta)
             ffb(2,7)=fb(1,2)*sin(beta)*d_sin(7)- &
                  (-Rb(2)*(-ddrb(7)/(two*drb**2*drab)-fr(3,1)/(two*drb*drab**2))+ &
                  Rab(2)*drr(7)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(7)/(two*drb**2*drab**3)-three*fr(3,1)/(two*drb*drab**4)))/sin(beta)
             ffb(2,8)=fb(1,2)*sin(beta)*d_sin(8)- &
                  (one/(two*drb*drab)-Rb(2)*(-ddrb(8)/(two*drb**2*drab)-fr(3,2)/(two*drb*drab**2))+ &
                  rr/(four*drb*drab**3)+Rab(2)*drr(8)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(8)/(two*drb**2*drab**3)-three*fr(3,2)/(two*drb*drab**4)))/sin(beta)
             ffb(2,9)=fb(1,2)*sin(beta)*d_sin(9)- &
                  (-Rb(2)*(-ddrb(9)/(two*drb**2*drab)-fr(3,3)/(two*drb*drab**2))+ &
                  Rab(2)*drr(9)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(9)/(two*drb**2*drab**3)-three*fr(3,3)/(two*drb*drab**4)))/sin(beta)
             ffb(2,10)=fb(1,2)*sin(beta)*d_sin(10)- &
                  (-Rb(2)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))+ &
                  Rab(2)*drr(10)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(2,11)=fb(1,2)*sin(beta)*d_sin(11)- &
                  (-one/(two*drb*drab)-Rb(2)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))+ &
                  rr/(four*drb*drab**3)+Rab(2)*drr(11)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(2,12)=fb(1,2)*sin(beta)*d_sin(12)- &
                  (-Rb(2)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))+ &
                  Rab(2)*drr(12)/(two*drb*drab**3)+ &
                  Rab(2)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(3,1)=fb(1,3)*sin(beta)*d_sin(1)- &
                  (-Rb(3)*(-ddrb(1)/(two*drb**2*drab)-fr(1,1)/(two*drb*drab**2))+ &
                  Rab(3)*drr(1)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(1)/(two*drb**2*drab**3)-three*fr(1,1)/(two*drb*drab**4)))/sin(beta)
             ffb(3,2)=fb(1,3)*sin(beta)*d_sin(2)- &
                  (-Rb(3)*(-ddrb(2)/(two*drb**2*drab)-fr(1,2)/(two*drb*drab**2))+ &
                  Rab(3)*drr(2)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(2)/(two*drb**2*drab**3)-three*fr(1,2)/(two*drb*drab**4)))/sin(beta)
             ffb(3,3)=fb(1,3)*sin(beta)*d_sin(3)- &
                  (-Rb(3)*(-ddrb(3)/(two*drb**2*drab)-fr(1,3)/(two*drb*drab**2))- &
                  rr/(four*drb*drab**3)+Rab(3)*drr(3)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(3)/(two*drb**2*drab**3)-three*fr(1,3)/(two*drb*drab**4)))/sin(beta)
             ffb(3,4)=fb(1,3)*sin(beta)*d_sin(4)- &
                  (-Rb(3)*(-ddrb(4)/(two*drb**2*drab)-fr(2,1)/(two*drb*drab**2))+ &
                  Rab(3)*drr(4)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(4)/(two*drb**2*drab**3)-three*fr(2,1)/(two*drb*drab**4)))/sin(beta)
             ffb(3,5)=fb(1,3)*sin(beta)*d_sin(5)- &
                  (-Rb(3)*(-ddrb(5)/(two*drb**2*drab)-fr(2,2)/(two*drb*drab**2))+ &
                  Rab(3)*drr(5)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(5)/(two*drb**2*drab**3)-three*fr(2,2)/(two*drb*drab**4)))/sin(beta)
             ffb(3,6)=fb(1,3)*sin(beta)*d_sin(6)- &
                  (-Rb(3)*(-ddrb(6)/(two*drb**2*drab)-fr(2,3)/(two*drb*drab**2))- &
                  rr/(four*drb*drab**3)+Rab(3)*drr(6)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(6)/(two*drb**2*drab**3)-three*fr(2,3)/(two*drb*drab**4)))/sin(beta)
             ffb(3,7)=fb(1,3)*sin(beta)*d_sin(7)- &
                  (-Rb(3)*(-ddrb(7)/(two*drb**2*drab)-fr(3,1)/(two*drb*drab**2))+ &
                  Rab(3)*drr(7)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(7)/(two*drb**2*drab**3)-three*fr(3,1)/(two*drb*drab**4)))/sin(beta)
             ffb(3,8)=fb(1,3)*sin(beta)*d_sin(8)- &
                  (-Rb(3)*(-ddrb(8)/(two*drb**2*drab)-fr(3,2)/(two*drb*drab**2))+ &
                  Rab(3)*drr(8)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(8)/(two*drb**2*drab**3)-three*fr(3,2)/(two*drb*drab**4)))/sin(beta)
             ffb(3,9)=fb(1,3)*sin(beta)*d_sin(9)- &
                  (one/(two*drb*drab)-Rb(3)*(-ddrb(9)/(two*drb**2*drab)-fr(3,3)/(two*drb*drab**2))+ &
                  rr/(four*drb*drab**3)+Rab(3)*drr(9)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(9)/(two*drb**2*drab**3)-three*fr(3,3)/(two*drb*drab**4)))/sin(beta)
             ffb(3,10)=fb(1,3)*sin(beta)*d_sin(10)- &
                  (-Rb(3)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))+ &
                  Rab(3)*drr(10)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(3,11)=fb(1,3)*sin(beta)*d_sin(11)- &
                  (-Rb(3)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))+ &
                  Rab(3)*drr(11)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(3,12)=fb(1,3)*sin(beta)*d_sin(12)- &
                  (-one/(two*drb*drab)-Rb(3)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))+ &
                  rr/(four*drb*drab**3)+Rab(3)*drr(12)/(two*drb*drab**3)+ &
                  Rab(3)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(4:6,:)=ffb(1:3,:)
             !...............................
             ffb(7,1:6)=ffb(1:6,7)
             ffb(7,7)=fb(3,1)*sin(beta)*d_sin(7)- &
                  (-one/(drb*drab)-rra(1)*(-ddrb(7)/(two*drb**2*drab)-fr(3,1)/(two*drb*drab**2))- &
                  rr/(drb**3*drab)+Rb(1)*drr(7)/(drb**3*drab)+ &
                  Rb(1)*rr*(-three*ddrb(7)/(drb**4*drab)-fr(3,1)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(1)*drr(7)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(7)/(two*drb**2*drab**3)-three*fr(3,1)/(two*drb*drab**4)))/sin(beta)
             ffb(7,8)=fb(3,1)*sin(beta)*d_sin(8)- &
                  (-rra(1)*(-ddrb(8)/(two*drb**2*drab)-fr(3,2)/(two*drb*drab**2))+ &
                  Rb(1)*drr(8)/(drb**3*drab)+ &
                  Rb(1)*rr*(-three*ddrb(8)/(drb**4*drab)-fr(3,2)/(drb**3*drab**2))- &
                  Rab(1)*drr(8)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(8)/(two*drb**2*drab**3)-three*fr(3,2)/(two*drb*drab**4)))/sin(beta)
             ffb(7,9)=fb(3,1)*sin(beta)*d_sin(9)- &
                  (-rra(1)*(-ddrb(9)/(two*drb**2*drab)-fr(3,3)/(two*drb*drab**2))+ &
                  Rb(1)*drr(9)/(drb**3*drab)+ &
                  Rb(1)*rr*(-three*ddrb(9)/(drb**4*drab)-fr(3,3)/(drb**3*drab**2))- &
                  Rab(1)*drr(9)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(9)/(two*drb**2*drab**3)-three*fr(3,3)/(two*drb*drab**4)))/sin(beta)
             ffb(7,10)=fb(3,1)*sin(beta)*d_sin(10)- &
                  (-rra(1)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))+ &
                  rr/(drb**3*drab)+Rb(1)*drr(10)/(drb**3*drab)+ &
                  Rb(1)*rr*(-three*ddrb(10)/(drb**4*drab)-fr(4,1)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(1)*drr(10)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(7,11)=fb(3,1)*sin(beta)*d_sin(11)- &
                  (-rra(1)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))+ &
                  Rb(1)*drr(11)/(drb**3*drab)+ &
                  Rb(1)*rr*(-three*ddrb(11)/(drb**4*drab)-fr(4,2)/(drb**3*drab**2))- &
                  Rab(1)*drr(11)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(7,12)=fb(3,1)*sin(beta)*d_sin(12)- &
                  (-rra(1)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))+ &
                  Rb(1)*drr(12)/(drb**3*drab)+ &
                  Rb(1)*rr*(-three*ddrb(12)/(drb**4*drab)-fr(4,3)/(drb**3*drab**2))- &
                  Rab(1)*drr(12)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(8,1:6)=ffb(1:6,8)
             ffb(8,7)=fb(3,2)*sin(beta)*d_sin(7)- &
                  (-rra(2)*(-ddrb(7)/(two*drb**2*drab)-fr(3,1)/(two*drb*drab**2))+ &
                  Rb(2)*drr(7)/(drb**3*drab)+ &
                  Rb(2)*rr*(-three*ddrb(7)/(drb**4*drab)-fr(3,1)/(drb**3*drab**2))- &
                  Rab(2)*drr(7)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(7)/(two*drb**2*drab**3)-three*fr(3,1)/(two*drb*drab**4)))/sin(beta)
             ffb(8,8)=fb(3,2)*sin(beta)*d_sin(8)- &
                  (-one/(drb*drab)-rra(2)*(-ddrb(8)/(two*drb**2*drab)-fr(3,2)/(two*drb*drab**2))- &
                  rr/(drb**3*drab)+Rb(2)*drr(8)/(drb**3*drab)+ &
                  Rb(2)*rr*(-three*ddrb(8)/(drb**4*drab)-fr(3,2)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(2)*drr(8)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(8)/(two*drb**2*drab**3)-three*fr(3,2)/(two*drb*drab**4)))/sin(beta)
             ffb(8,9)=fb(3,2)*sin(beta)*d_sin(9)- & 
                  (-rra(2)*(-ddrb(9)/(two*drb**2*drab)-fr(3,3)/(two*drb*drab**2))+ &
                  Rb(2)*drr(9)/(drb**3*drab)+ &
                  Rb(2)*rr*(-three*ddrb(9)/(drb**4*drab)-fr(3,3)/(drb**3*drab**2))- &
                  Rab(2)*drr(9)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(9)/(two*drb**2*drab**3)-three*fr(3,3)/(two*drb*drab**4)))/sin(beta)
             ffb(8,10)=fb(3,2)*sin(beta)*d_sin(10)- &
                  (-rra(2)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))+ &
                  Rb(2)*drr(10)/(drb**3*drab)+ &
                  Rb(2)*rr*(-three*ddrb(10)/(drb**4*drab)-fr(4,1)/(drb**3*drab**2))- &
                  Rab(2)*drr(10)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(8,11)=fb(3,2)*sin(beta)*d_sin(11)- &
                  (-rra(2)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))+ &
                  rr/(drb**3*drab)+Rb(2)*drr(11)/(drb**3*drab)+ &
                  Rb(2)*rr*(-three*ddrb(11)/(drb**4*drab)-fr(4,2)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(2)*drr(11)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(8,12)=fb(3,2)*sin(beta)*d_sin(12)- &
                  (-rra(2)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))+ &
                  Rb(2)*drr(12)/(drb**3*drab)+ &
                  Rb(2)*rr*(-three*ddrb(12)/(drb**4*drab)-fr(4,3)/(drb**3*drab**2))- &
                  Rab(2)*drr(12)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(9,1:6)=ffb(1:6,9)
             ffb(9,7)=fb(3,3)*sin(beta)*d_sin(7)- &
                  (-rra(3)*(-ddrb(7)/(two*drb**2*drab)-fr(3,1)/(two*drb*drab**2))+ &
                  Rb(3)*drr(7)/(drb**3*drab)+ &
                  Rb(3)*rr*(-three*ddrb(7)/(drb**4*drab)-fr(3,1)/(drb**3*drab**2))- &
                  Rab(3)*drr(7)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(7)/(two*drb**2*drab**3)-three*fr(3,1)/(two*drb*drab**4)))/sin(beta)
             ffb(9,8)=fb(3,3)*sin(beta)*d_sin(8)- &
                  (-rra(3)*(-ddrb(8)/(two*drb**2*drab)-fr(3,2)/(two*drb*drab**2))+ &
                  Rb(3)*drr(8)/(drb**3*drab)+ &
                  Rb(3)*rr*(-three*ddrb(8)/(drb**4*drab)-fr(3,2)/(drb**3*drab**2))- &
                  Rab(3)*drr(8)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(8)/(two*drb**2*drab**3)-three*fr(3,2)/(two*drb*drab**4)))/sin(beta)
             ffb(9,9)=fb(3,3)*sin(beta)*d_sin(9)- &
                  (-one/(drb*drab)-rra(3)*(-ddrb(9)/(two*drb**2*drab)-fr(3,3)/(two*drb*drab**2))- &
                  rr/(drb**3*drab)+Rb(3)*drr(9)/(drb**3*drab)+ &
                  Rb(3)*rr*(-three*ddrb(9)/(drb**4*drab)-fr(3,3)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(3)*drr(9)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(9)/(two*drb**2*drab**3)-three*fr(3,3)/(two*drb*drab**4)))/sin(beta)
             ffb(9,10)=fb(3,3)*sin(beta)*d_sin(10)- &
                  (-rra(3)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))+ &
                  Rb(3)*drr(10)/(drb**3*drab)+ &
                  Rb(3)*rr*(-three*ddrb(10)/(drb**4*drab)-fr(4,1)/(drb**3*drab**2))- &
                  Rab(3)*drr(10)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(9,11)=fb(3,3)*sin(beta)*d_sin(11)- &
                  (-rra(3)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))+ &
                  Rb(3)*drr(11)/(drb**3*drab)+ &
                  Rb(3)*rr*(-three*ddrb(11)/(drb**4*drab)-fr(4,2)/(drb**3*drab**2))- &
                  Rab(3)*drr(11)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(9,12)=fb(3,3)*sin(beta)*d_sin(12)- &
                  (-rra(3)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))+ &
                  rr/(drb**3*drab)+Rb(3)*drr(12)/(drb**3*drab)+ &
                  Rb(3)*rr*(-three*ddrb(12)/(drb**4*drab)-fr(4,3)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(3)*drr(12)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(10,1:9)=ffb(1:9,10)
             ffb(10,10)=fb(4,1)*sin(beta)*d_sin(10)- &
                  (one/(drb*drab)+rrb(1)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))- &
                  rr/(drb**3*drab)-Rb(1)*drr(10)/(drb**3*drab)- &
                  Rb(1)*rr*(-three*ddrb(10)/(drb**4*drab)-fr(4,1)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(1)*drr(10)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(10,11)=fb(4,1)*sin(beta)*d_sin(11)- &
                  (rrb(1)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))- &
                  Rb(1)*drr(11)/(drb**3*drab)- &
                  Rb(1)*rr*(-three*ddrb(11)/(drb**4*drab)-fr(4,2)/(drb**3*drab**2))- &
                  Rab(1)*drr(11)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(10,12)=fb(4,1)*sin(beta)*d_sin(12)- &
                  (rrb(1)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))- &
                  Rb(1)*drr(12)/(drb**3*drab)- &
                  Rb(1)*rr*(-three*ddrb(12)/(drb**4*drab)-fr(4,3)/(drb**3*drab**2))- &
                  Rab(1)*drr(12)/(two*drb*drab**3)- &
                  Rab(1)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(11,1:9)=ffb(1:9,11)
             ffb(11,10)=fb(4,2)*sin(beta)*d_sin(10)- &
                  (rrb(2)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))- &
                  Rb(2)*drr(10)/(drb**3*drab)- &
                  Rb(2)*rr*(-three*ddrb(10)/(drb**4*drab)-fr(4,1)/(drb**3*drab**2))- &
                  Rab(2)*drr(10)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(11,11)=fb(4,2)*sin(beta)*d_sin(11)- &
                  (one/(drb*drab)+rrb(2)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))- &
                  rr/(drb**3*drab)-Rb(2)*drr(11)/(drb**3*drab)- &
                  Rb(2)*rr*(-three*ddrb(11)/(drb**4*drab)-fr(4,2)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(2)*drr(11)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(11,12)=fb(4,2)*sin(beta)*d_sin(12)- &
                  (rrb(2)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))- &
                  Rb(2)*drr(12)/(drb**3*drab)- &
                  Rb(2)*rr*(-three*ddrb(12)/(drb**4*drab)-fr(4,3)/(drb**3*drab**2))- &
                  Rab(2)*drr(12)/(two*drb*drab**3)- &
                  Rab(2)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             ffb(12,1:9)=ffb(1:9,12)
             ffb(12,10)=fb(4,3)*sin(beta)*d_sin(10)- &
                  (rrb(3)*(-ddrb(10)/(two*drb**2*drab)-fr(4,1)/(two*drb*drab**2))- &
                  Rb(3)*drr(10)/(drb**3*drab)- &
                  Rb(3)*rr*(-three*ddrb(10)/(drb**4*drab)-fr(4,1)/(drb**3*drab**2))- &
                  Rab(3)*drr(10)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(10)/(two*drb**2*drab**3)-three*fr(4,1)/(two*drb*drab**4)))/sin(beta)
             ffb(12,11)=fb(4,3)*sin(beta)*d_sin(11)- &
                  (rrb(3)*(-ddrb(11)/(two*drb**2*drab)-fr(4,2)/(two*drb*drab**2))- &
                  Rb(3)*drr(11)/(drb**3*drab)- &
                  Rb(3)*rr*(-three*ddrb(11)/(drb**4*drab)-fr(4,2)/(drb**3*drab**2))- &
                  Rab(3)*drr(11)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(11)/(two*drb**2*drab**3)-three*fr(4,2)/(two*drb*drab**4)))/sin(beta)
             ffb(12,12)=fb(4,3)*sin(beta)*d_sin(12)- &
                  (one/(drb*drab)+rrb(3)*(-ddrb(12)/(two*drb**2*drab)-fr(4,3)/(two*drb*drab**2))- &
                  rr/(drb**3*drab)-Rb(3)*drr(12)/(drb**3*drab)- &
                  Rb(3)*rr*(-three*ddrb(12)/(drb**4*drab)-fr(4,3)/(drb**3*drab**2))- &
                  rr/(four*drb*drab**3)-Rab(3)*drr(12)/(two*drb*drab**3)- &
                  Rab(3)*rr*(-ddrb(12)/(two*drb**2*drab**3)-three*fr(4,3)/(two*drb*drab**4)))/sin(beta)
             !...............................
             !...............................
             k1=3*(ia1-1)
             k2=3*(ia2-1)
             k3=3*(ia3-1)
             k4=3*(ia4-1)
             do l=1,3
                do m=1,12
                   if(m <= 3) then
                      fst_r=fr(1,m)
                      fst_x=fx(1,m)
                      fst_a=fa(1,m)
                      fst_b=fb(1,m)
                      m1=k1+m
                   else if(m > 3 .and. m <= 6) then 
                      m1=k2+(m-3)
                      fst_r=fr(2,m-3)
                      fst_x=fx(2,m-3)
                      fst_a=fa(2,m-3)
                      fst_b=fb(2,m-3)
                   else if(m > 6 .and. m <= 9) then
                      m1=k3+(m-6)
                      fst_r=fr(3,m-6)
                      fst_x=fx(3,m-6)
                      fst_a=fa(3,m-6)
                      fst_b=fb(3,m-6)
                   else if(m > 9 .and. m <= 12) then
                      m1=k4+(m-9)
                      fst_r=fr(4,m-9)
                      fst_x=fx(4,m-9)
                      fst_a=fa(4,m-9)
                      fst_b=fb(4,m-9)
                   end if
                   H(k1+l,m1)=H(k1+l,m1)+(d2E_dr2*fst_r+d2E_drdx*fst_x+d2E_drda*fst_a+d2E_drdb*fst_b)*fr(1,l)+ &
                        (d2E_dxdr*fst_r+d2E_dx2*fst_x+d2E_dxda*fst_a+d2E_dxdb*fst_b)*fx(1,l)+ &
                        (d2E_dadr*fst_r+d2E_dadx*fst_x+d2E_da2*fst_a+d2E_dadb*fst_b)*fa(1,l)+ &
                        (d2E_dbdr*fst_r+d2E_dbdx*fst_x+d2E_dbda*fst_a+d2E_db2*fst_b)*fb(1,l)+ &
                        dE_dr*ffr(l,m)+dE_dx*ffx(l,m)+dE_da*ffa(l,m)+dE_db*ffb(l,m)
                   H(k2+l,m1)=H(k2+l,m1)+(d2E_dr2*fst_r+d2E_drdx*fst_x+d2E_drda*fst_a+d2E_drdb*fst_b)*fr(2,l)+ &
                        (d2E_dxdr*fst_r+d2E_dx2*fst_x+d2E_dxda*fst_a+d2E_dxdb*fst_b)*fx(2,l)+ &
                        (d2E_dadr*fst_r+d2E_dadx*fst_x+d2E_da2*fst_a+d2E_dadb*fst_b)*fa(2,l)+ &
                        (d2E_dbdr*fst_r+d2E_dbdx*fst_x+d2E_dbda*fst_a+d2E_db2*fst_b)*fb(2,l)+ &
                        dE_dr*ffr(l+3,m)+dE_dx*ffx(l+3,m)+dE_da*ffa(l+3,m)+dE_db*ffb(l+3,m)
                   H(k3+l,m1)=H(k3+l,m1)+(d2E_dr2*fst_r+d2E_drdx*fst_x+d2E_drda*fst_a+d2E_drdb*fst_b)*fr(3,l)+ &
                        (d2E_dxdr*fst_r+d2E_dx2*fst_x+d2E_dxda*fst_a+d2E_dxdb*fst_b)*fx(3,l)+ &
                        (d2E_dadr*fst_r+d2E_dadx*fst_x+d2E_da2*fst_a+d2E_dadb*fst_b)*fa(3,l)+ &
                        (d2E_dbdr*fst_r+d2E_dbdx*fst_x+d2E_dbda*fst_a+d2E_db2*fst_b)*fb(3,l)+ &
                        dE_dr*ffr(l+6,m)+dE_dx*ffx(l+6,m)+dE_da*ffa(l+6,m)+dE_db*ffb(l+6,m)
                   H(k4+l,m1)=H(k4+l,m1)+(d2E_dr2*fst_r+d2E_drdx*fst_x+d2E_drda*fst_a+d2E_drdb*fst_b)*fr(4,l)+ &
                        (d2E_dxdr*fst_r+d2E_dx2*fst_x+d2E_dxda*fst_a+d2E_dxdb*fst_b)*fx(4,l)+ &
                        (d2E_dadr*fst_r+d2E_dadx*fst_x+d2E_da2*fst_a+d2E_dadb*fst_b)*fa(4,l)+ &
                        (d2E_dbdr*fst_r+d2E_dbdx*fst_x+d2E_dbda*fst_a+d2E_db2*fst_b)*fb(4,l)+ &
                        dE_dr*ffr(l+9,m)+dE_dx*ffx(l+9,m)+dE_da*ffa(l+9,m)+dE_db*ffb(l+9,m)
                end do
             end do

          end if

       end do l_j
    end do l_i
!!$print*,"E_total=",E_coulomb
!!$do i=1,n_species
!!$print*,i,Grad(:,i)
!!$end do

  end subroutine dipol_dipol_E_and_F
  !****************************************************************

  !****************************************************************
  subroutine shutdown_electrostatic()
    
    integer(i4_kind) :: status

    if(allocated(dip_pair)) then
       deallocate(dip_pair,stat=status)
       if(status /= 0) call error_handler("MolMech: failed DIP_PAIR deallocation")
    end if

    if(allocated(N_interact)) then
       deallocate(N_interact,stat=status)
       if(status /= 0) call error_handler("MolMech: failed N_INTERACT deallocation")
    end if

  end subroutine shutdown_electrostatic
  !****************************************************************

  !****************************************************************
end module coulomb_module
