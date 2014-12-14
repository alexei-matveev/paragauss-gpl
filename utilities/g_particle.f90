program g_particle
  implicit none
  !Calculate free energy of free particle (atom, ion and so on)
  integer, parameter :: NELEM=103
  real*8 :: mass(NELEM)
  data mass/ &
         1.007825d0,  4.002603d0, 7.016005d0,   9.012182d0, &
        11.009305d0, 12.000000d0,14.003074d0 , 15.994915d0, &
        18.998403d0, 19.992439d0,22.989770d0 , 23.985045d0, &
        26.981541d0, 27.976928d0,30.973763d0 , 31.973909d0, &
        34.968853d0, 39.962383d0,38.963708d0 , 39.962591d0, &
        44.955914d0, 47.947947d0,50.943963d0 , 51.940510d0, &
        54.938046d0, 55.934939d0,58.933198d0 , 57.935347d0, &
        62.929599d0, 63.929145d0,68.925581d0 , 73.921179d0, &
        74.921596d0, 79.916521d0,78.918336d0 , 83.911506d0, &
        84.911800d0, 87.956250d0,88.905856d0 , 89.904708d0, &
        92.906378d0, 97.905405d0,97.907110d0 ,101.904348d0, &
       102.90550d0 ,106.903480d0,106.905095d0,113.903361d0, &
       114.90388d0 ,119.902199d0,120.903824d0,129.906230d0, &
       126.904477d0,131.90415d0 ,132.905770d0,137.905240d0, &
       138.906360d0,139.90544d0 ,140.907660d0,141.907730d0, &
       144.912691d0,151.91974d0 ,152.921240d0,157.924110d0, &
       158.92535d0 ,163.929180d0,164.930330d0,167.930310d0, &
       168.93423d0 ,173.938870d0,174.940790d0,177.943710d0, &
       180.94801d0 ,183.950950d0,186.955770d0,191.961490d0, &
       192.96294d0 ,194.964790d0,196.966560d0,201.970630d0, &
       204.97441d0 ,207.976640d0,208.980390d0,208.982420d0, &
       219.01130d0 ,222.017570d0,223.019734d0,226.025406d0, &
       227.02775d0 ,232.038050d0,231.035881d0,238.050786d0, &
       237.048169d0,244.06420d0 ,243.0614d0  ,247.0703d0,   &
       247.0703d0  ,251.0796d0  ,252.0829d0  ,257.0951d0,   &
       258.0986d0  ,259.1009d0  ,262.1100d0/
  real*8 :: total_mass
  real*8, parameter :: TWO_PI=6.2831853071795864769352867665588d0
  real*8, parameter :: h_Planck=6.62606896d-34
  real*8, parameter :: THREE=3.0d0
  real*8, parameter :: k_B=1.3806504d-23
!  real*8, parameter :: tmpT=298.15d0
  real*8, parameter :: FIVE_HALF=2.5d0
  real*8, parameter :: P=100000.0d0
  real*8, parameter :: THREE_HALF=1.5d0
  real*8, parameter :: R_gas=8.314472d0
  real*8, parameter :: E_h2Jmol=2.6254995d+06
  real*8, parameter :: u_Atom_Mass=1.660538782d-27
  real*8, parameter :: THOUSAND=1000.0d0
  real*8 :: q_t, E_t, S_t, U, S, entH, Gfe, Gfe_tot, tot_energy,tmpT
  integer :: atom_number

  write(6,*) 'Atomic number, electronic energy (au) and temperature(K)'
  read(5,*) atom_number, tot_energy, tmpT

  !*** Calculations of translational contributions !************************************
  total_mass=mass(atom_number)*u_Atom_Mass
  tot_energy=tot_energy*E_h2Jmol

  q_t = (sqrt(TWO_PI*total_mass)/h_Planck)**THREE*(k_B*tmpT)**FIVE_HALF/P
  
  E_t = THREE_HALF*R_gas*tmpT
  
  S_t = R_gas*(log(q_t)+FIVE_HALF)

  U = E_t
  S = S_t
  entH = U + R_gas*tmpT
  Gfe = entH - tmpT*S

  Gfe_tot = Gfe + tot_energy

  write(6,*)
  write(6,01)
01 format('  +---------------------------------------------------------------------------------------+')
  write(6,02)
02 format('  |',87X,'|')
  write(6,03)
03 format('  |     Thermodynamic properties section:',49X,'|')
  write(6,02)
  write(6,04) tmpT, P
04 format('  |',5X,'T   = ',F8.3,' K',21X,'p = ',F9.1,' Pa',29X,'|')
  write(6,02)
  write(6,01)
  write(6,02)
  write(6,05) total_mass/u_Atom_Mass
05 format('  |     Translational Contributions',10X,'total mass: ',F9.3,'  au',20X'|')
  write(6,06) E_t/THOUSAND, S_t/THOUSAND
06 format('  |',5X,'E_t  =  ',E12.6,'  kJ/mol',9X,'S_t  =  ',E12.6,'  kJ/K*mol',15X,'|')
  write(6,07) E_t/E_h2Jmol, S_t/E_h2Jmol
07 format('  |',12X,'(',E12.6,'      au)',15X,'(',E12.6,'      au/K)',14X,'|')
  write(6,01)
  write(6,02)
  write(6,18)
18 format('  |     Summary: THERMODYNAMICAL PROPERTIES',47X,'|')
  write(6,02)
  write(6,19) U/THOUSAND, S/THOUSAND
19 format('  |',5X,'U    =  ',E12.6,'  kJ/mol',9X,'S    =  ',E12.6,'  kJ/K*mol',15X,'|')
  write(6,07) U/E_h2Jmol, S/E_h2Jmol
  write(6,02)
  write(6,20) entH/THOUSAND, Gfe/THOUSAND
20 format('  |',5X,'H    =  ',E12.6,'  kJ/mol',9X,'G    =  ',E12.6,'  kJ/mol',17X,'|')
  write(6,07) entH/E_h2Jmol, Gfe/E_h2Jmol
  write(6,02)
  write(6,01)
  write(6,02)
  write(6,21)
21 format('  |     End of thermodynamic properties section',43X,'|')
  write(6,02)
  write(6,22)  Gfe_tot/THOUSAND
22 format('  |     G_tot = ',E17.10,' kJ/mol',13X,'(including "e_sum")',18X,'|')
  write(6,23)  Gfe_tot/E_h2Jmol
23 format('  |',12X,'(',E17.10,'     au)',49X,'|')
  write(6,02)
2201 format('  |',13X,E17.10,' kJ/mol',50X,'|')
  write(6,01)
  write(6,*)

end program g_particle
