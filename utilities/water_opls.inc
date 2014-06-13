 &task
   job_task="gradients"
 /task
 
 &main_options
   energy_unit="kcal/mol"
 /main_options
 
 &nb_lists_options
   automatic_nb_lists=true
 /nb_lists_options
 
 # OPLS FORCE FIELD PARAMETERS FOR NON BONDED INTERACTIONS
 # BETWEEN QM SOLUTE SUBSYSTEM AND MM SOLVENT (TIP4P WATER)
 #
 # ATOMIC TYPES AND NAMES
 # ---------------------------------------------------------
 # TYPE          NAME FOR MOLMECH       CHARGE       OPLS
 #                  INPUT FILE                       TYPES
 # ---------------------------------------------------------
 # H(TIP4P)              H               0.52
 # O(TIP4P)              O               0.00
 # Q(TIP4P)              OQ             -1.04
 #
 # H(WATER)              HW              0.00(0.52)
 # O(WATER)              OW              0.00(-1.04)
 #
 # H(H3O+)               H+              0.524
 # O(H3O+)               O+             -0.571
 #
 # U(UO2_2+)             U               2.50
 # O(UO2_2+)             OU             -0.25
 #
 # H(CH3,CH2,CH)         HC              0.00
 # CH2(C-CH2-C)          CH2             0.00          9
 # CH3(CH3-C)            CH3C            0.00          7
 # C(HCOO-)              COO             0.70          17
 #  (RCOO-)
 #  (HCOOH)                              0.50          1
 #  (RCOOH)
 # CH3(H3C-O)            CH3O            0.25          63
 # O(-COO-)              OO             -0.80          18
 #  (-C=O)                              -0.45          2
 # O(RCOOH)              OH             -0.58          ?
 #  (RCOOCH3)                           -0.40          62
 # H(RCOOH)              HO              0.435         24
 # C(C6H6)               CA              0.00          11
 #  (C6H6R)
 #
 # INTERNAL PARAMETERS FOR WATER HAS TO BE USED TO EXCLUDE
 # INTERACTION BETWEEN ATOMS WITHIN MOLECULE AND KEEP INTERNAL
 # ENERGY EQUAL TO ZERO.
 # ONLY VAN DER WAALS INTERACTIONS BETWEEN QM AND MM PARTS ARE
 # USED FOR QM/MM CALCULATIONS.
 # REAL MM CHARGES OF QM ATOMS ARE ZERO. BUT IF YOU WANT TO
 # OPTIMIZE BY MM TIP4P WATERS AROUND FIXED QM CLUSTER QM
 # ATOMS HAS TO CARY DEFINITE CHARGES WHICH CAN BE TAKEN
 # FROM HERE.
 #
 # JORGENCEN, W., TIRADO-RIVES, J. JACS 1988, V110,P1657
 # GUILBAUD, P., WIPFF, G. J MOL STRUC (THEOCHEM) 1996, V366, P55
 # SINGH, U. C., BROWN, F. K., BASH, P. A., KOLMANN, P. A.
 # JACS 1987, V109, P1607
 # JORGENSEN, W. L., CHANDRASEKHAR, J., MADURA, J. D.,
 # IMPEY, R. W., KLEIN. M. L. J CHEM PHYS 1983, V79, P926
 
 &species
   name="H"
   main_name="H_"
   charge= 0.52
 /species
 
 &species
   name="O"
   main_name="O_"
   charge= 0.00
 /species
 
 &species
   name="OQ"
   main_name="O_"
   charge=-1.04
 /species
 
 &species
   name="HW"
   main_name="H_"
   charge= 0.00
 /species
 
 &species
   name="OW"
   main_name="O_"
   charge= 0.00
 /species
 
 &species
   name="H+"
   main_name="H_"
   charge= 0.00
 /species
 
 &species
   name="O+"
   main_name="O_"
   charge= 0.00
 /species
 
 &species
   name="U"
   main_name="U_"
   charge= 0.00
 /species
 
 &species
   name="OU"
   main_name="O_"
   charge= 0.00
 /species
 
 &species
   name="HC"
   main_name="H_"
   charge= 0.00
 /species
 
 &species
   name="CH2"
   main_name="C_"
   charge= 0.00
 /species
 
 &species
   name="CH3C"
   main_name="C_"
   charge= 0.00
 /species
 
 &species
   name="COO"
   main_name="C_"
   charge= 0.00
 /species
 
 &species
   name="CH3O"
   main_name="C_"
   charge= 0.00
 /species
 
 &species
   name="OO"
   main_name="O_"
   charge= 0.00
 /species
 
 &species
   name="OH"
   main_name="O_"
   charge= 0.00
 /species
 
 &species
   name="HO"
   main_name="H_"
   charge= 0.00
 /species
 
 &species
   name="CA"
   main_name="C_"
   charge= 0.00
 /species
 
 &potential
   pot_name="HARM_STR"
   atom_name="H","O"
   ff_parameter=0.0, 0.9572
 /potential
 
 &potential
   pot_name="HARM_STR"
   atom_name="H","OQ"
   ff_parameter=0.0, 0.9572
 /potential
 
 &potential
   pot_name="HARM_BND"
   atom_name="H","O","H"
   ff_parameter=0.0, 104.52
 /potential
 
 &potential
   pot_name="HARM_BND"
   atom_name="H","OQ","H"
   ff_parameter=0.0, 104.52
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","O"
   ff_parameter= 600003.2, 609.9
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","OW"
   ff_parameter= 600003.2, 609.9
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","O+"
   ff_parameter= 463442.6, 531.7
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","U"
   ff_parameter= 487786.0, 697.0
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","OU"
   ff_parameter= 636795.6, 669.7
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","CH2"
   ff_parameter=1887009.2,1010.3
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","CH3C"
   ff_parameter=2214253.5,1181.0
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","COO"
   ff_parameter=1396011.7, 844.0
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","CH3O"
   ff_parameter=1923237.1,1117.5
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","OO"
   ff_parameter= 477492.7, 587.0
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","OH"
   ff_parameter= 465649.1, 549.9
 /potential
 
 &potential
   pot_name="L_J"
   atom_name="O","CA"
   ff_parameter=1428863.5, 863.9
 /potential
 
