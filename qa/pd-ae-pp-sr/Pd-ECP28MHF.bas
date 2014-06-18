  &UNIQUE_ATOM_BASISSET
    LMAX_OB     = 2,
    LMAX_CH     = 2,
    LMAX_XC     = 0,
    ZC          = 28,
    LMAX_PSEUDO = 4
  /UNIQUE_ATOM_BASISSET
  # ORBITAL BASIS:
# Pd s ECP28MHF : 8 1 1.3
#  (8s7p6d)/[6s5p3d]-Basissatz zu PP. von D. Andrae  
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 8,
    N_UNCONTRACTED_FCTS = 5,
    N_CONTRACTED_FCTS = 1
  /UNIQUE_ATOM_BASIS
  # exponents:
  0.0160000 0.0493600 0.1233030 0.6239210 1.4063570 %
  3.1821100 7.1657170 8.4756400
  # contractions:
  0.0 0.0 0.0 0.0 0.0  %
  0.4441355 2.2938092 -1.7692795  

# Pd p ECP28MHF : 7 2 1.2 3.4
#  (8s7p6d)/[6s5p3d]-Basissatz zu PP. von D. Andrae  
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 7,
    N_UNCONTRACTED_FCTS = 3,
    N_CONTRACTED_FCTS = 2
  /UNIQUE_ATOM_BASIS
  # exponents:
  0.0247540 0.0791750 0.3726570 %
  0.6396240 1.2161850 3.3925940 4.2460970
  # contractions:
  0.0 0.0 0.0 0.0 0.0  %
  5.6533768 -5.4567269  
  0.0 0.0 0.0  %
  0.2619366 0.7615573 %
  0.0 0.0 

# Pd d ECP28MHF : 6 1 1.4
#  (8s7p6d)/[6s5p3d]-Basissatz zu PP. von D. Andrae  
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 6,
    N_UNCONTRACTED_FCTS = 2,
    N_CONTRACTED_FCTS = 1
  /UNIQUE_ATOM_BASIS
  # exponents:
  0.0600000 0.1815250 %
  0.4622960 1.1096690 2.5612630 7.1696120
  # contractions:
  0.0 0.0  %
  0.3870177 0.4891212 0.2808141 -0.0216768  

  # PSEUDOPOTENTIAL:
# Pd ECP ECP28MHF : 28 4 0 32
#  Q=18., MEFIT, HF, Dirk Andrae, Diplomarbeit 1989.
  # TYPE: local
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 1 /
  # powers of the prefactors
  2
  # exponents
  1.000000
  # coefficients
  0.000000

  # TYPE: L = 0
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 2 /
  # powers of the prefactors
  2 2
  # exponents
  11.800000 5.880000
  # coefficients
  235.245962 34.682650

  # TYPE: L = 1
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 2 /
  # powers of the prefactors
  2 2
  # exponents
  10.820000 5.460000
  # coefficients
  170.430396 25.301872

  # TYPE: L = 2
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 2 /
  # powers of the prefactors
  2 2
  # exponents
  9.870000 4.500000
  # coefficients
  70.206324 14.777436

  # TYPE: L = 3
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 2 /
  # powers of the prefactors
  2 2
  # exponents
  13.070000 6.530000
  # coefficients
  -31.401554 -5.298965

  # CORE DENSITY: r2- and s-fakeed:
  &UNIQUE_ATOM_CORE_DENSITY N_EXPONENTS = 1 /
  -1.0
   0.0
  &UNIQUE_ATOM_CORE_DENSITY N_EXPONENTS = 1 /
  -1.0
   0.0

  # CHARHE FIT: r2-type, generated from orbitals p-basis:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 7,
    N_UNCONTRACTED_FCTS = 7,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  0.049508 0.15835 0.745314 1.279248 2.43237 6.785188 8.492194

  # CHARHE FIT:  s-type, generated from orbitals s-basis:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 8,
    N_UNCONTRACTED_FCTS = 8,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  0.032 0.09872 0.246606 1.247842 2.812714 6.36422 14.331434 16.95128

  # CHARHE FIT:  p-type:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 5,
    N_UNCONTRACTED_FCTS = 5,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  0.1 0.25 0.625 1.5625 3.90625

  # CHARHE FIT:  d-type:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 5,
    N_UNCONTRACTED_FCTS = 5,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  0.2 0.5 1.25 3.125 7.8125

  # XC FIT: r2- and s-fakeed:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 1,
    N_UNCONTRACTED_FCTS = 1,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  -1.0    

  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 1,
    N_UNCONTRACTED_FCTS = 1,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  -1.0    
