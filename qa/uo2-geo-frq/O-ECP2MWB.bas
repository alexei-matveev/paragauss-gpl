  &UNIQUE_ATOM_BASISSET
    LMAX_OB     = 1,
    LMAX_CH     = 2,
    LMAX_XC     = 0,
    ZC          = 2,
    LMAX_PSEUDO = 3
  /UNIQUE_ATOM_BASISSET
  # ORBITAL BASIS:
# O  s ECP2MWB : 4 1 1.3
#  (4s5p)/[2s3p]-Basissatz fuer PP. von Ref 17.             
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 4,
    N_UNCONTRACTED_FCTS = 1,
    N_CONTRACTED_FCTS = 1
  /UNIQUE_ATOM_BASIS
  # exponents:
  0.296070 %
  0.976483 5.911346 47.105518
  # contractions:
  0.0  %
  -0.563118 0.129568 -0.014408  

# O  p ECP2MWB : 5 1 1.3
#  (4s5p)/[2s3p]-Basissatz fuer PP. von Ref 17.               
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 5,
    N_UNCONTRACTED_FCTS = 2,
    N_CONTRACTED_FCTS = 1
  /UNIQUE_ATOM_BASIS
  # exponents:
  0.070200 0.284189 %
  1.078253 3.900702 16.692219
  # contractions:
  0.0 0.0  %
  0.500188 0.222613 0.044856  

  # PSEUDOPOTENTIAL:
# O  ECP ECP2MWB : 2 3 0 16
#  Q=6., MEFIT, WB, Ref 17.   
  # TYPE: local
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 1 /
  # powers of the prefactors
  2
  # exponents
  1.00000000
  # coefficients
  0.00000000

  # TYPE: L = 0
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 1 /
  # powers of the prefactors
  2
  # exponents
  10.44567000
  # coefficients
  50.77106900

  # TYPE: L = 1
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 1 /
  # powers of the prefactors
  2
  # exponents
  18.04517400
  # coefficients
  -4.90355100

  # TYPE: L = 2
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 1 /
  # powers of the prefactors
  2
  # exponents
  8.16479800
  # coefficients
  -3.31212400

  # CORE DENSITY: r2- and s-fakeed:
  &UNIQUE_ATOM_CORE_DENSITY N_EXPONENTS = 1 /
  -1.0
   0.0
  &UNIQUE_ATOM_CORE_DENSITY N_EXPONENTS = 1 /
  -1.0
   0.0

  # CHARHE FIT: r2-type, generated from orbitals p-basis:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 5,
    N_UNCONTRACTED_FCTS = 5,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  0.1404 0.568378 2.156506 7.801404 33.384438

  # CHARHE FIT:  s-type, generated from orbitals s-basis:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 4,
    N_UNCONTRACTED_FCTS = 4,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  0.59214 1.952966 11.822692 94.211036

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

