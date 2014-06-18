  &UNIQUE_ATOM_BASISSET
    LMAX_OB     = 1,
    LMAX_CH     = 2,
    LMAX_XC     = 0,
    ZC          = 2,
    LMAX_PSEUDO = 3
  /UNIQUE_ATOM_BASISSET
  # ORBITAL BASIS:
# F s ECP2MWB  : 4 1 1.3
#  (4s5p)/[2s3p]-Basissatz fuer PP. von Ref 17.            
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 4,
    N_UNCONTRACTED_FCTS = 1,
    N_CONTRACTED_FCTS = 1
  /UNIQUE_ATOM_BASIS
  # exponents:
  0.3700810 %
  1.2141230 9.4142770 51.6427630
  # contractions:
  0.0  %
  0.5898350 -0.1530090 0.0085660  

# F p ECP2MWB  : 5 1 1.3
#  (4s5p)/[2s3p]-Basissatz fuer PP. von Ref 17.         
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 5,
    N_UNCONTRACTED_FCTS = 2,
    N_CONTRACTED_FCTS = 1
  /UNIQUE_ATOM_BASIS
  # exponents:
  0.0847720 0.3465300 %
  1.3423110 4.9545730 22.3008300
  # contractions:
  0.0 0.0  %
  0.5076660 0.2372870 0.0518510  

  # PSEUDOPOTENTIAL:
# F  ECP  ECP2MWB :  2 3 0 22
#  Q=7., MEFIT, WB, Ref 17; CPP: alpha=0.0016;delta=6.5859;ncut=1.
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
  22.350400 11.175200
  # coefficients
  102.597952 19.049663

  # TYPE: L = 1
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 2 /
  # powers of the prefactors
  2 2
  # exponents
  26.476800 13.238400
  # coefficients
  -15.143960 2.802921

  # TYPE: L = 2
  &UNIQUE_ATOM_PSEUDOPOT N_EXPONENTS = 1 /
  # powers of the prefactors
  2
  # exponents
  0.031600
  # coefficients
  -0.001863

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
  0.169544 0.69306 2.684622 9.909146 44.60166

  # CHARHE FIT:  s-type, generated from orbitals s-basis:
  &UNIQUE_ATOM_BASIS
    N_EXPONENTS = 4,
    N_UNCONTRACTED_FCTS = 4,
    N_CONTRACTED_FCTS = 0
  /UNIQUE_ATOM_BASIS
  0.740162 2.428246 18.828554 103.285526

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

