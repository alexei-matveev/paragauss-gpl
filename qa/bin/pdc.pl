#!/usr/bin/perl
#
# Extracts Potential-Derived Charges part from output
#
#  Unique number    Atom       Charge
# -------------------------------------
#        1            O       -0.79542
#        2            H        0.39771
#        2            H        0.39771
# -------------------------------------
# ESP rrms:  0.12565
# Dipole moment of PDC:    0.0000   0.0000   0.9147

while( <> ){
	if( /^  Unique number    Atom       Charge/../^ Dipole moment of PDC:/ ){
		print;
	}
}
