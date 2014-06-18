#!/usr/bin/perl
#
# Extracts frequency part from output
#

while( <> ){
	if( /SD: FREQUENCIES:/../^SD:$/ ){
		print;
	}
	if( /SD: FREQUENCIES \(NUMERICAL CARTESIAN HESSIAN\):/../^SD:$/ ){
		print;
	}

}
