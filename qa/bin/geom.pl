#!/usr/bin/perl
#
# $Id: geom.pl,v 1.1 2001/07/11 13:30:29 matveev Exp matveev $
# $Log: geom.pl,v $
# Revision 1.1  2001/07/11  13:30:29  matveev
# Initial revision
#
#
# Extracts geometry part from output
#

$txt = "";
while( <> ){
	if( /Geometry in au/../^ -------/ ){
		# print "+$_";
		$txt .= $_;
		# print;
	}else{
		# print "-$_";
		if( $txt ){
			print $txt;
			$toprint = $txt;
			$txt = "";
		}
	}
}

#print $toprint;
