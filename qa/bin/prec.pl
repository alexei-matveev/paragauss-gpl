#!/usr/bin/perl

$D=10;
if($ARGV[0] =~ /-(\d+)/){
	$D=$1;
	shift(@ARGV);
}
$DD =$D+1;
#$pow = 10**$D;

#warn("## all numbers will be cut off to $D digits\n");

while(<>){
#	while(/(\d+\.\d{$DD,})/o){
#		$r = int( $1 * $pow + 0.5 )/$pow;
#		s/\d+\.\d{$DD,}/$r/;
#	}
        # cut figures beyond D after decimal point:
	s/(\d+\.\d{$D})[0-9]*/\1/og;
        # nomralize -0.00000 to 0.00000:
	s/-(0.0+)\b/ \1/og;
	print;
}
