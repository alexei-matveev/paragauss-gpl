#!/usr/bin/perl

$D=10;
if($ARGV[0] =~ /-(\d+)/){
	$D=$1;
	shift(@ARGV);
}
$DD =$D+1;

#warn("## all numbers will be cut off to $D digits\n");

while(<>){
#	while(/\d+\.\d{$DD,}/o){
#		s/(\d+)\.(\d{$D})(\d+)/\1.\2/;
#	}
	s/(\d+\.\d{$D})[0-9]*/\1/og;
	print;
}
