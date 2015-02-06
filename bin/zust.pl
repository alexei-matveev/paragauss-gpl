#!/usr/bin/perl
#
# This  script prints the  number of  electrons occupying  orbitals of
# different symmetry and spin as derived from the output.
#
# Copyright (C) 2000-2015 Alexei Matveev
#
# This program is free software; you can redistribute it and/or modify
# it under  the terms of the  GNU General Public License  version 2 as
# published by the Free Software Foundation [1].
#
# This program is distributed in the  hope that it will be useful, but
# WITHOUT  ANY   WARRANTY;  without  even  the   implied  warranty  of
# MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE. See  the GNU
# General Public License for more details.
#
# [1] http://www.gnu.org/licenses/gpl-2.0.html
#
# $Source: /home/matveev/bin.pl/RCS/zust.pl,v $
# $Revision: 1.8 $
# $Author: matveev $
#
# Changes after last checkout:
# Print usage information.
# Ignore comments in cmdmode.
# -spectrum N takes Nth spectrum (for GeoOpt)
#
# $Log: zust.pl,v $
# Revision 1.8  2004/10/28 11:36:29  matveev
# fix the polarized spectra with new fileds
#
# Revision 1.7  2004/10/28 10:49:38  matveev
# fixes for the last NR/SR spectra (new entry)
#
# Revision 1.6  2004/10/28 10:38:00  matveev
# -spectrum selects one of many in GeoOpt,
# print usage info
#
# Revision 1.5  2001/07/09  16:31:25  matveev
# Extension for polarized (NR) spectrum: summary only, DOS not tested
#
# Revision 1.4  2001/07/09  15:52:36  matveev
# check in before polarized (NR) extension
#
#
# 20-04-2000:
# Summarize electronic state
# Original name "zust" is due to Stefan Vent and
# comes apparantly from german "Zustand" (State)
#
# Usage: zust.pl output
#        takes spin-orbit and unpolarized nonrelativistic outputs

$au2ev = 27.211658;

$verbose = 0;

$dos   = 0;
$read_cmdfile =0;

$polarized = 0;
$which_spectrum = 1;

while( $ARGV[0] =~ /^-(\w+)/ ){
    $opt = $1;
    if( $opt eq "dos" ){
	$dos = 1;
	shift(@ARGV);
    }elsif( $opt eq "width" ){
	shift(@ARGV);
	$params{"WIDTH"} = $ARGV[0];
	$dos   = 1;
	shift(@ARGV);
    }elsif( $opt eq "occ" ){
	$params{"OCCUPIED"} = 1;
	shift(@ARGV);
    }elsif( $opt eq "spectrum" ){
	shift(@ARGV);
	$which_spectrum = $ARGV[0];
	shift(@ARGV);
    }else{
	&usage();
	exit(1);
    }
}

sub usage{
	print STDERR "Usage: $0 [-dos [-occ] [-width W]] <input> [<cmdfile>]\n";
	print STDERR "          <input> is an output of ParaGauss\n";
	print STDERR "          -dos for printing Desnity of States\n";
	print STDERR "          <cmdfile> may contain DOS settings\n";
	print STDERR "          -occ, -width W, provide command-line access\n";
	print STDERR "          to some of them\n";
}

LINE:
while(<>){
	next LINE if /^\s*#/;
	if(/^ This is a open shell calculation/){
		$polarized = 1;
	}

	if(/^\s*-{14,} Spectrum -{14,}/){
		$which_spectrum--;
		print "which_spectrum = $which_spectrum\n" if $verbose;
		next LINE if $which_spectrum != 0;
		print "zust: reading spectrum\n" if $verbose;
		LEVEL:
		while(<>){
			next LEVEL if /^\s*#/;
			last LEVEL if /^\s*-{30,}/i;
			#
			# Unpolarized spectrum:
			#
			#^   3     1   AU      2.0000  -18074.290317
			# latest version:
			#^   1     1   A1G   1    2.0000    -114936.577333
			# (shifted to the right ny spaces for "beta" electrons)
			if( $polarized ){
				if(/^( ){40,}/){
					$updn = "dn";
				}else{
					$updn = "up";
				}
			}else{
					$updn = "";
			}
			if(/^\s*\d+\s+(\d+)\s+(\S+)\s+(\d+)?\s+(\d+)?\s+(\S+)\s+(\S+)/){
				$n    = $1;
				$irr  = $2;
				$new_entry1 = $3;
				$new_entry2 = $4;
				$occ  = $5;
				$e_ev = $6;
				$e_au = $e_ev/$au2ev;
				print "$n $irr ($new_entry1, $new_entry2) $occ $e_ev\n" if $verbose;

				if( ! $updn ){
					&pushlevel($n,$irr,$occ,$e_au);
				}else{
					&pushlevel($n,$irr,$occ,$e_au,$updn);
				}
				next LEVEL;
			}
			#
			# SinOrbit spectrum:
			#^   1     1   E1/2G     2.0000  -82555.170465   -3033.816928 /* 1s1/2 */
			#
			if(/^\s*\d+\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
				$n    = $1;
				$irr  = $2;
				$occ  = $3;
				$e_ev = $4;
				$e_au = $5;

				&pushlevel($n,$irr,$occ,$e_au);
				next LEVEL;
			}
		}
		print "zust: stop reading spectrum\n" if $verbose;
	}
	if(/^\s*Group\s+(\S+)\s+eigenvalues of class operator and primitive characters/){
		if(! defined $group){
			$group = $1;
			print "# $group group\n" if $verbose;
		}
		$found = 0;
		ORDER:
		while(<>){
			if(/\s*(\S+)\s+([+-]?\d*\.\d+\s+)+\s*$/){
				$irr = $1;
				$found++;
				print "$found $irr\n" if $verbose;
				$irrep_order{$irr} = $found;
			}elsif($found){
				last ORDER;
			}
		}

		# now clean array:
		foreach $irr ( keys %irrep_order ){
			if( $irr =~ /[+-]?\d*\.\d+/){
				print "$irr is not an irrep\n" if $verbose;
				delete $irrep_order{$irr};
			}
			if( $irr =~ /Irr/){
				print "$irr is not an irrep\n" if $verbose;
				delete $irrep_order{$irr};
			}
		}
	}
	if(/^_PARAMETERS_/){
	    $dos = 1;
	    $read_cmdfile = 1;
	    &cmdmode();
	    print "exited cmdmode\n" if $verbose;
	}
}
print "exited getline cycle\n" if $verbose;

if( ! $dos ){
	print "# $group group\n";
	$total = 0;
	foreach $irr ( sort{$irrep_order{$a}<=>$irrep_order{$b}} keys %state ){
		if ( ! $polarized ){
			printf " %10.4f # %-10s\n",$state{$irr},$irr;
			$total += $state{$irr};
		}else{
			printf " %10.4f %10.4f # %-10s\n",
				$pol_state{$irr,"up"},$pol_state{$irr,"dn"},$irr;
			$pol_total{"up"} += $pol_state{$irr,"up"};
			$pol_total{"dn"} += $pol_state{$irr,"dn"};
		}
	}
	if ( ! $polarized ) {
		printf "#%10.4f  %-10s\n",$total,"TOTAL";
	}else{
		printf "#%10.4f %10.4f # %-10s\n",$pol_total{"up"},$pol_total{"dn"},"TOTAL U/D";
		printf "#      %10.4f      # %-10s\n",$pol_total{"up"}+$pol_total{"dn"},"TOTAL";
	}
}

    if( $dos && ! $read_cmdfile ){
	&dosplot(*energy,*occupation,*params,*graph);
	&printgr(*graph);
    }

#
# end of commands, below only subs
#

sub pushlevel{
	local($n,$irr,$occ,$e_au,$updn) = @_;

	local($level) = "$n $irr";
	$level =~ y/a-z/A-Z/;
	if( $updn ){ $level .= " ".$updn; }

	print "$level $occ $e_ev $e_au\n" if $verbose;

	die("zust: energy(1) already defined for $level\n")
		if defined $energy{$level};

	$energy{$level}     = $e_au;
	$occupation{$level} = $occ;

	$state{$irr}   += $occ;
	$energy_total  += $occ * $e_au;
	if ( $updn ){
		$pol_state{$irr,$updn} +=  $occ;
	}

	if( ! defined $degeneracy{$irr} ){
		$degeneracy{$irr} = $occ;
	}else{
		if( $occ > $degeneracy{$irr} ){
			$degeneracy{$irr} = $occ;
		}
	}
}

sub cmdmode{
    print "sub cmdmode entered\n" if $verbose;
    local($key,$val);
    while(<>){
	s/^#.*$//;
	if(/^\s*(\w+)\s+(\S+)/){
	    $key = $1;
	    $val = $2;
	    $params{$key} = $val;
	}
	if(/^_GO_/){
	    @graph = 0;
	    &dosplot(*energy,*occupation,*params,*graph);
	    &printgr(*graph);
	    last;
	}
    }
}

sub defaults{
    print "sub defaults entered\n" if $verbose;
    local(*par) = @_;

    if( ! defined $par{"WIDTH"} ){
	$par{"WIDTH"} = 0.05;
    }
    if( ! defined $par{"BOTTOM_ENERGY"} ){
	$par{"BOTTOM_ENERGY"} = -10.0;
    }
    if( ! defined $par{"TOP_ENERGY"} ){
	$par{"TOP_ENERGY"} = 0.0;
    }
    if( ! defined $par{"NPOINTS"} ){
	$par{"NPOINTS"} = 500;
    }
    if( ! defined $par{"OCCUPIED"} ){
	$par{"OCCUPIED"} = 0;
    }
    if( ! defined $par{"ENERGY_UNITS"} ){
	$par{"ENERGY_UNITS"} = "au";
    }
    local($lev);
    if( $par{"ENERGY_UNITS"} eq "au" ){
	foreach $lev ( keys %energy ){
	    $energy{$lev} = $energy{$lev} * $au2ev;
	}
	$par{"ENERGY_UNITS"} = "eV";
    }
    &checkdegen(*degeneracy);
}

sub checkdegen{
    local(*deg) = @_;
    local(%assume);

    $assume{"A"}=1;
    $assume{"B"}=1;
    $assume{"E"}=2;
    $assume{"F"}=3; # is there any in LCGTO ?
    $assume{"T"}=3;
    $assume{"G"}=4;
    $assume{"H"}=5;
    $assume{"I"}=6;

    $assume{"1E"}=1;
    $assume{"2E"}=1;

    local($pair) = 1.0;

    foreach $irr ( keys %deg ){
	if($verbose){
	    print "found degeneracy of $irr is $deg{$irr}\n";
	}
	if( $irr =~ m|^[EGI]\d/2| ){
	    $pair = 1.0;
	}
	if( $irr =~ m|[12]E\d/2| ){
	    $pair = 0.5;
	}
    }


    local($irr,$checked);
    foreach $irr ( keys %deg ){
	$checked = 0;
	foreach $lett ( keys %assume ){
	    if( $irr =~ /^$lett/ ){
		$checked = 1;
		if( $deg{$irr} != $assume{$lett}*$pair ){
		    warn "Check degeneracy of $irr! Is it really $deg{$irr}?\n";
		}
	    }
	}
	if( ! $checked ){
	    warn "Dont know irrep $irr\n";
	}
    }
}


sub printgr{
   print "sub printgr entered\n" if $verbose;
   local(*g) = @_;

   local($x);
   foreach $x ( sort{$a<=>$b} keys %g ){
       printf "%16.8f\t%16.8f\n",$x, $g{$x};
   }
}

sub printpar{
   print "sub printpar entered\n" if $verbose;
   local(*g) = @_;

   local($x);
   foreach $x ( sort keys %g ){
       printf "%-16s\t%-16s\n",$x, $g{$x};
   }
}

sub dosplot{
    print "sub dosplot entered\n" if $verbose;
    local(*eny,*occ,*par,*gr) = @_;
    local($i);

    &defaults(*par);
    &printpar(*par) if $verbose;

    local($width) = $par{"WIDTH"};
    local($bot)   = $par{"BOTTOM_ENERGY"};
    local($top)   = $par{"TOP_ENERGY"};
    local($e_bot) = $bot - 3*$width;
    local($e_top) = $top + 3*$width;
    local($np)    = $par{"NPOINTS"};

    local($occ_weight) = $par{"OCCUPIED"};

    local($dE) = ($top - $bot) / $np;

    local(@Xs,@Ys);
    @Xs=0; @Ys=0;

    for($i=0; $i<$np; $i++){
	$Xs[$i] = $bot + ( $i + 0.5 ) * $dE;
    }

    local($lev,$e,$o,$x,$weight);
    foreach $lev ( keys %eny ){
	$e = $eny{$lev};
	if( ($e > $e_bot) && ($e < $e_top)){
	    $o = $occ{$lev};
	    if($occ_weight){
		$weight = $o;
	    }else{
		$lev =~ /^\d+\s+(\S+)/;
		$irr = $1;
		$weight = $degeneracy{$irr};
	    }
	    print ">> $lev $e $o\n" if $verbose;
	    for($i=0; $i<$np; $i++){
		$x = ( $Xs[$i] - $e ) / $width;
		$Ys[$i] += $weight * &invexp2($x);
	    }
	}
    }
    for($i=0; $i<$np; $i++){
	$gr{$Xs[$i]} = $Ys[$i];
    }
}

sub invexp2{
    # exp( - x^2 )
    local($x) = @_;
    $x *= $x;
    ($x>50.0)?(0.0):(exp(-$x));
}




