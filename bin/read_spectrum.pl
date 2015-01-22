#
# read_spectrum.pl: exctracted from
# zust.pl to be used in dos.pl
# $Id: read_spectrum.pl,v 1.2 2008/11/06 12:25:16 matveev Exp matveev $
# $Log: read_spectrum.pl,v $
# Revision 1.2  2008/11/06 12:25:16  matveev
# fix reading newer spectra, non-rel and srel have preceedence
#
# Revision 1.1  2002/03/20 08:49:55  matveev
# Initial revision
#
#

$verbose = 0;

$au2ev = 27.211658;

sub read_spectrum{
	local($polarized,*occupation) = @_;

	local($n,$irr,$occ,$e_ev,$e_au);
	print "zust: reading spectrum\n" if $verbose;
	LEVEL:
	while(<>){
                print "-$_" if $verbose;
		next LEVEL if /^\s*#/;
		last LEVEL if /^\s*-{30,}/i;
		#
		# Unpolarized spectrum:
		#
		#^   3     1   AU      2.0000  -18074.290317
                #
		# Polarized spectrum: shifted to the right ny spaces for "beta" electrons
                #
                #^ 377    58   T2G   5  1    3.0000         -1.082981
                #^                                              379    58   T2G   5  2    3.0000         -1.077196
                #
		if( $polarized ){
			if(/^( ){40,}/){
				$updn = "dn";
			}else{
				$updn = "up";
			}
		}else{
				$updn = "";
		}
		if(/^\s*\d+\s+(\d+)\s+(\S+)\s+(\d+(\s+[12])?)?\s+(\S+)\s+(\S+)/){
                    #   indx    n      irr                        occ     eny
			$n    = $1; # principal quantum number
			$irr  = $2; # irrep symbol
			$ixs  = $3; # indices of irrep and, eventually spin (1|2)
                        $iab  = $4; # one or two with preceeding spaces
			$occ  = $5; # occupation
			$e_ev = $6; # energy in eV
			$e_au = $e_ev/$au2ev;

			if( ! $updn ){
				&pushlevel(*occupation,$n,$irr,$occ,$e_au);
			}else{
				&pushlevel(*occupation,$n,$irr,$occ,$e_au,$updn);
			}
                        print "(SR) $n,$irr,$updn ( $e_ev eV ) == $occ\n" if $verbose;
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
                        print "(SO) $n,$irr ( $e_ev eV ) == $occ\n" if $verbose;
			next LEVEL;
		}
	}
	print "zust: stop reading spectrum\n" if $verbose;
}

sub pushlevel{
	local(*occupation,$n,$irr,$occ,$e_au,$updn) = @_;

	local($level) = "$n,$irr";
	$level =~ y/a-z/A-Z/;
	if( $updn ){ $level .= ",$updn"; }

#       print "$level $occ $e_ev $e_au\n" if $verbose;

	#die("read_spectrum: energy(1) already defined for $level\n")
	#	if defined $energy{$level};

	# $energy{$level}     = $e_au;
	$occupation{$level} = $occ;

	# $state{$irr}   += $occ;
	# $energy_total  += $occ * $e_au;
	# if ( $updn ){
	# 	$pol_state{$irr,$updn} +=  $occ;
	# }

	if( ! defined $degeneracy{$irr} ){
		$degeneracy{$irr} = $occ;
	}else{
		if( $occ > $degeneracy{$irr} ){
			$degeneracy{$irr} = $occ;
		}
	}
}

1;
