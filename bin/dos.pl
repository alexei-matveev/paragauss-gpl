#!/usr/bin/perl
#
# dos.pl: utility to plot DOS using POAN
# of ParaGauss.
# $Id: dos.pl,v 1.6 2015/01/22 13:26:17 matveev Exp matveev $
# $Log: dos.pl,v $
# Revision 1.6  2015/01/22 13:26:17  matveev
# print debug only if verbose
#
# Revision 1.5  2010/11/05 15:08:43  matveev
# state as of 2008-11-12
#
# Revision 1.4  2008/11/06 12:39:27  matveev
# barplot added, accept command line arguments
#
# Revision 1.3  2008/11/05 17:45:27  matveev
# compute and print as a coment the weighted center of the dos plot
#
# Revision 1.2  2008/08/07 12:25:26  matveev
# print matching levels sorted by energy
#
# Revision 1.1  2002/03/20 08:48:48  matveev
# Initial revision
#
#
#

#usage:
#usage: Usage:
#usage:
#usage: dos.pl [options] output [cmdfile] > plot.xy
#usage:
#usage:	  output  -- output of ParaGauss with POPAN
#usage:	  cmdfile -- parameters file (see below),
#usage:	             comments ignored
#usage:	  graph will be the output graph in XY format
#usage:
#usage: Options:
#usage:
#usage:	--pattern  PATTERN
#usage:	--bottom   BOTTOM_ENERGY
#usage:	--top      TOM_ENERGY
#usage:	--occupied (0|1)
#usage:	--output   (gauss|bars)
#usage:
#usage:	--skip N   if you want to ignore N first
#usage:            geometry iterations
#usage:
#usage:
#usage:# This is an example parameter file,
#usage:# missing parameters default to specified values
#usage:# (except for PATTERN)
#usage:# It is enough to specify two lines to make it work:
#usage:#1:_PARAMETERS_
#usage:#2:_GO_
#usage:#
#usage:_PARAMETERS_
#usage:# this line MUST be present! See also the end of file!
#usage:
#usage:WIDTH   0.05
#usage:# width of the gaussian
#usage:
#usage:NPOINTS 500
#usage:# how smooth should be the graph
#usage:
#usage:BOTTOM_ENERGY -10.0
#usage:TOP_ENERGY      0.0
#usage:
#usage:OCCUPIED 0
#usage:# set it to one if you want only
#usage:# occupied levels.
#usage:
#usage:PATTERN  [SPD]<3>
#usage:# Perl-style pattern of what contributions should be
#usage:# counted, in this case S-,P-,D-contributions of 3rd unique
#usage:# atom will be printed. Deafult is ".*" i.e. everything matches.
#usage:# 
#usage:# In general, levels are indexed by complex structure:
#usage:# Spin restricted case:
#usage:#	N,IRR,L<A>
#usage:#  e.g. 2,A1G,P<5> is P-shell contribution of unique atom 5 to 2nd orbital of A1G
#usage:# Polarized case:
#usage:#	N,IRR,UD,L<A>
#usage:# where UD = "up" or "dn" for alpha and beta electrons
#usage:
#usage:_GO_
#usage:# this line MUST be present!
#usage:
#usage:
#usage:
#usage:

push(@INC,"/home/matveev/bin.pl");

require("read_spectrum.pl");

$verbose = 0;
$Verbose = 0;
$VErbose = 0;

print "@ARGV\n" if $verbose;

# output gauss-broadened plot:
$output = "gauss";

# show significant features in DOS comments
# that have magnitude higher than:
$feature = 0.1;

while( $ARGV[0] =~ /-.*/ ){
  if( $ARGV[0] eq "--skip" || $ARGV[0] eq "-s" ){
          shift(@ARGV);
          $skip_blocks = $ARGV[0];
          shift(@ARGV);
  }elsif( $ARGV[0] eq "--pattern" || $ARGV[0] eq "-p" ){
          shift(@ARGV);
          $params{'PATTERN'} = $ARGV[0];
          shift(@ARGV);
  }elsif( $ARGV[0] eq "--bottom" || $ARGV[0] eq "-b" ){
          shift(@ARGV);
          $params{'BOTTOM_ENERGY'} = $ARGV[0];
          shift(@ARGV);
  }elsif( $ARGV[0] eq "--top" || $ARGV[0] eq "-t" ){
          shift(@ARGV);
          $params{'TOP_ENERGY'} = $ARGV[0];
          shift(@ARGV);
  }elsif( $ARGV[0] eq "--occupied" || $ARGV[0] eq "-o" ){
          shift(@ARGV);
          $params{'OCCUPIED'} = $ARGV[0];
          shift(@ARGV);
  }elsif( $ARGV[0] eq "--output" ){
          shift(@ARGV);
          $output = $ARGV[0];
          shift(@ARGV);
  }else{
    &usage("no such option $ARGV[0]");
  }
}

$graph_center   = 0.0;
$graph_integral = 0.0;

$polarized = 0;

while(<>){
	$line_count++;
	print "-main: $_" if $Verbose;
	if( ! $go_ahead ){
	if(/^ This is a open shell calculation/){
		$polarized = 1;
		$go_ahead = 1;
	}
	}
	next if ! $process_rest && ! /POPULATION/ ;
	
	if(/\*\*\* POPULATION ANALYSIS\*\*\*/){
#		print "-found($line_count): $_";
		if( $skip_blocks > 0 ){
			$skip_blocks--;
			print "skip_blocks reduced to $skip_blocks\n" if $verbose;
			next;
		}
		$process_rest = 1;
		&read_popan();
	}
	if(/^\s*-{14,} Spectrum -{14,}/){
#		print "now read_spectrum (line $line_count)\n";
		&read_spectrum($polarized,*occupation);
	}
	if(/_PARAMETERS_/){
#		print "now cmdmode (line $line_count)\n";
		&cmdmode();
		last;
	}
}

# now filter and output the data:
&go();

sub go{
  local(%graph,%bars);
  @graph = 0;
  &defaults(*params); # set to defaults the *unset* params
  &barplot(*energy,*occupation,*contribution,*params,*bars);

  if( $output eq "gauss" ){

    print STDERR "output gauss\n" if $verbose;
    # turn bars into gaussians:
    &dosplot(*bars,*params,*graph);

    &print_graph(*graph);

  }elsif( $output eq "bars" ){

    print STDERR "output bars\n" if $verbose;
    # output bars:
    &print_graph(*bars);

  }else{

    &usage();

  }
}

sub read_popan{
	local($irr_sym,$irr_ind);
	local($parse_it,$updn);
	while(<>){
		print "-pop: $_" if $Verbose;
		if(/MAJORITY SPIN/){
			$updn = "up";
		}
		if(/MINORITY SPIN/){
			$updn = "dn";
		}
		if(/IRREDUCIBLE REPRESENTATION:\s*(\d+)\s+(\S+)/){
			$irr_ind = $1;
			$irr_sym = $2;
			# push(@irreps,$irr_sym);
			print "IRREP: $irr_sym ($irr_ind)\n" if $verbose;
		}
		if(/CONTRACTIONS IN COLUMNS/){
			$parse_it = 0;
			next;
		}
		if(/L VALUES IN COLUMNS/){
			$parse_it = 1;
		}
		if(s/^     EIGENVALUE/               /){
			if( ! $parse_it ){
				print "popan: ignoring this block!\n" if $verbose;
				next;
			}
			&pop_read_header($irr_sym,*header);
			print "HEAD: = @header\n" if $verbose;
			&pop_read_dat($irr_sym,$updn,*header,*energy,*contribution);
		}
		last if /SUMMARY OVER SPINS AND IRREPS/;
	}
}

sub pop_read_header{
	local($irr,*h) = @_;
	@h = ();
	local($count) = 0;
	do{
		chomp;
		print "pop_read_header: >$_<\n" if $verbose;
		while(s/(\d+)\s*([SPDFGH])//){
			push(@h,"$2<$1>");
			$count++;
			print "pop_read_header: += $2<$1>\n" if $verbose;
		}
		# empty line bug:
		do{ $_ = <> }while( /^$/ );
	}while(/^              /);
	$n_columns{$irr} = $count;
}

sub pop_read_dat{
	local($irr,$updn,*h,*e,*o) = @_;
	local($ind,$eny,$count,$pos,$cont,$lev,$ent);
	$count = 0;
	do{
		if(s/^\s*(\d+)\s+(\S+)/               /){
			print if $Verbose;
			$ind = $1;
			$eny = $2;
			$count++;
			$lev = "$ind,$irr";
			if($updn){
				$lev .= ",$updn";
			}
			$e{$lev} = $eny;
                        print "== $lev ( $eny ) ==\n" if $verbose;
			print STDERR "ind,count=$ind,$count\n" if $Verbose;
			die("pop_read_dat: count mismatch\n") if $ind != $count;
			$pos = 0;
			do{
				print "pop_read_dat: >$_<\n" if $verbose;
				while(s/(\S+)//){
					$ent = "$lev,$h[$pos]";
					if( defined $o{$ent} ){
						print "popan: $ent is already defined!=",
							$o{$ent},"\n";
					}
					$cont = $1;
					$o{$ent} = $cont if ! $cont == 0.0;
					print "+= $ent = $cont\n" if $verbose;
					$pos++;
				}
				# empty line bug in output:
				do{ $_ = <>; print "empty line!" if $verbose; }while( /^$/ );
			}while(/^          /);
		}else{
			die("pop_read_dat: wrong time to call me >$_<\n");
		}
	}while( ! /SUM/ );
	$n_levels{$irr} = $count;
	$degeneracy{$irr} = -1;
}



sub cmdmode{
    print "sub cmdmode entered\n" if $verbose;
    local($key,$val);
    do{
	print "cmdmode: >$_<\n" if $verbose;
	s/^#.*$//;
	if(/^\s*(\w+)\s+(\S+)/){
	    $key = $1;
	    $val = $2;
	    $params{$key} = $val;
	}
	if(/^_GO_/){
            &go();
	   #@graph = 0;
	   #&dosplot(*energy,*occupation,*contribution,*params,*graph);
	   #&print_graph(*graph);
	}
    }while(<>);
}

sub defaults{
    print "sub defaults entered\n" if $verbose;
    local(*par) = @_;

    if( ! defined $par{"PATTERN"} ){
	$par{"PATTERN"} = ".*";
    }
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
#    local($lev);
#    if( $par{"ENERGY_UNITS"} eq "au" ){
#	foreach $lev ( keys %energy ){
#	    $energy{$lev} = $energy{$lev} * $au2ev;
#	}
#	$par{"ENERGY_UNITS"} = "eV";
#    }
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

    local($pair) = 2.0;
    if( $polarized ){
	$pair = 1.0;
    }
    local($irr);

    print "checkdegen: entered %deg\n" if $verbose;

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


    local($checked);
    foreach $irr ( keys %deg ){
	$checked = 0;
	foreach $lett ( keys %assume ){
	    if( $irr =~ /^$lett/ ){
		$checked = 1;
		if( $deg{$irr} != $assume{$lett}*$pair ){
		    warn("Check degeneracy of $irr! Is it really $deg{$irr}?\n");
		    warn("   ... forcing it to ", $assume{$lett}*$pair,"\n");
                    $deg{$irr} = $assume{$lett}*$pair;
		}
	    }
	}
	if( ! $checked ){
	    warn "Dont know irrep $irr\n";
	}
    }
}


sub print_graph{
   print "sub print_graph entered\n" if $verbose;
   print "#GRAPH:	X		F(X)\n";
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

sub barplot{
#    print "sub barplot entered\n" if $VErbose;
    local(*eny,*occ,*cont,*par,*bars) = @_;

    local($i,$lev);
    if( $verbose ){
	foreach $lev ( keys %eny ){
	    print "e( $lev )= $eny{$lev}\n";
	}
	print "\n";
    }

#   &defaults(*par);
    &printpar(*par) if $verbose;

    local($width) = $par{"WIDTH"};
    local($bot)   = $par{"BOTTOM_ENERGY"};
    local($top)   = $par{"TOP_ENERGY"};
    local($e_bot) = $bot - 3*$width;
    local($e_top) = $top + 3*$width;
    local($patt)  = $par{"PATTERN"};
    print STDERR "PATTERN /$patt/o\n" if $verbose;

    local($occ_weight) = $par{"OCCUPIED"};

    local($e,$c,$o,$weight,$ent,$irr);
#	print "BEFORE: size= $#cont\n";
#	foreach $ent ( keys %cont ){
#		delete $cont{$ent} if $cont{$ent} == 0.0;
#	}
#	print "AFTER: size= $#cont\n";

    # energy of matching and non-vanishing contributions
    local(%bar_energy);
    # cumulative contribution of matching levels:
    local(%bar_height);

    foreach $ent ( keys %cont ){
	$c = $cont{$ent};
        # dont proceed with zero contributions:
        next if $c == 0.0;
#	if($c == 0.0){
#		$skipped_entries++;
#		next;
#	}else{
#		$counted_entries++;
#	}
	$ent =~ /^(\d+,[^,]+),(\w+)/;
	$lev = $1;
	if( $polarized ){
	    $lev .= ",$2";
	}
#	die("barplot: eny for $lev is not def!\n")
#	    if ! defined $eny{$lev};
	$e = $eny{$lev};
	print "e( $lev ( $ent ) ) = $e\n" if $verbose;
	next if $e =~ /\*\*\*\*/; # bad format
	next if ($e < $e_bot) || ($e > $e_top) ;
        next if $ent !~ /$patt/o;

	if($occ_weight){
            if( defined $occ{$lev} ){
              $o = $occ{$lev};
            }else{
              warn("barplot: occ for $lev is not def!\n");
              $o = 0.0;
	    }
	}else{
	    $lev =~ /^\d+,([^,]+)/;
	    $irr = $1;
	    die("barplot: bad irr $irr\n")
		if ! defined $degeneracy{$irr};
	    $o = $degeneracy{$irr};
	}
	$weight = $o * $c;
	print ">> $lev $e $o * $c\n" if $verbose;

        # for bar-plot collect the energies and accumulate matching contributions:
        $bar_energy{$lev}  = $e; # re-assigned the same value many times, FIXME by using $eny{$lev}
        $bar_height{$lev} += $weight;

        # Graph summary: compute the weight-center of the graph:
        $graph_center   += $weight * $e;
        $graph_integral += $weight;
    }
    # Graph summary: print weight-center of the graph as a comment:
    if( $graph_integral != 0 ){
      printf "# %10s: %10s %8.3f %8.3f\n", $patt, "Center", $graph_center / $graph_integral, $graph_integral;
    }else{
      printf "# %10s: %10s %8.3f %8.3f\n", $patt, "Empty?", $graph_center, $graph_integral;
    }

    # print significant matching contributions sorted by energy as #-comments:
    foreach $lev ( reverse sort{ $eny{$a} <=> $eny{$b} } keys %bar_height ){
      $weight = $bar_height{$lev};
      next if abs($weight) < $feature;
      printf "# %10s: %10s %8.3f %8.3f\n", $patt, $lev, $eny{$lev}, $weight;
    }

    # output bars by energy => weight hash:
    foreach $lev ( keys %bar_energy ){
      $e = $bar_energy{$lev};
      $bars{$e} = $bar_height{$lev};
    }
}

sub dosplot{ # turns bar plot into dos plot
#    print "sub dosplot entered\n" if $VErbose;
    local(*bars,*par,*graph) = @_;

    &printpar(*par) if $verbose;

    local($width) = $par{"WIDTH"};
    local($bot)   = $par{"BOTTOM_ENERGY"};
    local($top)   = $par{"TOP_ENERGY"};
    local($e_bot) = $bot - 3*$width;
    local($e_top) = $top + 3*$width;
    local($np)    = $par{"NPOINTS"};

    local($dE) = ($top - $bot) / $np;

    local(@Xs,@Ys);
    @Xs=0; @Ys=0;

    for($i=0; $i<$np; $i++){
	$Xs[$i] = $bot + ( $i + 0.5 ) * $dE;
    }

    local($e,$x,$weight);

    # here the bars are turned into gaussians:
    foreach $e ( reverse sort keys %bars ){
      next if $e > $e_top || $e < $e_bot;
      $weight = $bars{$e};
      next if $weight == 0.0;

      for($i=0; $i<$np; $i++){
          $x = ( $Xs[$i] - $e ) / $width;
          $Ys[$i] += $weight * &invexp2($x);
      }
    }

    # output hash x => y:
    for($i=0; $i<$np; $i++){
	$graph{$Xs[$i]} = $Ys[$i];
    }
}

sub invexp2{
    # exp( - x^2 )
    local($x) = @_;
    $x *= $x;
    ($x>50.0)?(0.0):(exp(-$x));
}

sub usage{
        local($msg) = @_;
        local($me) = $0;

        warn("$0: $msg\n");

        open(ME,$me) || die("$0: cannot read my source\n");
        while(<ME>){
                if(s/^#usage://){
                        print STDERR;
                }
        }
        close(ME);
        die();
}




