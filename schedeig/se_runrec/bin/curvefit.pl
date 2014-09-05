#!/usr/bin/perl -w
#
# ParaGauss,  a program package  for high-performance  computations of
# molecular systems
#
# Copyright (C) 2014     T. Belling,     T. Grauschopf,     S. Krüger,
# F. Nörtemann, M. Staufer,  M. Mayer, V. A. Nasluzov, U. Birkenheuer,
# A. Hu, A. V. Matveev, A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman,
# D. I. Ganyushin,   T. Kerdcharoen,   A. Woiterski,  A. B. Gordienko,
# S. Majumder,     M. H. i Rotllant,     R. Ramakrishnan,    G. Dixit,
# A. Nikodem, T. Soini, M. Roderus, N. Rösch
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
# Please see the accompanying LICENSE file for further information.
#

use warnings;
use strict;

use File::Copy;

use lib './lib/';
use Algorithm::CurveFit;


my $file_in = $ARGV[0];
my $file_out_plain = $ARGV[1];
my $file_out_f90 = $ARGV[2];

my %containerHash;

open( my $filehandle, "<", $file_in ) or die "Error opening file $file_in";



while(<$filehandle>){

    s/#.*//;
    next if /^(\s)*$/;
    chomp;

    my $line = $_;

    

    my @partial = split /\t/, $line;

    if( @partial == 3 ){

	push( @{ $containerHash{$partial[0]}{x} }, $partial[1]);
	push( @{ $containerHash{$partial[0]}{y} }, $partial[2]);
	
    }
}


close($filehandle) or die "Error closing file $!";

my $count = scalar(keys %containerHash);


#system("cp", "se_timefunction_module.f90.beg.template", "");
copy("se_timefunction_module.f90.beg.template", $file_out_f90) or die "Copy failed: $!";
open( my $filehandle_f90, ">>", $file_out_f90 ) or die "Open file failed: $!";
print $filehandle_f90 "!... <<START>> this section is automatically generated\n    if( len(filename) .eq. 0   .or.   ios .ne. 0 ) then\n";
#printf $filehandle_f90 "      allocate( timeFunctions($count) )\n\n      do i=1, $count\n        allocate( timeFunctions(i)%%coefficients(4) )\n      end do\n\n";
printf $filehandle_f90 "      call timefunction_allocate( timeFunctions, $count, 3 )\n\n";

# open file for plain coefficients
open($filehandle, ">", $file_out_plain) or die "Error opening file $!";

print $filehandle "# Polynomial of degree 3: c0 + c1 * x^1 + c2 * x^2 + c3 * x^3\n#\n";
print $filehandle "# The data are stored in the following maner:\n";
print $filehandle "## Number of Processors\n## c0    c1    c2    c3\n\n";
print $filehandle "N: $count\nDEG: 3\n\n";


my $i=1;
foreach my $key ( sort{$a <=> $b} keys %containerHash){
    #print "$key\n";

    my @xdata = @{ $containerHash{$key}{x} };
    my @ydata = @{ $containerHash{$key}{y} };
    
    my $formula = 'c0 + c1 * x^1 + c2 * x^2 + c3 * x^3';
    my $variable = 'x';

    my @parameters = (
	# Name    Guess   Accuracy
	['c0',     0.1,     0.0001],
	['c1',     0.1,     0.0005],
	['c2',     0.1,     0.0005],
	['c3',     0.1,     0.0005],
	);
    
    my $max_iter = 100;

    my $square_residual = Algorithm::CurveFit->curve_fit(
	formula            => $formula,
	params             => \@parameters,
	variable           => $variable,
	xdata              => \@xdata,
	ydata              => \@ydata,
	maximum_iterations => $max_iter,
	);

    printf $filehandle  "P: $key\nC: %+0.6E    %+0.6E    %+0.6E    %+0.6E\n\n", $parameters[0][1], $parameters[1][1], $parameters[2][1], $parameters[3][1];

    printf $filehandle_f90 "      timeFunctions($i)%%procCount = $key\n";
    printf $filehandle_f90 "      timeFunctions($i)%%coefficients(1) = %+0.6E\n", $parameters[0][1];
    printf $filehandle_f90 "      timeFunctions($i)%%coefficients(2) = %+0.6E\n", $parameters[1][1];
    printf $filehandle_f90 "      timeFunctions($i)%%coefficients(3) = %+0.6E\n", $parameters[2][1];
    printf $filehandle_f90 "      timeFunctions($i)%%coefficients(4) = %+0.6E\n", $parameters[3][1];
    printf $filehandle_f90 "      timeFunctions($i)%%scaleFactor = 1.0\n";
    printf $filehandle_f90 "      timeFunctions($i)%%polyDegr = 3\n";
    printf $filehandle_f90 "      timeFunctions($i)%%minSize = getMinSize( $key )\n\n";
    $i++;
}

close($filehandle) or die "Error closing file $!";

print $filehandle_f90 "    end if\n!... <<END>> this section is automatically generated\n";
close($filehandle_f90) or die "Error closing file $!";


system("cat se_timefunction_module.f90.end.template >> $file_out_f90");

print "Polynomials generated. Stored to files: $file_out_plain and $file_out_f90\n";
