#!/usr/bin/perl -w
#
# ParaGauss, a program package for high-performance computations
# of molecular systems
# Copyright (C) 2014
# T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
# M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
# A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
# T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
# M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
# M. Roderus, N. Rösch
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 as published
# by the Free Software Foundation [1].
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# [1] http://www.gnu.org/licenses/gpl-2.0.html
#
# Please see the accompanying LICENSE file for further information.
#

use warnings;
use strict;



my $inFile = $ARGV[0];
my $outFile = $ARGV[1];
my %matSizeHash;

open(my $fileHandle, "<", $inFile) or die "Error opening file $inFile";


while( <$fileHandle> ){
    my $lineBuffer = $_;

    # Next while-iteration unless line contains 'Matrix dimension'
    next unless $lineBuffer =~ /Code:/;
    next if $lineBuffer =~ /^(#)/;

    my @lineBuffer = split(/\s+/, $lineBuffer);

    $matSizeHash{$lineBuffer[2]}{$lineBuffer[3]} = $lineBuffer[4];
}


close( $fileHandle );



open($fileHandle, ">", $outFile) or die "Error opening file $inFile";

for my $procKey ( sort{$a <=> $b} keys %matSizeHash ){
    for my $matKey ( sort{$a <=> $b} keys %{$matSizeHash{$procKey}}){
	print $fileHandle "$procKey\t$matKey\t$matSizeHash{$procKey}{$matKey}\n";
    }
    print $fileHandle "\n";
}

close( $fileHandle );



sub log2 {
    my $n = shift;
    return log($n)/log(2);
}
sub log3 {
    my $n = shift;
    return log($n)/log(3);
}
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
