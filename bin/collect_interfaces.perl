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
# collect_interfaces TB 10/95
#
# usage :   help:           collect_interfaces
#           argument mode:  collect_interfaces [-o <outputfile>] <filename> [<filename> ...]
#           makefile mode:  collect_interfaces [-o <outputfile>] -m [<filename>]
#           listfile mode:  collect_interfaces [-o <outputfile>] -l [<filename>]
#
#
# Collects the public interfaces of f90 source files and writes them to
# "docomentation/interfaces" or any specified <outputfile>.
# The public interfaces must be extracted and put to .doc files before
# with the help of program make_interface (type without arguments to obtain help).
#
# 
# It is allowed to pass filenames ending with ".o" as arguments. In this case, 
# the programm looks for the correspondibg source files ending with ".f90" or ".f".
#
#
# Three modes of selecting the files to collect the interface from exist:
#
# Argument mode:
# Takes the files given as command line arguments.
#
# Listfile mode:
# Reads in file "listfile.collect_interfaces" or <filename> that should contain 
# the name of one fortran file in each line.
#
# Makefile mode:
# The file "Makefile" or <filename> must include the following line (beginning with #) :
#   # for collect_interfaces: <objectmakroname> [ <objectmakroname> ... ]
# collect_interfaces reads the files to process from the objectmakro definitions.
# 
# to obtain help statement, type "collect_interfaces".



###########  subroutines


sub shorten_list {
	local(@newlist);
	local($last) = "";
        local($elem);
	foreach $elem (@_) {
		if ($last ne $elem) {
			push(@newlist,$elem);
			$last = $elem;
		}
	}
	@newlist;
}



sub nostandardname {
	local($file,$subrt) = @_ ;
	local($nostandard) = 1;
	$subrt =~ tr/A-Z/a-z/;
	$file =~ tr/A-Z/a-z/;
	if ( $file =~ /(\w+_)module.doc$/ ) {
	    local ($modulename) = $1;
	    if ( $subrt=~/^$modulename/ ) {
		$nostandard  = 0;
	    }
	}
	$nostandard;
}



sub help_message {
    print  <<"ENDMESSAGE";

collect_interfaces TB 10/95

usage :   help:           collect_interfaces
          argument mode:  collect_interfaces [-o <outputfile>] <filename> [<filename> ...]
          makefile mode:  collect_interfaces [-o <outputfile>] -m [<filename>]
          listfile mode:  collect_interfaces [-o <outputfile>] -l [<filename>]

Collects the public interfaces of f90 source files and writes them to
"docomentation/interfaces" or any specified <outputfile>.
The public interfaces must be extracted and put to .doc files before
with the help of program make_interface (type without arguments to obtain help).


It is allowed to pass filenames ending with ".o" as arguments. In this case, 
the programm looks for the correspondibg source files ending with ".f90" or ".f".


Three modes of selecting the files to collect the interface from exist:

Argument mode:
Takes the files given as command line arguments.

Listfile mode:
Reads in file "listfile.collect_interfaces" or <filename> that should contain 
the name of one fortran file in each line.

Makefile mode:
The file "Makefile" or <filename> must include the following line (beginning with #) :
  # for collect_interfaces: <objectmakroname> [ <objectmakroname> ... ]
collect_interfaces reads the files to process from the objectmakro definitions.


ENDMESSAGE
}









########### executable part




# read command line arguments
$mode = "argmode";
$listfile = "listfile.collect_interfaces";
$makefile = "Makefile";
$outfile  = "documentation/interfaces";
if ( ! @ARGV ) {
    &help_message;
    exit;
}
READARGS: while ( @ARGV ) {
    $arg1 = shift(@ARGV);
    # -o is used to specify output file
    if ( $arg1 eq "-o" ) {
	if ( ! @ARGV ) {
	    &help_message;
	    die "Error: flag -o but no output file given as argument\n";
	}
	$outfile = shift(@ARGV);
    }
    # -l invokes list mode
    elsif ( $arg1 eq "-l" ) {
	$mode = "listmode";
	if ( @ARGV ) {
	    $arg1 = shift(@ARGV);
	    if ($arg1 ne "-o") {
		$listfile = $arg1;
		if ( @ARGV ) {
		    &help_message;
		    die "Error: to many arguments\n";
		}
	    }
	}
    }
    # -m invokes list mode
    elsif ( $arg1 eq "-m" ) {
	$mode = "makefilemode";
	if ( @ARGV ) {
	    $arg1 = shift(@ARGV);
	    if ($arg1 ne "-o") {
		$makefile = $arg1;
		if ( @ARGV ) {
		    &help_message;
		    die "Error: to many arguments\n";
		}
	    }
	}
    }
    else {
	@files = @ARGV;
	unshift(@files,$arg1);
	last READARGS;
    }
}


if ( $mode eq "argmode" ) {
    if ( ! @files ) {
	&help_message;
	die "Error: not enough arguments\n";
    }
}


# read list of object files @files
if ( $mode eq "makefilemode" ) {		# makefile mode
    local(@objmakros);

    # read line with list of object makros
    open(MAKEFILE,$makefile) || die "Error: $makefile not found: $!\n";
    while(<MAKEFILE>) {
	if ( /^# for collect_interfaces:\s+(.+)/ ) {
	    @objmakros = split(/\s+/,$1);
	}
    }
    close(MAKEFILE);
    if ( ! @objmakros ) {
	&help_message;
	die "Error: $makefile does not contain line listing object makros \n";
    }

    # read in object files from makro definitions
    foreach $objmakro ( @objmakros ) {
	open(MAKEFILE,$makefile) || die "Error: $makefile not found: $!\n";
	local( $found ) = 0;
	REODOBYMAKRO: while(<MAKEFILE>) {
	    if ( /^$objmakro\s*=\s+([^\\]*)$/ ) {
		splice(@files,0,0,split(/\s+/,$1));
		$found = 1;
		last REODOBYMAKRO;
	    }
	    if ( /^$objmakro\s*=\s+([^\\]*)[\\]/ ) {
		splice(@files,0,0,split(/\s+/,$1));
		$found = 1;
	    }
	    else {
		if ( $found ) {
		    if ( /^\s+([^\\]*)$/ ) {
			splice(@files,0,0,split(/\s+/,$1));
			last REODOBYMAKRO;
		    }
		    if ( /^\s+([^\\]*)[\\]/ ) {
			splice(@files,0,0,split(/\s+/,$1));
		    }
		    else {
			last REODOBYMAKRO;
		    }
		}
	    }
	}
	close(MAKEFILE);
	if ( ! $found ) {
	    die "Error: object makro definition of $objmakro not found in $makefile";
	}
    }
}

if ( $mode eq "listmode" ) {				# list mode
    open(LISTFILE,$listfile) || die "Error: $listfile not found: $!\n";
    while(<LISTFILE>) {
	chop;
	push(@files,$_);
    }
    close(LISTFILE);
}


# transform filenames ending with .o to filenames ending with .f90 or .f
foreach $file ( @files ) {
    if ( $file =~ /\.o$/ ) {
	local($ffile) = $file;
	$ffile =~ s/\.o$/.f/;
	local($fpfile) = $file;
	$fpfile =~ s/\.o$/.fp/;
	local($f90file) = $file;
	$f90file =~ s/\.o$/.f90/;
	local($f90pfile) = $file;
	$f90pfile =~ s/\.o$/.f90p/;
	if ( -e $f90pfile ) {
	    $file = $f90pfile;
	}
	elsif ( -e $f90file ) {
	    $file = $f90file;
	}
	elsif ( -e $fpfile ) {
	    $file = $fpfile;
	}
	elsif ( -e $ffile ) {
	    $file = $ffile;
	}
	else {
	    die "no source file for $file .\n";
	}
    }
    elsif ( ! -e $file ) {
	die "$file not found.\n";
    }
}


# take away directory part from @files and store complete name in $fullname{$file}
foreach $file ( @files ) {
    local($newname) = $file;
    $newname =~ s/^(.*\/)+//;
    $fullname{$newname} = $file;
    $file = $newname;
}


# sort and shorten @files 
@files = sort @files;
@files = &shorten_list(@files);



# extract names of public functions and subroutines from doc-files
foreach $file ( @files ) {
    # print "processing $file \n";
    local ( $docfile ) = $fullname{$file};
    $docfile =~ s/\.\w+$/.doc/;
    open(FILE,$docfile) || die "Error: $docfile not found: $!\n";
    local ( $subroutinesection ) = 0;
    while (<FILE>) {
	if ( $subroutinesection ) {
	    if ( /^\s*[sS][uU][bB][rR][oO][uU][tT][iI][nN][eE]\s+(\w+)/ ) {
		local ( $subrt ) = $1;
		if (  &nostandardname($docfile,$subrt) ) {
		    $subrt =~ tr/A-Z/a-z/;
		    push(@subroutines,$subrt);
		    $fullname{$subrt} = "contained in $file";
		}
	    }
	    elsif ( /^[^!]*[fF][uU][nN][cC][tT][iI][oO][nN]\s+(\w+)/ ) {
		local ( $subrt ) = $1;
		if (  &nostandardname($docfile,$subrt) ) {
		    $subrt =~ tr/A-Z/a-z/;
		    push(@subroutines,$subrt);
		    $fullname{$subrt} = "contained in $file";
		}
	    }
	}
	elsif ( /^=============== Public functions and subroutines: ===============/ ) {
	    $subroutinesection = 1;
	}
    }
    close (FILE);
}


# add those subroutines to files list
push(@files,@subroutines);
@files = sort @files;



# print @files in alphatic order to $outfile
open(OUTFILE,"> $outfile") || die "Error: $outfile can not be generated: $!\n";
require "ctime.pl";
print OUTFILE "\n\nThis file was automatically created by collect_interfaces\nat ",&ctime(time);
foreach $file ( @files ) {
    if ( $fullname{$file} =~ /^contained in / ) {
	print OUTFILE "\n\n\n\n\n         **********************************************\n";
	print OUTFILE "             $file\n";
	print OUTFILE "             $fullname{$file}\n";
	print OUTFILE "         **********************************************\n";
    }
    else {
	local ( $docfile ) = $fullname{$file};
	$docfile =~ s/\.\w+$/.doc/;
	open(FILE,$docfile) || die "Error: $docfile not found: $!\n";
	print OUTFILE "\n\n\n\n\n         **********************************************\n";
	print OUTFILE "             $file\n";
	print OUTFILE "         **********************************************\n\n\n";
	while (<FILE>) {
	    print OUTFILE;
	}
	close (FILE);
    }
}
print OUTFILE "\n\n\n";
close(OUTFILE);
