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

# make_interface TB 5/97
#
# usage :               make_interface <filename> [<filename> ... ]
# for this help mesage: make_interface
#
# Cuts intefrace information from Fortran source file in the form of 
# "documentation/standard_module.f90" or "documentation/standard_subroutine.f90".
#
# The f90 source files must contain comment lines indicating start and end
# of interface (beginning with !):
#   ! Public interface of module
#   ! End of public interface of module
# It can additionally contain pairs of lines indicating an interrupt in the
# public interface that will not be printed (lines beginning with !):
#   !== Interrupt of public interface of module =========
#   !== Interrupt end of public interface of module =====
# For fixed form sources, the "!" may be replaced by "C".
# 
# In case the fortran file contains a module with public subroutines or functions,
# the file must include the following lines (beginning with !) :
#   !------------ public functions and subroutines ------------------
# make_interface reads the names of the public subroutines from the 
# "public" statement following the first line
# Each of the subroutine declarations must include a line 
#   !** End of interface *****************************************
# The subroutine haeds are read in from the subroutine or function
# statement up to this line.
# 
# The following lines are removed from f90 sources:
#  implicit none
#  !------------ Modules used ------------
# 
# It is allowed to pass filenames ending with ".o" as arguments. In this case, 
# the programm looks for the correspondibg source files ending with ".f90" or ".f".
#




###########  subroutines


sub help_message {
    print  <<"ENDMESSAGE";

 make_interface TB 5/97

 usage :               make_interface <filename> [<filename> ... ]
 for this help mesage: make_interface

 Cuts intefrace information from Fortran source file in the form of 
 "documentation/standard_module.f90" or "documentation/standard_subroutine.f90".

 The f90 source files must contain comment lines indicating start and end
 of interface (beginning with !):
   ! Public interface of module
   ! End of public interface of module
 It can additionally contain pairs of lines indicating an interrupt in the
 public interface that will not be printed (lines beginning with !):
   !== Interrupt of public interface of module =========
   !== Interrupt end of public interface of module =====
 For fixed form sources, the "!" may be replaced by "C".
 
 In case the fortran file contains a module with public subroutines or functions,
 the file must include the following lines (beginning with !) :
   !------------ public functions and subroutines ------------------
 make_interface reads the names of the public subroutines from the 
 "public" statement following the first line
 Each of the subroutine declarations must include a line 
   !** End of interface *****************************************
 The subroutine haeds are read in from the subroutine or function
 statement up to this line.
 
 The following lines are removed from f90 sources:
  implicit none
  !------------ Modules used ------------
 
 It is allowed to pass filenames ending with ".o" as arguments. In this case, 
 the programm looks for the correspondibg source files ending with ".f90" or ".f".

ENDMESSAGE
}









########### executable part




# read command line arguments
if ( @ARGV ) {
    @files = @ARGV;
}
else {
    &help_message;
    exit;
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
	    die "make_interface: Error: no source file for $file .\n";
	}
    }
    elsif ( ! -e $file ) {
	die "make_interface: Error: $file not found.\n";
    }
}



# read list of object files @objfiles
foreach $file ( @files ) {

    local( $docfile ) = $file;
    $docfile =~ s/\.\w+$/.doc/;

    local(@subroutines);

    # read line with list of subroutines
    open(FILE,$file) || die "make_interface: Error: $file not found: $!\n";
    local( $found ) = 0;
    local( $publicfound ) = 0;
    local( $ismodule ) = 0;
    READSUBLIST: while(<FILE>) {
	if ( /^\s*[mM][oO][dD][uU][lL][eE]\s+\w+/ ) {
	    $ismodul = 1;
	}
	if ( /^\s*!-*\s*public\s+functions\s+and\s+subroutines\s*-*$/ ) {
	   $found = 1;
	}
	elsif ( $found ) {
	    if ( /^\s*public\s*:*\s+([^&!]*)&/ ) {
		splice(@subroutines,0,0,split(/\s*,\s*/,$1));
		$publicfound = 1;
	    }
	    elsif ( /^\s*PUBLIC\s*:*\s+([^&!]*)&/ ) {
		splice(@subroutines,0,0,split(/\s*,\s*/,$1));
		$publicfound = 1;
	    }
	    elsif ( /^\s*public\s*:*\s+&/ ) {
		$publicfound = 1;
	    }
	    elsif ( /^\s*PUBLIC\s*:*\s+&/ ) {
		$publicfound = 1;
	    }
	    elsif ( /^\s*public\s*:*\s+([^!]*)/ ) {
		splice(@subroutines,0,0,split(/\s*,\s*/,$1));
		last READSUBLIST;
	    }
	    elsif ( /^\s*PUBLIC\s*:*\s+([^!]*)/ ) {
		splice(@subroutines,0,0,split(/\s*,\s*/,$1));
		last READSUBLIST;
	    }
	    elsif ( $publicfound ) {
		if ( /^\s*([^&!]*)&/ ) {
		    splice(@subroutines,0,0,split(/\s*,\s*/,$1));
		}
		elsif ( /^\s*([^!]*)/ ) {
		    splice(@subroutines,0,0,split(/\s*,\s*/,$1));
		    last READSUBLIST;
		}
	    }
	}
    }
    close(FILE);


    if ( @subroutines ) {
	@subroutines = sort @subroutines;
    }
    elsif ( $ismodul ) {
	warn "make_interface: Warning: no declaration of public subroutines found in $file\n";
    }



    open(DOCFILE,"> $docfile") || die "make_interface: Error: $docfile can not be generated: $!\n";


# search interfaes in @files and print to $outfile
    open(FILE,$file) || die "make_interface: Error: $file not found: $!\n";
    $found = 0;
    local ( $interrupt ) = 0;
    local( $insertfound ) = 0;
    READFILE: while (<FILE>) {
	if ( /^\s*!\s*[Ee]nd\s+of\s+[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
	    last READFILE;
	}
	elsif ( /^[cC]\s*[Ee]nd\s+of\s+[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
	    last READFILE;
	}
	if ( $found ) {
	    if ( $interrupt ) {
		if ( /^\s*!=*\s*[Ii]nterrupt\s*[Ee]nd\s+of\s+[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
		    $interrupt = 0;
		}
		elsif ( /^[cC]=*\s*[Ii]nterrupt\s*[Ee]nd\s+of\s+[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
		    $interrupt = 0;
		}
	    }
	    elsif ( $insertfound ) {
		if ( /^\s*!\s*insert_subroutine_heads\s+end/ ) {
		    $insertfound = 0;
		}
	    }
	    elsif ( /^\s*!\s*insert_subroutine_heads\s+begin/ ) {
		$insertfound = 1;
	    }
	    elsif ( /^\s*!=*\s*[Ii]nterrupt\s+of\s+[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
		$interrupt = 1;
	    }
	    elsif ( /^[cC]=*\s*[Ii]nterrupt\s+of\s+[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
		$interrupt = 1;
	    }
	    elsif ( ! /^\s*!\s* [Ii][Mm][Pp][Ll][Ii][Cc][Ii][Tt]\s+[Nn][Oo][Nn][Ee]/ && 
                    ! /^\s*[Ii][Mm][Pp][Ll][Ii][Cc][Ii][Tt]\s+[Nn][Oo][Nn][Ee]/ && 
		    ! /^\s*!\s*!?-*\s*[Mm]odules\s+used\s*-*/ ) {
		print DOCFILE;
	    }
	}
	elsif ( /^\s*!=*\s*[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
	    $found = 1;
        }
	elsif ( /^[cC]=*\s*[Pp]ublic\s+[Ii]nterface\s+of\s+[Mm]odule/ ) {
	    $found = 1;
	}
    }
    if ( ! $found ) {
	warn "Warning: no interface found for $file \n";
    }
    elsif ( ! FILE ) {
	warn "Warning: no end of interface found for $file \n";
    }
    close (FILE);


    if ( @subroutines ) {
	print DOCFILE "\n\n=============== Public functions and subroutines: ===============\n\n";
    }


    # look for subroutine heads and write them to $docfile
    foreach $subroutine ( @subroutines ) {
	$subroutine =~ s/\n$//;
	print DOCFILE "\n";
	open(FILE,"$file") || die "Error: $file not found: $!\n";
	local ( $found ) = 0;
	local ( $firstline ) = 0;
        READSUBHAED: while (<FILE>) {
	    if ( /^\s*subroutine\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
		$firstline = 1;
	    }
	    elsif ( /^\s*SUBROUTINE\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
		$firstline = 1;
	    }
	    elsif ( /^[^!]*function\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
		$firstline = 1;
	    }
	    elsif ( /^[^!]*FUNCTION\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
		$firstline = 1;
	    }
	    if ( $found ) {
		if ( /\s*!\**\s*End\s+of\s+interface\s*\**/ ) {
		    last READSUBHAED;
		}
		elsif ( /^\s*subroutine/ && ! $firstline ) {
		    warn "Warning: no end of subroutine $subroutine interface found for $file \n";
		    last READSUBHAED;
		}
		elsif ( /^\s*SUBROUTINE/ && ! $firstline ) {
		    warn "Warning: no end of subroutine $subroutine interface found for $file \n";
		    last READSUBHAED;
		}
		elsif ( /^\[^!]*FUNCTION/ && ! $firstline ) {
		    warn "Warning: no end of function $subroutine interface found for $file \n";
		    last READSUBHAED;
		}
		elsif ( /^\[^!]*function/ && ! $firstline ) {
		    warn "Warning: no end of function $subroutine interface found for $file \n";
		    last READSUBHAED;
		}
		else {
		    print DOCFILE;
		    $firstline = 0;
		}
	    }
	}
	close(FILE);
	if ( ! $found ) {
	    warn "WARNING: declaration of $subroutine not found in $file \n";
	}
    }
    print DOCFILE "\n";

    close(DOCFILE);

}
