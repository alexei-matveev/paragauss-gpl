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

# insert_subroutine_haeds TB 10/95
#
# usage :   insert_subroutine_haeds <filename> [<filename> ... ]
#
# Edits fortran source file of module in the form of 
# "documentation/standard_module.f90" in order to include
# subroutine haeders in commented form to the definition of public
# interface of modules.
# Edits <filename> and saves original file in <filename>.old.
# The file must include the following lines (beginning with !) :
#   !------------ public functions and subroutines ------------------
#   ! insert_subroutine_heads begin
#   ! insert_subroutine_heads end
# Makedependent reads the names of the public subroutines from the 
# "public" statement following the first line
# and inserts the subroutine haeds between the start and end line,
# removing anything old between these lines.
# Each of the subroutine declarations must include a line 
#   !** End of interface *****************************************
# The subroutine haeds are read in from the subroutine or function
# statement up to this line.



###########  subroutines


sub help_message {
    print  <<"ENDMESSAGE";

insert_subroutine_haeds TB 10/95

usage :   insert_subroutine_haeds <filename> [<filename> ... ]

Edits fortran source file of module in the form of 
"documentation/standard_module.f90" in order to include
subroutine haeders in commented form to the definition of public
interface of modules.
Edits <filename> and saves original file in <filename>.old.
The file must include the following lines (beginning with !) :
  !------------ public functions and subroutines ------------------
  ! insert_subroutine_heads begin
  ! insert_subroutine_heads end
Makedependent reads the names of the public subroutines from the 
"public" statement following the first line
and inserts the subroutine haeds between the start and end line,
removing anything old between these lines.
Each of the subroutine declarations must include a line 
  !** End of interface *****************************************
The subroutine haeds are read in from the subroutine or function
statement up to this line.

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






# read list of object files @objfiles
foreach $file ( @files ) {

    local( $filetmp ) = $file.".tmp";
    # check if file "<file>.tmp" exists, if so die
    if ( -e $filetmp ) {
	die "Error: file $filetmp exists. This file would be overwritten.";
    }

    local(@subroutines);

    # read line with list of subroutines
    open(FILE,$file) || die "Error: $file not found: $!\n";
    local( $found ) = 0;
    local( $publicfound ) = 0;
    READSUBLIST: while(<FILE>) {
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

    if ( ! @subroutines ) {
	warn "Error: no declaration of public subroutines found in $file\n";
	&help_message;
	die "exiting";
    }


    @subroutines = sort @subroutines;


    open(FILETMP,"> $filetmp") || die "Error: $filetmp can not be generated: $!\n";


    # write beginning of file to $filetmp
    open(FILE,"$file") || die "Error: $file not found: $!\n";
    WRITEBEGINNING: while(<FILE>) {
	print FILETMP;
	if ( /^\s*!\s*insert_subroutine_heads\s+begin/ ) {
	    last WRITEBEGINNING;
	}
    }
    close(FILE);

    # look for subroutine heads and write them to $filetmp
    foreach $subroutine ( @subroutines ) {
	$subroutine =~ s/\n$//;
	print FILETMP "\n";
	open(FILE,"$file") || die "Error: $file not found: $!\n";
	local ( $found ) = 0;
        READSUBHAED: while (<FILE>) {
	    if ( /^\s*subroutine\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
	    }
	    elsif ( /^\s*SUBROUTINE\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
	    }
	    elsif ( /^[^!]*function\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
	    }
	    elsif ( /^[^!]*FUNCTION\s+$subroutine[\s\n\(]/ ) {
		$found = 1;
	    }
	    if ( $found ) {
		if ( /\s*!\**\s*End\s+of\s+interface\s*\**/ ) {
		    last READSUBHAED;
		}
		else {
		    print FILETMP "! $_";
		}
	    }
	}
	close(FILE);
	if ( ! $found ) {
	    warn "WARNING: declaration of $subroutine not found in $file \n";
	}
    }
    print FILETMP "\n";

    # write end of file to $filetmp
    open(FILE,"$file") || die "Error: $file not found: $!\n";
    $found = 0;
    while(<FILE>) {
	if ( /^\s*!\s*insert_subroutine_heads\s+end/ ) {
	    $found = 1;
	}
	if ( $found ) {
	    print FILETMP;
	}
    }
    close(FILE);


    close(FILETMP);


    rename($file,"$file.old") || warn "can not move $file to $file.old: $!\n";
    rename($filetmp,$file) || die  "can not move $filetmp to $file: $!\n";
    print "$file was edited and old file saved as $file.old\n";


}
