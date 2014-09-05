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

# removes_subroutine_haeds TB 5/97
#
# usage :   remove_subroutine_haeds <filename> [<filename> ... ]
#
# Edits fortran source file of module in the form of 
# "documentation/standard_module.f90" in order to remove
# subroutine haeders in commented form from the definition of public
# interface of modules.
# Edits <filename> and saves original file in <filename>.old.
# The file must include the following lines (beginning with !) :
#   ! insert_subroutine_heads begin
#   ! insert_subroutine_heads end
# everything between these two lines and the lines themselves are removed
# from the source file



###########  subroutines


sub help_message {
    print  <<"ENDMESSAGE";

removes_subroutine_haeds TB 5/97

usage :   remove_subroutine_haeds <filename> [<filename> ... ]

Edits fortran source file of module in the form of 
"documentation/standard_module.f90" in order to remove
subroutine haeders in commented form from the definition of public
interface of modules.
Edits <filename> and saves original file in <filename>.old.
The file must include the following lines (beginning with !) :
  ! insert_subroutine_heads begin
  ! insert_subroutine_heads end
everything between these two lines and the lines themselves are removed
from the source file

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
    open(FILETMP,"> $filetmp") || die "Error: $filetmp can not be generated: $!\n";


    # write beginning of file to $filetmp
    open(FILE,"$file") || die "Error: $file not found: $!\n";
    WRITEBEGINNING: while(<FILE>) {
	if ( /^\s*!\s*insert_subroutine_heads\s+begin/ ) {
	    last WRITEBEGINNING;
	}
	print FILETMP;
    }
    # write end of file to $filetmp
    local ( $found ) = 0;
    while(<FILE>) {
	if ( $found ) {
	    print FILETMP;
	}
	if ( /^\s*!\s*insert_subroutine_heads\s+end/ ) {
	    $found = 1;
	}
    }
    close(FILE);

    close(FILETMP);


    rename($file,"$file.old") || warn "can not move $file to $file.old: $!\n";
    rename($filetmp,$file) || die  "can not move $filetmp to $file: $!\n";
    print "$file was edited and old file saved as $file.old\n";


}
