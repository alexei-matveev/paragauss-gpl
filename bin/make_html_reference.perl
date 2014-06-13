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
# make_html_reference TB 2/98
#
# usage :   help:           make_html_reference
#           argument mode:  make_html_reference [-o <outputfile>] [-v <version>] <filename> [<filename> ...]
#           makefile mode:  make_html_reference [-o <outputfile>] [-v <version>] -m [<filename>]
#           listfile mode:  make_html_reference [-o <outputfile>] [-v <version>] -l [<filename>]
#
# Collects the public interfaces of f90 source files and writes them to
# "docomentation/html" or any specified <outputfile>.
# The public interfaces must be extracted and put to .doc files before
# with the help of program make_interface (type without arguments to obtain help).
#
# It is allowed to pass filenames ending with ".o" as arguments. In this case, 
# the programm looks for the correspondibg source files ending with ".f90" or ".f".
#
# The switch [-o <outputfile>] allows to specify the output file.
# The switch [-v <version>] allows to specify a Progrm Version written to haeder of output.
#
#
# Three modes of selecting the files to collect the interface from exist:
#
# Argument mode:
# Takes the files given as command line arguments.
#
# Listfile mode:
# Reads in file "listfile.make_html_reference" or <filename> that should contain 
# the name of one fortran file in each line.
#
# Makefile mode:
# The file "Makefile" or <filename> must include the following line (beginning with #) :
#   # for make_html_reference: <objectmakroname> [ <objectmakroname> ... ]
# make_html_reference reads the files to process from the objectmakro definitions.
# 
# to obtain help statement, type "make_html_reference".



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



sub help_message {
    print  <<"ENDMESSAGE";

make_html_reference TB 10/95

usage :   help:           make_html_reference
          argument mode:  make_html_reference [-o <outputfile>] [-v <version>] <filename> [<filename> ...]
          makefile mode:  make_html_reference [-o <outputfile>] [-v <version>] -m [<filename>]
          listfile mode:  make_html_reference [-o <outputfile>] [-v <version>] -l [<filename>]

Collects the public interfaces of f90 source files and writes them to
"docomentation/html" or any specified <outputfile>.
The public interfaces must be extracted and put to .doc files before
with the help of program make_interface (type without arguments to obtain help).

It is allowed to pass filenames ending with ".o" as arguments. In this case, 
the programm looks for the correspondibg source files ending with ".f90" or ".f".

The switch [-o <outputfile>] allows to specify the output file.
The switch [-v <version>] allows to specify a Progrm Version written to haeder of output.


Three modes of selecting the files to collect the interface from exist:

Argument mode:
Takes the files given as command line arguments.

Listfile mode:
Reads in file "listfile.make_html_reference" or <filename> that should contain 
the name of one fortran file in each line.

Makefile mode:
The file "Makefile" or <filename> must include the following line (beginning with #) :
  # for make_html_reference: <objectmakroname> [ <objectmakroname> ... ]
make_html_reference reads the files to process from the objectmakro definitions.


ENDMESSAGE
}









########### executable part




# read command line arguments
$mode = "argmode";
$listfile = "listfile.make_html_reference";
$makefile = "Makefile";
$outfile  = "documentation/reference.html";
$version = "";
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
    # -m invokes makefile mode
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
    # -v allows to spezify version written to haeder of output file
    elsif ( $arg1 eq "-v" ) {
	if ( ! @ARGV ) {
	    &help_message;
	    die "Error: flag -v but no version as argument\n";
	}
	$version = shift(@ARGV);
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
	if ( /^# for make_html_reference:\s+(.+)/ ) {
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
    $newname =~ s/\.\w+$//;
    $fullname{$newname} = $file;
    $file = $newname;
}


# sort and shorten @files 
@files = sort @files;
@files = &shorten_list(@files);





# extract names of public functions and subroutines from doc-files
foreach $file ( @files ) {
    # warn "processing $file\n";
    local ( $docfile ) = $fullname{$file};
    $docfile =~ s/\.\w+$/.doc/;
    open(FILE,$docfile) || die "Error: $docfile not found: $!\n";
    local ( $subroutinesection ) = 0;
    local ( $typesection ) = 0;
    local ( $interfacesection ) = 0;
    local ( $continueline ) = 0;
    local ( $string ) = "";
    local ( $linenumber ) = 0;
    local( $startlinenumber );
    while (<FILE>) {
	$linenumber++;
	if ( $subroutinesection ) {
	    if ( /^\s*[sS][uU][bB][rR][oO][uU][tT][iI][nN][eE]\s+(\w+)/ ) {
		local($name) = $1;
		$name =~ tr/A-Z/a-z/;
		local($finalname) = "$name,$file,$linenumber";
		push(@names,$finalname);
		$types{$finalname} = "module subroutine";
	    }
	    elsif ( /^[^!]*[fF][uU][nN][cC][tT][iI][oO][nN]\s+(\w+)/ ) {
		local($name) = $1;
		$name =~ tr/A-Z/a-z/;
		local($finalname) = "$name,$file,$linenumber";
		push(@names,$finalname);
		$types{$finalname} = "module function";
	    }
	}
	elsif ( /^=============== Public functions and subroutines: ===============/ ) {
	    $subroutinesection = 1;
	}
	else {
            # concatinate continue lines in $string
	    if ( /^\s*&(\s*([^\s&!]+\s*)*)&/ ) {
		if ( $continueline ) {
		    $string = $string . $1;
		}
		else {
		    die "input file syntax error with continuation line\n";
		}
	    }
	    elsif ( /^(\s*([^\s&!]+\s*)+)&/ ) {
		if ( $continueline ) {
		    $string = $string . $1;
		}
		else {
		    $string =  $1;
		    $continueline = 1;
		    $startlinenumber = $linenumber;
		}
	    }
	    elsif ( /^\s*&(\s*([^\s&!]+\s*)*)/ ) {
		if ( $continueline ) {
		    $string = $string . $1;
		    $continueline = 0;
		}
		else {
		    die "input file syntax error with continuation line\n";
		}
	    }
	    else {
		if ( $continueline ) {
		    $string = $string . $_;
		}
		else {
		    $string =  $_;
		    $startlinenumber = $linenumber;
		}
		$continueline = 0;
	    }
	    if ( ! $continueline ) {
		if ( $typesection ) {
		    if ( /^\s*[eE][nN][dD]\s*[tT][yY][pP][eE]/ ) {
			$typesection = 0;
		    }
		}
		elsif ( $interfacesection ) {
		    if ( /^\s*[eE][nN][dD]\s*[iI][nN][tT][eE][rR][fF][aA][cC][eE]/ ) {
			$interfacesection = 0;
		    }
		}
		else {
		    # type declaration
		    if ( $string =~ /^\s*[tT][yY][pP][eE]\s*(,\s*\w+\s*)?(::)?\s*(\w+)\s*(!.*)?$/ ) {
			local($name) = $3;
			$name =~ tr/A-Z/a-z/;
			local($finalname) = "$name,$file,$startlinenumber";
			push(@names,$finalname);
			$types{$finalname} = "type";
			$typesection = 1;
		    }
		    # interface declaration with name
		    elsif ( $string =~ /^\s*[iI][nN][tT][eE][rR][fF][aA][cC][eE]\s+(\S+)\s*(!.*)?$/ ) {
			local($name) = $1;
			$name =~ tr/A-Z/a-z/;
			local($finalname) = "$name,$file,$startlinenumber";
			push(@names,$finalname);
			$types{$finalname} = "interface";
			$interfacesection = 1;
		    }
		    # interface declaration without name
		    elsif ( $string =~ /^\s*[iI][nN][tT][eE][rR][fF][aA][cC][eE]\s*(!.*)?$/ ) {
			$interfacesection = 1;
		    }
		    # variable declaration
		    elsif ( $string =~ /^\s*(\w+(\s*\([^\)]*\))?(\s*,\s*\w+)*)(\s+|(\s*::))\s*(\w+(\s*\([^\)]*\))?(\s*=\s*[^\s,]+)?(\s*,\s*(\w+(\s*\([^\)]*\))?(\s*=\s*[^\s,]+)?))*)\s*(!.*)?$/ ) {
			local ($type) = $1;
			local($vars) = $6;
			$type =~ s/\s*,\s*/, /g;
			$type =~ tr/A-Z/a-z/;
			$vars =~ tr/A-Z/a-z/;
			if ( ! ($type =~ /^(public|private|module|use)/) ) {
			    if ( $type =~ /^(subroutine|function)/ ) {
				local($name) = $vars;
				if ( $name =~ /(\w+)\s*(\([^\)]*\))/ ) {
				    local($finalname) = "$1,$file,$startlinenumber";
				    push(@names,$finalname);
				    $types{$finalname} = "$type";
				}
				else {
				    local($finalname) = "$name,$file,$startlinenumber";
				    push(@names,$finalname);
				    $types{$finalname} = $type;
				}
			    }
			    else {
				local($nbr) = 1;
				while ( $vars ne "" ) {
				    if ( $nbr > 100 ) {
					die "Error: can not split $vars\n";
				    }
				    $vars =~ s/(\w+(\s*\([^\)]*\))?(\s*=\s*[^\s,]+)?)(\s*,\s*(\w+(\s*\([^\)]*\))?(\s*=\s*[^\s,]+)?(\s*,\s*\w+(\s*\([^\)]*\))?(\s*=\s*[^\s,]+)?)*))?/$5/;
				    local($name) = $1;
				    if ( $name =~ /^(\w+)\s*(\([^\)]*\))/ ) {
					local($finalname) = "$1,$file,$startlinenumber";
					push(@names,$finalname);
					$types{$finalname} = "$type, dimension$2";
				    }
				    elsif ( $name =~ /^(\w+)\s*=[^,]+/ ) {
					local($finalname) = "$1,$file,$startlinenumber";
					push(@names,$finalname);
					$types{$finalname} = $type;
				    }
				    else {
					local($finalname) = "$name,$file,$startlinenumber";
					push(@names,$finalname);
					$types{$finalname} = $type;
				    }
				    $nbr++;
				}
			    }
			}
		    }
		}
	    }
	}
    }
    close (FILE);
}


# sort  @names
@names = sort( @names );


# open output file $outfile
open(OUTFILE,"> $outfile") || die "Error: $outfile can not be generated: $!\n";

# print HTML haeder
print OUTFILE "<HTML>\n";
print OUTFILE "\n";
print OUTFILE "<HEAD>\n";
print OUTFILE "<TITLE>ParaGauss Programing Reference</TITLE>\n";
print OUTFILE "</HEAD>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<BODY BGCOLOR=\"#ffffff\">\n";
print OUTFILE "<FONT SIZE=3>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<B><FONT SIZE=8>\n";
print OUTFILE "<P ALIGN=\"CENTER\">\n";
print OUTFILE "ParaGauss Programing Reference\n";
print OUTFILE "</P>\n";
print OUTFILE "</FONT></B>\n";
print OUTFILE "<BR>\n";
print OUTFILE "<FONT SIZE=4><I>\n";
print OUTFILE "<P ALIGN=\"CENTER\">\n";
print OUTFILE "Th Belling, U Birkenheuer, S Kr&uuml;ger, M Mayer, F N&ouml;rtemann,\n";
print OUTFILE "M Staufer, N R&ouml;sch\n";
print OUTFILE "</P>\n";
print OUTFILE "<BR>\n";
print OUTFILE "<P ALIGN=\"CENTER\">\n";
print OUTFILE "ParaGauss Version $version\n";
print OUTFILE "</P>\n";
print OUTFILE "<BR>\n";
print OUTFILE "<P ALIGN=\"CENTER\">\n";
print OUTFILE "This file was automatically created by make_html_reference\n";
print OUTFILE "<BR>\n";
require "ctime.pl";
print OUTFILE "at ",&ctime(time);
print OUTFILE "</P>\n";
print OUTFILE "</FONT></I>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<BR><BR><HR SIZE=4 NOSHADE><BR><BR>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<H1>\n";
print OUTFILE "<A NAME=\"Contents\"></A>\n";
print OUTFILE "<B><FONT SIZE=7>Contents</FONT></B>\n";
print OUTFILE "</H1>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<OL type=\"A\"><FONT SIZE=4>\n";
print OUTFILE "<LI><A HREF=\"#List_of_Public_Names\">List of Public Names</A>\n";
print OUTFILE "<P>\n";

# print reference for variables beginning with char
foreach $char ( a .. z ) {
    print OUTFILE "<A HREF=\"#Variables_$char\">$char</A>\n";
}

print OUTFILE "</P>\n";
print OUTFILE "<LI><A HREF=\"#Haeders_of_Program_Files\">Haeders of Program Files</A>\n";
print OUTFILE "  <OL><FONT SIZE=3>\n";

# print @files (names only) in alphatic order for contents
foreach $file ( @files ) {
   print OUTFILE "  <LI><A HREF=\"#Files_$file\">$file</A>\n";
}

print OUTFILE "  </FONT></OL>\n";
print OUTFILE "</FONT></OL>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<BR><BR><HR SIZE=4 NOSHADE><BR><BR>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<H1>\n";
print OUTFILE "<A NAME=\"List_of_Public_Names\"></A>\n";
print OUTFILE "<B><FONT SIZE=7></FONT>List of Public Names</B>\n";
print OUTFILE "</H1>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<P>\n";

# prepare character array
@characters = ( a .. z );
# print reference for variables beginning with char
foreach $char ( @characters ) {
    print OUTFILE "<A HREF=\"#Variables_$char\">$char</A>\n";
}
unshift(@characters," ");

print OUTFILE "</P>\n";
print OUTFILE "\n";

$i = 0;
$tableopen = 0;
# print names, types and locations of variables etc
foreach $finalname ( @names ) {
    $finalname =~ /([^,]+),([^,]+),([^,]+)/;
    local($name) = $1;
    local($location) = $2;
    local($linenumber) = $3;
    $char = substr($name,0,1);
    $newtable = 0;
    while ( $char ne $characters[$i] && $i < $#characters ) {
	$i++;
	if ( $tableopen == 1 ) {
	    print OUTFILE "</TABLE>\n\n";
	    print OUTFILE "<P><I>\n";
	    print OUTFILE "Go to: <A HREF=\"#Haeders_of_Program_Files\">Haeders of Program Files</A>\n";
	    print OUTFILE "<A HREF=\"#List_of_Public_Names\">List of Public Names</A>\n";
	    print OUTFILE "</I></P><BR>\n\n";
	    $tableopen = 0;
	}
	print OUTFILE "<H3><B>\n<A NAME=\"Variables_$characters[$i]\">$characters[$i]</A>\n</B></H3>\n";
	$newtable = 1;
    }
    if ( $newtable == 1 ) {
	print OUTFILE "<TABLE BORDER=4>\n";
	print OUTFILE "<TR><TH>Name<TH>Type<TH>Contained in</TR>\n";
	$tableopen = 1;
    }
    print OUTFILE "<TR><TD><A HREF=\"#Files_$location","_$linenumber\">$name</A>";
    print OUTFILE "<TD>$types{$finalname}<TD><A HREF=\"#Files_$location\">$location</A></TR>\n\n";
}

if ( $tableopen == 1 ) {
    print OUTFILE "</TABLE>\n\n";
    print OUTFILE "<P><I>\n";
    print OUTFILE "Go to: <A HREF=\"#Haeders_of_Program_Files\">Haeders of Program Files</A>\n";
    print OUTFILE "<A HREF=\"#List_of_Public_Names\">List of Public Names</A>\n";
    print OUTFILE "</I></P>\n\n";
    $tableopen = 0;
}

print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<BR><BR><HR SIZE=4 NOSHADE><BR><BR>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<H1>\n";
print OUTFILE "<A NAME=\"Haeders_of_Program_Files\"></A>\n";
print OUTFILE "<B><FONT SIZE=7></FONT>Haeders of Program Files</B>\n";
print OUTFILE "</H1>\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "<P>\n";

# print @files (names only) in alphatic order for contents
foreach $file ( @files ) {
   print OUTFILE "<A HREF=\"#Files_$file\">$file</A>\n";
}

print OUTFILE "</P>\n";
print OUTFILE "\n";
print OUTFILE "\n";


# print @files in alphatic order to $outfile
foreach $file ( @files ) {
    print OUTFILE "<BR><BR><HR SIZE=2 NOSHADE><BR><BR>\n";
    print OUTFILE "\n";
    print OUTFILE "<H2>\n";
    print OUTFILE "<A NAME=\"Files_$file\"></A>\n";
    print OUTFILE "<B><FONT SIZE=5>$file</FONT></B>\n";
    print OUTFILE "</H2>\n";
    print OUTFILE "\n";
    print OUTFILE "<P>\n";
    print OUTFILE "Haeder contained in file \"$fullname{$file}\"\n";
    print OUTFILE "</P>\n";
    print OUTFILE "<BR>\n";
    print OUTFILE "<PRE>\n";

    local ( $docfile ) = $fullname{$file};
    $docfile =~ s/\.\w+$/.doc/;
    open(FILE,$docfile) || die "Error: $docfile not found: $!\n";
    local($linenumber) = 0;
    while (<FILE>) {
	$linenumber++;
	print OUTFILE "<A NAME=\"Files_$file","_$linenumber\"></A>";
	s/>/&gt\;/g;
	s/</&lt\;/g;
	s/&/&amp\;/g;
	s/"/&quot\;/g;
	print OUTFILE;
    }
    close (FILE);

    print OUTFILE "</PRE>\n";
    print OUTFILE "\n";
    print OUTFILE "<P><I>\n";
    print OUTFILE "Go to: <A HREF=\"#Haeders_of_Program_Files\">Haeders of Program Files</A>\n";
    print OUTFILE "<A HREF=\"#List_of_Public_Names\">List of Public Names</A>\n";
    print OUTFILE "</I></P>\n";
    print OUTFILE "\n";
    print OUTFILE "\n";
    print OUTFILE "\n";
}
print OUTFILE "\n\n\n";


print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "\n";
print OUTFILE "</FONT>\n";
print OUTFILE "</BODY>\n";
print OUTFILE "\n";
print OUTFILE "</HTML>\n";
close(OUTFILE);

