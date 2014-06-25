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

# makedependent TB 10/95

###########  subroutines

sub help_message {
    print  <<"ENDMESSAGE";

makedependent TB 10/95

usage :   help:           makedependent -h
          makefile mode:  makedependent [-d|-n] [-imd] [<filename>]
          listfile mode:  makedependent [-d|-n] -l [<filename>]


Produces description of fortran 90 file dependencies in makefile format.
Sorts the object files according to their dependencies: 
  files depending on others are put before them.
Detects cross-dependencies between modules.
Looks in the corresponding fortran 90 source files for use statements and
determines dependencies between the files in this way.

It is assumed that the modules are located in source files named according
to the module names with minor letters only.

The flag -d makes makedepent write a file *.dep listing the
dependencies *.mod for each object file *.o that can be used
by during the make procedure to minimize necessary work

The flag -n surpresses writing of dependencies in Makefile Format.
Only check for cross-dependencies is performed and *.dep files
are generated if Flag -d is given as well.

The flag -imd forces the same dependencies to be written
in the .dep files as those in Makefile.
i.e. only immediate dependencies



Two modes of operation exist:

Listfile mode:
Invoke with "makedependent -l" or "makedependent -l <filename>".
Reads in file "objlist.makedependent" or <filename> that should contain 
the name of one object file in each line.
Prints out the object files and their dependencies to STDOUT.

Makefile mode:
Invoke with "makedependent" or "makedependent <filename>".
Edits file "Makefile" or <filename> and saves original file in <filename>.old.
The file must include the following three lines (beginning with #) :
  # for makedependent: <objectmakroname> [ <objectmakroname> ... ]
  # makedependent start insert
  # makedependent end insert
Makedependent reads the object files to process from the objectmakro definitions
and inserts the necessary dependency descriptions between the start and end line,
removing anything old between these lines.

to obtain help statement, type "makedependent -h".

ENDMESSAGE
}

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

sub rm_f90_unix_list {
	local(@newlist);
        local($elem);
	foreach $elem (@_) {
		if ($elem ne "f90_unix") {
			push(@newlist,$elem);
		}
	}
	@newlist;
}

sub to_minuscles_list {
	local(@newlist);
        local($elem);
	foreach $elem (@_) {
		$elem =~ tr/A-Z/a-z/;
		push(@newlist,$elem);
	}
	@newlist;
}

sub add_o_list {
	local(@newlist);
        local($elem);
	foreach $elem (@_) {
		push(@newlist,$elem . ".o");
	}
	@newlist;
}

sub full_name {
	local($string,@array) = @_;
	local($i);
	SEARCHLOOP: for($i=0; $i < @array; $i++) {
	    if ($array[$i] =~ /\/$string/) { last SEARCHLOOP; }
	    if ($array[$i] =~ /^$string/) { last SEARCHLOOP; }
	}
	$array[$i];
}

sub list_to_string {
	local($string);
        local($elem);
	foreach $elem (@_) {
		$string .= $elem;
		$string .= " ";
	}
	$string =~ s/ $//;
	$string;
}

sub string_to_list {
	local($string) = @_;
	local(@newlist);
	foreach $elem ( split(/ /,$string) ) {
	    push(@newlist,$elem);
	}
	@newlist;
}

sub position {
	local($string,@array) = @_;
	local($i);
	for($i=0; $array[$i] ne $string && $i < @array; $i++) {}
	$i;
}

sub change_position {
	local($string1,$string2,@array) = @_;
	local($i1) = &position($string1,@array);
	local($i2) = &position($string2,@array);
	local($string3) = $array[$i1];
	$array[$i1] = $array[$i2];
	$array[$i2] = $string3;
	@array;
}









########### executable part

select(STDOUT); $|=1;



# read command line arguments
$immediate_dep=0;
$makefilemode = 1;
$listmode     = 0;
$objlistfile = "objlist.makedependent";
$makefile = "Makefile";
$makefiletmp = "Makefile.tmp";
$writedependencefiles = 0;
$nomakefile = 0;
while ( @ARGV ) {
    $arg1 = shift(@ARGV);
    # -h invokes help message
    if ( $arg1 eq "-h" ) {
	&help_message;
	exit;
    }
    # -l invokes list mode
    elsif ( $arg1 eq "-l" ) {
	$makefilemode = 0;
	$listmode     = 1;
    }
    # -d makes makedependent write a list file with
    # *.dep with its dependencies for each object file *.o
    elsif ( $arg1 eq "-d" ) {
	$writedependencefiles = 1;
    }
    # no Makefile output
    elsif ( $arg1 eq "-n" ) {
	$nomakefile = 1;
    }
    # only immediate dependencies to .dep file
    elsif ( $arg1 eq "-imd" ) {
	$immediate_dep=1;
    }
    # read other arguments as filenames
    else {
	if ( $makefilemode ) {
	    $makefile = $arg1;
	}
	else {
	    $objlistfile = $arg1;
	}
    }
}

$verbose = 1;
if( $listmode ){ $verbose = 1; }

if( $writedependencefiles ){
	if( $immediate_dep){
		print STDERR "only immediate dependencies to .dep files\n" if $verbose;
	}else{
		print STDERR "all dependencies to .dep files\n" if $verbose;
	}
}



# check if file "makefile.tmp" exists, if so die
if ( $makefilemode ) {		# makefile mode
    if ( -e $makefiletmp ) {
	die "makedependent: Error: file $makefiletmp exists. This file would be overwritten.";
    }
}



# read list of object files @objfiles
if ( $makefilemode ) {		# makefile mode
    local(@objmakros);

    # print "reading objecf files from Makefile $makefile\n";

    # read line with list of object makros
    open(MAKEFILE,$makefile) || die "makedependent: Error: $makefile not found: $!\n";
    while(<MAKEFILE>) {
	if ( /^# for makedependent:\s+(.+)/ ) {
	    @objmakros = split(/\s+/,$1);
	}
    }
    close(MAKEFILE);
    if ( ! @objmakros ) {
	die "makedependent: Error: $makefile does not contain line listing object makros \n";
    }

    # read in object files from makro definitions
    foreach $objmakro ( @objmakros ) {
	open(MAKEFILE,$makefile) || die "makedependent: Error: $makefile not found: $!\n";
	local( $found ) = 0;
	REODOBYMAKRO: while(<MAKEFILE>) {
	    if ( /^$objmakro\s*=\s+([^\\]*)$/ ) {
		splice(@objfiles,0,0,split(/\s+/,$1));
		$found = 1;
		last REODOBYMAKRO;
	    }
	    if ( /^$objmakro\s*=\s+([^\\]*)[\\]/ ) {
		splice(@objfiles,0,0,split(/\s+/,$1));
		$found = 1;
	    }
	    else {
		if ( $found ) {
		    if ( /^\s+([^\\]*)$/ ) {
			splice(@objfiles,0,0,split(/\s+/,$1));
			last REODOBYMAKRO;
		    }
		    if ( /^\s+([^\\]*)[\\]/ ) {
			splice(@objfiles,0,0,split(/\s+/,$1));
		    }
		    else {
			last REODOBYMAKRO;
		    }
		}
	    }
	}
	close(MAKEFILE);
	if ( ! $found ) {
	    die "makedependent: Error: object makro definition of $objmakro not found in $makefile\n";
	}
    }

}

elsif( $listmode ) {				# list mode
    # print "reading objecf files from list file $objlistfile\n";
    open(OBJLIST,$objlistfile) || die "makedependent: Error: $objlistfile not found: $!\n";
    while(<OBJLIST>) {
	chop;
	# push(@objfiles,$_);
	splice(@objfiles,0,0,split(/\s+/,$_));
    }
    close(OBJLIST);
} else {
	die "never should have happend\n";
}


# sort and shorten @objfiles to eliminate double references
@objfiles = sort @objfiles;
@objfiles = &shorten_list(@objfiles);


# initialise @used that stores information if objects used by others
foreach $objfile (@objfiles) {
    $used{$objfile} = 0;
}




# search for use statements in corresponding fortran files
# and store them in $dependencies{$objfile}.
# store information if objects are used by others in $used{$objfile}.
foreach $objfile (@objfiles) {

	# array to store dependencies of @objfile
	local(@dep);

	# file names of fortran sources
	$f90file = $objfile;
        $f90file =~ s/\.o$/\.f90/;
	$ffile = $objfile;
        $ffile =~ s/\.o$/\.f/;
	$f90pfile = $objfile;
        $f90pfile =~ s/\.o$/\.f90p/;
	$fpfile = $objfile;
        $fpfile =~ s/\.o$/\.fp/;

	# search for module names in use statements of fortran sources
	$foundfile = open(FFILE,$f90pfile);
	if ( $foundfile) {
	    $sourcefile{$objfile} = $f90pfile;
	}
	else {
	    $foundfile = open(FFILE,$f90file);
	    if ( $foundfile) {
		$sourcefile{$objfile} = $f90file;
	    }
	    else {
		$foundfile = open(FFILE,$fpfile);
		if ( $foundfile ) {
		    $sourcefile{$objfile} = $fpfile;
		}
		else {
		    open(FFILE,$ffile) || die "makedependent: Error no fortran file for $objfile, tried: $f90pfile $f90file $fpfile $ffile \n";
		    $sourcefile{$objfile} = $ffile;
		}
	    }
	}
	if( $verbose ){
	print STDERR "searching USE in $sourcefile{$objfile}                                      \r";
	}

	while(<FFILE>) {
        	push(@dep,$1) if /^\s*use\s+(\w+)/i;
	}
	close(FFILE);


	# module names to minuscles
	@dep = &to_minuscles_list(@dep);

	# sort @dep list
	@dep = sort @dep;

	# eliminate double occurences of elements in @dep list
	@dep = &shorten_list(@dep);

	# remove f90_unix module from @dep list
	@dep = &rm_f90_unix_list(@dep);

	# add ".o" to module names to obtain names of object files
	@dep = &add_o_list(@dep);

	# look for full names including directories in $objfile
	foreach $file ( @dep ) {
	    $file = &full_name($file,@objfiles);
	}

	# register for each file that it is used
	foreach $file ( @dep ) {
	    $used{$file} = 1;
	}

	# write @dep list to blank-seperated string and 
	# store the string in $dependencies{$objfile}
	$dependencies{$objfile} = &list_to_string(@dep);


        # check if file .flags exist  and store in $flagsfile{$objfile}
        # will be added to dependency list later on.
	$flfile = $objfile;
        $flfile =~ s/\.o$/\.flags/;
	$foundfile = open(FLFILE,$flfile);
	if ( $foundfile) {
	    $flagsfile{$objfile} = $flfile;
	    close(FLFILE);
	}
	else {
	    $flagsfile{$objfile} = "";
	}
}
	if( $verbose ){
	print STDERR "search done, checking dependencies                      \r";
	}

# generate a list of full recursive dependencies for each obj
# and check for cross dependencies at the same time

# initialising list with zero strings
foreach $objfile (@objfiles) {
    $fulldep{$objfile} = "";
}

sub check_fulldependency {
    local($intermediate_obj,@predecessor_objs) = @_;
    local($elem);
    if ( $dependencies{$intermediate_obj} ne "" ) {
	if ( $fulldep{$intermediate_obj} ne "" ) { # use $fulldep{$intermediate_obj} for check
	    foreach $elem ( split(/ /,$fulldep{$intermediate_obj}) ) {
		foreach $pred ( @predecessor_objs ) {
		    if ( $elem eq $pred ) {
			die "makedependent: Error: cross-dependency: $elem uses $pred via @predecessor_objs \n";
		    }
		}
	    }
	}
	else {			# build up $fulldep{$intermediate_obj} and check at one sweep
	    local($string) = "";
	    foreach $elem ( split(/ /,$dependencies{$intermediate_obj}) ) {
		foreach $pred ( @predecessor_objs ) {
		    if ( $elem eq $pred ) {
			die "makedependent: Error: cross-dependency: $elem uses $pred via @predecessor_objs \n";
		    }
		}
		$string .= &check_fulldependency($elem,(@predecessor_objs,$elem));
		$string .= " ";
	    }
	    $string .= $dependencies{$intermediate_obj};
	    $string =~ s/ $//;
            # sort and shorten dependency list
	    local(@list) = &string_to_list($string);
	    @list = sort @list;
	    @list = &shorten_list(@list);
	    $fulldep{$intermediate_obj} = &list_to_string(@list);
	}
    }
    $fulldep{$intermediate_obj};
}


# print "checking for cross-dependencies\n";
foreach $objfile (@objfiles) {
    if ( $fulldep{$objfile} eq "" ) {
	$fulldep{$objfile} = &check_fulldependency($objfile,($objfile));
    }
}

# print "print files containing all dependendencies\n";
if ( $writedependencefiles && ! $immediate_dep ) {
    foreach $objfile (@objfiles) {
	local($depfile) = $objfile;
	$depfile =~ s/\.o/.dep/g;
	if( $verbose ){
	print STDERR "writing full deps in $depfile                           \r";
	}
	@completedep = &string_to_list($fulldep{$objfile});
	open(DEPFILE,"> $depfile") || die "makedependent: Error: $depfile can not be generated: $!\n";
	if ( $#completedep ) {
	    local($elem);
	    foreach $elem ( @completedep ) {
		if ( $elem ne "" ) {
		    $elem =~ s/\.o/.mod/g;
		    print DEPFILE "$elem\n";
		}
	    }
	}
	close(DEPFILE);
    }
}			      


# stop here if no Makile output is requested
if ( $nomakefile ) {
    print "makedependent: No cross dependencies detected\n";
    if ( $makefilemode ) {
	print "makedependent: Stopping without modifying $makefile due to Flag -n\n";
    }
    exit;
}


# sort object files
$finished = 0;
while ( ! $finished) {
    $finished = 1;
  OBJ: for($j=0;  $j < @objfiles; $j++) {
      $objfile = @objfiles[$j];
      for($i = 0; $i < $j; $i ++) {
	  $objfile2 = @objfiles[$i];
	  if ( $fulldep{$objfile} =~ /$objfile2/ ) {
	      @objfiles = &change_position($objfile,$objfile2,@objfiles);
	      $finished = 0;
	      next OBJ;
	  }
      }
  }
}


# print out object files and dependencies
if ( $makefilemode ) {		# makefilemode
    local( $found ) = 0;
    open(MAKEFILE,"$makefile") || die "makedependent: Error: $makefile not found: $!\n";
    open(MAKEFILETMP,"> $makefiletmp") || die "makedependent: Error: $makefiletmp can not be generated: $!\n";
    while(<MAKEFILE>) {
	if ( /^# makedependent end insert/ ) {
	    $found = 0;
	}
	if ( ! $found ) {
	    print MAKEFILETMP $_;
	}
	if ( /^# makedependent start insert/ ) {
	    $found = 1;
	    foreach $objfile (@objfiles) {
		local($dep) = $dependencies{$objfile};
		$dep =~ s/\.o/.mod/g;
		if ( $used{$objfile} ) {
		    local($modtarget) = $objfile;
		    $modtarget =~ s/\.o/.mod/g;
		    print MAKEFILETMP "$objfile $modtarget",":\t","$dep $sourcefile{$objfile} $flagsfile{$objfile}\n";
		}
		else {
		    print MAKEFILETMP "$objfile",":\t","$dep $sourcefile{$objfile} $flagsfile{$objfile}\n";
		}
		print MAKEFILETMP "";
		#am:
		if( $writedependencefiles && $immediate_dep ) {
			local($dep_file) = $objfile;
			$dep_file =~ s/\.o/.dep/g;
			print "writing deps in $dep_file                           \r";
			$dep =~ s/ /\n/g;
			open( DEP_FILE, "> $dep_file") || die "makedependent: cannot open $dep_file";
			print DEP_FILE "$dep\n";
			close(DEP_FILE);
		}
	    }
	}
    }
    close(MAKEFILE);
    close(MAKEFILETMP);
    rename($makefile,"$makefile.old") || warn "makedependent: can not move $makefile to $makefile.old: $!\n";
    rename($makefiletmp,$makefile) || die  "makedependent: can not move $makefiletmp to $makefile: $!\n";
    print "makedependent: $makefile was edited and old file saved as $makefile.old\n";
}

elsif( $listmode ) {				# list mode
    foreach $objfile (@objfiles) {
	local($dep) = $dependencies{$objfile};
	$dep =~ s/\.o/.mod/g;
	if ( $used{$objfile} ) {
	    local($modtarget) = $objfile;
	    $modtarget =~ s/\.o/.mod/g;
	    print "$objfile $modtarget",":\t","$dep $sourcefile{$objfile} \n";
	}
	else {
	    print "$objfile",":\t","$dep $sourcefile{$objfile} \n";
	}
		#am:
		if( $writedependencefiles && $immediate_dep ) {
			local($dep_file) = $objfile;
			$dep_file =~ s/\.o/.dep/g;
			print STDERR "writing deps in $dep_file                           \r";
			$dep =~ s/ /\n/g;
			open( DEP_FILE, "> $dep_file") || die "makedependent: cannot open $dep_file";
			print DEP_FILE "$dep\n";
			close(DEP_FILE);
		}
	print "";
    }
			print STDERR "                                                    \r";
}else{
	die "should have never happend\n";
}
