#
#    \033          ascii ESCape
#    \033[<NUM>G   move to column <NUM> (linux console, xterm, not vt100)
#    \033[<NUM>C   move <NUM> columns forward but only upto last column
#    \033[<NUM>D   move <NUM> columns backward but only upto first column
#    \033[<NUM>A   move <NUM> rows up
#    \033[<NUM>B   move <NUM> rows down
#    \033[1m       switch on bold
#    \033[31m      switch on red
#    \033[32m      switch on green
#    \033[33m      switch on yellow
#    \033[m        switch off color/bold
#    \017          exit alternate mode (xterm, vt100, linux console)
#    \033[10m      exit alternate mode (linux console)
#    \015          carriage return (without newline)
#

if test -t 1 -a "$TERM" != "raw" -a "$TERM" != "dumb" && stty size > /dev/null 2>&1 ; then
	rc_tab="\033[30G"
	rc_ok="\033[32mok\033[m\017"
	rc_no="\033[31mNO\033[m\017"
	rc_failed="\033[31mFAILED\033[m\017"
else
	rc_tab="\t"
	rc_ok="ok"
	rc_no="NO"
	rc_failed="FAILED"
fi

make_summary(){
	local out
	out=$1
	# extract energies and gradients:
	#echo "========================= $out ============================="
	echo "GEOMETRY:"
	geom.pl $out/output
	echo "ENERGY:"
	Esum  $out
	echo "GRADIENTS:"
	grads $out
	echo "FREQUENCIES:"
	# frequencies from cartesian analytica hessian go to output:
	freq.pl $out/output
	# freqs from numerical differentiation go to flepo files:
	if [ -f $out/optimizer/flepo.1 ]; then
	freq.pl $out/optimizer/flepo.*
	fi
	# with newer version flepo files are located in output directory:
	if [ -f $out/flepo.1 ]; then
	freq.pl $out/flepo.*
	fi
	echo "POTENTIAL-DERIVED CHARGES:"
	pdc.pl $out/output
        echo "TDDFT EXCITATION ENERGIES AND OSCILLATOR STRENGTHS:"
        tddft.pl $out/output  
}

diff_prec(){
	local sum0
	local sum1
	local prec
	local retval
	sum0=$1
	sum1=$2
	prec=$3
	if [ -z "$prec" ]; then
		prec=6
	fi
	prec.pl -$prec $sum0 > tmp$sum0
	prec.pl -$prec $sum1 > tmp$sum1
	diff tmp$sum0 tmp$sum1
	retval=$?
	rm -f tmp$sum0 tmp$sum1
	return $retval
}

esum() {
	local out
	local cat
	if [ -d $1 ] ; then
		out=$1/output
	else
		out=$1
	fi
	if test ! -f $out -a -f $out.gz ; then
		out=$out.gz
		cat=zcat
	else
		cat=cat
	fi
	$cat $out | grep e_sum | tail -40
}
e_xc() {
	local out
	local cat
	if [ -d $1 ] ; then
		out=$1/output
	else
		out=$1
	fi
	if test ! -f $out -a -f $out.gz ; then
		out=$out.gz
		cat=zcat
	else
		cat=cat
	fi
	$cat $out | grep e_xc | tail -30
}

grads() {
	local out
	local cat
	if [ -d $1 ] ; then
		out=$1/output
	else
		out=$1
	fi
	if test ! -f $out -a -f $out.gz ; then
		out=$out.gz
		cat=zcat
	else
		cat=cat
	fi
	$cat $out | grep 'Equal C'
}

sders() {
	local out
	local cat
	if [ -d $1 ] ; then
		out=$1/output
	else
		out=$1
	fi
	if test ! -f $out -a -f $out.gz ; then
		out=$out.gz
		cat=zcat
	else
		cat=cat
	fi
	$cat $out | grep '^SD:'
}

function Esum() {
	esum $* | fgrep ' [' | tail -1
	esum $* | fgrep -v FINAL | fgrep '('
}

function efinal() {
	esum $* | fgrep FINAL
}

ewc() {
	grep e_sum $1 | fgrep ' [' | wc -l
}
