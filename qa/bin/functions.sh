
make_summary(){
	local out
	out=$1
	# extract energies and gradients:
	#echo "========================= $out ============================="
	echo "ENERGY:"
	Esum  $out
	echo "GRADIENTS:"
	grads $out
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
	$cat $out | grep e_sum | tail -30
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

egra() {
	esum $*
	grads $*
}

function Esum() {
	esum $* | fgrep -v FINAL | fgrep '('
}

function efinal() {
	esum $* | fgrep FINAL
}

ewc() {
	grep e_sum $1 | wc -l
}
