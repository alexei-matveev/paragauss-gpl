== Testsuite ==

The  testsuite  is  used  to  detect  regressions  in  the  developing
program. Passing a testsuite is a prereuisite for a "public" release.

To  get  the  sources  (inputs  and  scripts)  of  the  testsuite  use
[http://www.darcs.net Darcs] (available on theo2):

  darcs get ~matveev/darcs/testsuite
  cd testsuite/qa
  chmod +x bin/*

The <tt>chmod</tt> commad  is needed to set the  executable bit on the
scripts after initial download.

For update of the already downloaded testsuite issue:

  cd testsuite
  darcs pull ~matveev/darcs/testsuite

To compare two ParaGauss versions

* either edit the <tt>bin/script</tt> or set an evironment variable

  setenv versions "V3.0 V3.1" # tcsh syntax
  export versions="V3.0 V3.1" # bash syntax

* either edit the <tt>bin/hostfile</tt> or set an evironment variable

  setenv PE_HOSTFILE /some/other/location/hostfile # tcsh syntax
  export PE_HOSTFILE=/some/other/location/hostfile # bash syntax

* and run the testsuite by e.g.

  cd pd-ae-pp-sr
  ../bin/script i.*

In a (rare) case of success you will see only something like this:

  Running i.pd2    with V3.0 ... done. with V3.1 ... done. Compare to V3.0 -- ok
  Running i.pd2,ae         with V3.0 ... done. with V3.1 ... done. Compare to V3.0 -- ok
  Running i.pd2,ae,sr      with V3.0 ... done. with V3.1 ... done. Compare to V3.0 -- ok
  Running i.pd2,ap         with V3.0 ... done. with V3.1 ... done. Compare to V3.0 -- ok
  Running i.pd2,ap,sr      with V3.0 ... done. with V3.1 ... done. Compare to V3.0 -- ok
  Running i.pd2,pp         with V3.0 ... done. with V3.1 ... done. Compare to V3.0 -- ok

In case  there are  differences (remember that  real numbers  could be
compare only  up to  some precision)  you will see  the output  of the
<tt>diff</tt> command:

  ...
  Running i.pd2,pp         with V3.0 ... done. with V3.1 ... done. Compare to V3.03,4c3,4
  <   e_sum  =          -253.284047  (VWN)
  <   e_sum  =          -253.284047  (scf)
  ---
  >   e_sum  =          -253.284051  (VWN)
  >   e_sum  =          -253.284051  (scf)
  
  FAILED test pd2,pp at precision 6!
  FAILED test pd2,pp@6
  FAILED 6 tests at precision 6: pd2@6 pd2,ae@6 pd2,ae,sr@6 pd2,ap@6 pd2,ap,sr@6 pd2,pp@6
  FAILED 6 of 6 tests!

and the summary  at the end saying which tests  and at which precision
have  failed.  The  test fails  by default  if any  two  numbers still
differ after truncating them to 8 decimal figures.

You can set the required precisions for comparisons by

  setenv precisions "5 3" # tcsh syntax
  export precisions="5 3" # bash syntax

default is <tt>precisions="8 6 4"</tt>.

Currently, inputs are available in the following subdirectories:

  qa/au-tddft
  qa/uf6-geo-frq
  qa/so-ee-sr
  qa/uo2-geo-frq
  qa/h2o-sol-pdc
  qa/uf6-solv-frq
  qa/h2o-solv-frq
  qa/h3o-solv-frq
  qa/occh3ch3-solv-frq
  qa/ochh-solv-frq
  qa/ch3oh
  qa/pd-ae-pp-sr
  qa/ch4-geo-frq
  qa/nico-dftplus-u
  qa/h2o2-various-options
  qa/h2o2-functionals

Please populate the testsuite with more inputs to cover all corners of
the ParaGauss functionality.

