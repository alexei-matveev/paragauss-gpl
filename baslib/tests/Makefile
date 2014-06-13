scripts = $(wildcard *.sh)
tests = $(scripts:.sh=.o)

all: $(tests) checksums

#
# cat A | diff B - is the same as diff B A:
#
%.o: %.sh
	./$(<) | diff $(<:.sh=.stdout) -

checksums:
	cat MD5SUMS | while read sum dash cmd; do echo $$($$cmd | md5sum) $$cmd; done | diff MD5SUMS -

#
# "make  update" rebuilds  *.stdout files,  in case  a  format changes
# significantly:
#
update: $(scripts:.sh=.update)

%.update: %.sh
	./$(<) > $(<:.sh=.stdout)
