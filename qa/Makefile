#
# FIXME: this isnt a Makefile, but a shell script ..
#

SUBDIRS ?= $(shell cat manifest)

all:
	for dir in $(SUBDIRS); do \
		echo === Descending into $$dir ===; \
		cd $$dir; \
		../bin/script $(runopts) `cat manifest`; \
		cd ..; \
	done

# will not work because "versions" are now in "local.conf":
summaries:
	for dir in $(SUBDIRS); do \
		echo === Descending into $$dir ===; \
		cd $$dir; \
	for ver in  $(versions); do \
		for inp in `cat manifest`; do \
			echo Summary of $$inp with $$ver ...; \
			nam=`echo $$inp | sed s/^i.//`; \
			out=o.$$nam,$$ver; \
			smr=s.$$nam,$$ver; \
			if [ -d $$out ]; then \
			../bin/mk_summary $$out > $$smr; \
			fi; \
		done; \
	done; \
		cd ..; \
	done

clean:
	for dir in $(SUBDIRS); do \
		echo === Descending into $$dir ===; \
		cd $$dir; \
		for entry in o.* s.*; do \
			echo deleting $$entry; \
			rm -rf $$entry; \
		done; \
		cd ..; \
	done
