
ifndef versions
date = $(shell date +%Y%m%d)
versions = V3.0.40b13-32 V3.1cvs$(date)-32
endif

ifndef SUBDIRS
SUBDIRS = $(shell cat manifest)
endif

all:
	export versions="$(versions)"; \
	for dir in $(SUBDIRS); do \
		echo === Descending into $$dir ===; \
		cd $$dir; \
		../bin/script $(runopts) `cat manifest`; \
		cd ..; \
	done

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
