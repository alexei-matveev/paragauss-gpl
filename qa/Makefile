
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
		../bin/script `cat manifest`; \
		cd ..; \
	done
