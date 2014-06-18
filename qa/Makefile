
ifndef versions
versions = V3.1pre6-32 V3.1cvs20060905-32 V3.1am12-32
endif


SUBDIRS = 	$(shell cat manifest)
#SUBDIRS = 	pd-ae-pp-sr h2o-sol-pdc ch4-geo-frq

all:
	export versions="$(versions)"; \
	for dir in $(SUBDIRS); do \
		echo === Descending into $$dir ===; \
		cd $$dir; \
		../bin/script `cat manifest`; \
		cd ..; \
	done
