
SUBDIRS = clangO2  gcc48-O2  gccO0  gccO2

ALL = pointe5LETG.time  pointe5NONE.time  pointe5.time  point.time

.PHONY: subdirs $(SUBDIRS) ALL

%.time: ../../marx-test/speed/par/%.par install/bin/marx
	/usr/bin/time --output=$@ install/bin/marx @@$< > /dev/null 

subdirs: $(SUBDIRS)

$(SUBDIRS):
	cp Makefile $@/Makefile
	for timetest in $(ALL); do \
	  $(MAKE) -C $@ $$timetest; \
	done
