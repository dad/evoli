## Makefile.am

include $(top_srcdir)/config/Make-rules

## Cutee unit testing directives
# your test files here
test_files = t.selftest.h 

noinst_PROGRAMS = cutee

cutee_SOURCES = cutee.cc cutee.h

autocutee.mk: cutee Makefile.am $(test_files)
	$(CUTEE) -k -o autocutee.mk $(test_files)

include autocutee.mk