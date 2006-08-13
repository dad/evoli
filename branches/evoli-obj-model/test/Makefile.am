## Makefile.am

include $(top_srcdir)/config/Make-rules

MAINTAINERCLEANFILES += *.cutee.cc runtest.cc

## Cutee unit testing directives
# your test files here
test_files = t.genotype.h

noinst_PROGRAMS = cutee runtest

cutee_SOURCES = cutee.cc cutee.h

runtest_LDADD = $(top_builddir)/folder/libevolver.a

autocutee.mk: cutee Makefile.am $(test_files)
	$(CUTEE) -k -o autocutee.mk $(test_files)

include autocutee.mk