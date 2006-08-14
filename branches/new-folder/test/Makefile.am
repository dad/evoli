## Makefile.am

include $(top_srcdir)/config/Make-rules
INCLUDES = -I$(top_srcdir)/src/folder/ -I$(top_srcdir)/src/evolver/

MAINTAINERCLEANFILES += *.cutee.cc runtest.cc

## Cutee unit testing directives
# your test files here
test_files = t.protein.h

noinst_PROGRAMS = cutee runtest

cutee_SOURCES = cutee.cc cutee.h

runtest_LDADD = $(top_builddir)/src/evolver/libevolver.a $(top_builddir)/src/folder/libfolder.a

autocutee.mk: cutee Makefile.am $(test_files)
	$(CUTEE) -k -o autocutee.mk $(test_files)

include autocutee.mk