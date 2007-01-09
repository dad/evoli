## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/tools/ \
		-I$(top_srcdir)/src/gene/ \
		-I$(top_srcdir)/src/folder/ \
		-I$(top_srcdir)/src/evolver/

tr_driver_SOURCES = translation-experiment.cc tr-driver.cc
tr_driver_LDADD = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/gene/libgene.a \
		$(top_builddir)/src/tools/libtools.a

bin_PROGRAMS = tr-driver
