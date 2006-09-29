## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/tools -I$(top_srcdir)/src/folder/ -I$(top_srcdir)/src/evolver

driver_SOURCES = translation-experiment.cc tr-driver.cc
driver_LDADD = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/tools/libtools.a

cutoff_driver_SOURCES = translation-experiment.cc tr-cutoff-driver.cc
cutoff_driver_LDADD = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/tools/libtools.a

bin_PROGRAMS = driver cutoff-driver
