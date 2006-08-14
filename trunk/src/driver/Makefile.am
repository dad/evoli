## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/folder/ -I$(top_srcdir)/src/evolver

driver_SOURCES = translation-experiment.cc tr-driver.cc
driver_LDADD = $(top_builddir)/src/evolver/libevolver.a $(top_builddir)/src/folder/libfolder.a
bin_PROGRAMS = driver
