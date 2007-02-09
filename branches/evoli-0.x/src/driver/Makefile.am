## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/tools -I$(top_srcdir)/src/folder/ -I$(top_srcdir)/src/evolver

tr_driver_SOURCES = translation-experiment.cc tr-driver.cc
tr_driver_LDADD = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/tools/libtools.a

neut_driver_SOURCES = neut-driver.cc
neut_driver_LDADD = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/tools/libtools.a

aaseq_driver_SOURCES = aa-seq-driver.cc
aaseq_driver_LDADD = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/tools/libtools.a


struct_driver_SOURCES = struct-driver.cc
struct_driver_LDADD = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/tools/libtools.a


bin_PROGRAMS = tr-driver neut-driver aaseq-driver struct-driver
