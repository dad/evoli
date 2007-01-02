## Makefile.am

include $(top_srcdir)/config/Make-rules
INCLUDES = -I$(top_srcdir)/src/tools -I$(top_srcdir)/src/folder

lib_LIBRARIES = libevolver.a

##libevolver_a_LIBADD = $(top_builddir)/src/folder/libfolder.a
libevolver_a_SOURCES = fitness-evaluator.cc population.cc mutator.cc

test_SOURCES = test.cc

test_LDADD = $(top_builddir)/src/evolver/libevolver.a \
                $(top_builddir)/src/folder/libfolder.a \
                $(top_builddir)/src/tools/libtools.a

bin_PROGRAMS = test