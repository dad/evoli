## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/folder/ -I$(top_srcdir)/src/evolver/

sequence_generator_SOURCES = sequence-generator.cc
sequence_generator_LDADD = $(top_builddir)/src/folder/libfolder.a

structure_printer_SOURCES = structure-printer.cc
structure_printer_LDADD = $(top_builddir)/src//folder/libfolder.a

misfold_SOURCES = misfold.cc
misfold_LDADD = $(top_builddir)/src/evolver/libevolver.a $(top_builddir)/src/folder/libfolder.a

bin_PROGRAMS = sequence-generator structure-printer misfold
