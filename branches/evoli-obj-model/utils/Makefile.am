## Makefile.am

include $(top_srcdir)/config/Make-rules

sequence_generator_SOURCES = sequence-generator.cc
sequence_generator_LDADD = $(top_builddir)/src/libfolder.a

structure_printer_SOURCES = structure-printer.cc
structure_printer_LDADD = $(top_builddir)/src/libfolder.a

misfold_SOURCES = misfold.cc
misfold_LDADD = $(top_builddir)/src/libevolver.a

bin_PROGRAMS = sequence-generator structure-printer misfold
