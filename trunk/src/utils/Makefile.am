## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/tools/ -I$(top_srcdir)/src/folder/ -I$(top_srcdir)/src/evolver/

libraries = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/tools/libtools.a

sequence_generator_SOURCES = sequence-generator.cc
sequence_generator_LDADD = $(libraries)

decoy_sequence_generator_SOURCES = decoy-sequence-generator.cc
decoy_sequence_generator_LDADD = $(libraries)

decoy_sequence_analyzer_SOURCES = decoy-sequence-analyzer.cc
decoy_sequence_analyzer_LDADD = $(libraries)

structure_printer_SOURCES = structure-printer.cc
structure_printer_LDADD = $(libraries)

misfold_SOURCES = misfold.cc
misfold_LDADD = $(libraries)

bin_PROGRAMS = sequence-generator decoy-sequence-generator decoy-sequence-analyzer structure-printer misfold

CLEANFILES = pdbcontacts.pyc pdb.pyc

