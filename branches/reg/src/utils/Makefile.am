## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/tools/ \
		-I$(top_srcdir)/src/gene/ \
		-I$(top_srcdir)/src/folder/ \
		-I$(top_srcdir)/src/evolver/

libraries = $(top_builddir)/src/evolver/libevolver.a \
		$(top_builddir)/src/folder/libfolder.a \
		$(top_builddir)/src/gene/libgene.a \
		$(top_builddir)/src/tools/libtools.a

gb_analyzer_SOURCES = gb-analyzer.cc
gb_analyzer_LDADD   = $(libraries)

sequence_generator_SOURCES = sequence-generator.cc
sequence_generator_LDADD = $(libraries)

decoy_sequence_generator_SOURCES = decoy-sequence-generator.cc
decoy_sequence_generator_LDADD = $(libraries)

decoy_sequence_analyzer_SOURCES = decoy-sequence-analyzer.cc
decoy_sequence_analyzer_LDADD = $(libraries)

neutral_evolve_SOURCES = neutral-evolve.cc
neutral_evolve_LDADD = $(libraries)

structure_printer_SOURCES = structure-printer.cc
structure_printer_LDADD = $(libraries)

misfold_SOURCES = misfold.cc
misfold_LDADD = $(libraries)

get_weights_SOURCES = get-weights.cc
get_weights_LDADD = $(libraries)

evolved_dg_dist_SOURCES = evolved-dg-dist.cc
evolved_dg_dist_LDADD = $(libraries)

snp_mistrans_stability_SOURCES = snp-mistrans-stability.cc
snp_mistrans_stability_LDADD = $(libraries)

mistime_SOURCES = mistime.cc
mistime_LDADD = $(libraries)

bin_PROGRAMS = sequence-generator decoy-sequence-generator decoy-sequence-analyzer \
	structure-printer misfold get-weights gb-analyzer evolved-dg-dist snp-mistrans-stability \
	neutral-evolve mistime

CLEANFILES = pdbcontacts.pyc pdb.pyc

