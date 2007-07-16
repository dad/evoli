## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/tools/

lib_LIBRARIES = libgene.a

libgene_a_SOURCES = genetic-code.cc translator.cc protein.cc codon.cc mutator.cc upstream.cc expressible-gene.cc


mutator_test_SOURCES = mutator-test.cc

mutator_test_LDADD = $(top_builddir)/src/gene/libgene.a \
                $(top_builddir)/src/tools/libtools.a

bin_PROGRAMS = mutator-test
