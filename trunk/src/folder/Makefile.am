## Makefile.am

include $(top_srcdir)/config/Make-rules

INCLUDES = -I$(top_srcdir)/src/tools/ \
		-I$(top_srcdir)/src/gene/

lib_LIBRARIES = libfolder.a

libfolder_a_SOURCES = compact-lattice-folder.cc \
		decoy-contact-folder.cc \
		protein-contact-energies.cc
