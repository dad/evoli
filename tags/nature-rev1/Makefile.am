## Makefile.am

BASE_SUBDIRS = src test

# build apidoc if requested
if BUILD_APIDOC
  SUBDIRS = $(BASE_SUBDIRS) doc
else
  SUBDIRS = $(BASE_SUBDIRS)	
endif

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure configure.scan autoscan.log \
	config/missing config/depcomp config/install-sh
