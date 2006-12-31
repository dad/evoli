## Makefile.am

SUBDIRS = src test

# build apidoc if requested
if BUILD_APIDOC
  SUBDIRS += doc
endif

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure configure.scan autoscan.log \
	config/missing config/depcomp config/install-sh
