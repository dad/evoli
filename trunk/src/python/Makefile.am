## Makefile.am

all: $(top_builddir)/src/folder/libfolder.a
	python setup.py build $(top_builddir) $(top_srcdir)

install: $(top_builddir)/src/folder/libfolder.a
	python setup.py install $(top_builddir) $(top_srcdir)

clean:
	rm -rf build

maintainer-clean: clean
	rm Makefile Makefile.in
