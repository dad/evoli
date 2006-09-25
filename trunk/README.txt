E.voli 0.1
==========

To build this project from svn:
   ./bootstrap
   ./configure
   make

The build structure has been tested with automake v1.9, and may
not work with earlier versions of automake.

For a quick test of some components:
	utils/sequence-generator 25 -5 2 1 599 > tmp-str.txt
	utils/structure-printer tmp-str.txt
	rm tmp-str.txt

For unit tests:
	test/runtest

To use Python bindings:
	cd src/python
	make
	cd ../..
	python test/folder-test.py
