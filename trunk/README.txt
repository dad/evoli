E.voli 0.1
==========

To build this project:
   ./bootstrap
   ./configure
   make

For a quick test of some components:
	utils/sequence-generator -5 2 1 599 > tmp-str.txt
	utils/structure-printer tmp-str.txt
	rm tmp-str.txt

For unit tests:
	test/runtest

