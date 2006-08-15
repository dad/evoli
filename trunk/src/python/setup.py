from distutils.core import setup, Extension
import os, sys

srcdir = sys.argv.pop()
print "top_srcdir =", srcdir
builddir = sys.argv.pop()
print "top_builddir =", builddir

folder = Extension('folder',
				   sources = ['foldermodule.cc'],
				   include_dirs = [srcdir + '/src/folder', srcdir + '/src/evolver'],
				   libraries = ['folder', 'evolver'],
				   library_dirs = [srcdir + '/src/folder', srcdir + '/src/evolver'])

setup (name = 'Folder',
	          version = '0.1',
	          description = 'Module exposing lattice protein folding to Python',
	          ext_modules = [folder])

