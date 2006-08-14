from distutils.core import setup, Extension
import os

folder = Extension('folder',
				   sources = ['foldermodule.cc'],
				   include_dirs = ['../src/'],
				   libraries = ['folder'],
				   library_dirs = ['../src/'])

setup (name = 'Folder',
	          version = '0.1',
	          description = 'Module exposing lattice protein folding to Python',
	          ext_modules = [folder])

