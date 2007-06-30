# This file is part of the evoli project.
# Copyright (C) 2006 Claus Wilke <cwilke@mail.utexas.edu>, Allan Drummond <dadrummond@gmail.com>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1

from distutils.core import setup, Extension
import os, sys

srcdir = sys.argv.pop()
print "top_srcdir =", srcdir
builddir = sys.argv.pop()
print "top_builddir =", builddir


folder = Extension('folder',
				   sources = ['foldermodule.cc'],
				   include_dirs = [srcdir + '/src/tools', \
							srcdir + '/src/gene', \
							srcdir+'/src/folder'],
				   libraries = ['folder', 'gene', 'tools'],
				   library_dirs = [builddir + '/src/tools', \
	 						builddir + '/src/gene', \
							builddir + '/src/folder'])

setup (name = 'CompactLatticeFolder',
	          version = '0.1',
	          description = 'Module exposing lattice protein folding to Python',
	          ext_modules = [folder])

decoyfolder = Extension('decoyfolder',
				   sources = ['decoyfoldermodule.cc'],
				   include_dirs = [srcdir + '/src/tools', \
							srcdir + '/src/gene', \
							srcdir+'/src/folder'],
				   libraries = ['folder', 'gene', 'tools'],
				   library_dirs = [builddir + '/src/tools', \
							builddir + '/src/gene', \
							builddir + '/src/folder'])

setup (name = 'DecoyContactFolder',
	          version = '0.1',
	          description = 'Module exposing decoy-contact protein folding to Python',
	          ext_modules = [decoyfolder])

codon = Extension('codon',
				   sources = ['codonmodule.cc'],
				   include_dirs = [srcdir + '/src/tools',\
							srcdir + '/src/gene'],
				   libraries = ['gene', 'tools' ],
				   library_dirs = [builddir + '/src/tools', \
							builddir + '/src/gene'])

setup (name = 'Codon',
	          version = '0.1',
	          description = 'Module exposing functionality related to calculating synonymous and nonsynonymous substitutions and sites to Python',
	          ext_modules = [codon])

folder = Extension('misfold',
				   sources = ['misfoldmodule.cc'],
				   include_dirs = [srcdir + '/src/tools', \
							srcdir + '/src/gene', \
							srcdir + '/src/evolver', \
							srcdir+'/src/folder'],
				   libraries = ['folder', 'gene', 'tools', 'evolver'],
				   library_dirs = [builddir + '/src/tools', \
	 						builddir + '/src/gene', \
	 						builddir + '/src/evolver', \
							builddir + '/src/folder'])

setup (name = 'Misfold',
	          version = '0.1',
	          description = 'Module exposing misfolding assays of lattice protein folding to Python',
	          ext_modules = [folder])
