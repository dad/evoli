#! /usr/local/bin/python

# This file is part of the evoli project.
# Copyright (C) 2006 Allan Drummond <dadrummond@gmail.com>
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


import random, math, os, sys
# Import the modules
import folder, decoyfolder, misfold
#sys.path = [os.path.expanduser('~/research/lib/lib')] + sys.path
#import translate

# The 20 canonical amino acids
aas = 'ACDEFGHIKLMNPQRSTVWY'
# Initialize the protein folder -- very important!
# side_length refers to the protein.  E.g., in the 5x5
# model, side_length=5.

if False:
	side_length = 5
	prot_length = side_length*side_length
	folder.init(side_length)

	for i in range(10000):
		# Create a random polypeptide
		prot = ''.join([random.choice(aas) for xi in range(prot_length)])
		# Fold it and retrieve its lowest-free-energy conformation, sid, and its
		# free energy of folding, dg.
		(sid, dg) = folder.fold(prot)
		# Print them out
		if dg < -3:
			print "%d\t%d\t%1.3f\t%s" % (i, sid, dg, prot)
		'''
		print "sid g"
		for i in range(1):
			# Create a random polypeptide
			prot = ''.join([random.choice(aas) for i in range(prot_length)])
			# Fold it and retrieve its lowest-free-energy conformation, sid, and its
			# free energy of folding, dg.
			(sid, dg) = folder.fold(prot)
			for j in range(50):
				G = folder.getEnergy(prot, j)	
				# Print them out
				print "%d\t%1.3f\t%s" % (sid, dg, prot)	
				print "%d\t%1.3f" % (j, G)
		'''		

if False:
	prot_length = 300
	log_nconf = 10*math.log(10)
	map_file = os.path.abspath("test/data/rand_contact_maps/maps.txt")
	map_dir = os.path.abspath("test/data/rand_contact_maps/")+"/"

	decoyfolder.init(prot_length, log_nconf, map_file, map_dir)
	print "sid g"
	for i in range(1):
		# Create a random polypeptide
		prot = ''.join([random.choice(aas) for i in range(prot_length)])
		# Fold it and retrieve its lowest-free-energy conformation, sid, and its
		# free energy of folding, dg.
		(sid, dg) = decoyfolder.fold(prot)
		for j in range(50):
			G = decoyfolder.getEnergy(prot, j)	
			# Print them out
			print "%d\t%1.3f" % (j, G)
		
if False:
	prot_length = 300
	log_nconf = 10*math.log(10)
	map_file = os.path.abspath("test/data/rand_contact_maps/maps.txt")
	map_dir = os.path.abspath("test/data/rand_contact_maps/")+"/"

	decoyfolder.init(prot_length, log_nconf, map_file, map_dir)
	print "sid dg"
	p = "PRPEEEKKKREREEKRRKEDKLERIRDLPRKILKMIVEPKRRKKGETEDDDEKESKRREEMEKFKREFFTICIKLLECEEEMARRREKRREEEDIDSLRELMKDCRRFIDDPRRVEQQSQRLDFRSRRKLEDEKDDEDKRKPDFLFEFEMCEEDMRRRPLDRVKDICRVCCEMDEEEEIREEEEFFRPEEEDMKLKSFRESFKDVRRCILRKFEKSRREKSAEFLRHEIPMFSSEDEEDRKKKDRRRQRPMMRHFMKRIKEKEEERKKREFKEQEEPKPKSFKWKTEEEMEELGEQEKRV"
	(sid, dg) = decoyfolder.fold(p)
	print sid, dg

if True:
	side_length = 5
	struct_id = 599
	target_fraction_accurate = 0.1
	ca_cost = 4
	max_free_energy = -5
	
	misfold.init(side_length, struct_id, max_free_energy, ca_cost, target_fraction_accurate)
	#print misfold.getErrorRate()
	gene = "CCCCTGTACCGTACGACGAAATCTAACACTGGATCATGGCCTTCTGATTGGAAACCCCTACCTTATGAGTCAAAG"
	#(facc, frob, ftrunc, ffold) = misfold.calcOutcomes(gene);
	#print "%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % ("mis", facc, frob, ftrunc, ffold)
	if True:
		f = file(os.path.expanduser('~/research/trsim/data/trs599ca4nmut-genes.txt'),'r')
		for line in f.readlines()[2:]:
			flds = line.strip().split('\t')
			gene = flds[-1]
			# print gene
			#(facc, frob, ftrunc, ffold) = misfold.calcOutcomes(gene);
			#print "%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % (flds[0], facc, frob, ftrunc, ffold)
			(facc, frob, ftrunc, ffold) = misfold.countOutcomes(gene, 1000);
			print "%s\t%d\t%d\t%d\t%d" % (flds[0], facc, frob, ftrunc, ffold)
			#print line,
		f.close()
	
