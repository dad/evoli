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
sys.path = [os.path.expanduser('~/research/lib/'), os.path.expanduser('~/research/lib/src')] + sys.path
import folder, decoyfolder, misfold, mistrans
import translate

# The 20 canonical amino acids
aas = 'ACDEFGHIKLMNPQRSTVWY'
# Initialize the protein folder -- very important!
# side_length refers to the protein.  E.g., in the 5x5
# model, side_length=5.

if True:
	print "****\nTesting folder module"
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
	sid = random.randint(0,1080)
	dg = -2
	gene = folder.getSequenceForStructure(sid, dg)
	(new_sid, new_dg) = folder.fold(translate.Translate(gene))
	assert(sid == new_sid and new_dg <= dg)

if False:
	print "\n****\nTesting misfold module with compact lattice folder"
	side_length = 5
	struct_id = 599
	target_fraction_accurate = 0.85
	ca_cost = 5
	max_free_energy = -5

	prot_length = side_length*side_length
	misfold.init(folder, prot_length, struct_id, max_free_energy, ca_cost, target_fraction_accurate, 111)
	err_rate = misfold.getErrorRate()
	gene = "ATTATTGTCTCGAAGGGTGCTATCTCCGCCGTCAGTTCCTTCGCAAAGTACATCTTCTTGCTTCTAACTAAAGAC"
	(facc, frob, ftrunc, ffold) = misfold.calcOutcomes(gene);
	print "%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % ("calcOutcomes     ", facc, frob, ftrunc, ffold)
	(facc, frob, ftrunc, ffold) = misfold.countOutcomes(gene, 10000);
	print "%s\t%d\t%d\t%d\t%d" % ("countOutcomes      ", facc, frob, ftrunc, ffold)
	(facc, frob, ftrunc, ffold) = misfold.countOutcomeFractions(gene, 10000);
	print "%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % ("countOutcomeFractions", facc, frob, ftrunc, ffold)

if False:
	print "\n****\nTesting decoyfolder module"
	sid = 1
	prot_length = 300
	log_nconf = 10*math.log(10)
	map_file = os.path.abspath("test/data/williams_contact_maps/maps.txt")
	map_dir = os.path.abspath("test/data/williams_contact_maps/")+"/"

	decoyfolder.init(prot_length, log_nconf, map_file, map_dir)

	gene = decoyfolder.getSequenceForStructure(prot_length, sid, max_free_energy)
	print gene
	prot = translate.Translate(gene)
	(sid, dg) = decoyfolder.fold(prot)
	print sid, dg
	print "sid g"
	for i in range(10):
		# Create a random polypeptide
		prot = ''.join([random.choice(aas) for xi in range(prot_length)])
		# Fold it and retrieve its lowest-free-energy conformation, sid, and its
		# free energy of folding, dg.
		(sid, dg) = decoyfolder.fold(prot)
		# Print them out
		print "%d\t%1.3f" % (i, dg)

if False:
	print "\n****\nTesting misfold module with decoy contact folder"
	target_fraction_accurate = 0.85
	sid = 1
	ca_cost = 5
	max_free_energy = -5
	prot_length = 300

	gene = decoyfolder.getSequenceForStructure(prot_length, sid, max_free_energy)
	misfold.init(decoyfolder, prot_length, sid, max_free_energy, ca_cost, target_fraction_accurate, 111)
	#err_rate = misfold.getErrorRate()
	(facc, frob, ftrunc, ffold) = misfold.calcOutcomes(gene);
	print "%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % ("calcOutcomes     ", facc, frob, ftrunc, ffold)
	(facc, frob, ftrunc, ffold) = misfold.countOutcomes(gene, 1000);
	print "%s\t%d\t%d\t%d\t%d" % ("countOutcomes      ", facc, frob, ftrunc, ffold)
	(facc, frob, ftrunc, ffold) = misfold.countOutcomeFractions(gene, 1000);
	print "%s\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % ("countOutcomeFractions", facc, frob, ftrunc, ffold)

if False:
	print "\n****\nTesting mistrans module with decoy contact folder"
	target_fraction_accurate = 0.85
	sid = 1
	ca_cost = 5
	max_free_energy = -5
	prot_length = 300

	gene = decoyfolder.getSequenceForStructure(prot_length, sid, max_free_energy)
	misfold.init(decoyfolder, prot_length, sid, max_free_energy, ca_cost, target_fraction_accurate, 111)
	n = 100
	only_mistrans = True
	prots = mistrans.getTranslatedProteins(gene, n, only_mistrans)
	assert( len(prots) == n )
	for p in prots:
		print p

if True:
	print "\n****\nTesting mistrans module with compact lattice folder"
	side_length = 5
	struct_id = 599
	target_fraction_accurate = 0.85
	ca_cost = 5
	max_free_energy = -5

	prot_length = side_length*side_length
	gene = folder.getSequenceForStructure(struct_id, max_free_energy)
	mistrans.init(folder, prot_length, struct_id, max_free_energy, ca_cost, target_fraction_accurate, 111)
	n = 100
	only_mistrans = True
	prots = mistrans.getTranslatedProteins(gene, n, only_mistrans)
	assert( len(prots) == n )
	for p in prots:
		print p

