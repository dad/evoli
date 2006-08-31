#! /usr/local/bin/python

import random, math, os, sys
# Import the modules
import folder, decoyfolder
#sys.path = [os.path.expanduser('~/research/lib')] + sys.path
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

	for i in range(100000):
		# Create a random polypeptide
		prot = ''.join([random.choice(aas) for xi in range(prot_length)])
		# Fold it and retrieve its lowest-free-energy conformation, sid, and its
		# free energy of folding, dg.
		(sid, dg) = folder.fold(prot)
		# Print them out
		if dg < 0:
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
	log_nconf = 160*math.log(10)
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
		
if True:
	prot_length = 300
	log_nconf = 160*math.log(10)
	map_file = os.path.abspath("test/data/rand_contact_maps/maps.txt")
	map_dir = os.path.abspath("test/data/rand_contact_maps/")+"/"

	decoyfolder.init(prot_length, log_nconf, map_file, map_dir)
	print "sid dg"
	p = "PRPEEEKKKREREEKRRKEDKLERIRDLPRKILKMIVEPKRRKKGETEDDDEKESKRREEMEKFKREFFTICIKLLECEEEMARRREKRREEEDIDSLRELMKDCRRFIDDPRRVEQQSQRLDFRSRRKLEDEKDDEDKRKPDFLFEFEMCEEDMRRRPLDRVKDICRVCCEMDEEEEIREEEEFFRPEEEDMKLKSFRESFKDVRRCILRKFEKSRREKSAEFLRHEIPMFSSEDEEDRKKKDRRRQRPMMRHFMKRIKEKEEERKKREFKEQEEPKPKSFKWKTEEEMEELGEQEKRV"
	(sid, dg) = decoyfolder.fold(p)
	print sid, dg
