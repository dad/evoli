#! /usr/local/bin/python

import random
# Import the module
import folder

# The 20 canonical amino acids
aas = 'ACDEFGHIKLMNPQRSTVWY'
# Initialize the protein folder -- very important!
# side_length refers to the protein.  E.g., in the 5x5
# model, side_length=5.
side_length = 5
folder.init(side_length, "hello")

for i in range(100):
	# Create a random polypeptide
	prot = ''.join([random.choice(aas) for i in range(side_length*side_length)])
	# Fold it and retrieve its lowest-free-energy conformation, sid, and its
	# free energy of folding, dg.
	(sid, dg) = folder.foldProtein(prot)
	# Print them out
	print "%d\t%1.3f\t%s" % (sid, dg, prot)
