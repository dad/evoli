import random, math, os, sys
# Import the modules
sys.path = sys.path + [os.path.expanduser('~/research/lib/src'), os.path.expanduser('~/research/lib')]
import folder, decoyfolder, misfold, stats
import translate

# The 20 canonical amino acids
aas = 'ACDEFGHIKLMNPQRSTVWY'
# Initialize the protein folder -- very important!
# side_length refers to the protein.  E.g., in the 5x5
# model, side_length=5.

def main():
	prot = sys.argv[1]
	folder.init(int(math.sqrt(len(prot))))
	(sid, dg) = folder.fold(prot)
	print sid, dg


main()
