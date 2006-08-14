#! /usr/local/bin/python

import sys, os, math, string
sys.path = [os.path.expanduser("~/research/lib")]+sys.path
import folder
import translate

folder.init(5)
for line in file(os.path.expanduser("~/research/trsim/analysis/trs599ca6-genes.txt"), 'r').readlines()[2:]:
	flds = line.strip().split()
	gene = flds[2]
	prot = translate.Translate(gene)
	(sid, dg) = folder.foldProtein(prot)
	print sid, dg, prot
	

