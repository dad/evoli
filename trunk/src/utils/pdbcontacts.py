#! /usr/local/bin/python
"""Script for calculating PDB contacts.

Written by D. Allan Drummond, 2006.

"""

import sys, string, os
import pdb

ARG_PDB_FILE = 'pdb'
ARG_OUTPUT_FILE = 'o'
ARG_HELP = 'help'
ARG_CHAINS = 'chains'

def parse_arguments(args):
	# Turn linear arguments into a dictionary of (option, [values,...]) pairs
	arg_dict = {}
	key = None
	for arg in args[1:]:
		if arg[0] == '-':
			key = arg[1:]
			arg_dict[key] = None
		else:
			if arg_dict.has_key(key):
				if arg_dict[key]:
					if type(arg_dict[key]) is list:
						arg_dict[key] = arg_dict[key]+[arg]
					else:
						arg_dict[key] = [arg_dict[key],arg]
				else:
					arg_dict[key] = arg
			else:
				arg_dict[key] = arg
	return arg_dict

def print_usage(args):
	print 'Usage: python', args[0].split(os.path.sep)[-1], " [options]"
	print 'Options:\n',\
		"\t-%s <PDB file>\n" % ARG_PDB_FILE, \
		"\t[-%s <PDB chains list>]\n" % ARG_CHAINS,\
		"\t[-%s <contacts output file>]" % ARG_OUTPUT_FILE

def confirm_arguments(arg_dict):
	# Are arguments okay?
	res = True
	arg_keys = arg_dict.keys()

	try:
		if len(arg_keys) == 0:
			res = False
			return
			
		if not ARG_PDB_FILE in arg_keys:
			print "  You must provide a PDB file (-%s <file>)" % ARG_PDB_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_PDB_FILE]):
			print "  Can't find PDB file %s" % arg_dict[ARG_PDB_FILE]
			res = False
		
	except:
		res = False
	return res

def getPDBContacts(residues, contact_distance, chain_ids, nonbonded_only=False):
	# Get set of residues within contact_distance angstroms in a PDB
	contacts = []
	# Loop through all residue pairs
	for i in range(0,len(residues)-1):
		for j in range(i+1, len(residues)):
			if nonbonded_only and (j-i)<=1:
				continue
			resi = residues[i]
			resj = residues[j]
			# If both residues are present (gaps == None)
			if resi and resj: # and (resi.chain in chain_ids) and (resj.chain in chain_ids):
				# If residues are within contact_distance angstroms of each other, add as a contact
				contact = resi.isContactCBeta(resj, contact_distance)
				if contact:
					contacts.append((i, j, resi, resj))
	return contacts

def writeContactFile(contacts, outfile):
	for k in range(len(contacts)):
		(i, j, resi, resj) = contacts[k]
		#outfile.write("%d\t%s\t%d\t%s\t%f\n" % (i, pdb.three_to_one_map[resi.residue], j, pdb.three_to_one_map[resj.residue], resi.getDistanceCBeta(resj)))
		outfile.write("%d\t%s\t%d\t%s\n" % (i, pdb.three_to_one_map[resi.residue], j, pdb.three_to_one_map[resj.residue]))		

def main(args):
	arg_dict = parse_arguments(args)
	if not confirm_arguments(arg_dict):
		if args[0].split(os.path.sep)[-1] == "pdbcontacts.py":
			print_usage(args)
		return

	# Flags and values
	
	# Inputs:
	#	The PDB file name.
	pdb_file = arg_dict[ARG_PDB_FILE]

	# The PDB chains
	# Many PDB files include multiple chains.  The chain_identifier list includes those
	# chains which correspond to the protein whose contacts are being evaluated.
	# Most often, chain 'A' (in the case of multiple chains) or chain ' ' (only one chain)
	# will be the appropriate choice.
	if arg_dict.has_key(ARG_CHAINS):
		chains = arg_dict[ARG_CHAINS]
		if type(chains) is list:
			chain_identifiers = chains + [' ']
		else:
			chain_identifiers = [chains, ' ']
	else:
		chain_identifiers = ['A',' ']

	# Generate the contacts
	# Read in the PDB file to create a list of residues.
	raw_residues = pdb.File().read(file(pdb_file, 'r'))
	residues = [res for res in raw_residues if res.chain in chain_identifiers]
	#sequence = pdb.sequence(residues, chain_ids)
	
	# 	The contact file name for output.
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		contact_file = file(arg_dict[ARG_OUTPUT_FILE], 'w')
	else:
		contact_file = sys.stdout
		

	# With an aligned set of residues and parents, we can now compute the SCHEMA contacts.
	# Note that for more than two parents, some of these contacts may only be broken by 
	# specific chimera patterns.
	contact_distance = 6.0  # Residues closer than this distance, in angstroms, are in contact.
	pdb_contacts = getPDBContacts(residues, contact_distance, chain_identifiers, nonbonded_only=True)
	writeContactFile(pdb_contacts, contact_file)
	if False: # DAD: debugging
		print residues[294].getDistanceCBeta(residues[296])
		print len(pdb_contacts)
	
	if not contact_file == sys.stdout:
		contact_file.close()

		
def main_wrapper():
	main(sys.argv)

main_wrapper()
