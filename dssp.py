#!/usr/bin/python3.6

from optparse import OptionParser
import datetime as dt
import re
from Bio.SeqUtils import seq1
from Bio.PDB import *
import math as m

def readPDB():
	with open(opt.input,'r') as f:
		pdb = f.read()
	return(pdb)

def makeHeader():
	"""
	Make and return the header of the dssp output, using data from the pdb file
	"""
	header = "==== Secondary Structure Assignment using DSSP method ====\nDATE\t\t{}\n".format(dt.date.today())
	header += "REFERENCE\tW. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n"
	# regex : 
	header_regex = re.compile("HEADER[ |2-9]{1,}([^\n]*)")
	compnd_regex = re.compile("COMPND[ |2-9]{1,}([^\n]*)")
	source_regex = re.compile("SOURCE[ |2-9]{1,}([^\n]*)")
	author_regex = re.compile("AUTHOR {1,}([^\n]*)")
	header += "HEADER\t\t{}\n".format(' '.join(map(str.rstrip,header_regex.findall(pdb))))
	header += "COMPND\t\t{}\n".format(' '.join(map(str.rstrip,compnd_regex.findall(pdb)))) # all COMPND
	header += "SOURCE\t\t{}\n".format(' '.join(map(str.rstrip,source_regex.findall(pdb)))) # all SOURCE
	header += "AUTHOR\t\t{}\n".format(author_regex.findall(pdb)[0].rstrip())
	return(header)

def DSSPlines():
	p = PDBParser()
	structure = p.get_structure('A',opt.input)
	dssp = []
	for model in structure:
		for chain in model:
			for res in range(1,len(chain)+1):
				line = {}
				line['chain'] = chain.id
				line['index'] = res
				line['aa'] = seq1(chain[res].get_resname())
				line['x-ca'] = round(chain[res]['CA'].get_coord()[0],1)
				line['y-ca'] = round(chain[res]['CA'].get_coord()[1],1)
				line['z-ca'] = round(chain[res]['CA'].get_coord()[2],1)
				n = chain[res]['N'].get_vector() 
				ca = chain[res]['CA'].get_vector() 
				c = chain[res]['C'].get_vector()
				#PHI calculation
				try:
					cp = chain[res-1]['C'].get_vector() 
					line['phi'] = round((calc_dihedral(cp, n, ca, c)*180)/m.pi,1) # degree = (radian*180)/pi
					
				except:
					line['phi'] = 360.0
				#PSI calculation
				try:
					nn = chain[res+1]['N'].get_vector()
					line['psi'] = round((calc_dihedral(n, ca, c, nn)*180)/m.pi,1)
				except:
					line['psi'] = 360.0
				dssp.append(line)
	return(dssp)

if __name__ == "__main__":
	parser = OptionParser(usage="usage: %prog [options]",version="%prog 1.0")
	parser.add_option("-i", "--input",dest='input',
	                  help="The file name of a PDB formatted file containing the protein structure data.", metavar="FILE")
	parser.add_option("-o", "--output",dest='output',
	                  help="The  file  name  of  a  DSSP  file to create.", metavar="FILE")
	(opt, args) = parser.parse_args()
	if not opt.input:
		parser.error('Input pdb file not given.')
	pdb = readPDB()
	makeHeader()
	dssp = DSSPlines()
	print(dssp[0])
	print(dssp[10])
	print(dssp[-1])

	  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN

