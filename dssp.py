#!/home/helene/anaconda3/bin/python3.6

from optparse import OptionParser
import datetime as dt
import re
from Bio.SeqUtils import seq1
from Bio.PDB import *
import math as m

def makeHeader():
	"""
	Make and return the header of the dssp output, using data from the pdb file
	"""
	with open(opt.input,'r') as f:
		pdb = f.read()

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
	i = 0
	for chain in structure.get_chains():
		start = 9999
		start = next(res.id[1] for res in chain.get_residues() if (res.id[1] < start))
		#end = next(res.id[1] for res in chain.get_residues() if (res.id[0] == 'W'))
		for res in chain.get_residues():
			if (res.get_id()[0] == 'W'):
				break
			end = res.get_id()[1]

		for res in range(start,end):
			i += 1
			line = {}
			line['index'] = i
			line['res_nb'] = res
			line['aa'] = seq1(chain[res].get_resname())
			line['chain'] = chain.id
			line['x-ca'] = chain[res]['CA'].get_coord()[0]
			line['y-ca'] = chain[res]['CA'].get_coord()[1]
			line['z-ca'] = chain[res]['CA'].get_coord()[2]
			n = chain[res]['N'].get_vector() 
			ca = chain[res]['CA'].get_vector() 
			c = chain[res]['C'].get_vector()
			# PHI calculation
			try:
				cp = chain[res-1]['C'].get_vector() 
				line['phi'] = (calc_dihedral(cp, n, ca, c)*180)/m.pi # degree = (radian*180)/pi
				
			except:
				line['phi'] = 360.0
			# PSI calculation
			try:
				nn = chain[res+1]['N'].get_vector()
				line['psi'] = (calc_dihedral(n, ca, c, nn)*180)/m.pi
			except:
				line['psi'] = 360.0

			# ALPHA angle calculation
			try:
				cap = chain[res-1]['CA'].get_vector()
				can = chain[res+1]['CA'].get_vector()
				cann = chain[res+2]['CA'].get_vector()
				line['alpha'] = (calc_dihedral(cap,ca,can,cann)*180)/m.pi;
			except:
				line['alpha'] = 360.0
			if (line['alpha'] < 0):
				chirality = '-';
			else:
				chirality = '+';

			# KAPPA angle
			try:
				capp = chain[res-2]['CA'].get_vector()
				cann = chain[res+2]['CA'].get_vector()
				line['kappa'] = 0.0
			except:
				line['kappa'] = 360.0
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
	print(makeHeader())
	dssp = DSSPlines()

#	print("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN")
	print("  #  RESIDUE AA KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA")
	
	for i in range(0,len(dssp)):
		l = dssp[i]
		print("{:>5d}{:>5d}{:>2s}{:>2s}{:>7.1f}{:>6.1f}{:>6.1f}{:>6.1f}{:>7.1f}{:>7.1f}{:>7.1f}"\
			.format(l['index'],l['res_nb'],l['chain'],l['aa'],l['kappa'],l['alpha'],l['phi'],l['psi'],l['x-ca'],l['y-ca'],l['z-ca']))

	

