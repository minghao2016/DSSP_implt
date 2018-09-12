#!/home/helene/anaconda3/bin/python3.6

from optparse import OptionParser
from subprocess import call
from subprocess import Popen, PIPE
import datetime as dt
import re
from Bio.SeqUtils import seq1
from Bio.PDB import *
import numpy as np

q1 = 0.41
q2 = 0.20
f = 332

def calc_energy(i):
	try:
		n_i5 = chain[res+i]['N']
		h_i5 = chain[res+i]['H']
		o = chain[res]['O']
		c = chain[res]['C']
		rON = n_i5 - o
		rCH = h_i5 - c
		rOH = h_i5 - o
		rCN = n_i5 - c
		E = q1*q2*(1/rON + 1/rCH - 1/rOH - 1/rCN)*f
		line['Ei'+str(i)] = E
		if (E < -0.5):
			line['i'+str(i)] = 'Y'
		else:
			line['i'+str(i)] = '.'
	except:
		line['Ei'+str(i)] = 0
		line['i'+str(i)] = '.'

def lineHeader(dic):
	l = ""
	for mol_id,mol_items in dic.items():
		l += "MOL_ID: " + mol_id + "; "
		for key, item in mol_items.items():
			if (item != ""):
				l += key.upper() + ": " + item.upper() + "; "
	return(l)

def makeHeader():
	"""
	Make and return the header of the dssp output, using data from the pdb file
	"""
	header = "==== Secondary Structure Assignment using DSSP method ====\nDATE\t\t{}\n".format(dt.date.today())
	header += "REFERENCE\tW. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n"
	header += "HEADER\t\t{}{:>28}\n".format(structure.header["head"].upper(),structure.header["deposition_date"])
	header += "COMPND\t\t{}\n".format(lineHeader(structure.header["compound"])) # all COMPND
	header += "SOURCE\t\t{}\n".format(lineHeader(structure.header["source"])) # all SOURCE
	header += "AUTHOR\t\t{}".format(structure.header["author"].upper())
	return(header+"\n")

def displayResults():
	#	print("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN")
	header = makeHeader()
	#descp = "  #  RESIDUE AA    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
	descp = "  #  RESIDUE AA   E"
	if (opt.output):
		with open(opt.output,'w') as filout:
			filout.write(header+descp+'\n')
			for i in range(0,len(dssp)):
				l = dssp[i]
				filout.write("{:>5d}{:>5d}{:>2s}{:>2s}{:>9.3f}{:>6.1f}{:>6.1f}{:>6.1f}{:>6.1f}{:>7.1f}{:>7.1f}{:>7.1f}\n"\
					.format(l['index'],l['res_nb'],l['chain'],l['aa'],l['tco'],l['kappa'],l['alpha'],l['phi'],l['psi'],l['x-ca'],l['y-ca'],l['z-ca']))
	else:
		print(header,descp,sep="")
		for i in range(0,len(dssp)):
			l = dssp[i]
			print("{:>5d}{:>5d}{:>2s}{:>2s}{:>8.2f}{:>3s}{:>8.2f}{:>3s}{:>8.2f}{:>3s}"\
				.format(l['index'],l['res_nb'],l['chain'],l['aa'],l['Ei5'],l['i5'],l['Ei4'],l['i4'],l['Ei3'],l['i3']))	
			#print("{:>5d}{:>5d}{:>2s}{:>2s}{:>9.3f}{:>6.1f}{:>6.1f}{:>6.1f}{:>6.1f}{:>7.1f}{:>7.1f}{:>7.1f}"\
			#	.format(l['index'],l['res_nb'],l['chain'],l['aa'],l['tco'],l['kappa'],l['alpha'],l['phi'],l['psi'],l['x-ca'],l['y-ca'],l['z-ca']))

if __name__ == "__main__":
	parser = OptionParser(usage="usage: %prog [options]",version="%prog 1.0")
	parser.add_option("-i", "--input",dest='input',
	                  help="The file name of a PDB formatted file containing the protein structure data.", metavar="FILE")
	parser.add_option("-o", "--output",dest='output',
	                  help="The  file  name  of  a  DSSP  file to create.", metavar="FILE")
	(opt, args) = parser.parse_args()
	if not opt.input:
		parser.error('Input pdb file not given.')

	call(["./bin/reduce -NOFLIP "+opt.input+" 1>"+opt.input+".H"+" 2>"+opt.input+".H.log"],shell=True)

	p = PDBParser()
	structure = p.get_structure(opt.input,opt.input+".H")
	dssp = []
	index,	nb_chains = 0, 0
	for chain in structure.get_chains():
		nb_chains += 1
		start = 9999
		start = next(res.id[1] for res in chain.get_residues() if (res.id[1] < start))
		#end = next(res.id[1] for res in chain.get_residues() if (res.id[0] == 'W'))
		for res in chain.get_residues():
			if (res.get_id()[0] == 'W'):
				break
			end = res.get_id()[1]

		for res in range(start,end):
			index += 1
			line = {}
			line['index'] = index
			line['res_nb'] = res
			line['aa'] = seq1(chain[res].get_resname())
			line['chain'] = chain.id
			line['x-ca'] = chain[res]['CA'].get_coord()[0]
			line['y-ca'] = chain[res]['CA'].get_coord()[1]
			line['z-ca'] = chain[res]['CA'].get_coord()[2]
			n = chain[res]['N'].get_vector() 
			ca = chain[res]['CA'].get_vector() 
			c = chain[res]['C'].get_vector()
			o = chain[res]['O'].get_vector()


			#HBONDS i+5 :
			calc_energy(5)
			calc_energy(4)
			calc_energy(3)


			# PHI calculation
			try:
				cp = chain[res-1]['C'].get_vector() 
				line['phi'] = (calc_dihedral(cp, n, ca, c)*180)/np.pi # degree = (radian*180)/pi
				
			except:
				line['phi'] = 360
			# PSI calculation
			try:
				nn = chain[res+1]['N'].get_vector()
				line['psi'] = (calc_dihedral(n, ca, c, nn)*180)/np.pi
			except:
				line['psi'] = 360

			# ALPHA angle calculation
			try:
				cap = chain[res-1]['CA'].get_vector()
				can = chain[res+1]['CA'].get_vector()
				cann = chain[res+2]['CA'].get_vector()
				line['alpha'] = (calc_dihedral(cap,ca,can,cann)*180)/np.pi
			except:
				line['alpha'] = 360
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
				line['kappa'] = 360
			#TCO
			try:
				cp = chain[res-1]['C'].get_vector()
				op = chain[res-1]['O'].get_vector()
				p1 = c - o
				p2 = cp - op
				x = np.dot(p1,p1) * np.dot(p2,p2)
				if(x > 0):
					line['tco'] = np.dot(p1,p2) / np.sqrt(x)
			except:
				line['tco'] = 0
			dssp.append(line)
	displayResults()