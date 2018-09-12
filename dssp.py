#!/home/helene/anaconda3/bin/python3.6

import sys
import os
from Bio.SeqUtils import seq1
from Bio.PDB import *
import numpy as np
srcPath = 'src'
sys.path.append(os.path.abspath(srcPath))
from functions import *

################# FIXED VARIABLES #################

# For the calculation of the electrostatic
# interaction energy between two H-bonding groups
Q1, Q2 = 0.41, 0.20 # partial charges
F = 332 # factor

# Patterns
NONE = ' '
START = '>'
END = '<'
START_END = 'X'
###################################################

def residueDesc():
	"""
	Description of a residue
	index, residue number, aa name, chain number, coord (x,y,z)-Ca
	"""
	l['structure'] = NONE
	l['index'] = index
	l['res_nb'] = res
	l['aa'] = seq1(chain[res].get_resname())
	l['chain'] = chain.id
	l['x-ca'] = chain[res]['CA'].get_coord()[0]
	l['y-ca'] = chain[res]['CA'].get_coord()[1]
	l['z-ca'] = chain[res]['CA'].get_coord()[2]

def nTurnPattern():
	for n in range(3,6):
		if (testHbond(chain,res,n) == True):
			l[str(n)+'-turns'] = START
		else:
			l[str(n)+'-turns'] = NONE

def TCOCalc():
	"""
	TCO calculation
	"""	
	try:
		Cp = chain[res-1]['C'].get_vector()
		Op = chain[res-1]['O'].get_vector()
		p1 = C - O
		p2 = Cp - Op
		x = np.dot(p1,p1) * np.dot(p2,p2)
		if(x > 0):
			l['tco'] = np.dot(p1,p2) / np.sqrt(x)
	except:
		l['tco'] = 0

def kappaCalc():
	"""
	Kappa calculation
	"""	
	try:
		CApp = chain[res-2]['CA'].get_vector()
		CAnn = chain[res+2]['CA'].get_vector()
		l['kappa'] = 0.0
	except:
		l['kappa'] = 360
			#TCO

def alphaCalc():
	"""
	Alpha calculation (dihedral angle)
	"""
	try:
		CAp = chain[res-1]['CA'].get_vector()
		CAn = chain[res+1]['CA'].get_vector()
		CAnn = chain[res+2]['CA'].get_vector()
		l['alpha'] = (calc_dihedral(CAp,CA,CAn,CAnn)*180)/np.pi
	except:
		l['alpha'] = 360

	if (l['alpha'] < 0):
		l['chirality'] = '-'
	elif (l['alpha'] > 0 and l['alpha'] != 360):
		l['chirality'] = '+'
	else:
		l['chirality'] = ' '

def psiCalc():
	"""
	Phi calculation (dihedral angle)
	"""
	try:
		Nn = chain[res+1]['N'].get_vector()
		l['psi'] = (calc_dihedral(N, CA, C, Nn)*180)/np.pi
	except:
		l['psi'] = 360

def phiCalc():
	"""
	Phi calculation (dihedral angle)
	"""
	try:
		Cp = chain[res-1]['C'].get_vector() 
		l['phi'] = (calc_dihedral(Cp, N, CA, C)*180)/np.pi # degree = (radian*180)/pi
		
	except:
		l['phi'] = 360

if __name__ == "__main__":
	opt = argsParsing()
	hydrAddition(opt)

	p = PDBParser()
	structure = p.get_structure(opt.input,opt.input+".H")

	dssp = []
	index,	nb_chains = 0, 0

	for chain in structure.get_chains():
		nb_chains += 1
		# first residue id number of the current chain 
		first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))
		# last residue id number of the current chain 
		for res in chain.get_residues():
			if (res.get_id()[0] == 'W'):
				break
			last = res.get_id()[1]

		for res in range(first,last):
			index += 1
			l = {}
			residueDesc()

			N = chain[res]['N'].get_vector() 
			CA = chain[res]['CA'].get_vector() 
			C = chain[res]['C'].get_vector()
			O = chain[res]['O'].get_vector()

			nTurnPattern()
			TCOCalc()
			kappaCalc()
			alphaCalc()
			psiCalc()
			phiCalc()

			dssp.append(l)

	for i in range(0,len(dssp)-4):
		if ((dssp[i]['4-turns'] == '>' or dssp[i]['4-turns'] == 'X') and \
		((dssp[i+1]['4-turns'] == '>' or dssp[i+1]['4-turns'] == 'X') or (dssp[i-1]['4-turns'] == '>' or dssp[i-1]['4-turns'] == 'X'))):
			if (dssp[i-4]['4-turns'] == '>' or dssp[i-4]['4-turns'] == 'X'):
				dssp[i]['4-turns'] = 'X'
			if (dssp[i+2]['4-turns'] == ' ' and dssp[i+3]['4-turns'] == ' ' and dssp[i+4]['4-turns'] == ' ' and dssp[i+5]['4-turns'] == ' '):
				dssp[i+2]['4-turns'] = '<'
				dssp[i+3]['4-turns'] = '<'
				dssp[i+4]['4-turns'] = '<'
				dssp[i+5]['4-turns'] = '<'
	
	displayResults(opt,structure,dssp)