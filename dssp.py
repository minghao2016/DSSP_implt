#!/home/helene/anaconda3/bin/python3.6

from Bio.SeqUtils import seq1
from Bio.PDB import *

import sys
import os
srcPath = 'src'
sys.path.append(os.path.abspath(srcPath))

from functions import *
from calc import *

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
		if (testHbond(chain,res-1,n) == True):
			l[str(n)+'-turns'] = START
		else:
			l[str(n)+'-turns'] = NONE

def startHelixPattern():
	for n in range(3,6):
		if (isHelix(n) == True):
			l[str(n)+'-turns'] = START
			try:
				if (dssp[res-1][str(n)+'-turns'] == NONE):
					dssp[res-1][str(n)+'-turns'] = START
			except:
				print('no')
		else:
			l[str(n)+'-turns'] = NONE

def isHelix(n):
	try:
		if (testHbond(chain,res-1,n) == True and testHbond(chain,res,n) == True):
			return(True)
		return(False)
	except:
		return(False)

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


			startHelixPattern()

			l['tco'] = TCOCalc(chain,res)
			l['kappa'] = kappaCalc(chain,res)
			l['alpha'] = alphaCalc(chain,res)
			l['chirality'] = chirality(l['alpha'])
			l['psi'] = psiCalc(chain,res)
			l['phi'] = phiCalc(chain,res)

			dssp.append(l)
	displayResults(opt,structure,dssp)


	for n in range(3,6):
		for i in range(0,len(dssp)-4):
			if (dssp[i-4][str(n)+'-turns'] == '>' or dssp[i-4][str(n)+'-turns'] == 'X'):
				if (dssp[i][str(n)+'-turns'] == '>'):
					dssp[i][str(n)+'-turns'] = 'X'
			#elif (dssp[i-4][str(n)+'-turns'] == '>' or dssp[i-4][str(n)+'-turns'] == 'X'):

	
"""
	for i in range(0,len(dssp)-4):
		if (dssp[i-4]['4-turns'] == '>' or dssp[i-4]['4-turns'] == 'X'):
			if (dssp[i]['4-turns'] == '>'):
				dssp[i]['4-turns'] = 'X'
			elif  (dssp[i]['4-turns'] == ' '):
				for j in range(i,i+4):
					dssp[j]['4-turns'] = '<'
			
			


	displayResults(opt,structure,dssp)
"""
'''
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
'''
	
	