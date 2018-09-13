#!/home/helene/anaconda3/bin/python3.6

import sys
import os
# Implemented modules in /src foler :
sys.path.append(os.path.abspath('src'))
from calculation import *
from structure import *
from additional_functions import *

# Biopython
from Bio.SeqUtils import seq1
from Bio.PDB import *

NONE = ' '
HELIX_STRUC = {'3':'G','4':'H','5':'I'}
HELICES = [5,4,3]

NONE = ' '
START = '>'
END = '<'
MIDDLE = {'3':'3','4':'4','5':'5'}
START_END = 'X'


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
	l['tco'] = TCOCalc(chain,res,C,O)
	l['kappa'] = kappaCalc(chain,res)
	l['alpha'] = alphaCalc(chain,res,CA)
	l['chirality'] = chirality(l['alpha'])
	l['psi'] = psiCalc(chain,res,N,CA,C)
	l['phi'] = phiCalc(chain,res,N,CA,C)
	for n in HELICES:
		l[str(n)+'-hbonds'] = testHbond(chain,res,n)
		l[str(n)+'-turns'] = {'start':NONE,'middle':NONE,'end':NONE,'rlt':NONE}

if __name__ == "__main__":
	opt = argsParsing()
	hydrAddition(opt)

	p = PDBParser()
	pdb = p.get_structure(opt.input,opt.input+".H")

	dssp = []
	index,	nb_chains = 0, 0

	for chain in pdb.get_chains():
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

			N = chain[res]['N'].get_vector() 
			CA = chain[res]['CA'].get_vector() 
			C = chain[res]['C'].get_vector()
			O = chain[res]['O'].get_vector()
			
			residueDesc()
			dssp.append(l)

	dssp = nTurnPatterns(dssp)
	dssp = setStructure(dssp)

	displayResults(opt,pdb,dssp)
	
	"""
	for i in range(0,len(dssp)):
		l = dssp[i]
		print("{:>5d}{:>5d}{:>2s}{:>2s}{:>2s}{:>3s}{:>1s}{:>1s}{:>2s}"\
		.format(l['index'],l['res_nb'],l['chain'],l['aa'],l['structure'],l['4-turns']['rlt'],l['4-turns']['start'],l['4-turns']['middle'],l['4-turns']['end'],\
		l['chirality']))
	"""
	