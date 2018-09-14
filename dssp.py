#!/home/helene/anaconda3/bin/python3.6

__author__ = "Hélène Kabbech"

import sys
import os

# Implemented modules in /src foler :
sys.path.append(os.path.abspath('src'))
from calculation import *
from structure import *
from additional_functions import *

# Biopython modules used :
from Bio.SeqUtils import seq1
from Bio.PDB import *


NONE = ' '

if __name__ == "__main__":
    opt = argsParsing()
    hydrAddition(opt)

    p = PDBParser(QUIET=True) # QUIET=T : Warnings issued are suppressed
    pdb = p.get_structure(opt.input,opt.input+".H")

    dssp = []
    index,    nb_chains = 0, 0

    for chain in pdb.get_chains():
        nb_chains += 1
        # first residue id number of the current chain
        first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))
        # last residue id number of the current chain
        for res in chain.get_residues():
            if (res.get_id()[0] != NONE):
                break
            last = res.get_id()[1]

        for res in range(first,last+1):
            index += 1
            l = {}
            l['index'] = index
            l['res_nb'] = res
            l['chain'] = chain.id
            l['aa'] = seq1(chain[res].get_resname())
            l['structure'] = NONE
            for n in [5,4,3]:
                l[str(n)+'-turn'] = {'start':NONE,'middle':NONE,'end':NONE,'rlt':NONE}
            l['bridge'] = {'pi':'','pj':'','api':'','apj':''}
            l['BP1'], l['BP2'] = 0, 0
            l['sheet'] = NONE
            l['tco'] = TCOCalc(chain,res)
            l['kappa'] = kappaCalc(chain,res)
            l['bend'] = NONE
            l['alpha'] = alphaCalc(chain,res)
            l['chirality'] = chirality(l['alpha'])
            l['psi'] = psiCalc(chain,res)
            l['phi'] = phiCalc(chain,res)
            l['x-ca'] = chain[res]['CA'].get_coord()[0]
            l['y-ca'] = chain[res]['CA'].get_coord()[1]
            l['z-ca'] = chain[res]['CA'].get_coord()[2]
            dssp.append(l)

    dssp = foundHelices(dssp,chain)
    dssp = foundStrands(dssp,chain)
    displayResults(opt,pdb,dssp)