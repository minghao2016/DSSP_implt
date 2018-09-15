#!/home/helene/anaconda3/bin/python3.6

__author__ = "Hélène Kabbech"

import sys
import os

# Biopython modules used :
from Bio.SeqUtils import seq1
from Bio.PDB import *

# Implemented modules in /src foler :
sys.path.append(os.path.abspath('src'))
#from calculation import *
from structure import *
from additional_functions import *
import classes

if __name__ == "__main__":
    opt = argsParsing()
    hydrAddition(opt)

    p = PDBParser(QUIET=True) # QUIET=T : Warnings issued are suppressed
    pdb = p.get_structure(opt.input,opt.input+".H")
    #pdb = p.get_structure("1BTA","data/1bta.pdb.H")

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
        for resNum in range(first,last+1):
            index += 1
            r = classes.Residue(chain, index, chain.id, resNum)
            r.tco_calculation(chain)
            r.kappa_calculation(chain)
            r.alpha_calculation(chain)
            r.chirality_assignation()
            r.phi_calculation(chain)
            r.psi_calculation(chain)
            dssp.append(r)

    dssp = foundHelices(dssp)
    dssp = foundStrands(dssp)
    displayResults(opt,pdb,dssp)