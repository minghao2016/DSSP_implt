#!/usr/bin/python3.6

__author__ = "Hélène Kabbech"

import sys
import os

# Implemented modules in /src foler :
sys.path.append(os.path.abspath('src'))
from structure import *
from additional_functions import *
import classes

# Biopython :
from Bio.PDB import *

if __name__ == "__main__":
    opt = argsParsing() # Input and Output

    hydrAddition(opt.input) # Hydrogen addition (*.pdb.H file created)

    p = PDBParser(QUIET=True) # QUIET=T : Warnings issued are suppressed
    pdb = p.get_structure(opt.input,opt.input+".H")

    resList = [] # List of Residue instances
    index = 0

    for chain in pdb.get_chains():
        # First residue of the current chain
        first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))

        # Last residue of the current chain
        for res in chain.get_residues():
            if (res.get_id()[0] != ' '):
                break
            last = res.get_id()[1]

        for resNum in range(first,last+1):
            index += 1
            r = classes.Residue(chain, index, chain.id, resNum)
            r.tco_calculation(chain)
            r.kappa_calculation(chain)
            r.bend_assignation()
            r.alpha_calculation(chain)
            r.chirality_assignation()
            r.phi_calculation(chain)
            r.psi_calculation(chain)
            resList.append(r)

    resList = setSSE(resList)
    displayResults(opt,pdb,resList)