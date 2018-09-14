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

__author__ = "Hélène Kabbech"

NONE = ' '
HELIX_STRUC = {'3':'G','4':'H','5':'I'}
HELICES = [5,4,3]

NONE = ' '
START = '>'
END = '<'
MIDDLE = {'3':'3','4':'4','5':'5'}
START_END = 'X'

alphabet = []
for letter in range(97,123):
    alphabet.append(chr(letter))


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
    l['tco'] = TCOCalc(chain,res)
    l['kappa'] = kappaCalc(chain,res)
    l['alpha'] = alphaCalc(chain,res)
    l['chirality'] = chirality(l['alpha'])
    l['psi'] = psiCalc(chain,res)
    l['phi'] = phiCalc(chain,res)
    l['bridge'] = {'i':'','j':''}
    for n in HELICES:
        l[str(n)+'-hbonds'] = testHbond(chain,res,res+n)
        l[str(n)+'-turns'] = {'start':NONE,'middle':NONE,'end':NONE,'rlt':NONE}

def isParallelBridge(dssp,i,j):
    try:
        if ((testHbond(chain,i-1,j) == True and testHbond(chain,j,i+1) == True)\
        or (testHbond(chain,j-1,i) == True and testHbond(chain,i,j+1) == True)):
            return(True)
        return(False)
    except:
        return(False)

def isAntiparallelBridge(dssp,i,j):
    try:
        if ((testHbond(chain,i,j) == True and testHbond(chain,j,i) == True)\
        or (testHbond(chain,i-1,j+1) == True and testHbond(chain,j-1,i+1) == True)):
            return(True)
        return(False)
    except:
        return(False)

if __name__ == "__main__":
    opt = argsParsing()
    hydrAddition(opt)

    p = PDBParser()
    pdb = p.get_structure(opt.input,opt.input+".H")

    dssp = []
    index,    nb_chains = 0, 0

    for chain in pdb.get_chains():
        nb_chains += 1
        # first residue id number of the current chain
        first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))
        # last residue id number of the current chain
        for res in chain.get_residues():
            if (res.get_id()[0] != ' '):
                break
            last = res.get_id()[1]

        for res in range(first,last+1):
            index += 1
            l = {}

            residueDesc()
            dssp.append(l)

    dssp = nTurnPatterns(dssp)


    a, A = -1, -1
    for res_i in range(len(dssp)):
        for res_j in range(res_i+2,len(dssp)+1):
                if (isParallelBridge(dssp,res_i,res_j) == True):
                    if (dssp[res_i-2]['bridge']['i'] == ''):
                        a += 1
                    dssp[res_i-1]['bridge']['i'] += alphabet[a]
                    dssp[res_j-1]['bridge']['j'] += alphabet[a]
                #elif (isAntiparallelBridge(dssp,res_i,res_j) == True):
                #    dssp[res_i-1]['bridge']['i'] += alphabet[a].upper()
                #    dssp[res_j-1]['bridge']['j'] += alphabet[a].upper()
                #    if (dssp[res_i-2]['bridge']['i'] == ''):
                #        A += 1
            
    dssp = setStructure(dssp)

    displayResults(opt,pdb,dssp)