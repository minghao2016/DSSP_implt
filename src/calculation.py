"""@package calculation

Small functions to generate variables that describe a residue.
"""

import numpy as np
from Bio.PDB import *

def TCOCalc(chain,res):
    """TCO calculation for a given residue.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        Given chain from the pdb file.
    res : int
        Residue number.

    Returns
    -------
    tco : int/float
        Cosine of angle TCO of the given residue.
    """
    try:
        C = chain[res]['C'].get_vector()
        O = chain[res]['O'].get_vector()
        Cp = chain[res-1]['C'].get_vector()
        Op = chain[res-1]['O'].get_vector()
        p1 = C - O
        p2 = Cp - Op
        x = np.dot(p1,p1) * np.dot(p2,p2)
        if(x > 0):
            tco = np.dot(p1,p2) / np.sqrt(x)
    except:
        tco = 0
    return(tco)

def kappaCalc(chain,res):
    """KAPPA calculation for a given residue

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        Given chain from the pdb file
    res : int
        Residue number

    Returns
    -------
    kappa : int/float
        KAPPA angle of the given residue
    """
    try:
        CApp = chain[res-2]['CA'].get_vector()
        CAnn = chain[res+2]['CA'].get_vector()
        kappa = 0.0
    except:
        kappa = 360
    return(kappa)

def alphaCalc(chain,res):
    """Alpha calculation for a given residue.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        Given chain from the pdb file.
    res : int
        Residue number.

    Returns
    -------
    alpha : int/float
        Dihedral angle alpha of the given residue.
    """
    try:
        CA =  chain[res]['CA'].get_vector()
        CAp = chain[res-1]['CA'].get_vector()
        CAn = chain[res+1]['CA'].get_vector()
        CAnn = chain[res+2]['CA'].get_vector()
        alpha = (calc_dihedral(CAp,CA,CAn,CAnn)*180)/np.pi
    except:
        alpha = 360
    return(alpha)

def chirality(alpha):
    """Chirality assignation of a residue

    Parameters
    ----------
    alpha : int/float
        Alpha dihedral angle of a given residue.

    Returns
    -------
    chirality : str
        '+' , '-' or ' '.
    """
    if (alpha < 0):
        chirality = '-'
    elif (alpha > 0 and alpha != 360):
        chirality = '+'
    else:
        chirality = ' '
    return(chirality)

def psiCalc(chain,res):
    """Psi calculation for a given residue.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        Given chain from the pdb file.
    res : int
        Residue number.

    Returns
    -------
    psi : int/float
        Dihedral angle psi of the given residue.
    """
    try:
        N = chain[res]['N'].get_vector()
        CA = chain[res]['CA'].get_vector()
        C = chain[res]['C'].get_vector()
        Nn = chain[res+1]['N'].get_vector()
        # degree = (radian*180)/pi
        psi = (calc_dihedral(N, CA, C, Nn)*180) / np.pi
    except:
        psi = 360
    return(psi)

def phiCalc(chain,res):
    """Phi calculation for a given residue.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        Given chain from the pdb file.
    res : int
        Residue number.

    Returns
    -------
    phi : int/float
        Dihedral angle phi of the given residue.
    """
    try:
        N = chain[res]['N'].get_vector()
        CA = chain[res]['CA'].get_vector()
        C = chain[res]['C'].get_vector()
        Cp = chain[res-1]['C'].get_vector()
        phi = (calc_dihedral(Cp, N, CA, C)*180) / np.pi
    except:
        phi = 360
    return(phi)
