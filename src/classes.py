"""@package classes

Residue class and its methods to calculate supplement parameters.
"""
import numpy as np
# Biopython :
from Bio.PDB import *
from Bio.SeqUtils import seq1

NONE = ' '

class Nturn:
    """Class that regroups patterns of an nturn"""
    def __init__(self):
        """The constructor of a Residue class."""
        self.start = NONE
        self.middle = NONE
        self.end = NONE
        self.result = NONE

class Residue:
    """Class that groups informations about a residue."""
    def __init__(self, chain, index, chainID, resNum):
        """The constructor of a Residue class."""
        self.index = index
        self.chainID = chainID
        self.resNum = resNum
        self.resName = seq1(chain[resNum].get_resname())
        self.structure = NONE # 'H', 'B', 'E', 'G', 'I', 'T', 'S' or ' '
        self.nturns = { 3: Nturn(), 4: Nturn(), 5: Nturn() }
        self.bend = NONE
        self.chirality = NONE 
        self.bridge_1 = NONE
        self.bridge_2 = NONE
        self.bp1 = 0
        self.bp2 = 0
        self.sheet = NONE
        self.tco = 0
        self.kappa = 360
        self.alpha = 360
        self.phi = 360
        self.psi = 360
        self.CA = chain[self.resNum]['CA'].get_vector()
        self.C = chain[self.resNum]['C']
        self.O = chain[self.resNum]['O']
        self.N = chain[self.resNum]['N']
        try:
            self.H = chain[self.resNum]['H']
        except:
            pass

    def tco_calculation(self,chain):
        """Cosine of angle TCO of the residue."""
        try:
            Cp = chain[self.resNum-1]['C'].get_vector()
            Op = chain[self.resNum-1]['O'].get_vector()
            p1 = self.C.get_vector() - self.O.get_vector()
            p2 = Cp - Op
            x = np.dot(p1,p1) * np.dot(p2,p2)
            if(x > 0):
                self.tco = np.dot(p1,p2) / np.sqrt(x)
        except:
            pass

    def kappa_calculation(self,chain):
        """KAPPA angle of the residue."""
        try:
            ppCA = chain[self.resNum-2]['CA'].get_vector()
            nnCA = chain[self.resNum+2]['CA'].get_vector()
            p1 = self.CA - ppCA
            p2 = nnCA - self.CA
            x = np.dot(p1,p1) * np.dot(p2,p2)
            if (x > 0):
                a = np.dot(p1,p2) / np.sqrt(x)
                b = np.sqrt(1 - a * a)
                self.kappa = np.arctan2(b,a) * (180 / (4 * np.arctan(1.0)))
        except:
            pass
    def bend_assignation(self):
        if (self.kappa != 360 and self.kappa > 70):
            self.bend = 'S'

    def alpha_calculation(self,chain):
        """Dihedral angle alpha of the residue."""
        try:
            CAp = chain[self.resNum-1]['CA'].get_vector()
            CAn = chain[self.resNum+1]['CA'].get_vector()
            CAnn = chain[self.resNum+2]['CA'].get_vector()
            self.alpha = (calc_dihedral(CAp,self.CA,CAn,CAnn)*180) / np.pi
        except:
            pass

    def chirality_assignation(self):
        """Chirality assignation of the residue."""
        if (self.alpha < 0):
            self.chirality = '-'
        elif (self.alpha > 0 and self.alpha != 360):
            self.chirality = '+'

    def psi_calculation(self,chain):
        """Dihedral angle psi of the residue."""
        try:
            Nn = chain[self.resNum+1]['N'].get_vector()
            # degree = (radian*180)/pi
            self.psi = (calc_dihedral(self.N.get_vector(), self.CA, self.C.get_vector(), Nn)*180) / np.pi
        except:
            pass

    def phi_calculation(self,chain):
        """Dihedral angle phi of the residue."""
        try:
            Cp = chain[self.resNum-1]['C'].get_vector()
            # degree = (radian*180)/pi
            self.phi = (calc_dihedral(Cp, self.N.get_vector(), self.CA, self.C.get_vector())*180) / np.pi
        except:
            pass

