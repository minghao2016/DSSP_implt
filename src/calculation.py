import numpy as np
from Bio.PDB import *

def TCOCalc(chain,res,C,O):
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
			tco = np.dot(p1,p2) / np.sqrt(x)
	except:
		tco = 0
	return(tco)

def kappaCalc(chain,res):
	"""
	Kappa calculation
	"""	
	try:
		CApp = chain[res-2]['CA'].get_vector()
		CAnn = chain[res+2]['CA'].get_vector()
		kappa = 0.0
	except:
		kappa = 360
			#TCO
	return(kappa)

def alphaCalc(chain,res,CA):
	"""
	Alpha calculation (dihedral angle)
	"""
	try:
		CAp = chain[res-1]['CA'].get_vector()
		CAn = chain[res+1]['CA'].get_vector()
		CAnn = chain[res+2]['CA'].get_vector()
		alpha = (calc_dihedral(CAp,CA,CAn,CAnn)*180)/np.pi
	except:
		alpha = 360
	return(alpha)

def chirality(alpha):
	"""
	Chirality assignation
	"""
	if (alpha < 0):
		chirality = '-'
	elif (alpha > 0 and alpha != 360):
		chirality = '+'
	else:
		chirality= ' '
	return(chirality)

def psiCalc(chain,res,N,CA,C):
	"""
	Phi calculation (dihedral angle)
	"""
	try:
		Nn = chain[res+1]['N'].get_vector()
		psi = (calc_dihedral(N, CA, C, Nn)*180)/np.pi
	except:
		psi = 360
	return(psi)

def phiCalc(chain,res,N,CA,C):
	"""
	Phi calculation (dihedral angle)
	"""
	try:
		Cp = chain[res-1]['C'].get_vector() 
		phi = (calc_dihedral(Cp, N, CA, C)*180)/np.pi # degree = (radian*180)/pi
		
	except:
		phi = 360
	return(phi)

def testHbond(chain,res,n):
	"""
	Test if there is an H-bond between two a residue i
	and a residue i+n (n = [3,4,5])
	"""
	# Fixed variables :
	Q1, Q2 = 0.41, 0.20 # partial charges
	F = 332 # factor
	try:
		# Residue i
		O = chain[res]['O'] # Oxygen
		C = chain[res]['C'] # Carbon
		# Residue i+n
		Nn = chain[res+n]['N'] # Azote
		Hn = chain[res+n]['H'] # Hydrogen

		# r(ABn) = Interatomic distance from atom A (residue i) to B (residue i+n)
		rONn = Nn - O
		rCHn = Hn - C
		rOHn = Hn - O
		rCNn = Nn - C

		# Electrostatic interaction energy between two H-bonding groups.
		# E in kcal/mol
		E = Q1 * Q2 * (1/rONn + 1/rCHn - 1/rOHn - 1/rCNn) * F
		if (E < -0.5):
			return(True)
		else:
			return(False)
	except:
		# There is no residue i+n
		return(False)