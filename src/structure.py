"""@package structure
Documentation for this module.

More details.
"""
from Bio.PDB import *
# helix structural patterns

STRUC = {'none':' ','3-helix':'G','4-helix':'H','5-helix':'I','strand':'E'}

# n-turn patterns
NTURN = {'none':' ','start':'>','middle':{3:'3',4:'4',5:'5'},'end':'<','start end':'X'}

ABC_LOWER, ABC_UPPER = [], []
for letter in range(65,91):
    ABC_UPPER.append(chr(letter))
    ABC_LOWER.append(chr(letter+32))

def testHbond(chain,i,j):
    """    Test if there is an H-bond between residue i and j.

    Parameters
    ----------
    chain : Bio.PDB.Chain.Chain
        Given chain from the pdb file.
    i : int
        Residue number.
    j : int
        Residue number.

    Returns
    -------
    Boolean
        True if there is an H-bond, False if not.
    """
    # Fixed variables :
    Q1, Q2 = 0.41, 0.20 # partial charges
    F = 332 # factor
    try:
        # Residue i
        O = chain[i]['O'] # Oxygen
        C = chain[i]['C'] # Carbon
        # Residue i+n
        Nn = chain[j]['N'] # Azote
        Hn = chain[j]['H'] # Hydrogen

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

def isHelix(chain,res,n):
    """
    n-helix(i,i+n-1) = [n-turn(i-1) and n-turn(i)]
    """
    try:
        if (testHbond(chain,res,res+n) == True and testHbond(chain,res+1,res+1+n) == True):
            return(True)
        return(False)
    except:
        return(False)

def setHelicesStructure(dssp):
    for n in [5,4,3]:
        for res in range(len(dssp)):
            if (dssp[res][str(n)+'-turn']['rlt'] == NTURN['start']):
                while (dssp[res+2][str(n)+'-turn']['rlt'] != NTURN['none']):
                    res += 1
                    if (dssp[res]['structure'] == STRUC['none']):
                        dssp[res]['structure'] = STRUC[str(n)+'-helix']
            
    return(dssp)

def foundHelices(dssp,chain):
    """
    Helix patterns :
    3-helix : >>3<<
    4-helix : >>44<<
    5-helix : >>555<<
    """
    for res in range(len(dssp)):
        # n-Helix, with n = [3,4,5]
        for n in [5,4,3]:
            if (isHelix(chain,res,n)):
                for i in range(0,2):
                    # helix start : res(i):'>' & res(i+1):'>'
                    dssp[res+i][str(n)+'-turn']['start'] = NTURN['start']
                    # helix end : res(i+n):'<' & res(i+n+1):'<'
                    dssp[res+n+i][str(n)+'-turn']['end'] =  NTURN['end']
                for i in range(2,n):
                    # helix middle : res(i+2):'[3,4,5]' & res(i+3):'[4,5]' & res(i+4):'[5]'
                    dssp[res+i][str(n)+'-turn']['middle'] =  NTURN['middle'][n]

            if (dssp[res][str(n)+'-turn']['end'] == NTURN['end']):
                if (dssp[res][str(n)+'-turn']['start'] == NTURN['start']):
                    # '>' + '[3,4,5]' + '<' = 'X' (start & end)
                    dssp[res][str(n)+'-turn']['rlt'] = NTURN['start end']
                else:
                    # ' ' + '<' = '<' (end)
                    dssp[res][str(n)+'-turn']['rlt'] = NTURN['end']
            else:
                if (dssp[res][str(n)+'-turn']['start'] == NTURN['start']):
                    # '>' + ' ' = '>' (start)
                    dssp[res][str(n)+'-turn']['rlt'] = NTURN['start']
                else:
                    if (dssp[res][str(n)+'-turn']['middle'] != NTURN['none']):
                        # (' ' + ' ') and (middle != ' ') =  [3,4,5] (middle)
                        dssp[res][str(n)+'-turn']['rlt'] = NTURN['middle'][n]

    setHelicesStructure(dssp)
    return(dssp)

def isParallelBridge(chain,i,j):
    try:
        if ((testHbond(chain,i-1,j) == True and testHbond(chain,j,i+1) == True)\
        or (testHbond(chain,j-1,i) == True and testHbond(chain,i,j+1) == True)):
            return(True)
        return(False)
    except:
        return(False)

def isAntiparallelBridge(chain,i,j):
    try:
        if ((testHbond(chain,i,j) == True and testHbond(chain,j,i) == True)\
        or (testHbond(chain,i-1,j+1) == True and testHbond(chain,j-1,i+1) == True)):
            return(True)
        return(False)
    except:
        return(False)

def setStrandsStructure(dssp):
    for res in range(len(dssp)):
        if (dssp[res]['bridge']['pi']+dssp[res]['bridge']['pj'] != ''):
            dssp[res]['structure'] = STRUC['strand']
    return(dssp)

def foundStrands(dssp,chain):
    n, N = -1, -1
    for res_i in range(len(dssp)):
        if (dssp[res_i]['structure'] != STRUC['none']):
            continue
        for res_j in range(res_i+2,len(dssp)+1):
                if (isParallelBridge(chain,res_i,res_j) == True):
                    if (dssp[res_i-2]['bridge']['pi'] == ''):
                        n += 1
                        if (n == 26): n = 0
                    dssp[res_i-1]['bridge']['pi'] += ABC_LOWER[n]
                    dssp[res_i-1]['BP1'] = res_j 
                    dssp[res_j-1]['bridge']['pj'] += ABC_LOWER[n]
                    dssp[res_j-1]['BP2'] = res_i 
    setStrandsStructure(dssp)                    
    return(dssp)

"""
                elif (isAntiparallelBridge(chain,res_i,res_j) == True):
                    if (dssp[res_i-2]['bridge']['api'] == ''):
                        N += 1
                        if (N == 26): N = 0
                    dssp[res_i-1]['bridge']['api'] += ABC_UPPER[N]
                    dssp[res_i-1]['BP2'] = res_j
                    dssp[res_j-1]['bridge']['apj'] += ABC_UPPER[N]
                    dssp[res_j-1]['BP1'] = res_i 
"""
