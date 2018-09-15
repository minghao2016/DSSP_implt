"""@package structure
Documentation for this module.

More details.
"""

NONE = ' '

# n-turn patterns
START = '>'
END = '<'
START_END = 'X'
MIDDLE = { 3: '3', 4: '4', 5: '5' }

# Structural patterns
HELIX = { 3: 'G', 4: 'H', 5: 'I' }
STRAND = 'E'




ABC_LOWER, ABC_UPPER = [], []
for letter in range(65,91):
    ABC_UPPER.append(chr(letter))
    ABC_LOWER.append(chr(letter+32))

def testHbond(ri,rj):
    """Test if there is an H-bond between residue i and j. Return a boolean."""
    # Fixed variables :
    Q1, Q2 = 0.41, 0.20 # partial charges
    F = 332 # factor
    try:
        # dAB = Interatomic distance from atom A (residue i) to B (residue j)
        dON = rj.N - ri.O
        dCH = rj.H - ri.C
        dOH = rj.H - ri.O
        dCN = rj.N - ri.C

        # Electrostatic interaction energy between two H-bonding groups.
        # E in kcal/mol
        E = Q1 * Q2 * (1/dON + 1/dCH - 1/dOH - 1/dCN) * F
        if (E < -0.5):
            return(True)
        return(False)
    except:
        # There is no residue i+n
        return(False)

def isHelix(dssp,res,n):
    """n-helix(i,i+n-1) = [n-turn(i-1) and n-turn(i)]"""
    try:
        if (testHbond(dssp[res],dssp[res+n]) == True and testHbond(dssp[res+1],dssp[res+1+n]) == True):
            return(True)
        return(False)
    except:
        return(False)

def setHelicesStructure(dssp):
    """ """
    for n in [5,4,3]:
        for res in range(len(dssp)):
            if (dssp[res].nturns[n].result == START):
                while (dssp[res+2].nturns[n].result != NONE):
                    res += 1
                    if (dssp[res].structure == NONE):
                        dssp[res].structure = HELIX[n]
    return(dssp)

def foundHelices(dssp):
    """
    Helix patterns :
    3-helix : >>3<<
    4-helix : >>44<<
    5-helix : >>555<<
    """
    for res in range(len(dssp)):
        # n-Helix, with n = [3,4,5]
        for n in [5,4,3]:
            if (isHelix(dssp,res,n)):
                for i in range(0,2):
                    # helix start : res(i):'>' & res(i+1):'>'
                    dssp[res+i].nturns[n].start = START
                    # helix end : res(i+n):'<' & res(i+n+1):'<'
                    dssp[res+n+i].nturns[n].end =  END
                for i in range(2,n):
                    # helix middle : res(i+2):'[3,4,5]' & res(i+3):'[4,5]' & res(i+4):'[5]'
                    dssp[res+i].nturns[n].middle =  MIDDLE[n]

            if (dssp[res].nturns[n].end == END):
                if (dssp[res].nturns[n].start == START):
                    # '>' + '[3,4,5]' + '<' = 'X' (start & end)
                    dssp[res].nturns[n].result = START_END
                else:
                    # ' ' + '<' = '<' (end)
                    dssp[res].nturns[n].result = END
            else:
                if (dssp[res].nturns[n].start == START):
                    # '>' + ' ' = '>' (start)
                    dssp[res].nturns[n].result = START
                else:
                    if (dssp[res].nturns[n].middle != NONE):
                        # (' ' + ' ') and (middle != ' ') =  [3,4,5] (middle)
                        dssp[res].nturns[n].result = MIDDLE[n]

    setHelicesStructure(dssp)
    return(dssp)

def isParallelBridge(i,j):
    try:
        if ((testHbond(i-1,j) == True and testHbond(j,i+1) == True)\
        or (testHbond(j-1,i) == True and testHbond(i,j+1) == True)):
            return(True)
        return(False)
    except:
        return(False)

def isAntiparallelBridge(i,j):
    try:
        if ((testHbond(i,j) == True and testHbond(j,i) == True)\
        or (testHbond(i-1,j+1) == True and testHbond(j-1,i+1) == True)):
            return(True)
        return(False)
    except:
        return(False)

def setStrandsStructure(dssp):
    for res in range(1,len(dssp)-1):
        if ((dssp[res-1].bridge_1 != NONE or dssp[res+1].bridge_1 != NONE) \
            and dssp[res].bridge_1!= NONE):
            dssp[res].structure = STRAND
    return(dssp)

def foundStrands(dssp):
    n = -1
    s = 0
    left = True
    for res_i in range(len(dssp)):
        if (dssp[res_i].structure != NONE):
            continue
        for res_j in range(res_i+2,len(dssp)+1):
            if (isParallelBridge(res_i,res_j) == True):
                #dssp[res_i-1]['sheet'] = ABC_UPPER[s]
                #dssp[res_j-1]['sheet'] = ABC_UPPER[s]
                if (dssp[res_i-2].bridge_1 == NONE):
                    n += 1
                    if (n == 26): n = 0
                if (dssp[res_i-1].bridge_1 == NONE):
                    dssp[res_i-1].bridge_1 = ABC_LOWER[n]
                    dssp[res_i-1].bp1 = res_j
                else:
                    dssp[res_i-1].bridge_2 = ABC_LOWER[n]
                    dssp[res_i-1].bp2 = res_j

                if (dssp[res_j-1].bridge_1 == NONE):
                    dssp[res_j-1].bridge_1 = ABC_LOWER[n]
                    dssp[res_j-1].bp1 = res_i
                else:
                    dssp[res_j-1].bridge_2 = ABC_LOWER[n]
                    dssp[res_j-1].bp2 = res_i
                if (dssp[res_i-1].bridge_2 != NONE):
                    s += 1
    setStrandsStructure(dssp)                    
    return(dssp)
"""
def foundStrands(dssp,chain):
    n = -1
    s = 0
    for res_i in range(len(dssp)):
        if (dssp[res_i]['structure'] != STRUC['none']):
            continue
        for res_j in range(res_i+2,len(dssp)+1):
                if (isParallelBridge(chain,res_i,res_j) == True):
                    #dssp[res_i-1]['sheet'] = ABC_UPPER[s]
                    #dssp[res_j-1]['sheet'] = ABC_UPPER[s]
                    if (dssp[res_i-2]['bridge 1'] == NONE):
                        n += 1
                        if (n == 26): n = 0
                    if (dssp[res_i-1]['bridge 1'] == NONE):
                        dssp[res_i-1]['bridge 1'] = ABC_LOWER[n]
                        dssp[res_i-1]['bp1'] = res_j
                    else:
                        dssp[res_i-1]['bridge 2'] = ABC_LOWER[n]
                        dssp[res_i-1]['bp2'] = res_j

                    if (dssp[res_j-1]['bridge 1'] == NONE):
                        dssp[res_j-1]['bridge 1'] = ABC_LOWER[n]
                        dssp[res_j-1]['bp1'] = res_i
                    else:
                        dssp[res_j-1]['bridge 2'] = ABC_LOWER[n]
                        dssp[res_j-1]['bp2'] = res_i
                    if (dssp[res_i-1]['bridge 2'] != NONE):
                        s += 1
    setStrandsStructure(dssp)                    
    return(dssp)
"""

"""
                elif (isAntiparallelBridge(chain,res_i,res_j) == True):
                    if (dssp[res_i-2]['bridge']['api'] == ''):
                        N += 1
                        if (N == 26): N = 0
                    dssp[res_i-1]['bridge']['api'] += ABC_UPPER[N]
                    dssp[res_i-1]['bp2'] = res_j
                    dssp[res_j-1]['bridge']['apj'] += ABC_UPPER[N]
                    dssp[res_j-1]['bp1'] = res_i 
"""
