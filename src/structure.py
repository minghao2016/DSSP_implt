"""@package structure
Documentation for this module.

Secondary structure elements
Helices
Strands
n-Turns
Bridges
"""

import classes

################ VARIABLES ################
NONE = ' '

# n-turn patterns:
START = '>'
END = '<'
START_END = 'X'
MIDDLE = { 3: '3', 4: '4', 5: '5' }

# Secondary structural patterns :
HELIX = { 3: 'G', 4: 'H', 5: 'I' }
STRAND = 'E'
TURN = 'T'
BRIDGE = 'B'
BEND = 'S'

# Alphabet :
ABC_LOWER, ABC_UPPER = [], []
for letter in range(65,91):
    ABC_UPPER.append(chr(letter))
    ABC_LOWER.append(chr(letter+32))

##########################################

def isHbond(ri,rj):
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

def isHelix(dssp,i,n):
    """n-helix(i,i+n-1) = [n-turn(i-1) and n-turn(i)]"""
    try:
        if (isHbond(dssp[i],dssp[i+n]) == True and isHbond(dssp[i+1],dssp[i+1+n]) == True):
            return(True)
        return(False)
    except:
        return(False)

def isParallelBridge(dssp,i,j):
    try:
        if ((isHbond(dssp[i-1],dssp[j]) == True and isHbond(dssp[j],dssp[i+1]) == True)\
        or (isHbond(dssp[j-1],dssp[i]) == True and isHbond(dssp[i],dssp[j+1]) == True)):
            return(True)
        return(False)
    except:
        return(False)

def isAntiparallelBridge(dssp,i,j):
    try:
        if ((isHbond(dssp[i],dssp[j]) == True and isHbond(dssp[j],dssp[i]) == True)\
        or (isHbond(dssp[i-1],dssp[j+1]) == True and isHbond(dssp[j-1],dssp[i+1]) == True)):
            return(True)
        return(False)
    except:
        return(False)

def setHelixStruct(dssp,n):
    """ """
    for res in range(len(dssp)):
        if (dssp[res].structure != NONE): continue
        if (dssp[res].nturns[n].result == START):
            while (dssp[res+2].nturns[n].result != NONE):
                res += 1
                if (dssp[res].structure == NONE):
                    dssp[res].structure = HELIX[n]
    return(dssp)

def setIsolatedBridgeStruct(dssp):
    for res in range(1,len(dssp)-1):
        if (dssp[res].structure != NONE): continue
        if ((dssp[res].bp1+dssp[res].bp2) != 0):
            if (dssp[res-1].bp1+dssp[res-1].bp2 == 0 and dssp[res+1].bp1+dssp[res+1].bp2 == 0):
                dssp[res].structure = BRIDGE
                dssp[dssp[res].bp1-1].structure = BRIDGE
    return(dssp)

def setStrandStruct(dssp):
    for res in range(1,len(dssp)-1):
        if (dssp[res].structure != NONE): continue
        if ((dssp[res].bp1+dssp[res].bp2) != 0):
            if (dssp[res-1].bp1+dssp[res-1].bp2 != 0 or dssp[res+1].bp1+dssp[res+1].bp2 != 0):
                dssp[res].structure = STRAND
        else:
            if (dssp[res-1].bp1+dssp[res-1].bp2 != 0 and dssp[res+1].bp1+dssp[res+1].bp2 != 0\
                and dssp[res-2].bp1+dssp[res-2].bp2 != 0 and dssp[res+2].bp1+dssp[res+2].bp2 != 0):
                dssp[res].structure = STRAND
    return(dssp)

def setNturnsStruct(dssp):
    """ """
    for n in [5,4,3]:
        for res in range(len(dssp)):
            if (dssp[res].nturns[n].result == START):
                while (dssp[res+1].structure != NONE): res += 1
                while (dssp[res+2].nturns[n].result != NONE):
                    res += 1
                    if (dssp[res].structure == NONE):
                        dssp[res].structure = TURN
    return(dssp)

def foundHelices(dssp,n):
    """
    Helix patterns :
    3-helix : >>3<<
    4-helix : >>44<<
    5-helix : >>555<<
    """
    for res in range(0,len(dssp)):
        if (dssp[res].structure != NONE): continue
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

    setHelixStruct(dssp,n)
    return(dssp)

def foundStrands(dssp):
    n = -1
    newStrand = True
    for res_i in range(len(dssp)):
        if (dssp[res_i].structure != NONE): continue
        for res_j in range(res_i+2,len(dssp)):
            if (dssp[res_j].structure != NONE): continue
            if (isParallelBridge(dssp,res_i,res_j) == True or isAntiparallelBridge(dssp,res_i,res_j) == True):
                if (isParallelBridge(dssp,res_i,res_j) == True):
                    alphabet = ABC_LOWER
                    if (dssp[res_i-1].bp1 == 0 or dssp[res_i-1].bp1+1 != res_j+1):
                        newStrand = True
                    else:
                        newStrand = False
                elif (isAntiparallelBridge(dssp,res_i,res_j) == True):
                    alphabet = ABC_UPPER
                    if (dssp[res_i-1].bp1 == 0 or dssp[res_i-1].bp1-1 != res_j+1):
                        newStrand = True
                    else:
                        newStrand = False
                if (dssp[res_i+1].bp1 == 0 and dssp[res_i-1].bp2 == 0):
                    if (newStrand == True):
                        n += 1
                        if (n == 26): n = 0
                    dssp[res_i].bp1 = res_j+1
                    dssp[res_i].bridge_1 = alphabet[n]
                else:
                    dssp[res_i].bp2 = res_j+1
                    if (dssp[res_i-1].bp2 == 0 and dssp[res_i-2].bp2 == 0): # New strand
                        n += 1
                        if (n == 26): n = 0
                    dssp[res_i].bridge_2 = alphabet[n]
                dssp[res_j].bp1 = res_i+1
                dssp[res_j].bridge_1 = alphabet[n]
    dssp=setIsolatedBridgeStruct(dssp)
    dssp=setStrandStruct(dssp)
    return(dssp)

def foundNturns(dssp):
    for res in range(len(dssp)-5):
        for n in [3,4,5]:
            if (isHbond(dssp[res],dssp[res+n]) == True):
                dssp[res].nturns[n].start = START
                dssp[res+n].nturns[n].end =  END
                for i in range(1,n):
                    dssp[res+i].nturns[n].middle =  MIDDLE[n]
    for res in range(len(dssp)):
        for n in [3,4,5]:
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
    setNturnsStruct(dssp)
    return(dssp)

def foundAndSetBend(dssp):
    for res in range(len(dssp)):
        if (dssp[res].kappa != 360 and dssp[res].kappa > 70):
            dssp[res].bend = BEND
            if (dssp[res].structure != NONE):
                continue
            else:
                dssp[res].structure = BEND
    return(dssp)

def setSSE(dssp):
    """ """
    foundHelices(dssp,4)
    foundStrands(dssp)
    foundHelices(dssp,3)
    foundHelices(dssp,5)
    foundNturns(dssp)
    foundAndSetBend(dssp)
    return(dssp)