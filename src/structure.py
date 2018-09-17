"""@package structure

Secondary structures assignment
    1. 4-helices (H)
    2. isolated beta-bridge (B)
    3. beta-strands (E)
    4. 3-helices (G)
    5. 5-helices (I)
    6. n-turns (T)
    7.bends (S)
"""

############################## V A R I A B L E S ###############################
NONE = ' '

# n-turn patterns:
START = '>'
END = '<'
START_END = 'X'
MIDDLE = { 3: '3', 4: '4', 5: '5' }

# Secondary structure patterns :
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

################################ H - B O N D S #################################

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
        return(False)

######################## N T U R N S  &  H E L I C E S #########################

def isHelix(resList,i,n):
    """Test if residue i and i+1 are each bonded to an Hydrogen.
    -> Residues i to i+1+n form an n-helix. Return a boolean."""
    try:
        if (isHbond(resList[i],resList[i+n]) == True and isHbond(resList[i+1],resList[i+1+n]) == True):
            return(True)
        return(False)
    except:
        return(False)

def setNturnPatternResult(resList,res,n):
    """n-turn pattern determination of a residue"""
    if (resList[res].nturns[n].end == END):
        # '>' + '[3,4,5]' + '<' = 'X' (start & end)
        if (resList[res].nturns[n].start == START):
            resList[res].nturns[n].result = START_END
        # ' ' + '<' = '<' (end)
        else:
            resList[res].nturns[n].result = END
    else:
        # '>' + ' ' = '>' (start)
        if (resList[res].nturns[n].start == START):
            resList[res].nturns[n].result = START
        else:
            # (' ' + ' ') and (middle != ' ') =  [3,4,5] (middle)
            if (resList[res].nturns[n].middle != NONE):
                resList[res].nturns[n].result = MIDDLE[n]
    return(resList)

def setHelixStruct(resList,n):
    """Assignment of n-helix to all concerned residues"""
    for res in range(len(resList)):
        # Secondary structure already determined. There is no overwrite
        if (resList[res].structure != NONE): continue
        if (resList[res].nturns[n].result == START):
            # Start and end patterns are not include in the helix structure
            while (resList[res+2].nturns[n].result != NONE):
                res += 1
                if (resList[res].structure == NONE): # no overwrite
                    resList[res].structure = HELIX[n]
    return(resList)

def setNturnsStruct(resList):
    """Assignment of turns (3, 4 and 5 types) to all concerned residues"""
    for n in [5,4,3]:
        for res in range(len(resList)):
            # Start and end patterns are not include
            if (resList[res].nturns[n].result == START):
                while (resList[res+1].structure != NONE): res += 1 # no overwrite
                while (resList[res+2].nturns[n].result != NONE):
                    res += 1
                    if (resList[res].structure == NONE):
                        resList[res].structure = TURN
    return(resList)

def foundHelices(resList,n):
    """Found all n-helices (3, 4 or 5 type) and assign the structural element (G, H or I)"""
    for res in range(0,len(resList)):
        if (resList[res].structure != NONE): continue # no overwrite
        if (isHelix(resList,res,n)):
            # start (>>) and end (<<) pattern positioning :
            for i in range(0,2):
                resList[res+i].nturns[n].start = START
                resList[res+n+i].nturns[n].end =  END
            # middle patterns positioning [3-helix (3), 4-helix (44) or 5-helix (555)] :
            for i in range(2,n):
                resList[res+i].nturns[n].middle =  MIDDLE[n]
        setNturnPatternResult(resList,res,n)
    setHelixStruct(resList,n)
    return(resList)

def foundNturns(resList):
    """Found all n-turns (3, 4 and 5 types) and assign the structural element (T)"""
    for res in range(len(resList)):
        for n in [3,4,5]:
            try:
                if (isHbond(resList[res],resList[res+n]) == True):
                    # start (>) and end (<) pattern positioning :
                    resList[res].nturns[n].start = START
                    resList[res+n].nturns[n].end =  END
                    # middle pattern positioning [3-helix (33), 4-helix (444) or 5-helix (5555)] :
                    for i in range(1,n):
                        resList[res+i].nturns[n].middle =  MIDDLE[n]
            except:
                pass
            setNturnPatternResult(resList,res,n)
    setNturnsStruct(resList)
    return(resList)

###################### B R I D G E S  &  S T R A N D S #########################

def isParallelBridge(resList,i,j):
    """Test if residues i and j form a parallel bridge. Return a boolean."""
    try:
        if ((isHbond(resList[i-1],resList[j]) == True and isHbond(resList[j],resList[i+1]) == True)\
        or (isHbond(resList[j-1],resList[i]) == True and isHbond(resList[i],resList[j+1]) == True)):
            return(True)
        return(False)
    except:
        return(False)

def isAntiparallelBridge(resList,i,j):
    """Test if residues i and j form an antiparallel bridge. Return a boolean."""
    try:
        if ((isHbond(resList[i],resList[j]) == True and isHbond(resList[j],resList[i]) == True)\
        or (isHbond(resList[i-1],resList[j+1]) == True and isHbond(resList[j-1],resList[i+1]) == True)):
            return(True)
        return(False)
    except:
        return(False)

def setBridgesStruct(resList):
    """Assignment of isolated beta-bridges and then beta-strands to all concerned residues."""
    for res in range(1,len(resList)-1):
        if (resList[res].structure != NONE): continue  # no overwrite
        if ((resList[res].bp1+resList[res].bp2) != 0):
            # Isolated Beta-bridge :
            if (resList[res-1].bp1+resList[res-1].bp2 == 0 and resList[res+1].bp1+resList[res+1].bp2 == 0):
                resList[res].structure = BRIDGE
                resList[resList[res].bp1-1].structure = BRIDGE
            # Extanded strand :
            elif (resList[res-1].bp1+resList[res-1].bp2 != 0 or resList[res+1].bp1+resList[res+1].bp2 != 0):
                resList[res].structure = STRAND
        else:
            # Surrounded by strands = Considered as a strand :
            if (resList[res-1].bp1+resList[res-1].bp2 != 0 and resList[res+1].bp1+resList[res+1].bp2 != 0\
                and resList[res-2].bp1+resList[res-2].bp2 != 0 and resList[res+2].bp1+resList[res+2].bp2 != 0):
                resList[res].structure = STRAND
    return(resList)

def foundBridges(resList):
    """Found all parallel and antiparallel bridges and assign the structural element
    for isolated beta-bridges (B) and extended strands (E)."""
    n = -1
    newStrand = True
    for res_i in range(len(resList)):
        if (resList[res_i].structure != NONE): continue  # no overwrite
        for res_j in range(res_i+2,len(resList)):
            if (resList[res_j].structure != NONE): continue  # no overwrite
            if (isParallelBridge(resList,res_i,res_j) == True or isAntiparallelBridge(resList,res_i,res_j) == True):
                if (isParallelBridge(resList,res_i,res_j) == True):
                    alphabet = ABC_LOWER
                    if (resList[res_i-1].bp1 == 0 or resList[res_i-1].bp1+1 != res_j+1):
                        newStrand = True
                    else:
                        newStrand = False
                elif (isAntiparallelBridge(resList,res_i,res_j) == True):
                    alphabet = ABC_UPPER
                    if (resList[res_i-1].bp1 == 0 or resList[res_i-1].bp1-1 != res_j+1):
                        newStrand = True
                    else:
                        newStrand = False
                if (resList[res_i+1].bp1 == 0 and resList[res_i-1].bp2 == 0):
                    if (newStrand == True):
                        n += 1
                        if (n == 26): n = 0
                    resList[res_i].bp1 = res_j+1
                    resList[res_i].bridge_1 = alphabet[n]
                else:
                    resList[res_i].bp2 = res_j+1
                    if (resList[res_i-1].bp2 == 0 and resList[res_i-2].bp2 == 0): # New strand
                        n += 1
                        if (n == 26): n = 0
                    resList[res_i].bridge_2 = alphabet[n]
                resList[res_j].bp1 = res_i+1
                resList[res_j].bridge_1 = alphabet[n]
    setBridgesStruct(resList)
    return(resList)

################################## B E N D S ###################################

def setBendStruct(resList):
    """Assignment of bends to all concerned residues"""
    for res in range(len(resList)):
        if (resList[res].structure != NONE): continue # no overwrite
        if (resList[res].bend != NONE):
            resList[res].structure = BEND
    return(resList)

######## S E C O N D A R Y   S T R U C T U R E S    A S S I G N M E N T ########


def setSSE(resList):
    """Assignment of secondary structure elements in this order :
        1. 4-helices (H)
        2. isolated beta-bridge (B)
        3. beta-strands (E)
        4. 3-helices (G)
        5. 5-helices (I)
        6. n-turns (T)
        7.bends (S)"""
    foundHelices(resList,4) # alpha-helices
    foundBridges(resList)
    foundHelices(resList,3) # 3_10-helices
    foundHelices(resList,5) # pi-helices
    foundNturns(resList)
    setBendStruct(resList)
    return(resList)
