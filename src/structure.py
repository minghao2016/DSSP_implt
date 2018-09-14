"""@package structure
Documentation for this module.

More details.
"""

# helix structural patterns
HELIX_STRUC = {'3':'G','4':'H','5':'I'}
HELICES = [5,4,3]

# n-turn patterns
NONE = ' '
START = '>'
END = '<'
MIDDLE = {'3':'3','4':'4','5':'5'}
START_END = 'X'

def isHelix(dssp,res,n):
    """
    n-helix(i,i+n-1) = [n-turn(i-1) and n-turn(i)]
    """
    try:
        if (dssp[res][str(n)+'-hbonds'] == True and dssp[res+1][str(n)+'-hbonds'] == True):
            return(True)
        return(False)
    except:
        return(False)

def nTurnPatterns(dssp):
    """
    Helix patterns :
    3-helix : >>3<<
    4-helix : >>44<<
    5-helix : >>555<<
    """

    for res in range(len(dssp)):
        # n-Helix, with n = [3,4,5]
        for n in HELICES:
            if (isHelix(dssp,res,n)):
                for i in range(0,2):
                    # helix start : res(i):'>' & res(i+1):'>'
                    dssp[res+i][str(n)+'-turns']['start'] = START
                    # helix end : res(i+n):'<' & res(i+n+1):'<'
                    dssp[res+n+i][str(n)+'-turns']['end'] = END
                for i in range(2,n):
                    # helix middle : res(i+2):'[3,4,5]' & res(i+3):'[4,5]' & res(i+4):'[5]'
                    dssp[res+i][str(n)+'-turns']['middle'] = MIDDLE[str(n)]

            if (dssp[res][str(n)+'-turns']['end'] == END):
                if (dssp[res][str(n)+'-turns']['start'] == START):
                    # '>' + '[3,4,5]' + '<' = 'X' (start & end)
                    dssp[res][str(n)+'-turns']['rlt'] = START_END
                else:
                    # ' ' + '<' = '<' (end)
                    dssp[res][str(n)+'-turns']['rlt'] = END
            else:
                if (dssp[res][str(n)+'-turns']['start'] == START):
                    # '>' + ' ' = '>' (start)
                    dssp[res][str(n)+'-turns']['rlt'] = START
                else:
                    if (dssp[res][str(n)+'-turns']['middle'] != NONE):
                        # (' ' + ' ') and (middle != ' ') =  [3,4,5] (middle)
                        dssp[res][str(n)+'-turns']['rlt'] = MIDDLE[str(n)]
    return(dssp)

def setStructure(dssp):
    for res in range(len(dssp)):
        # Helices
        for n in HELICES:
            if (dssp[res][str(n)+'-turns']['rlt'] == START):
                while (dssp[res+2][str(n)+'-turns']['rlt'] != NONE):
                    res += 1
                    if (dssp[res]['structure'] == NONE):
                        dssp[res]['structure'] = HELIX_STRUC[str(n)]

        # Strands
        if (dssp[res]['bridge']['i']+dssp[res]['bridge']['j'] != ''):
             dssp[res]['structure'] = 'E'
            
    return(dssp)
