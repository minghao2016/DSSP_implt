"""@package additional_functions

Functions that improve the use of the program :
Creation of a help and Displaying results
"""
import os
from optparse import OptionParser
from subprocess import call
import datetime as dt

def argsParsing():
    """Automatic creation of the help and parsing of arguments. Return the input and output files."""
    # Help creation :
    parser = OptionParser(usage="usage: %prog [options]", version="%prog 1.0")
    parser.add_option("-i", "--input", dest='input', metavar="FILE",
        help="the file name of a PDB formatted file containing the protein structure data")
    parser.add_option("-o", "--output", dest='output', metavar="FILE",
        help="the  file  name  of  a  DSSP  file to create")

    (opt, args) = parser.parse_args() # Input and Output files are  attribute of opt
    if not opt.input: parser.error('Input pdb file not given. Use the [-i FILE] option.')
    if not os.path.isfile(opt.input): parser.error('Input file does not exist.')
    return(opt)

def hydrAddition(pdbfile):
    """Addition of hydrogens in the pdb file using a shell command line to run the reduce program"""
    call(["./bin/reduce -NOFLIP " + pdbfile + " 1>" + pdbfile + ".H" + " 2>" + pdbfile + ".H.log"],shell=True)

def lineHeader(dic):
    """Access to all items of a dictionary. Return a string of all existing items"""
    l = ''
    for mol_id,mol_items in dic.items():
        l += 'MOL_ID: ' + mol_id + '; '
        for key, item in mol_items.items():
            if (item != ''):
                l += key.upper() + ': ' + item.upper() + '; '
    return(l)

def makeHeader(pdb):
    """Make and return the header of the resList output using data from the pdb file"""
    header = "==== Secondary Structure Assignment using DSSP method ====\nDATE\t\t{}\n".format(dt.date.today())
    header += "REFERENCE\tW. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n"
    header += "HEADER\t\t{}{:>28}\n".format(pdb.header["head"].upper(),pdb.header["deposition_date"])
    header += "COMPND\t\t{}\n".format(lineHeader(pdb.header["compound"])) # all COMPND lines
    header += "SOURCE\t\t{}\n".format(lineHeader(pdb.header["source"])) # all SOURCE lines
    header += "AUTHOR\t\t{}".format(pdb.header["author"].upper())
    return(header+"\n")

def displayResults(opt,pdb,resList):
    """ """
    header = makeHeader(pdb)
    descp = "  #  RESIDUE AA STRUCTURE BP1 BP2    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"

    # Results save in an output file :
    if (opt.output):
        f = open(opt.output,'w')
        f.write(header+descp+'\n')
    # Or they are displayed on the terminal :
    else:
        print(header,descp,sep='')

    for i in range(len(resList)): 
        r = resList[i] # A Residue (instance of the Residue class)
        line = "{:>5d}{:>5d}{:>2s}{:>2s}{:>3s}{:>2s}{:>1s}{:>1s}{:>1s}{:>1s}{:>1s}{:>1s}{:>4d}{:>4d}{:>1s}{:>7.3f}{:>6.1f}{:>6.1f}{:>6.1f}{:>6.1f}{:>7.1f}{:>7.1f}{:>7.1f}\n"\
                .format(r.index,r.resNum,r.chainID,r.resName,r.structure,r.nturns[3].result,r.nturns[4].result,r.nturns[5].result,
                r.bend,r.chirality,r.bridge_1,r.bridge_2,r.bp1,r.bp2,r.sheet,r.tco,r.kappa,r.alpha,r.phi,r.psi,r.CA[0],r.CA[1],r.CA[2])
        if (opt.output):
            f.write(line)
        else:
            print(line,end='')
