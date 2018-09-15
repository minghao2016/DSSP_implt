"""@package additional_functions
Documentation for this module.

More details.
"""
import os
from optparse import OptionParser
from subprocess import call
import datetime as dt

def argsParsing():
    """Short summary.

    Returns
    -------
    type
        Description of returned object.

    """
    """
    Automatic creation of the help and parsing of arguments
    Options added : -i (input file), -o (output file)
    Return the input and output files
    """
    parser = OptionParser(usage="usage: %prog [options]", version="%prog 1.0")
    parser.add_option("-i", "--input", dest='input', metavar="FILE",
        help="the file name of a PDB formatted file containing the protein structure data")
    parser.add_option("-o", "--output", dest='output', metavar="FILE",
        help="the  file  name  of  a  DSSP  file to create")
    (opt, args) = parser.parse_args()
    if not opt.input: parser.error('Input pdb file not given. Use the [-i FILE] option.')
    if not os.path.isfile(opt.input): parser.error('Input file does not exist.')
    return(opt)

def hydrAddition(opt):
    """Short summary.

    Parameters
    ----------
    opt : type
        Description of parameter `opt`.

    Returns
    -------
    type
        Description of returned object.

    """
    """
    Addition of hydrogens in the pdb file using a
    shell command line to run the reduce program
    """
    call(["./bin/reduce -NOFLIP "+opt.input+" 1>"+opt.input+".H"+" 2>"+opt.input+".H.log"],shell=True)

def lineHeader(dic):
    """
    Access to all items of a dictionary and return a string
    of all existing items
    """
    l = ""
    for mol_id,mol_items in dic.items():
        l += "MOL_ID: " + mol_id + "; "
        for key, item in mol_items.items():
            if (item != ""):
                l += key.upper() + ": " + item.upper() + "; "
    return(l)

def makeHeader(pdb):
    """
    Make and return the header of the dssp output using data from the pdb file
    """
    header = "==== Secondary Structure Assignment using DSSP method ====\nDATE\t\t{}\n".format(dt.date.today())
    header += "REFERENCE\tW. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n"
    header += "HEADER\t\t{}{:>28}\n".format(pdb.header["head"].upper(),pdb.header["deposition_date"])
    header += "COMPND\t\t{}\n".format(lineHeader(pdb.header["compound"])) # all COMPND
    header += "SOURCE\t\t{}\n".format(lineHeader(pdb.header["source"])) # all SOURCE
    header += "AUTHOR\t\t{}".format(pdb.header["author"].upper())
    return(header+"\n")

def displayResults(opt,pdb,dssp):
    #    print("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN")
    header = makeHeader(pdb)
    descp = "  #  RESIDUE AA STRUCTURE BP1 BP2    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
    if (opt.output):
        with open(opt.output,'w') as filout:
            filout.write(header+descp+'\n')
            for i in range(len(dssp)):
                l = dssp[i]
                filout.write("{:>5d}{:>5d}{:>2s}{:>2s}{:>3s}{:>2s}{:>1s}{:>1s}{:>1s}{:>1s}{:>1s}{:>1s}{:>4d}{:>4d}{:>1s}{:>7.3f}{:>6.1f}{:>6.1f}{:>6.1f}{:>6.1f}{:>7.1f}{:>7.1f}{:>7.1f}\n"\
                .format(r.index,r.resNum,r.chainID,r.resName,r.structure,r.turn3['rlt'],r.turn4['rlt'],r.turn5['rlt'],r.bend,r.chirality,r.bridge_1,r.bridge_2,r.bp1,r.bp2,r.sheet,\
                    r.tco,r.kappa,r.alpha,r.phi,r.psi,r.CA[0],r.CA[1],r.CA[2]))

    else:
        print(header,descp,sep="")
        for i in range(len(dssp)):
            r = dssp[i]
            print("{:>5d}{:>5d}{:>2s}{:>2s}{:>3s}{:>2s}{:>1s}{:>1s}{:>1s}{:>1s}{:>1s}{:>1s}{:>4d}{:>4d}{:>1s}{:>7.3f}{:>6.1f}{:>6.1f}{:>6.1f}{:>6.1f}{:>7.1f}{:>7.1f}{:>7.1f}"\
                .format(r.index,r.resNum,r.chainID,r.resName,r.structure,r.nturns[3].result,r.nturns[4].result,r.nturns[5].result,r.bend,r.chirality,r.bridge_1,r.bridge_2,r.bp1,r.bp2,r.sheet,\
                    r.tco,r.kappa,r.alpha,r.phi,r.psi,r.CA[0],r.CA[1],r.CA[2]))
