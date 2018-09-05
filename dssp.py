#!/usr/bin/python3

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options]",version="%prog 1.0")
parser.add_option("-i", "--input",dest='input',
                  help="The file name of a PDB formatted file containing the protein structure data.", metavar="FILE")
parser.add_option("-o", "--output",dest='output',
                  help="The  file  name  of  a  DSSP  file to create.", metavar="FILE")


(options, args) = parser.parse_args()

if not options.input:
    parser.error('Input pdb file not given.')

if not options.output:
    options.output = options.input+".dssp"


