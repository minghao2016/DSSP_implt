#!/usr/bin/python3

from optparse import OptionParser
import datetime
import re

def readPDB():
	with open(opt.input,'r') as f:
		pdb = f.read()
	return(pdb)

def makeHeader():
	title = "==== Secondary Structure Definition using DSSP method ====\nDATE\t\t{}\n".format(datetime.date.today())
	reference = "REFERENCE\tW. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637\n"

	# regex : 
	header_regex = re.compile("HEADER[ |2-9]{1,}([^\n]*)")
	compnd_regex = re.compile("COMPND[ |2-9]{1,}([^\n]*)")
	source_regex = re.compile("SOURCE[ |2-9]{1,}([^\n]*)")
	author_regex = re.compile("AUTHOR {1,}([^\n]*)")

	header = "HEADER\t\t{}\n".format(' '.join(map(str.rstrip,header_regex.findall(pdb))))
	compnd = "COMPND\t\t{}\n".format(' '.join(map(str.rstrip,compnd_regex.findall(pdb)))) # all COMPND
	source = "SOURCE\t\t{}\n".format(' '.join(map(str.rstrip,source_regex.findall(pdb)))) # all SOURCE
	author = "AUTHOR\t\t{}\n".format(author_regex.findall(pdb)[0].rstrip())
	
	print(title,reference,header,compnd,source,author,sep='')

if __name__ == "__main__":
	parser = OptionParser(usage="usage: %prog [options]",version="%prog 1.0")
	parser.add_option("-i", "--input",dest='input',
	                  help="The file name of a PDB formatted file containing the protein structure data.", metavar="FILE")
	parser.add_option("-o", "--output",dest='output',
	                  help="The  file  name  of  a  DSSP  file to create.", metavar="FILE")
	(opt, args) = parser.parse_args()

	if not opt.input:
	    parser.error('Input pdb file not given.')

	pdb = readPDB()

	makeHeader()
