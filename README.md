![Python version](https://img.shields.io/badge/python-3.6-brightgreen.svg)

# Secondary Structure Assignment using DSSP method

## Requirements
- Python 3.6
- Packages : math, re, optparse, Bio.SeqUtils, Bio.PDBmath, datetime

## Installation
```
git clone https://github.com/kabhel/DSSP_implmt.git
cd DSSP_implmt
chmod +x dssp.py
```

## Examples
```
./dssp.py -i data/1bta.pdb
./dssp.py -i data/1pnk.pdb
```

## Help
```
$ ./dssp.py -h
Usage: dssp.py [options]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i FILE, --input=FILE
                        The file name of a PDB formatted file containing the
                        protein structure data.
  -o FILE, --output=FILE
                        The  file  name  of  a  DSSP  file to create.
```

## Reference
W. Kabsch and C.Sander, Biopolymers 22 (1983) 2577-2637 
