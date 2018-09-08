# DSSP : Dictionary of protein secondary structure

## Requirements
- Python3

## Installation
```
git clone https://github.com/kabhel/DSSP_implmt.git
cd DSSP_implmt
chmod +x dssp.py
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
