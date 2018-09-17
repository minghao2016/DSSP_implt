# DSSP implementation
Secondary Structure Assignment using the DSSP method.

## Installation
### Requirements
- Python 3.6
- Biopython package (for the Bio.SeqUtils and Bio.PDBmath modules)

### Clone the repository
```shell
git clone https://github.com/kabhel/DSSP_implmt.git
cd DSSP_implmt
```
## Run the program
### Examples
```shell
./dssp.py -i data/1bta.pdb
./dssp.py -i data/1itv.pdb -o res/1itv.dssp
```

### Get help
```
$ ./dssp.py -h
Usage: dssp.py [options]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i FILE, --input=FILE
                        the file name of a PDB formatted file containing the
                        protein structure data
  -o FILE, --output=FILE
                        the  file  name  of  a  DSSP  file to create
```

## Reference
W. Kabsch and C.Sander, Biopolymers 22 (1983) 2577-2637 
