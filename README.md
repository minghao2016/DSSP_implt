![Python version](https://img.shields.io/badge/python-3.6-brightgreen.svg)

# DSSP implementation
Secondary Structure Assignment using the DSSP method.

## Installation
### Requirements
- Python 3.6
- Biopython package (Bio.SeqUtils and Bio.PDBmath)

### Clone the repository
```shell
git clone https://github.com/kabhel/DSSP_implmt.git
cd DSSP_implmt
chmod +x dssp.py
```
## Run the program
### Examples
- Protein with one domain : 
```shell
./dssp.py -i data/1bta.pdb
```
- Protein with more than one domain :
```shell
./dssp.py -i data/1pnk.pdb
```

### Get help
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

## Output parameters
- **RESIDUE :** two columns of residue numbers. First column is DSSP's sequential residue number, starting at the first residue actually in the data set and including chain breaks; this number is used to refer to residues throughout. Second column gives crystallographers' 'residue sequence number','insertion code' and 'chain identifier' (see protein data bank file record format manual), given for reference only.
- **AA :** One letter amino acid code, lower case for SS-bridge CYS.
- **S (first column in STRUCTURE block) :** Compromise summary of secondary structure, intended to approximate crystallographers' intuition, based on columns 19-38, which are the principal result of DSSP analysis of the atomic coordinates.
- **BP1 BP2 :** Residue number of first and second bridge partner followed by one letter sheet label
- **ACC :** Number of water molecules in contact with this residue *10. or residue water exposed surface in Angstrom**2.
- **N-H-->O etc. :** Hydrogen bonds; e.g. -3,-1.4 means: if this residue is residue i then N-H of i is h-bonded to C=O of i-3 with an electrostatic H-bond energy of -1.4 kcal/mol. There are two columns for each type of H-bond, to allow for bifurcated H-bonds.
- **TCO :** Cosine of angle between C=O of residue i and C=O of residue i-1. For alpha-helices, TCO is near +1, for beta-sheets TCO is near -1. Not used for structure definition.
- **KAPPA :** Virtual bond angle (bend angle) defined by the three C-alpha atoms of residues i-2, i, i+2. Used to define bend (structure code 'S').
- **ALPHA :** Virtual torsion angle (dihedral angle) defined by the four C-alpha atoms of residues i-1, i, i+1, i+2. Used to define chirality (structure code '+' or '-').
- **PHI PSI :** IUPAC peptide backbone torsion angles
- **X-CA Y-CA Z-CA :** C-alpha atom coordinates

## Reference
W. Kabsch and C.Sander, Biopolymers 22 (1983) 2577-2637 
