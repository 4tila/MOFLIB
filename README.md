# MOFLIB

Algorithm for generating metal organic framework semiconductors with a metropolis algorithm, where the the mutation is done directly at the molecular structure. The semiconductors generated until now have all been narrow gap semiconductors

## MOFLIB

MOFLIB  is a library developed  that contain some special functions used by this project. For example it converts a smile to two lists containing atomic numbers and bonds and functions to convert those atoms and bonds back to the original molecule.

```
>>> # smiles for benzene
>>> mol = molecule(smile="c1ccccc1")
>>> atoms, bonds = mol_info(mol)
>>> atoms
[6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1]
>>> bonds
[(1, 1, 2), (2, 2, 3), (1, 3, 4), (2, 4, 5), (1, 5, 6), (2, 1, 6), (1, 1, 7), (1, 2, 8), (1, 3, 9), (1, 4, 10), (1, 5, 11), (1, 6, 12)]
>>> mol = molecule(atoms=atoms, bonds=bonds)
>>> print(mol)
c1ccccc1

```

It contains also functions to generate arbitrary mutations of molecules

```
>>> from MOFLIB import *
>>> mol = molecule(smile="c1ccccc1")
>>> mol2 = mutation(mol)
>>> print(mol2)
c1ccncc1	

>>> mol2 = mutation(mol)
>>> print(mol2)
C1=CC=CC1	

>>> mol2 = mutation(mol)
>>> print(mol2)
c1(ccccc1)F	

```
Here are some examples of possible mutations of benzene done by the algorithm together with the unmutated molecule for comparison (first image)

<img src="https://github.com/4tila/MOFLIB/blob/main/imgs/transparent_img01.png" width="250" height="250" /><img src="https://github.com/4tila/MOFLIB/blob/main/imgs/transparent_img02.png" width="250" height="250" /><img src="https://github.com/4tila/MOFLIB/blob/main/imgs/transparent_img03.png" width="250" height="250" />

<img src="https://github.com/4tila/MOFLIB/blob/main/imgs/transparent_img04.png" width="250" height="250" /><img src="https://github.com/4tila/MOFLIB/blob/main/imgs/transparent_img05.png" width="250" height="250" /><img src="https://github.com/4tila/MOFLIB/blob/main/imgs/transparent_img06.png" width="250" height="250" />


A function to generate a adjacency matrix with information on the atomic numbers

``` 
>>> mol = molecule(smile="C")
>>> M = mol_matrix(mol)
>>> for i in range(len(M)): print("atomic number of index %d is %d"%(i, M[i][0]))
...
atomic number of index 0 is 6
atomic number of index 1 is 1
atomic number of index 2 is 1
atomic number of index 3 is 1
atomic number of index 4 is 1
>>> #adjcency matrix
>>> for i in range(len(M)): print(*M[i][1])
... 
0 1 1 1 1
1 0 0 0 0
1 0 0 0 0
1 0 0 0 0
1 0 0 0 0

```

A function that generate a dycarboxilate linker using a molecule as a seed

```
>>> mol = molecule(smile="c1ccccc1")
>>> M = dicarboxylate_linker(mol)
>>> mol2 = matrix_mol(M)
>>> print(mol2)
c1ccc(c(c1)c1c(cccc1)C(=O)[O])C(=O)[O]	
```
Here's the result of applying the function to benzene

<img src="https://github.com/4tila/MOFLIB/blob/main/imgs/transparent_img07.png" width="250" height="250" />


To find out more about those functions use

```
import MOFLIB
help(MOFLIB)
```

## Encoder

The "encoder.py" file is going to take nodes, edges and topology of MOF's in the QMOF database and base on openbabel's utility for molecular fingerprint generate a encoding of the MOF that is going to be interpreted by the machine learning model (called X).

It also going to produce a column Y that is 1 if the the bandgap computed by PBE is between 1.8 eV and 2.8eV and 0 if it's outside the interval.

## Classifier

The classifier is a simple XGBoost with objective set to "binary logistic" that is going to classify the MOF's (saved in X) to 0 or 1 depending if it's inside or outside the interval of energy mentioned. The accuracy is obtained using 1/10 of the training data as validation and saved on "model01.acc" file. The classifier model is saved in "model01.mod"

## Metropolis Algorithm 

The algorithm starts with a seed and causes mutations at each iteration to produce new candidates. If the mutation increases the objective function, then the new candidate replaces the old one, otherwise it has a chance to replace with certain probability. This makes the algorithm avoids local optima.

The metropolis algorithm in this work is going to be used to maximize the probability the MOF has PBE bandgap between 1.8eV and 2.8 eV. That is given a topology and a SBU, the algorithm will search for new linkers that maximize this probability,

In the case of this work the objective function is going be the classifier model previously mentioned, the seed is a smile given by the GDB-11 database and the mutation is done by the function "mutation" that directly changes the structure of the molecule. I talked about that function in the section "MOFLIB" that describes the library. You can also find out more using the python help function.

To generate the linker will use the function "dicarboxylate_linker" from the library MOFLIB that was developed.

Candidates generated with more than 50% of probability are saved to "./molecules" directory. The first 4 characters of the file are the probability  and the rest is a random string to uniquely identify each linker

## Generate Cif

The "generate_cif.py" makes a cif from all the linkers in the molecules directory. The "generate_cif_file.py" allows you to optimize only one of the molecules where a linker can be chosen by using the pyfzf interface.

```
python3 generate_cif.py
```

## Proof Of Concept

Molecule 5020_mbPXoIZAGm (in the molecules file) with corresponding CIF 5020_mbPXoIZAGm.cif (At cifs file) has a band gap of 0.391 eV by GGA-PBE, with calculations done on materials studio.

<img src="https://github.com/4tila/MOFLIB/blob/main/imgs/5020.jpg" width="500" height="400" />

## Acknowledgement

I'd like to thank professor David Azevedo that came with the idea of the project and guided the research. This work was done at universidade de Brasília and funded by CNPQ and it's free to use. My email is loboatila@gmail.com in case you would like to contact me
