from MOFLIB import *
from pickle import load
import pormake as pm
from random import random
mol = read_mol('N622.mol')
atoms, bonds=mol_info(mol)
atoms = [a if a!=0 else 9 for a in atoms]
pos = mol_pos(mol)
O = [27, 28,29, 30]
N=len(atoms)
pos+=[(pos[O[i]][0]+random()*1.0,pos[O[i]][1]+random()*1.0, pos[O[i]][2]+random()*1.0 ) for i in range(4)]
atoms+=[1]*4
bonds+=[(1, O[i]+1, N+i+1) for i in range(4)]
mol = molecule(atoms=atoms, bonds=bonds, pos=pos)
mol = optimizedMol(mol)

atoms, bonds=mol_info(mol)
atoms = [a if a!=9 else 0 for a in atoms]
pos = mol_pos(mol)

mol = molecule(atoms=atoms, bonds=bonds, pos=pos)
write_xyz(mol, 'sbu.xyz')
new = pm.BuildingBlock('sbu.xyz')
new.view()
