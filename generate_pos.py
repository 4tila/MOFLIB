from os import listdir
from MOFLIB import *
from pickle import load, dump

dir_ = listdir('./molecules/')
for mol_file in dir_:
    try:
        with open('./molecules/'+mol_file, 'rb') as f: atoms, bonds, pos = load(f)

        mol = molecule(atoms=atoms, bonds=bonds)
        mol.make3D(steps=2)
        #ol.localopt()
        
        pos = mol_pos(mol)
        with open('./molecules/'+mol_file, 'wb+') as f: dump([atoms, bonds, pos], f)
    except EOFError: pass

