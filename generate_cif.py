from os import listdir
from MOFLIB import *
import pormake as pm
from pickle import load
import numpy as np

db = pm.Database()
pcu = db.get_topo('pcu')
pcu.describe()
sbu = pm.BuildingBlock('sbu.xyz')
#sbu.view()

dir_ = listdir('./molecules/')
for mol_file in dir_:
    try:
        with open('./molecules/'+mol_file, 'rb') as f: atoms, bonds, pos = load(f)
        mol = molecule(atoms=atoms, bonds=bonds, pos=pos)

        M = mol_matrix(mol)
        pos = np.array(pos)
        x =-1
        for i in range(len(M)-3):
            if M[len(M)-3][1][i]>0: x=i
        y =-1
        for i in range(len(M)//2-3):
            if M[len(M)//2-3][1][i]>0: y=i

        alpha=0.9
        pos[len(M)-3]+=alpha*(pos[x]-pos[len(M)-3])
        pos[len(M)//2-3]+=alpha*(pos[y]-pos[len(M)//2-3])

        M = [[(M[i][0] if i!=len(M)-3 and i!=len(M)//2-3 else 0), [M[i][1][j] for j in range(len(M))]] for i in range(len(M))]
        delete = [len(M)-1, len(M)-2, len(M)//2-1, len(M)//2-2]
        M = [[(M[i][0] if i!=len(M)-3 and i!=len(M)//2-3 else 0), [M[i][1][j] for j in range(len(M)) if j not in delete]] for i in range(len(M)) if i not in delete]
        pos = [pos[i] for i in range(len(pos)) if i not in delete]

        mol = matrix_mol(M, pos)
        write_xyz(mol, 'linker.xyz')
        linker = pm.BuildingBlock('linker.xyz')
        bbs = [sbu]
        for i in range(3): bbs.append(linker.copy())
        builder = pm.Builder()
        MOF = builder.build(topology=pcu, bbs=bbs)
        MOF.write_cif('cifs/%s.cif'%(mol_file))
    except EOFError: pass
