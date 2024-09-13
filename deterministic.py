import sys
sys.path.insert(1, '../')
from MOFLIB import *
import numpy as np
from pickle import load

#    The issue https://github.com/openbabel/openbabel/issues/1934 in github
#    shows that the openbabel make3D is not deterministic which makes the current
#    algorithm for generating molecules undeterministic. The following code shows
#    proof of that

mol_file = '5020_mbPXoIZAGm'
pos1, pos2 = list(), list()
with open('../molecules/'+mol_file, 'rb') as f:
    atoms, bonds = load(f)
    mol = molecule(atoms=atoms, bonds=bonds)
    mol.make3D(steps=2)
    mol.localopt()
    pos1 = mol_pos(mol)

with open('../molecules/'+mol_file, 'rb') as f:
    atoms, bonds = load(f)
    mol = molecule(atoms=atoms, bonds=bonds)
    mol.make3D(steps=2)
    mol.localopt()
    pos2 = mol_pos(mol)

pos1 = np.array(pos1)
pos2 = np.array(pos2)
diff = pos1-pos2
a=np.sqrt(sum([diff[i][0]*diff[i][0] for i in range(len(diff))]))/len(diff)
b=np.sqrt(sum([diff[i][1]*diff[i][1] for i in range(len(diff))]))/len(diff)
c=np.sqrt(sum([diff[i][2]*diff[i][2] for i in range(len(diff))]))/len(diff)
print('RMSE in x axis: %.6f'%(a))
print('RMSE in y axis: %.6f'%(b))
print('RMSE in z axis: %.6f'%(c))
print("RMSE: %.6f"%(np.sqrt(a*a+b*b+c*c)))
