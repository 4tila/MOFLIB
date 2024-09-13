import os
import numpy as np
from openbabel import openbabel, pybel
from random import random, randint as ri, shuffle

# TODO cite openbabel in the article
# TODO use benzene as linker in the demonstration of the article
# TODO cite qmof and materials project
# TODO cite pormake

CHARS = 'q w e r t y u i o p a s d f g h j k l z x c v b n m Q W E R T Y U I O P A S D F G H J K L Z X C V B N M'.split(' ')
subs = {1:[9], 9:[1], 8:[16], 16:[8]}

def write_xyz(mol, fname='test.xyz'):
    """
    def write_xyz(mol, fname='test.xyz')

    writes molecule to .xyz with bonds to "fname". Atom set to 0 is saved as element X.
    Element X can be interpreted by PORMAKE as connection point.

    """
    atoms, bonds= mol_info(mol)
    mol.write('xyz', fname, overwrite=True)
    with open(fname,'r') as f: s=f.read()
    for b in bonds: s+='{}{}{}{}\n'.format(''.ljust(4), str(b[1]-1).ljust(4), str(b[2]-1).ljust(4), bond_type(b[0]))
    with open(fname, 'w+') as f: f.write(s.replace('*', 'X'))
    return ;
def bond_type(b):
    """
    def bond_type(b):

    returns a bond type (simple, double, triple) given number of bonds.
    This is used by write_xyz function.
    """
    if b==1: return 'S'
    if b==2: return 'D'
    if b==3: return 'T'
def energy(mol):
    """
    def energy(mol):

    returns energy of molecule given openbabel molecule.
    
    The function saves molecule to auxiliar file "test.mol" and uses the obenergy utility to compute energy.
    """
    mol.write('mol', 'test.mol', overwrite=True)
    os.popen('obenergy -ff UFF test.mol > aux').read()
    with open('aux', 'r') as f: E=float(f.read().split('\n')[-2].replace('kJ/mol', '').split('=')[1])
    return E
def read_mol(fname='test.mol'):
    """
    def read_mol(fname='test.mol'):

    return openbabel molecule given a molecule's file. General file types are .mol and .pdb
    default filename is "test.mol"
    """
    tt = fname.split('.')[-1]
    return [m for m in pybel.readfile(tt, fname)][0]
def optimizedMol(mol=None, smile=None, verbose=False, iters=1000):
    """
    def optimizedMol(mol=None, smile=None, verbose=False, iters=1000)
    
    given openbabel molecule or smiles of molecule, returns openbabel optimized molecule using UFF.

    It uses the obenergy utility to perform optimization. The maximum number of iterations is given by
    iters variable. default is 1000

    If smile and molecule is passed, then it will use the smile.

    If verbose is set to True then it will print obminimize messages.
    """
    if smile!=None:
        mol = molecule(smile=smile)
        mol.make3D()
    fname = random_string()+'.mol'
    mol.write('mol', fname, overwrite=True)
    msg = os.popen('obminimize -h -ff GAFF -n %d %s > test.pdb'%(iters, fname)).read()
    if verbose: print(msg)
    os.popen('rm %s'%(fname))
    mol = read_mol('test.pdb')
    return mol
def optimizedMolInfo(mol=None, smile=None, verbose=False):
    """
    def optimizedMolInfo(mol=None, smile=None, verbose=False):

    given openbabel molecule or smile, returns tuple with atoms, pos and bonds in this order.

    If smile and molecule is passed, then it will use the smile.
    """
    mol = optimizedMol(smile=smile, verbose=verbose) if smile!=None else optimizedMol(mol=mol, verbose=verbose)
    atoms, bonds = mol_info(mol)
    pos = mol_pos(mol)
    return atoms, pos, bonds
def random_string(n=10):
    """
    def random_string(n=10):

    given a number "n", returns a string of that size
    """

    return ''.join([CHARS[ri(0, len(CHARS)-1)] for i in range(n)])
def rnd(a, b):
    """
    def rnd(a, b):

    given two numbers "a" and "b", returns a number between a and b
    """
    return a+(b-a)*random()
def convert_to_bit(ls, N):
    """
    def convert_to_bit(ls, N):

    given a list of positions in a 1D array of length N, return a 1D array
    of length N which is set to 1 at the specified positions and 0 otherwise.
    """

    x = [0]*N
    for e in ls:
        if e<N: x[e]=1
        else: print('OVERFLOW --- Attempt to put at position %d in list of size %d'%(e, N))
    return x
def padding(ls, N):
    """
    def padding(ls, N):

    given a list "ls", if the list has length less than N then it returns
    onlythe first N elements, otherwise returns the list padded with 0's
    so that it has length N
    """
    return [ls[i] if i<len(ls) else 0 for i in range(N)]

def rdLs(ls):
    """
    def rdLs(ls):

    given a list "ls", returns randomly one of its elements
    
    if the list has length 0, then it returns None
    """
    return ls[ri(0, len(ls)-1)] if len(ls)>0 else None
def molecule(atoms=None, bonds=None, pos=list(), smile=None):
    """
    def molecule(atoms=None, bonds=None, pos=list(), smile=None):

    given a list of atoms, bonds between those atoms and atomics
    positions or just its smiles, returns a openbabel molecule.

    atoms is a list of integers containing the atomic numbers.

    bonds is a list of 3D tuples. For each tuple, the first element 
    is the number os bonds (eg.: 2 is double bond), and the second and
    third are the atoms that are going to be connected by that bond.

    If list of positions is not passed, then the function generates random
    positions.

    example: 
        
        >>> atoms = [6, 1, 1, 1, 1]
        >>> bonds = [(1, 1, 2), (1, 1, 3), (1, 1, 4), (1, 1, 5)]
        >>> print(molecule(atoms=atoms, bonds=bonds))
        C

        >>> atoms = [6, 7, 1]
        >>> bonds = [(3, 1, 2), (1, 1, 3)]
        >>> print(molecule(atoms=atoms, bonds=bonds))
        C#N
        >>> mol = molecule(smile="C1=CC=CC=C1")
    """
    if smile!=None:
        mol = pybel.readstring('smi', smile)
        mol.addh()
        return mol
    if len(pos)==0: pos = [[random()*5 for j in range(3)] for i in range(len(atoms))]
    mol = openbabel.OBMol()
    for i in range(len(atoms)):
        a = mol.NewAtom()
        a.SetAtomicNum(atoms[i])
        a.SetVector(pos[i][0], pos[i][1], pos[i][2])
    for b in bonds: mol.AddBond(b[1], b[2], b[0])
    return pybel.Molecule(mol)
def mol_pos(mol):
    """
    def mol_pos(mol):

    given openbabel molecule, returns the positions of each of its atoms
    """
    return [a.coords for a in mol]
def mol_info(mol_input):
    """
    def mol_info(mol_input):

    given openbabel molecule or openbabel.openbabel.OBMol, returns tuple of atoms and bonds.
    
    example:
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
    """
    flag = type(mol_input)==pybel.Molecule
    if flag: mol = mol_input.__dict__['OBMol']
    else: mol = mol_input.copy()
    A, B = mol.NumAtoms(), mol.NumBonds()
    bonds = [(mol.GetBond(i).GetBondOrder(), mol.GetBond(i).GetBeginAtomIdx(), mol.GetBond(i).GetEndAtomIdx()) for i in range(B)]
    atoms = [mol.GetAtom(i+1).GetAtomicNum() for i in range(A)]
    return atoms, bonds
def dicarboxylate_linker(mol):
    """
    def dicarboxylate_linker(mol):

    given openbabel molecule, returns a dicarboxylate a adjacency matrix based on it.

    The function makes a list of carbons that has at least one hydrogen and selects two
    distinct carbons, call it C1 and C2 and removes one hydrogen of each of them. It then
    duplicates the molecule. Call it C1' and C2' carbons corresponding to C1 and C2 in the
    duplicate. It joins by a simple bond C1 to C1' and connects C2 and C2' to two CO_2^{-}
    molecules.

    *** ATTENTION ***

    The adjacency matrix is different from the usual adjacency matrix because it also contains
    the atomic numbers embedded in it. To find out about the adjacency matrix format type help(mol_matrix)

    example:
        >>> from MOFLIB import *
        >>> mol = molecule(smile="c1ccccc1")
        >>> linker = dicarboxylate_linker(mol)
        >>> linker = matrix_mol(linker)# converts matrix to openbabel molecule
        >>> print(linker)
        c1ccc(cc1c1cccc(c1)C(=O)[O])C(=O)[O]
    """
    M = mol_matrix(mol)
    cs = [i for i in range(len(M)) if M[i][0]==6 and len([j for j in range(len(M)) if M[j][0]==1 and M[i][1][j]==1])>0]
    if len(cs)<2: return M
    shuffle(cs)
    c1, c2 = cs[0], cs[1]
    h1, h2=-1, -1
    for i in range(len(M)):
        if M[i][0]==1 and M[c1][1][i]==1: h1=i
        if M[i][0]==1 and M[c2][1][i]==1: h2=i
    M[c1][1][h1]=0
    M[c2][1][h2]=0
    for i in range(len(M)): M[i][1]+=[0, 0, 0]
    M.append([6, [0 for j in range(len(M[0][1]))]])
    for i in range(2): M.append([8, [0 for j in range(len(M[0][1]))]])
    M[-3][1][-2]=M[-2][1][-3]=2
    M[-3][1][-1]=M[-1][1][-3]=1
    M[c1][1][-3]=M[-3][1][c1]=1
    N=len(M)
    for i in range(N):
        M.append([M[i][0], [0 for j in range(N)]+M[i][1]])
        M[i][1]=M[i][1]+[0 for j in range(N)]
    M[c2][1][c2+N]=M[c2+N][1][c2]=1
    delete = [h1, h2, h1+N, h2+N]
    M = [[M[i][0], [M[i][1][j] for j in range(len(M[i][1])) if j not in delete]] for i in range(len(M)) if i not in delete]
    return M
def matrix_delete(M, delete):
    """
    def matrix_delete(M, delete):

    given adjacency matrix, deletes rows and columns at delete. Use delete as set instead of list,
    otherwise it will run slower.

    *** ATTENTION ***

    The adjacency matrix is different from the usual adjacency matrix because it also contains
    the atomic numbers embedded in it. To find out about the adjacency matrix format type help(mol_matrix)

    """
    return [[M[i][0], [M[i][1][j] for j in range(len(M[i][1])) if j not in delete]] for i in range(len(M)) if i not in delete]
def mol_matrix(mol):
    """
    def mol_matrix(mol):

    given a openbabel mol, returns a adjacency matrix.

    The adjacency matris is different from the usual adjacency matrix because it also has the atomic
    numbers embedded which are the first element of each of the matrix elements. The second element would
    be the corresponding row of the adjacency matrix.

    example:

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

    """
    atoms, bonds = mol_info(mol)
    bonds = [(b[0], b[1]-1, b[2]-1) for b in bonds]
    matrix = [[atoms[i], [0 for j in range(len(atoms))]] for i in range(len(atoms))]
    for b in bonds:
        matrix[b[1]][1][b[2]]=b[0]
        matrix[b[2]][1][b[1]]=b[0]
    return matrix
def matrix_mol(M, pos=list()):
    """
    def matrix_mol(M, pos=list()):

    given adjcency matrix, returns openbabel molecule

    *** ATTENTION ***

    The adjacency matrix is different from the usual adjacency matrix because it also contains
    the atomic numbers embedded in it. To find out about the adjacency matrix format type help(mol_matrix)


    """
    atoms = [M[i][0] for i in range(len(M))]
    bonds=list()
    for i in range(len(M)):
        for j in range(i+1, len(M)):
            if M[i][1][j]>0: bonds.append((M[i][1][j], i+1, j+1))
    if pos==list(): return molecule(atoms=atoms, bonds=bonds)
    else: return molecule(atoms=atoms, bonds=bonds, pos=pos)

def mutation(mol):
    """
    def mutation(mol):

    given a openbabel molecule, returns a random mutation using the mutation functions
    """
    return rdLs([mut01, mut02, mut03, mut04, mut05])(mol)
def mut01(mol):
    """
    def mut01(mol):

    given a molecule, returns a new molecule based on it

    It substitutes one atom by another valid one. The "subs" dictionary contains valid
    substitutions. Change the dictionary to change the substitions done by mut01.
    """
    atoms, bonds = mol_info(mol)
    dc = dict()
    for i in range(len(atoms)):
        if dc.get(atoms[i], None)==None: dc[atoms[i]] = list()
        dc[atoms[i]].append(i)

    x = rdLs(list(set(subs.keys())&set(dc.keys())))
    y = rdLs(subs[x])
    z = rdLs(dc[x])
    if x==None or y==None or z==None: return mol
    atoms[z] = y
    return molecule(atoms=atoms, bonds=bonds)
def mut02(mol):
    """
    def mut02(mol):

    given a molecule, returns a new molecule based on it

    It picks one carbon that has at least one hydrogen, deletes the
    hydrogen, replaces the carbon by nitrogen and returns new molecule.
    """
    if type(mol)!=pybel.Molecule: mol = molecule(smile=mol)
    M = mol_matrix(mol)
    CS = [i for i in range(len(M)) if M[i][0]==6]
    HS = [list() for i in range(len(M))]
    for x in CS: HS[x] = [i for i in range(len(M)) if M[x][1][i]>0 and M[i][0]==1]
    CS = [x for x in CS if len(HS[x])>0]
    x=rdLs(CS)
    if x==None: return mol
    H = HS[x][0]
    M[x][0]=7
    M = [[M[i][0], [M[i][1][j] for j in range(len(M)) if j!=H]] for i in range(len(M)) if i!=H]
    return matrix_mol(M)
def mut03(mol):
    """
    def mut03(mol):

    given a molecule, returns a new molecule based on it

    It replaces one hydrogen of a random carbon by a CH3^{-}
    """
    if type(mol)!=pybel.Molecule: mol = molecule(smile=mol)
    M = mol_matrix(mol)
    H = rdLs([i for i in range(len(M)) if M[i][0]==1])
    if H==None: return mol
    M[H][0]=6
    for i in range(len(M)): M[i][1]+=[0]*3
    for i in range(3):M.append([1, [0]*len(M[0][1])])
    M[H][1][-1],M[H][1][-2],M[H][1][-3],M[-1][1][H],M[-2][1][H],M[-3][1][H]=1,1,1,1,1,1
    return matrix_mol(M)
def mut04(mol):
    """
    def mut04(mol):

    given a molecule, returns a new molecule based on it

    It picks two random carbons that have at least one hydrogen, removes them and join them by a bond.
    If the carbons already have a bond, the bond type increases by one
    """
    M = mol_matrix(mol)
    CS = [i for i in range(len(M))]
    HS = [list() for i in range(len(M))]
    for x in CS: HS[x] = [i for i in range(len(M)) if M[x][1][i]>0 and M[i][0]==1]
    CS = [x for x in CS if len(HS[x])>0]
    if len(CS)<2: return mol
    shuffle(CS)
    x1, x2 = CS[0], CS[1]
    H1, H2 = HS[x1][0], HS[x2][0]
    M[x1][1][x2]+=1
    M[x2][1][x1]+=1
    M = [[M[i][0], [M[i][1][j] for j in range(len(M)) if j!=H1 and j!=H2]] for i in range(len(M)) if i!=H1 and i!=H2]
    return matrix_mol(M)
def mut05(mol):
    """
    def mut05(mol):

    given a molecule, returns a new molecule based on it

    It picks a random atom of carbon, nitrogen or oxygen, deletes it and bonds atoms that were together
    with it. If the atoms still need bonds to not be ions, then hydrogens are bonded to them.
    """
    M = mol_matrix(mol)
    N = len(M)
    S = [sum([M[i][1][j] for j in range(N)]) for i in range(N)]
    ls = [6, 7, 8]
    x = rdLs([i for i in range(N) if M[i][0] in ls])
    if x==None: return mol
    con = [i for i in range(len(M)) if M[i][0] in ls and M[x][1][i]>0]
    for i in range(len(con)):
        for j in range(i+1, len(con)):
            m = min(M[con[i]][1][x], M[con[j]][1][x])
            M[con[i]][1][con[j]]=M[con[j]][1][con[i]]=M[con[i]][1][con[j]]+m
            M[con[i]][1][x]=M[x][1][con[i]]=M[x][1][con[i]]-m
            M[con[j]][1][x]=M[x][1][con[j]]=M[x][1][con[j]]-m
    delete = [x]+[i for i in range(len(M)) if M[i][1][x]>0 and M[i][0]==1]
    M = [[M[i][0], [M[i][1][j] for j in range(len(M)) if j not in delete]] for i in range(len(M)) if i not in delete]
    hs = [0 for i in range(len(M))]
    for i in range(len(M)):
        if M[i][0]==6: hs[i]=4-sum([M[i][1][j] for j in range(len(M))])
        if M[i][0]==7: hs[i]=3-sum([M[i][1][j] for j in range(len(M))])
        if M[i][0]==8: hs[i]=sum([M[i][1][j] for j in range(len(M))])%2
    H = sum(hs)
    pos = len(M)
    for i in range(len(M)): M[i][1]+=[0]*H
    M+=[[1, [0]*(len(M)+H)]]*H
    for i in range(len(hs)):
        for j in range(hs[i]):
            M[i][1][pos]=M[pos][1][i]=1
            pos+=1
    return matrix_mol(M)
