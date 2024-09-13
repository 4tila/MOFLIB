from numpy import exp
from MOFLIB import *
from pickle import load, dump

def save_file(P, edge):
    if P>0.5 and edge.charge==0 and '.' not in str(edge):
        edge.make3D(steps=2)
        edge.localopt()
        pos = mol_pos(edge)
        atoms, bonds = mol_info(edge)
        fname = './molecules/'+('%.4f'%(P)).replace('0.', '')+'_'+random_string()
        print(fname)
        with open(fname,'wb+' ) as f: dump([atoms, bonds, pos], f)
def encode(mol): return padding([a.atomicnum for a in mol.atoms], 140)+convert_to_bit(mol.calcfp('FP3').bits, 60)+convert_to_bit(mol.calcfp('FP4').bits, 310)+convert_to_bit(mol.calcfp('MACCS').bits, 170)
with open('./model01.mod', 'rb') as f: clf = load(f)
with open('./GDB/gdb11_size07.smi', 'r') as f: smile = rdLs(f.read().split('\n'))
topologies = ['fcu', 'pcu']
topologies = {topologies[i]:i for i in range(len(topologies))}
top_str = 'pcu'
topology = [0]*2
topology[topologies[top_str]]=1
N=680

node = molecule(smile='[Mn]')
aux_mol = molecule(smile=smile)
edgeM = dicarboxylate_linker(aux_mol)
edge = matrix_mol(edgeM)
code_node, code_edge = encode(node), encode(edge)
P = clf.predict_proba([code_node+[0]*N*3+code_edge+[0]*N*7+topology])[0][1]
save_file(P, edge)

KT = 0.1

for _ in range(20):
    aux_mol2 = mutation(aux_mol)
    edgeM2 = dicarboxylate_linker(aux_mol2)
    edge2 = matrix_mol(edgeM2)
    if edge2==None: exit(0)
    code_edge2 = encode(edge2)
    P2 = clf.predict_proba([code_node+[0]*N*3+code_edge2+[0]*N*7+topology])[0][1]
    if exp((P2-P)/KT)>random():
        edge, P,  aux_mol, edgeM = edge2, P2,aux_mol2, edgeM2
        print('='*100)
        print(edge)
        print(P)
        save_file(P, edge)
