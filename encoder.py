from pickle import dump as pdump
from openbabel import pybel
from pandas import read_csv
from sklearn.feature_extraction.text import CountVectorizer
from MOFLIB import *

df = read_csv('qmof_database/qmof.csv')
df = df[df['info.mofid.smiles_nodes'].notna()]
df = df[df['info.mofid.smiles_linkers'].notna()]
df = df[df['info.mofid.topology'].notna()]
df = df[df['info.mofid.topology']!='TIMEOUT']
df = df[df['info.mofid.topology']!='ERROR,UNKNOWN']
df = df[df['info.mofid.topology']!='MISMATCH']
df = df[df['outputs.pbe.bandgap'].notna()]

nodes = df['info.mofid.smiles_nodes'].tolist()
nodes = [eval(n) for n in nodes]
linkers = df['info.mofid.smiles_linkers'].tolist()
linkers = [eval(n) for n in linkers]
topology = df['info.mofid.topology'].tolist()
bandgap = df['outputs.pbe.bandgap'].tolist()
bandgap = [bandgap[i] for i in range(len(bandgap)) if topology[i] in ['pcu', 'fcu']]
linkers= [linkers[i] for i in range(len(linkers)) if topology[i] in ['pcu', 'fcu']]
nodes = [nodes[i] for i in range(len(nodes)) if topology[i] in ['pcu', 'fcu']]
topology = [topology[i] for i in range(len(topology)) if topology[i] in ['pcu', 'fcu']]
vec = CountVectorizer()
topology = vec.fit_transform(topology).toarray().tolist()
print(vec.get_feature_names_out())
del df

# There are at most 4 different types of nodes and 8 types of linkers
# Any linker has at most 136 atoms
# Any node at most 36
mols = dict()
for n in nodes:
    for e in n:
        if mols.get(e, None)==None: mols[e] = pybel.readstring('smi', e)

for n in linkers:
    for e in n:
        if mols.get(e, None)==None: mols[e] = pybel.readstring('smi', e)
for m in mols: mols[m].addh()

#clf = XGBClassifier(n_estimators=10, max_depth=5, learning_rate=1e-1, objective='binary:logistic')
fingerprint, atoms= dict(), dict()
for m in mols: fingerprint[m] =convert_to_bit(mols[m].calcfp('FP3').bits, 60)+convert_to_bit(mols[m].calcfp('FP4').bits, 310)+convert_to_bit(mols[m].calcfp('MACCS').bits, 170)
for m in mols: atoms[m] = [a.atomicnum for a in mols[m].atoms]
for m in mols: atoms[m] = padding(atoms[m], 140)

encoding = dict()
for m in mols: encoding[m] = atoms[m]+fingerprint[m] # 669 entries
del mols, fingerprint, atoms
N = len(encoding[nodes[0][0]]) # 680 
print(N)
X, Y= list(), list()
for i in range(len(nodes)):
    exs = [encoding.get(n, None)!=None for n in nodes[i]]+[encoding.get(n, None)!=None for n in linkers[i]] #exists
    if False not in exs:
        x = list()
        for n in nodes[i]: x+=encoding[n]
        x+=[0]*(N*(4-len(nodes[i])))
        for n in linkers[i]: x+=encoding[n]
        x+=[0]*(N*(8-len(linkers[i])))
        x+=topology[i]
        X.append(x)
        Y.append(int(bandgap[i]>1.8 and bandgap[i]<2.8))
with open('encoding.bin', 'wb+') as f: pdump([X, Y], f)
