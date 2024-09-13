from xgboost import XGBClassifier
from pickle import load, dump
clf = XGBClassifier(n_estimators=30, max_depth=5, learning_rate=1e-1, objective='binary:logistic')
with open('encoding.bin', 'rb') as f: X, Y = load(f)
N = len(X)*9//10
clf.fit(X[0:N], Y[0:N])
X, Y = X[N:], Y[N:]
Y2 = clf.predict(X).tolist()
accuracy = sum([int(Y[i]==Y2[i]) for i in range(len(Y))])/len(Y)
with open('model01.mod', 'wb+') as f: dump(clf, f)
with open('model01.acc', 'w+') as f: f.write("%.4f\n"%(accuracy))

