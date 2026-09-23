import numpy as np
import pandas as pd

# read the dataset
data = pd.read_csv("hollins.dat", header=None, sep=r"\s+", engine="python", dtype=str)

nr_pages = np.array(data.iloc[0,0], dtype=int)
nr_edges = np.array(data.iloc[0,1], dtype=int)

pages = data.iloc[1:nr_pages+1,:]
edges = data.iloc[1+nr_pages:1+nr_pages+nr_edges,:]

pageslink = np.array(pages.iloc[:,1], dtype=str)

edgesnp = np.array([edges.iloc[:,0], edges.iloc[:,1]], dtype=int).T
edges0 = np.vstack([edgesnp[:,0]-1, edgesnp[:,1]-1 ]).T

adj_list = [[] for _ in range(nr_pages)]
for u, v in edges0:
    adj_list[u].append(v)

backlinks = [[] for _ in range(nr_pages)]
for u, v in edges0:
    backlinks[v].append(u)

# check for the dangling nodes
dangling_nodes_counter = 0
e_dangling_nodes = np.zeros(nr_pages)
for i in range (nr_pages):
    if len(adj_list[i]) == 0:
        dangling_nodes_counter += 1
        e_dangling_nodes[i]=1

print(dangling_nodes_counter, e_dangling_nodes)

Nj = np.zeros(nr_pages, dtype=int)
for i in range(nr_pages):
    Nj[i] = len(adj_list[i])

#CSR format
AA = []
JA = []
IA = [0]

for i in range(nr_pages):
    nz = 0
    for k in backlinks[i]:
        if Nj[k]!=0:
            val = 1.0/Nj[k]
            AA.append(val)
            JA.append(k)
            nz += 1
    IA.append(IA[-1] + nz)

def CSR_matvec(AA, JA, IA, v):
    n = len(IA)-1
    w = np.zeros(n)
    for i in range(n):
        for k in range(IA[i],IA[i+1]):
            w[i] += AA[k]*v[JA[k]]

    return w

def PowerMethod(AA, JA, IA, m, maxIter, relTol):
    inv_n = 1.0/nr_pages
    teleport = m*inv_n
    v0 = np.full(nr_pages, inv_n)

    for k in range(maxIter):
        v = (1-m)*CSR_matvec(AA,JA,IA,v0)
        v += (1-m)*inv_n*(np.dot(e_dangling_nodes,v0))
        v += teleport

        v = v/v.sum()

        if np.linalg.norm(v-v0,1)<relTol:
            return v, k
        else:
            v0 = v

    return v, k

import time
start = time.time()

m=0.15
maxIter=100
relTol=1e-4

v, k = PowerMethod(AA,JA,IA,m,maxIter,relTol)
end = time.time()
print("Executing time {:.6f} seconds".format(end - start))

print("Number of iterations: ", k)

rankingord = np.argsort(v)[::-1]

vfinord= v[rankingord]
pageslinkord = pageslink[rankingord]

print(f"Top 10 PageRank scores:\n {vfinord[:10]}")
print(f"Corresponding page links:\n {pageslinkord[:10]}")