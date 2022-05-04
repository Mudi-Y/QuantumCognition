import networkx as nx
from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system import LeapHybridSampler
from dwave.system.composites import EmbeddingComposite
import math
import matplotlib.pyplot as plt
import numpy as np


m = 4 #as defined by paper, number of evidence states

num_reads = 1000

#paramaters given by Gunnar
sigma = math.sqrt(1/3)  #sigma squared is diffusion
delta = 1 #drift
m = 4 #number of qubits

#Hamiltonian mxm grid
h = [[0]*m for _ in range(m)]

for i in range(m):
    h[i][i] = (delta*(i+1)/m)
    if i+1 < m:
        h[i][i+1] = sigma**2
    if i-1 >= 0: 
        h[i][i-1] = sigma**2


linear = {}
qudratic = {('q1','q2'):2, ('q2','q3'):2, ('q3','q4'):2}

qubo = {**linear, **qudratic}

# sampler = LeapHybridSampler()
# qsampleset = np.array([0]*m)
# for i in range(2): 
#     out = sampler.sample_qubo(qubo)
#     qsampleset += out.record.sample[0]

# print(qsampleset)

sampler = EmbeddingComposite(DWaveSampler())
qsampleset = sampler.sample_qubo(qubo, num_reads=num_reads)



#plot results
def calcy (set):
    ret = np.array([0]*len(set.record.sample[0]))
    for entry in set.record:
        ret += entry[0]*entry[2]
    return ret

x = list(range(1, m+1))
y = calcy(qsampleset)

plt.bar(x, y, align='center')
plt.xticks(x, x)
plt.xlabel("Qubit")
plt.ylabel("1 Count")
plt.savefig('/workspace/QuantumCognition/fig1')

