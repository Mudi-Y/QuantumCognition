import networkx as nx
from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system import LeapHybridSampler, LeapHybridCQMSampler
from dwave.system.composites import EmbeddingComposite
from dimod import ConstrainedQuadraticModel, Binary
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

import pickle


m = 4 #as defined by paper, number of evidence states
num_reads = 1000

#paramaters given by Gunnar
sigma = math.sqrt(2)  #sigma squared is diffusion
delta = 4 #drift
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

#####################sample with Qudratic Sampler
sampler = EmbeddingComposite(DWaveSampler())
out = []
tmax = 15

count = 0
for time in range(0, tmax+1, tmax//15):
    if time == 0:
        time = 1

    qsampleset = sampler.sample_qubo(qubo, num_reads=num_reads, annealing_time = time)

    def get_counts(set):
        ret = {}
        num_feasible = 0
        for entry in set.record:
            if sum(entry[0]) == 1: #filter for only valid solutions
                ret[str(entry[0])] = entry[2]
                num_feasible += entry[2]
        return ret, num_feasible

    results_count, num_feasible = get_counts(qsampleset)

    lists = sorted(results_count.items())

    x, y = zip(*lists)
    res = (list(x),(list(y)))
    out.append(res)

with open(f'stepped_anneal_{tmax}_drift{delta}_diff{int(sigma**2)}.pickle', 'wb') as handle:
    pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)


# with open(f'stepped_anneal_{tmax}_drift{delta}_diff{int(sigma**2)}.pickle', 'rb') as handle:
#     b = pickle.load(handle)

# print(len(b))