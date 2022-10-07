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
import neal

num_reads = 100
time = 1000 #range between 0.5 and 2000 micro seconds

m = 4 #as defined by paper, number of evidence states

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

# linear = {('q1','q1'):0.25, ('q2','q2'):0.5, ('q3','q3'):0.75, ('q4','q4'):1}
# qudratic = {('q1','q2'):2.0/3.0, ('q2','q3'):2.0/3.0, ('q3','q4'):2.0/3.0}

qubo = {**linear, **qudratic}

#####################sample with Qudratic Sampler
sampler = neal.SimulatedAnnealingSampler()
#range between 0.5 and 2000.0 micro seconds
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

plt.bar(x, y)
plt.xlabel("Final State")
plt.ylabel("Count")
plt.title(f"Final States after {num_feasible} Feasible Samples on Simulated Annealer")
plt.xticks(rotation=0)
plt.savefig('QuantumWalk/SimulatedAnneal/Plots/test.png')