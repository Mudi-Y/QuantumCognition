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

#####################sample with Qudratic Sampler
sampler = EmbeddingComposite(DWaveSampler())
qsampleset = sampler.sample_qubo(qubo, num_reads=num_reads)

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
plt.title(f"Final States after {num_feasible} Feasible Samples on QPU")
plt.xticks(rotation=0)
plt.savefig('/workspace/QuantumCognition/fig3_qpu')



# ######################sample with CQM Model using hybrid sampler
# cqm = ConstrainedQuadraticModel()
# system_state = [Binary(f'qubit_{i}_state') for i in range(m)]

# #obj: [ 4 x_0*x_1 + 4 x_1*x_2 + 4 x_2*x_3 ]/2 -2
# cqm.set_objective((4*system_state[0]*system_state[1] + 4*system_state[1]*system_state[2] + 4*system_state[2]*system_state[3])/2 - 2)
# sum_to_one = cqm.add_constraint(sum(system_state) == 1, label='sum to one')

# sampler = LeapHybridCQMSampler()
# sampleset = sampler.sample_cqm(cqm, time_limit=5)

# feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)
# num_feasible = len(feasible_sampleset.record)
# results = [str(sampleset.record[i][0]) for i in range(num_feasible)]
# results.sort()

# #plot histogram of results
# results_count = Counter(results)
# plt.bar(results_count.keys(), results_count.values())
# plt.xlabel("Final State")
# plt.ylabel("Count")
# plt.title(f"Final States after {num_feasible} Feasible Samples on LeapHybridCQMSampler")
# plt.xticks(rotation=0)
# plt.savefig('/workspace/QuantumCognition/fig3_hybrid')



