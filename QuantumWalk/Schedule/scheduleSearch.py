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
import pandas as pd

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

def get_counts(set):
            ret = {}
            num_feasible = 0
            for entry in set.record:
                if sum(entry[0]) == 1: #filter for only valid solutions
                    ret[str(entry[0])] = entry[2]
                    num_feasible += entry[2]
            return ret, num_feasible

percents = [0.30, 0.35, 0.40, 0.45, 0.50]
pauses = [5, 40, 160, 320]

data = pd.DataFrame(columns=['[0 0 0 1]', '[0 0 1 0]', '[0 1 0 0]', '[1 0 0 0]'])


for i in range(3):
    for percent in percents:
        for pause in pauses:
            
            schedule = [[0.0, 0.0], [10, percent], [10+pause, percent], [20+pause, 1.0]]

            qsampleset = sampler.sample_qubo(qubo, num_reads=num_reads, anneal_schedule= schedule)
            results_count, num_feasible = get_counts(qsampleset)

            lists = sorted(results_count.items())
            x, y = zip(*lists)

            data.loc[f'P{pause}_Per{percent}_T{i+1}'] = y

            plt.bar(x, y)
            plt.xlabel("Final State")
            plt.ylabel("Count")
            plt.title(f"Final States: P {pause}, Per {percent}, T {i+1}")
            plt.xticks(rotation=0)
            plt.savefig(f'/workspace/QuantumCognition/QuantumWalk/Schedule/ScheduleSearchPlots/T2/P{pause}_Per{percent}_T{i+1}.png')
            plt.clf()

data.to_pickle('/workspace/QuantumCognition/QuantumWalk/Schedule/ScheduleSearck.pkl')


