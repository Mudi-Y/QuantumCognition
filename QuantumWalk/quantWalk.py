import networkx as nx
from collections import defaultdict
from itertools import combinations
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import math

m = 103 #as defined by paper, number of evidence states

num_reads = 1000

#paramaters given by Gunnar
sigma = math.sqrt(1/3)
mu = 1
gamma = 8 #lagrange paramater


#unsure how to define delta as specified in paper, made it linear with respect to i
def delta(i):
    return mu*(i)


#Hamiltonian mxm grid
h = [[0]*m for _ in range(m)]

for i in range(m):
    h[i][i] = delta(i+1)*((i+1)/m)
    if i+1 < m:
        h[i][i+1] = sigma**2
    if i-1 >= 0: 
        h[i][i-1] = sigma**2


#create Q matrix for Ising
Q = defaultdict(int)

for i in range(m):
    Q[(i,i)] = mu*sigma*h[i][i]
    if i+1 < m:
        Q[(i,i+1)] = sigma*h[i][i+1]
    if i-1 >= 0:
        Q[(i,i-1)] = sigma*h[i][i-1]

chain_strength = gamma*m

sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_ising(Q, {}, chain_strength=chain_strength, num_reads=num_reads, label='QRW Test')

#get result
sample = response.record.sample[0]
print(sample)