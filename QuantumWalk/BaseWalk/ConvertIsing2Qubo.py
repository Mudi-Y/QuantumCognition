import dimod
from dwave.system import DWaveSampler, LeapHybridSampler, EmbeddingComposite
import numpy as np
import matplotlib.pyplot as plt

m=4
num_reads = 1000

linear = {('q1','q1'):0.25, ('q2','q2'):0.5, ('q3','q3'):0.75, ('q4','q4'):1}
linear = {'q1':0.25, 'q2':0.5, 'q3':0.75, 'q4':1}
qudratic = {('q1','q2'):2.0/3.0, ('q2','q3'):2.0/3.0, ('q3','q4'):2.0/3.0}

# qubo = dimod.ising_to_qubo(linear, qudratic) #converting ising to qubo
# qubo = {**linear, **qudratic}

sampler = LeapHybridSampler()
isampleset = sampler.sample_ising(linear, qudratic)
# qsampleset = sampler.sample_qubo(qubo)

# sampler = EmbeddingComposite(DWaveSampler())
# isampleset = sampler.sample_ising(linear, qudratic, num_reads=num_reads)
# qsampleset = sampler.sample_qubo(qubo, num_reads=num_reads)



qubo =({(('q1', 'q1'), ('q1', 'q1')): 0.5, (('q2', 'q2'), ('q2', 'q2')): 1.0, (('q3', 'q3'), ('q3', 'q3')): 1.5, (('q4', 'q4'), ('q4', 'q4')): 2.0, ('q1', 'q2'): 2.6666666666666665, ('q1', 'q1'): -1.3333333333333333, ('q2', 'q2'): -2.6666666666666665, ('q2', 'q3'): 2.6666666666666665, ('q3', 'q3'): -2.6666666666666665, ('q3', 'q4'): 2.6666666666666665, ('q4', 'q4'): -1.3333333333333333}, -0.5)
