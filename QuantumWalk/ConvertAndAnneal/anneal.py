import numpy as np
import sympy as sp
import time
import sys
import pickle
from sympy.parsing.sympy_parser import parse_expr

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

from dimod import BinaryQuadraticModel

import neal
import dimod

from misc import HamiltonianToPauli, PauliToString

file = open('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/bqm', 'rb')
bqm = pickle.load(file)
file.close()

file = open('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/vals', 'rb')
vals = pickle.load(file)
Sign_vals, num_new_groups, num_orig_vars, Cp_num_exp, orig_vars, L = vals
file.close()


#D-Wave
if sys.argv[1] == "simulated":
    sampler = neal.SimulatedAnnealingSampler()
    print("##########\n","Simulated Annealing")
else:
    sampler = EmbeddingComposite(DWaveSampler())
    print("##########\n","QPU")

#Solve the BQM
num_reads = 100
start = time.time()
anneal_schedule = []
# anneal_schedule = [(0.0,0.0),(10.0,0.25),(20,0.25),(30.0,0.5),(40.0,0.5),(50.0,0.75),(60.0,0.75),(70.0,1.0)]
# anneal_schedule = [(0.0,0.0),(10.0,0.4),(115.0,0.4),(125.0,1.0)]
# anneal_schedule = [(0.0,0.0),(5.0,0.4),(35.0,0.4),(40.0,0.8),(50.0,0.8),(55.0,1.0)]
# anneal_schedule = [(0.0,0.0),(10.0,0.4),(15.0,0.4),(30.0,0.8),(35.0,0.8),(65.0,1.0)]
anneal_schedule = [(0.0, 0.0), (5.0, 0.4), (15.0, 0.4), (25.0, 0.8), (35.0, 0.8), (55.0, 1.0)]
sampleset = sampler.sample(bqm, num_reads=num_reads, num_sweeps = 100) if sys.argv[1] == "simulated" else (sampler.sample(bqm, num_reads=num_reads, anneal_schedule = anneal_schedule) if anneal_schedule else sampler.sample(bqm, num_reads=num_reads))
print("Anneal Schedule: ", anneal_schedule)

end = time.time()
print("Number of Qubits: ", len(sampleset.first.sample))
print("Time: ", end-start)
print("Num Reads: ", num_reads)
print("Info: ", sampleset.info)
#Get minimum energy record
solution = sampleset.first
# print(solution.energy, solution.sample)
print("Energy: ", sampleset.first.energy)
print("Sample: ", sampleset.first.sample)

#Get the data
bin_vals_from_var = solution.sample
spin_vals_from_var = {}
qubit_vals = []
for bin_var in bin_vals_from_var:
    spin_vals_from_var['z'+bin_var[1:]] = 2*bin_vals_from_var[bin_var] -1
    if bin_vals_from_var[bin_var]  == 1:
        qubit_vals.append('0')
    else:
        qubit_vals.append('1')
print("Spin Values: ",spin_vals_from_var)
print("New Eigenvector: ", "".join([str(q) for q in qubit_vals[:num_orig_vars]]))
print("New Cost: ", solution.energy)

Cp_num = Cp_num_exp
for spin_var in orig_vars:
    Cp_num = Cp_num.subs(spin_var, spin_vals_from_var[str(spin_var)])
    #print(spin_var,spin_vals_from_var[str(spin_var)])
print("Sum of bi^2: ", Cp_num) #This is sum of bi^2

experimental_eigenvalue = float(str(L)) + solution.energy / Cp_num
print("Starting Lambda: ", L)
print("Estimated Eigenvalue: ", experimental_eigenvalue)
#print(np.linalg.eig(Hamiltonian))
I = np.eye(2, 2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0,-1]], dtype=complex)
Y = complex(0,-1)*np.matmul(Z,X)
v,d = np.linalg.eig(3*np.kron(I,I)+6*np.kron(I,X)-1*np.kron(I,Z)+3*np.kron(X,X) +3*np.kron(Y,Y)-2*np.kron(Z,I))
print("Eigen Values: ",v)
print("Eigen Vectors: ",d)
print("Absolute Error: ", abs(min(v) - experimental_eigenvalue))
print("Groups: ", num_new_groups)
print("Signs: ", Sign_vals)
