import numpy as np
import time
import os
import pickle

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

import neal
import dimod

def AnnealBQM(BQMPath, valsPath, method, schedule, outpath):

    #create log file
    logfile = outpath
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    logfile = open(logfile, 'w')

    file = open(BQMPath, 'rb')
    bqm = pickle.load(file)
    file.close()

    file = open(valsPath, 'rb')
    vals = pickle.load(file)
    Sign_vals, num_new_groups, num_orig_vars, Cp_num_exp, orig_vars, L = vals
    file.close()


    #D-Wave
    if method == "simulated":
        sampler = neal.SimulatedAnnealingSampler()
        logfile.write("##########\n")
        logfile.write("Simulated Annealing")
    else:
        sampler = EmbeddingComposite(DWaveSampler())
        logfile.write("##########\n")
        logfile.write("QPU")

    #Solve the BQM
    num_reads = 100
    start = time.time()
    # anneal_schedule = [(0.0,0.0),(10.0,0.25),(20,0.25),(30.0,0.5),(40.0,0.5),(50.0,0.75),(60.0,0.75),(70.0,1.0)]
    # anneal_schedule = [(0.0,0.0),(10.0,0.4),(115.0,0.4),(125.0,1.0)]
    # anneal_schedule = [(0.0,0.0),(5.0,0.4),(35.0,0.4),(40.0,0.8),(50.0,0.8),(55.0,1.0)]
    # anneal_schedule = [(0.0,0.0),(5.0,0.4),(35.0,0.4),(40.0,1.0)] #original suggested by Raghav
    # anneal_schedule = [(0.0,0.0),(10.0,0.4),(15.0,0.4),(30.0,0.8),(35.0,0.8),(65.0,1.0)]
    # anneal_schedule = [(0.0, 0.0), (5.0, 0.4), (15.0, 0.4), (25.0, 0.8), (35.0, 0.8), (55.0, 1.0)]
    anneal_schedule = schedule
    sampleset = sampler.sample(bqm, num_reads=num_reads, num_sweeps = 100) if method == "simulated" else (sampler.sample(bqm, num_reads=num_reads, anneal_schedule = anneal_schedule) if anneal_schedule else sampler.sample(bqm, num_reads=num_reads))
    logfile.write("Anneal Schedule:  ")
    logfile.write(str(anneal_schedule)+'\n')

    end = time.time()
    logfile.write("Number of Qubits:  ")
    logfile.write(str(len(sampleset.first.sample))+'\n')
    
    logfile.write("Time:  ")
    logfile.write(str(end-start)+'\n')

    logfile.write("Num Reads:  ")
    logfile.write(str(num_reads)+'\n')

    logfile.write("Info:  ")
    logfile.write(str(sampleset.info)+'\n')
    
    #Get minimum energy record
    solution = sampleset.first
    
    # print(solution.energy, solution.sample)
    logfile.write("Energy:  ")
    logfile.write(str(sampleset.first.energy)+'\n')


    logfile.write("Sample:  ")
    logfile.write(str(sampleset.first.sample)+'\n')

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
    logfile.write("Spin Values:  ")
    logfile.write(str(spin_vals_from_var)+'\n')
    logfile.write("New Eigenvector:  ")
    logfile.write(str("".join([str(q) for q in qubit_vals[:num_orig_vars]]))+'\n')
    logfile.write("New Cost:  ")
    logfile.write(str(solution.energy)+'\n')

    Cp_num = Cp_num_exp
    for spin_var in orig_vars:
        Cp_num = Cp_num.subs(spin_var, spin_vals_from_var[str(spin_var)])
        #print(spin_var,spin_vals_from_var[str(spin_var)])
    logfile.write("Sum of bi^2:  ")
    logfile.write(str(Cp_num)+'\n')

    experimental_eigenvalue = float(str(L)) + solution.energy / Cp_num
    logfile.write("Starting Lambda:  ")
    logfile.write(str(L)+'\n')

    logfile.write("Estimated Eigenvalue:  ")
    logfile.write(str(experimental_eigenvalue)+'\n')

    #print(np.linalg.eig(Hamiltonian))
    I = np.eye(2, 2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.array([[1, 0], [0,-1]], dtype=complex)
    Y = complex(0,-1)*np.matmul(Z,X)
    v,d = np.linalg.eig(3*np.kron(I,I)+6*np.kron(I,X)-1*np.kron(I,Z)+3*np.kron(X,X) +3*np.kron(Y,Y)-2*np.kron(Z,I))
    logfile.write("Eigen Values:  ")
    logfile.write(str(v)+'\n')

    logfile.write("Eigen Vectors:  ")
    logfile.write(str(d)+'\n')

    logfile.write("Absolute Error:  ")
    logfile.write(str(abs(min(v) - experimental_eigenvalue))+'\n')

    logfile.write("Groups:  ")
    logfile.write(str(num_new_groups)+'\n')

    logfile.write("Signs:  ")
    logfile.write(str(Sign_vals)+'\n')
    
    logfile.close()

