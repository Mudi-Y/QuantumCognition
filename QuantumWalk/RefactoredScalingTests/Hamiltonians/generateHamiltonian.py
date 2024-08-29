import math
import pandas as pd
import numpy as np


def genHamiltonian(m=4, integer=True):
    if integer: #integer value hamiltonians based on orig H used in our paper (emailed to me by Raghav)
        #Hamiltonian mxm grid
        h = [[0]*m for _ in range(m)]

        start = 0
        for i in range(m):
            h[i][i] = start 
            start += 2
            if i+1 < m:
                h[i][i+1] = (m-1)*2
            if i-1 >= 0: 
                h[i][i-1] = (m-1)*2
        
        return pd.DataFrame(h)
    
    else: #hamiltonain based on Gunnar's paramaters
        sigma = math.sqrt(1/3)  #sigma squared is diffusion
        delta = 1 #drift
        # m = 4 #number of qubits

        #Hamiltonian mxm grid
        h = [[0]*m for _ in range(m)]

        for i in range(m):
            h[i][i] = (delta*(i+1)/m)
            if i+1 < m:
                h[i][i+1] = sigma**2
            if i-1 >= 0: 
                h[i][i-1] = sigma**2
            
        return pd.DataFrame(h)

for m in [4, 8, 16, 32, 64]:
    h = genHamiltonian(m, integer=True)
    h.to_csv(f'/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Hamiltonians/hamiltonian_{m}.txt', sep=" ", index=False, header=False)
    v,_ = np.linalg.eig(h)
    print(f"HSize {m}, min eig {min(v)}")
