import math
import pandas as pd

def genHamiltonian(m=4):

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

for m in [4, 8, 16, 32, 64]:
    h = genHamiltonian(m)
    h.to_csv(f'/workspaces/QuantumCognition/QuantumWalk/ScalingTests/hamiltonians/hamiltonian_{m}.txt', sep=" ", index=False, header=False)