import math
import pandas as pd
import numpy as np
import scipy as sp

def getHamiltonianMatrix(numQubits,drift,diffusion):
    numStates = 2**numQubits

    diag = np.array([drift * stateNum for stateNum in range(numStates)],dtype="float64")
    diag_1 = diffusion*np.ones(numStates-1)
    hamiltonian = np.diagflat(diag) + np.diagflat(diag_1,k=1) + np.diagflat(diag_1,k=-1)

    if not sp.linalg.ishermitian(hamiltonian):
        raise Exception("Hamiltonian is not Hermitian")

    return pd.DataFrame(hamiltonian).astype(int)


drift = 2
diffusion = 6

for m in [2,3,4]:
    h = getHamiltonianMatrix(m, drift, diffusion)
    h.to_csv(f'/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Hamiltonians/hamiltonian_{2**m}.txt', sep=" ", index=False, header=False)