# /usr/local/bin/python -i /workspace/QuantumCognition/QuantumWalk/IsingAnnealer/ising_annealer.py

from ising_encoded_model import IsingModel
a = IsingModel(q=2, t=1, sigmaSq=2, delta=4)

res = a.annealQPU()
print(res, end = "\n")
a.eigs()


# 1/2|00> + 1/4|01> + 1/2|10> + 1/4|11>