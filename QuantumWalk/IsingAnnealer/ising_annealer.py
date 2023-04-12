# /usr/local/bin/python -i /workspace/QuantumCognition/QuantumWalk/IsingAnnealer/ising_annealer.py

from ising_model import IsingModel
a = IsingModel(q=2, t=1, sigmaSq=2, delta=4)

# res = a.annealQPU()
res = a.calculateEnergy([1,-1,1,-1])
print(res, end = "\n")