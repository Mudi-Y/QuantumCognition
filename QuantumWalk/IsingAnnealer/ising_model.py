import numpy as np
import math as math
import matplotlib.pyplot as plt
from scipy import linalg
from numpy.core.numeric import array_equal

class IsingModel:
  "create a quantum circuit with initial state Init, q qubits, time=t (for generating U), diffusion=sigmaSq, drift=delta"
  def __init__(self, q=4, Init=None, t=1, sigmaSq=(1/3), delta=1):
    self.q = q
    self.t = t
    self.sigma = np.sqrt(sigmaSq)
    self.sigmaSq = sigmaSq
    self.delta = delta
    self.init_H()
    self.init_U()
    self.states = self.genStates()
    self.isingModel = self.gen_ising()
  


  def init_H(self):
    #State Hamiltonian
    m = 2**self.q #number of possible states
    #default paramaters given by Gunnar
    delta = self.delta #drift
    sigmaSq = self.sigmaSq #diffusion
    sigma = self.sigma  #sigma, not really used

    #Hamiltonian mxm grid
    h = [[0]*m for _ in range(m)]
    for i in range(m):
        h[i][i] = (delta*(i+1)/m)
        if i+1 < m:
            h[i][i+1] = sigmaSq
        if i-1 >= 0: 
            h[i][i-1] = sigmaSq
    self.H = np.array(h)


  def init_U(self):
    self.U = linalg.expm(-(1j)*self.t*self.H)


  def gen_ising(self):
    linear = self.linearTerms()
    quad = self.quadraticTerms()
    return (linear, quad)

  
  def linearTerms(self):
    linear = {}

    for state in self.states:
      diag = int(state, 2)

      key = (state, state)
      value = self.U[diag][diag]
      linear[key] = value
    return linear
  
  def quadraticTerms(self):
    quad = {}

    for state1i in range(2**self.q):
      for state2i in range(state1i+1, 2**self.q):
        state1str = self.states[state1i]
        state2str = self.states[state2i]

        key = (state1str, state2str)
        value = self.U[state1i][state2i]
        quad[key] = value
    return quad


  def annealQPU(self):
    from dwave.system import DWaveSampler, EmbeddingComposite
    sampler = EmbeddingComposite(DWaveSampler())
    lin, quad = self.isingModel
    sampleset = sampler.sample_ising(lin, quad, num_reads=1000)
    print(sampleset)
  

  def annealHybrid(self):
    return


  def eig(self):
    l, v = np.linalg.eig(self.H)
    self.elambdas = l
    self.evectors = v
    print(f"Lambdas: {l}")
    print(f"Vectors: {v}")
    return (l, v)

  
  def genStates(self):
    "generate binary string representation of quantum states"
    n = 2**self.q
    ret = []
    for i in range(0,n):
      ret.append(bin(i)[2:].zfill(self.q))
    return ret