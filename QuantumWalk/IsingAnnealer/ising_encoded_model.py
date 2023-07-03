import numpy as np
import math as math
import matplotlib.pyplot as plt
from scipy import linalg
from numpy.core.numeric import array_equal
from dwave.system import DWaveSampler, EmbeddingComposite

#based on ising_mode.py. Differences are in LinearTerms and QudraticTerms. Encoding based on binary expansion.

class IsingModel:
  "create a quantum circuit with initial state Init, q qubits, time=t (for generating U), diffusion=sigmaSq, drift=delta"
  def __init__(self, q=4, Init=None, t=1, sigmaSq=(1/3), delta=1):
    self.q = q
    self.t = t
    self.sigma = np.sqrt(sigmaSq)
    self.sigmaSq = sigmaSq
    self.delta = delta
    self.init_H()
    self.states = self.genStates()
    self.isingModel = self.gen_ising()
  

  def init_H(self):
    #State Hamiltonian
    m = 2**self.q #number of possible states
    #default paramaters given by Gunnar
    delta = self.delta #drift
    sigmaSq = self.sigmaSq #diffusion

    #Hamiltonian mxm grid
    h = [[0]*m for _ in range(m)]
    for i in range(m):
        h[i][i] = (delta*(i+1)/m)
        if i+1 < m:
            h[i][i+1] = sigmaSq
        if i-1 >= 0: 
            h[i][i-1] = sigmaSq
    self.H = np.array(h)


  def gen_ising(self):
    "create Insing formulation in accordance with definiton discussed with Thi Ha"
    linear = self.linearTerms()
    quad = self.quadraticTerms()
    return (linear, quad)

  
  def linearTerms(self):
    "create a dictionary of the linear terms"
    linear = {}

    for state in self.states:
      diag = int(state, 2)

      # key = (state, state)
      key = state
      # multiply by two because in qudratic terms we count both (a,b) and (b,a) interactions
      value = 2 * self.H[diag][diag]
      if value != 0:
        linear[key+"a"] = value*(1/2)
        linear[key+"a"] = value*(1/4)
    return linear
  

  def quadraticTerms(self):
    "create a dictionary of qudradic terms"
    quad = {}

    for state1i in range(2**self.q):
      for state2i in range(2**self.q):
        state1str = self.states[state1i]
        state2str = self.states[state2i]

        key1 = (state1str+"a", state2str+"a")
        key2 = (state1str+"b", state2str+"b")
        value = self.H[state1i][state2i]
        if value != 0:
          quad[key1] = value*(1/2)
          quad[key2] = value*(1/4)
    return quad


  def annealQPU(self):
    "anneal on QPU"
    sampler = EmbeddingComposite(DWaveSampler())
    lin, quad = self.isingModel
    sampleset = sampler.sample_ising(lin, quad, num_reads=1000)
    return sampleset


  def eig(self):
    l, v = np.linalg.eig(self.H)
    self.elambdas = l
    self.evectors = v
    print(f"Lambdas: \n{l}", end="\n\n")
    print(f"Vectors: \n{v}",end="\n\n")
    return (l, v)

  
  def genStates(self):
    "generate binary string representation of quantum states"
    n = 2**self.q
    ret = []
    for i in range(0,n):
      ret.append(bin(i)[2:].zfill(self.q))
    return ret

  
  def calculateEnergy(self, spins):
    "given state as a list of spins for each state, calculate energy of Ising formulation"
    spins = {state : energy for state, energy in zip(self.states, spins)}
    ret = 0

    #linear terms
    for state, coef in self.isingModel[0].items():
      spin = spins[state]
      ret += spin*coef

    #qudratic terms
    for state, coef in self.isingModel[1].items():
      spin1 = spins[state[0]]
      spin2 = spins[state[1]]
      mult = spin1*spin2
      ret += mult*coef
    
    return ret

