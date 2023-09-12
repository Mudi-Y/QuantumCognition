import numpy as np
import sys
import sympy as sp
import pickle

from ConstructBQM import buildBQM
from anneal import AnnealBQM

class AnnealConfig:
    """Annealing paramaters. 
        Default values: method = 'simulated', schedule = [(0.0,0.0),(10.0,0.4),(15.0,0.4),(30.0,0.8),(35.0,0.8),(65.0,1.0)], groups = [2,3,4,5,6], runs = 3, trial = 1, output = '/Results/{method}{trial}."""
    def __init__(self, method = "simulated", 
        schedule = [(0.0,0.0),(10.0,0.4),(15.0,0.4),(30.0,0.8),(35.0,0.8),(65.0,1.0)], 
        groups = [2,3,4,5,6],
        runs = 3,
        trial = 1,
        output= None):
        self.m_method = method
        self.m_schedule = schedule
        self.m_groups = groups
        self.m_runs = runs
        self.m_trial = trial
        self.m_output = output if output else f'Results/{self.m_method}{self.m_trial}/'


class AnnealHandler:
    """Class to handle all information relevant to annealing. Stores problem Hamiltonian and AnnealConfig. Constructs BQM, anneals, and plots results.
        Default values: hamiltonian = [[0, 6, 0, 0], [6, 2, 6., 0], [0, 6, 4, 6], [0, 0, 6, 6.0]], AnnealConfig = AnnealConfig class defaults."""
    def __init__(self, hamiltonian = [[0, 6, 0, 0], [6, 2, 6., 0], [0, 6, 4, 6], [0, 0, 6, 6.0]],
        AnnealConfig = AnnealConfig()):
        self.m_hamiltonian = np.array(hamiltonian)
        self.m_AnnealConfig = AnnealConfig
    
    def ConstructBQMs(self):
        """Construct BQMs from provided hamiltonian and AnnealConfig paramaters. Save in the /BQM/ directory in both .pkl and human readable formats."""
        for group in self.m_AnnealConfig.m_groups:
            BQMPath = f'Results/BQM/BQM_group{group}'
            valsPath = f'Results/vals/vals_group{group}'
            buildBQM(group, self.m_hamiltonian, BQMPath, valsPath)


    def Anneal(self):
        for group in self.m_AnnealConfig.m_groups:
            BQMPath = f'Results/BQM/BQM_group{group}'
            valsPath = f'Results/vals/vals_group{group}'
            for run in range(1, self.m_AnnealConfig.m_runs+1):
                outpath = self.m_AnnealConfig.m_output+f'{self.m_AnnealConfig.m_method}_group{group}_attempt{run}'
                AnnealBQM(BQMPath, valsPath, self.m_AnnealConfig.m_method, self.m_AnnealConfig.m_schedule, outpath)


    def Plot(self):
        """Generate plots form annealing results. Save in /Plots/ directory."""
        pass