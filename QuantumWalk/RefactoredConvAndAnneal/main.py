import numpy as np
import pickle

from ConstructBQM import buildBQM
from anneal import AnnealBQM
from plots import GeneratePlots

class AnnealConfig:
    """Annealing paramaters. 
        Default values: method = 'simulated', schedule = [(0.0,0.0),(10.0,0.4),(15.0,0.4),(30.0,0.8),(35.0,0.8),(65.0,1.0)], 
        groups = [2,3,4], num_runs = 3, trial = 1, lambdaVal = -3.35, output = {working directory}
        results in {output}/Results/Trial_{trial}."""
    def __init__(self, 
        method = "simulated", 
        schedule = [(0.0,0.0),(10.0,0.4),(15.0,0.4),(30.0,0.8),(35.0,0.8),(65.0,1.0)], 
        groups = [2,3,4],
        num_runs = 3,
        trial = 1,
        lambdaVal = -3.35,
        output = None):

        self.m_method = method
        self.m_schedule = schedule
        self.m_groups = groups
        self.m_runs = num_runs
        self.m_trial = trial
        self.m_lambdaVal = lambdaVal

        #output file path processing
        if output:
            self.m_output = output if output[-1] =="/" else output + "/"
        else: 
            self.m_output = ""
        self.m_output += f'Results/Trial_{self.m_trial}/'


    def ReadConfig(self, file):
        import configparser, ast
        configParser = configparser.RawConfigParser()   
        configParser.read(file)

        self.m_method = configParser.get("configuration","method")
        self.m_method = ast.literal_eval(self.m_method)
        self.m_schedule = configParser.get("configuration","schedule")
        self.m_schedule = ast.literal_eval(self.m_schedule)
        self.m_groups = configParser.get("configuration","groups")
        self.m_groups = ast.literal_eval(self.m_groups)
        self.m_runs = configParser.getint("configuration", "num_runs")
        self.m_trial = configParser.getint("configuration","trial")
        self.m_lambdaVal = configParser.getfloat("configuration","lambdaVal")
        self.m_output = f'Results/Trial_{self.m_trial}/'


class AnnealHandler:
    """Class to handle all information relevant to annealing. Stores problem Hamiltonian and AnnealConfig. Constructs BQM, anneals, and plots results.
        Default values: hamiltonian = [[0, 6., 0, 0], [6., 2., 6., 0], [0, 6., 4., 6.], [0, 0, 6., 6.]], AnnealConfig = AnnealConfig class defaults."""
    def __init__(self, hamiltonian = [[0, 6., 0, 0], [6., 2., 6., 0], [0, 6., 4., 6.], [0, 0, 6., 6.]],
        annealConfig = None):
        self.m_hamiltonian = np.array(hamiltonian)
        self.m_AnnealConfig = annealConfig if annealConfig else AnnealConfig()


    def ReadHamiltonian(self, file):
        self.m_hamiltonian = np.loadtxt(file, dtype=float)

    def ReadConfig(self, file): #wrapper, just calls corresponding method in AnnealConfig class
        self.m_AnnealConfig.ReadConfig(file)

    def GetAnnealMethod(self):
        return self.m_AnnealConfig.m_method


    def SetAnnealSimulated(self):
        self.m_AnnealConfig.m_method = "simulated"


    def SetAnnealQPU(self):
        self.m_AnnealConfig.m_method = "QPU"


    def ConstructBQMs(self):
        """Construct BQMs from provided hamiltonian and AnnealConfig paramaters.
        Save in the /BQM/ directory in both .pkl and human readable formats."""
        for group in self.m_AnnealConfig.m_groups:
            BQMPath = self.m_AnnealConfig.m_output + f'BQM/BQM_group{group}'
            valsPath = self.m_AnnealConfig.m_output + f'BQM/BQM_group{group}_vals'
            buildBQM(group, self.m_hamiltonian, self.m_AnnealConfig.m_lambdaVal, BQMPath, valsPath)


    def Anneal(self):
        for group in self.m_AnnealConfig.m_groups:
            BQMPath = self.m_AnnealConfig.m_output + f'BQM/BQM_group{group}'
            valsPath = self.m_AnnealConfig.m_output + f'BQM/BQM_group{group}_vals'
            for run in range(1, self.m_AnnealConfig.m_runs+1):
                outpath = self.m_AnnealConfig.m_output + f'{self.m_AnnealConfig.m_method}/group{group}/attempt{run}'
                AnnealBQM(BQMPath, valsPath, self.m_AnnealConfig.m_method, self.m_AnnealConfig.m_schedule, outpath)


    def Plot(self, plots = ["acc", "time", "size"]):
        """Generate plots form annealing results. Running both both simulated and QPU first for same trial number and number of runs. Save in /Plots/ directory.
            plots size vs accuracy, time, and number of qubits in three seperate plots. Specify which to plot by providing plots = ["acc", "time", "size"]. Default is all three"""
        datapath = self.m_AnnealConfig.m_output #anneal method, group, and run are added to the path in plot.py
        outpath = self.m_AnnealConfig.m_output + f'Plots/'
        GeneratePlots(self.m_AnnealConfig.m_trial, self.m_AnnealConfig.m_groups, self.m_AnnealConfig.m_runs, datapath, outpath, plots)



if __name__ == "__main__":
    # A = AnnealHandler()
    # A.ReadConfig("/workspaces/QuantumCognition/QuantumWalk/RefactoredConvAndAnneal/config.txt")
    # A.ReadHamiltonian("/workspaces/QuantumCognition/QuantumWalk/RefactoredConvAndAnneal/hamiltonian.txt")
    # A.ConstructBQMs()
    # A.Anneal()
    # A.SetAnnealQPU()
    # A.Anneal()
    # A.Plot()


<<<<<<< HEAD
    # A = AnnealHandler()
    # for i in [4, 8, 16, 32, 64]:
    #     A.ReadConfig("/workspaces/QuantumCognition/QuantumWalk/RefactoredConvAndAnneal/config.txt")
    #     A.m_AnnealConfig.m_trial = i
    #     A.m_AnnealConfig.m_output = f'Results/Trial_{i}/'
    #     A.ReadHamiltonian(f"/workspaces/QuantumCognition/QuantumWalk/ScalingTests/hamiltonians/hamiltonian_{i}.txt")
    #     A.ConstructBQMs()
    #     A.Anneal()
    #     A.SetAnnealQPU()
    #     A.Anneal()
    #     # A.Plot()
=======
    A = AnnealHandler()
    for i in [4, 8, 16, 32, 64]:
        A.ReadConfig("/workspaces/QuantumCognition/QuantumWalk/RefactoredConvAndAnneal/config.txt")
        A.m_AnnealConfig.m_trial = i
        A.m_AnnealConfig.m_output = f'Results/Trial_{i}/'
        A.ReadHamiltonian(f"/workspaces/QuantumCognition/QuantumWalk/ScalingTests/hamiltonians/hamiltonian_{i}.txt")
        A.ConstructBQMs()
        A.Anneal()
        A.SetAnnealQPU()
        A.Anneal()
        # A.Plot()
>>>>>>> refs/remotes/origin/main
