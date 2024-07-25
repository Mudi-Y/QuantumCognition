import sys, os

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../RefactoredConvAndAnneal/")
from RefactoredConvAndAnneal import main as RCA


#path for results
outPath = "/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/ImportantResults/"

# #find appropriate hyperparamaters for each size of hamiltonian
# #paramas are group sizes, anneal schedule, and lambda
# hSizes = [4, 8, 16, 32, 64]
# groups = [2, 3, 4, 5, 6, 7, 8, 9]
# lambdas = [0, -1,-2,-3, -4, -5, -6]

# annealSchedules = []
# times = [5.0, 10.0, 15.0, 20.0]
# for time in times:
#     schedule = [(0.0, 0.0), (time, 0.4), (time + 5.0, 0.4), (time + 5.0 + time, 0.8), (time + 5.0 + time + 5.0, 0.8), (time + 5.0 + time + 5.0 + time, 1.0)]
#     annealSchedules.append(schedule)


#find appropriate hyperparamaters for each size of hamiltonian
#paramas are group sizes, anneal schedule, and lambda
hSizes = [16]
groups = [5, 6, 7]
lambdas = [-2,-3, -4]

annealSchedules = []
times = [5.0, 10.0, 15.0]
for time in times:
    schedule = [(0.0, 0.0), (time, 0.4), (time + 5.0, 0.4), (time + 5.0 + time, 0.8), (time + 5.0 + time + 5.0, 0.8), (time + 5.0 + time + 5.0 + time, 1.0)]
    annealSchedules.append(schedule)

#Create an anneal handler for each experiment. Save annealing results
# #SA
# counter = 1
# for size in hSizes:
#     for group in groups:
#         for lamb in lambdas:
#             for schedule in annealSchedules:
#                 #create anneal configuration
#                 config = RCA.AnnealConfig(
#                     method="simulated",
#                     schedule=schedule,
#                     groups=[group],
#                     num_runs=3,
#                     trial=1,
#                     lambdaVal= lamb,
#                     output=outPath + f"HSize-{size}/Group-{group}/Lambda-{lamb}/Schedule-{schedule}"
#                 )
                
#                 #create anneal handler
#                 handler = RCA.AnnealHandler(annealConfig=config)
#                 handler.ReadHamiltonian(f"/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Hamiltonians/hamiltonian_{size}.txt")
#                 handler.ConstructBQMs()
#                 handler.Anneal()

#                 print("H Size: ", size)
#                 print("Iter: ", counter)
#                 counter += 1



#QA
counter = 1
for size in hSizes:
    for group in groups:
        for lamb in lambdas:
            for schedule in annealSchedules:
                #create anneal configuration
                config = RCA.AnnealConfig(
                    method="QPU",
                    schedule=schedule,
                    groups=[group],
                    num_runs=3,
                    trial=1,
                    lambdaVal= lamb,
                    output=outPath + f"HSize-{size}/Group-{group}/Lambda-{lamb}/Schedule-{schedule}"
                )
                
                #create anneal handler
                handler = RCA.AnnealHandler(annealConfig=config)
                handler.ReadHamiltonian(f"/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Hamiltonians/hamiltonian_{size}.txt")
                handler.ConstructBQMs()
                handler.Anneal()

                print("H Size: ", size)
                print("Iter: ", counter)
                counter += 1