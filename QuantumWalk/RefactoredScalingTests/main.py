import sys, os

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../RefactoredConvAndAnneal/")
from RefactoredConvAndAnneal import main as RCA


#path for results
outPath = "/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Outputs/"

#Trial number. VERY IMPORTANT! Increase each time to prevent over-writing old results
trial = 3

# #find appropriate hyperparamaters for each size of hamiltonian
# #paramas are group sizes, anneal schedule, and lambda
# hSizes = [4, 8, 16, 32, 64]
# groups = [5, 6, 7]
# lambdas = [-2,-3, -4]

# annealSchedules = []
# times = [5.0, 10.0, 15.0]
# for time in times:
#     schedule = [(0.0, 0.0), (time, 0.4), (time + 5.0, 0.4), (time + 5.0 + time, 0.8), (time + 5.0 + time + 5.0, 0.8), (time + 5.0 + time + 5.0 + time, 1.0)]
#     annealSchedules.append(schedule)


#find appropriate hyperparamaters for each size of hamiltonian
#paramas are group sizes, anneal schedule, and lambda
hSizes = [8]
groups = [5, 6, 7]
lambdas = [-15, -20, -25]
annealSchedules = []
times = [5.0, 10.0, 15.0]
for time in times:
    schedule = [(0.0, 0.0), (time, 0.4), (time + 5.0, 0.4), (time + 5.0 + time, 0.8), (time + 5.0 + time + 5.0, 0.8), (time + 5.0 + time + 5.0 + time, 1.0)]
    annealSchedules.append(schedule)

# # Create an anneal handler for each experiment. Save annealing results
# #SA
# for size in hSizes:
#     counter = 1
#     print("########## H Size: ", size, " ##########")
#     for group in groups:
#         for lamb in lambdas:
#             for schedule in annealSchedules:
#                 #create anneal configuration
#                 config = RCA.AnnealConfig(
#                     method="simulated",
#                     schedule=schedule,
#                     groups=[group],
#                     num_runs=3,
#                     trial=trial,
#                     lambdaVal= lamb,
#                     output=outPath + f"HSize-{size}/Group-{group}/Lambda-{lamb}/Schedule-{schedule}"
#                 )
                
#                 #create anneal handler
#                 handler = RCA.AnnealHandler(annealConfig=config)
#                 handler.ReadHamiltonian(f"/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Hamiltonians/hamiltonian_{size}.txt")
#                 handler.ConstructBQMs()
#                 handler.Anneal()

#                 print("Iter: ", counter, "/", (len(groups) * len(lambdas) * len(annealSchedules)))
#                 counter += 1



# QA
for size in hSizes:
    counter = 1
    print("########## H Size: ", size, " ##########")
    for group in groups:
        for lamb in lambdas:
            for schedule in annealSchedules:
                #create anneal configuration
                config = RCA.AnnealConfig(
                    method="QPU",
                    schedule=schedule,
                    groups=[group],
                    num_runs=3,
                    trial=trial,
                    lambdaVal= lamb,
                    output=outPath + f"HSize-{size}/Group-{group}/Lambda-{lamb}/Schedule-{schedule}"
                )
                
                #create anneal handler
                handler = RCA.AnnealHandler(annealConfig=config)
                handler.ReadHamiltonian(f"/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Hamiltonians/hamiltonian_{size}.txt")
                handler.ConstructBQMs()
                handler.Anneal()

                print("Iter: ", counter, "/", (len(groups) * len(lambdas) * len(annealSchedules)))
                counter += 1