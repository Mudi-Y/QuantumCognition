import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os

path = "/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Plotting/Plots/" #where plots go

#read in data
simulatedDF = pd.read_pickle("/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Plotting/simulated.pkl")
quantumDF = pd.read_pickle("/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Plotting/quantum.pkl")

#Data to display
sa_accuracy = simulatedDF['acc_avg'].tolist()
qa_accuracy = quantumDF['acc_avg'].tolist()
sa_time = simulatedDF['time_avg'].tolist()
qa_time = quantumDF['time_avg'].tolist()

sa_acc_std = simulatedDF['acc_std'].tolist()
qa_acc_std = quantumDF['acc_std'].tolist()
sa_time_std = simulatedDF['time_std'].tolist()
qa_time_std = quantumDF['time_std'].tolist()

#create plot output directory
os.makedirs(os.path.dirname(path), exist_ok=True)


#font specifications for plots
plt.rcParams['font.size'] = 12
plt.rcParams["font.family"] = "Serif"

#plot two seperate plots

#Params from experiments:
sizes = [4, 8, 16, 32, 64]
groups = [5, 6, 7]
lambdas = [-2, -3, -4]

x = sizes

# basic plots
# accuracy of anneals
fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(left=0.1, right=0.95, top=0.97, bottom=0.12)
for i in range(len(sizes), 82, len(sizes)):
    ax.errorbar(x, sa_accuracy[i-len(sizes) : i], yerr = sa_acc_std[i-len(sizes) : i], ls='none', capsize=6, color='k') 
    ax.errorbar(x, qa_accuracy[i-len(sizes) : i], yerr = qa_acc_std[i-len(sizes) : i], ls='none', capsize=6, color='k')
    ax.plot(x, sa_accuracy[i-len(sizes) : i], '-v', color='#f8a652', label="SA" if i == len(sizes) else "")
    ax.plot(x, qa_accuracy[i-len(sizes) : i], '-o', color='#5f7d41',label="QA-P" if i == len(sizes) else "")
major_ticks = x
ymajor_ticks = [0,1,2,3,4]
yminor_ticks = [0,0.5,1,1.5,2,2.5,3.5,4]
ax.set_xscale('log', base=2)
ax.set_xticks(major_ticks)
ax.set_yticks(ymajor_ticks)
ax.set_yticks(yminor_ticks, minor = True)
ax.set_xlabel("Hamiltonian Size")
ax.set_ylabel("Absolute Error")
ax.legend(borderpad=.1, loc='upper left', labelspacing=.15)
ax.grid(visible=True, which='major', color='gray', linestyle='-')
plt.grid(visible=True, which='minor', color='gray', linestyle='--')
# plt.savefig(path+'Accuracy.png')
# plt.savefig(path+'annealing_accuracy.pdf')
plt.close(fig)


#time taken for each anneal
fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(left=0.12, right=0.95, top=0.97, bottom=0.12)
for i in range(len(sizes), 82, len(sizes)):
    ax.errorbar(x, sa_time[i-len(sizes) : i], yerr = sa_time_std[i-len(sizes) : i], ls='none', capsize=6, color='k') 
    ax.errorbar(x, qa_time[i-len(sizes) : i], yerr = qa_time_std[i-len(sizes) : i], ls='none', capsize=6, color='k')
    ax.plot(x, sa_time[i-len(sizes) : i], '-v', color='#f8a652', label="SA" if i == len(sizes) else "")
    ax.plot(x, qa_time[i-len(sizes) : i], '-o', color='#5f7d41',label="QA-P" if i == len(sizes) else "")
major_ticks = x
ymajor_ticks = [0,10,20,30]
yminor_ticks = [0,5,10,15,20,25,30]
ax.set_xscale('log', base=2)
ax.set_xticks(major_ticks)
ax.set_yticks(ymajor_ticks)
ax.set_yticks(yminor_ticks, minor = True)
ax.set_xlabel("Hamiltonian Size") 
ax.set_ylabel("Time (ms)")
ax.legend(borderpad=.1, loc='upper left', labelspacing=.15)
ax.grid(visible=True, which='major', color='gray', linestyle='-')
ax.grid(visible=True, which='minor', color='gray', linestyle='--')
# plt.savefig(path+'Timing.png')
# plt.savefig(path+'annealing_timing.pdf')
plt.close(fig)


#find setting that has lowest value of 'param' for 'dataframe' and hamiltonian 'size' fitting the 'condition':
#dataframe['param2'] <= 'value'
def findBestConfig(dataframe, size = 4, param="acc_avg", param2 = "acc_std", value = float('inf')):
    size_string = f"hSize_{size}"
    rows = dataframe.loc[(dataframe['hSize'] == size_string) & 
                         (dataframe[param2] <= value)]
    rows = rows.reset_index(drop=True)
    
    if not rows.empty:
        best = rows.iloc[0]
        for i, row in rows.iterrows():
            if row[param] < best[param]:
                best = row
        return best
    return "empty dataframe"


#average time taken to reach absolute error e
e = 1
#get relevant data for SA and QA for each size hamiltonian
saData = []
saStd = []
qaData = []
qaStd = []
for size in sizes:
    saResult = findBestConfig(simulatedDF, size=size, param="acc_avg", param2="acc_avg", value=e)
    qaResult = findBestConfig(quantumDF, size=size, param="acc_avg", param2="acc_avg", value=e)
    saData.append(saResult["time_avg"])
    saStd.append(saResult["time_std"])
    qaData.append(qaResult["time_avg"])
    qaStd.append(qaResult["time_std"])
#plot time scaling
fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(1, 1, 1)
ax.errorbar(x, saData, yerr = saStd, ls='none', capsize=6, color='k') 
ax.errorbar(x, qaData, yerr = qaStd, ls='none', capsize=6, color='k')
ax.plot(x, saData, '-v', color='#f8a652', label="SA")
ax.plot(x, qaData, '-o', color='#5f7d41',label="QA-P")
major_ticks = x
ymajor_ticks = [0,10,20,30]
yminor_ticks = [0,5,10,15,20,25,30]
ax.set_xscale('log', base=2)
ax.set_xticks(major_ticks)
ax.set_yticks(ymajor_ticks)
ax.set_yticks(yminor_ticks, minor = True)
ax.set_xlabel("Hamiltonian Size") 
ax.set_ylabel("Time (ms)")
ax.legend(borderpad=.1, loc='upper left', labelspacing=.15)
ax.grid(visible=True, which='major', color='gray', linestyle='-')
ax.grid(visible=True, which='minor', color='gray', linestyle='--')
# plt.savefig(path+'TimeScaling.png')
# plt.savefig(path+'time_scaling.pdf')
plt.close(fig)



#average absolute error after time t
t = 30
#get relevant data for SA and QA for each size hamiltonian
saData = []
saStd = []
qaData = []
qaStd = []
for size in sizes:
    saResult = findBestConfig(simulatedDF, size=size, param="acc_avg", param2="time_avg", value=t)
    qaResult = findBestConfig(quantumDF, size=size, param="acc_avg", param2="time_avg", value=t)
    saData.append(saResult["acc_avg"])
    saStd.append(saResult["acc_std"])
    qaData.append(qaResult["acc_avg"])
    qaStd.append(qaResult["acc_std"])
#plot accuracy scaling
fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(left=0.12, right=0.95, top=0.97, bottom=0.12)
ax.errorbar(x, saData, yerr = saStd, ls='none', capsize=6, color='k') 
ax.errorbar(x, qaData, yerr = qaStd, ls='none', capsize=6, color='k')
ax.plot(x, saData, '-v', color='#f8a652', label="SA")
ax.plot(x, qaData, '-o', color='#5f7d41',label="QA-P")
major_ticks = x
ymajor_ticks = [0,0.2, 0.4, 0.6, 0.8, 1.0]
yminor_ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
ax.set_xscale('log', base=2)
ax.set_xticks(major_ticks)
ax.set_yticks(ymajor_ticks)
ax.set_yticks(yminor_ticks, minor = True)
ax.set_xlabel("Hamiltonian Size")
ax.set_ylabel("Absolute Error")
ax.legend(borderpad=.1, loc='upper left', labelspacing=.15)
ax.grid(visible=True, which='major', color='gray', linestyle='-')
plt.grid(visible=True, which='minor', color='gray', linestyle='--')
# plt.savefig(path+'AccuracyScaling.png')
# plt.savefig(path+'accuracy_scaling.pdf')
plt.close(fig)


# Scatter plots of accuracy vs time for each hamiltonian size
for size in sizes:
    fig = plt.figure(figsize = (5, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.9, bottom=0.12)
    xdata = simulatedDF.loc[simulatedDF["hSize"] == f"hSize_{size}"]["time_avg"]
    ydata = simulatedDF.loc[simulatedDF["hSize"] == f"hSize_{size}"]["acc_avg"]
    sdata = simulatedDF.loc[simulatedDF["hSize"] == f"hSize_{size}"]["n_qubits"]
    ax.scatter(xdata, ydata, s=sdata)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Absolute Error")
    ax.set_title(f"Simulated Annealing Hamiltonian Size {size}")
    ax.grid(visible=True, which='major', color='gray', linestyle='-')
    plt.grid(visible=True, which='minor', color='gray', linestyle='--')
    plt.savefig(path+'acc_vs_time/SA/'+f'SA_accuracy_vs_time_{size}.png')
    plt.savefig(path+'acc_vs_time/SA/'+f'SA_accuracy_vs_time_{size}.pdf')
    plt.close(fig)

for size in sizes:
    fig = plt.figure(figsize = (5, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.9, bottom=0.12)
    xdata = quantumDF.loc[quantumDF["hSize"] == f"hSize_{size}"]["time_avg"]
    ydata = quantumDF.loc[quantumDF["hSize"] == f"hSize_{size}"]["acc_avg"]
    sdata = quantumDF.loc[quantumDF["hSize"] == f"hSize_{size}"]["n_qubits"]
    ax.scatter(xdata, ydata, s=sdata)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Absolute Error")
    ax.set_title(f"Quantum Annealing Hamiltonian Size {size}")
    ax.grid(visible=True, which='major', color='gray', linestyle='-')
    plt.grid(visible=True, which='minor', color='gray', linestyle='--')
    plt.savefig(path+'acc_vs_time/QA/'+f'QA_accuracy_vs_time_{size}.png')
    plt.savefig(path+'acc_vs_time/QA/'+f'QA_accuracy_vs_time_{size}.pdf')
    plt.close(fig)



# #draw embedding graphs
# #for each BQM, embed into the solver and draw the topology graph
# from dwave.system.samplers import DWaveSampler
# from dwave.system.composites import EmbeddingComposite
# import pickle
# import dwave_networkx as dnx
# import matplotlib.pyplot as plt
# import networkx as nx
# from pathlib import Path

# sampler = DWaveSampler(token="DEV-bc8e40b9fdc4a711ebf6cdcf71747114f89d8baf")
# composite = EmbeddingComposite(sampler)
# print(sampler.properties["topology"])

# for size in [4, 8, 16, 32, 64]:
#     for group in [5,6,7]:
#         fig = plt.figure(figsize = (100, 100))
#         bqm_path =  Path(__file__).parent / f"../Outputs/HSize-{size}/Group-{group}/Lambda--2/Schedule-[(0.0, 0.0), (10.0, 0.4), (15.0, 0.4), (25.0, 0.8), (30.0, 0.8), (40.0, 1.0)]/Results/Trial_1/BQM/BQM_group{group}"
#         file = open(bqm_path, 'rb')
#         bqm = pickle.load(file)
#         file.close()

#         res = composite.sample(bqm, return_embedding=True)
#         embedding = res.info["embedding_context"]["embedding"]

#         G = dnx.pegasus_graph(16) #16 because the topology of the solver s "pegasus" with size "16"
#         dnx.draw_pegasus_embedding(G, embedding, show_labels=True)
#         plt.savefig(path+'embeddings/'+f'hSize_{size}_group{group}.png')
#         plt.savefig(path+'embeddings/'+f'hSize_{size}_group{group}.pdf')
#         plt.close(fig)