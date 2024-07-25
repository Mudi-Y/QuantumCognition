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

#x axis is hamiltonian sizes
x = [4, 8, 16]

fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(left=0.1, right=0.95, top=0.97, bottom=0.12)
for i in range(3, 82, 3):
    ax.errorbar(x, sa_accuracy[i-3 : i], yerr = sa_acc_std[i-3 : i], ls='none', capsize=6, color='k') 
    ax.errorbar(x, qa_accuracy[i-3 : i], yerr = sa_acc_std[i-3 : i], ls='none', capsize=6, color='k')
    ax.plot(x, sa_accuracy[i-3 : i], '-v', color='#f8a652', label="SA" if i == 3 else "")
    ax.plot(x, qa_accuracy[i-3 : i], '-o', color='#5f7d41',label="QA-P" if i == 3 else "")
major_ticks = x
ymajor_ticks = [0,1,2,3,4]
yminor_ticks = [0,0.5,1,1.5,2,2.5,3.5,4]
ax.set_xscale('log', base=2)
ax.set_xticks(major_ticks)
ax.set_yticks(ymajor_ticks)
ax.set_yticks(yminor_ticks, minor = True)
ax.set_xlabel("Hamiltonian Size")
ax.set_ylabel("Absolute error")
ax.legend(borderpad=.1, loc='upper left', labelspacing=.15)
ax.grid(visible=True, which='major', color='gray', linestyle='-')
plt.grid(visible=True, which='minor', color='gray', linestyle='--')
plt.savefig(path+'Accuracy.png')
plt.savefig(path+'annealing_accuracy.pdf')


fig = plt.figure(figsize = (5, 5))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(left=0.12, right=0.95, top=0.97, bottom=0.12)
for i in range(3, 82, 3):
    ax.errorbar(x, sa_time[i-3 : i], yerr = sa_time_std[i-3 : i], ls='none', capsize=6, color='k') 
    ax.errorbar(x, qa_time[i-3 : i], yerr = sa_time_std[i-3 : i], ls='none', capsize=6, color='k')
    ax.plot(x, sa_time[i-3 : i], '-v', color='#f8a652', label="SA" if i == 3 else "")
    ax.plot(x, qa_time[i-3 : i], '-o', color='#5f7d41',label="QA-P" if i == 3 else "")
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
plt.savefig(path+'Timing.png')
plt.savefig(path+'annealing_timing.pdf')