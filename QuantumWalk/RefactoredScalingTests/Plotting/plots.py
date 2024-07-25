import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os
    
def extract_value(filename, p):
    with open(filename, mode="rt", encoding="utf-8") as docFile:
        doc = docFile.read()
        val = re.findall(p, doc)
        if "{" in val[0]:
            val = float(re.findall(r"[-+]?(?:\d*\.*\d+)", val[0])[0])
        else:
            val = float(val[0].split("  ")[1])
    return val


#Extract Data from Logs
#setup code for extraction. IMPORTANT!
path = "/workspaces/QuantumCognition/QuantumWalk/ScalingTests/Plotting/Plots/"
numbers = [4, 8, 16, 32, 64]
annealMethods = ["simulated", "QPU"]
group = 6
attempts = 3

#make dataframes to store data
SA = pd.DataFrame(index = numbers, columns = ["acc", "time", "acc_avg", "time_avg", "acc_std", "time_std"])
QA = pd.DataFrame(index = numbers, columns = ["acc", "time", "acc_avg", "time_avg", "acc_std", "time_std"])

#make relevant fields arrays
for number in numbers:
    for data in ["acc", "time"]:
        SA[data][number] = []
        QA[data][number] = []


datapath = "/workspaces/QuantumCognition/QuantumWalk/ScalingTests/Results/"

#extraction script
for annealMethod in annealMethods:
        for number in numbers:
            for i in range(attempts):
                fname = datapath + f"Trial_{number}/{annealMethod}/group{group}/attempt{i+1}"
                if annealMethod == "simulated":
                    #extract and append data for simulated
                    p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                    SA["acc"][number].append(extract_value(fname, p))

                    p = re.compile("Time:  [-+]?(?:\d*\.*\d+)")
                    SA["time"][number].append(extract_value(fname, p)*1000) #convert from seconds to miliseconds
                else:
                    #extract and append data for QPU
                    p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                    QA["acc"][number].append(extract_value(fname, p))

                    p = re.compile(".*qpu_sampling_time.* [-+]?(?:\d*\.*\d+)")
                    QA["time"][number].append(extract_value(fname, p)*0.001) #convert from microseconds to miliseconds



#Fill in Avg, STDEV
for number in numbers:
    for data in ["acc_avg", "time_avg"]:
        SA[data][number] = np.mean(SA[data[:-4]][number])
        QA[data][number] = np.mean(QA[data[:-4]][number])
for number in numbers:
    for data in ["acc_std", "time_std"]:
        SA[data][number] = np.std(SA[data[:-4]][number])
        QA[data][number] = np.std(QA[data[:-4]][number])

#font specifications
plt.rcParams['font.size'] = 12
plt.rcParams["font.family"] = "Serif"


#X values are hamiltonian sizes
x = numbers


#Data to display
sa_accuracy = SA['acc_avg'].tolist()
qa_accuracy = QA['acc_avg'].tolist()
sa_time = SA['time_avg'].tolist()
qa_time = QA['time_avg'].tolist()

sa_acc_std = SA['acc_std'].tolist()
qa_acc_std = QA['acc_std'].tolist()
sa_time_std = SA['time_std'].tolist()
qa_time_std = QA['time_std'].tolist()


#create plot output directory
os.makedirs(os.path.dirname(path), exist_ok=True)

#plot two seperate plots

fig = plt.figure(figsize = (2.5, 2.5))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(left=0.225, right=0.959, top=0.97, bottom=0.2)
ax.errorbar(x, sa_accuracy, yerr = sa_acc_std, ls='none', capsize=6, color='k') 
ax.plot(x, sa_accuracy, '-v', color='#f8a652',label="SA")
ax.errorbar(x, qa_accuracy, yerr = sa_acc_std, ls='none', capsize=6, color='k') 
ax.plot(x, qa_accuracy, '-o', color='#5f7d41',label="QA-P")
major_ticks = numbers
ymajor_ticks = [0,1,2,3,4]
yminor_ticks = [0,0.5,1,1.5,2,2.5,3.5,4]
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


fig = plt.figure(figsize = (2.5, 2.5))
ax = fig.add_subplot(1, 1, 1)
plt.subplots_adjust(left=0.225, right=0.959, top=0.975, bottom=0.2)
ax.errorbar(x, sa_time, yerr = sa_time_std, ls='none', capsize=6, color='k') 
ax.plot(x, sa_time, '-v', color='#f8a652', label="SA")
ax.errorbar(x, qa_time, yerr = qa_time_std, ls='none', capsize=6, color='k') 
ax.plot(x, qa_time, '-o', color='#5f7d41', label="QA-P")
major_ticks = numbers
ymajor_ticks = [0,10,20,30]
yminor_ticks = [0,5,10,15,20,25,30]
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