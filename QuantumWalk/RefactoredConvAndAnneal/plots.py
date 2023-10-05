import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os


def GeneratePlots(trial, groups, runs, datapath, outpath, plots):
    #extraction function
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
    path = outpath
    number = trial
    annealMethods = ["simulated", "QPU"]
    groups = groups
    attempts = runs

    #make dataframes to store data
    SA = pd.DataFrame(index = groups, columns = ["acc", "time", "size", "acc_avg", "time_avg", "size_avg", "acc_std", "time_std", "size_std"])
    QA = pd.DataFrame(index = groups, columns = ["acc", "time", "size", "acc_avg", "time_avg", "size_avg", "acc_std", "time_std", "size_std"])

    #make relevant fields arrays
    for group in groups:
        for data in ["acc", "time", "size"]:
            SA[data][group] = []
            QA[data][group] = []


    #extraction script
    for annealMethod in annealMethods:
            for group in groups:
                for i in range(attempts):
                    fname = datapath + f"{annealMethod}/group{group}/attempt{i+1}"
                    if annealMethod == "simulated":
                        #extract and append data for simulated
                        p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                        SA["acc"][group].append(extract_value(fname, p))

                        p = re.compile("Number of Qubits:  [-+]?(?:\d*\.*\d+)")
                        SA["size"][group].append(extract_value(fname, p))

                        p = re.compile("Time:  [-+]?(?:\d*\.*\d+)")
                        SA["time"][group].append(extract_value(fname, p)*1000) #convert from seconds to miliseconds
                    else:
                        #extract and append data for QPU
                        p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                        QA["acc"][group].append(extract_value(fname, p))

                        p = re.compile("Number of Qubits:  [-+]?(?:\d*\.*\d+)")
                        QA["size"][group].append(extract_value(fname, p))

                        p = re.compile(".*qpu_sampling_time.* [-+]?(?:\d*\.*\d+)")
                        QA["time"][group].append(extract_value(fname, p)*0.001) #convert from microseconds to miliseconds

    #Fill in Avg, STDEV
    for group in groups:
        for data in ["acc_avg", "time_avg", "size_avg"]:
            SA[data][group] = np.mean(SA[data[:-4]][group])
            QA[data][group] = np.mean(QA[data[:-4]][group])
    for group in groups:
        for data in ["acc_std", "time_std", "size_std"]:
            SA[data][group] = np.std(SA[data[:-4]][group])
            QA[data][group] = np.std(QA[data[:-4]][group])

    #font specifications
    plt.rcParams['font.size'] = 12
    plt.rcParams["font.family"] = "Serif"


    #X values are group sizes
    x = groups

    #Dummy data to display
    # sa_accuracy = np.sin(x ** 2)
    # qa_accuracy = np.cos(x ** 2)
    # sa_time = -np.sin(x ** 2)
    # qa_time = -np.cos(x ** 2)
    # sa_size = np.sin(x ** 2)
    # qa_size = np.cos(x ** 2)

    #Data to display
    sa_accuracy = SA['acc_avg'].tolist()
    qa_accuracy = QA['acc_avg'].tolist()
    sa_time = SA['time_avg'].tolist()
    qa_time = QA['time_avg'].tolist()
    sa_size = SA['size_avg'].tolist()
    qa_size = QA['size_avg'].tolist()

    sa_acc_std = SA['acc_std'].tolist()
    qa_acc_std = QA['acc_std'].tolist()
    sa_time_std = SA['time_std'].tolist()
    qa_time_std = QA['time_std'].tolist()


    #create plot output directory
    os.makedirs(os.path.dirname(path), exist_ok=True)

    #plot three seperate plots
    if "acc" in plots:
        fig = plt.figure(figsize = (2.5, 2.5))
        ax = fig.add_subplot(1, 1, 1)
        plt.subplots_adjust(left=0.225, right=0.959, top=0.97, bottom=0.2)
        ax.errorbar(x, sa_accuracy, yerr = sa_acc_std, ls='none', capsize=6, color='k') 
        ax.plot(x, sa_accuracy, '-v', color='#f8a652',label="SA")
        ax.errorbar(x, qa_accuracy, yerr = sa_acc_std, ls='none', capsize=6, color='k') 
        ax.plot(x, qa_accuracy, '-o', color='#5f7d41',label="QA-P")
        major_ticks = [2,4,6,8,10]
        minor_ticks = [2,3,4,5,6,7,8,9,10]
        ymajor_ticks = [0,1,2,3,4]
        yminor_ticks = [0,0.5,1,1.5,2,2.5,3.5,4]
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(ymajor_ticks)
        ax.set_yticks(yminor_ticks, minor = True)
        ax.set_xlabel("Number of groups (r)")
        ax.set_ylabel("Absolute error")
        ax.legend(borderpad=.1, loc='upper left', labelspacing=.15)
        ax.grid(visible=True, which='major', color='gray', linestyle='-')
        plt.grid(visible=True, which='minor', color='gray', linestyle='--')
        plt.savefig(path+'Accuracy.png')
        plt.savefig(path+'annealing_accuracy.pdf')

    if "time" in plots:
        fig = plt.figure(figsize = (2.5, 2.5))
        ax = fig.add_subplot(1, 1, 1)
        plt.subplots_adjust(left=0.225, right=0.959, top=0.975, bottom=0.2)
        ax.errorbar(x, sa_time, yerr = sa_time_std, ls='none', capsize=6, color='k') 
        ax.plot(x, sa_time, '-v', color='#f8a652', label="SA")
        ax.errorbar(x, qa_time, yerr = qa_time_std, ls='none', capsize=6, color='k') 
        ax.plot(x, qa_time, '-o', color='#5f7d41', label="QA-P")
        major_ticks = [2,4,6,8,10]
        minor_ticks = [2,3,4,5,6,7,8,9,10]
        ymajor_ticks = [0,10,20,30]
        yminor_ticks = [0,5,10,15,20,25,30]
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(ymajor_ticks)
        ax.set_yticks(yminor_ticks, minor = True)
        ax.set_xlabel("Number of groups (r)") 
        ax.set_ylabel("Time (ms)")
        ax.legend(borderpad=.1, loc='upper left', labelspacing=.15)
        ax.grid(visible=True, which='major', color='gray', linestyle='-')
        ax.grid(visible=True, which='minor', color='gray', linestyle='--')
        plt.savefig(path+'Timing.png')
        plt.savefig(path+'annealing_timing.pdf')

    if "size" in plots:
        fig = plt.figure(figsize = (2.5, 2.5))
        ax = fig.add_subplot(1, 1, 1)
        plt.subplots_adjust(left=0.225, right=0.959, top=0.975, bottom=0.2)
        ax.plot(x, qa_size, '-o', color='#5f7d41')
        major_ticks = [2,4,6,8,10]
        minor_ticks = [2,3,4,5,6,7,8,9,10]
        ymajor_ticks = [0,20,40,60,80]
        yminor_ticks = [0,10,20,30,40,50,60,70,80]
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(ymajor_ticks)
        ax.set_yticks(yminor_ticks, minor = True)
        ax.set_xlabel("Number of groups (r)")
        ax.set_ylabel("Number of qubits")
        ax.grid(visible=True, which='major', color='gray', linestyle='-')
        ax.grid(visible=True, which='minor', color='gray', linestyle='--')
        plt.savefig(path+'Size.png')
        plt.savefig(path+'annealing_size.pdf')
