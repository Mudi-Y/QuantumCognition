import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re


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
number = "4"
path = "/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/"
annealMethods = ["simulated", "QPU"]
groups=[2,3,4,5,6,7,8,9]
attempts = 3

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
                fname = path+f"Results/{annealMethod}{number}/{annealMethod}_group{group}_attempt{i+1}"
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


# #Ploting 3x1 plots (3 x 1" x 9")
# fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize = (3, 6))
# fig.suptitle('Vertically stacked subplots')

# ax1.scatter(x, sa_accuracy, color='#f8a652')
# ax1.scatter(x, qa_accuracy, color='#5f7d41')
# ax1.set_title("Accuracy")
# ax1.set(ylabel="Absolute Error")

# ax2.scatter(x, sa_time, color='#f8a652')
# ax2.scatter(x, qa_time, color='#5f7d41')
# ax2.set_title("Timing")
# ax2.set(ylabel="Micro Seconds")

# ax3.scatter(x, qa_size, color='r')
# ax3.set_title("Size")
# ax3.set(xlabel="Number of Groups", ylabel="Number of Qubits")

# fig.tight_layout()
# fig.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/vertical')
# fig.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/vertical.pdf')


#plot three seperate plots
plt.figure()
plt.errorbar(x, sa_accuracy, yerr = sa_acc_std, ls='none', capsize=6, color='k') 
plt.plot(x, sa_accuracy, '-o', color='#f8a652',label="SA")
plt.errorbar(x, qa_accuracy, yerr = sa_acc_std, ls='none', capsize=6, color='k') 
plt.plot(x, qa_accuracy, '-o', color='#5f7d41',label="QA-P")
plt.xlabel("Number of Groups (r)")
plt.ylabel("Absolute Error")
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/Accuracy')
plt.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/Accuracy.pdf')


plt.figure()
plt.errorbar(x, sa_time, yerr = sa_time_std, ls='none', capsize=6, color='k') 
plt.plot(x, sa_time, '-o', color='#f8a652', label="SA")
plt.errorbar(x, qa_time, yerr = qa_time_std, ls='none', capsize=6, color='k') 
plt.plot(x, qa_time, '-o', color='#5f7d41', label="QA-P")
plt.xlabel("Number of Groups (r)") 
plt.ylabel("Milliseconds")
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/Timing')
plt.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/Timing.pdf')


plt.figure()
plt.plot(x, qa_size, '-o', color='k')
plt.xlabel("Number of Groups (r)")
plt.ylabel("Number of Qubits")
plt.tight_layout()
plt.grid()
plt.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/Size')
plt.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/Size.pdf')

