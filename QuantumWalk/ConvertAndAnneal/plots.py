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
x = np.array([x for x in range(2,18,2)])
SA = pd.DataFrame(index = [2,4,6,8,10,12,14,16], columns = ["acc", "time", "size", "acc_avg", "time_avg", "size_avg"])
QA = pd.DataFrame(index = [2,4,6,8,10,12,14,16], columns = ["acc", "time", "size", "acc_avg", "time_avg", "size_avg"])

#make relevant fields arrays
for group in [2,4,6,8,10,12,14,16]:
    for data in ["acc", "time", "size"]:
        SA[data][group] = []
        QA[data][group] = []


number = "2"
path = "/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/"
annealMethods = ["simulated", "QPU"]
groups=[2,4,6,8,10,12,14,16]


for annealMethod in annealMethods:
        for group in groups:
            for i in range(10):
                fname = path+f"Results/{annealMethod}{number}/{annealMethod}_group{group}_attempt{i+1}"
                if annealMethod == "simulated":
                    #extract and append data for simulated
                    p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                    SA["acc"][group].append(extract_value(fname, p))

                    p = re.compile("Number of Qubits:  [-+]?(?:\d*\.*\d+)")
                    SA["size"][group].append(extract_value(fname, p))

                    p = re.compile("Time:  [-+]?(?:\d*\.*\d+)")
                    SA["time"][group].append(extract_value(fname, p)*1000000) #convert from seconds to microseconds
                else:
                    #extract and append data for QPU
                    p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                    QA["acc"][group].append(extract_value(fname, p))

                    p = re.compile("Number of Qubits:  [-+]?(?:\d*\.*\d+)")
                    QA["size"][group].append(extract_value(fname, p))

                    p = re.compile(".*qpu_sampling_time.* [-+]?(?:\d*\.*\d+)")
                    QA["time"][group].append(extract_value(fname, p))

#make relevant fields arrays
for group in [2,4,6,8,10,12,14,16]:
    for data in ["acc_avg", "time_avg", "size_avg"]:
        SA[data][group] = np.mean(SA[data[:-4]][group])
        QA[data][group] = np.mean(QA[data[:-4]][group])

#font specifications
plt.rcParams['font.size'] = 12
plt.rcParams["font.family"] = "Serif"


#X values are group sizes
x = np.array([x for x in range(2,18,2)])

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

#Ploting 3x1 plots (3 x 1" x 9")
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize = (3, 6))
fig.suptitle('Vertically stacked subplots')

ax1.scatter(x, sa_accuracy, color='#f8a652')
ax1.scatter(x, qa_accuracy, color='#5f7d41')
ax1.set_title("Accuracy")
ax1.set(ylabel="Absolute Error")
ax1.set_yscale('log')

ax2.scatter(x, sa_time, color='#f8a652')
ax2.scatter(x, qa_time, color='#5f7d41')
ax2.set_title("Timing")
ax2.set(ylabel="Micro Seconds")

ax3.scatter(x, qa_size, color='r')
ax3.set_title("Size")
ax3.set(xlabel="Number of Groups", ylabel="Number of Qubits")

fig.tight_layout()
fig.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/vertical')
fig.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/vertical.pdf')



#Plotting 1x3 plots (3 x 3" x 3")
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (9,3))
fig.suptitle('Horizontally stacked subplots')

ax1.scatter(x, sa_accuracy, color='#f8a652')
ax1.scatter(x, qa_accuracy, color='#5f7d41')
ax1.set_title("Accuracy")
ax1.set(xlabel="Number of Groups", ylabel="Absolute Error")
ax1.set_yscale('log')

ax2.scatter(x, sa_time, color='#f8a652')
ax2.scatter(x, qa_time, color='#5f7d41')
ax2.set_title("Timing")
ax2.set(xlabel="Number of Groups", ylabel="Micro Seconds")

ax3.scatter(x, qa_size, color='r')
ax3.set_title("Size")
ax3.set(xlabel="Number of Groups", ylabel="Number of Qubits")

fig.tight_layout()
fig.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/horizontal')
fig.savefig('/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/plots/horizontal.pdf')