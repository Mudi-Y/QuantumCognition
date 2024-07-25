import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import re
import os

#helper function to extract values
def extract_value(filename, p):
    with open(filename, mode="rt", encoding="utf-8") as docFile:
        doc = docFile.read()
        val = re.findall(p, doc)
        
        #potentially failed anneal
        if val == []:
            return -1.0
        
        if "{" in val[0]:
            val = float(re.findall(r"[-+]?(?:\d*\.*\d+)", val[0])[0])
        else:
            val = float(val[0].split("  ")[1])
    return val

#Extract Data from Logs
path = "/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Plotting/Plots/" #where plots go
datapath = "/workspaces/QuantumCognition/QuantumWalk/RefactoredScalingTests/Outputs/"

#dataframes
simulatedDF = pd.DataFrame(columns=["hSize", "group", "lambda", "schedule", "data"])
quantumDF = pd.DataFrame(columns=["hSize", "group", "lambda", "schedule", "data"])

#values to look through data and fill in dataframe
attempts = [1,2,3] #number of attempts for each aneal configuration
hSizes = [4, 8, 16]
groups = [5, 6, 7]
lambdas = [-2,-3, -4]
annealSchedules = []
times = [5.0, 10.0, 15.0]
for time in times:
    schedule = [(0.0, 0.0), (time, 0.4), (time + 5.0, 0.4), (time + 5.0 + time, 0.8), (time + 5.0 + time + 5.0, 0.8), (time + 5.0 + time + 5.0 + time, 1.0)]
    annealSchedules.append(schedule)

for size in hSizes:
    for group in groups:
        for lamb in lambdas:
            for schedule in annealSchedules:
                #make dataframes to store data
                SA = pd.DataFrame(index=attempts, columns = ["acc", "time", "acc_avg", "time_avg", "acc_std", "time_std"])
                QA = pd.DataFrame(index=attempts, columns = ["acc", "time", "acc_avg", "time_avg", "acc_std", "time_std"])

                # #make acc and time fields arrays
                # for data in ["acc", "time"]:
                #     for attempt in attempts:
                #         SA[data][attempt] = []
                #         QA[data][attempt] = []

                for method in ["simulated", "QPU"]:
                    for attempt in attempts:
                        fname = datapath + f"HSize-{size}/Group-{group}/Lambda-{lamb}/Schedule-{schedule}/Results/Trial_1/{method}/group{group}/attempt{attempt}"
                        if method == "simulated":
                            #extract and append data for simulated
                            p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                            SA["acc"][attempt] = extract_value(fname, p)

                            p = re.compile("Time:  [-+]?(?:\d*\.*\d+)")
                            SA["time"][attempt] = extract_value(fname, p)*1000 #convert from seconds to miliseconds
                        else:
                            #extract and append data for QPU
                            p = re.compile("Absolute Error:  [-+]?(?:\d*\.*\d+)")
                            QA["acc"][attempt] = extract_value(fname, p)

                            p = re.compile(".*qpu_sampling_time.* [-+]?(?:\d*\.*\d+)")
                            QA["time"][attempt] = extract_value(fname, p)*0.001 #convert from microseconds to miliseconds

                #Fill in Avg, STDEV
                for df in [SA, QA]:
                    df["acc_avg"] = np.mean(df["acc"])
                    df["time_avg"] = np.mean(df["time"])
                    df["acc_std"] = np.std(df["acc"])
                    df["time_std"] = np.std(df["time"])
                
                #create meta dataframes to be appended to main dataframes
                SA_meta = pd.DataFrame({"hSize": f"hSize_{size}",
                                        "group": f"group_{group}",
                                        "lambda": f"lambda_{lamb}",
                                        "schedule": f"schedule_{schedule}"}, index = [0])
                QA_meta = pd.DataFrame({"hSize": f"hSize_{size}",
                                        "group": f"group_{group}",
                                        "lambda": f"lambda_{lamb}",
                                        "schedule": f"schedule_{schedule}"}, index = [0])
                
                simulatedDF = pd.concat([simulatedDF, SA_meta], ignore_index=True)
                quantumDF = pd.concat([quantumDF, QA_meta], ignore_index=True)
