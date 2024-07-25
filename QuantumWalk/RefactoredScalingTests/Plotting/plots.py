import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os

#import extracted data values
from extract import numbers, path, SA, QA

#font specifications for plots
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