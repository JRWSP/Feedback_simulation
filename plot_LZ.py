# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 22:57:57 2020

@author: jiraw
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams.update({'font.size': 30})

def read_data(filename):
     data = np.genfromtxt(filename, delimiter = ",")
     data = data[1:]
     data = data.transpose()
     return data

LZ100 = np.load("LZ_100.npy")
LZ300 = np.load("LZ_300.npy")
LZ500 = np.load("LZ_500.npy")

LZ100_ANA = read_data("LZ_Analytic_Com_100.csv")
LZ300_ANA = read_data("LZ_Analytic_Com_300.csv")
LZ500_ANA = read_data("LZ_Analytic_Com_500.csv")
LZ = read_data("LZ.csv")
avg = []
avg.append(np.mean(LZ100, axis=0))
avg.append(np.mean(LZ300, axis=0))
avg.append(np.mean(LZ500, axis=0))

t = 2
step = 2000
dt = t/step
x = np.linspace(0,t,step)
omega = 2*np.pi

s = 40
w = 1
fig = plt.figure(figsize=(14, 8))
plt.plot(LZ100_ANA[0], LZ100_ANA[1], c='g', linewidth=w)
plt.scatter(x[::50], avg[0][0][::50], c='g', label=r"$\tau = 0.1T_{\Omega}$", s=s)
plt.plot(LZ300_ANA[0], LZ300_ANA[1], c='y', linewidth=w)
plt.scatter(x[::50], avg[1][0][::50], c='y', label=r"$\tau = 0.3T_{\Omega}$", s=s)
plt.plot(LZ500_ANA[0], LZ500_ANA[1], c='r', linewidth=w)
plt.scatter(x[::50], avg[2][0][::50], c='r', label=r"$\tau = 0.5T_{\Omega}$", s=s)
plt.plot(LZ[0], LZ[1], c='b', label=r"$Lindblad$", linewidth=w)
plt.legend(loc=0, fontsize=24)
plt.xticks([0, 0.5, 1.0, 1.5, 2.0], ["0", r"$0.5T_{\Omega}$", r"$1.0T_{\Omega}$", r"$1.5T_{\Omega}$", r"$2.0T_{\Omega}$"])
plt.xlabel(r"$Time$")
plt.ylabel(r"$S_x$")
plt.tight_layout()
plt.savefig("LZ_n", dpi=150)