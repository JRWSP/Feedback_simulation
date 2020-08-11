# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 22:57:57 2020

@author: jiraw
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams.update({'font.size': 24})

def read_data(filename):
     data = np.genfromtxt(filename+".csv", delimiter = ",")
     data = data[1:]
     data = data.transpose()
     return data

LZ100 = read_data("HZLZ_100")
LZ300 = read_data("HZLZ_300")
LZ500 = read_data("HZLZ_500")

LZ100_ANA = read_data("Analytic_Com_100")
LZ300_ANA = read_data("Analytic_Com_300")
LZ500_ANA = read_data("Analytic_Com_500")
LZ = np.load("HZLZ.npy")

t = 2
step = 2000
dt = t/step
x = np.linspace(0,t,step)
omega = 2*np.pi

s = 20
w = 1
ter = 1200

fig = plt.figure(figsize=(12, 8))
ax1 = plt.subplot(211)

plt.plot(LZ100_ANA[0][:ter], LZ100_ANA[3][:ter], c='g', linewidth=w)
plt.scatter(x[:ter:50], LZ100[3][:ter:50], c='g', label=r"$\tau = 0.1T_{\Omega}$", s=s)
plt.plot(LZ300_ANA[0][:ter], LZ300_ANA[3][:ter], c='y', linewidth=w)
plt.scatter(x[:ter:50], LZ300[3][:ter:50], c='y', label=r"$\tau = 0.3T_{\Omega}$", s=s)
plt.plot(LZ500_ANA[0][:ter], LZ500_ANA[3][:ter], c='r', linewidth=w)
plt.scatter(x[:ter:50], LZ500[3][:ter:50], c='r', label=r"$\tau = 0.5T_{\Omega}$", s=s)
plt.plot(x[:ter], LZ[2][:ter], c='b', label=r"$Lindblad$", linewidth=w)
plt.legend()
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)

plt.plot(x[:ter], LZ[0][:ter], c='b', label="Lindblad", linewidth=w)
plt.plot(LZ100_ANA[0][:ter], LZ100_ANA[1][:ter], c='g', linewidth=w)
plt.scatter(x[:ter:50], LZ100[1][:ter:50], c='g', s=s)
plt.plot(LZ300_ANA[0][:ter], LZ300_ANA[1][:ter], c='y', linewidth=w)
plt.scatter(x[:ter:50], LZ300[1][:ter:50], c='y', s=s)
plt.plot(LZ500_ANA[0][:ter], LZ500_ANA[1][:ter], c='r', linewidth=w)
plt.scatter(x[:ter:50], LZ500[1][:ter:50], c='r', s=s)



plt.xticks([0, 0.25, 0.5, 0.75, 1.0], ["0", r"$0.25T_{\Omega}$", r"$0.50T_{\Omega}$", r"$0.75T_{\Omega}$", r"$1.00T_{\Omega}$"])
plt.xlabel(r"$Time$")
ax1.set_ylabel(r"$S_x$")
ax2.set_ylabel(r"$S_z$")
plt.tight_layout()
plt.savefig("HZLZ_n", dpi=150)