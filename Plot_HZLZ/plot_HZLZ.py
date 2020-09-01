# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 22:57:57 2020

@author: jiraw
"""

import numpy as np
from qutip import *
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
Omega = 2*np.pi

H = 0.5 * Omega * sigmaz() #Hamiltonian
L = sigmaz() #Coupling Operator
c_ops = []
e_ops = [sigmax(), sigmay(), sigmaz()]
rho0 = 0.5 * (1 + (sigmax() + sigmay())/np.sqrt(2) ) # Initial state
result = mesolve(H, rho0, x, c_ops, e_ops)
Mx, My, Mz = result.expect[0], result.expect[1], result.expect[2]


s = 20
w = 1
ter = 1200

fig = plt.figure(figsize=(12, 8))
palette = plt.get_cmap('Set1')
"""
ax1 = plt.subplot(211)
plt.plot(LZ100_ANA[0][:ter], LZ100_ANA[3][:ter], c='g', linewidth=w)
plt.scatter(x[:ter:50], LZ100[3][:ter:50], c='g', s=s)
plt.plot(LZ300_ANA[0][:ter], LZ300_ANA[3][:ter], c='y', linewidth=w)
plt.scatter(x[:ter:50], LZ300[3][:ter:50], c='y', s=s)
plt.plot(LZ500_ANA[0][:ter], LZ500_ANA[3][:ter], c='r', linewidth=w)
plt.scatter(x[:ter:50], LZ500[3][:ter:50], c='r', s=s)
plt.plot(x[:ter], LZ[2][:ter], c='b', label=r"$Lindblad$", linewidth=w)
plt.legend()
plt.setp(ax1.get_xticklabels(), visible=False)
ax2 = plt.subplot(212, sharex=ax1, sharey=ax1)
"""

plt.plot(LZ100_ANA[0][:ter], LZ100_ANA[1][:ter], color=palette(0))
plt.scatter(x[:ter:50], LZ100[1][:ter:50], color=palette(0), s=s, label=r"$\tau = 0.1T_{\Omega}$")
plt.plot(LZ300_ANA[0][:ter], LZ300_ANA[1][:ter], color=palette(1), linewidth=w)
plt.scatter(x[:ter:50], LZ300[1][:ter:50], color=palette(1), s=s, label=r"$\tau = 0.3T_{\Omega}$")
plt.plot(LZ500_ANA[0][:ter], LZ500_ANA[1][:ter], color=palette(2), linewidth=w)
plt.scatter(x[:ter:50], LZ500[1][:ter:50], color=palette(2), s=s, label=r"$\tau = 0.5T_{\Omega}$")
plt.plot(x[:ter], LZ[0][:ter], '--', color=palette(3), label="Lindblad", linewidth=w)
plt.plot(x[:ter], Mx[:ter], '--', color=palette(4), label=r"$\tau = 0$")

plt.xticks([0, 0.25, 0.5, 0.75, 1.0], ["0", r"$0.25$", r"$0.50$", r"$0.75$", r"$1.00$"])
plt.xlabel(r"Time (units of $T_{\Omega})$")
#ax1.set_ylabel(r"$S_z$")
plt.ylabel(r"$S_x$")
plt.legend(loc=0, fontsize=24)
plt.tight_layout()
#plt.savefig("HZLZ_n", dpi=150)