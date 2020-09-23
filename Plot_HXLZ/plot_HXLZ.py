#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:30:00 2019

@author: quantuminw
"""

import numpy as np
from qutip import *
import matplotlib.pyplot as plt
import string
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams.update({'font.size': 24})

def read_data(filename):
     data = np.genfromtxt(filename, delimiter = ',', dtype=str)
     data = data.tolist()
     del data[0]
     for ii in range(4):
          for jj in range(len(data)):
               data[jj][ii] = complex(data[jj][ii])
               data[jj][ii] = np.real(data[jj][ii])
     data2 = list(map(list, zip(*data)))
     return np.array(data2)

t = 1.25
step = 1250
dt = t/step
x = np.linspace(0,t,step)
Omega = 2 * np.pi *10
H = 0.5 * Omega * sigmax() #Hamiltonian
L = sigmaz() #Coupling Operator
c_ops = [L]
e_ops = [sigmax(), sigmay(), sigmaz()]
rho0 = 0.5 * (1 + (sigmax() + sigmay())/np.sqrt(2) ) # Initial state
n = 5000


result = mesolve(H, rho0, x, c_ops, e_ops)
Lx, Ly, Lz = result.expect[0], result.expect[1], result.expect[2]
result = mesolve(H, rho0, x, [], e_ops)
Mx, My, Mz = result.expect[0], result.expect[1], result.expect[2]

Sim100 = np.load("w10_100dt.npy")
Sim100 = np.sum(Sim100, axis=0)/n
Sim125 = np.load("w10_125dt.npy")
Sim125 = np.sum(Sim125, axis=0)/n
Sim150 = np.load("w10_150dt.npy")
Sim150 = np.sum(Sim150, axis=0)/n
Sim175 = np.load("w10_175dt.npy")
Sim175 = np.sum(Sim175, axis=0)/n
Sim200 = np.load("w10_200dt.npy")
Sim200 = np.sum(Sim200, axis=0)/n
Sim250 = np.load("w10_250dt.npy")
Sim250 = np.sum(Sim250, axis=0)/n
#Sim0 = np.load("HXLZ0_Traj_10w.npy")

l=1
j = 500
k = 10
size = 50
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(20, 14))
#plt.title("H = σ x")
palette = plt.get_cmap('Set1')

ax1.plot(x[:j], Lx[:j], '--', linewidth=l, label=r"$Lindblad$",  color=palette(0))
ax1.plot(x[:j], Mx[:j], '--', linewidth=l, label=r"$\tau = 0$", color=palette(1))

ax1.scatter(x[:j:k], Sim100[0][:j:k], marker='.', s=size, label=r"$\tau = 1.0T_{\Omega}$", color=palette(2))
#ax1.plot(x[:j], Sim100[0][:j], linewidth=1, color=palette(2))
ax1.scatter(x[:j:k], Sim150[0][:j:k], marker='*', s=size, label=r"$\tau = 1.5T_{\Omega}$", color=palette(3))
#ax1.plot(x[:j], Sim150[0][:j], linewidth=1, color=palette(3))
ax1.scatter(x[:j:k], Sim175[0][:j:k], marker='x', s=size, label=r"$\tau = 1.75T_{\Omega}$", color=palette(4))
#ax1.plot(x[:j], Sim175[0][:j], linewidth=1, color=palette(4))
ax1.scatter(x[:j:k], Sim200[0][:j:k], marker='.', s=size, label=r"$\tau = 2.0T_{\Omega}$", color=palette(6))
#ax1.plot(x[:j], Sim200[0][:j], linewidth=1, color=palette(6))
ax1.scatter(x[:j:k], Sim250[0][:j:k], marker='*', s=size, label=r"$\tau = 2.5T_{\Omega}$", color=palette(7))
#ax1.plot(x[:j], Sim250[0][:j], linewidth=1, color=palette(7))

ax1.axvline(x[100], linewidth=1, color=palette(2))
ax1.axvline(x[150], linewidth=1, color=palette(3))
ax1.axvline(x[200], linewidth=1, color=palette(6))
ax1.axvline(x[250], linewidth=1, color=palette(7))
ax1.axvline(x[175], linewidth=1, color=palette(4))
plt.setp(ax1.get_xticklabels(), visible=False)

ax1.set_ylabel(r"$S_x$")
lgd=ax1.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
ax1.set_ylim([-0.1, 1.0])
ax1.text(-0.09, 0.9, "(a)", transform=ax1.transAxes, size=20)

ax2.plot(x[:j], Lz[:j], '--', linewidth=l, label=r"$Lindblad$", color=palette(0))
ax2.plot(x[:j], Mz[:j], '--', linewidth=l, label=r"$\tau = 0$", color=palette(1))
ax2.scatter(x[:j:k], Sim100[2][:j:k], marker='.', s=size, color=palette(2))
ax2.plot(x[:j], Sim100[2][:j], linewidth=1, color=palette(2))

ax2.scatter(x[:j:k], Sim200[2][:j:k], marker='.', s=size, color=palette(6))
ax2.plot(x[:j], Sim200[2][:j], linewidth=1, color=palette(6))
ax2.scatter(x[:j:k], Sim150[2][:j:k], marker='*', s=size, color=palette(3))
ax2.plot(x[:j], Sim150[2][:j], linewidth=1, color=palette(3))
ax2.scatter(x[:j:k], Sim250[2][:j:k], marker='*', s=size, color=palette(7))
ax2.plot(x[:j], Sim250[2][:j], linewidth=1, color=palette(7))

plt.setp(ax2.get_xticklabels(), visible=False)
ax2.axvline(x[100], linewidth=1, color=palette(2))
ax2.axvline(x[200], linewidth=1, color=palette(6))
ax2.axvline(x[150], linewidth=1, color=palette(3))
ax2.axvline(x[250], linewidth=1, color=palette(7))
ax2.set_ylabel(r"$S_z$")
ax2.set_ylim([-1.0, 1.0])
ax2.text(-0.09, 0.9, "(b)", transform=ax2.transAxes, size=20)

ax3.scatter(x[:j:k], Sim175[2][:j:k], marker='x', s=size, color=palette(4), label=r"$\tau = 1.75T_{\Omega}$")
ax3.plot(x[:j], Sim175[2][:j], linewidth=1, color=palette(4), label='_nolegend_')
ax3.axvline(x[175], linewidth=1, color=palette(4), label='_nolegend_')
ax3.plot(x[:j], Lz[:j], '--', linewidth=l, label='_nolegend_', color=palette(0))
ax3.plot(x[:j], Mz[:j], '--', linewidth=l, label='_nolegend_', color=palette(1))
ax3.legend()
plt.setp(ax1.get_xticklabels(), fontsize=6)

ax3.set_ylabel(r"$S_z$")
ax3.set_ylim([-1.0, 1.0])
ax3.text(-0.09, 0.9, "(c)", transform=ax3.transAxes, size=20)

plt.xlabel(r"Time (units of $T_{\Omega})$")

plt.xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5], ["0", r"$1.0$", r"$2.0$", r"$3.0$", r"$4.0$", r"$5.0$"])
plt.tight_layout(pad=1.0)
#plt.show()
plt.savefig("HXLZ_n", dpi=150, bbox_extra_artists=(lgd,), bbox_inches='tight')
