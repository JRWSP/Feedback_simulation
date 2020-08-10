#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:30:00 2019

@author: quantuminw
"""

import numpy as np
from qutip import *
import matplotlib.pyplot as plt

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

t = 2
step = 2000
dt = t/step
x = np.linspace(0,t,step)
Omega = 2 * np.pi
H = 0.5 * Omega * sigmax() #Hamiltonian
L = sigmaz() #Coupling Operator
c_ops = []
e_ops = [sigmax(), sigmay(), sigmaz()]
rho0 = 0.5 * (1 + (sigmax() + sigmay())/np.sqrt(2) ) # Initial state

result = mesolve(H, rho0, x, c_ops, e_ops)
Mx, My, Mz = result.expect[0], result.expect[1], result.expect[2]

#Sim100 = read_data("HXLZ_100.csv")
#Sim300 = read_data("HXLZ_300.csv")
#Sim500 = read_data("HXLZ_500.csv")
Sim0 = np.load("HXLZ0.npy")
Sim1 = np.load("HXLZ1_Traj.npy")
Sim100 = np.load("HXLZ100_Traj.npy")
#Sim300 = np.load("HZLZ300_Traj.npy")
#Sim500 = np.load("HZLZ500_Traj.npy")
Tj = np.random.randint(0, 19)

k = 500
b = Bloch()

ini = np.array([Mx[0], My[0], Mz[0]])
b.add_points(ini)
point = [Sim1[Tj][0][:k], Sim1[Tj][1][:k], Sim1[Tj][2][:k]]
b.add_points(point ,'l')
point = [Sim100[Tj][0][:k], Sim100[Tj][1][:k], Sim100[Tj][2][:k]]
b.add_points(point ,'l')
point = [Sim0[Tj][0][:k], Sim0[Tj][1][:k], Sim0[Tj][2][:k]]
b.add_points(point ,'l')
b.point_color = ['b', 'r', 'g', 'b']
#b.point_size  = [15, 25, 5, 5]

b.render(b.fig, b.axes)
b.fig.savefig("Bloch_Intro2.png", dpi=300, transparent=True)