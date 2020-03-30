#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 23:46:46 2018

@author: jirawat
"""

from qutip import *
from qutip.expect import expect_rho_vec
from math import *
import os
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Pool
from tqdm import tqdm, trange
import csv


def save(t, x, y, z, filename):
     data = zip(t, x, y, z)

     with open(filename, "w") as output:
          writer = csv.writer(output)
          writer.writerow(['t', 'sx', 'sy', 'sz'])
          for row in data:
               writer.writerow(row)
          
def Operate(U_op, initial_state):
    return U_op * initial_state * U_op.dag()

def Wiener(dt, step):
    local_rand = np.random.RandomState(None)
    s = local_rand.normal(0.0, np.sqrt(dt), step)  # * np.sqrt(dt)
    s = s.tolist()
    return s
        
def traj_delay(delay, exp_traj_x, exp_traj_y, exp_traj_z):
            #Time delay added
    stock = Wiener(dt, step)
    for ii in range(1, step):
        if ii <= delay:
            if ii == 1:
                new_rho = rho0
            Y = ( (L * e**(1j*theta) + L.dag() * e**(-1j*theta) )*new_rho).tr()/new_rho.tr()
            Y2 = stock[ii]/dt + Y
            M = 1 + ( 1j*L*Y2*dt ) - (0.5 * dt * L**2)
            
            #Evolution
            new_rho = Operate(M, new_rho)
            #new_rho = Operate(U, new_rho)
            new_rho = new_rho.unit()
        else :
            Y = ( (L * e**(1j*theta) + L.dag() * e**(-1j*theta) )*new_rho).tr()/new_rho.tr()
            Y2 = stock[ii]/dt + Y
            Y1 = stock[ii-delay]/dt + Y
            M = 1 + ( 1j*L*Y2*dt ) - (0.5 * dt * L**2)
            F = 1 - ( 1j*L*Y1*dt ) - (0.5 * dt * L**2)
            
            new_rho = Operate(F*M, new_rho)
            #new_rho = Operate(U, new_rho)
            #new_rho = Operate(F, new_rho)
            new_rho = new_rho.unit()
            
        exp_traj_x[ii] =  (new_rho * sigmax()).tr() / new_rho.tr()
        exp_traj_y[ii] =  (new_rho * sigmay()).tr() / new_rho.tr()
        exp_traj_z[ii] =  (new_rho * sigmaz()).tr() / new_rho.tr()

def process(delay):
    # List for expectation results
    exp_traj_x = np.zeros(step)
    exp_traj_y = np.zeros(step)
    exp_traj_z = np.zeros(step)
    #Execute simulation function and print progression    
    for ii in range(nsub):
        traj_delay(delay, exp_traj_x, exp_traj_y, exp_traj_z)
    exp_traj_x[0] = (rho0 * sigmax()).tr()
    exp_traj_y[0] = (rho0 * sigmay()).tr()
    exp_traj_z[0] = (rho0 * sigmaz()).tr()
    return exp_traj_x, exp_traj_y, exp_traj_z

#Define Parameter
#Time, timestep
n_proc = 1#multiprocessing.cpu_count()
list_delay = [100]#, [100], [150]]
list_filename = ["LZ_100.csv"] 
""",
                 "purity_100dt.csv", 
                 "purity_500dt.csv", 
                 "purity_1000dt.csv", 
                 "purity_5000dt.csv", 
                 "purity_10000dt.csv"] 
     """
t = 2
step = 200
dt = t/step
x = np.linspace(0,t,step)
nsub = 1    #Number of trajectories/core
Epsilon = 0 * np.pi
Omega = 2 * np.pi
H = 0.5 * Epsilon * sigmaz() + 0.5 * Omega * sigmax() #Hamiltonian
theta = 0.5 * np.pi #Measurement angle
U = 1 - 1j*H*dt  #Unitary
L = sigmaz() #Coupling Operator
rho0 = 0.5 * (1 + (sigmax() + sigmay())/np.sqrt(2) ) # Initial state

if __name__ == "__main__":
    dat = []
    for indx in tqdm(range(5)): #Iterate list_delay
         t1_U = process(100)
         dat.append(t1_U)
    dat = np.array(dat)
    
    fig = plt.plot()
    plt.ylim([-1, 1])
    plt.ylabel({-1, 1})
    for ii in range(len(dat)):
        plt.plot(x, np.real(dat[ii][0]), linewidth=0.2, c='b')
    plt.ylabel('Sx')
    plt.xlabel('Time')
    plt.title("test")
    plt.legend()
    plt.show()
    #np.save("LZ_100.csv", dat)
    #plt.savefig("MUF", dpi=500)
"""
         t1_total = [[]]*3
         for jj in  range(3):
              t1_total[jj] = t1_U[0][jj]/n_proc
              for ii in range(1, n_proc):
                   t1_total[jj] += t1_U[ii][jj]/n_proc
"""
