from __future__ import division

# Importing various packages
from math import exp, sqrt
from random import random, seed
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys
#Trial wave function for quantum dots in two dims
def WaveFunction(r):
    r1 = r[0,0]**2 + r[0,1]**2 + r[0, 2]**2
    r2 = r[1,0]**2 + r[1,1]**2 + r[1,2]**2
    r12 = sqrt((r[0,0]-r[1,0])**2 + (r[0,1]-r[1,1])**2 + (r[0,2]-r[1,2])**2)
    return exp(-2*(r1+r2)+r12/2)

#Local energy  for quantum dots in two dims, using analytical local energy
def LocalEnergy(r):
    
    r1 = (r[0,0]**2 + r[0,1]**2 + r[0,2]**2)
    r2 = (r[1,0]**2 + r[1,1]**2 + r[1,2]**2)
    r12 = sqrt((r[0,0]-r[1,0])**2 + (r[0,1]-r[1,1])**2 + (r[0,2]-r[1,2])**2)
    x21 = r[1,0] - r[0,0]
    y21 = r[1,1] - r[0,1]
    z21 = r[1,2] - r[0,2]
    x1 = r[1,0]/r2 - r[0,0]/r1
    y1 = r[1,1]/r2 - r[0,1]/r1
    z1 = r[1,2]/r2 - r[0,2]/r1
    return ((x21 * x1 + y21 * y1 + z21*z1)/r12 - 11/4)

def metropolis():
    R0 = np.zeros((2,3))
    R1 = np.zeros((2,3))
    E = []
    for i in range(2):
        for j in range(3):
            R0[i][j] = random()
    wv0 = WaveFunction(R0)**2
    le0 = LocalEnergy(R0)
    delta = 0.1
    for i in range(2):
        for j in range(3):
            R1[i][j] = R0[i][j] + delta
    wv1 = WaveFunction(R1)**2
    x = wv0/wv1
    E.append(le0)
    for counter in range(1000):
        t = random()
        if x > t:
            for i in range(2):
                for j in range(3):
                    R0[i][j] = R1[i][j]
            if x > 0.5:
                y = delta * 0.01
                delta = delta + y
                for i in range(2):
                    for j in range(3):
                        R1[i][j] = R1[i][j] + y
            if x < 0.5:
                y = delta * 0.01
                delta = delta - y
                for i in range(2):
                    for j in range(3):
                        R1[i][j] = R1[i][j] - y
               
            le0 = LocalEnergy(R1)
            E.append(le0)
            wv00 = WaveFunction(R0)**2
            wv11 = WaveFunction(R1)**2
            if wv11 == 0:
                continue
            x = wv00/wv11
    return E
    
EE = metropolis()
print(EE)            
mean = sum(EE)/len(EE)
s = 0
for v in range(len(EE)):
    s = s + EE[v]**2 - mean**2
s =s/len(EE)
print(mean)
print(s)
plt.plot(EE)
plt.title("Energy vesus number of Iterations")
plt.legend()
plt.show()

