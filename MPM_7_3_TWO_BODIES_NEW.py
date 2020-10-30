# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 00:51:00 2020

@author: Артур
"""

import matplotlib.pyplot as plt
import numpy as np
from sympy import sin

# constants
rho = 1000


numelem = 2
pCount = 2 # 400 # number of particles

#массив размера pCount, состоящий из единиц
Mp = np.ones(pCount,dtype = np.float) # mass
Vp = np.ones(pCount,dtype = np.float) # volume
Fp = np.ones(pCount,dtype = np.float) # gradient deformation
s = np.zeros((pCount,4),dtype = np.float) # stress
eps = np.zeros((pCount,3),dtype = np.float) # strain
vp = np.zeros((pCount,2),dtype = np.float) # velocity
xp = np.zeros((pCount,2),dtype = np.float) # position

element_1 = [[0,1,0],[0,0,1]] # Треугольник (0,0)-(1,0)-(0,1)
element_2 = [[5,5,4],[5,4,5]]  # Треугольник (5,5)-(5,4)-(4,5)
elements = [element_1,element_2]

# initialize particle position, mass ,volume, velocity
for e in range(0,numelem): # если numelem = 2, то e = 0, e = 1
    coord = elements[e] # элементы - треугольники
    x = coord[0] # абсциссы узлов данного элемента
    y = coord[1] # ординаты узлов данного элемента
    a = 1/2*np.absolute((x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]))
    print('e=',e)
    print('coord=',coord)
    print('x=',x,' y=',y)
    print('a=',a)
    #coord = node1(element1(e,:),:)
    #a = det([coord,[1;1;1]])/2
    Vp[e] = a
    Mp[e] = a*rho
    #xp[e] = mean(coord)
    
    # this is partcular for the two impact discs example
    
    #if xp(e,1) < 0.5:
    #    vp(e,:) = [v v]
    #else:
    #    vp(e,:) = [-v -v]
    #Fp(e,:) = [1 0 0 1]
Vp0 = Vp

nodeCount = 10
# build the structural grid
# initialize nodal data
nmass = np.zeros((nodeCount,1),dtype = np.float) # nodal mass vector
nmomentum = np.zeros((nodeCount,2),dtype = np.float) #nodal momentum vector
niforce = np.zeros((nodeCount,2),dtype = np.float) #nodal internal force vector
neforce = np.zeros((nodeCount,2),dtype = np.float)



    









