# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 00:19:24 2020

@author: Артур
"""

###Vibration of a continuum bar

import matplotlib.pyplot as plt
import numpy as np
from sympy import sin

# Изучаем осевую (axial) вибрацию стержня
# с модулем Юнга E = 100 и длиной L = 1

elemCount = 13  # Число элементов
nodeCount = 14  # Число узлов

L = 25

# координаты всех 14-ти узлов сетки
nodes = np.linspace(0, L, elemCount + 1)
print('nodes = ', nodes)

# массив размера (elemCount x 2), состоящий из номеров элементов
elements = np.zeros((elemCount, 2), dtype=np.int16)

for ie in range(0, elemCount):  # заполняем этот массив
    elements[ie, 0] = ie  # левый узел элемента
    elements[ie, 1] = ie + 1  # правый узел
print('elements = ', elements)

# size - размер массива
nodeCount = np.size(nodes)
print(nodeCount)

velo = np.zeros(nodeCount)
print('velo = ', velo)

# номер частицы - центра масс
mld = np.floor(elemCount / 2) + 1
print(mld)

E = 100
L = 25
rho = 1

c = (E / rho) ** (1 / 2)
beta1 = np.pi / 2 / L
omega1 = beta1 * c

deltax = L / elemCount  # расстояние между узлами

# material point at the center of elements
xp = np.zeros(elemCount)
for p in range(0, elemCount):
    xp[p] = 0.5 * (nodes[p] + nodes[p + 1])

print('xp = ', xp)

# Длина массива xp
pCount = np.size(xp)
print(pCount)

Mp = deltax * np.ones(pCount)  # mass
print('Mp = ', Mp)

Vp = deltax * np.ones(pCount)  # volume

Fp = np.ones(pCount)  # gradient deformation

Vp0 = Vp

sp = np.zeros(pCount)  # stress
vp = np.zeros(pCount)  # velocity

# initial velocities
p = 1
while p < pCount:
    vp[p] = 0.1 * sin(beta1 * xp[p])
    p = p + 1

print('vp = ', vp)

dtime = 0.1 * deltax / c
time = 10  # time = 100
t = 0

ta = []
va = []
xa = []

nmass = np.zeros(nodeCount)  # nodal mass vector
nmomentum = np.zeros(nodeCount)  # nodal momentum vector
niforce = np.zeros(nodeCount)  # nodal internal force vector
neforce = np.zeros(nodeCount)  # nodal external force vector

tol = 0.0001

while t < time:

    nmass[:] = 0
    nmomentum[:] = 0
    niforce[:] = 0

    # loop over computational cells of elements
    for e in range(0, elemCount):
        esctr = elements[e, :]  # e-тая строка массива elements
        enode = nodes[esctr]  # координаты двух соседних узлов
        Le = enode[1] - enode[0]  # расстояние между данными узлами
        mpts = [e]  # ??? предположим, что это массив [e]
        # либо, как вариант: (1,nodeCount, nodeCount)
        # !!! в оригинале mpts = mpoints{e} # particles inside element e!!!

        # loop over particles (inside element e)
        for p in range(0, len(mpts)):
            pid = mpts[p]
            N1 = 1 - np.absolute(xp[pid] - enode[0]) / Le
            N2 = 1 - np.absolute(xp[pid] - enode[1]) / Le
            dN1 = -1 / Le
            dN2 = 1 / Le
            # particle mass and momentum to node
            nmass[esctr[0]] = nmass[esctr[0]] + N1 * Mp[pid]
            nmass[esctr[1]] = nmass[esctr[1]] + N2 * Mp[pid]
            nmomentum[esctr[0]] = nmomentum[esctr[0]] + N1 * Mp[pid] * vp[pid]
            nmomentum[esctr[1]] = nmomentum[esctr[1]] + N2 * Mp[pid] * vp[pid]
            # internal force
            niforce[esctr[0]] = niforce[esctr[0]] - Vp[pid] * sp[pid] * dN1
            niforce[esctr[1]] = niforce[esctr[1]] - Vp[pid] * sp[pid] * dN2

    # update nodal momenta
    nmomentum[0] = 0  # boundary conditions f1 = m1*a1 = 0, a1 = 0
    niforce[0] = 0

    for i in range(0, nodeCount):
        nmomentum[i] = nmomentum[i] + niforce[i] * dtime

    # update particle velocity and position and stresses
    for e in range(0, elemCount):
        esctr = elements[e, :]  # e-тая строка массива elements
        enode = nodes[esctr]  # координаты двух соседних узлов
        Le = enode[1] - enode[0]  # расстояние между данными узлами
        mpts = [e]  # ???
        # loop over particles
        for p in range(0, len(mpts)):
            pid = mpts[p]
            N1 = 1 - np.absolute(xp[pid] - enode[0]) / Le
            N2 = 1 - np.absolute(xp[pid] - enode[1]) / Le
            dN1 = -1 / Le
            dN2 = 1 / Le
            v1 = 0
            v2 = 0
            if nmass[esctr[0]] > tol:
                vp[pid] = vp[pid] + dtime * niforce[esctr[0]] / nmass[esctr[0]]
                v1 = nmomentum[esctr[0]] / nmass[esctr[0]]
                xp[pid] = xp[pid] + dtime * N1 * v1
            if nmass[esctr[1]] > tol:
                vp[pid] = vp[pid] + dtime * niforce[esctr[1]] / nmass[esctr[1]]
                v2 = nmomentum[esctr[1]] / nmass[esctr[1]]
                xp[pid] = xp[pid] + dtime * N2 * v2
            # if (esctr[0] == 1) v1 = 0
            Lp = dN1 * v1 + dN2 * v2  # gradient velocity
            Fp[pid] = (1 + Lp * dtime) * Fp[pid]  # gradient deformation
            Vp[pid] = Fp[pid] * Vp0[pid]  # volume
            dEps = dtime * Lp  # strain increment
            sp[pid] = sp[pid] + E * dEps  # stress update
        # ====== end of the iteration (for p) ======
    # ====== end of the iteration (for e) ======

    # store time, velocity for plotting
    ta.append(t)
    va.append(vp[e])
    xa.append(xp[e])

    # advance for the next time step
    t = t + dtime

# ===== end of the iteration (while t < time) =======

plt.plot(ta, xa)
plt.show()
