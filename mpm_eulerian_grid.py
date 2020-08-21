import matplotlib.pyplot as plt
import numpy as np

L = 1
E = 10

nodes = [0, L]
elements = [1, 2]

xp = 0.5 * L
mass_p = 1
vol_p = 1
vel_p = 0.1

s = 0
q = mass_p * vel_p

time = 10
dtime = 0.01
t = 0

ta = []
va = []
xa = []

while t < time:
    N1 = 1 - np.absolute(xp - nodes[0])  ###
    N2 = 1 - np.absolute(xp - nodes[1])

    dN1 = - 1 / L
    dN2 = 1 / L

    m1 = N1 * mass_p
    m2 = N2 * mass_p

    mv1 = N1 * q
    mv2 = N2 * q

    mv1 = 0

    fint1 = - vol_p * s * dN1
    fint2 = - vol_p * s * dN2

    f1 = fint1
    f2 = fint2

    f1 = 0

    mv1 = mv1 + f1 * dtime
    mv2 = mv2 + f2 * dtime

    vel_p = vel_p + dtime * (N1 * f1 / m1 + N2 * f2 / m2)
    xp = xp + dtime * (N1 * mv1 / m1 + N2 * mv2 / m2)

    q = mass_p * vel_p

    v1 = N1 * mass_p * vel_p / m1
    v2 = N2 * mass_p * vel_p / m2

    v1 = 0
    Lp = dN1 * v1 + dN2 * v2
    dEps = dtime * Lp

    s = s + E * dEps

    ta.append(t)
    va.append(vel_p)
    xa.append(xp)

    t = t + dtime

plt.plot(ta, va)
plt.show()
