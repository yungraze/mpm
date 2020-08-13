import matplotlib.pyplot as plt
import numpy as np

L = 2.0
E = 5 * np.pi**2

nodes = [0, L, 2*L]

xp = 3.0
mass_p = 15.0
vol_p = 1.0
vel_p = 1.0

s = 0
q = mass_p * vel_p

time = 50
dtime = 0.05
t = 0

ta = []
va = []
xa = []

while t < time:

    N1 = 0
    N2 = 0
    N3 = 0

    if (xp >= -L)&(xp <= L):
        N1 = 1 - np.absolute(xp - nodes[0])/L

    if (xp >= 0)&(xp <= 2*L):
        N2 = 1 - np.absolute(xp - nodes[1])/L

    if (xp >= L)&(xp <= 3*L):
        N3 = 1 - np.absolute(xp - nodes[2])/L

    dN1 = 0
    dN2 = 0
    dN3 = 0

    if (xp >= 0)&(xp <= L):
        dN1 = -1/L
        dN2 = 1/L

    if (xp > L)&(xp <= 2*L):
        dN2 = -1/L
        dN3 = 1/L

    if (xp > 2*L)&(xp <= 3*L):
        dN3 = -1/L

    #Разбиение массы точки по узлам
    m1 = N1 * mass_p
    m2 = N2 * mass_p
    m3 = N3 * mass_p

    #Разбиение импульса по узлам
    mv1 = N1 * q
    mv2 = N2 * q
    mv3 = N3 * q

    mv1 = 0

    f1 = - vol_p * s * dN1
    f2 = - vol_p * s * dN2
    f3 = - vol_p * s * dN3

    f1 = 0
    #f3 = 0

    mv1 = mv1 + f1 * dtime
    mv2 = mv2 + f2 * dtime
    mv3 = mv3 + f3 * dtime

    Nm1 = 0
    Nm2 = 0
    Nm3 = 0

    if m1 > 0 :
        Nm1 = N1/m1

    if m2 > 0 :
        Nm2 = N2/m2

    if m3 > 0 :
        Nm3 = N3/m3

    vel_p = vel_p + dtime * (Nm1 * f1 + Nm2 * f2 + Nm3 * f3)  #измененное
    #vel_p = vel_p + dtime * (Nm1 * f1 + Nm2 * f2 + Nm3 * f3) #исходное
    #xp = xp + dtime * vel_p
    xp = xp + dtime * (Nm1 * mv1 + Nm2 * mv2 + Nm3 * mv3)

    q = mass_p * vel_p

    v1 = Nm1 * q
    v2 = Nm2 * q
    v3 = Nm3 * q

    v1 = 0
    #v3 = 0

    Lp = dN1 * v1 + dN2 * v2 + dN3 * v3
    dEps = dtime * Lp

    s = s + E * dEps

    ta.append(t)
    va.append(vel_p)
    xa.append(xp)

    #print('t=',t,'N1=',N1,'N2=',N2,'N3=',N3,'xp=',xp)
    #print('dN1=',dN1,'dN2=',dN2,'dN3=',dN3)
    #print('Nm1=',Nm1,'Nm2=',Nm2,'Nm3=',Nm3)
    print('t=',t)
    print('f1=',f1,'f2=',f2,'f3=',f3)
    print('mv1=',mv1,'mv2=',mv2,'mv3=',mv3)
    print('Nm1=',Nm1,'Nm2=',Nm2,'Nm3=',Nm3)
    print('delta_v (with N)=',N1 * f1 + N2 * f2 + N3 * f3)
    #print('delta_v (with Nm)=',Nm1 * f1 + Nm2 * f2 + Nm3 * f3)
    #print('v1=',v1,'v2=',v2,'v3=',v3)
    print('Lp=',Lp,'s=',s)
    #print('t=',t,' N1+N2+N3=',N1+N2+N3)
    print()
    t = t + dtime

plt.plot(ta, xa)
plt.show()
