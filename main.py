from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as plb

# #-------------Create DOLPHIN mesh and define function space----------
ae = 15
edge = 1 * ae
tol = 1E-14

# Define domain
domain = Rectangle(Point(-edge, -edge), Point(edge, edge))

# generate mesh and determination of the Function Space
mesh = generate_mesh(domain, 6)

V = FunctionSpace(mesh, 'P', 1)

# -------------Define boundary condition------------
G = 1
m_sun = 1.9
rad_sun = 3

n = V.dim()
print('number_of_cels = ', n)

# размерность пространства (число измерений)
d = mesh.geometry().dim()
print('dim_domain =', d)

# массив координат узлов
F_dof_coordinates = V.tabulate_dof_coordinates().reshape(n, d)
F_dof_coordinates.resize((n, d))

# массив x-координат узлов и y-ооординат узлов
F_dof_x = F_dof_coordinates[:, 0]
F_dof_y = F_dof_coordinates[:, 1]

Phi = [[0] * 5 for k in range(n)]

for i in range(0, n, 1):
    for j in range(0, n, 1):
        phi = - G * m_sun / (F_dof_x[i] ** 2 + F_dof_y[j] ** 2) ** (0.5)
        d_phi_x = G * m_sun * F_dof_x[i] / (F_dof_x[i] ** 2 + F_dof_y[j] ** 2) ** (1.5)
        d_phi_y = G * m_sun * F_dof_y[j] / (F_dof_x[i] ** 2 + F_dof_y[j] ** 2) ** (1.5)
        Phi[i][0] = F_dof_x[i]
        Phi[j][1] = F_dof_y[j]
        Phi[j][2] = phi
        Phi[j][3] = d_phi_x
        Phi[j][4] = d_phi_y

for i in range(n):
    print(Phi[i])

# --------------------------Ploting pylab------------------------
# get array componets and triangulation :
u = Function(V)

v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:, 0]
y = mesh.coordinates()[:, 1]
t = mesh.cells()

# print('x =', x)

ax = plb.axes()
cm = plb.get_cmap('viridis')

ax.tricontourf(x, y, t, v, 10, cmap=cm)
ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

plb.ylim(-edge, edge)
plb.xlim(-edge, edge)

# Output in the file
print("gen_mesh.png obtained")
plb.savefig("gen_mesh.%s" % "png", bbox_inches="tight")

# --------------------------Ploting------------------------
# get array componets and triangulation :
# v = phi.compute_vertex_values(mesh)
# x = mesh.coordinates()[:,0]
# y = mesh.coordinates()[:,1]
# t = mesh.cells()
#
# from pylab import *
# from mpl_toolkits.axes_grid1 import make_axes_locatable
#
# ax = plt.axes()
#
# cm = get_cmap('viridis')
# c = ax.tricontourf(x, y, t, v, 10, cmap = cm)
# p = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)
#
# #-----------Output in the file-------------------
# for ext in ["png"]:
#     plt.savefig("results/potential_func.%s" % (ext,), bbox_inches="tight")

# import matplotlib.pyplot as plt
#
# data_x=[]
# y0=0
# data_phi=[]
#
# for i in range(100):
#     x1=0.0001+0.01*i*(ae-0.0001)
#     p=Point(x1,y0)
#     F1=phi(p.x(),p.y())
#     data_x.append(x1)
#     data_phi.append(F1)
#
# fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
# plt.plot(data_x,data_phi)
# plt.xlim(-ae,ae)
# plt.savefig('results/Phi.png')
