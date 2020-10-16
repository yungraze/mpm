import matplotlib.pyplot as plt
import numpy as np
import fenics as fen
import mshr
from fenics.dynamic_solver import dynamic_solver_func

# -------------Create DOLPHIN mesh and define function space----------
ae = 1.0
edge = 2 * ae  # размер области
G = 1.0
m_sun = 1.0
rad_sun = 0.05
ncells = 256

tol = 1E-14

domain = mshr.Rectangle(fen.Point(-edge, -edge), fen.Point(edge, edge))

# generate mesh and determination of the Function Space
mesh = mshr.generate_mesh(domain, ncells)

P1 = fen.FiniteElement('CG', fen.triangle, 1)

V = fen.FunctionSpace(mesh, P1)
Q = fen.FunctionSpace(mesh, 'DG', 0)

# -------------Define boundary condition------------
phi_L = fen.Expression('-G*m_sun/pow(pow(x[1], 2) + pow(edge, 2),0.5)', degree=2, G=G, m_sun=m_sun, edge=edge)
phi_R = fen.Expression('-G*m_sun/pow(pow(x[1], 2) + pow(edge, 2),0.5)', degree=2, G=G, m_sun=m_sun, edge=edge)
phi_H = fen.Expression('-G*m_sun/pow(pow(x[0], 2) + pow(edge, 2),0.5)', degree=2, G=G, m_sun=m_sun, edge=edge)
phi_B = fen.Expression('-G*m_sun/pow(pow(x[0], 2) + pow(edge, 2),0.5)', degree=2, G=G, m_sun=m_sun, edge=edge)


def boundary_L(x, on_boundary):
    return on_boundary and fen.near(x[0], -edge, tol)


def boundary_R(x, on_boundary):
    return on_boundary and fen.near(x[0], edge, tol)


def boundary_H(x, on_boundary):
    return on_boundary and fen.near(x[1], edge, tol)


def boundary_B(x, on_boundary):
    return on_boundary and fen.near(x[1], -edge, tol)


bc_L = fen.DirichletBC(V, phi_L, boundary_L)
bc_R = fen.DirichletBC(V, phi_R, boundary_R)
bc_H = fen.DirichletBC(V, phi_H, boundary_H)
bc_B = fen.DirichletBC(V, phi_B, boundary_B)

bc = [bc_L, bc_R, bc_H, bc_B]

# -------------Define variational problem------------
phi = fen.Function(V)
v = fen.TestFunction(V)

pi = np.pi

rho_sun = fen.Expression('pow(pow(x[0],2)+pow(x[1],2),0.5)<=rad_sun ? 3*m_sun/(4*pi*pow(rad_sun,3)): 0', degree=2,
                         rad_sun=rad_sun, m_sun=m_sun, pi=pi)

J = fen.Expression('pow(x[0]*x[0] + x[1]*x[1],0.5)', degree=2)

Func = (- J * phi.dx(0) * v.dx(0) - J * phi.dx(1) * v.dx(1) - J * 4 * np.pi * rho_sun * v) * fen.dx

# ---------------------Compute solution---------------------
fen.solve(Func == 0, phi, bc)

# --------------------------Ploting------------------------
# get array componets and triangulation :
v = phi.compute_vertex_values(mesh)
x = mesh.coordinates()[:, 0]
y = mesh.coordinates()[:, 1]
t = mesh.cells()

ax = plt.axes()

cm = plt.get_cmap('viridis')
c = ax.tricontourf(x, y, t, v, 10, cmap=cm)
p = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)
# plt.grid()
plt.xlim(-edge, edge)
plt.ylim(-edge, edge)

# -----------Output in the file-------------------
plt.savefig("results/potential_func_one_body.png", bbox_inches="tight")

x_0 = ae
y_0 = 0.0
v_x_0 = 0.0
v_y_0 = 1.0

data_time = []
data_x = []
data_y = []
data_vx = []
data_vy = []

for i in range(1, 4*5, 1):
    print()
    print('**********************')
    print('ITERATION =', i)
    print()

    time_step = 0.001

    data_time.append(time_step * (i - 1))
    data_x.append(x_0)
    data_y.append(y_0)
    data_vx.append(v_x_0)
    data_vy.append(v_y_0)

    # Where to evaluate
    x = np.array([x_0, y_0])
    vel = np.array([v_x_0, v_y_0])

    uu = fen.project(phi, V)
    phidx = fen.project(uu.dx(0), Q)
    phidy = fen.project(uu.dx(1), Q)

    p = fen.Point(x[0], x[1])

    dfdx = phidx(p.x(), p.y())
    dfdy = phidy(p.x(), p.y())

    sol = dynamic_solver_func(1000, time_step * (i - 1), time_step * i, dfdx, dfdy, x[0], vel[0], x[1], vel[1])
    x_0 = sol[0]
    y_0 = sol[1]
    v_x_0 = sol[2]
    v_y_0 = sol[3]

fig = plt.figure(figsize=(8, 8), facecolor='pink', frameon=True)
plt.plot(data_time, data_x)
plt.grid()
plt.savefig('results/X_from_time.png')

fig = plt.figure(figsize=(8, 8), facecolor='pink', frameon=True)
plt.plot(data_time, data_y)
plt.grid()
plt.savefig('results/Y_from_time.png')

fig = plt.figure(figsize=(8, 8), facecolor='pink', frameon=True)
plt.plot(data_x, data_y)
plt.grid()
plt.savefig('results/trajectory{}.png'.format(ncells))
