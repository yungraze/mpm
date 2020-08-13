import numpy as np
import matplotlib.pyplot as plt
import fenics as fen
import mshr

# #-------------Create DOLPHIN mesh and define function space----------
ae = 1.0
edge = 2 * ae  # размер области
G = 1.0
m_sun = 1.0
rad_sun = 0.1
rad_earth = 0.2
m_earth = 1.5

tol = 1E-14

domain = mshr.Rectangle(fen.Point(-edge, -edge), fen.Point(edge, edge))

# generate mesh and determination of the Function Space
mesh = mshr.generate_mesh(domain, 32)

P1 = fen.FiniteElement('CG', fen.triangle, 1)

V = fen.FunctionSpace(mesh, P1)
Q = fen.FunctionSpace(mesh, 'DG', 0)

# -------------Define boundary condition------------
phi_L = fen.Expression('-G*m_sun/pow(pow(x[1],2) + pow(2*ae,2),0.5) - G*m_earth/pow(pow(x[1],2) + pow(3*ae,2),0.5)',
                       degree=3, G=G, m_sun=m_sun, m_earth=m_earth, ae=ae)
phi_R = fen.Expression('-G*m_sun/pow(pow(x[1],2) + pow(2*ae,2),0.5) - G*m_earth/pow(pow(x[1],2) + pow(ae,2),0.5)',
                       degree=3, G=G, m_sun=m_sun, m_earth=m_earth, ae=ae)
phi_H = fen.Expression(
    '-G*m_sun/pow(pow(x[0],2) + pow(2*ae,2),0.5) - G*m_earth/pow(pow(x[0] - ae,2) + pow(2*ae,2),0.5)', degree=3, G=G,
    m_sun=m_sun, m_earth=m_earth, ae=ae)
phi_B = fen.Expression(
    '-G*m_sun/pow(pow(x[0],2) + pow(2*ae,2),0.5) - G*m_earth/pow(pow(x[0] - ae,2) + pow(2*ae,2),0.5)', degree=3, G=G,
    m_sun=m_sun, m_earth=m_earth, ae=ae)


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

Pi = np.pi

rho_sun = fen.Expression('pow(pow(x[0],2) + pow(x[1],2),0.5) <= rad_sun ? 3*m_sun/(4*Pi*pow(rad_sun,3)) : 0',
                         degree=3,
                         rad_sun=rad_sun,
                         m_sun=m_sun,
                         Pi=Pi)
rho_earth = fen.Expression('pow(pow(x[0]-ae,2) + pow(x[1],2),0.5) <= rad_earth ? 3*m_earth/(4*Pi*pow(rad_earth,3)) : 0',
                           degree=3,
                           rad_earth=rad_earth,
                           m_earth=m_earth,
                           Pi=Pi,
                           ae=ae)

J = fen.Expression('pow(x[0]*x[0] + x[1]*x[1],0.5)', degree=2)

Func = (- J * phi.dx(0) * v.dx(0)
        - J * phi.dx(1) * v.dx(1)
        + J * 4 * np.pi * rho_sun * v
        + 4 * np.pi * rho_earth * v) * fen.dx

# ---------------------Compute solution---------------------
fen.solve(Func == 0, phi, bc)

# --------------------------Plotting------------------------
# get array components and triangulation :
v = phi.compute_vertex_values(mesh)
x = mesh.coordinates()[:, 0]
y = mesh.coordinates()[:, 1]
t = mesh.cells()

ax = plt.axes()

cm = plt.get_cmap('viridis')
c = ax.tricontourf(x, y, t, v, 10, cmap=cm)
p = ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

# -----------Output in the file-------------------
plt.savefig('results/potential_func.png', bbox_inches="tight")
