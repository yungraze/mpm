import fenics as fen
import mshr
import matplotlib.pyplot as plt

# -------------Create DOLPHIN mesh and define function space----------
ae = 15
edge = 1 * ae
tol = 1E-14

# Define domain
domain = mshr.Rectangle(fen.Point(-edge, -edge), fen.Point(edge, edge))

# generate mesh and determination of the Function Space
mesh = mshr.generate_mesh(domain, 6)

V = fen.FunctionSpace(mesh, 'P', 1)

# -------------Define boundary condition------------
G = 1
m_sun = 1.9
rad_sun = 3

n = V.dim()
print('number_of_cels = ', n)

# размерность пространства (число измерений)
d = mesh.geometry().dim()
print('dim_domain = ', d)

# массив координат узлов
F_dof_coordinates = V.tabulate_dof_coordinates().reshape(n, d)
F_dof_coordinates.resize((n, d))

# массив x-координат узлов и y-ооординат узлов
F_dof_x = F_dof_coordinates[:, 0]
F_dof_y = F_dof_coordinates[:, 1]

Phi = [[0] * 5 for k in range(n)]

for i in range(0, n, 1):
    for j in range(0, n, 1):
        phi = - G * m_sun / (F_dof_x[i] ** 2 + F_dof_y[j] ** 2) ** 0.5
        d_phi_x = G * m_sun * F_dof_x[i] / (F_dof_x[i] ** 2 + F_dof_y[j] ** 2) ** 1.5
        d_phi_y = G * m_sun * F_dof_y[j] / (F_dof_x[i] ** 2 + F_dof_y[j] ** 2) ** 1.5
        Phi[i][0] = F_dof_x[i]
        Phi[j][1] = F_dof_y[j]
        Phi[j][2] = phi
        Phi[j][3] = d_phi_x
        Phi[j][4] = d_phi_y

for i in range(n):
    print(Phi[i])

# --------------------------Plotting------------------------
# get array components and triangulation :
u = fen.Function(V)

v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:, 0]
y = mesh.coordinates()[:, 1]
t = mesh.cells()

ax = plt.axes()
cm = plt.get_cmap('viridis')

ax.tricontourf(x, y, t, v, 10, cmap=cm)
ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

plt.ylim(-edge, edge)
plt.xlim(-edge, edge)

# Output in the file
plt.savefig("results/gen_mesh.png", bbox_inches="tight")

data_x = []
y0 = 0
data_phi = []

for i in range(100):
    x1 = 0.0001 + 0.01 * i * (ae - 0.0001)
    p = fen.Point(x1, y0)
    F1 = u(p.x(), p.y())
    data_x.append(x1)
    data_phi.append(F1)

fig = plt.figure(figsize=(8, 6), facecolor='pink', frameon=True)
plt.plot(data_x, data_phi)
plt.xlim(-ae, ae)
plt.savefig('results/Phi.png')
