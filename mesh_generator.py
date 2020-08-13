import fenics as fen
import mshr
import matplotlib.pyplot as plt

# determination of the domain
domain = mshr.Rectangle(fen.Point(0., 0.), fen.Point(5., 5.)) \
         - mshr.Rectangle(fen.Point(2., 1.25), fen.Point(3., 1.75)) \
         - mshr.Circle(fen.Point(1, 4), .25) \
         + mshr.Circle(fen.Point(4, 2), 1.25)

# generate mesh and determination of the Function Space
mesh = mshr.generate_mesh(domain, 20)

V = fen.FunctionSpace(mesh, 'P', 1)

# determination element of the Function Space
u = fen.Function(V)

# --------------------------Plotting------------------------
# get array components and triangulation :
v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:, 0]
y = mesh.coordinates()[:, 1]
t = mesh.cells()

ax = plt.axes()
cm = plt.get_cmap('viridis')

ax.tricontourf(x, y, t, v, 10, cmap=cm)
ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

plt.ylim(0, 5)
plt.xlim(0, 5)

# Output in the file
plt.savefig("results/gen_mesh.png", bbox_inches="tight")

# #--------------------------Plotting PVD------------------------
# #Plotting the solution using the plot command
# plot(u)
# plot(mesh)
#
# #Save solution to file in VTK format
# vtkfile = File('The heat equation/The heat equation.pvd')
# vtkfile << u
