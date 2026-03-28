"""Plot the navigation field -grad(u) as streamlines."""

import numpy as np
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.tri import LinearTriInterpolator

BUILD = "../../build/developers/RobotNavigation"

# Read the Gmsh mesh file
mesh_data = meshio.read("meshes/room.msh")

# Extract (x, y) coordinates of each node
# we drop the z-column because it is always 0 in 2D
node_coords = mesh_data.points[:, :2]

# Extract triangle connectivity
# each row has 3 node indices forming a triangle
triangles = mesh_data.cells_dict["triangle"]

# Load the FEM solution from the csv
# columns are [x, y, u, grad_u_x, grad_u_y]
solution = np.genfromtxt(f"{BUILD}/solution_data.csv", delimiter=",", skip_header=1)
potential = solution[:, 2]
neg_grad_u_x = solution[:, 3]  # -∂u/∂x at each node
neg_grad_u_y = solution[:, 4]  # -∂u/∂y at each node

# Build a matplotlib triangulation object
triangulation = mtri.Triangulation(node_coords[:, 0], node_coords[:, 1], triangles)

fig, ax = plt.subplots(figsize=(12, 6))

# Filled contour plot of the potential field
contour_filled = ax.tricontourf(triangulation, potential, levels=30, cmap="viridis")
fig.colorbar(contour_filled, ax=ax, label=r"$u$", shrink=0.6)

# matplotlib's streamplot requires data on a regular grid, but our FEM
# solution lives on an unstructured triangular mesh. We use
# LinearTriInterpolator to evaluate the fields at regular grid points
# via barycentric interpolation within each triangle.
grid_x = np.linspace(0.01, 9.99, 80)
grid_y = np.linspace(0.01, 5.99, 48)
grid_X, grid_Y = np.meshgrid(grid_x, grid_y)

flux_x = LinearTriInterpolator(triangulation, neg_grad_u_x)(grid_X, grid_Y)
flux_y = LinearTriInterpolator(triangulation, neg_grad_u_y)(grid_X, grid_Y)

ax.streamplot(
    grid_x, grid_y, flux_x, flux_y, color="k", density=1.8, linewidth=1.0, arrowsize=1.5
)

ax.plot([0.5, 2.0], [0, 0], "r-", linewidth=4, label="doors")
ax.plot([7.0, 9.0], [6, 6], "r-", linewidth=4)
ax.legend(loc="lower right")
ax.set_aspect("equal")
ax.set_title(r"Navigation field $-\nabla u$")
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
plt.tight_layout()
plt.savefig("navigation.png", dpi=150)
print("Saved navigation.png")
