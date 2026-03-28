"""Plot the robot path overlaid on the potential field."""

import numpy as np
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

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

# Load the robot path computed by RK4 integration of dx/dt = -grad(u)
# columns are [x, y]
path = np.genfromtxt(f"{BUILD}/path_data.csv", delimiter=",", skip_header=1)

# Build a matplotlib triangulation object for plotting
triangulation = mtri.Triangulation(node_coords[:, 0], node_coords[:, 1], triangles)

fig, ax = plt.subplots(figsize=(12, 6))

# Semi-transparent potential field as background
contour_filled = ax.tricontourf(
    triangulation, potential, levels=30, cmap="viridis", alpha=0.5
)
ax.tricontour(triangulation, potential, levels=15, colors="grey", linewidths=0.3)
fig.colorbar(contour_filled, ax=ax, label=r"$u$", shrink=0.6)

# Robot path from start to exit
ax.plot(path[:, 0], path[:, 1], "r-", linewidth=3, zorder=10, label="Robot path (RK4)")
ax.plot(path[0, 0], path[0, 1], "ko", markersize=10, zorder=11, label="Start")
ax.plot(path[-1, 0], path[-1, 1], "r*", markersize=15, zorder=11, label="Exit")

ax.plot([0.5, 2.0], [0, 0], "r-", linewidth=4, label="doors")
ax.plot([7.0, 9.0], [6, 6], "r-", linewidth=4)
ax.legend(loc="upper right", fontsize=10)
ax.set_aspect("equal")
ax.set_title(r"Robot path: $\dot{x} = -\nabla u(x)$ integrated with RK4")
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
plt.tight_layout()
plt.savefig("robot.png", dpi=150)
print("Saved robot.png")
