"""Plot the potential field u from the Poisson BVP.

NPDE homework RobotNavigation
Author: Erick Schulz
Date: 28.03.2026
Developed at ETH Zurich
"""

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
potential = solution[:, 2]  # u values at each node

# Build a matplotlib triangulation object for plotting
triangulation = mtri.Triangulation(node_coords[:, 0], node_coords[:, 1], triangles)

fig, ax = plt.subplots(figsize=(12, 6))

# Filled contour plot of the potential field
contour_filled = ax.tricontourf(triangulation, potential, levels=30, cmap="viridis")
# Overlay contour lines for readability
ax.tricontour(triangulation, potential, levels=15, colors="k", linewidths=0.3)

fig.colorbar(contour_filled, ax=ax, label=r"$u$", shrink=0.6)

ax.plot([0.5, 2.0], [0, 0], "r-", linewidth=4, label="doors")
ax.plot([7.0, 9.0], [6, 6], "r-", linewidth=4)
ax.legend(loc="lower right")
ax.set_aspect("equal")
ax.set_title(r"Potential field $u$:   $-\Delta u = 1$ in $\Omega$")
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
plt.tight_layout()
plt.savefig("potential.png", dpi=150)
print("Saved potential.png")
