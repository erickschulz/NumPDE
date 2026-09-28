"""Plot the room mesh with doors highlighted.

NPDE homework RobotNavigation
Author: Erick Schulz
Date: 28.03.2026
Developed at ETH Zurich
"""

import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# Read the Gmsh mesh file
# meshio parses it into node coordinates and cells
mesh_data = meshio.read("meshes/room.msh")

# Extract (x, y) coordinates of each node
# we drop the z-column because it is always 0 in 2D
node_coords = mesh_data.points[:, :2]

# Extract triangle connectivity
# each row has 3 node indices forming a triangle
triangles = mesh_data.cells_dict["triangle"]

# Build a matplotlib triangulation object for plotting the mesh
triangulation = mtri.Triangulation(node_coords[:, 0], node_coords[:, 1], triangles)

fig, ax = plt.subplots(figsize=(10, 6))
ax.triplot(triangulation, linewidth=0.3, color="k")
ax.plot([0.5, 2.0], [0, 0], "r-", linewidth=4, label="doors")
ax.plot([7.0, 9.0], [6, 6], "r-", linewidth=4)
ax.legend(loc="lower right")
ax.set_aspect("equal")
ax.set_title("Room mesh")
plt.tight_layout()
plt.savefig("mesh.png", dpi=150)
print("Saved mesh.png")
