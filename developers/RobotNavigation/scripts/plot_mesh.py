"""Quick mesh visualization to validate the gmsh geometry."""
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def read_msh2(path):
    """Minimal reader for gmsh MSH2 ASCII format."""
    with open(path) as f:
        lines = f.readlines()

    # Parse nodes
    i = lines.index("$Nodes\n") + 1
    n_nodes = int(lines[i])
    nodes = np.zeros((n_nodes, 2))
    for j in range(n_nodes):
        parts = lines[i + 1 + j].split()
        nodes[j] = [float(parts[1]), float(parts[2])]

    # Parse elements — keep only triangles (type 2)
    i = lines.index("$Elements\n") + 1
    n_elems = int(lines[i])
    triangles = []
    for j in range(n_elems):
        parts = lines[i + 1 + j].split()
        elem_type = int(parts[1])
        if elem_type == 2:  # 3-node triangle
            n_tags = int(parts[2])
            idx = 3 + n_tags
            tri = [int(parts[idx]) - 1,
                   int(parts[idx + 1]) - 1,
                   int(parts[idx + 2]) - 1]
            triangles.append(tri)

    return nodes, np.array(triangles)


nodes, tris = read_msh2("meshes/room.msh")
triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tris)

fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.triplot(triang, linewidth=0.3, color="k")
ax.set_aspect("equal")
ax.set_title("Room mesh")

# Mark doors in red
door1_x = [0.5, 2.0]
door2_x = [7.0, 9.0]
ax.plot(door1_x, [0, 0], "r-", linewidth=3, label="doors")
ax.plot(door2_x, [6, 6], "r-", linewidth=3)
ax.legend()
plt.tight_layout()
plt.savefig("mesh.png", dpi=150)
print("Saved mesh.png")
