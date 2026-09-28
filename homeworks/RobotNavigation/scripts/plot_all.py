"""Run all plotting scripts.

NPDE homework RobotNavigation
Author: Erick Schulz
Date: 28.03.2026
Developed at ETH Zurich
"""

import subprocess
import sys

for script in [
    "plot_mesh.py",
    "plot_potential.py",
    "plot_navigation.py",
    "plot_path.py",
]:
    subprocess.run([sys.executable, f"scripts/{script}"], check=True)
