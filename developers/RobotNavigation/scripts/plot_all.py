"""Run all plotting scripts."""

import subprocess
import sys

for script in [
    "plot_mesh.py",
    "plot_potential.py",
    "plot_navigation.py",
    "plot_path.py",
]:
    subprocess.run([sys.executable, f"scripts/{script}"], check=True)
