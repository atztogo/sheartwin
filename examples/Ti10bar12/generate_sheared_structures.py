"""Generate sheared structures of LAMMPS for Ti 10-12 twin."""
import numpy as np
from sheartwin import ShearTwin
import spglib
from phonopy.structure.cells import get_cell_parameters
from phonopy.interface.calculator import write_crystal_structure

Q_10bar12 = [[2, 0, -2], [1, 1, -1], [1, 0, 1]]


def T_10bar12(t):
    """Return change of basis matrix for {10-12}."""
    return [[1, 0, t], [0, 1, 0], [0, 0, 1]]


sheartwin = ShearTwin(cell_yaml_filename="phonopy_symcells_Ti.yaml")
sheartwin.Q = Q_10bar12
sheartwin.T = T_10bar12
cell = sheartwin.unitcell
Q = sheartwin.Q
lattice = np.dot(cell.cell.T, Q)  # (a', b', c') as column vectors of extended unit cell
a, b, c = get_cell_parameters(lattice.T)
tt = -2 * np.dot(lattice[:, 2], lattice[:, 0]) / a**2  # twinning shear -2a.c/|a|^2
print(tt)
sheared_cells = []
for i in range(21):
    t = tt * i / 20.0
    sheared_cells.append(sheartwin.run(t))
print(spglib.get_spacegroup(sheared_cells[0]))
print(sheared_cells[0])
print(spglib.get_spacegroup(sheared_cells[-1]))
print(sheared_cells[-1])
print()

# Lattices of extend unit cells are same although orientations are different.
abc_matrix = get_cell_parameters(np.dot(sheared_cells[0].cell.T, Q).T)
abc_twin = get_cell_parameters(np.dot(sheared_cells[-1].cell.T, Q).T)
np.testing.assert_allclose(abc_matrix, abc_twin)
for i, sh_cell in enumerate(sheared_cells):
    write_crystal_structure(f"structure_{i}", sh_cell, interface_mode="lammps")
