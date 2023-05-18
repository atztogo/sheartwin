"""Phonon calculation under shaer for Ti.

Based on

Atsushi Togo, Yuta Inoue, and Isao Tanaka
Phys. Rev. B 102, 024106 (2020)

"""
import os
import numpy as np
from typing import Optional, Union
import phonopy
from phonopy.interface.phonopy_yaml import read_cell_yaml
from phonopy.structure.atoms import PhonopyAtoms
from collections.abc import Sequence, Callable


class ShearTwin:
    """Class to generate sheared unit cells."""

    def __init__(
        self,
        unitcell: Optional[PhonopyAtoms] = None,
        cell_yaml_filename: Optional[Union[str, bytes, os.PathLike]] = None,
        Q: Optional[Union[Sequence, np.ndarray]] = None,
        T: Optional[Callable] = None,
    ):
        """Init method."""
        self._Q = Q
        self._T = T
        self._unitcell = unitcell
        if cell_yaml_filename is not None:
            self._set_cell_from_yaml(cell_yaml_filename)

    def run(self, t: float) -> PhonopyAtoms:
        """Return deformed unit cell."""
        if self._unitcell is None:
            raise RuntimeError("unitcell is not set.")
        if self._Q is None:
            raise RuntimeError("Q is not set.")
        if self._T is None:
            raise RuntimeError("T is not set.")

        cell = self._unitcell.copy()
        cell.cell = np.dot(cell.cell.T, self._get_QTQinv(t)).T
        return cell

    @property
    def Q(self):
        """Setter and gettter of Q matrix."""
        return self._Q

    @Q.setter
    def Q(self, matrix: Union[Sequence, np.ndarray]):
        self._Q = matrix

    @property
    def T(self) -> Optional[Callable]:
        """Setter and gettter of T matrix function."""
        return self._T

    @T.setter
    def T(self, matrix_function: Callable):
        self._T = matrix_function

    @property
    def unitcell(self) -> Optional[PhonopyAtoms]:
        """Return unit cell."""
        return self._unitcell

    def _set_cell_from_yaml(
        self, cell_yaml_filename: Union[str, bytes, os.PathLike]
    ) -> phonopy.Phonopy:
        """Set unit cell from cell yaml file."""
        self._unitcell = read_cell_yaml(cell_yaml_filename)

    def _get_QTQinv(self, t: float) -> np.ndarray:
        Q = self._Q
        T = self._T
        return np.dot(Q, np.dot(T(t), np.linalg.inv(Q)))
