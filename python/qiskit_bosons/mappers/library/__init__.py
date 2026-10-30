# This code is a Qiskit project.
#
# (C) Copyright IBM 2026.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at https://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# ruff: noqa: D205,D212,D415
r"""
==============
Mapper Library
==============

.. currentmodule:: qiskit_bosons.mappers.library

This module provides the boson-to-qubit mappers.

:func:`occupation_encoding` maps a :class:`~qiskit_bosons.operators.BosonOperator` in the occupation
(Fock) basis, truncating each mode and laying its levels out on qubits as a
:class:`~qiskit_bosons.mappers.LevelSpec` describes.

A bosonic encoding is **not** an algebra isomorphism: a bosonic mode's Fock space is
infinite-dimensional and must be truncated first, which breaks the canonical commutation relation on
the highest retained level. See :func:`occupation_encoding` for the exact statement and its
consequences.

.. autosummary::
   :toctree: ../stubs/

   occupation_encoding
"""

from __future__ import annotations

from qiskit_bosons._lib.mappers.mappers_library.occupation import occupation_encoding

__all__ = [
    "occupation_encoding",
]
