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
"""
======================
Representation Mappers
======================

.. currentmodule:: qiskit_bosons.mappers

This module provides a framework for implementing custom boson-to-qubit mapper routines.

Level encodings
===============

A :class:`.LevelEncoding` decides how the :math:`d` levels of a single system are laid out on qubits.
It says nothing about what the levels are, so it is shared by every mapper that works in a basis of
discrete levels: the occupation-number basis of a truncated bosonic mode, vibrational modals, or the
points of a position grid. The mappers in :mod:`~qiskit_bosons.mappers.library` are given one as
part of a :class:`.LevelSpec` (see below); for example
:func:`~qiskit_bosons.mappers.library.occupation_encoding` maps a
:class:`~qiskit_bosons.operators.BosonOperator` in the occupation basis.

.. autosummary::
   :toctree: ../stubs/

   LevelEncoding

A :class:`.LevelSpec` combines a level encoding with the truncation of every mode of a register. It
holds one cutoff per mode (or one cutoff shared by all modes) and from them derives how many qubits
each mode occupies and where its qubits start, so everything that has to agree about that is handed
the same spec.

.. autosummary::
   :toctree: ../stubs/

   LevelSpec

Generator-based mappers
=======================

The :mod:`~qiskit_bosons.mappers.library` submodule ships the ready-made mappers. The generator
function below is the extension point for writing your own: supply the qubit operator that a single
:class:`~qiskit_bosons.operators.boson_action.BosonAction` maps to, and it lifts that into a mapper
over whole operators.

.. autosummary::
   :toctree: ../stubs/

   map_boson_action_generators
"""

from qiskit_bosons._lib.mappers.levels import LevelEncoding, LevelSpec

from .boson_generators import map_boson_action_generators

__all__ = [
    "LevelEncoding",
    "LevelSpec",
    "map_boson_action_generators",
]
