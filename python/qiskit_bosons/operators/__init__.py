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
========================
Operator Representations
========================

.. currentmodule:: qiskit_bosons.operators

This module provides the data structures for representing bosonic operators.

.. autosummary::
   :toctree: ../stubs/

   OperatorTrait

Boson Operator
--------------

This operator represents bosonic operators in terms of the second-quantization creation and
annihilation operators, which obey the canonical commutation relations

.. math::

    \left[b_i, b^\dagger_j\right] = \delta_{ij}, \qquad
    \left[b_i, b_j\right] = \left[b^\dagger_i, b^\dagger_j\right] = 0 \, .

.. autosummary::
   :toctree: ../stubs/
   :template: autosummary/typealias.rst

   ~boson_action.BosonAction

.. autosummary::
   :toctree: ../stubs/

   BosonOperator
   cre
   ann

.. note::
   A bosonic mode has an **infinite-dimensional** Fock space. :class:`BosonOperator` is exact and
   carries no occupation cutoff; a truncation is introduced only when mapping to qubits (see
   :mod:`qiskit_bosons.mappers`). Two consequences are visible directly in this API: a mode index may
   repeat within a single term, and :meth:`~BosonOperator.max_occupation` and
   :meth:`~BosonOperator.mode_occupations` report the cutoff such a term demands.

Operator Library
----------------

A library of the model Hamiltonians that recur in bosonic simulation is provided in
:mod:`qiskit_bosons.operators.library`.
"""

from __future__ import annotations

from qiskit_bosons._lib.operators.boson_operator import BosonOperator

from .boson_action import BosonAction, ann, cre
from .operator_trait import OperatorTrait

__all__ = [
    "BosonAction",
    "BosonOperator",
    "OperatorTrait",
    "ann",
    "cre",
]
