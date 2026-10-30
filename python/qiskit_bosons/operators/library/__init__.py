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
================
Operator Library
================

.. currentmodule:: qiskit_bosons.operators.library

Building on top of the operator representations provided by :mod:`~qiskit_bosons.operators`, this
module provides tools to simplify the construction of commonly used operators.

Model Hamiltonians are built directly from your own lattice rather than from helpers here; the
:ref:`Bose-Hubbard guide <bose_hubbard_getting_started>` works one through end to end.

Commutator Generators
---------------------

The brackets through which the bosonic algebra is expressed. The commutator is the defining one:
:math:`[b_i, b_j^\dagger] = \delta_{ij}` *is* the bosonic algebra, so these are how that structure
gets interrogated in practice (checking whether two terms commute, whether an operator conserves
boson number, and so on).

.. autosummary::
   :toctree: ../stubs/

   commutator
   anti_commutator
   double_commutator

Single-mode building blocks
---------------------------

The dimensionless quadrature operators, which are the pieces most other bosonic operators are
written in terms of and whose normalization is easy to get wrong.

These operators are Hermitian and satisfy the commutation relation:

.. math::

       [\hat{q}_i, \hat{p}_j] = i\delta_{ij}

.. doctest::

    >>> from qiskit_bosons.operators import BosonOperator
    >>> from qiskit_bosons.operators.library import commutator, momentum, position
    >>> commutator(position(0), momentum(0)).equiv(1j * BosonOperator.one())
    True
    >>> commutator(position(0), momentum(1)).equiv(BosonOperator.zero())
    True

.. autosummary::
   :toctree: ../stubs/

   position
   momentum
"""

from .commutators import anti_commutator, commutator, double_commutator
from .quadratures import momentum, position

__all__ = [
    "anti_commutator",
    "commutator",
    "double_commutator",
    "momentum",
    "position",
]
