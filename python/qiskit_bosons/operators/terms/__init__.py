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
==============
Operator Terms
==============

.. currentmodule:: qiskit_bosons.operators.terms

This module collects routines that operate on the individual terms of an operator, partitioning
them based on their structure. The submodules below group related functionality.

Grouping
--------

Terms may be tagged with group indices via
:attr:`~qiskit_bosons.operators.BosonOperator.groups`, which an
:class:`~qiskit_bosons.circuit.library.Evolution` gate then honors by decomposing into one
sub-evolution per group. Grouping terms *before* the boson-to-qubit mapping is what lets a grouping
chosen from the problem's structure survive into the circuit.

No automatic grouping strategy ships yet; assign the indices directly. What this module does
provide is the *analysis* of an existing grouping: :func:`group_coeff_means` reduces it to one
sampling weight per group, :func:`groups_are_hermitian` reports whether each group is separately
Hermitian (which a randomized product formula, and Qiskit's ``PauliEvolutionGate``, require), and
:func:`groups_have_uniform_coeffs` whether its coefficients are numerically equal.

.. autosummary::
   :toctree: ../stubs/

   group_coeff_means
   groups_are_hermitian
   groups_have_uniform_coeffs

Filtering
---------

This module provides convenience functions for removing terms from an operator that do not
contribute meaningfully to a downstream computation.

Members
^^^^^^^

.. autosummary::
   :toctree: ../stubs/

   filter_diagonal_terms

Ordering
--------

This module provides functions to reorder the terms of an operator. While the
resulting operators are mathematically identical, the term order can have
significant implications on algorithm behavior at runtime, due to the term
iteration order.

Members
^^^^^^^

.. autosummary::
   :toctree: ../stubs/

   canonical_order
   group_order
   order_terms

Group-ordered operators
^^^^^^^^^^^^^^^^^^^^^^^

:func:`group_order` makes each group one contiguous run of terms. That layout is what lets
:meth:`~qiskit_bosons.operators.OperatorTrait.split_out_groups` binary-search the requested groups
instead of scanning every term, so ordering once up front pays off before repeatedly transpiling with
:class:`~qiskit_bosons.circuit.library.synthesis.BosonicQDrift`, which looks up its sampled groups on
every call.
"""

from .filtering import filter_diagonal_terms
from .grouping import group_coeff_means, groups_are_hermitian, groups_have_uniform_coeffs
from .ordering import canonical_order, group_order, order_terms

__all__ = [
    "canonical_order",
    "filter_diagonal_terms",
    "group_coeff_means",
    "group_order",
    "groups_are_hermitian",
    "groups_have_uniform_coeffs",
    "order_terms",
]
