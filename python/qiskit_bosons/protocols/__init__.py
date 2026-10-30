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
=========
Protocols
=========

.. currentmodule:: qiskit_bosons.protocols

This module provides the protocols used throughout the package.

A :class:`typing.Protocol` is an interface that a class can fulfill structurally, without explicitly
inheriting from it. For example, any class that implements a ``__len__`` method satisfies
:class:`~typing.Sized`, whether or not it says so anywhere in its class hierarchy.

The protocols here follow the same idea, each centered on a small set of dunder-style methods. A
class implementing them thereby satisfies the protocol and becomes usable wherever it is expected,
without any explicit registration.

.. list-table::
   :header-rows: 1

   * - Protocol
     - Methods
     - Helper functions
   * - :class:`.SupportsCommutators`
     - ``_commutator_``, ``_anti_commutator_``, ``_double_commutator_``
     - :func:`.commutator`, :func:`.anti_commutator`, :func:`.double_commutator`

.. autosummary::
   :toctree: ../stubs/
   :class: sd-d-none

   SupportsCommutators
"""

from .commutators import SupportsCommutators

__all__ = [
    "SupportsCommutators",
]
