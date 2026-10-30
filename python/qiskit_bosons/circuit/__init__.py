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
Circuit Representation
======================

.. currentmodule:: qiskit_bosons.circuit

With this module you can build :class:`~qiskit.circuit.QuantumCircuit` objects that are
interpreted as acting on `bosonic` modes rather than `qubits`.

The following objects are simple type aliases to make this re-interpretation more clear:

.. autosummary::
   :toctree: ../stubs/
   :template: autosummary/typealias.rst

   BosonicDAGCircuit
   BosonicMode
   BosonicRegister
   BosonicSpecifier

For actually expressing your circuits, you should use the :class:`.BosonicCircuit` class along with
the :mod:`~qiskit_bosons.circuit.library` of :class:`.BosonicGate` implementations.

.. autosummary::
   :toctree: ../stubs/

   BosonicCircuit
   BosonicGate
"""

from collections.abc import Sequence
from typing import TypeAlias

from qiskit.circuit import QuantumRegister, Qubit
from qiskit.dagcircuit import DAGCircuit

BosonicMode: TypeAlias = Qubit
"""A type alias of :class:`~qiskit.circuit.Qubit`.

Although this does not really give us any functional guarantees, it serves to better document the
API contract that instances of this type are interpreted as bosonic modes.
"""

BosonicRegister: TypeAlias = QuantumRegister
"""A type alias of :class:`~qiskit.circuit.QubitRegister`.

Although this does not really give us any functional guarantees, it serves to better document the
API contract that instances of this type are interpreted as registers of bosonic modes.
"""

BosonicSpecifier: TypeAlias = (
    BosonicMode | BosonicRegister | int | slice | Sequence[BosonicMode | int]
)
"""A type alias equivalent to Qiskit's ``QubitSpecifier`` but for bosonic modes."""

BosonicDAGCircuit: TypeAlias = DAGCircuit
"""A type alias of :class:`~qiskit.dagcircuit.DAGCircuit`.

Although this does not really give us any functional guarantees, it serves to better document the
API contract that instances of this type only contain circuit instructions of type
:class:`.BosonicGate` or :class:`~qiskit.circuit.library.Barrier`.
"""

# NOTE: we must explicitly define the type aliases _before_ the following imports to ensure that
# they can actually use those type aliases themselves.
# ruff: noqa: E402
from .bosonic_circuit import BosonicCircuit
from .bosonic_gate import BosonicGate

__all__ = [
    "BosonicCircuit",
    "BosonicDAGCircuit",
    "BosonicGate",
    "BosonicMode",
    "BosonicRegister",
    "BosonicSpecifier",
]
