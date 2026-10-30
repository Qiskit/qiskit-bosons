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
=========================
Transpiler Pass Utilities
=========================

.. currentmodule:: qiskit_bosons.transpiler.passes.utils

Several convenience functions are provided to simplify the implementation of common operations
during the implementation of custom transpiler passes.

.. autosummary::
   :toctree: ../stubs/

   map_node_single_register
"""

from __future__ import annotations

from qiskit.circuit import QuantumRegister
from qiskit.dagcircuit import DAGOpNode

from .. import B2QLayout


def map_node_single_register(
    in_node: DAGOpNode, b2q_layout: B2QLayout
) -> tuple[list[int], QuantumRegister]:
    """Map the single bosonic mode register acted upon by the circuit operation.

    Given a :external:class:`~qiskit.dagcircuit.DAGOpNode` and a
    :type:`~qiskit_bosons.transpiler.B2QLayout` this function will:

    1. ensure that all acted-upon bosonic modes belong to the same register
    2. map the node's indices to their global indices with respect to said bosonic mode register
    3. map this bosonic mode register to the corresponding qubit register

    .. caution::
       The returned global bosonic mode register indices are knowingly dependent on the order of
       the bosonic modes during the construction of their original
       :type:`~qiskit_bosons.circuit.BosonicRegister`. This is unlikely to cause any issues but
       makes the usage of :class:`.RelabelModes` possibly more generally.

    Args:
        in_node: the node whose indices to parse.
        b2q_layout: the mapping of bosonic to qubit registers.

    Returns:
        The list of the global bosonic mode indices and mapped qubit register instance.

    Raises:
        NotImplementedError: when the provided node acts on bosonic modes that do not all belong
            to the same :class:`.BosonicRegister`.
    """
    encountered_bosonic_registers: set[QuantumRegister] = set()
    breg_indices = []
    for boson in in_node.qargs:
        for breg in b2q_layout:
            if boson in breg:
                encountered_bosonic_registers.add(breg)
                # HACK: we explicitly use the private _index attribute here to enable the general
                # usage of the RelabelModes pass. This _should_ have no adversarial effect on other
                # scenarios unless the user tampered with the BosonicRegister themselves, too.
                breg_indices.append(boson._index)
                break

    if len(encountered_bosonic_registers) > 1:
        raise NotImplementedError(
            "Cannot map a BosonicGate acting on bosonic modes that are spread across "
            "multiple BosonicRegister instances."
        )

    breg = encountered_bosonic_registers.pop()
    qreg = b2q_layout[breg]

    return breg_indices, qreg
