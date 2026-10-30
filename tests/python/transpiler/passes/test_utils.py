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

"""Tests for the transpiler pass utility helpers."""

from __future__ import annotations

import pytest
from qiskit.circuit import QuantumRegister
from qiskit.circuit.library import XXPlusYYGate
from qiskit.dagcircuit import DAGOpNode
from qiskit_bosons.transpiler.passes.utils import map_node_single_register


def test_map_node_single_register_maps_indices():
    """A node acting on one bosonic register maps to global indices and the qubit register."""
    breg = QuantumRegister(4, "b")
    qreg = QuantumRegister(8, "q")
    b2q_layout = {breg: qreg}

    # A two-mode gate acting on modes 1 and 3 of the single register.
    node = DAGOpNode(XXPlusYYGate(0.5), qargs=(breg[1], breg[3]), cargs=())

    breg_indices, mapped_qreg = map_node_single_register(node, b2q_layout)
    assert breg_indices == [1, 3]
    assert mapped_qreg is qreg


def test_map_node_single_register_rejects_multiple_registers():
    """A node whose modes span two bosonic registers is unsupported."""
    breg_a = QuantumRegister(2, "a")
    breg_b = QuantumRegister(2, "b")
    b2q_layout = {breg_a: QuantumRegister(2), breg_b: QuantumRegister(2)}

    # A two-mode gate straddling the two registers.
    node = DAGOpNode(XXPlusYYGate(0.5), qargs=(breg_a[0], breg_b[1]), cargs=())

    with pytest.raises(NotImplementedError, match="spread across multiple BosonicRegister"):
        map_node_single_register(node, b2q_layout)
