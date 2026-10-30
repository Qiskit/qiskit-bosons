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

"""OccupationB2QLayout tests.

Unlike the fermionic trivial layout, a boson-to-qubit layout must widen the register: a mode needs
more than one qubit. These tests pin the width per encoding and truncation.
"""

from __future__ import annotations

import unittest

import pytest
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import InitializeModes
from qiskit_bosons.mappers import LevelEncoding, LevelSpec
from qiskit_bosons.transpiler import BosonicCircuitToDAG
from qiskit_bosons.transpiler.passes import OccupationB2QLayout


def run_layout(num_modes: int, max_occupation: int, encoding: str):
    """Runs the layout pass over an ``num_modes``-mode circuit and returns its ``b2q_layout``."""
    circuit = BosonicCircuit(num_modes)
    circuit.append(InitializeModes.vacuum(num_modes), circuit.modes)
    dag = BosonicCircuitToDAG().run(circuit)

    pass_ = OccupationB2QLayout(LevelSpec(max_occupation, encoding))
    pass_.property_set = {}
    returned = pass_.run(dag)

    assert returned is dag, "this is an analysis pass and must not modify the DAG"
    return pass_, pass_.property_set["b2q_layout"]


def num_qubits_per_mode(max_occupation: int, encoding: str) -> int:
    """The number of qubits one mode occupies, written out independently of the layout."""
    levels = max_occupation + 1
    if encoding in ("unary_projected", "unary_unprojected"):
        return levels
    return (levels - 1).bit_length()


class TestNumQubitsPerMode(unittest.TestCase):
    def test_matches_the_mapper(self):
        for encoding in ("unary_projected", "standard_binary", "gray"):
            for max_occupation in (1, 2, 4, 7):
                spec = LevelSpec(max_occupation, encoding)
                assert spec.resolve(1).widths == [num_qubits_per_mode(max_occupation, encoding)]

    def test_stores_its_arguments(self):
        pass_ = OccupationB2QLayout(LevelSpec(5, "gray"))
        assert pass_.spec.max_occupation == 5
        assert pass_.spec.encoding == LevelEncoding.GRAY


class TestRegisterWidth(unittest.TestCase):
    def test_unary_uses_d_qubits_per_mode(self):
        _, layout = run_layout(3, 2, "unary_projected")
        assert [len(qreg) for qreg in layout.values()] == [9]

    def test_compact_encodings_use_the_logarithm(self):
        for encoding in ("standard_binary", "gray"):
            _, layout = run_layout(3, 2, encoding)
            assert [len(qreg) for qreg in layout.values()] == [6], encoding

    def test_width_scales_with_the_truncation(self):
        # d = 5 needs 3 qubits per mode under a compact encoding
        _, layout = run_layout(3, 4, "gray")
        assert [len(qreg) for qreg in layout.values()] == [9]

    def test_width_is_num_modes_times_per_mode(self):
        for encoding in ("unary_projected", "standard_binary", "gray"):
            for num_modes in (1, 2, 4):
                for max_occupation in (1, 2, 4):
                    _, layout = run_layout(num_modes, max_occupation, encoding)
                    (qreg,) = layout.values()
                    assert len(qreg) == num_modes * num_qubits_per_mode(max_occupation, encoding)

    def test_layout_is_keyed_by_the_mode_register(self):
        circuit = BosonicCircuit(2)
        circuit.append(InitializeModes.vacuum(2), circuit.modes)
        dag = BosonicCircuitToDAG().run(circuit)
        pass_ = OccupationB2QLayout(LevelSpec(2, "gray"))
        pass_.property_set = {}
        pass_.run(dag)

        layout = pass_.property_set["b2q_layout"]
        assert set(layout) == set(dag.qregs.values())

    def test_qubit_register_is_wider_than_the_mode_register(self):
        """This is exactly why the trivial fermionic layout is unusable here."""
        for encoding in ("unary_projected", "standard_binary", "gray"):
            _, layout = run_layout(3, 2, encoding)
            for breg, qreg in layout.items():
                assert len(qreg) > len(breg), encoding

    def test_one_pass_serves_circuits_of_different_widths(self):
        """The registers are derived from the circuit, so the same pass handles any width."""
        pass_ = OccupationB2QLayout(LevelSpec(2, "gray"))
        widths = []
        for num_modes in (1, 3, 5):
            circuit = BosonicCircuit(num_modes)
            circuit.append(InitializeModes.vacuum(num_modes), circuit.modes)
            pass_.property_set = {}
            pass_.run(BosonicCircuitToDAG().run(circuit))
            (qreg,) = pass_.property_set["b2q_layout"].values()
            widths.append(len(qreg))
        assert widths == [2, 6, 10]


class TestErrorPaths(unittest.TestCase):
    def test_unknown_encoding_raises(self):
        # The name is resolved once, when the spec is built, so the error surfaces there.
        with pytest.raises(ValueError):
            OccupationB2QLayout(LevelSpec(2, "not_an_encoding"))
