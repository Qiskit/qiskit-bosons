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

"""OccupationInitializeModesSynthesis tests.

The gate is encoding-agnostic; the plugin is where the truncation and the encoding enter. The most
important behavior here is the transpile-time check that the declared occupation actually fits the
configured cutoff -- there is no such check on the gate itself, and no fermionic analogue of it.
"""

from __future__ import annotations

import unittest

import numpy as np
import pytest
from qiskit.quantum_info import Statevector
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import InitializeModes
from qiskit_bosons.mappers import LevelEncoding, LevelSpec
from qiskit_bosons.transpiler.passes import OccupationInitializeModesSynthesis
from qiskit_bosons.transpiler.presets import generate_preset_occupation_pass_manager

ENCODINGS = ("unary_projected", "standard_binary", "gray")


def level_code(level: int, encoding: str) -> int:
    """The bit pattern an occupation level is stored as, written out independently of the plugin."""
    match encoding:
        case "unary_projected" | "unary_unprojected":
            return 1 << level
        case "standard_binary":
            return level
        case "gray":
            return level ^ (level >> 1)
        case _:  # pragma: no cover
            raise AssertionError(encoding)


def expected_index(occupation, encoding: str, max_occupation: int) -> int:
    """The qubit basis index the prepared Fock state must land on."""
    per_mode = LevelSpec(max_occupation, encoding).resolve(1).total_width
    index = 0
    for mode, level in enumerate(occupation):
        index |= level_code(level, encoding) << (mode * per_mode)
    return index


def lower(occupation, encoding: str, max_occupation: int):
    """Transpiles a bare ``InitializeModes`` under the occupation preset."""
    num_modes = len(occupation)
    circuit = BosonicCircuit(num_modes)
    circuit.append(InitializeModes(occupation), circuit.modes)
    pass_manager = generate_preset_occupation_pass_manager(LevelSpec(max_occupation, encoding))
    return pass_manager.run(circuit)


class TestConstruction(unittest.TestCase):
    def test_stores_its_arguments(self):
        spec = LevelSpec(4, "gray")
        plugin = OccupationInitializeModesSynthesis(spec)
        assert plugin.spec is spec
        assert plugin.spec.max_occupation == 4
        assert plugin.spec.encoding == LevelEncoding.GRAY


class TestPreparedState(unittest.TestCase):
    """The Fock state must land on exactly the basis index the encoding assigns to it."""

    def test_vacuum_needs_no_gates_for_compact_encodings(self):
        for encoding in ("standard_binary", "gray"):
            circuit = lower([0, 0], encoding, 2)
            assert "x" not in circuit.count_ops(), encoding

    def test_vacuum_still_needs_one_gate_per_mode_for_unary(self):
        """Unary stores level 0 as the one-hot pattern on qubit 0, not as all-zeros."""
        circuit = lower([0, 0], "unary_projected", 2)
        assert circuit.count_ops()["x"] == 2

    def test_prepares_the_declared_fock_state(self):
        for encoding in ENCODINGS:
            for occupation in ([0, 0], [1, 0], [0, 2], [2, 1], [1, 1]):
                max_occupation = 2
                circuit = lower(occupation, encoding, max_occupation)
                amplitudes = np.asarray(Statevector(circuit).data)
                index = expected_index(occupation, encoding, max_occupation)
                assert abs(amplitudes[index]) == pytest.approx(1.0), (encoding, occupation)

    def test_a_single_mode_at_every_level(self):
        for encoding in ENCODINGS:
            max_occupation = 3
            for level in range(max_occupation + 1):
                circuit = lower([level], encoding, max_occupation)
                amplitudes = np.asarray(Statevector(circuit).data)
                index = expected_index([level], encoding, max_occupation)
                assert abs(amplitudes[index]) == pytest.approx(1.0), (encoding, level)

    def test_only_x_gates_are_emitted(self):
        """A product Fock state from the all-zero state needs no entangling gate at all."""
        for encoding in ENCODINGS:
            circuit = lower([2, 1, 0], encoding, 2)
            assert set(circuit.count_ops()) <= {"x"}

    def test_gate_count_matches_the_code_popcount(self):
        for encoding in ENCODINGS:
            occupation = [3, 1, 2]
            circuit = lower(occupation, encoding, 3)
            expected = sum(level_code(level, encoding).bit_count() for level in occupation)
            assert circuit.count_ops().get("x", 0) == expected, encoding

    def test_register_width(self):
        for encoding in ENCODINGS:
            circuit = lower([1, 0, 2], encoding, 2)
            assert circuit.num_qubits == LevelSpec(2, encoding).resolve(3).total_width

    def test_a_larger_cutoff_widens_but_does_not_change_the_state_labels(self):
        """Raising the cutoff pads with unoccupied levels; the state is still |1, 0>."""
        for encoding in ("standard_binary", "gray"):
            narrow = np.asarray(Statevector(lower([1, 0], encoding, 1)).data)
            wide = np.asarray(Statevector(lower([1, 0], encoding, 3)).data)
            assert abs(narrow[expected_index([1, 0], encoding, 1)]) == pytest.approx(1.0)
            assert abs(wide[expected_index([1, 0], encoding, 3)]) == pytest.approx(1.0)


class TestTruncationCheck(unittest.TestCase):
    """The transpile-time guard, which has no fermionic counterpart."""

    def test_occupation_above_the_cutoff_raises(self):
        with pytest.raises(ValueError, match="exceeds the truncation"):
            lower([3, 0], "unary_projected", 1)

    def test_error_names_both_numbers(self):
        with pytest.raises(ValueError, match=r"occupation of 5.*max_occupation=2"):
            lower([5], "gray", 2)

    def test_occupation_exactly_at_the_cutoff_is_accepted(self):
        for encoding in ENCODINGS:
            circuit = lower([2, 2], encoding, 2)
            amplitudes = np.asarray(Statevector(circuit).data)
            assert abs(amplitudes[expected_index([2, 2], encoding, 2)]) == pytest.approx(1.0)

    def test_only_the_largest_entry_matters(self):
        """A large *total* boson number is fine as long as no single mode exceeds the cutoff."""
        circuit = lower([2, 2, 2], "gray", 2)
        assert sum([2, 2, 2]) > 2
        amplitudes = np.asarray(Statevector(circuit).data)
        assert abs(amplitudes[expected_index([2, 2, 2], "gray", 2)]) == pytest.approx(1.0)

    def test_the_check_fires_for_every_encoding(self):
        for encoding in ENCODINGS:
            with pytest.raises(ValueError, match="exceeds the truncation"):
                lower([4], encoding, 2)

    def test_vacuum_is_fine_at_the_smallest_cutoff(self):
        # A single-level mode needs no qubit under a compact encoding, so a spec made only of such
        # modes is rejected (see `test_a_spec_without_qubits_is_rejected`); the unary encoding
        # still gives it one.
        circuit = lower([0, 0], "unary_projected", 0)
        assert circuit.num_qubits == 2
        assert circuit.count_ops()["x"] == 2


class TestLocalPlacement(unittest.TestCase):
    def test_a_gate_on_a_subset_of_modes(self):
        """A gate placed on modes 1..2 of a 4-mode register must set only those modes' qubits."""
        encoding, max_occupation = "gray", 2
        per_mode = LevelSpec(max_occupation, encoding).resolve(1).total_width

        circuit = BosonicCircuit(4)
        circuit.append(InitializeModes([1, 2]), circuit.modes[1:3])
        lowered = generate_preset_occupation_pass_manager(LevelSpec(max_occupation, encoding)).run(
            circuit
        )

        index = 0
        for mode, level in ((1, 1), (2, 2)):
            index |= level_code(level, encoding) << (mode * per_mode)
        amplitudes = np.asarray(Statevector(lowered).data)
        assert abs(amplitudes[index]) == pytest.approx(1.0)

    def test_the_truncation_check_applies_to_a_local_gate_too(self):
        circuit = BosonicCircuit(4)
        circuit.append(InitializeModes([5]), circuit.modes[2:3])
        with pytest.raises(ValueError, match="exceeds the truncation"):
            generate_preset_occupation_pass_manager(LevelSpec(2, "gray")).run(circuit)
