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

"""Bosonic mode rotation synthesis tests.

A ``ModeRotation`` is broken down in bosonic space into atomic, single-group ``Evolution`` gates.
These tests pin that structure and that its product reproduces the rotation on the single-particle
sector, independently of any encoding.
"""

from __future__ import annotations

import unittest

import numpy as np
import pytest
import scipy.linalg
from qiskit.quantum_info import random_unitary
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, ModeRotation
from qiskit_bosons.circuit.library.synthesis import (
    B2BSynthesisPluginManager,
    BosonicClementsDecomposition,
)

CLEMENTS = BosonicClementsDecomposition()


def single_particle_block(op, num_modes: int) -> np.ndarray:
    """Reads the matrix ``K_ij`` back out of ``sum_ij K_ij b^dagger_i b_j``."""
    block = np.zeros((num_modes, num_modes), dtype=complex)
    for term, coeff in op.iter_terms():
        (_, row), (_, col) = term
        block[row, col] += coeff
    return block


def embedded_unitary(circuit: BosonicCircuit, num_modes: int) -> np.ndarray:
    """The single-particle unitary of a circuit of ``Evolution`` gates, exp(-i H) in order."""
    total = np.eye(num_modes, dtype=complex)
    for instruction in circuit._inner.data:
        gate = instruction.operation
        modes = [circuit.register.index(mode) for mode in instruction.qubits]
        local = single_particle_block(gate.operator, len(modes))
        factor = np.eye(num_modes, dtype=complex)
        factor[np.ix_(modes, modes)] = scipy.linalg.expm(-1j * gate.params[0] * local)
        total = factor @ total
    return total


class TestBosonicClementsDecomposition(unittest.TestCase):
    def test_is_registered_as_an_entry_point(self):
        manager = B2BSynthesisPluginManager()
        assert manager.method("ModeRotation", "Clements") is BosonicClementsDecomposition

    def test_a_gate_has_no_decomposition_by_default(self):
        assert ModeRotation(np.eye(2)).synthesis is None

    def test_a_gate_without_a_synthesis_is_left_whole_by_decompose(self):
        circuit = BosonicCircuit(3)
        circuit.append(ModeRotation(random_unitary(3, seed=7).data), circuit.modes)
        assert dict(circuit.decompose().count_ops()) == {"ModeRotation": 1}

    def test_a_generic_unitary_gives_a_mesh_of_evolutions_and_phases(self):
        """``M(M-1)/2`` two-mode rotations and ``M`` phases."""
        gate = ModeRotation(random_unitary(3, seed=7).data, synthesis=CLEMENTS)
        circuit = BosonicCircuit(3)
        circuit.append(gate, circuit.modes)
        decomposed = circuit.decompose()
        assert dict(decomposed.count_ops()) == {"Evolution": 6}
        widths = sorted(instruction.operation.num_modes for instruction in decomposed._inner.data)
        assert widths == [1, 1, 1, 2, 2, 2]

    def test_every_factor_is_atomic_and_a_single_group(self):
        gate = ModeRotation(random_unitary(4, seed=2).data, synthesis=CLEMENTS)
        circuit = BosonicCircuit(4)
        circuit.append(gate, circuit.modes)
        for instruction in circuit.decompose()._inner.data:
            factor = instruction.operation
            assert isinstance(factor, Evolution)
            assert factor.atomic
            assert factor.operator.num_groups() == 1
            assert factor.params[0] == 1.0

    def test_every_factor_is_hermitian(self):
        gate = ModeRotation(random_unitary(4, seed=3).data, synthesis=CLEMENTS)
        circuit = BosonicCircuit(4)
        circuit.append(gate, circuit.modes)
        for instruction in circuit.decompose()._inner.data:
            assert instruction.operation.operator.is_hermitian(0.0)

    def test_atomic_factors_are_not_decomposed_further(self):
        gate = ModeRotation(random_unitary(3, seed=7).data, synthesis=CLEMENTS)
        circuit = BosonicCircuit(3)
        circuit.append(gate, circuit.modes)
        once = circuit.decompose()
        assert dict(once.decompose().count_ops()) == dict(once.count_ops())

    def test_the_product_of_the_factors_is_the_rotation(self):
        """On the single-particle sector, which the homomorphism makes exact."""
        for num_modes, seed in ((2, 1), (3, 2), (4, 3), (5, 4)):
            unitary = random_unitary(num_modes, seed=seed).data
            circuit = BosonicCircuit(num_modes)
            circuit.append(ModeRotation(unitary, synthesis=CLEMENTS), circuit.modes)
            rebuilt = embedded_unitary(circuit.decompose(), num_modes)
            np.testing.assert_allclose(rebuilt, unitary, atol=1e-9)

    def test_a_beamsplitter_is_a_single_two_mode_factor(self):
        circuit = BosonicCircuit(2)
        circuit.append(
            ModeRotation(ModeRotation.beamsplitter(np.pi / 4).rotation_unitary, synthesis=CLEMENTS),
            circuit.modes,
        )
        decomposed = circuit.decompose()
        assert dict(decomposed.count_ops()) == {"Evolution": 1}
        assert next(iter(decomposed._inner.data)).operation.num_modes == 2

    def test_the_identity_decomposes_to_nothing(self):
        circuit = BosonicCircuit(3)
        circuit.append(ModeRotation(np.eye(3), synthesis=CLEMENTS), circuit.modes)
        assert dict(circuit.decompose().count_ops()) == {}

    def test_a_pure_phase_shift_gives_only_single_mode_factors(self):
        circuit = BosonicCircuit(3)
        circuit.append(
            ModeRotation(
                ModeRotation.phase_shift([0.3, 0.0, -0.4]).rotation_unitary, synthesis=CLEMENTS
            ),
            circuit.modes,
        )
        decomposed = circuit.decompose()
        assert dict(decomposed.count_ops()) == {"Evolution": 2}
        assert all(i.operation.num_modes == 1 for i in decomposed._inner.data)

    def test_a_loose_tolerance_drops_small_phases(self):
        circuit = BosonicCircuit(2)
        circuit.append(
            ModeRotation(
                np.diag(np.exp(1j * np.array([1e-6, 0.5]))),
                synthesis=BosonicClementsDecomposition(tol=1e-3),
            ),
            circuit.modes,
        )
        assert dict(circuit.decompose().count_ops()) == {"Evolution": 1}


class TestModeRotationGate(unittest.TestCase):
    def test_equal_unitaries_give_equal_gates_and_different_ones_do_not(self):
        first = ModeRotation(random_unitary(2, seed=1).data)
        assert first == ModeRotation(first.rotation_unitary.copy())
        assert first != ModeRotation(random_unitary(2, seed=2).data)

    def test_equal_gates_hash_equally(self):
        unitary = random_unitary(3, seed=1).data
        assert hash(ModeRotation(unitary)) == hash(ModeRotation(unitary.copy()))

    def test_the_inverse_is_the_adjoint_rotation(self):
        gate = ModeRotation(random_unitary(3, seed=4).data)
        np.testing.assert_allclose(gate.inverse().rotation_unitary, gate.rotation_unitary.conj().T)
        assert isinstance(gate.inverse(), ModeRotation)

    def test_the_inverse_keeps_the_synthesis(self):
        synthesis = BosonicClementsDecomposition(tol=1e-9)
        gate = ModeRotation(random_unitary(3, seed=4).data, synthesis=synthesis)
        assert gate.inverse().synthesis is synthesis

    def test_the_synthesis_is_read_only(self):
        gate = ModeRotation(np.eye(2))
        with pytest.raises(AttributeError):
            gate.synthesis = BosonicClementsDecomposition()  # type: ignore[misc]
