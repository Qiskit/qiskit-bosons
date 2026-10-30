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

"""ModeRotation tests.

A ``ModeRotation`` is a passive linear-optical transformation. These tests pin the matrix
conventions of its constructors and the round trip between the unitary, its anti-Hermitian
generator, and the Hermitian :class:`~qiskit_bosons.operators.BosonOperator` that a synthesis
exponentiates.
"""

from __future__ import annotations

import unittest

import numpy as np
import scipy.linalg
from qiskit_bosons.circuit import BosonicGate
from qiskit_bosons.circuit.library import ModeRotation
from qiskit_bosons.operators import BosonOperator

RNG = np.random.default_rng(1234)


def random_unitary(dim: int) -> np.ndarray:
    """A Haar-ish random unitary from a seeded QR decomposition."""
    matrix = RNG.normal(size=(dim, dim)) + 1j * RNG.normal(size=(dim, dim))
    unitary, upper = np.linalg.qr(matrix)
    return unitary * (np.diag(upper) / np.abs(np.diag(upper)))


def single_particle_block(op: BosonOperator, num_modes: int) -> np.ndarray:
    """Reads the matrix ``K_ij`` back out of ``sum_ij K_ij b^dagger_i b_j``."""
    block = np.zeros((num_modes, num_modes), dtype=complex)
    for term, coeff in op.iter_terms():
        assert len(term) == 2, "a mode-rotation generator is quadratic"
        (raise_action, row), (lower_action, col) = term
        assert raise_action and not lower_action
        block[row, col] += coeff
    return block


class TestBeamsplitter(unittest.TestCase):
    def test_is_unitary(self):
        for theta in (0.0, 0.3, np.pi / 4, np.pi / 2, 1.7):
            for phi in (0.0, 0.5, np.pi):
                unitary = ModeRotation.beamsplitter(theta, phi).rotation_unitary
                np.testing.assert_allclose(unitary @ unitary.conj().T, np.eye(2), atol=1e-12)

    def test_clements_transfer_matrix(self):
        theta, phi = 0.7, 0.4
        cos, sin = np.cos(theta), np.sin(theta)
        phase = np.exp(1j * phi)
        expected = np.array([[phase * cos, -sin], [phase * sin, cos]], dtype=complex)
        np.testing.assert_allclose(
            ModeRotation.beamsplitter(theta, phi).rotation_unitary, expected, atol=1e-12
        )

    def test_balanced_at_quarter_pi(self):
        """``theta = pi/4`` is the 50:50 beamsplitter: every entry has magnitude ``1/sqrt(2)``."""
        unitary = ModeRotation.beamsplitter(np.pi / 4).rotation_unitary
        np.testing.assert_allclose(np.abs(unitary), np.full((2, 2), 2**-0.5), atol=1e-12)

    def test_zero_angle_is_the_identity(self):
        np.testing.assert_allclose(
            ModeRotation.beamsplitter(0.0).rotation_unitary, np.eye(2), atol=1e-12
        )

    def test_default_phase_is_zero(self):
        np.testing.assert_allclose(
            ModeRotation.beamsplitter(0.6).rotation_unitary,
            ModeRotation.beamsplitter(0.6, 0.0).rotation_unitary,
            atol=1e-12,
        )

    def test_acts_on_two_modes(self):
        assert ModeRotation.beamsplitter(0.3).num_modes == 2

    def test_is_a_bosonic_gate(self):
        gate = ModeRotation.beamsplitter(0.3)
        assert isinstance(gate, BosonicGate)
        assert gate.name == "ModeRotation"


class TestPhaseShift(unittest.TestCase):
    def test_is_diagonal_and_unitary(self):
        phases = [0.1, -0.7, 2.5]
        unitary = ModeRotation.phase_shift(phases).rotation_unitary
        np.testing.assert_allclose(unitary, np.diag(np.exp(1j * np.array(phases))), atol=1e-12)
        np.testing.assert_allclose(unitary @ unitary.conj().T, np.eye(3), atol=1e-12)

    def test_zero_phases_are_the_identity(self):
        np.testing.assert_allclose(
            ModeRotation.phase_shift([0.0, 0.0]).rotation_unitary, np.eye(2), atol=1e-12
        )

    def test_num_modes_is_the_phase_count(self):
        assert ModeRotation.phase_shift([0.1] * 5).num_modes == 5

    def test_composes_additively(self):
        first = ModeRotation.phase_shift([0.3, 0.9]).rotation_unitary
        second = ModeRotation.phase_shift([0.4, -0.2]).rotation_unitary
        combined = ModeRotation.phase_shift([0.7, 0.7]).rotation_unitary
        np.testing.assert_allclose(second @ first, combined, atol=1e-12)


class TestGeneratorRoundTrip(unittest.TestCase):
    def test_generator_is_anti_hermitian(self):
        for gate in (
            ModeRotation.beamsplitter(0.6, 0.4),
            ModeRotation.phase_shift([0.2, -0.5, 1.1]),
            ModeRotation(random_unitary(3)),
        ):
            generator = gate.generator()
            np.testing.assert_allclose(generator.conj().T, -generator, atol=1e-10)

    def test_exponentiating_the_generator_recovers_the_unitary(self):
        for gate in (
            ModeRotation.beamsplitter(0.6, 0.4),
            ModeRotation.phase_shift([0.2, -0.5, 1.1]),
            ModeRotation(random_unitary(4)),
        ):
            np.testing.assert_allclose(
                scipy.linalg.expm(gate.generator()), gate.rotation_unitary, atol=1e-10
            )

    def test_from_generator_round_trip(self):
        unitary = random_unitary(3)
        gate = ModeRotation(unitary)
        rebuilt = ModeRotation.from_generator(gate.generator())
        np.testing.assert_allclose(rebuilt.rotation_unitary, unitary, atol=1e-10)

    def test_from_generator_of_an_anti_hermitian_matrix_is_unitary(self):
        hermitian = RNG.normal(size=(3, 3))
        hermitian = hermitian + hermitian.T
        unitary = ModeRotation.from_generator(1j * hermitian).rotation_unitary
        np.testing.assert_allclose(unitary @ unitary.conj().T, np.eye(3), atol=1e-10)

    def test_from_generator_of_zero_is_the_identity(self):
        np.testing.assert_allclose(
            ModeRotation.from_generator(np.zeros((3, 3))).rotation_unitary,
            np.eye(3),
            atol=1e-12,
        )


class TestToOperator(unittest.TestCase):
    """The bridge from the matrix form to the operator algebra."""

    def test_single_particle_block_is_the_hermitian_generator(self):
        gate = ModeRotation(random_unitary(3))
        operator = gate.to_operator()
        np.testing.assert_allclose(
            single_particle_block(operator, 3), 1j * gate.generator(), atol=1e-10
        )

    def test_operator_is_number_conserving(self):
        """A mode rotation is *passive*: it is generated by ``b^dagger_i b_j`` terms only."""
        operator = ModeRotation(random_unitary(3)).to_operator()
        assert operator.conserves_particle_number()

    def test_operator_is_exactly_hermitian(self):
        """``H = i log(U)`` is Hermitian to the last bit, not merely up to the round-off of ``logm``.

        A ``PauliEvolutionGate`` rejects the imaginary part that any residual would leave behind.
        """
        for num_modes in (2, 3, 4):
            operator = ModeRotation(random_unitary(num_modes)).to_operator()
            assert operator.is_hermitian(0.0)

    def test_the_antihermitian_generator_is_minus_i_times_the_operator(self):
        operator = ModeRotation(random_unitary(3)).to_operator()
        assert not (1j * operator).is_hermitian(1e-9)

    def test_exp_of_the_hermitian_generator_reproduces_the_unitary(self):
        r"""``exp(-i H) = U`` on the single-particle block, with ``H = i log(U)``."""
        for num_modes in (2, 3):
            gate = ModeRotation(random_unitary(num_modes))
            hermitian = single_particle_block(gate.to_operator(), num_modes)
            np.testing.assert_allclose(hermitian, hermitian.conj().T, atol=1e-10)
            np.testing.assert_allclose(
                scipy.linalg.expm(-1j * hermitian), gate.rotation_unitary, atol=1e-10
            )

    def test_beamsplitter_generator(self):
        gate = ModeRotation.beamsplitter(np.pi / 4, 0.0)
        hermitian = single_particle_block(gate.to_operator(), 2)
        np.testing.assert_allclose(
            scipy.linalg.expm(-1j * hermitian), gate.rotation_unitary, atol=1e-10
        )

    def test_tolerance_drops_round_off_terms(self):
        """``logm`` of a block-diagonal unitary leaves noise in nominally-zero entries."""
        gate = ModeRotation.phase_shift([0.3, 0.7, -0.2])
        # a diagonal rotation has a diagonal generator, i.e. exactly `num_modes` terms
        assert len(gate.to_operator()) == 3
        # a wide-open tolerance drops everything
        assert len(gate.to_operator(tol=10.0)) == 0

    def test_identity_rotation_yields_the_zero_operator(self):
        gate = ModeRotation(np.eye(3))
        assert len(gate.to_operator()) == 0
        assert gate.to_operator().equiv(BosonOperator.zero())

    def test_operator_support_is_the_acted_modes(self):
        gate = ModeRotation(random_unitary(3))
        assert gate.to_operator().get_support() == {0, 1, 2}

    def test_max_occupation_of_the_generator_is_one(self):
        """Each generator term raises one mode once, so it demands only ``d >= 2``."""
        assert ModeRotation(random_unitary(3)).to_operator().max_occupation() == 1


class TestArbitraryUnitary(unittest.TestCase):
    def test_num_modes_is_the_matrix_dimension(self):
        for dim in (1, 2, 3, 5):
            assert ModeRotation(random_unitary(dim)).num_modes == dim

    def test_input_is_cast_to_complex(self):
        gate = ModeRotation(np.eye(2, dtype=int))
        assert gate.rotation_unitary.dtype == complex

    def test_accepts_a_nested_list(self):
        gate = ModeRotation([[0.0, 1.0], [1.0, 0.0]])
        np.testing.assert_allclose(gate.rotation_unitary, np.array([[0.0, 1.0], [1.0, 0.0]]))

    def test_unitarity_is_not_verified(self):
        """Documented as the caller's responsibility, so a non-unitary matrix constructs fine."""
        gate = ModeRotation(np.array([[2.0, 0.0], [0.0, 2.0]]))
        assert gate.num_modes == 2
