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

"""Occupation-encoding preset pipeline tests.

The payoff of keeping the truncation out of the circuit is that one bosonic circuit lowers under any
encoding by changing one argument. These tests exercise that, and close with an end-to-end physics
check: a small Bose-Hubbard evolution transpiled and simulated against ``expm_multiply`` on the
mapped Hamiltonian.
"""

from __future__ import annotations

import itertools
import unittest
from itertools import pairwise

import numpy as np
import pytest
from qiskit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import Operator, SparsePauliOp, Statevector, random_unitary
from qiskit.synthesis import LieTrotter
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes, ModeRotation
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler.presets import generate_preset_occupation_pass_manager
from scipy.linalg import expm
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply

from ..utils import preset_without_qubit_stage

ENCODINGS = ("unary_projected", "standard_binary", "gray")

NUM_SITES = 3
MAX_OCCUPATION = 2  # d = 3 levels per mode
NUM_LEVELS = MAX_OCCUPATION + 1
TOTAL_TIME = 0.4
OCCUPATION = [1, 0, 1]  # one boson on each end of the chain


def bose_hubbard(
    num_sites: int = NUM_SITES, tunneling: float = 1.0, interaction: float = 2.0
) -> BosonOperator:
    """The normal-ordered Bose-Hubbard chain, grouped so that every group is Hermitian.

    Each bond's two hopping terms share a group -- a lone ``b^dagger_i b_j`` is **not** Hermitian and
    would map to complex coefficients that ``PauliEvolutionGate`` rejects -- and all the mutually
    commuting diagonal on-site terms share one further group.
    """
    edges = [(site, site + 1) for site in range(num_sites - 1)]

    terms: dict = {}
    for left, right in edges:
        terms[(cre(left), ann(right))] = -tunneling
        terms[(cre(right), ann(left))] = -tunneling
    for site in range(num_sites):
        terms[(cre(site), cre(site), ann(site), ann(site))] = 0.5 * interaction

    operator = BosonOperator.from_dict(terms)

    groups = []
    for term, _ in operator.iter_terms():
        modes = {mode for _, mode in term}
        if len(modes) == 1:
            groups.append(len(edges))
        else:
            groups.append(edges.index(tuple(sorted(modes))))
    operator.groups = groups

    return operator


def level_codes(num_levels: int, encoding: str) -> list[int]:
    match encoding:
        case "unary_projected" | "unary_unprojected":
            return [1 << level for level in range(num_levels)]
        case "standard_binary":
            return list(range(num_levels))
        case "gray":
            return [level ^ (level >> 1) for level in range(num_levels)]
        case _:  # pragma: no cover
            raise AssertionError(encoding)


def fock_index(occupation, encoding: str) -> int:
    """The qubit basis index of a product Fock state."""
    codes = level_codes(NUM_LEVELS, encoding)
    per_mode = LevelSpec(MAX_OCCUPATION, encoding).resolve(1).total_width
    index = 0
    for mode, level in enumerate(occupation):
        index |= codes[level] << (mode * per_mode)
    return index


def sparse_hamiltonian(operator: BosonOperator, encoding: str) -> csc_matrix:
    """The mapped Hamiltonian as a sparse matrix, for the exact reference evolution."""
    observable = occupation_encoding(operator, LevelSpec(MAX_OCCUPATION, encoding)).simplify()
    return csc_matrix(SparsePauliOp.from_sparse_observable(observable).to_matrix(sparse=True))


def site_densities(state: np.ndarray, encoding: str) -> list[float]:
    """The expected boson number on each site, read off the physical basis states only."""
    densities = []
    for site in range(NUM_SITES):
        total = 0.0
        for levels in itertools.product(range(NUM_LEVELS), repeat=NUM_SITES):
            total += levels[site] * abs(state[fock_index(levels, encoding)]) ** 2
        densities.append(total)
    return densities


class TestPipelineShape(unittest.TestCase):
    def test_lowers_under_every_encoding(self):
        hamiltonian = bose_hubbard()
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(InitializeModes(OCCUPATION), circuit.modes)
        circuit.append(Evolution(NUM_SITES, hamiltonian, time=TOTAL_TIME), circuit.modes)
        decomposed = circuit.decompose()

        widths = {}
        for encoding in ENCODINGS:
            pass_manager = preset_without_qubit_stage(LevelSpec(MAX_OCCUPATION, encoding))
            lowered = pass_manager.run(decomposed)
            widths[encoding] = lowered.num_qubits
            # one PauliEvolutionGate per group, plus the initialization
            assert lowered.count_ops()["PauliEvolution"] == hamiltonian.num_groups()

        assert widths == {
            "unary_projected": NUM_SITES * NUM_LEVELS,
            "standard_binary": NUM_SITES * 2,
            "gray": NUM_SITES * 2,
        }

    def test_qubit_count_is_num_modes_times_per_mode(self):
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(InitializeModes.vacuum(NUM_SITES), circuit.modes)
        for encoding in ENCODINGS:
            for max_occupation in (1, 2, 4):
                lowered = generate_preset_occupation_pass_manager(
                    LevelSpec(max_occupation, encoding)
                ).run(circuit)
                assert lowered.num_qubits == (
                    LevelSpec(max_occupation, encoding).resolve(NUM_SITES).total_width
                )

    def test_compact_encodings_need_fewer_qubits(self):
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(InitializeModes(OCCUPATION), circuit.modes)
        unary = generate_preset_occupation_pass_manager(
            LevelSpec(MAX_OCCUPATION, "unary_projected")
        ).run(circuit)
        gray = generate_preset_occupation_pass_manager(LevelSpec(MAX_OCCUPATION, "gray")).run(
            circuit
        )
        assert gray.num_qubits < unary.num_qubits

    def test_the_same_circuit_serves_several_truncations(self):
        """Nothing about the circuit changes; only the pass manager's arguments do."""
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(InitializeModes(OCCUPATION), circuit.modes)
        circuit.append(Evolution(NUM_SITES, bose_hubbard(), time=TOTAL_TIME), circuit.modes)
        decomposed = circuit.decompose()

        widths = [
            generate_preset_occupation_pass_manager(LevelSpec(cutoff, "gray"))
            .run(decomposed)
            .num_qubits
            for cutoff in (1, 2, 4, 8)
        ]
        assert widths == [3, 6, 9, 12]
        # the source circuit was never rebuilt
        assert dict(decomposed.count_ops()) == {"InitializeModes": 1, "Evolution": 3}

    def test_ungrouped_evolution_lowers_per_term(self):
        """Without groups the decomposition is term-by-term, so each term gets its own gate."""
        operator = BosonOperator.from_dict(
            {(cre(0), ann(0)): 1.0, (cre(1), ann(1)): 1.0, (cre(0), ann(1)): 0.0}
        )
        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, operator, time=0.1), circuit.modes)
        lowered = preset_without_qubit_stage(LevelSpec(1, "gray")).run(circuit.decompose())
        assert lowered.count_ops()["PauliEvolution"] == len(operator)

    def test_identity_mode_rotation_lowers_to_nothing(self):
        """The mesh of an identity rotation has no beamsplitters and no phases to emit."""
        circuit = BosonicCircuit(2)
        circuit.append(ModeRotation(np.eye(2)), circuit.modes)
        lowered = generate_preset_occupation_pass_manager(LevelSpec(MAX_OCCUPATION, "gray")).run(
            circuit
        )
        assert lowered.num_qubits == 4
        assert "PauliEvolution" not in lowered.count_ops()

    def test_mode_rotation_is_a_mesh_of_local_evolutions(self):
        """A generic unitary gives M(M-1)/2 beamsplitters and M phases."""
        unitary = random_unitary(3, seed=7).data
        circuit = BosonicCircuit(3)
        circuit.append(ModeRotation(unitary), circuit.modes)
        lowered = preset_without_qubit_stage(LevelSpec(MAX_OCCUPATION, "gray")).run(circuit)
        assert lowered.count_ops()["PauliEvolution"] == 3 + 3

    def test_mode_rotation_circuit_converges_with_product_formula_repetitions(self):
        """The factors are exact only as exponentials; the *circuit* needs Trotter repetitions.

        Each beamsplitter generator has non-commuting Pauli terms, so a single Lie-Trotter step
        leaves a large error that shrinks as the repetition count grows.
        """
        circuit = BosonicCircuit(3)
        circuit.append(ModeRotation(random_unitary(3, seed=3).data), circuit.modes)
        lowered = preset_without_qubit_stage(LevelSpec(MAX_OCCUPATION, "gray")).run(circuit)

        def ideal_factor(gate) -> np.ndarray:
            observable = gate.operator
            return expm(
                -1j * gate.time * SparsePauliOp.from_sparse_observable(observable).to_matrix()
            )

        ideal = np.eye(2**lowered.num_qubits, dtype=complex)
        for instruction in lowered.data:
            ideal = ideal_factor(instruction.operation) @ ideal

        errors = []
        for reps in (1, 4, 16):
            # `Operator(gate)` of a PauliEvolutionGate is its exact matrix, so the product formula
            # only takes effect once the gates are decomposed
            synthesized = QuantumCircuit(lowered.num_qubits)
            for instruction in lowered.data:
                gate = instruction.operation
                synthesized.append(
                    PauliEvolutionGate(
                        gate.operator, time=gate.time, synthesis=LieTrotter(reps=reps)
                    ),
                    range(lowered.num_qubits),
                )
            matrix = Operator(synthesized.decompose(reps=3)).data
            errors.append(np.abs(matrix - ideal).max())

        assert errors[0] > 0.1
        assert errors[0] > errors[1] > errors[2]
        assert errors[2] < 0.1

    def test_mode_rotation_matches_the_exact_rotation(self):
        """The mesh reproduces ``U`` on the single-particle sector, for every encoding."""
        unitary = random_unitary(3, seed=3).data
        gate = ModeRotation(unitary)
        circuit = BosonicCircuit(3)
        circuit.append(gate, circuit.modes)
        for encoding in ENCODINGS:
            lowered = preset_without_qubit_stage(LevelSpec(MAX_OCCUPATION, encoding)).run(circuit)
            matrix = np.eye(2**lowered.num_qubits, dtype=complex)
            for instruction in lowered.data:
                observable = instruction.operation.operator
                hamiltonian = SparsePauliOp.from_sparse_observable(observable).to_matrix()
                matrix = expm(-1j * instruction.operation.time * hamiltonian) @ matrix
            for source in range(3):
                for target in range(3):
                    occupation_in = [int(m == source) for m in range(3)]
                    occupation_out = [int(m == target) for m in range(3)]
                    amplitude = matrix[
                        fock_index(occupation_out, encoding), fock_index(occupation_in, encoding)
                    ]
                    assert amplitude == pytest.approx(unitary[target, source], abs=1e-9)

    def test_mode_rotation_generator_maps_to_real_coefficients(self):
        """The generator is exactly Hermitian, so its image has no imaginary part to reject."""
        gate = ModeRotation.beamsplitter(np.pi / 4)
        hermitian = occupation_encoding(gate.to_operator(), LevelSpec(MAX_OCCUPATION, "gray"))

        coefficients = np.array([complex(term.coeff) for term in hermitian])
        assert len(coefficients) > 0
        assert np.abs(coefficients.imag).max() < 1e-14

    def test_layout_mismatch_is_reported(self):
        """A spec for three modes handed a two-mode circuit is reported, not silently padded."""
        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, bose_hubbard(2), time=0.1), circuit.modes)
        # a pass manager built for 3 modes handed a 2-mode circuit
        pass_manager = generate_preset_occupation_pass_manager(
            LevelSpec([MAX_OCCUPATION] * 3, "gray")
        )
        with pytest.raises(ValueError, match="covers 3 modes"):
            pass_manager.run(circuit)

    def test_the_qubit_stage_is_always_there_and_runs_before_the_output_stage(self):
        """The qubit stage works on the DAG, so it has to run before ``output`` makes a circuit."""
        pass_manager = generate_preset_occupation_pass_manager(LevelSpec(1, "gray"))
        stages = list(pass_manager._stages)
        assert stages.index("qubit") < stages.index("output")
        assert stages[-1] == "output"

    def test_the_optimization_stage_decomposes_mode_rotations(self):
        """The stage sits between ``input`` and ``layout`` and holds the mode rotation pass."""
        pass_manager = generate_preset_occupation_pass_manager(LevelSpec(1, "gray"))
        stages = list(pass_manager._stages)
        assert stages.index("input") < stages.index("optimization") < stages.index("layout")
        tasks = pass_manager.optimization.to_flow_controller().tasks
        assert [type(task).__name__ for task in tasks] == ["B2BSynthesis"]

    def test_without_arguments_the_qubit_stage_synthesizes_the_evolution(self):
        """The default qubit stage replaces each ``PauliEvolutionGate`` by basis gates."""
        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, bose_hubbard(2), time=0.1), circuit.modes)
        lowered = generate_preset_occupation_pass_manager(LevelSpec(1, "gray")).run(circuit)
        assert "PauliEvolution" not in lowered.count_ops()
        assert lowered.num_qubits == 2

    def test_keyword_arguments_reach_the_qubit_stage(self):
        """Arguments are forwarded to Qiskit's preset pass manager, so they set the basis gates."""
        circuit = BosonicCircuit(2)
        circuit.append(InitializeModes([1, 0]), circuit.modes)
        circuit.append(Evolution(2, bose_hubbard(2), time=0.1), circuit.modes)
        basis = ["rz", "sx", "x", "cx"]
        lowered = generate_preset_occupation_pass_manager(
            LevelSpec(1, "gray"), basis_gates=basis, optimization_level=0
        ).run(circuit)
        assert set(lowered.count_ops()) <= set(basis)

    def test_an_emptied_qubit_stage_keeps_the_evolution_gates(self):
        """Replacing the stage with an empty pass manager leaves the exact ``PauliEvolutionGate``."""
        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, bose_hubbard(2), time=0.1), circuit.modes)
        lowered = preset_without_qubit_stage(LevelSpec(1, "gray")).run(circuit)
        # the circuit was not decomposed in bosonic space, so the groups share one gate
        assert lowered.count_ops()["PauliEvolution"] == 1


class TestGroupingIsHermitian(unittest.TestCase):
    def test_every_group_is_hermitian(self):
        for group in bose_hubbard().split_out_groups():
            assert group.is_hermitian()

    def test_a_lone_hopping_term_is_rejected_by_pauli_evolution(self):
        """Grouping is not cosmetic: an ungrouped bond term is non-Hermitian and fails to lower."""
        lone = BosonOperator.from_dict({(cre(0), ann(1)): -1.0})
        assert not lone.is_hermitian()

        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, lone, time=0.1), circuit.modes)
        with pytest.raises(Exception):  # noqa: B017 - Qiskit's own complex-coefficient rejection
            generate_preset_occupation_pass_manager(LevelSpec(1, "gray")).run(circuit)


class TestBoseHubbardDynamics(unittest.TestCase):
    """End-to-end physics: transpiled Trotter dynamics against an exact sparse reference."""

    ENCODING = "gray"

    def setUp(self):
        self.hamiltonian = bose_hubbard()
        self.initial = self._initial_for(self.ENCODING)
        self.reference = self._reference_for(self.ENCODING)

    def _initial_for(self, encoding: str) -> np.ndarray:
        """The initial Fock state as a qubit state vector under ``encoding``."""
        width = LevelSpec(MAX_OCCUPATION, encoding).resolve(NUM_SITES).total_width
        initial = np.zeros(2**width, dtype=complex)
        initial[fock_index(OCCUPATION, encoding)] = 1.0
        return initial

    def _reference_for(self, encoding: str) -> list[float]:
        """The *exact* truncated dynamics, from ``expm_multiply`` on the mapped Hamiltonian."""
        exact = expm_multiply(
            -1j * TOTAL_TIME * sparse_hamiltonian(bose_hubbard(), encoding),
            self._initial_for(encoding),
        )
        return site_densities(exact, encoding)

    def _trotter_densities(self, steps: int, encoding: str | None = None) -> list[float]:
        """Densities after ``steps`` repetitions of the whole grouped evolution.

        Convergence comes from repeating the entire grouped ``Evolution`` at ``time / steps``, *not*
        from raising a product formula's repetition count -- the latter cannot interleave separate
        gates; see :meth:`test_raising_lie_trotter_reps_does_not_converge`.
        """
        encoding = encoding or self.ENCODING
        circuit = BosonicCircuit(NUM_SITES)
        for _ in range(steps):
            circuit.append(
                Evolution(NUM_SITES, self.hamiltonian, time=TOTAL_TIME / steps), circuit.modes
            )

        pass_manager = generate_preset_occupation_pass_manager(LevelSpec(MAX_OCCUPATION, encoding))
        lowered = pass_manager.run(circuit.decompose())
        # decompose BEFORE evolving: an un-decomposed PauliEvolutionGate carries its *exact*
        # definition, which bypasses the product formula and makes this test pass vacuously
        basis = lowered.decompose(reps=6)

        state = np.asarray(Statevector(self._initial_for(encoding)).evolve(basis).data)
        return site_densities(state, encoding)

    def test_reference_is_symmetric_and_conserves_the_boson_number(self):
        assert sum(self.reference) == pytest.approx(2.0)
        assert self.reference[0] == pytest.approx(self.reference[2])
        # the middle site starts empty and fills up
        assert 0.0 < self.reference[1] < self.reference[0]

    def test_error_decreases_with_more_steps(self):
        errors = []
        for steps in (1, 2, 4, 8):
            densities = self._trotter_densities(steps)
            errors.append(
                max(abs(got - want) for got, want in zip(densities, self.reference, strict=False))
            )

        for coarse, fine in pairwise(errors):
            assert fine < coarse, f"error did not decrease: {errors}"
        # a first-order product formula reaches a small error quickly at this time
        assert errors[-1] < 0.01
        assert errors[0] < 0.05

    def _physical_weight(self, steps: int, encoding: str) -> float:
        """The probability remaining on representable (in-code) basis states."""
        circuit = BosonicCircuit(NUM_SITES)
        for _ in range(steps):
            circuit.append(
                Evolution(NUM_SITES, self.hamiltonian, time=TOTAL_TIME / steps), circuit.modes
            )
        lowered = generate_preset_occupation_pass_manager(LevelSpec(MAX_OCCUPATION, encoding)).run(
            circuit.decompose()
        )
        state = np.asarray(
            Statevector(self._initial_for(encoding)).evolve(lowered.decompose(reps=6)).data
        )
        return sum(
            abs(state[fock_index(levels, encoding)]) ** 2
            for levels in itertools.product(range(NUM_LEVELS), repeat=NUM_SITES)
        )

    def test_the_mapped_hamiltonian_is_block_diagonal(self):
        """The exact evolution can never leave the physical subspace, under any encoding.

        The mapped Hamiltonian has *exactly* zero matrix elements between representable and
        unrepresentable basis states. This is what makes any observed leakage below attributable to
        the product formula rather than to the encoding.
        """
        for encoding in ENCODINGS:
            dense = np.asarray(sparse_hamiltonian(self.hamiltonian, encoding).todense())
            physical = {
                fock_index(levels, encoding)
                for levels in itertools.product(range(NUM_LEVELS), repeat=NUM_SITES)
            }
            unphysical = [index for index in range(dense.shape[0]) if index not in physical]
            if not unphysical:
                continue
            coupling = np.abs(dense[np.ix_(unphysical, sorted(physical))]).max()
            assert coupling == pytest.approx(0.0, abs=1e-12), encoding

    def test_unary_and_gray_circuits_stay_in_the_physical_subspace_exactly(self):
        """For these two encodings the *Trotterized* circuit also preserves the subspace exactly."""
        for encoding in ("unary_projected", "gray"):
            for steps in (1, 4):
                assert self._physical_weight(steps, encoding) == pytest.approx(1.0, abs=1e-12), (
                    encoding,
                    steps,
                )

    def test_standard_binary_leaks_into_its_unused_code_word(self):
        """Standard binary leaks into the unused ``0b11`` code word -- as Trotter error, so it decays.

        At ``d = 3`` a compact encoding spends 2 qubits on 3 levels, leaving one code word unused. The
        *sum* of the mapped Hamiltonian's terms does not couple to it (see
        :meth:`test_the_mapped_hamiltonian_is_block_diagonal`), but the individual Pauli exponentials
        a product formula emits do, so a finite Trotter step parks a little amplitude there.

        The leakage is therefore first-order Trotter error and must shrink with the step count -- it
        falls off roughly as ``1/steps`` here. Worth knowing when reading occupation measurements off
        a standard-binary circuit at a non-power-of-two ``d``: renormalize onto the physical subspace,
        or use enough steps. The Gray code happens not to couple to its unused word at this ``d``.
        """
        leakage = [1.0 - self._physical_weight(steps, "standard_binary") for steps in (1, 2, 4, 8)]

        assert leakage[0] > 1e-3, leakage
        for coarse, fine in pairwise(leakage):
            assert fine < coarse, leakage
        # roughly first order: each doubling of the step count roughly quarters it here
        assert leakage[-1] < leakage[0] / 10, leakage

    def test_boson_number_sector_mixing_is_trotter_error_and_shrinks(self):
        """Trotterization mixes boson-number sectors, but the leak vanishes with more steps.

        The exact evolution conserves the boson number exactly (the Hamiltonian does), so the total
        density is exactly 2. A first-order product formula does not, because the individual groups
        are only number-conserving *together*. The residual is therefore Trotter error, not a
        representation defect -- it must decrease with the step count and it does.
        """
        totals = [sum(self._trotter_densities(steps)) for steps in (1, 4, 16)]
        deviations = [abs(total - 2.0) for total in totals]

        for coarse, fine in pairwise(deviations):
            assert fine < coarse, deviations
        assert deviations[0] < 0.05
        assert deviations[-1] < 1e-3

        # the exact reference, by contrast, conserves it to machine precision
        assert sum(self.reference) == pytest.approx(2.0, abs=1e-10)

    def test_all_encodings_approach_the_same_physics(self):
        """Every encoding converges to the same exact dynamics, differing only in Trotter error.

        The encodings do **not** agree to machine precision at a fixed step count: each produces a
        different Pauli decomposition, hence a different term ordering in the product formula and a
        different Trotter error (~1e-3 apart at 4 steps). What must agree is the limit, so the test
        checks each encoding against the *exact* reference rather than against each other.
        """
        for encoding in ENCODINGS:
            coarse = self._trotter_densities(1, encoding)
            fine = self._trotter_densities(8, encoding)

            reference = self._reference_for(encoding)
            coarse_error = max(abs(a - b) for a, b in zip(coarse, reference, strict=False))
            fine_error = max(abs(a - b) for a, b in zip(fine, reference, strict=False))

            assert fine_error < coarse_error, encoding
            assert fine_error < 0.01, (encoding, fine_error)

    def test_the_exact_reference_is_encoding_independent(self):
        """The *exact* truncated dynamics is the same under all three encodings, to machine precision.

        This is the invariant the previous test's limit rests on: the encoding is a representation
        choice, so once the product-formula error is removed the physics coincides exactly.
        """
        baseline = self._reference_for(self.ENCODING)
        for encoding in ENCODINGS:
            np.testing.assert_allclose(
                self._reference_for(encoding), baseline, atol=1e-10, err_msg=encoding
            )

    def test_raising_lie_trotter_reps_does_not_converge(self):
        """A product formula cannot interleave *separate* gates, so the outer error plateaus.

        This is the documented flip side of grouping being an upstream decision: the partition
        survives into the circuit as one gate per group, and ``LieTrotter(reps=n)`` only subdivides
        the terms *within* one of them. The remaining error between the non-commuting groups is a
        single outer Trotter step and does not vanish as ``reps`` grows -- unlike repeating the whole
        evolution, which :meth:`test_error_decreases_with_more_steps` shows converging.
        """
        from qiskit.passmanager import MultiStagePassManager
        from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
        from qiskit_bosons.transpiler.passes import (
            B2QSynthesis,
            MapperFnEvolutionSynthesis,
            OccupationB2QLayout,
        )

        def mapper_fn(operator, num_qubits):
            return occupation_encoding(
                operator, LevelSpec(MAX_OCCUPATION, self.ENCODING), num_qubits
            )

        def error_at(reps: int) -> float:
            synthesis = B2QSynthesis()
            synthesis.methods["Evolution"] = MapperFnEvolutionSynthesis(
                mapper_fn, LieTrotter(reps=reps)
            )
            pass_manager = MultiStagePassManager(
                input=BosonicCircuitToDAG(),
                layout=OccupationB2QLayout(LevelSpec(MAX_OCCUPATION, self.ENCODING)),
                synthesis=synthesis,
                output=QuantumDAGToCircuit(),
            )
            circuit = BosonicCircuit(NUM_SITES)
            circuit.append(Evolution(NUM_SITES, self.hamiltonian, time=TOTAL_TIME), circuit.modes)
            basis = pass_manager.run(circuit.decompose()).decompose(reps=6)
            state = np.asarray(Statevector(self.initial).evolve(basis).data)
            densities = site_densities(state, self.ENCODING)
            return max(
                abs(got - want) for got, want in zip(densities, self.reference, strict=False)
            )

        plateau = [error_at(reps) for reps in (1, 4, 16)]

        # it improves only marginally and stalls well above the achievable error ...
        assert plateau[-1] > 0.01, plateau
        assert plateau[-1] > plateau[0] / 3, plateau

        # ... whereas repeating the whole grouped evolution at the same total time converges
        stepped = max(
            abs(got - want)
            for got, want in zip(self._trotter_densities(8), self.reference, strict=False)
        )
        assert stepped < plateau[-1] / 2, (stepped, plateau)

    def test_decomposing_before_evolving_is_load_bearing(self):
        """Skipping ``decompose`` makes a Trotter test vacuous: the exact gate definition is used.

        This guards the methodology of every test above. An un-decomposed
        :class:`~qiskit.circuit.library.PauliEvolutionGate` simulates its *exact* ``exp(-itH)``,
        bypassing the product formula entirely. Put the whole Hamiltonian in **one** group and the
        point is unmissable: the un-decomposed circuit reproduces the reference to machine precision
        -- a synthesis test against it would pass no matter how broken the synthesis was -- while the
        decomposed one shows the real first-order Trotter error.
        """
        single_group = bose_hubbard()
        single_group.groups = [0] * len(single_group)

        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, single_group, time=TOTAL_TIME), circuit.modes)
        lowered = preset_without_qubit_stage(LevelSpec(MAX_OCCUPATION, self.ENCODING)).run(
            circuit.decompose()
        )
        assert lowered.count_ops()["PauliEvolution"] == 1

        undecomposed = site_densities(
            np.asarray(Statevector(self.initial).evolve(lowered).data), self.ENCODING
        )
        decomposed = site_densities(
            np.asarray(Statevector(self.initial).evolve(lowered.decompose(reps=6)).data),
            self.ENCODING,
        )

        exact_error = max(abs(a - b) for a, b in zip(undecomposed, self.reference, strict=False))
        trotter_error = max(abs(a - b) for a, b in zip(decomposed, self.reference, strict=False))

        assert exact_error < 1e-12, exact_error
        assert trotter_error > 0.01, trotter_error

    def test_grouped_evolution_is_not_exact_even_undecomposed(self):
        """With several groups, even the un-decomposed circuit carries an outer Trotter error.

        Each group's gate is exact on its own, but the groups do not commute, so the product across
        them is still a single first-order step. This is why convergence needs repeating the whole
        evolution rather than refining anything inside it.
        """
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, self.hamiltonian, time=TOTAL_TIME), circuit.modes)
        lowered = preset_without_qubit_stage(LevelSpec(MAX_OCCUPATION, self.ENCODING)).run(
            circuit.decompose()
        )
        assert lowered.count_ops()["PauliEvolution"] == self.hamiltonian.num_groups()

        densities = site_densities(
            np.asarray(Statevector(self.initial).evolve(lowered).data), self.ENCODING
        )
        error = max(abs(a - b) for a, b in zip(densities, self.reference, strict=False))
        assert error > 1e-3, error
