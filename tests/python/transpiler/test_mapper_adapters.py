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

"""Tests for ``simplify``, ``group_wise`` and the term order that the Evolution plugin preserves.

Never convert a mapped operator with ``SparsePauliOp.from_sparse_observable`` for anything but a
tiny one: it expands every projector into two Paulis, so a unary diagonal term of ``w`` qubits
becomes ``2**w`` terms. The helpers below read ``SparseObservable.coeffs`` directly instead.
"""

from __future__ import annotations

import numpy as np
import pytest
from qiskit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.passmanager import MultiStagePassManager
from qiskit.quantum_info import Operator, SparseObservable, SparsePauliOp, Statevector
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
from qiskit_bosons.transpiler.passes import (
    B2QSynthesis,
    MapperFnEvolutionSynthesis,
    OccupationB2QLayout,
    OccupationInitializeModesSynthesis,
    expand_projectors,
    group_wise,
    simplify,
)
from scipy.linalg import expm

from ..trotter_support import (
    MAX_OCCUPATION,
    NUM_SITES,
    bose_hubbard,
    distance,
    evolution_circuit,
    exact_state,
)

TIME = 0.8
NUM_QUBITS = NUM_SITES * (MAX_OCCUPATION + 1)  # unary: one qubit per level


def _map(operator, num_qubits=NUM_QUBITS, *, encoding="unary_projected"):
    return occupation_encoding(operator, LevelSpec(MAX_OCCUPATION, encoding), num_qubits)


def _max_imag(observable: SparseObservable) -> float:
    coeffs = np.asarray(observable.coeffs)
    return float(np.max(np.abs(coeffs.imag))) if coeffs.size else 0.0


def _labels(observable: SparseObservable) -> list[str]:
    return [str(term) for term in observable]


def _pipeline(mapper_fn) -> MultiStagePassManager:
    synthesis = B2QSynthesis()
    synthesis.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)
    synthesis.methods["InitializeModes"] = OccupationInitializeModesSynthesis(
        LevelSpec(MAX_OCCUPATION, "unary_projected")
    )
    return MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        layout=OccupationB2QLayout(LevelSpec(MAX_OCCUPATION, "unary_projected")),
        synthesis=synthesis,
        output=QuantumDAGToCircuit(),
    )


def _gate_operator(circuit) -> SparseObservable:
    (gate,) = [i.operation for i in circuit.data if isinstance(i.operation, PauliEvolutionGate)]
    return gate.operator


class TestTheMapperMergesItsOwnResult:
    """The occupation encoding returns a merged operator, as the fermionic mappers do."""

    @pytest.mark.parametrize("encoding", ["unary_projected", "standard_binary", "gray"])
    def test_nothing_is_left_to_merge(self, encoding):
        image = _map(bose_hubbard(), encoding=encoding)
        assert len(image) == len(image.simplify())

    @pytest.mark.parametrize("encoding", ["unary_projected", "standard_binary", "gray"])
    @pytest.mark.parametrize("scale", [1e-3, 1.0, 1e3, 1e6])
    def test_a_hermitian_operator_maps_to_real_coefficients(self, encoding, scale):
        """The reason the final merge is relative to the largest coefficient.

        The cancellation of a hopping pair leaves rounding noise of about 1e-16 times the size of the
        coefficients (Gray in particular). It is imaginary, and PauliEvolutionGate rejects any
        imaginary part, so an absolute tolerance would fail from some scale on.
        """
        operator = BosonOperator.from_dict(
            {
                (cre(0), ann(1)): -scale,
                (cre(1), ann(0)): -scale,
                (cre(0), cre(0), ann(1), ann(1)): 0.7 * scale,
                (cre(1), cre(1), ann(0), ann(0)): 0.7 * scale,
            }
        )
        image = occupation_encoding(operator, LevelSpec(3, encoding))
        assert _max_imag(image) == 0.0
        PauliEvolutionGate(image, time=0.1)  # must not raise

    def test_a_non_hermitian_term_keeps_its_imaginary_part(self):
        """Only noise is removed. A lone hopping term is genuinely complex, and stays so."""
        lone = BosonOperator.from_dict({(cre(0), ann(1)): -1.0})
        image = occupation_encoding(lone, LevelSpec(2, "unary_projected"))
        assert _max_imag(image) > 0.0
        with pytest.raises(ValueError, match="complex coefficients"):
            PauliEvolutionGate(image, time=0.1)

    def test_tiny_real_coefficients_are_not_lost(self):
        """The relative tolerance must not eat real terms of a small-scale operator."""
        tiny = BosonOperator.from_dict({(cre(0), ann(0)): 1e-12, (cre(1), ann(1)): 1e-12})
        reference = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (cre(1), ann(1)): 1.0})
        assert len(occupation_encoding(tiny, LevelSpec(2, "unary_projected"))) == len(
            occupation_encoding(reference, LevelSpec(2, "unary_projected"))
        )

    def test_the_result_is_reproducible_across_runs(self):
        """The mapping is parallel, but the merged result must not depend on how work was split."""
        operator = bose_hubbard()
        runs = [_labels(_map(operator, encoding="gray")) for _ in range(10)]
        assert all(run == runs[0] for run in runs)


class TestSimplify:
    def test_it_yields_a_canonical_term_order(self):
        operator = bose_hubbard()
        assert _labels(simplify(_map)(operator, NUM_QUBITS)) == _labels(_map(operator).simplify())

    def test_the_value_is_unchanged(self):
        operator = bose_hubbard()
        diff = (simplify(_map)(operator, NUM_QUBITS) - _map(operator)).simplify()
        assert len(diff) == 0

    def test_it_is_idempotent(self):
        once = simplify(_map)(bose_hubbard(), NUM_QUBITS)
        twice = simplify(simplify(_map))(bose_hubbard(), NUM_QUBITS)
        assert _labels(twice) == _labels(once)


class TestGroupWise:
    def test_the_sum_is_the_same_operator(self):
        operator = bose_hubbard()
        diff = (group_wise(_map)(operator, NUM_QUBITS) - _map(operator)).simplify()
        assert len(diff) == 0

    def test_an_ungrouped_operator_is_mapped_in_one_go(self):
        operator = bose_hubbard()
        operator.groups = None
        assert _labels(group_wise(_map)(operator, NUM_QUBITS)) == _labels(_map(operator))

    def test_terms_of_one_group_come_out_adjacent(self):
        """Mapping group by group concatenates the images, so no group's terms are interleaved."""
        operator = bose_hubbard()
        groups = operator.split_out_groups()
        images = [_labels(_map(group)) for group in groups]
        out = _labels(group_wise(_map)(operator, NUM_QUBITS))
        start = 0
        for image in images:
            assert out[start : start + len(image)] == image
            start += len(image)
        assert start == len(out)

    def test_every_group_maps_to_real_coefficients(self):
        for group in bose_hubbard().split_out_groups():
            assert _max_imag(_map(group)) == 0.0
        assert _max_imag(group_wise(_map)(bose_hubbard(), NUM_QUBITS)) == 0.0

    def test_a_non_hermitian_group_leaves_imaginary_coefficients(self):
        """The documented caution: grouping must keep each group Hermitian."""
        operator = bose_hubbard()
        operator.groups = list(range(len(operator)))  # one term per group
        assert _max_imag(group_wise(_map)(operator, NUM_QUBITS)) > 0.0


class TestThePluginPreservesTheMapperOrder:
    def test_the_gate_holds_exactly_what_the_mapper_returned(self):
        """Neither simplified nor re-sorted: the order is part of the mapper's output."""
        operator = bose_hubbard()
        lowered = _pipeline(group_wise(_map)).run(evolution_circuit(operator, TIME))
        assert _labels(_gate_operator(lowered)) == _labels(group_wise(_map)(operator, NUM_QUBITS))

    def test_the_default_mapper_order_is_not_the_canonical_one(self):
        """Guards the test above: if both orders coincided it would prove nothing."""
        grouped = group_wise(_map)(bose_hubbard(), NUM_QUBITS)
        assert _labels(grouped) != _labels(grouped.simplify())

    def test_wrapping_in_simplify_pins_the_canonical_order(self):
        operator = bose_hubbard()
        lowered = _pipeline(simplify(group_wise(_map))).run(evolution_circuit(operator, TIME))
        assert _labels(_gate_operator(lowered)) == _labels(
            group_wise(_map)(operator, NUM_QUBITS).simplify()
        )

    def test_a_mapper_returning_complex_coefficients_is_rejected_by_qiskit(self):
        """The plugin does not clean up after a mapper; PauliEvolutionGate objects instead."""
        operator = bose_hubbard()
        operator.groups = list(range(len(operator)))  # splits every hopping pair
        with pytest.raises(ValueError, match="complex coefficients"):
            _pipeline(group_wise(_map)).run(evolution_circuit(operator, TIME))

    def test_every_ordering_describes_the_same_evolution(self):
        operator = bose_hubbard()
        exact = exact_state(operator, TIME)
        for mapper in (_map, group_wise(_map), simplify(_map)):
            lowered = _pipeline(mapper).run(evolution_circuit(operator, TIME))
            state = Statevector.from_int(0, 2**lowered.num_qubits).evolve(lowered).data
            assert distance(np.asarray(state), exact) < 1e-6


def _xxz_bond_group() -> BosonOperator:
    """One bond of an XXZ chain: hopping, a density-density term and number terms, in one group.

    Hard-core bosons (cutoff one), so a boson is a spin 1/2. The three pieces commute on a bond when
    written as Paulis (XX + YY + ZZ), but not when the number terms are kept apart.
    """
    return BosonOperator.from_terms_with_groups(
        [
            ([cre(0), ann(1)], 1.0, 0),
            ([cre(1), ann(0)], 1.0, 0),
            ([cre(0), ann(0), cre(1), ann(1)], 2.8, 0),
            ([cre(0), ann(0)], -1.4, 0),
            ([cre(1), ann(1)], -1.4, 0),
        ]
    )


def _hard_core_map(operator, num_qubits=2):
    return occupation_encoding(operator, LevelSpec(1, "gray").resolve(2), num_qubits)


class TestExpandProjectors:
    def test_the_gate_of_an_unexpanded_group_is_not_the_exact_evolution(self):
        """The reason the wrapper exists: the cancelling Z pieces are rotated one by one."""
        terms = []
        for left, right in ((0, 1), (2, 3)):
            terms += [
                ([cre(left), ann(right)], 1.0, 0),
                ([cre(right), ann(left)], 1.0, 0),
                ([cre(left), ann(left), cre(right), ann(right)], 2.8, 0),
                ([cre(left), ann(left)], -1.4, 0),
                ([cre(right), ann(right)], -1.4, 0),
            ]
        operator = BosonOperator.from_terms_with_groups(terms)
        spec = LevelSpec(1, "gray").resolve(4)
        raw = occupation_encoding(operator, spec, 4)
        expanded = expand_projectors(lambda op, n: occupation_encoding(op, spec, n))(operator, 4)
        target = Operator(expm(-0.8j * SparsePauliOp.from_sparse_observable(raw).to_matrix()))
        for observable, should_match in ((raw, False), (expanded, True)):
            circuit = QuantumCircuit(4)
            circuit.append(PauliEvolutionGate(observable, time=0.8), range(4))
            assert Operator(circuit.decompose(reps=3)).equiv(target) is should_match

    def test_the_expansion_equals_the_operator(self):
        operator = _xxz_bond_group()
        before = _hard_core_map(operator)
        after = expand_projectors(_hard_core_map)(operator, 2)
        # only Pauli factors are left (`contains_projectors` is newer than the minimum Qiskit)
        projectors = {SparseObservable.BitTerm.ZERO, SparseObservable.BitTerm.ONE}
        assert not projectors & set(after.bit_terms)
        np.testing.assert_allclose(
            SparsePauliOp.from_sparse_observable(after).to_matrix(),
            SparsePauliOp.from_sparse_observable(before).to_matrix(),
            atol=1e-12,
        )

    def test_the_expanded_group_is_the_commuting_pauli_form(self):
        """XX + YY + ZZ per bond: every pair of terms commutes, so one Lie step is exact."""
        image = expand_projectors(_hard_core_map)(_xxz_bond_group(), 2)
        assert sorted(SparsePauliOp.from_sparse_observable(image).paulis.to_labels()) == [
            "II",
            "XX",
            "YY",
            "ZZ",
        ]
        assert _max_imag(image) == 0.0

    def test_each_string_takes_the_place_of_its_first_occurrence(self):
        """Merging must not sort: the order is what a user chose, through the group indices.

        Checked against a plain dictionary merge on random operators with many exact cancellations.
        """
        rng = np.random.default_rng(7)
        for _ in range(300):
            num_qubits = int(rng.integers(1, 5))
            count = int(rng.integers(2, 24))
            labels = ["".join(rng.choice(list("IXYZ"), num_qubits)) for _ in range(count)]
            coeffs = rng.integers(-2, 3, count).astype(float)
            expected = {}
            for label, coeff in zip(labels, coeffs, strict=True):
                expected[label] = expected.get(label, 0.0) + coeff
            expected = {label: coeff for label, coeff in expected.items() if coeff != 0.0}
            if not expected:
                continue

            observable = SparseObservable.from_sparse_pauli_op(SparsePauliOp(labels, coeffs))
            image = expand_projectors(lambda _op, _n, obs=observable: obs)(None, num_qubits)
            result = SparsePauliOp.from_sparse_observable(image)
            assert result.paulis.to_labels() == list(expected)
            np.testing.assert_allclose(result.coeffs.real, list(expected.values()))

    def test_a_fully_cancelling_operator_leaves_a_zero_identity(self):
        observable = SparseObservable.from_sparse_list([("Z", (0,), 1.0), ("Z", (0,), -1.0)], 2)
        image = expand_projectors(lambda _op, _n: observable)(None, 2)
        assert np.allclose(np.asarray(image.coeffs), 0.0)

    def test_group_wise_inside_keeps_a_shared_string_in_each_group(self):
        """A Pauli string that two groups share stays in both; outside, it moves to the first."""
        operator = BosonOperator.from_terms_with_groups(
            [
                ([cre(0), ann(1)], 1.0, 0),
                ([cre(1), ann(0)], 1.0, 0),
                ([cre(0), ann(0)], 0.5, 0),
                ([cre(0), ann(0)], 0.25, 1),
                ([cre(1), ann(1)], 0.25, 1),
            ]
        )
        inside = group_wise(expand_projectors(_hard_core_map))(operator, 2)
        outside = expand_projectors(group_wise(_hard_core_map))(operator, 2)

        def labels(observable):
            return SparsePauliOp.from_sparse_observable(observable).paulis.to_labels()

        assert labels(inside).count("IZ") == 2  # Z on mode 0, once per group
        assert labels(outside).count("IZ") == 1  # merged across the group boundary
        np.testing.assert_allclose(  # the same operator either way
            SparsePauliOp.from_sparse_observable(inside).to_matrix(),
            SparsePauliOp.from_sparse_observable(outside).to_matrix(),
            atol=1e-12,
        )

    def test_a_non_hermitian_group_is_still_rejected_by_qiskit(self):
        lone = BosonOperator.from_terms_with_groups([([cre(0), ann(1)], 1.0, 0)])
        image = expand_projectors(_hard_core_map)(lone, 2)
        with pytest.raises(ValueError, match="complex coefficients"):
            PauliEvolutionGate(image, time=0.1)

    def test_an_operator_without_terms_passes_through(self):
        empty = BosonOperator.from_terms_with_groups([])
        assert np.allclose(np.asarray(expand_projectors(_hard_core_map)(empty, 2).coeffs), 0.0)


def test_helpers_are_importable_from_the_package_root():
    from qiskit_bosons.transpiler import passes

    assert passes.simplify is simplify
    assert passes.group_wise is group_wise
    assert passes.expand_projectors is expand_projectors
