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
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.passmanager import MultiStagePassManager
from qiskit.quantum_info import SparseObservable, Statevector
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
from qiskit_bosons.transpiler.passes import (
    B2QSynthesis,
    MapperFnEvolutionSynthesis,
    OccupationB2QLayout,
    OccupationInitializeModesSynthesis,
    group_wise,
    simplify,
)

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


def test_helpers_are_importable_from_the_package_root():
    from qiskit_bosons.transpiler import passes

    assert passes.simplify is simplify
    assert passes.group_wise is group_wise
