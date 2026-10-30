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

"""Tests for the MergeEvolutions optimization pass.

Two evolutions under the same operator compose exactly. The pass must merge those and nothing else,
which for bosons includes operators that are only ``equiv``: the truncation does not commute with
normal ordering, so such a merge would change the circuit.
"""

from __future__ import annotations

import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution
from qiskit_bosons.circuit.library.synthesis import (
    BosonicLieTrotter,
    BosonicQDrift,
    BosonicSuzukiTrotter,
)
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicCircuitToDAG
from qiskit_bosons.transpiler.passes import MergeEvolutions
from scipy.linalg import expm


def _number(mode: int, coeff: float = 1.0) -> BosonOperator:
    return BosonOperator.from_dict({(cre(mode), ann(mode)): coeff})


def _hopping() -> BosonOperator:
    return BosonOperator.from_dict({(cre(0), ann(1)): -1.0, (cre(1), ann(0)): -1.0})


def _merge(circuit: BosonicCircuit, **kwargs):
    """Runs the merge pass on ``circuit`` and returns the resulting DAG."""
    return MergeEvolutions(**kwargs).run(BosonicCircuitToDAG().run(circuit))


def _times(dag) -> list[float]:
    return [node.op.params[0] for node in dag.topological_op_nodes()]


# --- positive cases -----------------------------------------------------------------------------


def test_merge_two_evolutions_sums_the_time():
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, _hopping(), time=0.25), circuit.modes)
    circuit.append(Evolution(2, _hopping(), time=0.5), circuit.modes)

    assert _times(_merge(circuit)) == pytest.approx([0.75])


def test_merge_run_of_three():
    circuit = BosonicCircuit(2)
    for time in (0.1, 0.2, 0.3):
        circuit.append(Evolution(2, _hopping(), time=time), circuit.modes)

    assert _times(_merge(circuit)) == pytest.approx([0.6])


def test_merge_keeps_groups_of_the_first_operator():
    operator = BosonOperator.from_terms_with_groups(
        [(((True, 0), (False, 1)), -1.0, 0), (((True, 1), (False, 0)), -1.0, 0)]
    )
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, operator, time=0.1), circuit.modes)
    circuit.append(Evolution(2, operator, time=0.2), circuit.modes)

    (node,) = _merge(circuit).op_nodes()
    assert node.op.operator.num_groups() == 1


def test_the_same_terms_in_another_order_are_the_same_operator():
    forward = BosonOperator.from_dict({(cre(0), ann(1)): -1.0, (cre(1), ann(0)): -1.0})
    backward = BosonOperator.from_dict({(cre(1), ann(0)): -1.0, (cre(0), ann(1)): -1.0})
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, forward, time=0.1), circuit.modes)
    circuit.append(Evolution(2, backward, time=0.2), circuit.modes)

    assert _times(_merge(circuit)) == pytest.approx([0.3])


def test_operators_equal_within_the_tolerance_merge():
    base = _hopping()
    nudged = BosonOperator.from_dict({(cre(0), ann(1)): -1.0 + 1e-14, (cre(1), ann(0)): -1.0})
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, base, time=0.1), circuit.modes)
    circuit.append(Evolution(2, nudged, time=0.2), circuit.modes)

    assert len(_merge(circuit).op_nodes()) == 1
    assert len(_merge(circuit, atol=1e-18).op_nodes()) == 2


def test_gates_on_different_modes_merge_independently():
    circuit = BosonicCircuit(3)
    for _ in range(2):
        circuit.append(Evolution(1, _number(0), time=0.1), [circuit.modes[0]])
        circuit.append(Evolution(1, _number(0), time=0.1), [circuit.modes[2]])

    assert len(_merge(circuit).op_nodes()) == 2


# --- negative cases -----------------------------------------------------------------------------


def test_different_operators_do_not_merge():
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, _hopping(), time=0.1), circuit.modes)
    circuit.append(Evolution(2, _number(0), time=0.1), circuit.modes)

    assert len(_merge(circuit).op_nodes()) == 2


def test_a_different_wire_order_does_not_merge():
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, _hopping(), time=0.1), circuit.modes)
    circuit.append(Evolution(2, _hopping(), time=0.1), circuit.modes[::-1])

    assert len(_merge(circuit).op_nodes()) == 2


def test_a_barrier_interrupts_a_run():
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, _hopping(), time=0.1), circuit.modes)
    circuit.barrier()
    circuit.append(Evolution(2, _hopping(), time=0.1), circuit.modes)

    assert len([n for n in _merge(circuit).op_nodes() if n.name == "Evolution"]) == 2


def test_an_atomic_and_a_decomposable_gate_do_not_merge():
    circuit = BosonicCircuit(2)
    circuit.append(Evolution(2, _hopping(), time=0.1, atomic=True), circuit.modes)
    circuit.append(Evolution(2, _hopping(), time=0.1), circuit.modes)

    assert len(_merge(circuit).op_nodes()) == 2


def test_two_atomic_gates_merge_and_stay_atomic():
    circuit = BosonicCircuit(2)
    for _ in range(2):
        circuit.append(Evolution(2, _hopping(), time=0.1, atomic=True), circuit.modes)

    (node,) = _merge(circuit).op_nodes()
    assert node.op.atomic


def test_different_synthesis_settings_do_not_merge():
    circuit = BosonicCircuit(2)
    circuit.append(
        Evolution(2, _hopping(), time=0.1, synthesis=BosonicLieTrotter(reps=1)), circuit.modes
    )
    circuit.append(
        Evolution(2, _hopping(), time=0.1, synthesis=BosonicLieTrotter(reps=2)), circuit.modes
    )

    assert len(_merge(circuit).op_nodes()) == 2


class _WeightedSynthesis:
    """A synthesis method whose settings include an array, which ``==`` cannot compare."""

    def __init__(self, weights: list[float]) -> None:
        self.weights = np.asarray(weights, dtype=float)

    def synthesize(self, gate):  # pragma: no cover - never called, only compared
        raise NotImplementedError


def test_no_merge_of_different_array_settings():
    """Array-valued settings are compared by value: different arrays must not merge, or raise."""
    circuit = BosonicCircuit(2)
    for time, weights in ((0.1, [1.0, 2.0]), (0.2, [1.0, 3.0])):
        synthesis = _WeightedSynthesis(weights)
        circuit.append(Evolution(2, _hopping(), time=time, synthesis=synthesis), circuit.modes)

    assert _times(_merge(circuit)) == pytest.approx([0.1, 0.2])


def test_merge_of_equal_array_settings():
    """Equal arrays held in distinct objects describe the same method, so the gates merge."""
    circuit = BosonicCircuit(2)
    for time in (0.1, 0.2):
        synthesis = _WeightedSynthesis([1.0, 2.0])
        circuit.append(Evolution(2, _hopping(), time=time, synthesis=synthesis), circuit.modes)

    assert _times(_merge(circuit)) == pytest.approx([0.3])


def test_qdrift_gates_with_different_weights_do_not_merge():
    """The weights of a ``BosonicQDrift`` are an array; its gates used to make the pass raise."""
    circuit = BosonicCircuit(2)
    for time, weights in ((0.1, [1.0, 2.0]), (0.2, [1.0, 3.0])):
        synthesis = BosonicQDrift(3, weights=weights)
        circuit.append(Evolution(2, _hopping(), time=time, synthesis=synthesis), circuit.modes)

    assert _times(_merge(circuit)) == pytest.approx([0.1, 0.2])


def test_the_synthesis_of_the_first_gate_is_kept():
    synthesis = BosonicSuzukiTrotter(order=2, reps=3)
    circuit = BosonicCircuit(2)
    for _ in range(2):
        circuit.append(Evolution(2, _hopping(), time=0.1, synthesis=synthesis), circuit.modes)

    (node,) = _merge(circuit).op_nodes()
    assert node.op.synthesis is not synthesis or node.op.synthesis.reps == 3
    assert node.op.synthesis.reps == 3


# --- what is bosonic about it -------------------------------------------------------------------


def _ordered_pair() -> tuple[BosonOperator, BosonOperator]:
    """``b b^dagger`` and ``b^dagger b + 1``: ``equiv``, but different once truncated."""
    return (
        BosonOperator.from_dict({(ann(0), cre(0)): 1.0}),
        BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (): 1.0}),
    )


def test_operators_that_are_only_equiv_are_not_merged():
    first, second = _ordered_pair()
    assert first.equiv(second, 1e-12)

    circuit = BosonicCircuit(1)
    circuit.append(Evolution(1, first, time=0.4), circuit.modes)
    circuit.append(Evolution(1, second, time=0.7), circuit.modes)

    assert len(_merge(circuit).op_nodes()) == 2


def test_merging_equiv_operators_would_have_changed_the_circuit():
    """This is why the pass does not use ``equiv``: the truncated unitaries really differ."""
    first, second = _ordered_pair()
    spec = LevelSpec(2, "gray").resolve(1)

    def unitary(operator: BosonOperator, time: float) -> np.ndarray:
        observable = occupation_encoding(operator, spec, spec.total_width)
        return expm(-1j * time * SparsePauliOp.from_sparse_observable(observable).to_matrix())

    separate = unitary(second, 0.7) @ unitary(first, 0.4)
    merged = unitary(first, 0.4 + 0.7)
    assert np.abs(separate - merged).max() > 0.1


def test_merging_identical_operators_does_not_change_the_truncated_circuit():
    operator = _hopping()
    spec = LevelSpec(2, "gray").resolve(2)
    observable = occupation_encoding(operator, spec, spec.total_width)
    hamiltonian = SparsePauliOp.from_sparse_observable(observable).to_matrix()

    separate = expm(-1j * 0.7 * hamiltonian) @ expm(-1j * 0.4 * hamiltonian)
    merged = expm(-1j * (0.4 + 0.7) * hamiltonian)
    np.testing.assert_allclose(separate, merged, atol=1e-12)


# --- interplay with the Trotterization ----------------------------------------------------------


def _grouped_hopping() -> BosonOperator:
    return BosonOperator.from_terms_with_groups(
        [
            (((True, 0), (False, 1)), -1.0, 0),
            (((True, 1), (False, 0)), -1.0, 0),
            (((True, 1), (False, 2)), -1.0, 1),
            (((True, 2), (False, 1)), -1.0, 1),
            (((True, 0), (True, 0), (False, 0), (False, 0)), 0.5, 2),
        ]
    )


def _lowered_unitary(circuit: BosonicCircuit, spec: LevelSpec) -> np.ndarray:
    """The unitary of a circuit of ``Evolution`` gates, each mapped and exponentiated exactly."""
    resolved = spec.resolve(len(circuit.register))
    matrix = np.eye(2**resolved.total_width, dtype=complex)
    for instruction in circuit._inner.data:
        gate = instruction.operation
        modes = [circuit.register.index(mode) for mode in instruction.qubits]
        local = gate.operator.relabel_modes(modes)
        observable = occupation_encoding(local, spec, resolved.total_width)
        hamiltonian = SparsePauliOp.from_sparse_observable(observable).to_matrix()
        matrix = expm(-1j * gate.params[0] * hamiltonian) @ matrix
    return matrix


def test_the_pass_merges_across_gate_boundaries_that_the_synthesis_cannot_see():
    """Two separately synthesized gates meet on the same factor at the seam of the circuits.

    The pass also fuses the two ``(1, 2)`` hops around the quartic term on mode 0, because the two
    act on disjoint modes and are therefore neighbors on their own wires. Disjoint operators
    commute, so that is exact as well.
    """
    operator = _grouped_hopping()
    synthesis = BosonicSuzukiTrotter(order=2)
    circuit = BosonicCircuit(3)
    for _ in range(2):
        circuit.append(Evolution(3, operator, time=0.1, synthesis=synthesis), circuit.modes)
    decomposed = circuit.decompose()
    # [0 1 2 1 0] [0 1 2 1 0]
    assert decomposed.count_ops()["Evolution"] == 10

    merged = _merge(decomposed)
    # the seam 0|0 fuses, and so do both `1 2 1` triples across the mode-0 term
    assert len(merged.op_nodes()) == 7


def test_merging_is_exact_on_the_truncated_circuit():
    operator = _grouped_hopping()
    circuit = BosonicCircuit(3)
    for _ in range(3):
        circuit.append(
            Evolution(3, operator, time=0.1, synthesis=BosonicSuzukiTrotter(order=2)), circuit.modes
        )
    decomposed = circuit.decompose()
    merged = BosonicCircuit(3)
    merged._inner = _merged_inner(decomposed)

    spec = LevelSpec(2, "gray")
    assert len(merged._inner.data) < len(decomposed._inner.data)
    np.testing.assert_allclose(
        _lowered_unitary(merged, spec), _lowered_unitary(decomposed, spec), atol=1e-10
    )


def _merged_inner(circuit: BosonicCircuit):
    """The circuit after the pass, converted back through the DAG."""
    from qiskit.converters import dag_to_circuit

    return dag_to_circuit(_merge(circuit))


def test_the_pass_agrees_with_a_synthesis_that_repeats_the_formula():
    """``repeat`` of a step, merged by the pass, is no longer than one synthesis with ``reps``."""
    operator = _grouped_hopping()
    step = BosonicCircuit(3)
    step.append(
        Evolution(3, operator, time=0.1, synthesis=BosonicSuzukiTrotter(order=2)), step.modes
    )
    repeated = step.repeat(3).decompose()
    assert repeated.count_ops()["Evolution"] == 15

    single = BosonicCircuit(3)
    single.append(
        Evolution(3, operator, time=0.3, synthesis=BosonicSuzukiTrotter(order=2, reps=3)),
        single.modes,
    )
    assert single.decompose().count_ops()["Evolution"] == 13
    assert len(_merge(repeated).op_nodes()) <= 13
