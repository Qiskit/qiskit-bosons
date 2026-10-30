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

"""Tests for the order in which passes visit the nodes of a DAG.

A pass that substitutes a node by a sub-circuit leaves the order in which the DAG *stores* its nodes
out of sync with the dependencies between them. A pass that rebuilds its output node by node has to
iterate in circuit order, and not in storage order, or the factors of successive substitutions end
up interleaved: a well-formed circuit, but for a different unitary, because the gates of two
substitutions share qubits and do not commute.
"""

from __future__ import annotations

import numpy as np
from qiskit.passmanager import MultiStagePassManager
from qiskit.quantum_info import Operator, SparsePauliOp
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution
from qiskit_bosons.circuit.library.synthesis import BosonicLieTrotter, BosonicQDrift
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicCircuitToDAG, BosonicPassManager, QuantumDAGToCircuit
from qiskit_bosons.transpiler.passes import (
    B2BSynthesis,
    B2QSynthesis,
    MapperFnEvolutionSynthesis,
    OccupationB2QLayout,
    RelabelModes,
)
from scipy.linalg import expm

NUM_MODES = 3
SPEC = LevelSpec(1, "gray")
"""One qubit per mode keeps the matrices small."""


def _hop_groups(*hops: tuple[int, int, float]) -> BosonOperator:
    """A Hermitian operator holding one adjoint-paired group per hop between two modes."""
    terms = {}
    for i, j, coeff in hops:
        terms[(cre(i), ann(j))] = coeff
        terms[(cre(j), ann(i))] = coeff
    operator = BosonOperator.from_dict(terms)
    # the operator does not preserve the order in which its terms were added, so the group of a term
    # (the pair of modes it hops between) is derived from the terms as stored
    pairs = [frozenset(mode for _, mode in term) for term, _ in operator.iter_terms()]
    ids = {pair: group for group, pair in enumerate(sorted(set(pairs), key=sorted))}
    operator.groups = [ids[pair] for pair in pairs]
    return operator


def _evolve(operator: BosonOperator, time: float) -> np.ndarray:
    observable = occupation_encoding(operator, SPEC, NUM_MODES)
    return expm(-1j * time * SparsePauliOp.from_sparse_observable(observable).to_matrix())


def _lie_trotter(operator: BosonOperator, time: float) -> np.ndarray:
    """The first-order product over the groups of ``operator``, in group order."""
    unitary = np.eye(2**NUM_MODES, dtype=complex)
    for group in sorted(set(operator.groups)):
        unitary = _evolve(operator.split_out_groups([group])[0], time) @ unitary
    return unitary


def _overlapping_evolutions() -> tuple[BosonicCircuit, list[tuple[BosonOperator, float]]]:
    """Three successive evolutions on overlapping modes, which therefore do not commute."""
    gates = [
        (_hop_groups((0, 1, 0.7), (1, 2, 0.3)), 0.6),
        (_hop_groups((1, 2, 0.4), (0, 1, 0.2)), 0.5),
        (_hop_groups((0, 1, 0.9), (1, 2, 0.5)), 0.8),
    ]
    circuit = BosonicCircuit(NUM_MODES)
    for operator, time in gates:
        circuit.append(Evolution(NUM_MODES, operator, time=time), circuit.modes)
    return circuit, gates


def _signature(nodes) -> list[tuple[str, tuple[int, ...], float]]:
    """Name, modes and evolution time of each node, which identify a gate of the examples above."""
    return [
        (
            node.name,
            tuple(mode._index for mode in node.qargs),
            round(float(node.op.params[0]), 9) if node.op.params else 0.0,
        )
        for node in nodes
    ]


def _substituted_dag():
    """A DAG after a pass replaced each of its nodes by a sub-circuit."""
    circuit, _ = _overlapping_evolutions()
    dag = BosonicCircuitToDAG().run(circuit)
    B2BSynthesis({"Evolution": BosonicLieTrotter()}).run(dag)
    return dag


def test_the_fixture_has_storage_out_of_sync_with_the_circuit():
    """If this fails, the tests below no longer exercise the bug they are named after."""
    dag = _substituted_dag()
    assert _signature(dag.op_nodes()) != _signature(dag.topological_op_nodes())


def test_b2b_synthesis_visits_the_gates_in_circuit_order():
    """The pass substitutes nodes itself, so it must not iterate in storage order.

    A first pass splits each gate into plain (not atomic) factors, which leaves the storage of the
    DAG out of sync with its wires. A second pass, with a plugin that records the order in which it
    is handed the gates, then has to see them in circuit order.
    """

    class Split:
        def synthesize(self, gate):
            circuit = BosonicCircuit(gate.num_modes)
            for factor in gate.operator.split_out_groups():
                circuit.append(
                    Evolution(gate.num_modes, factor, time=gate.params[0]), circuit.modes
                )
            return circuit

    class Recorder:
        def __init__(self):
            self.seen = []

        def synthesize(self, gate):
            self.seen.append(round(float(gate.params[0]), 9))
            return BosonicCircuit(gate.num_modes)

    circuit, _ = _overlapping_evolutions()
    dag = BosonicCircuitToDAG().run(circuit)
    B2BSynthesis({"Evolution": Split()}).run(dag)
    assert _signature(dag.op_nodes()) != _signature(dag.topological_op_nodes())

    recorder = Recorder()
    B2BSynthesis({"Evolution": recorder}).run(dag)

    assert recorder.seen == [0.6, 0.6, 0.5, 0.5, 0.8, 0.8]


def test_synthesis_lowers_successive_substitutions_in_circuit_order():
    """The factors of one substituted gate must not be interleaved with those of the next."""
    circuit, gates = _overlapping_evolutions()

    synthesis = B2QSynthesis()
    synthesis.methods["Evolution"] = MapperFnEvolutionSynthesis(
        lambda operator, num_qubits: occupation_encoding(operator, SPEC, num_qubits)
    )
    pass_manager = MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        optimization=BosonicPassManager(B2BSynthesis({"Evolution": BosonicLieTrotter()})),
        layout=OccupationB2QLayout(SPEC),
        synthesis=synthesis,
        output=QuantumDAGToCircuit(),
    )

    lowered = Operator(pass_manager.run(circuit)).data

    expected = np.eye(2**NUM_MODES, dtype=complex)
    for operator, time in gates:
        expected = _lie_trotter(operator, time) @ expected
    np.testing.assert_allclose(lowered, expected, atol=1e-10)


def test_relabel_modes_keeps_the_circuit_order_of_a_substituted_dag():
    """A permutation that changes nothing must leave the circuit exactly as it was."""
    dag = _substituted_dag()
    before = _signature(dag.topological_op_nodes())

    relabeled = RelabelModes(list(range(NUM_MODES))).run(dag)

    assert _signature(relabeled.topological_op_nodes()) == before


def test_qdrift_keeps_the_circuit_order_of_a_substituted_dag():
    """With one group per gate every sample is deterministic, so only the order can differ.

    The sampled gates sit on the whole register, so a gate is identified by its evolution time
    instead: it is the time of the input gate times a constant, so the order of the times is the
    order of the circuit.
    """
    gates = [
        (_hop_groups((0, 1, 1.0)), 0.1),
        (_hop_groups((1, 2, 1.0)), 0.2),
        (_hop_groups((0, 1, 1.0)), 0.3),
    ]
    circuit = BosonicCircuit(NUM_MODES)
    for operator, time in gates:
        circuit.append(Evolution(NUM_MODES, operator, time=time), circuit.modes)
    dag = BosonicCircuitToDAG().run(circuit)
    B2BSynthesis({"Evolution": BosonicLieTrotter()}).run(dag)
    assert _signature(dag.op_nodes()) != _signature(dag.topological_op_nodes())

    sampled = B2BSynthesis({"Evolution": BosonicQDrift(1, rng=0)}).run(dag)

    times = [float(node.op.params[0]) for node in sampled.topological_op_nodes()]
    assert times == sorted(times)
    assert len(set(times)) == len(times)
