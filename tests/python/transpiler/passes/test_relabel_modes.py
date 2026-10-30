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

"""Mode relabeling tests.

Unlike the fermionic pass, ``RelabelModes`` only takes an explicit permutation: there is no
excitation-span model to derive one from (an occupation encoding emits no parity strings), so the
pyomo-based and ``find_permutation`` cases have no analogue here.
"""

from __future__ import annotations

import pytest
from qiskit.circuit import QuantumRegister
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.dagcircuit import DAGCircuit
from qiskit.passmanager import MultiStagePassManager
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator
from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
from qiskit_bosons.transpiler.passes import (
    B2QSynthesis,
    MapperFnEvolutionSynthesis,
    OccupationB2QLayout,
    OccupationInitializeModesSynthesis,
    RelabelModes,
)

SPEC = LevelSpec(1, "standard_binary")
"""One qubit per mode: a two-level truncation under the binary encoding."""


def mapper_fn(operator, num_qubits):
    """Maps with the one spec that the layout of these tests also uses."""
    return occupation_encoding(operator, SPEC, num_qubits)


def _hopping(a: int, b: int, weight: float = 1.0) -> dict:
    return {((True, a), (False, b)): weight, ((True, b), (False, a)): weight}


def _qubit_support(gate: PauliEvolutionGate) -> set[int]:
    """The qubits the evolved operator acts on non-trivially."""
    return {qubit for _, indices, _ in gate.operator.to_sparse_list() for qubit in indices}


def _pipeline(relabel: RelabelModes) -> MultiStagePassManager:
    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)
    synth.methods["InitializeModes"] = OccupationInitializeModesSynthesis(SPEC)
    return MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        optimization=relabel,
        layout=OccupationB2QLayout(SPEC),
        synthesis=synth,
        output=QuantumDAGToCircuit(),
    )


def _evolution_support(circ: BosonicCircuit, relabel: RelabelModes) -> set[int]:
    qu_circ = _pipeline(relabel).run(circ)
    (gate,) = [i.operation for i in qu_circ.data if isinstance(i.operation, PauliEvolutionGate)]
    return _qubit_support(gate)


def test_relabel_modes_fixed_permutation():
    """Original mode ``i`` ends up at ``permutation[i]``, which moves the evolved qubits."""
    num_modes = 4
    circ = BosonicCircuit(num_modes)
    circ.append(InitializeModes([1, 0, 1, 0]), circ.modes)
    circ.append(Evolution(num_modes, BosonOperator.from_dict(_hopping(0, 1)), time=1.5), circ.modes)

    # Without a relabeling the hopping acts on qubits 0 and 1.
    assert _evolution_support(circ, RelabelModes([0, 1, 2, 3])) == {0, 1}

    # Both modes are moved: the pass places the register as ``[orig[p] for p in permutation]``.
    relabeled = _evolution_support(circ, RelabelModes([2, 0, 3, 1]))
    assert relabeled != {0, 1}
    assert len(relabeled) == 2


def test_relabel_modes_records_the_permutation():
    permutation = [2, 0, 3, 1]
    circ = BosonicCircuit(4)
    out_dag = RelabelModes(permutation).run(BosonicCircuitToDAG().run(circ))

    assert out_dag.metadata["permutation"] == permutation
    # the metadata is a copy, so mutating it does not reach back into the pass
    out_dag.metadata["permutation"].append(99)
    assert permutation == [2, 0, 3, 1]


def test_relabel_modes_local_indices():
    """An evolution acting on a sub-register keeps acting on exactly its own modes' qubits."""
    num_modes = 4
    circ = BosonicCircuit(2 * num_modes)
    circ.append(InitializeModes([1, 0, 1, 0]), circ.modes[:num_modes])
    circ.append(
        Evolution(num_modes, BosonOperator.from_dict(_hopping(0, 2)), time=1.5),
        circ.modes[:num_modes],
    )

    # Swapping modes 1 and 2 moves the hopping from qubits {0, 2} to qubits {0, 1}.
    permutation = [0, 2, 1, 3, 4, 5, 6, 7]
    assert _evolution_support(circ, RelabelModes([0, 1, 2, 3, 4, 5, 6, 7])) == {0, 2}
    assert _evolution_support(circ, RelabelModes(permutation)) == {0, 1}


def test_relabel_modes_keeps_repeated_mode_terms():
    """A repeated mode (no exclusion principle) survives relabeling as a single, still-local term.

    The term is ``b†_0 b_0 b†_0 b_0 = n_0^2``: with the two-level truncation used here ``b†b†bb``
    would vanish, which would make for an empty image and a vacuous test.
    """
    num_modes = 3
    circ = BosonicCircuit(num_modes)
    operator = BosonOperator.from_dict({((True, 0), (False, 0), (True, 0), (False, 0)): 1.0})
    circ.append(Evolution(num_modes, operator, time=0.5), circ.modes)

    assert _evolution_support(circ, RelabelModes([1, 2, 0])) == {1}


def test_relabel_modes_relabels_barrier():
    """A barrier is carried onto the relabeled register along with the gates.

    The pass rebuilds the register as ``[orig[p] for p in permutation]`` and re-emits each node on
    the same *positions* of that new register, so a barrier on positions 0 and 1 stays on positions
    0 and 1 while the modes sitting underneath them change. The permutation is asymmetric, so a
    barrier that was dropped, or placed on the wrong indices, would not match either assertion.
    """
    num_modes = 4
    circ = BosonicCircuit(num_modes)
    circ.barrier(circ.register[0], circ.register[1])

    dag = BosonicCircuitToDAG().run(circ)
    out_dag = RelabelModes([2, 0, 3, 1]).run(dag)

    nodes = list(out_dag.topological_op_nodes())
    assert [node.op.name for node in nodes] == ["barrier"]

    (relabeled_register,) = out_dag.qregs.values()
    assert sorted(relabeled_register.index(qubit) for qubit in nodes[0].qargs) == [0, 1]
    # those two positions now hold the modes originally at indices 2 and 0
    assert sorted(circ.register.index(qubit) for qubit in nodes[0].qargs) == [0, 2]


def test_run_rejects_multiple_registers():
    """``run`` only supports a single bosonic register."""
    dag = DAGCircuit()
    dag.add_qreg(QuantumRegister(2, "a"))
    dag.add_qreg(QuantumRegister(2, "b"))

    with pytest.raises(NotImplementedError, match="more than a single register"):
        RelabelModes([0, 1, 2, 3]).run(dag)
