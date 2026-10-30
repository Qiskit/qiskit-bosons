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

"""Evolution synthesis tests.

The ``simplify``/``group_wise`` adapters and the order-preservation of the plugin are covered in
``tests/python/transpiler/test_mapper_adapters.py``. This module covers the plugin's basic
synthesis path and its ``product_formula`` forwarding.
"""

from __future__ import annotations

from qiskit.circuit.library import PauliEvolutionGate
from qiskit.passmanager import MultiStagePassManager
from qiskit.synthesis import LieTrotter, SuzukiTrotter
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator
from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
from qiskit_bosons.transpiler.passes import (
    B2QSynthesis,
    MapperFnEvolutionSynthesis,
    OccupationB2QLayout,
    group_wise,
)

NUM_MODES = 4


SPEC = LevelSpec(1, "standard_binary")
"""One qubit per mode: a two-level truncation under the binary encoding."""


def mapper_fn(operator, num_qubits):
    """Maps with the one spec that the layout of these tests also uses."""
    return occupation_encoding(operator, SPEC, num_qubits)


def _single_evolution_circuit(time: float = 1.5) -> BosonicCircuit:
    hamil = BosonOperator.from_dict(
        {
            ((True, 0), (False, 2)): 2.0,
            ((True, 2), (False, 0)): 2.0,
            ((True, 1), (False, 3)): -2.0,
            ((True, 3), (False, 1)): -2.0,
        }
    )
    circ = BosonicCircuit(NUM_MODES)
    circ.append(Evolution(NUM_MODES, hamil, time=time), circ.modes)
    return circ


def _run(synth: B2QSynthesis, circ: BosonicCircuit | None = None):
    pm = MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        layout=OccupationB2QLayout(SPEC),
        synthesis=synth,
        output=QuantumDAGToCircuit(),
    )
    return pm.run(circ if circ is not None else _single_evolution_circuit())


def _run_with_product_formula(product_formula):
    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn, product_formula)
    return _run(synth)


def test_evolution_gate_synthesis():
    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)

    qu_circ = _run(synth)

    gates = qu_circ.data
    assert len(gates) == 1
    assert isinstance(gates[0].operation, PauliEvolutionGate)
    # the evolution time is carried over verbatim
    assert gates[0].operation.time == 1.5


def test_custom_qubit_ordering():
    """A mapper may relabel modes first, which moves the evolved qubits accordingly."""
    ordering = [0, 2, 1, 3]

    def reordering_mapper(op, num_qubits):
        return mapper_fn(op.relabel_modes(ordering), num_qubits)

    def support(circ):
        (instr,) = circ.data
        return {q for _, idx, _ in instr.operation.operator.to_sparse_list() for q in idx}

    plain = MapperFnEvolutionSynthesis(mapper_fn)
    reordered = MapperFnEvolutionSynthesis(reordering_mapper)

    # a single hop between modes 0 and 1, so the relabeling is visible in the evolved qubits
    hop = BosonOperator.from_dict({((True, 0), (False, 1)): 1.0, ((True, 1), (False, 0)): 1.0})
    circ = BosonicCircuit(NUM_MODES)
    circ.append(Evolution(NUM_MODES, hop, time=1.5), circ.modes)

    synth = B2QSynthesis()
    synth.methods["Evolution"] = plain
    assert support(_run(synth, circ)) == {0, 1}

    # mode 1 is relabeled to 2, so the same hop now acts on qubits 0 and 2
    synth.methods["Evolution"] = reordered
    qu_circ = _run(synth, circ)
    assert isinstance(qu_circ.data[0].operation, PauliEvolutionGate)
    assert support(qu_circ) == {0, 2}


def test_default_product_formula_is_none():
    """Omitting ``product_formula`` leaves the emitted gate on its default synthesis."""
    plugin = MapperFnEvolutionSynthesis(mapper_fn)
    assert plugin.product_formula is None

    (instruction,) = _run_with_product_formula(None).data
    assert isinstance(instruction.operation, PauliEvolutionGate)
    # a PauliEvolutionGate with no explicit synthesis falls back to a single-rep LieTrotter
    assert isinstance(instruction.operation.synthesis, LieTrotter)
    assert instruction.operation.synthesis.reps == 1


def test_product_formula_is_forwarded_to_pauli_evolution_gate():
    """An explicit product formula is attached to the emitted PauliEvolutionGate verbatim."""
    product_formula = SuzukiTrotter(order=2, reps=3)
    (instruction,) = _run_with_product_formula(product_formula).data
    assert isinstance(instruction.operation, PauliEvolutionGate)
    assert instruction.operation.synthesis is product_formula


def test_higher_order_product_formula_deepens_synthesis():
    """A higher-order / repeated product formula yields a deeper decomposition."""
    default_circ = _run_with_product_formula(None).decompose(reps=2)
    suzuki_circ = _run_with_product_formula(SuzukiTrotter(order=2, reps=3)).decompose(reps=2)

    assert suzuki_circ.size() > default_circ.size()


def test_product_formula_via_config_kwargs():
    """The product formula flows through the ``B2QSynthesisConfig`` keyword-argument form."""
    product_formula = SuzukiTrotter(order=2, reps=2)
    config = {"Evolution": ("MapperFn", (mapper_fn,), {"product_formula": product_formula})}
    synth = B2QSynthesis(config)
    assert synth.methods["Evolution"].product_formula is product_formula

    (instruction,) = _run(synth).data
    assert instruction.operation.synthesis is product_formula


def test_adapters_compose_with_a_product_formula():
    """The adapters are plain mapper functions, so they coexist with ``product_formula``."""
    product_formula = SuzukiTrotter(order=2, reps=2)
    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(group_wise(mapper_fn), product_formula)

    (instruction,) = _run(synth).data
    assert instruction.operation.synthesis is product_formula
