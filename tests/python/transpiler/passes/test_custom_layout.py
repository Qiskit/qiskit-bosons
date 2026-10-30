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

"""CustomB2QLayout tests.

The fermionic counterpart pins a Derby-Klassen encoding, which is Majorana-specific. What carries
over is the contract of the pass: it installs the user's layout verbatim, and the synthesis stage
then honors that layout, including a qubit register that is wider than the mode register.
"""

from __future__ import annotations

from qiskit.circuit import QuantumRegister
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.passmanager import MultiStagePassManager
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator
from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
from qiskit_bosons.transpiler.passes import (
    B2QSynthesis,
    CustomB2QLayout,
    MapperFnEvolutionSynthesis,
)

MAX_OCCUPATION = 2
ENCODING = "unary_projected"
SPEC = LevelSpec(MAX_OCCUPATION, ENCODING)


def mapper_fn(operator, num_qubits):
    """Maps the operator with the unary encoding at ``d = 3`` onto ``num_qubits`` qubits."""
    return occupation_encoding(operator, SPEC, num_qubits)


def test_custom_layout_is_installed_verbatim():
    circ = BosonicCircuit(2)
    layout = {circ.register: QuantumRegister(6, "q")}

    pass_ = CustomB2QLayout(layout)
    pass_.property_set = {}
    dag = BosonicCircuitToDAG().run(circ)

    assert pass_.layout is layout
    assert pass_.run(dag) is dag, "this is an analysis pass and must not modify the DAG"
    assert pass_.property_set["b2q_layout"] is layout


def test_custom_layout_drives_the_synthesis_width():
    """The synthesized circuit lives on exactly the qubit register the user supplied."""
    num_modes = 2
    per_mode = SPEC.resolve(1).total_width
    circ = BosonicCircuit(num_modes)
    operator = BosonOperator.from_dict({((True, 0), (False, 1)): 1.0, ((True, 1), (False, 0)): 1.0})
    circ.append(Evolution(num_modes, operator, time=0.5), circ.modes)

    qreg = QuantumRegister(num_modes * per_mode, "custom")
    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)
    pm = MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        layout=CustomB2QLayout({circ.register: qreg}),
        synthesis=synth,
        output=QuantumDAGToCircuit(),
    )

    qu_circ = pm.run(circ)

    assert qu_circ.num_qubits == num_modes * per_mode
    assert [reg.name for reg in qu_circ.qregs] == ["custom"]
    (instruction,) = qu_circ.data
    assert isinstance(instruction.operation, PauliEvolutionGate)
    assert instruction.operation.num_qubits == num_modes * per_mode
