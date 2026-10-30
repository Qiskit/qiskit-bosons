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

"""B2BSynthesis tests for ``ModeRotation`` gates.

The pass breaks every ``ModeRotation`` of a circuit down in bosonic space, ahead of the layout and
the boson-to-qubit synthesis. It must not depend on the truncation or the encoding, and the method
it is given takes precedence over the one of the gate. The same pipeline must also accept a gate that
is left whole when a plugin for it is registered.
"""

from __future__ import annotations

import numpy as np
import pytest
from qiskit.passmanager import MultiStagePassManager
from qiskit.quantum_info import SparsePauliOp, Statevector, random_unitary
from qiskit.transpiler import PassManager
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes, ModeRotation
from qiskit_bosons.circuit.library.synthesis import BosonicClementsDecomposition
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.transpiler import (
    BosonicCircuitToDAG,
    BosonicPassManager,
    QuantumDAGToCircuit,
)
from qiskit_bosons.transpiler.passes import (
    B2BSynthesis,
    B2QSynthesis,
    MapperFnEvolutionSynthesis,
    OccupationB2QLayout,
    OccupationInitializeModesSynthesis,
)
from qiskit_bosons.transpiler.presets import generate_preset_occupation_pass_manager
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply

from ...mode_rotation_support import BosonicReckDecomposition, DirectModeRotationSynthesis

CLEMENTS = BosonicClementsDecomposition()


def _rotation_circuit(num_modes: int = 3, seed: int = 7, **kwargs) -> BosonicCircuit:
    circuit = BosonicCircuit(num_modes)
    circuit.append(ModeRotation(random_unitary(num_modes, seed=seed).data, **kwargs), circuit.modes)
    return circuit


def _decompose(circuit: BosonicCircuit, synthesis=CLEMENTS, **kwargs) -> BosonicCircuit:
    return BosonicPassManager([B2BSynthesis({"ModeRotation": synthesis}, **kwargs)]).run(circuit)


def test_a_mode_rotation_is_replaced_by_atomic_evolutions():
    decomposed = _decompose(_rotation_circuit())
    assert dict(decomposed.count_ops()) == {"Evolution": 6}
    assert all(
        isinstance(i.operation, Evolution) and i.operation.atomic for i in decomposed._inner.data
    )


def test_the_pass_does_nothing_without_a_method_for_the_gate():
    circuit = _rotation_circuit()
    assert dict(BosonicPassManager([B2BSynthesis()]).run(circuit).count_ops()) == {
        "ModeRotation": 1
    }
    assert dict(
        BosonicPassManager([B2BSynthesis({"Evolution": "LieTrotter"})]).run(circuit).count_ops()
    ) == {"ModeRotation": 1}


def test_the_method_can_be_chosen_by_name():
    decomposed = BosonicPassManager([B2BSynthesis({"ModeRotation": "Clements"})]).run(
        _rotation_circuit()
    )
    assert dict(decomposed.count_ops()) == {"Evolution": 6}


def test_other_gates_are_left_untouched():
    circuit = BosonicCircuit(3)
    circuit.append(InitializeModes([1, 0, 1]), circuit.modes)
    circuit.append(ModeRotation(random_unitary(3, seed=1).data), circuit.modes)
    decomposed = _decompose(circuit)
    assert decomposed.count_ops()["InitializeModes"] == 1
    assert "ModeRotation" not in decomposed.count_ops()


def test_the_pass_is_idempotent():
    once = _decompose(_rotation_circuit())
    assert dict(_decompose(once).count_ops()) == dict(once.count_ops())


def test_a_filter_restricts_the_pass():
    circuit = BosonicCircuit(3)
    circuit.append(ModeRotation(random_unitary(3, seed=1).data), circuit.modes)
    circuit.append(ModeRotation.beamsplitter(0.3), circuit.modes[:2])
    decomposed = _decompose(circuit, filter=lambda node: node.op.num_modes == 2)
    assert decomposed.count_ops()["ModeRotation"] == 1


def test_the_synthesis_of_the_pass_wins_over_the_one_of_the_gate():
    """As for an ``Evolution``: a gate that carries a method is decomposed by the pass."""
    unitary = np.diag(np.exp(1j * np.array([1e-6, 0.5])))
    circuit = BosonicCircuit(2)
    circuit.append(
        ModeRotation(unitary, synthesis=BosonicClementsDecomposition(tol=1e-3)), circuit.modes
    )
    assert _decompose(circuit).count_ops()["Evolution"] == 2
    assert _decompose(circuit, BosonicClementsDecomposition(tol=1e-3)).count_ops()["Evolution"] == 1


def test_the_gates_of_the_circuit_are_not_mutated():
    circuit = _rotation_circuit()
    gate = circuit._inner.data[0].operation
    _decompose(circuit, BosonicClementsDecomposition(tol=1e-3))
    assert gate.synthesis is None


def test_the_same_decomposition_serves_every_encoding():
    """The bosonic circuit does not depend on the encoding; only its lowering does."""
    decomposed = _decompose(_rotation_circuit())
    for encoding in ("unary_projected", "unary_unprojected", "standard_binary", "gray"):
        pass_manager = generate_preset_occupation_pass_manager(LevelSpec(3, encoding))
        pass_manager.qubit = PassManager()
        lowered = pass_manager.run(decomposed)
        assert lowered.count_ops()["PauliEvolution"] == 6


def _pipeline(spec: LevelSpec, *, decomposition=None, direct=False) -> MultiStagePassManager:
    """A hand-built pipeline, with the mode rotation either decomposed or mapped directly."""

    def mapper_fn(operator, num_qubits):
        return occupation_encoding(operator, spec, num_qubits)

    synthesis = B2QSynthesis()
    synthesis.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)
    synthesis.methods["InitializeModes"] = OccupationInitializeModesSynthesis(spec)
    if direct:
        synthesis.methods["ModeRotation"] = DirectModeRotationSynthesis(mapper_fn)
    stages = {}
    if decomposition is not None:
        stages["optimization"] = BosonicPassManager([B2BSynthesis({"ModeRotation": decomposition})])
    return MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        **stages,
        layout=OccupationB2QLayout(spec),
        synthesis=synthesis,
        output=QuantumDAGToCircuit(),
    )


def test_a_rotation_without_a_decomposition_or_a_plugin_names_both_remedies():
    with pytest.raises(TypeError) as excinfo:
        _pipeline(LevelSpec(1, "gray")).run(_rotation_circuit(2))
    assert "B2BSynthesis" in str(excinfo.value)
    assert "register a plugin" in str(excinfo.value)


def _amplitudes(pass_manager, spec: LevelSpec, unitary: np.ndarray) -> float:
    """The worst deviation of the single-boson amplitudes from ``U[target, source]``."""
    num_modes = len(unitary)
    encoding = spec.encoding
    width = spec.resolve(num_modes).widths[0]
    worst = 0.0
    for source in range(num_modes):
        circuit = BosonicCircuit(num_modes)
        circuit.append(InitializeModes([int(m == source) for m in range(num_modes)]), circuit.modes)
        circuit.append(ModeRotation(unitary), circuit.modes)
        lowered = pass_manager.run(circuit)
        state = Statevector.from_int(0, 2**lowered.num_qubits).evolve(lowered).data
        for target in range(num_modes):
            occupation = [int(m == target) for m in range(num_modes)]
            index = sum(
                encoding.level_to_bits(n) << (width * mode) for mode, n in enumerate(occupation)
            )
            worst = max(worst, abs(state[index] - unitary[target, source]))
    return worst


ENCODINGS = ["unary_projected", "unary_unprojected", "standard_binary", "gray"]


@pytest.mark.parametrize("num_modes", [2, 3, 4])
@pytest.mark.parametrize("encoding", ENCODINGS)
@pytest.mark.parametrize(
    "decomposition", [CLEMENTS, BosonicReckDecomposition()], ids=["clements", "reck"]
)
def test_every_decomposition_reproduces_the_rotation(num_modes, encoding, decomposition):
    """The amplitudes ``<target| R |source>`` are ``U[target, source]``, checked analytically.

    The interface is not tied to one decomposition, and neither is the pipeline or the encoding.
    """
    unitary = random_unitary(num_modes, seed=7 + num_modes).data
    spec = LevelSpec(1, encoding)
    assert _amplitudes(_pipeline(spec, decomposition=decomposition), spec, unitary) < 1e-9


@pytest.mark.parametrize("num_modes", [2, 3])
@pytest.mark.parametrize("encoding", ENCODINGS)
def test_a_plugin_for_the_whole_rotation_reproduces_it_without_any_decomposition(
    num_modes, encoding
):
    """The framework lowers a ``ModeRotation`` that is left whole when a plugin for it exists."""
    unitary = random_unitary(num_modes, seed=7 + num_modes).data
    spec = LevelSpec(1, encoding)
    assert _amplitudes(_pipeline(spec, direct=True), spec, unitary) < 1e-9


def test_the_direct_plugin_emits_one_evolution_for_the_whole_rotation():
    spec = LevelSpec(1, "gray")
    lowered = _pipeline(spec, direct=True).run(_rotation_circuit())
    assert dict(lowered.count_ops()) == {"PauliEvolution": 1}


def test_reck_needs_more_depth_than_clements_for_the_same_rotation():
    """It is only a test of the interface, and not a competitor to the library mesh."""
    unitary = random_unitary(5, seed=3).data

    def depth(decomposition) -> int:
        circuit = BosonicCircuit(5)
        circuit.append(ModeRotation(unitary, synthesis=decomposition), circuit.modes)
        return circuit.decompose()._inner.depth()

    assert depth(BosonicReckDecomposition()) > depth(CLEMENTS)


def _final_state(circuit: BosonicCircuit, spec: LevelSpec, pass_manager=None) -> np.ndarray:
    pass_manager = pass_manager or generate_preset_occupation_pass_manager(spec)
    pass_manager.qubit = PassManager()  # keep each factor as one exact PauliEvolutionGate
    lowered = pass_manager.run(circuit)
    return Statevector.from_int(0, 2**lowered.num_qubits).evolve(lowered).data


def _whole_generator_state(unitary: np.ndarray, initial: np.ndarray, spec: LevelSpec) -> np.ndarray:
    """Exponentiates the generator of the whole rotation in one go, truncated once."""
    resolved = spec.resolve(len(unitary))
    observable = occupation_encoding(
        ModeRotation(unitary).to_operator(), spec, resolved.total_width
    )
    matrix = csc_matrix(SparsePauliOp.from_sparse_observable(observable).to_matrix(sparse=True))
    return expm_multiply(-1j * matrix, initial)


@pytest.mark.parametrize("encoding", ["unary_unprojected", "standard_binary", "gray"])
def test_a_cutoff_of_the_boson_number_makes_mesh_and_whole_rotation_agree(encoding):
    """Each factor conserves the boson number, so it cannot take a mode past ``N``."""
    unitary = random_unitary(3, seed=7).data
    occupation = [1, 1, 1]
    spec = LevelSpec(sum(occupation), encoding)

    prepare = BosonicCircuit(3)
    prepare.append(InitializeModes(occupation), prepare.modes)
    initial = _final_state(prepare, spec)

    mesh = BosonicCircuit(3)
    mesh.append(InitializeModes(occupation), mesh.modes)
    mesh.append(ModeRotation(unitary), mesh.modes)

    np.testing.assert_allclose(
        _final_state(mesh, spec), _whole_generator_state(unitary, initial, spec), atol=1e-9
    )


def test_a_cutoff_below_the_boson_number_makes_them_differ():
    """The truncated factors do not multiply to the truncated rotation (see the class docs)."""
    unitary = random_unitary(3, seed=7).data
    occupation = [1, 1, 1]
    spec = LevelSpec(2, "gray")

    prepare = BosonicCircuit(3)
    prepare.append(InitializeModes(occupation), prepare.modes)
    initial = _final_state(prepare, spec)

    mesh = BosonicCircuit(3)
    mesh.append(InitializeModes(occupation), mesh.modes)
    mesh.append(ModeRotation(unitary), mesh.modes)

    difference = np.abs(_final_state(mesh, spec) - _whole_generator_state(unitary, initial, spec))
    assert difference.max() > 0.1


def test_the_direct_plugin_is_the_whole_generator_truncated_once():
    """It agrees with the whole-generator path at *any* cutoff, unlike the mesh."""
    unitary = random_unitary(3, seed=7).data
    occupation = [1, 1, 1]
    spec = LevelSpec(2, "gray")

    prepare = BosonicCircuit(3)
    prepare.append(InitializeModes(occupation), prepare.modes)
    initial = _final_state(prepare, spec)

    direct = BosonicCircuit(3)
    direct.append(InitializeModes(occupation), direct.modes)
    direct.append(ModeRotation(unitary), direct.modes)

    np.testing.assert_allclose(
        _final_state(direct, spec, _pipeline(spec, direct=True)),
        _whole_generator_state(unitary, initial, spec),
        atol=1e-9,
    )
