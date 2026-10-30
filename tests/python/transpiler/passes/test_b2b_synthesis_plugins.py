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

"""Tests for the plugin framework of the boson-to-boson synthesis pass."""

from __future__ import annotations

import numpy as np
import pytest
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes, ModeRotation
from qiskit_bosons.circuit.library.synthesis import (
    B2BSynthesisPluginManager,
    BosonicClementsDecomposition,
    BosonicLieTrotter,
    BosonicQDrift,
    BosonicSuzukiTrotter,
)
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicCircuitToDAG, BosonicPassManager
from qiskit_bosons.transpiler.passes import B2BSynthesis

NUM_MODES = 3


def _operator() -> BosonOperator:
    operator = BosonOperator.from_dict(
        {
            (cre(0), ann(1)): -1.0,
            (cre(1), ann(0)): -1.0,
            (cre(1), ann(2)): -0.8,
            (cre(2), ann(1)): -0.8,
        }
    )
    operator.groups = [0, 0, 1, 1]
    return operator


def _circuit(*extra_gates) -> BosonicCircuit:
    circuit = BosonicCircuit(NUM_MODES)
    for gate in extra_gates:
        circuit.append(gate, circuit.modes)
    circuit.append(Evolution(NUM_MODES, _operator(), time=0.5), circuit.modes)
    return circuit


class _Marker:
    """A plugin that is registered nowhere, and counts the gates it is given."""

    def __init__(self):
        self.seen = []

    def synthesize(self, gate):
        self.seen.append(gate)
        return BosonicCircuit(gate.num_modes)


def test_the_builtin_plugins_are_installed():
    manager = B2BSynthesisPluginManager()

    assert {"LieTrotter", "SuzukiTrotter", "QDrift"} <= set(manager.method_names("Evolution"))
    assert manager.method("Evolution", "QDrift") is BosonicQDrift
    assert manager.method("Evolution", "LieTrotter") is BosonicLieTrotter
    assert manager.method("Evolution", "SuzukiTrotter") is BosonicSuzukiTrotter
    assert manager.method("ModeRotation", "Clements") is BosonicClementsDecomposition
    assert {"Evolution", "ModeRotation"} <= set(manager.op_names())


def test_an_unknown_gate_has_no_methods():
    assert B2BSynthesisPluginManager().method_names("NotAGate") == []


def test_config_by_name():
    synth = B2BSynthesis({"Evolution": "LieTrotter"})

    assert isinstance(synth.methods["Evolution"], BosonicLieTrotter)


def test_config_by_name_and_positional_arguments():
    synth = B2BSynthesis({"Evolution": ("LieTrotter", (3,))})

    assert synth.methods["Evolution"].reps == 3


def test_config_by_name_arguments_and_keyword_arguments():
    synth = B2BSynthesis({"Evolution": ("SuzukiTrotter", (), {"order": 4, "reps": 2})})

    method = synth.methods["Evolution"]
    assert isinstance(method, BosonicSuzukiTrotter)
    assert (method.order, method.reps) == (4, 2)


def test_config_by_instance_needs_no_entry_point():
    marker = _Marker()

    synth = B2BSynthesis({"Evolution": marker})

    assert synth.methods["Evolution"] is marker


def test_methods_can_be_populated_by_hand():
    synth = B2BSynthesis()
    assert synth.methods == {}

    synth.methods["Evolution"] = BosonicLieTrotter()

    out = BosonicPassManager(synth).run(_circuit())
    assert out.count_ops()["Evolution"] > 1


def test_an_empty_pass_changes_nothing():
    out = BosonicPassManager(B2BSynthesis()).run(_circuit())

    assert dict(out.count_ops()) == {"Evolution": 1}


@pytest.mark.parametrize(
    ("config", "message"),
    [
        ({"Evolution": "NoSuchMethod"}, "No synthesis method 'NoSuchMethod'"),
        ({"Evolution": ("NoSuchMethod", ())}, "No synthesis method 'NoSuchMethod'"),
        ({"NotAGate": "Anything"}, "No boson-to-boson synthesis is installed for 'NotAGate'"),
    ],
)
def test_an_unknown_gate_or_method_is_rejected_when_the_pass_is_built(config, message):
    with pytest.raises(KeyError, match=message):
        B2BSynthesis(config)


def test_a_gate_without_a_synthesis_cannot_be_configured():
    """``InitializeModes`` breaks down into nothing, so naming a method for it is a mistake."""
    with pytest.raises(KeyError, match="No boson-to-boson synthesis is installed"):
        B2BSynthesis({"InitializeModes": "Anything"})


@pytest.mark.parametrize("value", [42, None, ("LieTrotter",), ["LieTrotter"], object()])
def test_a_value_of_no_known_form_is_rejected(value):
    with pytest.raises(TypeError, match="must be a method name"):
        B2BSynthesis({"Evolution": value})  # type: ignore[dict-item]


def test_a_gate_without_a_synthesis_passes_through_whatever_else_is_configured():
    marker = _Marker()
    circuit = _circuit(InitializeModes([1, 0, 1]))

    out = BosonicPassManager(B2BSynthesis({"Evolution": marker})).run(circuit)

    assert [type(g).__name__ for g in marker.seen] == ["Evolution"]
    assert dict(out.count_ops()) == {"InitializeModes": 1}


def test_only_the_configured_gate_types_are_touched():
    marker = _Marker()

    BosonicPassManager(B2BSynthesis({"InitializeModes": marker})).run(_circuit())

    assert marker.seen == []


def test_the_plugin_is_matched_by_gate_type_name():
    first, second = _Marker(), _Marker()
    synth = B2BSynthesis({"Evolution": first})

    synth.methods["Evolution"] = second
    BosonicPassManager(synth).run(_circuit())

    assert first.seen == []
    assert len(second.seen) == 1


def test_one_pass_serves_several_gate_types():
    unitary = np.array([[0, 1], [1, 0]], dtype=complex)
    circuit = BosonicCircuit(2)
    circuit.append(ModeRotation(unitary), circuit.modes)
    circuit.append(
        Evolution(2, BosonOperator.from_dict({(cre(0), ann(0)): 1.0}), time=0.3), circuit.modes
    )

    out = BosonicPassManager(
        B2BSynthesis({"ModeRotation": "Clements", "Evolution": "LieTrotter"})
    ).run(circuit)

    assert "ModeRotation" not in out.count_ops()
    assert all(
        inst.operation.atomic for inst in out._inner.data if isinstance(inst.operation, Evolution)
    )


def test_a_registered_name_and_an_instance_behave_the_same():
    by_name = BosonicPassManager(B2BSynthesis({"Evolution": ("SuzukiTrotter", (), {"order": 2})}))
    by_instance = BosonicPassManager(B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=2)}))

    a, b = by_name.run(_circuit()), by_instance.run(_circuit())

    assert dict(a.count_ops()) == dict(b.count_ops())


def test_the_pass_runs_on_a_dag_directly():
    dag = BosonicCircuitToDAG().run(_circuit())

    out = B2BSynthesis({"Evolution": BosonicLieTrotter()}).run(dag)

    assert len(out.op_nodes()) == 2


class TestGateSynthesis:
    """The ``synthesis`` and ``atomic`` attributes now belong to every bosonic gate."""

    def test_a_gate_without_a_synthesis_has_no_definition(self):
        gate = InitializeModes([1, 0])
        assert gate.synthesis is None
        assert gate.atomic is False
        circuit = BosonicCircuit(2)
        circuit.append(gate, circuit.modes)
        assert dict(circuit.decompose().count_ops()) == {"InitializeModes": 1}

    def test_with_synthesis_returns_a_copy_and_leaves_the_gate_alone(self):
        gate = Evolution(NUM_MODES, _operator(), time=0.5)
        original = gate.synthesis
        method = BosonicSuzukiTrotter(order=2)

        tagged = gate.with_synthesis(method)

        assert tagged is not gate
        assert tagged.synthesis is method
        assert gate.synthesis is original
        assert isinstance(tagged, Evolution)
        assert tagged.operator is gate.operator

    def test_with_synthesis_discards_the_cached_definition(self):
        gate = Evolution(NUM_MODES, _operator(), time=0.5)
        first = gate.definition
        assert first is not None

        tagged = gate.with_synthesis(BosonicSuzukiTrotter(order=2))

        assert len(tagged.definition.data) > len(first.data)

    def test_an_atomic_gate_has_no_definition_whatever_its_synthesis(self):
        gate = Evolution(NUM_MODES, _operator(), time=0.5, atomic=True)
        assert gate.definition is None

    def test_the_mode_rotation_default_is_still_no_synthesis(self):
        assert ModeRotation(np.eye(2)).synthesis is None

    def test_a_mode_rotation_tagged_with_a_method_decomposes(self):
        unitary = np.array([[0, 1], [1, 0]], dtype=complex)
        tagged = ModeRotation(unitary).with_synthesis(BosonicClementsDecomposition())
        circuit = BosonicCircuit(2)
        circuit.append(tagged, circuit.modes)
        assert "ModeRotation" not in circuit.decompose().count_ops()

    def test_apply_false_tags_a_mode_rotation(self):
        unitary = np.array([[0, 1], [1, 0]], dtype=complex)
        circuit = BosonicCircuit(2)
        circuit.append(ModeRotation(unitary), circuit.modes)

        out = BosonicPassManager(B2BSynthesis({"ModeRotation": "Clements"}, apply=False)).run(
            circuit
        )

        assert dict(out.count_ops()) == {"ModeRotation": 1}
        gate = out._inner.data[0].operation
        assert isinstance(gate.synthesis, BosonicClementsDecomposition)
        assert "ModeRotation" not in out.decompose().count_ops()
