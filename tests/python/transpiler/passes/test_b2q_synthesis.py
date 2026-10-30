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

"""Boson-to-qubit synthesis tests.

Measurement is not part of the bosonic circuit yet, so the measurement-related cases of the
fermionic counterpart (clbit copying, ``FermionicMeasure``) have no analogue here.
"""

from __future__ import annotations

import pytest
from qiskit import ClassicalRegister
from qiskit.circuit import Barrier, QuantumRegister
from qiskit.circuit.library import XGate
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
    OccupationB2QLayout,
)

SPEC = LevelSpec(1, "standard_binary")
"""One qubit per mode: a two-level truncation under the binary encoding."""


def mapper_fn(operator, num_qubits):
    """Maps with the one spec that the layout of these tests also uses."""
    return occupation_encoding(operator, SPEC, num_qubits)


def _hopping() -> BosonOperator:
    return BosonOperator.from_dict({((True, 0), (False, 1)): 1.0, ((True, 1), (False, 0)): 1.0})


def _pass_manager(synth: B2QSynthesis, layout=None) -> MultiStagePassManager:
    return MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        layout=layout if layout is not None else OccupationB2QLayout(SPEC),
        synthesis=synth,
        output=QuantumDAGToCircuit(),
    )


@pytest.mark.parametrize(
    ("config", "message"),
    [
        ({"Evolution": "NoSuchMethod"}, "No synthesis method 'NoSuchMethod'"),
        ({"Evolution": ("NoSuchMethod", ())}, "No synthesis method 'NoSuchMethod'"),
        ({"NotAnInstruction": "Anything"}, "No boson-to-qubit synthesis is installed"),
    ],
)
def test_an_unknown_instruction_or_method_is_rejected_when_the_pass_is_built(config, message):
    with pytest.raises(KeyError, match=message):
        B2QSynthesis(config)


@pytest.mark.parametrize("value", [42, None, ("MapperFn",), ["MapperFn"], object()])
def test_a_value_of_no_known_form_is_rejected(value):
    with pytest.raises(TypeError, match="must be a method name"):
        B2QSynthesis({"Evolution": value})  # type: ignore[dict-item]


def test_config_by_instance_needs_no_entry_point():
    plugin = MapperFnEvolutionSynthesis(mapper_fn)

    synth = B2QSynthesis({"Evolution": plugin})

    assert synth.methods["Evolution"] is plugin


def test_missing_plugin():
    """Test the handling of a missing boson-to-qubit plugin."""
    num_modes = 2
    circ = BosonicCircuit(num_modes)
    circ.append(Evolution(num_modes, _hopping(), time=1.5), circ.modes)

    with pytest.raises(TypeError, match="No plugin registered"):
        _ = _pass_manager(B2QSynthesis()).run(circ)


def test_config_tuple_with_positional_args():
    """The ``(name, args)`` config form forwards positional args to the plugin ``__init__``."""
    synth = B2QSynthesis({"Evolution": ("MapperFn", (mapper_fn,))})
    assert isinstance(synth.methods["Evolution"], MapperFnEvolutionSynthesis)


def test_config_tuple_with_positional_and_keyword_args():
    """The ``(name, args, kwargs)`` config form forwards both positional and keyword args."""
    synth = B2QSynthesis({"Evolution": ("MapperFn", (), {"mapper_fn": mapper_fn})})
    assert isinstance(synth.methods["Evolution"], MapperFnEvolutionSynthesis)


def test_config_bare_name_is_rejected_when_the_plugin_needs_arguments():
    """The bare-name form requires an ``__init__`` that takes no arguments."""
    with pytest.raises(TypeError):
        B2QSynthesis({"Evolution": "MapperFn"})


def test_run_rejects_non_bosonic_gate():
    """A circuit instruction that is not a BosonicGate raises during synthesis."""
    circ = BosonicCircuit(2)
    # Append a plain (non-bosonic) gate directly on the mode register.
    circ._inner.append(XGate(), [circ.register[0]])

    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)

    with pytest.raises(ValueError, match="unsupported circuit instruction type"):
        _ = _pass_manager(synth).run(circ)


def test_run_rejects_reset():
    """``Reset`` is non-unitary and not a bosonic gate, so it stays unsupported."""
    num_modes = 2
    circ = BosonicCircuit(num_modes)
    circ.append(Evolution(num_modes, _hopping(), time=0.5), circ.modes)
    # there is deliberately no public route for this, hence reaching through ``_inner``
    circ._inner.reset(circ.modes[0])

    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)

    with pytest.raises(ValueError, match="unsupported circuit instruction type"):
        _pass_manager(synth).run(circ)


def _barrier_circuit(num_modes: int) -> BosonicCircuit:
    """A circuit whose barrier covers only a subset of its modes."""
    circ = BosonicCircuit(num_modes)
    circ.append(Evolution(2, _hopping(), time=0.5), circ.modes[:2])
    circ.barrier(circ.register[0], circ.register[1])
    return circ


def _barrier_widths(circ: BosonicCircuit, layout) -> tuple[int, ...]:
    """Runs the synthesis stage and returns the width of every barrier in the output."""
    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)
    out = _pass_manager(synth, layout).run(circ)
    return tuple(
        instr.operation.num_qubits for instr in out.data if isinstance(instr.operation, Barrier)
    )


def test_run_lowers_barrier_to_mapped_register():
    """A bosonic barrier is lowered onto every qubit of the register it was mapped to.

    The barrier covers only modes 0 and 1, but the lowering deliberately widens it to the whole
    mapped register: a ``B2QLayout`` relates whole registers, so the modes it covers have no
    per-qubit image in general.
    """
    num_modes = 4
    assert _barrier_widths(_barrier_circuit(num_modes), OccupationB2QLayout(SPEC)) == (num_modes,)


def test_run_lowers_barrier_under_non_trivial_layout():
    """The barrier spans the mapped register even when it is sized differently from the modes.

    This pins the register-wide lowering: an implementation mapping mode ``i`` onto qubit ``i``
    would emit a width-2 barrier here (or index out of the register) instead of covering all six
    qubits.
    """
    num_modes, num_qubits = 4, 6
    circ = _barrier_circuit(num_modes)
    layout = CustomB2QLayout({circ.register: QuantumRegister(num_qubits)})

    assert _barrier_widths(circ, layout) == (num_qubits,)


def test_run_preserves_barrier_label():
    """A labeled bosonic barrier keeps its label through the synthesis stage."""
    circ = BosonicCircuit(2)
    circ.barrier(label="sync")

    out = _pass_manager(B2QSynthesis()).run(circ)

    (barrier,) = [instr for instr in out.data if isinstance(instr.operation, Barrier)]
    assert barrier.operation.label == "sync"


def test_classical_registers_are_preserved():
    """Classical registers on the input circuit are carried over to the synthesized output."""
    num_modes = 2
    circ = BosonicCircuit(num_modes)
    circ.append(Evolution(num_modes, _hopping(), time=0.5), circ.modes)
    circ._inner.add_register(ClassicalRegister(2, "c"))

    synth = B2QSynthesis()
    synth.methods["Evolution"] = MapperFnEvolutionSynthesis(mapper_fn)

    qu_circ = _pass_manager(synth).run(circ)
    assert any(reg.name == "c" and len(reg) == 2 for reg in qu_circ.cregs)
