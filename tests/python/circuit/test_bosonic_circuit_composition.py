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

"""Tests for ``BosonicCircuit.barrier``, ``compose`` and ``repeat``."""

from __future__ import annotations

import pytest
from qiskit.circuit import Barrier, CircuitError, QuantumCircuit
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes
from qiskit_bosons.circuit.library.synthesis import BosonicQDrift, BosonicSuzukiTrotter
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicPassManager
from qiskit_bosons.transpiler.passes import B2BSynthesis, RelabelModes
from qiskit_bosons.transpiler.presets import generate_preset_occupation_pass_manager

from ..trotter_support import (
    MAX_OCCUPATION,
    NUM_SITES,
    bose_hubbard,
    circuit_state,
    distance,
    exact_state,
)

TIME = 0.4


def _number(mode: int) -> BosonOperator:
    return BosonOperator.from_dict({(cre(mode), ann(mode)): 1.0})


def _names(circuit: BosonicCircuit) -> list[str]:
    return [inst.operation.name for inst in circuit._inner.data]


class TestBarrier:
    def test_default_spans_every_mode(self):
        circuit = BosonicCircuit(3)
        circuit.barrier()
        (inst,) = circuit._inner.data
        assert isinstance(inst.operation, Barrier)
        assert len(inst.qubits) == 3

    def test_spans_only_the_modes_given(self):
        circuit = BosonicCircuit(4)
        circuit.barrier(circuit.register[1], circuit.register[2], label="mid")
        (inst,) = circuit._inner.data
        assert inst.operation.label == "mid"
        assert list(inst.qubits) == [circuit.register[1], circuit.register[2]]

    def test_append_still_refuses_non_bosonic_instructions(self):
        """A barrier is added through its own method; ``append`` stays a type guard."""
        circuit = BosonicCircuit(2)
        with pytest.raises(ValueError, match="Unsupported instruction type"):
            circuit.append(Barrier(2), circuit.modes)  # type: ignore[arg-type]

    def test_the_error_message_is_formatted(self):
        """It used to be passed ``%s`` and an argument as separate exception arguments."""
        circuit = BosonicCircuit(1)
        with pytest.raises(ValueError, match="Unsupported instruction type: <class 'str'>"):
            circuit.append("not a gate", [0])  # type: ignore[arg-type]

    def test_barrier_does_not_count_as_a_gate_but_is_counted(self):
        circuit = BosonicCircuit(2)
        circuit.append(InitializeModes([1, 0]), circuit.modes)
        circuit.barrier()
        assert dict(circuit.count_ops()) == {"InitializeModes": 1, "barrier": 1}

    def test_decompose_keeps_the_barrier(self):
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, bose_hubbard(), time=TIME), circuit.modes)
        circuit.barrier()
        assert circuit.decompose().count_ops()["barrier"] == 1


class TestCompose:
    def test_returns_a_new_circuit_by_default(self):
        left, right = BosonicCircuit(2), BosonicCircuit(2)
        left.append(InitializeModes([1, 0]), left.modes)
        right.append(Evolution(2, _number(0), time=TIME), right.modes)

        out = left.compose(right)

        assert isinstance(out, BosonicCircuit)
        assert out is not left
        assert _names(out) == ["InitializeModes", "Evolution"]
        # the operands are untouched
        assert _names(left) == ["InitializeModes"]
        assert _names(right) == ["Evolution"]

    def test_inplace_modifies_and_returns_none(self):
        left, right = BosonicCircuit(2), BosonicCircuit(2)
        left.append(InitializeModes([1, 0]), left.modes)
        right.append(Evolution(2, _number(0), time=TIME), right.modes)
        assert left.compose(right, inplace=True) is None
        assert _names(left) == ["InitializeModes", "Evolution"]

    def test_front_composes_before_the_existing_instructions(self):
        left, right = BosonicCircuit(2), BosonicCircuit(2)
        left.append(Evolution(2, _number(0), time=TIME), left.modes)
        right.append(InitializeModes([1, 0]), right.modes)
        assert _names(left.compose(right, front=True)) == ["InitializeModes", "Evolution"]

    def test_places_a_narrower_circuit_onto_chosen_modes(self):
        wide, narrow = BosonicCircuit(3), BosonicCircuit(1)
        narrow.append(Evolution(1, _number(0), time=TIME), narrow.modes)
        out = wide.compose(narrow, [2])
        (inst,) = out._inner.data
        assert list(inst.qubits) == [wide.register[2]]

    def test_the_result_keeps_this_circuits_register(self):
        left, right = BosonicCircuit(2), BosonicCircuit(2)
        out = left.compose(right)
        assert out.register is left.register
        assert out.modes == left.modes

    def test_inlines_instead_of_wrapping(self):
        left, right = BosonicCircuit(2), BosonicCircuit(2)
        right.append(InitializeModes([1, 0]), right.modes)
        right.append(Evolution(2, _number(0), time=TIME), right.modes)
        assert dict(left.compose(right).count_ops()) == {"InitializeModes": 1, "Evolution": 1}

    def test_a_barrier_is_carried_over(self):
        left, right = BosonicCircuit(2), BosonicCircuit(2)
        right.barrier()
        assert _names(left.compose(right)) == ["barrier"]

    def test_rejects_a_plain_quantum_circuit(self):
        """Inlining qubit gates would bypass the type guard of ``append``."""
        circuit = BosonicCircuit(2)
        qubit_circuit = QuantumCircuit(2)
        qubit_circuit.h(0)
        with pytest.raises(ValueError, match="Unsupported circuit type"):
            circuit.compose(qubit_circuit)  # type: ignore[arg-type]

    def test_rejects_a_wider_operand(self):
        with pytest.raises(CircuitError):
            BosonicCircuit(1).compose(BosonicCircuit(2))


class TestRepeat:
    @staticmethod
    def _step() -> BosonicCircuit:
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, bose_hubbard(), time=TIME), circuit.modes)
        return circuit

    def test_flattens_the_repetitions(self):
        out = self._step().repeat(4)
        assert isinstance(out, BosonicCircuit)
        assert dict(out.count_ops()) == {"Evolution": 4}

    def test_decompose_acts_on_the_repeated_body(self):
        """A wrapped repetition would peel off one layer; a flattened one decomposes the body."""
        groups = bose_hubbard().num_groups()
        assert self._step().repeat(3).decompose().count_ops()["Evolution"] == 3 * groups

    def test_zero_repetitions_is_the_identity(self):
        out = self._step().repeat(0)
        assert len(out._inner.data) == 0
        assert len(out.modes) == NUM_SITES

    def test_one_repetition_copies(self):
        step = self._step()
        out = step.repeat(1)
        assert out is not step
        assert _names(out) == _names(step)

    def test_negative_repetitions_are_rejected(self):
        with pytest.raises(ValueError, match="non-negative"):
            self._step().repeat(-1)

    def test_the_original_is_untouched(self):
        step = self._step()
        step.repeat(5)
        assert dict(step.count_ops()) == {"Evolution": 1}

    def test_barriers_go_between_repetitions_only(self):
        out = self._step().repeat(3, insert_barriers=True)
        assert _names(out) == ["Evolution", "barrier", "Evolution", "barrier", "Evolution"]

    def test_no_barrier_for_a_single_repetition(self):
        assert _names(self._step().repeat(1, insert_barriers=True)) == ["Evolution"]

    def test_insert_barriers_is_keyword_only(self):
        with pytest.raises(TypeError):
            self._step().repeat(2, True)  # type: ignore[misc]

    def test_metadata_and_register_are_kept(self):
        step = self._step()
        step.metadata = {"label": "step"}
        out = step.repeat(2)
        assert out.metadata == {"label": "step"}
        assert out.register is step.register


class TestBarrierThroughTheTranspiler:
    """A barrier must survive every pass and land on all qubits of the registers it touched."""

    def _circuit(self) -> BosonicCircuit:
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(InitializeModes([1, 0, 1]), circuit.modes)
        circuit.barrier()
        circuit.append(Evolution(NUM_SITES, bose_hubbard(), time=TIME), circuit.modes)
        return circuit

    def test_the_barrier_spans_every_qubit_of_the_mapped_register(self):
        spec = LevelSpec(MAX_OCCUPATION, "unary_projected")
        lowered = generate_preset_occupation_pass_manager(spec).run(self._circuit())
        barriers = [i for i in lowered.data if isinstance(i.operation, Barrier)]
        assert len(barriers) == 1
        # not narrowed to the qubits of the modes it was placed on: the whole encoded register
        assert len(barriers[0].qubits) == lowered.num_qubits
        assert lowered.num_qubits == spec.resolve(NUM_SITES).total_width

    def test_a_partial_barrier_still_covers_the_whole_register(self):
        """A layout relates whole registers, so there is no per-qubit image of a single mode."""
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(InitializeModes([1, 0, 1]), circuit.modes)
        circuit.barrier(circuit.register[0])
        lowered = generate_preset_occupation_pass_manager(
            LevelSpec(MAX_OCCUPATION, "unary_projected")
        ).run(circuit)
        (barrier,) = [i for i in lowered.data if isinstance(i.operation, Barrier)]
        assert len(barrier.qubits) == lowered.num_qubits

    def test_the_label_survives(self):
        circuit = BosonicCircuit(NUM_SITES)
        circuit.barrier(label="sync")
        lowered = generate_preset_occupation_pass_manager(
            LevelSpec(MAX_OCCUPATION, "unary_projected")
        ).run(circuit)
        (barrier,) = [i for i in lowered.data if isinstance(i.operation, Barrier)]
        assert barrier.operation.label == "sync"

    def test_a_barrier_does_not_change_the_physics(self):
        operator = bose_hubbard()
        with_barrier = self._circuit()
        without = BosonicCircuit(NUM_SITES)
        without.append(InitializeModes([1, 0, 1]), without.modes)
        without.append(Evolution(NUM_SITES, operator, time=TIME), without.modes)
        assert distance(circuit_state(with_barrier), circuit_state(without)) < 1e-12
        assert distance(circuit_state(with_barrier), exact_state(operator, TIME)) < 1e-6

    def test_relabel_modes_carries_the_barrier_along(self):
        out = BosonicPassManager([RelabelModes([2, 0, 1])]).run(self._circuit())
        assert out.count_ops()["barrier"] == 1

    def test_qdrift_copies_the_barrier_through(self):
        out = BosonicPassManager([B2BSynthesis({"Evolution": BosonicQDrift(5, rng=0)})]).run(
            self._circuit()
        )
        assert out.count_ops()["barrier"] == 1
        assert out.count_ops()["Evolution"] == 5

    def test_trotterization_leaves_the_barrier_alone(self):
        out = BosonicPassManager([B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=2)})]).run(
            self._circuit()
        )
        assert out.count_ops()["barrier"] == 1
        assert out.count_ops()["Evolution"] == 2 * bose_hubbard().num_groups() - 1

    def test_a_repeated_trotter_step_converges_to_the_exact_evolution(self):
        """The circuit-level way to reduce the outer error: repeat one short step."""
        operator = bose_hubbard()
        total = 0.8
        errors = []
        for steps in (1, 2, 4):
            repeated = BosonicCircuit(NUM_SITES)
            repeated.append(InitializeModes([1, 0, 1]), repeated.modes)
            body = BosonicCircuit(NUM_SITES)
            body.append(Evolution(NUM_SITES, operator, time=total / steps), body.modes)
            repeated.compose(body.repeat(steps), inplace=True)
            errors.append(
                distance(circuit_state(repeated.decompose()), exact_state(operator, total))
            )
        assert errors == sorted(errors, reverse=True)
        assert errors[0] / errors[-1] > 2.5
