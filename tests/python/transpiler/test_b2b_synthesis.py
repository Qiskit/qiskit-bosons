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

"""Tests for the ``B2BSynthesis`` pass and the ``Evolution`` synthesis machinery."""

from __future__ import annotations

import pytest
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes
from qiskit_bosons.circuit.library.synthesis import (
    BosonicLieTrotter,
    BosonicQDrift,
    BosonicSuzukiTrotter,
)
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicPassManager
from qiskit_bosons.transpiler.passes import B2BSynthesis

from ..trotter_support import (
    NUM_SITES,
    bose_hubbard,
    circuit_state,
    distance,
    evolution_circuit,
    exact_state,
)

TIME = 0.8


def _run(circuit: BosonicCircuit, *passes) -> BosonicCircuit:
    return BosonicPassManager(list(passes)).run(circuit)


class TestEvolutionGate:
    def test_defaults(self):
        gate = Evolution(NUM_SITES, bose_hubbard(), time=TIME)
        assert isinstance(gate.synthesis, BosonicLieTrotter)
        assert gate.atomic is False

    def test_synthesis_and_atomic_are_keyword_only_and_read_only(self):
        with pytest.raises(TypeError):
            Evolution(NUM_SITES, bose_hubbard(), TIME, BosonicLieTrotter())  # type: ignore[misc]
        gate = Evolution(NUM_SITES, bose_hubbard(), time=TIME)
        with pytest.raises(AttributeError):
            gate.synthesis = BosonicLieTrotter()  # type: ignore[misc]
        with pytest.raises(AttributeError):
            gate.atomic = True  # type: ignore[misc]

    def test_the_chosen_synthesis_decides_the_decomposition(self):
        operator = bose_hubbard()
        groups = operator.num_groups()
        plain = evolution_circuit(operator, TIME).decompose()
        suzuki = evolution_circuit(
            operator, TIME, synthesis=BosonicSuzukiTrotter(order=2)
        ).decompose()
        assert plain.count_ops()["Evolution"] == groups
        assert suzuki.count_ops()["Evolution"] == 2 * groups - 1

    def test_atomic_gates_are_left_in_place_by_decompose(self):
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, bose_hubbard(), time=TIME, atomic=True), circuit.modes)
        assert dict(circuit.decompose().count_ops()) == {"Evolution": 1}
        # ... however often it is repeated: an atomic gate has no definition to expand
        assert dict(circuit.decompose(reps=4).count_ops()) == {"Evolution": 1}

    def test_decomposing_repeatedly_reaches_a_fixed_point(self):
        """The factors a synthesis emits are atomic, so a further decomposition changes nothing."""
        circuit = evolution_circuit(bose_hubbard(), TIME)
        once = circuit.decompose(gates_to_decompose="Evolution")
        assert once.decompose(gates_to_decompose="Evolution").count_ops() == once.count_ops()

    def test_a_single_term_evolution_does_not_decompose_into_itself(self):
        """Used to recurse without end, e.g. through ``inverse()``."""
        operator = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        gate = Evolution(1, operator, time=TIME)
        gate.inverse().inverse()  # must not raise RecursionError
        circuit = BosonicCircuit(1)
        circuit.append(gate, circuit.modes)
        assert circuit.decompose(reps=6).count_ops()["Evolution"] == 1

    def test_inverse_keeps_operator_synthesis_and_atomicity(self):
        synthesis = BosonicSuzukiTrotter(order=2, reps=3)
        gate = Evolution(NUM_SITES, bose_hubbard(), time=TIME, synthesis=synthesis, atomic=True)
        inverse = gate.inverse()
        assert isinstance(inverse, Evolution)
        assert inverse.params[0] == pytest.approx(-TIME)
        assert inverse.synthesis is synthesis
        assert inverse.atomic is True
        assert inverse.operator is gate.operator
        assert inverse.inverse().params[0] == pytest.approx(TIME)

    def test_inverse_undoes_the_evolution(self):
        operator = bose_hubbard()
        circuit = evolution_circuit(operator, TIME)
        circuit.append(Evolution(NUM_SITES, operator, time=TIME).inverse(), circuit.modes)
        # evolve and un-evolve: back at the initial occupation (compared via the exact evolution)
        assert distance(circuit_state(circuit), exact_state(operator, 0.0)) < 1e-6


class TestB2BSynthesis:
    def test_expands_every_evolution_with_the_chosen_formula(self):
        operator = bose_hubbard()
        circuit = evolution_circuit(operator, TIME)
        out = _run(circuit, B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=2)}))
        assert out.count_ops()["Evolution"] == 2 * operator.num_groups() - 1
        assert out.count_ops()["InitializeModes"] == 1

    def test_the_expanded_factors_are_atomic(self):
        out = _run(
            evolution_circuit(bose_hubbard(), TIME),
            B2BSynthesis({"Evolution": BosonicLieTrotter()}),
        )
        gates = [
            inst.operation for inst in out._inner.data if isinstance(inst.operation, Evolution)
        ]
        assert gates
        assert all(gate.atomic for gate in gates)

    def test_matches_setting_the_synthesis_on_the_gate(self):
        operator = bose_hubbard()
        synthesis = BosonicSuzukiTrotter(order=2, reps=2)
        via_pass = _run(evolution_circuit(operator, TIME), B2BSynthesis({"Evolution": synthesis}))
        via_gate = evolution_circuit(operator, TIME, synthesis=synthesis).decompose()
        assert distance(circuit_state(via_pass), circuit_state(via_gate)) < 1e-12

    def test_a_higher_order_is_more_accurate_than_the_default(self):
        operator = bose_hubbard()
        exact = exact_state(operator, TIME)
        first = _run(
            evolution_circuit(operator, TIME), B2BSynthesis({"Evolution": BosonicLieTrotter()})
        )
        second = _run(
            evolution_circuit(operator, TIME),
            B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=4)}),
        )
        assert distance(circuit_state(second), exact) < distance(circuit_state(first), exact)

    def test_running_it_again_is_a_no_op(self):
        pass_ = B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=2)})
        once = _run(evolution_circuit(bose_hubbard(), TIME), pass_)
        twice = _run(once, pass_)
        assert twice.count_ops() == once.count_ops()

    def test_apply_false_only_selects_the_method(self):
        operator = bose_hubbard()
        synthesis = BosonicSuzukiTrotter(order=2)
        out = _run(
            evolution_circuit(operator, TIME), B2BSynthesis({"Evolution": synthesis}, apply=False)
        )
        # still one gate, but now carrying the chosen method ...
        assert out.count_ops()["Evolution"] == 1
        gate = next(i.operation for i in out._inner.data if isinstance(i.operation, Evolution))
        assert gate.synthesis is synthesis
        # ... which a later decomposition then honors
        assert out.decompose().count_ops()["Evolution"] == 2 * operator.num_groups() - 1

    def test_the_pass_leaves_atomic_gates_alone(self):
        synthesis = BosonicSuzukiTrotter(order=2)
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, bose_hubbard(), time=TIME, atomic=True), circuit.modes)
        out = _run(circuit, B2BSynthesis({"Evolution": synthesis}))
        gate = next(i.operation for i in out._inner.data if isinstance(i.operation, Evolution))
        assert gate.atomic
        # it passes through unchanged: it is not even tagged with the method
        assert gate.synthesis is not synthesis

    @pytest.mark.parametrize("apply", [True, False])
    def test_the_pass_does_not_tag_atomic_gates(self, apply):
        atomic = Evolution(NUM_SITES, bose_hubbard(), time=TIME, atomic=True)
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(atomic, circuit.modes)
        out = _run(circuit, B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=2)}, apply=apply))
        gate = next(i.operation for i in out._inner.data if isinstance(i.operation, Evolution))
        assert gate.synthesis is atomic.synthesis

    def test_the_filter_is_not_consulted_for_atomic_gates(self):
        seen = []
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, bose_hubbard(), time=TIME, atomic=True), circuit.modes)
        _run(
            circuit,
            B2BSynthesis(
                {"Evolution": BosonicSuzukiTrotter(order=2)}, filter=lambda node: seen.append(node)
            ),
        )
        assert seen == []

    def test_the_filter_restricts_the_gates_that_are_expanded(self):
        grouped, diagonal = bose_hubbard(), BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, grouped, time=TIME), circuit.modes)
        circuit.append(Evolution(NUM_SITES, diagonal, time=TIME), circuit.modes)

        pass_ = B2BSynthesis(
            {"Evolution": BosonicSuzukiTrotter(order=2)},
            filter=lambda node: node.op.operator.has_groups(),
        )
        out = _run(circuit, pass_)
        # the grouped gate expands, the other one is untouched
        assert out.count_ops()["Evolution"] == (2 * grouped.num_groups() - 1) + 1

    def test_non_evolution_instructions_are_untouched(self):
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(InitializeModes([1, 0, 1]), circuit.modes)
        out = _run(circuit, B2BSynthesis({"Evolution": BosonicLieTrotter()}))
        assert dict(out.count_ops()) == {"InitializeModes": 1}

    def test_the_input_gate_is_not_mutated(self):
        """A gate instance may be shared between circuits, so the pass must not edit it in place."""
        gate = Evolution(NUM_SITES, bose_hubbard(), time=TIME)
        original = gate.synthesis
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(gate, circuit.modes, copy=False)
        _run(circuit, B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=2)}, apply=False))
        assert gate.synthesis is original

    def test_ungrouped_operators_with_the_number_operator_are_expandable(self):
        """A lone ``b†b`` is Hermitian, so even a term-by-term split is a valid factor."""
        operator = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (cre(1), ann(1)): 2.0})
        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, operator, time=TIME), circuit.modes)
        out = _run(circuit, B2BSynthesis({"Evolution": BosonicLieTrotter()}))
        assert out.count_ops()["Evolution"] == 2


class TestQDriftIntegration:
    def test_sampled_gates_are_atomic_and_keep_the_synthesis_method(self):
        synthesis = BosonicSuzukiTrotter(order=2, reps=2)
        circuit = evolution_circuit(bose_hubbard(), TIME, synthesis=synthesis)
        out = _run(circuit, B2BSynthesis({"Evolution": BosonicQDrift(7, rng=1)}))
        gates = [i.operation for i in out._inner.data if isinstance(i.operation, Evolution)]
        assert len(gates) == 7
        assert all(gate.atomic for gate in gates)
        assert all(gate.synthesis is synthesis for gate in gates)

    def test_trotterization_leaves_sampled_gates_alone(self):
        circuit = evolution_circuit(bose_hubbard(), TIME)
        sampled = _run(circuit, B2BSynthesis({"Evolution": BosonicQDrift(7, rng=1)}))
        again = _run(sampled, B2BSynthesis({"Evolution": BosonicSuzukiTrotter(order=2)}))
        assert again.count_ops() == sampled.count_ops()

    def test_decompose_does_not_split_sampled_gates(self):
        sampled = _run(
            evolution_circuit(bose_hubbard(), TIME),
            B2BSynthesis({"Evolution": BosonicQDrift(7, rng=1)}),
        )
        assert sampled.decompose(reps=3).count_ops() == sampled.count_ops()
