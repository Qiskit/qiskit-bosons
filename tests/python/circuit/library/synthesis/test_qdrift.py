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

"""Tests for the ``weights`` argument of the qDRIFT synthesis method."""

from __future__ import annotations

import numpy as np
import pytest
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution
from qiskit_bosons.circuit.library.synthesis import BosonicQDrift
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.operators.terms import group_coeff_means, group_order
from qiskit_bosons.transpiler import BosonicPassManager
from qiskit_bosons.transpiler.passes import B2BSynthesis

from ....trotter_support import NUM_SITES, bose_hubbard, evolution_circuit

TIME = 0.8


def _sampled(circuit: BosonicCircuit, **kwargs) -> list[Evolution]:
    out = BosonicPassManager([B2BSynthesis({"Evolution": BosonicQDrift(**kwargs)})]).run(circuit)
    return [i.operation for i in out._inner.data if isinstance(i.operation, Evolution)]


def _terms(gates: list[Evolution]) -> list:
    """The terms of every sampled gate, as a comparable structure (term order is not guaranteed)."""
    return [sorted(map(repr, gate.operator.iter_terms())) for gate in gates]


def _ungrouped() -> BosonOperator:
    return BosonOperator.from_dict(
        {(cre(0), ann(0)): 1.0, (cre(1), ann(1)): 3.0, (cre(2), cre(2), ann(2), ann(2)): 6.0}
    )


class TestValidation:
    @pytest.mark.parametrize(
        ("weights", "message"),
        [
            ([[1.0, 2.0]], "one-dimensional"),
            ([], "must not be empty"),
            ([1.0, np.nan], "finite"),
            ([1.0, np.inf], "finite"),
            ([1.0, -1.0], "Negative sampling weights"),
            ([0.0, 0.0], "sum to zero"),
        ],
    )
    def test_rejects_malformed_weights(self, weights, message):
        with pytest.raises(ValueError, match=message):
            BosonicQDrift(4, weights=weights)

    def test_weights_default_to_none_and_are_stored_as_an_array(self):
        assert BosonicQDrift(4).weights is None
        stored = BosonicQDrift(4, weights=[1, 2, 3]).weights
        assert isinstance(stored, np.ndarray)
        assert stored.dtype == float

    def test_weights_is_keyword_only(self):
        with pytest.raises(TypeError):
            BosonicQDrift(4, [1.0, 2.0])  # type: ignore[misc]

    def test_length_must_match_the_number_of_terms(self):
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, _ungrouped(), time=TIME), circuit.modes)
        with pytest.raises(ValueError, match="3 entries for an operator with 3 terms, but got 2"):
            _sampled(circuit, num_terms=4, weights=[1.0, 2.0])

    def test_length_must_match_the_number_of_groups_when_grouped(self):
        operator = bose_hubbard()
        wrong = [1.0] * (operator.num_groups() + 1)
        with pytest.raises(ValueError, match=f"{operator.num_groups()} groups"):
            _sampled(evolution_circuit(operator, TIME), num_terms=4, weights=wrong)

    def test_a_weights_array_is_checked_against_every_gate(self):
        """The method is stateless, so it cannot count the gates it is given; it validates the
        length against each of them instead, which rejects a second gate of another size."""
        operator = _ungrouped()
        other = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, operator, time=TIME), circuit.modes)
        circuit.append(Evolution(NUM_SITES, other, time=TIME), circuit.modes)
        with pytest.raises(ValueError, match="exactly one entry per sampled piece"):
            _sampled(circuit, num_terms=4, weights=[1.0, 3.0, 6.0])

    def test_without_weights_several_gates_are_supported(self):
        """Each gate derives its own weights."""
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, _ungrouped(), time=TIME), circuit.modes)
        circuit.append(Evolution(NUM_SITES, _ungrouped(), time=TIME), circuit.modes)
        assert len(_sampled(circuit, num_terms=4, rng=1)) == 8


class TestWeightsOnTheGate:
    """The method can be the ``synthesis`` of a gate, which keeps weights next to their operator."""

    def test_each_gate_carries_the_weights_of_its_own_operator(self):
        first = _ungrouped()
        second = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(
            Evolution(
                NUM_SITES,
                first,
                time=TIME,
                synthesis=BosonicQDrift(4, weights=[1.0, 3.0, 6.0], rng=1),
            ),
            circuit.modes,
        )
        circuit.append(
            Evolution(
                NUM_SITES, second, time=TIME, synthesis=BosonicQDrift(5, weights=[1.0, 1.0], rng=1)
            ),
            circuit.modes,
        )

        gates = [
            i.operation
            for i in circuit.decompose()._inner.data
            if isinstance(i.operation, Evolution)
        ]

        assert len(gates) == 4 + 5
        assert all(gate.atomic for gate in gates)

    def test_a_gate_rejects_weights_of_another_operator(self):
        gate = Evolution(
            NUM_SITES, _ungrouped(), time=TIME, synthesis=BosonicQDrift(4, weights=[1.0, 2.0])
        )
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(gate, circuit.modes)
        with pytest.raises(ValueError, match="exactly one entry per sampled piece"):
            circuit.decompose()


class TestSampling:
    def test_default_weights_equal_the_coefficient_magnitudes(self):
        """Passing the magnitudes explicitly reproduces the default draw exactly."""
        operator = _ungrouped()
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, operator, time=TIME), circuit.modes)

        default = _sampled(circuit, num_terms=12, rng=7)
        explicit = _sampled(circuit, num_terms=12, weights=np.abs(operator.get_coeffs()), rng=7)
        assert _terms(default) == _terms(explicit)
        assert [g.params[0] for g in default] == pytest.approx([g.params[0] for g in explicit])

    def test_grouped_default_equals_group_coeff_means(self):
        operator = group_order(bose_hubbard())
        circuit = evolution_circuit(operator, TIME)
        default = _sampled(circuit, num_terms=12, rng=3)
        explicit = _sampled(circuit, num_terms=12, weights=group_coeff_means(operator), rng=3)
        assert [len(g.operator) for g in default] == [len(g.operator) for g in explicit]
        assert [g.params[0] for g in default] == pytest.approx([g.params[0] for g in explicit])

    def test_zero_weight_entries_are_never_sampled(self):
        operator = _ungrouped()
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, operator, time=TIME), circuit.modes)
        # give weight only to the number operator of mode 1, found via the *stored* term order
        keep = next(
            i for i, (term, _) in enumerate(operator.iter_terms()) if {m for _, m in term} == {1}
        )
        weights = [1.0 if i == keep else 0.0 for i in range(len(operator))]
        gates = _sampled(circuit, num_terms=40, weights=weights, rng=0)
        # every sampled gate is that single term (its unit-magnitude, sign-preserving version)
        assert len(gates) == 40
        assert len({tuple(map(repr, g.operator.iter_terms())) for g in gates}) == 1

    def test_the_scale_of_the_weights_sets_the_evolution_time(self):
        """Rescaling every weight by gamma samples identically but evolves for gamma * t."""
        operator = _ungrouped()
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, operator, time=TIME), circuit.modes)
        base = np.abs(operator.get_coeffs())

        times = _sampled(circuit, num_terms=5, weights=base, rng=11)
        scaled = _sampled(circuit, num_terms=5, weights=3.0 * base, rng=11)

        assert _terms(times) == _terms(scaled)
        assert [g.params[0] for g in scaled] == pytest.approx([3.0 * g.params[0] for g in times])

    def test_the_total_evolution_time_matches_lambda_t(self):
        """With ``w_j = |c_j|`` the sampled gates together evolve for ``lambda * t``."""
        operator = _ungrouped()
        circuit = BosonicCircuit(NUM_SITES)
        circuit.append(Evolution(NUM_SITES, operator, time=TIME), circuit.modes)
        gates = _sampled(circuit, num_terms=10, rng=5)
        assert sum(g.params[0] for g in gates) == pytest.approx(10.0 * TIME)  # lambda = 10

    def test_signs_are_kept_on_the_sampled_operator(self):
        operator = BosonOperator.from_dict({(cre(0), ann(0)): -2.0, (cre(1), ann(1)): 2.0})
        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, operator, time=TIME), circuit.modes)
        for gate in _sampled(circuit, num_terms=20, rng=2):
            ((_, coeff),) = list(gate.operator.iter_terms())
            assert abs(coeff) == pytest.approx(1.0)

    def test_a_fixed_seed_reproduces_the_sequence_of_circuits(self):
        circuit = evolution_circuit(bose_hubbard(), TIME)
        runs = []
        for _ in range(2):
            pass_ = BosonicPassManager([B2BSynthesis({"Evolution": BosonicQDrift(6, rng=42)})])
            runs.append([pass_.run(circuit).count_ops() for _ in range(3)])
        assert runs[0] == runs[1]


class TestSampler:
    def _gate(self) -> Evolution:
        return Evolution(3, _ungrouped(), time=TIME)

    def test_synthesize_equals_the_batched_sample(self):
        gate = self._gate()
        circuit = BosonicQDrift(9, rng=5).synthesize(gate)
        sampled = BosonicQDrift(9, rng=5).sampler(gate).sample()

        assert _terms([i.operation for i in circuit._inner.data]) == _terms(sampled)

    def test_sample_size_defaults_to_num_terms_and_can_be_overridden(self):
        sampler = BosonicQDrift(4, rng=1).sampler(self._gate())

        assert len(sampler.sample()) == 4
        assert len(sampler.sample(7)) == 7

    def test_sampled_gates_are_atomic_with_the_shared_time(self):
        sampler = BosonicQDrift(4, rng=1).sampler(self._gate())

        for gate in sampler.sample():
            assert gate.atomic
            assert float(gate.params[0]) == pytest.approx(sampler.time)
        assert sampler.time == pytest.approx(10.0 * TIME / 4)

    def test_draw_follows_the_distribution(self):
        gate = self._gate()
        # the term order of `from_dict` is not guaranteed, so find the term acting on mode 1
        weights = [
            0.0 if {action[1] for action in actions} == {1} else 1.0
            for actions, _ in gate.operator.iter_terms()
        ]
        assert sorted(weights) == [0.0, 1.0, 1.0]
        sampler = BosonicQDrift(5, rng=3, weights=weights).sampler(gate)

        supports = {tuple(sorted(sampler.draw().operator.get_support())) for _ in range(50)}

        # the zero-weight term (mode 1) is never drawn
        assert supports == {(0,), (2,)}

    def test_draw_works_for_a_grouped_operator(self):
        operator = group_order(
            BosonOperator.from_terms_with_groups(
                [((cre(0), ann(1)), 0.5, 0), ((cre(1), ann(0)), 0.5, 0), ((cre(2), ann(2)), 1.0, 1)]
            )
        )
        sampler = BosonicQDrift(3, rng=0).sampler(Evolution(3, operator, time=TIME))

        gate = sampler.draw()

        assert gate.atomic
        assert len(gate.operator) in (1, 2)

    def test_sampler_rejects_weights_of_the_wrong_length(self):
        with pytest.raises(ValueError, match="exactly one entry per sampled piece"):
            BosonicQDrift(3, weights=[1.0, 1.0]).sampler(self._gate())
