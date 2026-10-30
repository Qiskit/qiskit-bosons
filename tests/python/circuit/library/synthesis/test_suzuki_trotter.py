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

"""Tests for the bosonic Lie- and Suzuki-Trotter product formulas."""

from __future__ import annotations

from collections import defaultdict
from itertools import pairwise

import pytest
from qiskit_bosons.circuit.library import Evolution
from qiskit_bosons.circuit.library.synthesis import (
    B2BSynthesisPlugin,
    B2BSynthesisPluginManager,
    BosonicLieTrotter,
    BosonicSuzukiTrotter,
)
from qiskit_bosons.circuit.library.synthesis.suzuki_trotter import _suzuki_schedule
from qiskit_bosons.operators import BosonOperator, ann, cre

from ....trotter_support import (
    bose_hubbard,
    circuit_state,
    distance,
    evolution_circuit,
    exact_state,
)

TIME = 0.8


def _factors(synthesis: B2BSynthesisPlugin, operator: BosonOperator, num_modes: int = 3):
    """The ``(operator, time, qargs)`` of every gate a synthesis emits."""
    gate = Evolution(num_modes, operator, time=TIME)
    circuit = synthesis.synthesize(gate)
    return [(inst.operation, inst.qubits) for inst in circuit._inner.data]


class TestValidation:
    @pytest.mark.parametrize("order", [0, 3, 5, -2])
    def test_rejects_odd_and_non_positive_orders(self, order):
        with pytest.raises(ValueError, match="order of 1 or a positive even order"):
            BosonicSuzukiTrotter(order=order)

    @pytest.mark.parametrize("reps", [0, -1])
    def test_rejects_non_positive_reps(self, reps):
        with pytest.raises(ValueError, match="must be positive"):
            BosonicSuzukiTrotter(order=2, reps=reps)
        with pytest.raises(ValueError, match="must be positive"):
            BosonicLieTrotter(reps=reps)

    def test_attributes(self):
        synthesis = BosonicSuzukiTrotter(order=4, reps=3)
        assert (synthesis.order, synthesis.reps) == (4, 3)

    def test_lie_trotter_is_the_first_order_member(self):
        lie = BosonicLieTrotter(reps=2)
        assert isinstance(lie, BosonicSuzukiTrotter)
        assert (lie.order, lie.reps) == (1, 2)

    def test_is_registered_as_an_entry_point(self):
        manager = B2BSynthesisPluginManager()
        assert manager.method("Evolution", "SuzukiTrotter") is BosonicSuzukiTrotter
        assert manager.method("Evolution", "LieTrotter") is BosonicLieTrotter


class TestSchedule:
    def test_first_order_is_one_sweep_in_order(self):
        assert _suzuki_schedule(1, 3, 1) == [(0, 1.0), (1, 1.0), (2, 1.0)]

    def test_second_order_is_a_palindrome_with_halved_outer_factors(self):
        assert _suzuki_schedule(2, 3, 1) == [(0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)]

    def test_repetitions_divide_the_time(self):
        schedule = _suzuki_schedule(1, 2, 4, merge_adjacent=False)
        assert len(schedule) == 8
        assert all(scale == pytest.approx(0.25) for _, scale in schedule)

    def test_first_order_repetitions_alternate_and_never_touch(self):
        """The sweep ends on the factor the next one starts with only when there is one factor."""
        merged = _suzuki_schedule(1, 2, 4)
        assert [index for index, _ in merged] == [0, 1] * 4

    def test_repetition_boundaries_of_a_palindrome_are_merged(self):
        """The palindromes of neighboring repetitions end and begin on the same factor."""
        doubled = _suzuki_schedule(2, 3, 2)
        assert [index for index, _ in doubled] == [0, 1, 2, 1, 0, 1, 2, 1, 0]
        assert [scale for _, scale in doubled] == [
            0.25,
            0.25,
            0.5,
            0.25,
            0.5,
            0.25,
            0.5,
            0.25,
            0.25,
        ]

    def test_unmerged_repetitions_repeat_the_sweep(self):
        single = _suzuki_schedule(2, 3, 1, merge_adjacent=False)
        doubled = _suzuki_schedule(2, 3, 2, merge_adjacent=False)
        assert len(doubled) == 2 * len(single)
        assert [index for index, _ in doubled] == [index for index, _ in single] * 2
        assert [scale for _, scale in doubled] == [scale / 2 for _, scale in single] * 2

    @pytest.mark.parametrize("order", [1, 2, 4, 6])
    @pytest.mark.parametrize("reps", [1, 3])
    def test_merging_never_changes_what_each_factor_accumulates(self, order, reps):
        merged = _suzuki_schedule(order, 3, reps)
        unmerged = _suzuki_schedule(order, 3, reps, merge_adjacent=False)
        assert len(merged) <= len(unmerged)
        for factor in range(3):
            assert sum(s for i, s in merged if i == factor) == pytest.approx(
                sum(s for i, s in unmerged if i == factor)
            )

    @pytest.mark.parametrize("order", [1, 2, 4, 6])
    @pytest.mark.parametrize("reps", [1, 3])
    def test_a_merged_schedule_has_no_neighbors_of_the_same_factor(self, order, reps):
        indices = [index for index, _ in _suzuki_schedule(order, 3, reps)]
        assert all(a != b for a, b in pairwise(indices))

    def test_merge_adjacent_defaults_to_enabled_and_is_stored(self):
        assert BosonicSuzukiTrotter().merge_adjacent is True
        assert BosonicSuzukiTrotter(merge_adjacent=False).merge_adjacent is False

    @pytest.mark.parametrize("order", [1, 2, 4, 6])
    @pytest.mark.parametrize("reps", [1, 3])
    @pytest.mark.parametrize("num_factors", [2, 3, 5])
    def test_each_factor_always_evolves_for_the_whole_time(self, order, reps, num_factors):
        """The scales of any one factor sum to one, whatever the order."""
        totals: dict[int, float] = defaultdict(float)
        for index, scale in _suzuki_schedule(order, num_factors, reps):
            totals[index] += scale
        assert set(totals) == set(range(num_factors))
        for total in totals.values():
            assert total == pytest.approx(1.0)

    def test_a_single_factor_is_exact_at_any_order(self):
        assert _suzuki_schedule(8, 1, 5) == [(0, 1.0)]

    def test_factor_count_grows_by_five_per_two_orders(self):
        sizes = [len(_suzuki_schedule(order, 3, 1, merge_adjacent=False)) for order in (2, 4, 6)]
        # order 2 has 2n - 1 = 5 entries; each +2 orders composes five sweeps of the one below
        assert sizes == [5, 25, 125]

    def test_merging_removes_the_repeated_factor_at_the_sub_sweep_seams(self):
        assert [len(_suzuki_schedule(order, 3, 1)) for order in (2, 4)] == [5, 21]


class TestSynthesize:
    def test_one_factor_per_group_with_the_whole_time(self):
        operator = bose_hubbard()
        factors = _factors(BosonicLieTrotter(), operator)
        assert len(factors) == operator.num_groups()
        assert all(gate.params[0] == pytest.approx(TIME) for gate, _ in factors)

    def test_emitted_factors_are_atomic_and_carry_their_group(self):
        operator = bose_hubbard()
        groups = operator.split_out_groups()
        factors = _factors(BosonicLieTrotter(), operator)
        for (gate, _), group in zip(factors, groups, strict=True):
            assert isinstance(gate, Evolution)
            assert gate.atomic
            # relabeling onto the active modes preserves the value of the factor's *support*
            assert gate.operator.get_support() == set(range(len(gate.operator.get_support())))
            assert len(gate.operator) == len(group)

    def test_each_factor_acts_only_on_the_modes_it_touches(self):
        operator = bose_hubbard()
        for (gate, qargs), group in zip(
            _factors(BosonicLieTrotter(), operator), operator.split_out_groups(), strict=True
        ):
            assert gate.num_modes == len(group.get_support()) == len(qargs)

    def test_ungrouped_operators_are_split_term_by_term(self):
        operator = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        assert len(_factors(BosonicLieTrotter(), operator)) == len(operator)

    def test_repeated_modes_are_not_dropped(self):
        """There is no exclusion principle: a term repeating a mode is a perfectly good factor."""
        operator = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        ((gate, qargs),) = _factors(BosonicLieTrotter(), operator)
        assert len(qargs) == 1
        assert max(len(term) for term, _ in gate.operator.iter_terms()) == 4

    def test_a_single_group_is_one_gate_at_any_order(self):
        operator = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        assert len(_factors(BosonicSuzukiTrotter(order=6, reps=4), operator)) == 1

    def test_second_order_emits_the_palindrome(self):
        operator = bose_hubbard()
        factors = _factors(BosonicSuzukiTrotter(order=2), operator)
        assert len(factors) == 2 * operator.num_groups() - 1
        times = [gate.params[0] for gate, _ in factors]
        assert times == pytest.approx(times[::-1])

    def test_total_time_per_factor_is_preserved(self):
        operator = bose_hubbard()
        groups = operator.split_out_groups()
        by_support: dict = defaultdict(float)
        for gate, qargs in _factors(BosonicSuzukiTrotter(order=4, reps=2), operator):
            by_support[tuple(q._index for q in qargs)] += gate.params[0]
        assert len(by_support) == len(groups)
        for total in by_support.values():
            assert total == pytest.approx(TIME)


class TestAccuracy:
    """The outer Trotter error between non-commuting groups, against the exact evolution."""

    @staticmethod
    def _error(**synthesis_kwargs) -> float:
        operator = bose_hubbard()
        circuit = evolution_circuit(
            operator, TIME, synthesis=BosonicSuzukiTrotter(**synthesis_kwargs)
        )
        return distance(circuit_state(circuit.decompose()), exact_state(operator, TIME))

    def test_first_order_has_a_genuine_error(self):
        # guards the methodology: if this were ~1e-16 the comparisons below would be vacuous
        assert self._error(order=1) > 1e-3

    def test_more_repetitions_reduce_the_error(self):
        errors = [self._error(order=1, reps=reps) for reps in (1, 2, 4, 8)]
        assert errors == sorted(errors, reverse=True)
        # first order converges linearly in the step size
        assert errors[0] / errors[-1] == pytest.approx(8, rel=0.5)

    def test_second_order_beats_first_order(self):
        assert self._error(order=2) < self._error(order=1) / 3

    def test_higher_order_beats_lower_order(self):
        assert self._error(order=4) < self._error(order=2)

    def test_second_order_converges_quadratically(self):
        coarse, fine = self._error(order=2, reps=2), self._error(order=2, reps=4)
        assert coarse / fine == pytest.approx(4, rel=0.5)

    def test_matches_the_lie_trotter_class_at_equal_reps(self):
        operator = bose_hubbard()
        lie = evolution_circuit(operator, TIME, synthesis=BosonicLieTrotter(reps=3))
        suzuki = evolution_circuit(operator, TIME, synthesis=BosonicSuzukiTrotter(order=1, reps=3))
        assert distance(
            circuit_state(lie.decompose()), circuit_state(suzuki.decompose())
        ) == pytest.approx(0.0, abs=1e-12)

    def test_the_default_synthesis_is_first_order_lie_trotter(self):
        operator = bose_hubbard()
        default = evolution_circuit(operator, TIME)
        explicit = evolution_circuit(operator, TIME, synthesis=BosonicLieTrotter())
        assert distance(
            circuit_state(default.decompose()), circuit_state(explicit.decompose())
        ) == pytest.approx(0.0, abs=1e-12)

    def test_undecomposed_evolution_is_exact_at_this_level(self):
        """Handing the whole gate to the boson-to-qubit stage incurs no bosonic Trotter error."""
        operator = bose_hubbard()
        whole = evolution_circuit(operator, TIME)
        assert distance(circuit_state(whole), exact_state(operator, TIME)) < 1e-6
