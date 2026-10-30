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

"""Commutator and quadrature library tests.

For bosons the *commutator* is the bracket that defines the algebra, which is the reverse of the
fermionic situation. These tests cover the algebraic identities the brackets must satisfy and the
normalization of the quadrature operators, which is easy to get wrong.
"""

from __future__ import annotations

import unittest

import numpy as np
import pytest
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.operators.library import (
    anti_commutator,
    commutator,
    double_commutator,
    momentum,
    position,
)


def ladder(mode: int, *, raising: bool) -> BosonOperator:
    action = cre(mode) if raising else ann(mode)
    return BosonOperator.from_dict({(action,): 1.0})


def number(mode: int) -> BosonOperator:
    return BosonOperator.from_dict({(cre(mode), ann(mode)): 1.0})


class TestCommutator(unittest.TestCase):
    def test_anti_symmetry(self):
        first = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
        second = BosonOperator.from_dict({(cre(1), ann(0)): 2.0})
        assert commutator(first, second).equiv(-commutator(second, first))

    def test_self_commutator_vanishes(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(1)): 1.0 - 2.0j})
        assert commutator(op, op).equiv(BosonOperator.zero())

    def test_linearity(self):
        first = ladder(0, raising=True)
        second = ladder(0, raising=False)
        third = number(1)
        assert commutator(first + third, second).equiv(
            commutator(first, second) + commutator(third, second)
        )
        assert commutator(3.0 * first, second).equiv(3.0 * commutator(first, second))

    def test_jacobi_identity(self):
        a = ladder(0, raising=True)
        b = ladder(0, raising=False)
        c = number(1) + ladder(1, raising=True)
        total = (
            commutator(a, commutator(b, c))
            + commutator(b, commutator(c, a))
            + commutator(c, commutator(a, b))
        )
        assert total.equiv(BosonOperator.zero())

    def test_identity_commutes_with_everything(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(1)): 1.0})
        assert commutator(op, BosonOperator.one()).equiv(BosonOperator.zero())

    def test_number_operators_commute(self):
        assert commutator(number(0), number(1)).equiv(BosonOperator.zero())
        assert commutator(number(0), number(0)).equiv(BosonOperator.zero())

    def test_number_counts_repeated_creations(self):
        """``[n_0, b^dag_0 b^dag_0] = 2 b^dag_0 b^dag_0`` -- the eigenvalue is the boson count."""
        pair = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        assert commutator(number(0), pair).equiv(2.0 * pair)

    def test_hopping_terms_on_disjoint_bonds_commute(self):
        first = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        second = BosonOperator.from_dict({(cre(2), ann(3)): 1.0, (cre(3), ann(2)): 1.0})
        assert commutator(first, second).equiv(BosonOperator.zero())

    def test_hopping_terms_sharing_a_site_do_not_commute(self):
        first = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        second = BosonOperator.from_dict({(cre(1), ann(2)): 1.0, (cre(2), ann(1)): 1.0})
        assert not commutator(first, second).equiv(BosonOperator.zero())

    def test_total_number_commutes_with_a_hopping_term(self):
        """A passive (number-conserving) term commutes with the total boson number."""
        total = number(0) + number(1)
        hopping = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        assert commutator(total, hopping).equiv(BosonOperator.zero())

    def test_total_number_does_not_commute_with_a_bare_creation(self):
        total = number(0) + number(1)
        assert not commutator(total, ladder(0, raising=True)).equiv(BosonOperator.zero())


class TestAntiCommutator(unittest.TestCase):
    def test_symmetry(self):
        first = ladder(0, raising=True)
        second = ladder(0, raising=False)
        assert anti_commutator(first, second).equiv(anti_commutator(second, first))

    def test_ladder_pair(self):
        expected = BosonOperator.from_dict({(cre(0), ann(0)): 2.0, (): 1.0})
        assert anti_commutator(ladder(0, raising=False), ladder(0, raising=True)).equiv(expected)

    def test_self_anti_commutator_is_twice_the_square(self):
        op = BosonOperator.from_dict({(cre(0),): 1.0})
        assert anti_commutator(op, op).equiv(2.0 * (op @ op))

    def test_bosonic_ladder_anti_commutator_is_not_a_delta(self):
        """The fermionic identity ``{a, a^dag} = 1`` has no bosonic analogue."""
        got = anti_commutator(ladder(0, raising=False), ladder(0, raising=True))
        assert not got.equiv(BosonOperator.one())

    def test_repeated_creation_does_not_vanish(self):
        """A fermion would give ``{a^dag, a^dag} = 2 a^dag a^dag = 0``; a boson does not."""
        raising = ladder(0, raising=True)
        assert not anti_commutator(raising, raising).equiv(BosonOperator.zero())


class TestDoubleCommutator(unittest.TestCase):
    @staticmethod
    def _operands():
        return (
            BosonOperator.from_dict({(cre(0), ann(1)): 1.0}),
            BosonOperator.from_dict({(cre(1), ann(0)): 1.0}),
            number(0),
        )

    def test_sign_variants_differ(self):
        first, second, third = self._operands()
        assert not double_commutator(first, second, third, False).equiv(
            double_commutator(first, second, third, True)
        )

    def test_linear_in_each_argument(self):
        first, second, third = self._operands()
        for sign in (False, True):
            base = double_commutator(first, second, third, sign)
            assert double_commutator(3.0 * first, second, third, sign).equiv(3.0 * base)
            assert double_commutator(first, 3.0 * second, third, sign).equiv(3.0 * base)
            assert double_commutator(first, second, 3.0 * third, sign).equiv(3.0 * base)

    def test_matches_the_explicit_definition(self):
        """``sign=False`` is ``[[A,B],C]/2 + [A,[B,C]]/2``."""
        first, second, third = self._operands()
        expected = (
            commutator(commutator(first, second), third)
            + commutator(first, commutator(second, third))
        ) / 2
        assert double_commutator(first, second, third, False).equiv(expected)

    def test_vanishes_when_everything_commutes(self):
        for sign in (False, True):
            got = double_commutator(number(0), number(1), number(2), sign)
            assert got.equiv(BosonOperator.zero())


class TestQuadratures(unittest.TestCase):
    def test_both_are_hermitian(self):
        assert position(0).is_hermitian()
        assert momentum(0).is_hermitian()

    def test_position_normalization(self):
        """``q = (b^dag + b)/sqrt(2)``."""
        scale = 2.0**-0.5
        expected = BosonOperator.from_dict({(cre(1),): scale, (ann(1),): scale})
        assert position(1).equiv(expected)

    def test_momentum_normalization(self):
        """``p = i(b^dag - b)/sqrt(2)``."""
        scale = 2.0**-0.5
        expected = BosonOperator.from_dict({(cre(2),): 1j * scale, (ann(2),): -1j * scale})
        assert momentum(2).equiv(expected)

    def test_canonical_commutator(self):
        """``[q, p] = i``, which is the whole point of the ``1/sqrt(2)`` normalization."""
        got = commutator(position(0), momentum(0))
        assert got.equiv(1j * BosonOperator.one())

    def test_quadratures_on_distinct_modes_commute(self):
        assert commutator(position(0), momentum(1)).equiv(BosonOperator.zero())
        assert commutator(position(0), position(1)).equiv(BosonOperator.zero())

    def test_harmonic_oscillator_hamiltonian(self):
        """``(q^2 + p^2)/2 = n + 1/2``, the standard oscillator spectrum offset."""
        squares = (position(0) @ position(0) + momentum(0) @ momentum(0)) / 2
        expected = number(0) + 0.5 * BosonOperator.one()
        assert squares.equiv(expected)

    def test_quadrature_squares_demand_two_levels(self):
        assert (position(0) @ position(0)).max_occupation() == 2

    def test_position_and_momentum_are_not_equal(self):
        assert not position(0).equiv(momentum(0))

    def test_mode_index_is_honored(self):
        for mode in (0, 3, 7):
            assert position(mode).get_support() == {mode}
            assert momentum(mode).get_support() == {mode}


class TestTypeContract(unittest.TestCase):
    def test_output_type_follows_the_input(self):
        first = ladder(0, raising=True)
        second = ladder(0, raising=False)
        assert isinstance(commutator(first, second), BosonOperator)
        assert isinstance(anti_commutator(first, second), BosonOperator)
        assert isinstance(double_commutator(first, second, first, False), BosonOperator)

    def test_commutator_does_not_mutate_its_inputs(self):
        first = ladder(0, raising=True)
        second = ladder(0, raising=False)
        before = (len(first), len(second), first.get_coeffs(), second.get_coeffs())
        commutator(first, second)
        assert (len(first), len(second), first.get_coeffs(), second.get_coeffs()) == before

    def test_mismatched_type_is_rejected(self):
        with pytest.raises((TypeError, AttributeError)):
            commutator(ladder(0, raising=True), np.eye(2))  # type: ignore[arg-type]
