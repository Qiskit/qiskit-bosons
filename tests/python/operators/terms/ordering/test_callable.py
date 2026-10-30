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

"""Tests for ordering the terms of an operator by a user-provided key."""

from __future__ import annotations

from qiskit_bosons.operators import BosonOperator
from qiskit_bosons.operators.terms.ordering import order_terms


def test_order_terms_sorts_by_custom_key():
    # Two terms whose coefficient magnitudes fix the requested order.
    op = BosonOperator.from_dict(
        {
            ((True, 1), (False, 0)): 1.0,  # |coeff| = 1
            ((True, 0), (False, 1)): 3.0,  # |coeff| = 3
        }
    )

    ordered = order_terms(op, key=lambda term: abs(term[1]))

    # A new operator is returned; the input is untouched.
    assert ordered is not op
    assert op.equiv(ordered)
    # Ascending |coeff|: the magnitude-1 term comes first, carrying its own coefficient.
    assert list(ordered.iter_terms()) == [
        ([(True, 1), (False, 0)], 1 + 0j),
        ([(True, 0), (False, 1)], 3 + 0j),
    ]


def test_order_terms_reverse():
    op = BosonOperator.from_dict(
        {
            ((True, 1), (False, 0)): 1.0,
            ((True, 0), (False, 1)): 3.0,
        }
    )

    ordered = order_terms(op, key=lambda term: abs(term[1]), reverse=True)

    # Descending |coeff|: the magnitude-3 term comes first.
    assert op.equiv(ordered)
    assert list(ordered.iter_terms()) == [
        ([(True, 0), (False, 1)], 3 + 0j),
        ([(True, 1), (False, 0)], 1 + 0j),
    ]


def test_order_terms_preserves_groups():
    op = BosonOperator.from_dict(
        {
            ((True, 1), (False, 0)): 1.0,
            ((True, 0), (False, 1)): 3.0,
        }
    )
    # Give the two terms distinct group tags; the term<->tag association must survive reordering.
    op.groups = [0, 1]
    before = {tuple(term): group for term, _, group in op.iter_terms_with_groups()}

    ordered = order_terms(op, key=lambda term: abs(term[1]), reverse=True)

    # Groups are preserved, and each term still carries the same group tag it had before.
    assert ordered.groups is not None
    after = {tuple(term): group for term, _, group in ordered.iter_terms_with_groups()}
    assert after == before


def test_order_terms_key_on_group():
    op = BosonOperator.from_dict(
        {
            ((True, 0), (False, 1)): 3.0,
            ((True, 1), (False, 0)): 1.0,
        }
    )
    # `from_dict` does not preserve dict order, so tag by the actual stored positions: term 0 gets
    # group 1 and term 1 gets group 0, i.e. the stored order starts out group-descending.
    op.groups = [1, 0]
    before = {tuple(term): group for term, _, group in op.iter_terms_with_groups()}

    ordered = order_terms(op, key=lambda term: term[2])

    # Sorted by ascending group index; each term keeps the tag it started with.
    assert op.equiv(ordered)
    assert [group for _, _, group in ordered.iter_terms_with_groups()] == [0, 1]
    after = {tuple(term): group for term, _, group in ordered.iter_terms_with_groups()}
    assert after == before


def test_order_terms_keeps_repeated_modes():
    # b†_0 b†_0 is a legitimate bosonic term; reordering must neither drop nor merge it.
    op = BosonOperator.from_dict({((True, 1),): 1.0, ((True, 0), (True, 0)): 2.0})

    ordered = order_terms(op, key=lambda term: len(term[0]), reverse=True)

    assert isinstance(ordered, BosonOperator)
    assert op.equiv(ordered)
    assert list(ordered.iter_terms()) == [
        ([(True, 0), (True, 0)], 2 + 0j),
        ([(True, 1)], 1 + 0j),
    ]
