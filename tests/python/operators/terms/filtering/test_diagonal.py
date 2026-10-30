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

"""Tests for filtering out the diagonal terms of an operator."""

from __future__ import annotations

import copy

from qiskit_bosons.operators import BosonOperator
from qiskit_bosons.operators.terms.filtering import filter_diagonal_terms


def _is_diagonal(term):
    creations = sorted(idx for action, idx in term if action)
    annihilations = sorted(idx for action, idx in term if not action)
    return creations == annihilations


def test_filter_drops_constant_number_and_products():
    op = BosonOperator.from_dict(
        {
            (): 1.0,  # constant (diagonal)
            ((True, 0), (False, 0)): 2.0,  # n_0 (diagonal)
            ((True, 0), (True, 1), (False, 1), (False, 0)): 3.0,  # n_0 n_1 (diagonal)
            ((True, 0), (False, 1)): 4.0,  # b†_0 b_1 (off-diagonal)
        }
    )

    res = filter_diagonal_terms(op)
    assert res is None, "The function mutates in place and returns None."

    terms = [term for term, _ in op.iter_terms()]
    assert terms == [[(True, 0), (False, 1)]], "Only the off-diagonal term should survive."


def test_filter_drops_repeated_mode_number_products():
    # Unlike for fermions, a mode may repeat: b†_0 b†_0 b_0 b_0 = n_0 (n_0 - 1) is diagonal too.
    op = BosonOperator.from_dict(
        {
            ((True, 0), (True, 0), (False, 0), (False, 0)): 1.0,  # diagonal
            ((True, 0), (True, 0), (False, 0), (False, 1)): 2.0,  # off-diagonal
            ((True, 0), (True, 0)): 3.0,  # off-diagonal: changes the occupation by two
        }
    )

    filter_diagonal_terms(op)

    terms = sorted(tuple(term) for term, _ in op.iter_terms())
    assert terms == sorted(
        [
            ((True, 0), (True, 0), (False, 0), (False, 1)),
            ((True, 0), (True, 0)),
        ]
    )


def test_filter_keeps_off_diagonal_terms():
    data = {((True, 0), (False, 1)): 1.0, ((True, 1), (False, 0)): 2.0}
    op = BosonOperator.from_dict(data)
    filter_diagonal_terms(op)
    assert op.equiv(BosonOperator.from_dict(data))


def test_filter_does_not_mutate_a_deepcopy_source():
    op = BosonOperator.from_dict(
        {(): 1.0, ((True, 0), (False, 1)): 2.0},
    )
    other = copy.deepcopy(op)
    filter_diagonal_terms(other)
    # the original retains both terms (incl. the diagonal constant)
    assert len(list(op.iter_terms())) == 2
    assert len(list(other.iter_terms())) == 1


def test_filter_renumbers_groups_contiguously():
    op = BosonOperator.from_terms_with_groups(
        [
            (((True, 0), (False, 0)), 1.0, 0),  # diagonal, whole group 0 dropped
            (((True, 0), (False, 1)), 1.0, 1),
            (((True, 1), (False, 0)), 1.0, 2),
        ]
    )

    filter_diagonal_terms(op)

    assert sorted(group for _, _, group in op.iter_terms_with_groups()) == [0, 1]


def test_filter_bose_hubbard_hamiltonian():
    # A normal-ordered Bose-Hubbard chain: hopping is off-diagonal, on-site terms are diagonal.
    num_sites = 3
    terms = {}
    for site in range(num_sites - 1):
        terms[((True, site), (False, site + 1))] = -1.0
        terms[((True, site + 1), (False, site))] = -1.0
    for site in range(num_sites):
        terms[((True, site), (False, site))] = 0.5
        terms[((True, site), (True, site), (False, site), (False, site))] = 2.0
    op = BosonOperator.from_dict(terms)
    num_terms_before = len(list(op.iter_terms()))

    filter_diagonal_terms(op)

    surviving = [term for term, _ in op.iter_terms()]
    # every surviving term is off-diagonal, and exactly the hopping terms survive
    assert all(not _is_diagonal(term) for term in surviving)
    assert len(surviving) == 2 * (num_sites - 1)
    assert len(surviving) < num_terms_before
