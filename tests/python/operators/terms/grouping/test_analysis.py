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

"""Tests for the analysis of an existing grouping."""

from __future__ import annotations

import pytest
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.operators.terms.grouping import (
    group_coeff_means,
    groups_are_hermitian,
    groups_have_uniform_coeffs,
)


def _op(coeffs, groups=None):
    """``c0 b†_0 b_1 + c1 b†_1 b_0 + c2 b†_2 b_3``, in this exact term order.

    The first two terms are each other's adjoint, so they are a Hermitian pair when they share a
    group and carry equal coefficients. The third is unpaired and is never Hermitian on its own.
    """
    op = BosonOperator(
        coeffs, [True, False, True, False, True, False], [0, 1, 1, 0, 2, 3], [0, 2, 4, 6]
    )
    if groups is not None:
        op.groups = groups
    return op


class TestGroupCoeffMeans:
    def test_mean_not_sum_and_magnitude_not_real_part(self):
        op = _op([1.0, -3.0 + 4.0j, 2.0], [0, 0, 1])
        # group 0 averages |1| and |-3+4j| = 5, i.e. 3; group 1 holds the single |2|
        assert group_coeff_means(op) == pytest.approx([3.0, 2.0])

    def test_ungrouped_is_none(self):
        assert group_coeff_means(_op([1.0, 1.0, 1.0])) is None

    def test_unused_index_weighs_zero(self):
        # a weight of 0.0 (rather than NaN) keeps the index out of a sample
        op = _op([1.0, 1.0, 2.0], [0, 0, 2])
        assert group_coeff_means(op) == pytest.approx([1.0, 0.0, 2.0])

    def test_grouped_operator_without_terms(self):
        op = BosonOperator.zero()
        op.groups = []
        assert group_coeff_means(op) == []

    def test_rejects_unsupported_type(self):
        with pytest.raises(TypeError, match="expects a BosonOperator"):
            group_coeff_means("not an operator")


class TestGroupsAreHermitian:
    def test_pair_and_unpaired_term(self):
        op = _op([1.0, 1.0, 1.0], [0, 0, 1])
        assert groups_are_hermitian(op) == [True, False]

    def test_ungrouped_is_none(self):
        assert groups_are_hermitian(_op([1.0, 1.0, 1.0])) is None

    def test_unused_index_is_hermitian(self):
        # the zero operator is Hermitian
        op = _op([1.0, 1.0, 1.0], [0, 0, 2])
        assert groups_are_hermitian(op) == [True, True, False]

    def test_mismatched_pair_coefficients_are_not_hermitian(self):
        # b†_0 b_1 + 2 b†_1 b_0 is not its own adjoint
        op = _op([1.0, 2.0, 1.0], [0, 0, 1])
        assert groups_are_hermitian(op) == [False, False]

    def test_atol_is_keyword_only_and_respected(self):
        op = _op([1.0, 1.0 + 1e-6, 1.0], [0, 0, 1])
        assert groups_are_hermitian(op) == [False, False]
        assert groups_are_hermitian(op, atol=1e-4) == [True, False]
        with pytest.raises(TypeError):
            groups_are_hermitian(op, 1e-4)  # type: ignore[misc]

    def test_repeated_modes_are_not_collapsed(self):
        """There is no exclusion principle: a mode repeats within a term and is still analyzed."""
        op = BosonOperator.from_dict(
            {
                (cre(0), ann(0)): 1.0,  # the number operator, Hermitian by itself
                (cre(0), cre(0), ann(0)): 1.0,  # b†b†b alone is not
            }
        )
        # assign by the actual stored term order rather than assuming one
        op.groups = [0 if len(term) == 2 else 1 for term, _ in op.iter_terms()]
        assert groups_are_hermitian(op) == [True, False]

    def test_group_formed_by_a_term_and_its_adjoint(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(1)): 1.0})
        op = op + op.adjoint()
        op.groups = [0] * len(op.get_coeffs())
        assert groups_are_hermitian(op) == [True]

    def test_rejects_unsupported_type(self):
        with pytest.raises(TypeError, match="expects a BosonOperator"):
            groups_are_hermitian("not an operator")


class TestGroupsHaveUniformCoeffs:
    def test_uniform_and_not(self):
        op = _op([1.0, 1.0, 2.0], [0, 0, 1])
        assert groups_have_uniform_coeffs(op) == [True, True]
        assert groups_have_uniform_coeffs(_op([1.0, 2.0, 2.0], [0, 0, 1])) == [False, True]

    def test_ungrouped_is_none(self):
        assert groups_have_uniform_coeffs(_op([1.0, 1.0, 1.0])) is None

    def test_abs_compares_magnitudes(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0j, (cre(1), ann(0)): -1.0j})
        op.groups = [0, 0]
        assert groups_have_uniform_coeffs(op) == [True]  # equal magnitudes
        assert groups_have_uniform_coeffs(op, abs=False) == [False]  # but not equal coefficients

    def test_atol(self):
        op = _op([1.0, 1.0 + 1e-10, 1.0], [0, 0, 1])
        assert groups_have_uniform_coeffs(op) == [True, True]
        assert groups_have_uniform_coeffs(op, atol=1e-12) == [False, True]

    def test_unused_and_single_term_groups_are_trivially_uniform(self):
        op = _op([1.0, 2.0, 3.0], [0, 1, 3])
        assert groups_have_uniform_coeffs(op) == [True, True, True, True]

    def test_rejects_unsupported_type(self):
        with pytest.raises(TypeError, match="expects a BosonOperator"):
            groups_have_uniform_coeffs("not an operator")


def test_hermiticity_and_uniformity_are_independent():
    """Neither check implies the other, in either direction."""
    # uniform coefficients, but not Hermitian: the terms are not each other's adjoint
    uniform_not_hermitian = BosonOperator(
        [1.0, 1.0], [True, False, True, False], [0, 1, 2, 3], [0, 2, 4]
    )
    uniform_not_hermitian.groups = [0, 0]
    assert groups_are_hermitian(uniform_not_hermitian) == [False]
    assert groups_have_uniform_coeffs(uniform_not_hermitian) == [True]

    # Hermitian, but with mixed magnitudes: two conjugate pairs at different scales
    hermitian_mixed = BosonOperator(
        [1.0, 1.0, 5.0, 5.0],
        [True, False, True, False, True, False, True, False],
        [0, 1, 1, 0, 2, 3, 3, 2],
        [0, 2, 4, 6, 8],
    )
    hermitian_mixed.groups = [0, 0, 0, 0]
    assert groups_are_hermitian(hermitian_mixed) == [True]
    assert groups_have_uniform_coeffs(hermitian_mixed) == [False]
