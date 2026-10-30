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

"""BosonOperator tests.

These tests deliberately concentrate on what is *bosonic* about this operator, i.e. on the two
consequences of the canonical commutation relations that separate it from its fermionic
counterpart: there is no exclusion principle, and a contraction during reordering enters with a
plus sign.
"""

from __future__ import annotations

import copy
import pickle
import unittest

import numpy as np
import pytest
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.operators.library import anti_commutator, commutator
from qiskit_bosons.operators.terms.grouping import group_coeff_means
from qiskit_bosons.operators.terms.ordering import group_order


def _num(mode: int) -> BosonOperator:
    """The number operator ``n_mode = b^dagger_mode b_mode``."""
    return BosonOperator.from_dict({(cre(mode), ann(mode)): 1.0})


class TestCanonicalCommutationRelations(unittest.TestCase):
    """The defining algebra: ``[b_i, b_j^dagger] = delta_ij``."""

    def test_commutator_is_kronecker_delta(self):
        num_modes = 3
        for i in range(num_modes):
            for j in range(num_modes):
                b_i = BosonOperator.from_dict({(ann(i),): 1.0})
                b_dag_j = BosonOperator.from_dict({(cre(j),): 1.0})
                expected = BosonOperator.one() if i == j else BosonOperator.zero()
                assert commutator(b_i, b_dag_j).equiv(expected), f"[b_{i}, b^dag_{j}]"

    def test_like_operators_commute(self):
        """``[b_i, b_j] = [b_i^dagger, b_j^dagger] = 0``, including for ``i == j``."""
        for i, j in [(0, 0), (0, 1), (1, 0), (2, 1)]:
            anns = (
                BosonOperator.from_dict({(ann(i),): 1.0}),
                BosonOperator.from_dict({(ann(j),): 1.0}),
            )
            cres = (
                BosonOperator.from_dict({(cre(i),): 1.0}),
                BosonOperator.from_dict({(cre(j),): 1.0}),
            )
            assert commutator(*anns).equiv(BosonOperator.zero())
            assert commutator(*cres).equiv(BosonOperator.zero())

    def test_anti_commutator_is_not_the_defining_bracket(self):
        """``{b, b^dagger} = 2 n + 1``, which is emphatically not the delta a fermion would give."""
        b = BosonOperator.from_dict({(ann(0),): 1.0})
        b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
        expected = BosonOperator.from_dict({(cre(0), ann(0)): 2.0, (): 1.0})
        assert anti_commutator(b, b_dag).equiv(expected)
        # in particular it is *not* the identity
        assert not anti_commutator(b, b_dag).equiv(BosonOperator.one())

    def test_number_operator_ladder(self):
        """``[n, b^dagger] = b^dagger`` and ``[n, b] = -b``."""
        b = BosonOperator.from_dict({(ann(0),): 1.0})
        b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
        assert commutator(_num(0), b_dag).equiv(b_dag)
        assert commutator(_num(0), b).equiv(-b)


class TestNoExclusionPrinciple(unittest.TestCase):
    """The single most important bosonic invariant: ``b^dagger_0 b^dagger_0 != 0``."""

    def test_repeated_creation_is_not_zero(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        assert len(op) == 1
        assert not op.equiv(BosonOperator.zero())

    def test_repeated_creation_survives_normal_ordering(self):
        """Normal ordering must not silently drop a repeated-index term."""
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        ordered = op.normal_ordered().simplify()
        assert not ordered.equiv(BosonOperator.zero())
        assert ordered.equiv(op)
        # the term is stored verbatim, i.e. the mode index really does repeat
        assert sorted(ordered.get_modes()) == [0, 0]
        assert ordered.get_actions() == [True, True]

    def test_repeated_creation_survives_simplify(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        assert op.simplify().equiv(op)

    def test_power_equals_repeated_term(self):
        """``(b^dagger)**2 == b^dagger b^dagger``, not zero."""
        b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
        squared = (b_dag**2).simplify()
        assert squared.equiv(BosonOperator.from_dict({(cre(0), cre(0)): 1.0}))
        assert not squared.equiv(BosonOperator.zero())

    def test_repeated_annihilation_is_not_zero(self):
        b = BosonOperator.from_dict({(ann(0),): 1.0})
        assert not (b**2).simplify().equiv(BosonOperator.zero())

    def test_higher_powers_are_not_zero(self):
        """Arbitrarily many bosons fit in one mode, so no power ever vanishes."""
        b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
        for power in range(1, 6):
            op = (b_dag**power).simplify()
            assert not op.equiv(BosonOperator.zero()), f"(b^dag)**{power}"
            assert op.max_rank() == power


class TestNormalOrdered(unittest.TestCase):
    """The contraction sign, which is the other half of what makes this class bosonic."""

    def test_contraction_adds_with_a_plus(self):
        """``b b^dagger -> b^dagger b + 1``."""
        op = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        ordered = op.normal_ordered().simplify()
        expected = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (): 1.0})
        assert ordered.equiv(expected)

    def test_contraction_is_not_the_fermionic_answer(self):
        """Explicitly *not* ``-b^dagger b + 1``."""
        op = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        ordered = op.normal_ordered().simplify()
        fermionic = BosonOperator.from_dict({(cre(0), ann(0)): -1.0, (): 1.0})
        assert not ordered.equiv(fermionic)

    def test_no_contraction_across_distinct_modes(self):
        """``b_0 b^dagger_1 -> b^dagger_1 b_0`` with no constant and no sign."""
        op = BosonOperator.from_dict({(ann(0), cre(1)): 1.0})
        ordered = op.normal_ordered().simplify()
        expected = BosonOperator.from_dict({(cre(1), ann(0)): 1.0})
        assert ordered.equiv(expected)

    def test_double_contraction_constant(self):
        """``b b b^dagger b^dagger = b^dagger b^dagger b b + 4 b^dagger b + 2``.

        The constant is what a fermionic reordering could never produce; it is pinned below by
        acting on the vacuum, where every ladder term dies and only the ``2`` survives.
        """
        op = BosonOperator.from_dict({(ann(0), ann(0), cre(0), cre(0)): 1.0})
        ordered = op.normal_ordered().simplify()
        expected = BosonOperator.from_dict(
            {
                (cre(0), cre(0), ann(0), ann(0)): 1.0,
                (cre(0), ann(0)): 4.0,
                (): 2.0,
            }
        )
        assert ordered.equiv(expected)

    def test_double_contraction_constant_on_vacuum(self):
        """Pins the ``+2`` above: ``<vac| b b b^dag b^dag |vac> = 2``.

        ``b b^dagger b b^dagger`` applied to the vacuum gives ``b b^dagger |0> = |0>`` twice over,
        but the operator here is ``b^2 (b^dagger)^2``, whose vacuum expectation value is
        ``|| (b^dagger)^2 |0> ||^2 = || sqrt(2) |2> ||^2 = 2``.
        """
        ordered = (
            BosonOperator.from_dict({(ann(0), ann(0), cre(0), cre(0)): 1.0})
            .normal_ordered()
            .simplify()
        )
        # every non-constant normal-ordered term ends in an annihilation, hence kills the vacuum
        constants = [coeff for term, coeff in ordered.iter_terms() if not term]
        assert len(constants) == 1
        assert constants[0] == pytest.approx(2.0)

    def test_sandwich_variants_are_equivalent(self):
        """The three orderings differ structurally but describe the same operator."""
        op = BosonOperator.from_dict({(ann(1), cre(1), ann(0), cre(0)): 1.0})
        default = op.normal_ordered().simplify()
        inner = op.normal_ordered(sandwich=True).simplify()
        outer = op.normal_ordered(sandwich=False).simplify()
        assert default.equiv(inner)
        assert default.equiv(outer)
        # ... while being stored differently
        assert default != inner

    def test_normal_ordered_is_idempotent(self):
        op = BosonOperator.from_dict({(ann(0), cre(0), ann(1), cre(1)): 1.0})
        once = op.normal_ordered().simplify()
        twice = once.normal_ordered().simplify()
        assert once.equiv(twice)

    def test_normal_ordering_grows_the_term_count(self):
        """Reordering a single term can produce several: the contraction is not free."""
        op = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        assert len(op) == 1
        assert len(op.normal_ordered().simplify()) == 2


class TestEquivVersusSimplify(unittest.TestCase):
    """``equiv`` normal-orders; ``simplify`` is purely structural. This locks the distinction."""

    def test_equiv_identifies_bb_dagger_with_number_plus_one(self):
        lhs = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        rhs = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (): 1.0})
        assert lhs.equiv(rhs)

    def test_simplify_does_not_identify_them(self):
        """``simplify`` only merges duplicate *stored* terms, so it cannot see the identity."""
        lhs = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        rhs = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (): 1.0})
        assert lhs.simplify() != rhs.simplify()
        assert len(lhs.simplify()) == 1
        assert len(rhs.simplify()) == 2

    def test_equality_operator_is_structural(self):
        lhs = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        rhs = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (): 1.0})
        assert lhs != rhs
        assert lhs.equiv(rhs)

    def test_equiv_tolerance(self):
        op = BosonOperator.from_dict({(): 1e-7})
        zero = BosonOperator.zero()
        assert not op.equiv(zero)
        assert op.equiv(zero, 1e-6)
        assert not op.equiv(zero, 1e-9)

    def test_simplify_merges_duplicates_before_chopping(self):
        """A pile of tiny duplicates sums to something that must not be chopped away."""
        count = 1000
        coeffs = [1e-5] * count
        op = BosonOperator(coeffs, [], [], [0] * (count + 1))
        assert op.simplify(1e-4).equiv(BosonOperator.one() * (count * 1e-5), 1e-9)
        # whereas ichop truncates each term greedily
        greedy = copy.deepcopy(op)
        greedy.ichop(1e-4)
        assert greedy.equiv(BosonOperator.zero())


class TestMaxOccupation(unittest.TestCase):
    """The quantity with no fermionic counterpart."""

    def test_annihilations_only_is_zero(self):
        op = BosonOperator.from_dict({(ann(0), ann(0), ann(1)): 1.0})
        assert op.max_occupation() == 0

    def test_counts_per_mode_not_per_term(self):
        """A rank-4 term spread over four modes still only drives each to occupation 1."""
        op = BosonOperator.from_dict({(cre(0), cre(1), cre(2), cre(3)): 1.0})
        assert op.max_rank() == 4
        assert op.max_occupation() == 1

    def test_repeated_mode_raises_the_bound(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), cre(0)): 1.0})
        assert op.max_occupation() == 3

    def test_max_across_terms(self):
        op = BosonOperator.from_dict(
            {
                (cre(0), ann(1)): 1.0,
                (cre(2), cre(2), cre(2), ann(2)): 1.0,
                (cre(1), cre(1)): 1.0,
            }
        )
        assert op.max_occupation() == 3

    def test_annihilations_do_not_offset_creations(self):
        """It counts creations, not the net change: ``b^dag b^dag b b`` needs ``d >= 3``."""
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        assert op.max_occupation() == 2

    def test_empty_and_identity(self):
        assert BosonOperator.zero().max_occupation() == 0
        assert BosonOperator.one().max_occupation() == 0


class TestConservation(unittest.TestCase):
    """Particle-number and sector conservation, with the repeated indices bosons allow."""

    def test_conserves_particle_number_with_repeated_indices(self):
        # b^dag_0 b^dag_0 b_0 b_0: two in, two out
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        assert op.conserves_particle_number()

    def test_repeated_creation_does_not_conserve(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        assert not op.conserves_particle_number()

    def test_unbalanced_repeated_indices(self):
        # three creations, one annihilation, all on the same mode
        op = BosonOperator.from_dict({(cre(0), cre(0), cre(0), ann(0)): 1.0})
        assert not op.conserves_particle_number()

    def test_conserves_sector_with_repeated_indices(self):
        # both bosons created and destroyed inside block 0
        op = BosonOperator.from_dict({(cre(0), cre(1), ann(0), ann(1)): 1.0})
        assert op.conserves_sector([2, 2])
        assert op.conserves_sector([4])

    def test_conserves_sector_repeated_index_across_blocks(self):
        # b^dag_0 b^dag_0 b_2 b_2 moves two bosons from block 1 to block 0
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(2), ann(2)): 1.0})
        assert op.conserves_particle_number()
        assert not op.conserves_sector([2, 2])

    def test_empty_block_sizes_is_particle_number(self):
        op = BosonOperator.from_dict({(cre(0), ann(2)): 1.0})
        assert op.conserves_sector([]) == op.conserves_particle_number()

    def test_mode_beyond_final_block(self):
        op = BosonOperator.from_dict({(cre(0), ann(5)): 1.0})
        assert not op.conserves_sector([2, 2])

    def test_all_terms_must_conserve(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (cre(1),): 1.0})
        assert not op.conserves_particle_number()


class TestRelabelModes(unittest.TestCase):
    """Relabeling, including the two error paths."""

    def test_relabel_repeated_indices(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(1)): 1.0})
        relabeled = op.relabel_modes([2, 0])
        assert relabeled.equiv(BosonOperator.from_dict({(cre(2), cre(2), ann(0)): 1.0}))

    def test_relabel_round_trip(self):
        op = BosonOperator.from_dict({(cre(0), ann(1), cre(2)): 1.0})
        permutation = [2, 0, 1]
        inverse = [permutation.index(i) for i in range(len(permutation))]
        assert op.relabel_modes(permutation).relabel_modes(inverse).equiv(op)

    def test_duplicate_indices_raises(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
        with pytest.raises(ValueError, match=r"[Dd]uplicate"):
            op.relabel_modes([0, 0])

    def test_index_map_too_small_raises(self):
        op = BosonOperator.from_dict({(cre(0), ann(5)): 1.0})
        with pytest.raises(ValueError):
            op.relabel_modes([1, 0])

    def test_permutation_longer_than_needed_is_fine(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
        relabeled = op.relabel_modes([3, 2, 1, 0])
        assert relabeled.equiv(BosonOperator.from_dict({(cre(3), ann(2)): 1.0}))


class TestAdjoint(unittest.TestCase):
    def test_adjoint_reverses_and_flips(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 2.0j})
        adj = op.adjoint()
        assert adj.equiv(BosonOperator.from_dict({(cre(1), ann(0)): -2.0j}))

    def test_adjoint_of_repeated_creation(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        assert op.adjoint().equiv(BosonOperator.from_dict({(ann(0), ann(0)): 1.0}))

    def test_adjoint_is_involutive(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(1)): 1.0 - 2.0j, (): 3.0j})
        assert op.adjoint().adjoint().equiv(op)

    def test_adjoint_conjugates_coefficients(self):
        op = BosonOperator.from_dict({(): 1.0 + 2.0j})
        assert op.adjoint().get_coeffs() == [1.0 - 2.0j]


class TestIsHermitian(unittest.TestCase):
    def test_number_operator_is_hermitian(self):
        assert _num(0).is_hermitian()

    def test_repeated_creation_is_not_hermitian(self):
        assert not BosonOperator.from_dict({(cre(0), cre(0)): 1.0}).is_hermitian()

    def test_interaction_term_is_hermitian(self):
        """The Bose-Hubbard interaction ``b^dag b^dag b b`` is Hermitian despite the repeats."""
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        assert op.is_hermitian()

    def test_single_hopping_term_is_not_hermitian(self):
        """The lone hopping term is not Hermitian; only the bond pair is."""
        single = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
        assert not single.is_hermitian()
        paired = single + single.adjoint()
        assert paired.is_hermitian()

    def test_hermitian_requires_the_contraction(self):
        """``b b^dagger`` is Hermitian only because ``normal_ordered`` gets the plus sign right."""
        assert BosonOperator.from_dict({(ann(0), cre(0)): 1.0}).is_hermitian()

    def test_tolerance(self):
        op = BosonOperator.from_dict(
            {(cre(0), ann(1)): 1.00001j, (cre(1), ann(0)): -1j},
        )
        assert not op.is_hermitian()
        assert op.is_hermitian(1e-4)


class TestIchop(unittest.TestCase):
    def test_removes_small_terms_in_place(self):
        op = BosonOperator.from_dict({(): 1e-4, (cre(0),): 1e-10})
        assert len(op) == 2
        assert op.ichop() is None
        assert len(op) == 1

    def test_does_not_conjugate_retained_coefficients(self):
        """``ichop`` filters; it must not touch the coefficients it keeps."""
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0 + 2.0j, (): 1e-12})
        op.ichop()
        assert op.get_coeffs() == [1.0 + 2.0j]

    def test_uses_magnitude_not_real_part(self):
        """A purely imaginary coefficient of large magnitude must survive."""
        op = BosonOperator.from_dict({(cre(0),): 5.0j})
        op.ichop(1e-2)
        assert len(op) == 1
        assert op.get_coeffs() == [5.0j]

    def test_chopping_everything_leaves_zero(self):
        op = BosonOperator.from_dict({(): 1e-12, (cre(0),): 1e-12})
        op.ichop()
        assert len(op) == 0
        assert op.equiv(BosonOperator.zero())


class TestPickle(unittest.TestCase):
    def test_round_trip(self):
        op = BosonOperator.from_dict(
            {(cre(0), cre(0), ann(1)): 1.0 - 2.0j, (): 3.0, (ann(2),): -1.0j}
        )
        restored = pickle.loads(pickle.dumps(op))
        assert restored == op
        assert restored.groups is None

    def test_round_trip_preserves_groups(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        op.groups = [0, 1]
        restored = pickle.loads(pickle.dumps(op))
        assert restored == op
        assert restored.groups == [0, 1]

    def test_getnewargs_reconstructs(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 2.0})
        assert BosonOperator(*op.__getnewargs__()) == op

    def test_getstate_setstate(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        assert op.__getstate__() is None
        op.__setstate__([3])
        assert op.groups == [3]

    def test_deepcopy(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        clone = copy.deepcopy(op)
        assert clone == op
        clone += BosonOperator.one()
        assert len(op) == 1


class TestGroups(unittest.TestCase):
    def test_groups_round_trip(self):
        op = BosonOperator.from_dict(
            {(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0, (cre(0), cre(0), ann(0), ann(0)): 0.5}
        )
        op.groups = [0, 0, 1]
        assert op.has_groups()
        assert op.num_groups() == 2
        reconstructed = BosonOperator.from_terms_with_groups(op.iter_terms_with_groups())
        assert reconstructed.equiv(op)
        assert reconstructed.groups == op.groups

    def test_no_groups_by_default(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        assert not op.has_groups()
        assert op.groups is None
        assert op.num_groups() is None
        assert group_coeff_means(op) is None
        assert op.split_out_groups() is None

    def test_iter_terms_with_groups_fails_without_groups(self):
        """Iterating with groups on an ungrouped operator is an error.

        The native implementation makes this a ``panic!`` rather than a Python exception, so what
        surfaces is a ``pyo3_runtime.PanicException`` (a ``BaseException``, *not* a ``ValueError``).
        That is the documented contract of the underlying ``iter_with_groups`` -- the fermionic core
        pins it with a ``#[should_panic]`` test -- so it is locked here as-is rather than as the
        ``ValueError`` one might expect. Guard with :meth:`has_groups` instead of catching this.
        """
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        assert not op.has_groups()
        with pytest.raises(BaseException, match="only be called when groups are present"):
            list(op.iter_terms_with_groups())

    def test_groups_can_be_cleared(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        op.groups = [0]
        assert op.has_groups()
        op.groups = None
        assert not op.has_groups()

    def test_group_coeff_means_is_the_mean_absolute_magnitude(self):
        """Note the *mean*, not the sum, and that it uses ``abs`` of a complex coefficient."""
        op = BosonOperator(
            [1.0, 3.0j, -2.0],
            [True, False, True, False, True, False],
            [0, 1, 1, 0, 0, 0],
            [0, 2, 4, 6],
        )
        op.groups = [0, 0, 1]
        # group 0 holds |1| and |3j| -> mean 2.0; group 1 holds |-2| -> 2.0
        assert group_coeff_means(op) == pytest.approx([2.0, 2.0])

    def test_empty_group_weighs_zero(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (cre(1), ann(1)): 1.0})
        op.groups = [0, 2]
        assert op.num_groups() == 3
        assert group_coeff_means(op) == pytest.approx([1.0, 0.0, 1.0])

    def test_group_weights_was_removed(self):
        """The method became the free function :func:`group_coeff_means`."""
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        assert not hasattr(op, "group_weights")

    def test_groups_length_is_validated(self):
        """Too few indices would drop trailing terms, too many would invent groups."""
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        for wrong in ([], [0], [0, 0, 0]):
            with pytest.raises(ValueError, match="one group index per term"):
                op.groups = wrong
            # a rejected assignment leaves the operator untouched
            assert not op.has_groups()

        op.groups = [0, 1]
        with pytest.raises(ValueError, match="one group index per term"):
            op.groups = [0]
        assert op.groups == [0, 1]

    def test_groups_can_be_cleared_even_when_set(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        op.groups = [0, 1]
        op.groups = None
        assert op.groups is None

    def test_operations_changing_the_term_count_drop_the_groups(self):
        base = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        for mutate in (
            lambda op: op.__iadd__(BosonOperator.one()),
            lambda op: op.__isub__(BosonOperator.one()),
        ):
            op = copy.deepcopy(base)
            op.groups = [0, 0]
            mutate(op)
            assert not op.has_groups()
            assert len(op.get_coeffs()) == 3

    def test_unpickling_validates_groups(self):
        """The pickled groups are the one external route around the ``groups`` setter."""
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        with pytest.raises(ValueError, match="one group index per term"):
            op.__setstate__([0, 1, 2])
        assert not op.has_groups()
        op.__setstate__([0, 1])
        assert op.groups == [0, 1]

    def test_repeated_modes_survive_grouping(self):
        """A mode may repeat within a term, so group handling must not collapse such terms."""
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 0.5, (cre(0), ann(0)): 1.0})
        op.groups = [1, 0]
        ordered = group_order(op)
        assert ordered.equiv(op)
        assert sorted(ordered.groups) == ordered.groups
        assert ordered.has_groups()


class TestSplitOutGroups(unittest.TestCase):
    """Splitting a Bose-Hubbard-shaped operator into per-bond and diagonal groups."""

    @staticmethod
    def _grouped() -> BosonOperator:
        """A 3-site chain: one group per bond plus one group for the on-site interaction.

        The groups are derived from the *actual* term order rather than being written out
        positionally, because ``from_dict`` gives no ordering guarantee.
        """
        op = BosonOperator.from_dict(
            {
                (cre(0), ann(1)): 1.0,
                (cre(1), ann(0)): 1.0,
                (cre(1), ann(2)): 2.0,
                (cre(2), ann(1)): 2.0,
                (cre(0), cre(0), ann(0), ann(0)): 0.5,
            }
        )
        bonds = [(0, 1), (1, 2)]
        groups = []
        for term, _ in op.iter_terms():
            modes = {mode for _, mode in term}
            groups.append(len(bonds) if len(modes) == 1 else bonds.index(tuple(sorted(modes))))
        op.groups = groups
        return op

    def test_splits_into_all_groups(self):
        op = self._grouped()
        groups = op.split_out_groups()
        assert len(groups) == 3
        assert [len(g) for g in groups] == [2, 2, 1]
        # every group is separately Hermitian, which is exactly why one groups this way
        assert all(g.is_hermitian() for g in groups)

    def test_selection(self):
        op = self._grouped()
        selected = op.split_out_groups(group_indices=[2])
        assert len(selected) == 1
        assert selected[0].equiv(BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 0.5}))

    def test_selection_order_is_honored(self):
        op = self._grouped()
        selected = op.split_out_groups(group_indices=[2, 0])
        assert [len(g) for g in selected] == [1, 2]

    def test_duplicate_request_is_returned_once_per_occurrence(self):
        op = self._grouped()
        selected = op.split_out_groups(group_indices=[1, 1, 1])
        assert len(selected) == 3
        assert all(g.equiv(selected[0]) for g in selected)

    def test_empty_request(self):
        op = self._grouped()
        selected = op.split_out_groups(group_indices=[])
        assert selected == []

    def test_groups_reassemble_to_the_whole(self):
        op = self._grouped()
        total = BosonOperator.zero()
        for group in op.split_out_groups():
            total += group
        assert total.simplify().equiv(op)

    def test_out_of_range_index_yields_an_empty_operator(self):
        """Requesting a group index no term carries is not an error; it yields the zero operator.

        This is the same convention as :func:`group_coeff_means` weighing an unused index ``0.0``:
        an absent group is empty, not invalid.
        """
        op = self._grouped()
        selected = op.split_out_groups(group_indices=[7])
        assert len(selected) == 1
        assert len(selected[0]) == 0
        assert selected[0].equiv(BosonOperator.zero())


class TestArithmeticInvariants(unittest.TestCase):
    """Arithmetic sanity, with an eye on the bosonic composition rules."""

    def test_composition_is_left_multiplication(self):
        b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
        b = BosonOperator.from_dict({(ann(0),): 1.0})
        # `a & b` is B A, so `b_dag & b` == b b_dag
        assert (b_dag & b).equiv(BosonOperator.from_dict({(ann(0), cre(0)): 1.0}))
        assert (b_dag @ b).equiv(BosonOperator.from_dict({(cre(0), ann(0)): 1.0}))

    def test_composition_order_matters_on_the_same_mode(self):
        """``b b^dagger != b^dagger b`` -- they differ by exactly the identity."""
        b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
        b = BosonOperator.from_dict({(ann(0),): 1.0})
        forward = (b_dag @ b).simplify()
        backward = (b @ b_dag).simplify()
        assert not forward.equiv(backward)
        assert (backward - forward).simplify().equiv(BosonOperator.one())

    def test_number_operator_squared(self):
        """``n^2 = b^dag b^dag b b + n``, the identity behind ``n(n-1)``."""
        n_squared = (_num(0) @ _num(0)).simplify()
        expected = BosonOperator.from_dict(
            {(cre(0), cre(0), ann(0), ann(0)): 1.0, (cre(0), ann(0)): 1.0}
        )
        assert n_squared.equiv(expected)

    def test_max_rank_depends_on_form(self):
        op = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        assert op.max_rank() == 2
        # normal ordering emits a constant term, but the longest term is unchanged
        assert op.normal_ordered().simplify().max_rank() == 2

    def test_scalar_arithmetic(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 2.0})
        assert (op / 2).simplify().equiv(BosonOperator.from_dict({(cre(0), cre(0)): 1.0}))
        assert (-op).simplify().equiv(BosonOperator.from_dict({(cre(0), cre(0)): -2.0}))
        assert (op - op).simplify().equiv(BosonOperator.zero())

    def test_getters_round_trip(self):
        coeffs = [1.0, -2.0j]
        actions = [True, True, False]
        modes = [0, 0, 1]
        boundaries = [0, 2, 3]
        op = BosonOperator(coeffs, actions, modes, boundaries)
        assert np.allclose(op.get_coeffs(), coeffs)
        assert op.get_actions() == actions
        assert op.get_modes() == modes
        assert op.get_boundaries() == boundaries
        assert op.get_support() == {0, 1}
        assert BosonOperator(*op.__getnewargs__()) == op

    def test_repr_round_trip(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(1)): 1.5, (): -0.5j})
        assert op.equiv(eval(repr(op)))

    def test_not_iterable_directly(self):
        op = BosonOperator.one()
        with pytest.raises(TypeError):
            iter(op)
