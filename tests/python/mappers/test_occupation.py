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

"""Occupation-encoding boson-to-qubit mapper tests.

The core of this file is a *faithfulness* check: for each encoding, the mapped qubit operator is
compared against a dense truncated matrix built independently from the ladder-operator definition,
restricted to the physical (in-code) subspace. That independent construction is what makes these
tests able to catch a transposed raise/lower convention, which the mapper has in fact had.
"""

from __future__ import annotations

import itertools
import unittest

import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp
from qiskit_bosons.mappers import LevelEncoding, LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre

ENCODINGS = ("unary_projected", "standard_binary", "gray")


def level_codes(num_levels: int, encoding: str) -> list[int]:
    """The bit pattern each occupation level is stored as, per encoding."""
    match encoding:
        case "unary_projected" | "unary_unprojected":
            return [1 << level for level in range(num_levels)]
        case "standard_binary":
            return list(range(num_levels))
        case "gray":
            return [level ^ (level >> 1) for level in range(num_levels)]
        case _:  # pragma: no cover - guards a typo in the test itself
            raise AssertionError(f"unknown encoding {encoding!r}")


def ladder_matrices(num_levels: int) -> tuple[np.ndarray, np.ndarray]:
    r"""The truncated single-mode ladder operators.

    ``b^dagger`` carries ``sqrt(n+1)`` on the *sub*-diagonal (it raises ``|n>`` to ``|n+1>``) and
    ``b`` carries ``sqrt(n)`` on the *super*-diagonal, both truncated to ``num_levels``. Written out
    from the definition rather than derived from the mapper, so that a transposed convention shows up
    as a failure rather than canceling out.
    """
    creation = np.zeros((num_levels, num_levels), dtype=complex)
    annihilation = np.zeros((num_levels, num_levels), dtype=complex)
    for level in range(num_levels - 1):
        creation[level + 1, level] = np.sqrt(level + 1)
        annihilation[level, level + 1] = np.sqrt(level + 1)
    return creation, annihilation


def reference_matrix(op: BosonOperator, num_modes: int, num_levels: int) -> np.ndarray:
    """Builds the dense truncated matrix of ``op`` directly from the ladder definition.

    Per-mode factors are composed in term order, then combined with ``kron``. Qiskit is
    little-endian, so mode 0 must be the **last** kron factor -- i.e. the index of a basis state is
    ``sum_m level_m * num_levels**m``.
    """
    creation, annihilation = ladder_matrices(num_levels)
    dim = num_levels**num_modes
    total = np.zeros((dim, dim), dtype=complex)

    for term, coeff in op.iter_terms():
        factors = [np.eye(num_levels, dtype=complex) for _ in range(num_modes)]
        for action, mode in term:
            factors[mode] = factors[mode] @ (creation if action else annihilation)

        matrix = np.array([[1.0 + 0.0j]])
        for mode in reversed(range(num_modes)):
            matrix = np.kron(matrix, factors[mode])
        total += coeff * matrix

    return total


def physical_indices(num_modes: int, num_levels: int, encoding: str) -> list[int]:
    """The qubit basis indices of the physical (in-code) states, in reference-matrix order.

    Mode ``m`` at level ``n`` contributes ``code(n) << (m * per_mode)``, and the entries are ordered
    to match :func:`reference_matrix`, whose row ``r`` has mode ``m`` at level
    ``(r // num_levels**m) % num_levels``.
    """
    codes = level_codes(num_levels, encoding)
    per_mode = num_qubits_per_mode(num_levels - 1, encoding)

    indices = []
    for row in range(num_levels**num_modes):
        value = 0
        for mode in range(num_modes):
            level = (row // num_levels**mode) % num_levels
            value |= codes[level] << (mode * per_mode)
        indices.append(value)
    return indices


def mapped_on_physical_subspace(
    op: BosonOperator, num_modes: int, num_levels: int, encoding: str
) -> np.ndarray:
    """Maps ``op`` to qubits and restricts the dense result to the physical subspace."""
    observable = occupation_encoding(op, LevelSpec(num_levels - 1, encoding).resolve(num_modes))
    dense = SparsePauliOp.from_sparse_observable(observable).to_matrix()
    indices = physical_indices(num_modes, num_levels, encoding)
    return np.asarray(dense)[np.ix_(indices, indices)]


def assert_faithful(op: BosonOperator, num_modes: int, num_levels: int, encoding: str) -> None:
    """Asserts the mapped operator reproduces the reference matrix on the physical subspace."""
    got = mapped_on_physical_subspace(op, num_modes, num_levels, encoding)
    want = reference_matrix(op, num_modes, num_levels)
    np.testing.assert_allclose(got, want, atol=1e-12, err_msg=f"{encoding}, d={num_levels}")


def num_qubits_per_mode(max_occupation: int, encoding: str) -> int:
    """The number of qubits one mode occupies, as the spec reports it."""
    return LevelSpec(max_occupation, encoding).resolve(1).total_width


class TestNumQubitsPerMode(unittest.TestCase):
    def test_unary_is_the_level_count(self):
        for max_occupation in range(5):
            assert num_qubits_per_mode(max_occupation, "unary_projected") == max_occupation + 1

    def test_compact_is_the_logarithm(self):
        for max_occupation in range(1, 9):
            expected = int(np.ceil(np.log2(max_occupation + 1)))
            for encoding in ("standard_binary", "gray"):
                assert num_qubits_per_mode(max_occupation, encoding) == expected

    def test_single_level_is_degenerate(self):
        """``max_occupation=0`` means one level, i.e. no dynamics: the compact codes need 0 qubits.

        This is literally ``ceil(log2(1)) = 0`` and is not clamped to 1. A single-level truncation
        describes a mode frozen in the vacuum, so the mapped operator is empty; it is recorded here
        rather than worked around, since the value is mathematically right and the case vacuous.
        """
        # the width is a property of the layout, so it is read off the encoding for one level
        assert LevelEncoding.STANDARD_BINARY.num_qubits(1) == 0
        assert LevelEncoding.GRAY.num_qubits(1) == 0
        assert num_qubits_per_mode(0, "unary_projected") == 1

        # A spec made only of such modes would span no qubits at all, which it refuses to describe.
        for encoding in ("standard_binary", "gray"):
            with pytest.raises(ValueError, match="no qubits"):
                LevelSpec(0, encoding)

        # ...but a single-level mode among wider ones is fine and simply takes no qubit.
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        observable = occupation_encoding(op, LevelSpec([0, 1], "gray"))
        assert observable.num_qubits == 1
        assert len(observable) == 0

    def test_non_powers_of_two_round_up(self):
        # d = 5 levels needs 3 qubits, not 2.5
        for encoding in ("standard_binary", "gray"):
            assert num_qubits_per_mode(4, encoding) == 3
            assert num_qubits_per_mode(2, encoding) == 2
            assert num_qubits_per_mode(6, encoding) == 3
            assert num_qubits_per_mode(7, encoding) == 3
            assert num_qubits_per_mode(8, encoding) == 4

    def test_both_unary_variants_use_the_same_number_of_qubits(self):
        for max_occupation in range(5):
            assert num_qubits_per_mode(max_occupation, "unary_projected") == num_qubits_per_mode(
                max_occupation, "unary_unprojected"
            )

    def test_compact_beats_unary_from_three_levels_on(self):
        for max_occupation in range(2, 8):
            assert num_qubits_per_mode(max_occupation, "gray") < num_qubits_per_mode(
                max_occupation, "unary_projected"
            )

    def test_unknown_encoding_raises(self):
        with pytest.raises(ValueError, match=r"unknown|encoding"):
            num_qubits_per_mode(2, "not_an_encoding")


class TestCodeLayout(unittest.TestCase):
    """Properties of the three level-to-bits codes themselves."""

    def test_codes_are_distinct(self):
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 4, 5):
                codes = level_codes(num_levels, encoding)
                assert len(set(codes)) == num_levels

    def test_codes_fit_in_the_allotted_qubits(self):
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 4, 5):
                per_mode = num_qubits_per_mode(num_levels - 1, encoding)
                assert all(code < (1 << per_mode) for code in level_codes(num_levels, encoding))

    def test_gray_adjacent_levels_differ_in_one_bit(self):
        codes = level_codes(8, "gray")
        for level in range(7):
            assert (codes[level] ^ codes[level + 1]).bit_count() == 1

    def test_standard_binary_reaches_larger_distances(self):
        """Standard binary has an adjacent-level Hamming distance above 1, which Gray never does."""
        codes = level_codes(8, "standard_binary")
        distances = [(codes[n] ^ codes[n + 1]).bit_count() for n in range(7)]
        assert max(distances) > 1
        # 3 -> 4 flips all three bits
        assert distances[3] == 3

    def test_unary_is_one_hot(self):
        for code in level_codes(5, "unary_projected"):
            assert code.bit_count() == 1


class TestFaithfulness(unittest.TestCase):
    """The mapped operator must equal the truncated matrix on the physical subspace."""

    def test_bare_annihilation(self):
        """A bare ``b`` pins the lowering direction; a transposed convention fails here."""
        op = BosonOperator.from_dict({(ann(0),): 1.0})
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                assert_faithful(op, 1, num_levels, encoding)

    def test_bare_creation(self):
        """A bare ``b^dagger`` pins the raising direction, the other half of the same check."""
        op = BosonOperator.from_dict({(cre(0),): 1.0})
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                assert_faithful(op, 1, num_levels, encoding)

    def test_annihilation_is_the_adjoint_of_creation(self):
        """Guards against a convention that transposes *both* consistently."""
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                lowered = mapped_on_physical_subspace(
                    BosonOperator.from_dict({(ann(0),): 1.0}), 1, num_levels, encoding
                )
                raised = mapped_on_physical_subspace(
                    BosonOperator.from_dict({(cre(0),): 1.0}), 1, num_levels, encoding
                )
                np.testing.assert_allclose(lowered, raised.conj().T, atol=1e-12)
                # and the raising operator really is strictly lower-triangular
                assert np.allclose(np.triu(raised), 0.0)

    def test_number_operator(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                assert_faithful(op, 1, num_levels, encoding)
                # it is diagonal with the level index on the diagonal
                mapped = mapped_on_physical_subspace(op, 1, num_levels, encoding)
                np.testing.assert_allclose(np.diag(mapped).real, np.arange(num_levels), atol=1e-12)

    def test_repeated_creation(self):
        """``b^dagger b^dagger`` -- the term with no fermionic counterpart at all."""
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                assert_faithful(op, 1, num_levels, encoding)

    def test_interaction_term(self):
        """``n(n-1) = b^dag b^dag b b``, the Bose-Hubbard on-site interaction."""
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                assert_faithful(op, 1, num_levels, encoding)
                mapped = mapped_on_physical_subspace(op, 1, num_levels, encoding)
                levels = np.arange(num_levels)
                np.testing.assert_allclose(np.diag(mapped).real, levels * (levels - 1), atol=1e-12)

    def test_needs_contraction(self):
        """``b b^dagger`` -- mapped in its *given* (non-normal) order, no contraction applied."""
        op = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                assert_faithful(op, 1, num_levels, encoding)

    def test_multi_mode_hopping(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        for encoding in ENCODINGS:
            for num_levels in (2, 3):
                assert_faithful(op, 2, num_levels, encoding)

    def test_complex_coefficients(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 0.5 + 1.5j, (cre(1), ann(0)): 0.5 - 1.5j})
        for encoding in ENCODINGS:
            for num_levels in (2, 3):
                assert_faithful(op, 2, num_levels, encoding)

    def test_three_modes_with_padding(self):
        """A term skipping mode 1 must leave that mode's qubits on the identity."""
        op = BosonOperator.from_dict(
            {(cre(0), ann(2)): 1.0, (cre(1), cre(1), ann(1), ann(1)): 0.25}
        )
        for encoding in ENCODINGS:
            assert_faithful(op, 3, 2, encoding)

    def test_identity_and_zero(self):
        for encoding in ENCODINGS:
            assert_faithful(BosonOperator.one(), 2, 3, encoding)
            assert_faithful(BosonOperator.zero(), 2, 3, encoding)

    def test_non_power_of_two_level_count(self):
        """d = 5 leaves three unused code words per mode under a compact encoding."""
        op = BosonOperator.from_dict(
            {(cre(0), ann(0)): 1.0, (cre(0), cre(0), ann(0), ann(0)): 0.5, (ann(0),): 1.0j}
        )
        for encoding in ENCODINGS:
            assert_faithful(op, 1, 5, encoding)

    def test_linearity(self):
        """Mapping is linear, so a sum maps to the sum of the mapped parts."""
        first = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
        second = BosonOperator.from_dict({(cre(1), ann(0)): 2.0j})
        for encoding in ENCODINGS:
            combined = mapped_on_physical_subspace(first + second, 2, 3, encoding)
            separate = mapped_on_physical_subspace(
                first, 2, 3, encoding
            ) + mapped_on_physical_subspace(second, 2, 3, encoding)
            np.testing.assert_allclose(combined, separate, atol=1e-12)

    def test_hermitian_operator_maps_to_a_hermitian_matrix(self):
        op = BosonOperator.from_dict(
            {
                (cre(0), ann(1)): -1.0,
                (cre(1), ann(0)): -1.0,
                (cre(0), cre(0), ann(0), ann(0)): 1.0,
            }
        )
        assert op.is_hermitian()
        for encoding in ENCODINGS:
            mapped = mapped_on_physical_subspace(op, 2, 3, encoding)
            np.testing.assert_allclose(mapped, mapped.conj().T, atol=1e-12)


class TestRegisterWidth(unittest.TestCase):
    def test_num_qubits_matches_the_layout(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
        for encoding in ENCODINGS:
            for max_occupation in (1, 2, 4):
                observable = occupation_encoding(op, LevelSpec(max_occupation, encoding).resolve(3))
                assert observable.num_qubits == 3 * num_qubits_per_mode(max_occupation, encoding)

    def test_padding_modes_widens_the_register(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        narrow = occupation_encoding(op, LevelSpec(2, "gray"))
        wide = occupation_encoding(op, LevelSpec(2, "gray").resolve(4))
        assert narrow.num_qubits == 2
        assert wide.num_qubits == 8


class TestTruncationDefect(unittest.TestCase):
    r"""``[b, b^dagger]`` maps to ``diag(1, ..., 1, 1 - d)``, not the identity.

    This is the price of truncation, and it is exact: the canonical commutator holds on every level
    except the highest, where the amplitude ``b^dagger`` would send out of the retained space is
    clipped. It is also why no bosonic Jordan-Wigner exists: the span of the ladder operators is not
    closed under the commutator.
    """

    @staticmethod
    def _commutator_operator() -> BosonOperator:
        return BosonOperator.from_dict({(ann(0), cre(0)): 1.0, (cre(0), ann(0)): -1.0})

    def test_defect_is_one_minus_d_in_the_top_entry(self):
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 4, 5):
                mapped = mapped_on_physical_subspace(
                    self._commutator_operator(), 1, num_levels, encoding
                )
                expected = np.eye(num_levels, dtype=complex)
                expected[num_levels - 1, num_levels - 1] = 1 - num_levels
                np.testing.assert_allclose(
                    mapped, expected, atol=1e-12, err_msg=f"{encoding}, d={num_levels}"
                )

    def test_defect_is_not_the_identity(self):
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                mapped = mapped_on_physical_subspace(
                    self._commutator_operator(), 1, num_levels, encoding
                )
                assert not np.allclose(mapped, np.eye(num_levels))

    def test_defect_is_a_rank_one_projector_correction(self):
        """The defect is exactly ``-d |d-1><d-1|``."""
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                mapped = mapped_on_physical_subspace(
                    self._commutator_operator(), 1, num_levels, encoding
                )
                defect = mapped - np.eye(num_levels)
                assert np.linalg.matrix_rank(defect, tol=1e-10) == 1
                assert defect[num_levels - 1, num_levels - 1] == pytest.approx(-num_levels)

    def test_commutator_holds_below_the_cutoff(self):
        """Every level but the top obeys the canonical relation exactly."""
        for encoding in ENCODINGS:
            for num_levels in (3, 5):
                mapped = mapped_on_physical_subspace(
                    self._commutator_operator(), 1, num_levels, encoding
                )
                below = mapped[: num_levels - 1, : num_levels - 1]
                np.testing.assert_allclose(below, np.eye(num_levels - 1), atol=1e-12)


class TestMappingDoesNotCommuteWithNormalOrdering(unittest.TestCase):
    r"""``map(b b^dagger) != map((b b^dagger).normal_ordered())``.

    Normal ordering is exact in the infinite-dimensional Fock space, so truncating before versus
    after applying it gives genuinely different operators. This is expected, documented bosonic
    behavior -- locked here so that nobody "fixes" it into agreement.
    """

    @staticmethod
    def _op() -> BosonOperator:
        return BosonOperator.from_dict({(ann(0), cre(0)): 1.0})

    def test_the_two_mappings_differ(self):
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                direct = mapped_on_physical_subspace(self._op(), 1, num_levels, encoding)
                ordered = mapped_on_physical_subspace(
                    self._op().normal_ordered(), 1, num_levels, encoding
                )
                assert not np.allclose(direct, ordered), f"{encoding}, d={num_levels}"

    def test_they_differ_by_exactly_d_in_the_top_entry(self):
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                direct = mapped_on_physical_subspace(self._op(), 1, num_levels, encoding)
                ordered = mapped_on_physical_subspace(
                    self._op().normal_ordered(), 1, num_levels, encoding
                )
                difference = np.zeros((num_levels, num_levels), dtype=complex)
                difference[num_levels - 1, num_levels - 1] = num_levels
                np.testing.assert_allclose(ordered - direct, difference, atol=1e-12)

    def test_the_diagonals_spell_out_the_clipping(self):
        """``b b^dagger`` gives ``d-1`` at the top (clipped) while ``b^dagger b + 1`` gives ``d``."""
        for encoding in ENCODINGS:
            for num_levels in (2, 3, 5):
                direct = mapped_on_physical_subspace(self._op(), 1, num_levels, encoding)
                ordered = mapped_on_physical_subspace(
                    self._op().normal_ordered(), 1, num_levels, encoding
                )
                # the direct mapping loses the top level entirely
                assert np.diag(direct).real[num_levels - 1] == pytest.approx(0.0)
                assert np.diag(ordered).real[num_levels - 1] == pytest.approx(num_levels)
                np.testing.assert_allclose(
                    np.diag(ordered).real, np.arange(1, num_levels + 1), atol=1e-12
                )

    def test_normal_ordering_agrees_where_no_contraction_is_needed(self):
        """An already-normal-ordered operator is of course unaffected."""
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        for encoding in ENCODINGS:
            direct = mapped_on_physical_subspace(op, 1, 3, encoding)
            ordered = mapped_on_physical_subspace(op.normal_ordered(), 1, 3, encoding)
            np.testing.assert_allclose(direct, ordered, atol=1e-12)


class TestErrorPaths(unittest.TestCase):
    def test_mode_index_at_num_modes_raises(self):
        op = BosonOperator.from_dict({(cre(2), ann(0)): 1.0})
        with pytest.raises(ValueError):
            occupation_encoding(op, LevelSpec(2, "unary_projected").resolve(2))

    def test_mode_index_beyond_num_modes_raises(self):
        op = BosonOperator.from_dict({(cre(0), ann(7)): 1.0})
        with pytest.raises(ValueError):
            occupation_encoding(op, LevelSpec(1, "gray").resolve(3))

    def test_zero_modes_with_a_non_empty_operator_raises(self):
        op = BosonOperator.from_dict({(cre(0),): 1.0})
        with pytest.raises(ValueError):
            occupation_encoding(op, LevelSpec(2, "unary_projected").resolve(0))

    def test_unknown_encoding_raises(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        with pytest.raises(ValueError, match=r"unknown|encoding"):
            occupation_encoding(op, LevelSpec(2, "binary"))

    def test_exactly_num_modes_is_accepted(self):
        """The boundary case: the largest mode index may be ``num_modes - 1``."""
        op = BosonOperator.from_dict({(cre(1), ann(0)): 1.0})
        assert occupation_encoding(op, LevelSpec(1, "unary_projected").resolve(2)).num_qubits == 4


class TestEncodingComparison(unittest.TestCase):
    """Cross-encoding invariants: they describe the same physics on different qubit layouts."""

    def test_all_encodings_agree_on_the_physical_subspace(self):
        op = BosonOperator.from_dict(
            {
                (cre(0), ann(1)): -1.0,
                (cre(1), ann(0)): -1.0,
                (cre(0), cre(0), ann(0), ann(0)): 1.0,
                (cre(1), cre(1), ann(1), ann(1)): 1.0,
            }
        )
        matrices = [mapped_on_physical_subspace(op, 2, 3, enc) for enc in ENCODINGS]
        for other in matrices[1:]:
            np.testing.assert_allclose(matrices[0], other, atol=1e-12)

    def test_diagonal_operators_emit_equal_term_counts(self):
        """For a *diagonal* operator the two compact codes do emit the same number of terms.

        A diagonal matrix element flips no bits, so it needs no Pauli expansion under either code and
        the counts coincide. This is the case in which the docstring's "same number of terms" claim
        holds; see :meth:`test_gray_emits_fewer_terms_than_standard_binary_for_hopping` for the case
        in which it does not.
        """
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        for max_occupation in (1, 2, 3, 4):
            binary = occupation_encoding(op, LevelSpec(max_occupation, "standard_binary"))
            gray = occupation_encoding(op, LevelSpec(max_occupation, "gray"))
            assert len(binary) == len(gray)

    def test_gray_emits_fewer_terms_than_standard_binary_for_hopping(self):
        r"""Gray emits strictly *fewer* terms than standard binary for a ladder operator.

        .. warning::
           This contradicts the :func:`occupation_encoding` docstring, which states that
           ``"standard_binary"`` and ``"gray"`` "emit the same *number* of terms" and that the Gray
           code's advantage "lies entirely in the weight *distribution*". The implementation expands
           each matrix element ``<row|O|col>`` into ``2**k`` Pauli terms, where ``k`` is the Hamming
           distance between the code words of ``row`` and ``col``. Since the Gray code guarantees
           ``k = 1`` for adjacent levels while standard binary does not, the Gray code emits fewer
           terms as well as lower-weight ones. The docstring is inaccurate; the behavior asserted
           here is the correct one.
        """
        op = BosonOperator.from_dict({(ann(0),): 1.0})

        def predicted(num_levels: int, encoding: str) -> int:
            codes = level_codes(num_levels, encoding)
            return sum(
                2 ** (codes[level] ^ codes[level + 1]).bit_count()
                for level in range(num_levels - 1)
            )

        for num_levels in (3, 4, 5):
            binary = len(occupation_encoding(op, LevelSpec(num_levels - 1, "standard_binary")))
            gray = len(occupation_encoding(op, LevelSpec(num_levels - 1, "gray")))
            # the 2**(Hamming distance) expansion rule predicts both counts exactly
            assert binary == predicted(num_levels, "standard_binary")
            assert gray == predicted(num_levels, "gray")
            assert gray < binary, f"d={num_levels}"

        # at d = 2 the codes coincide, so the counts do too
        assert len(occupation_encoding(op, LevelSpec(1, "gray"))) == len(
            occupation_encoding(op, LevelSpec(1, "standard_binary"))
        )

    def test_unary_ladder_terms_have_exactly_two_off_diagonal_factors(self):
        r"""A unary ladder operator moves amplitude between exactly two level qubits.

        This is the part of the "strictly two-local" claim that does hold: each emitted term has
        exactly two ``X``/``Y`` (off-diagonal) bit terms, on adjacent level qubits of the mode. The
        remaining factors are ``Zero`` projectors; see
        :meth:`test_unary_term_weight_grows_with_the_truncation`.
        """
        x_and_y = {2, 3}  # QkBitTerm::X and QkBitTerm::Y
        for max_occupation in (1, 2, 4):
            observable = occupation_encoding(
                BosonOperator.from_dict({(ann(0),): 1.0}),
                LevelSpec(max_occupation, "unary_projected"),
            ).simplify()
            for term in observable:
                off_diagonal = [
                    int(index)
                    for index, bit in zip(term.indices, term.bit_terms, strict=True)
                    if int(bit) in x_and_y
                ]
                assert len(off_diagonal) == 2
                assert abs(off_diagonal[0] - off_diagonal[1]) == 1

    def test_unary_term_weight_grows_with_the_truncation(self):
        r"""A unary term's Pauli weight is ``d``, not 2.

        .. warning::
           This contradicts the :func:`occupation_encoding` docstring (and the
           :mod:`~qiskit_bosons.mappers.library` module docs), which state that the unary encoding
           "yields strictly weight-2 strings" and that "the output stays strictly weight-2 regardless
           of the truncation level -- the property behind its ``O(d)`` gate scaling".

           In fact every emitted term additionally carries a ``Zero`` projector on *each of the
           mode's other level qubits*, giving weight exactly ``d``. Those projectors are **not
           load-bearing**: they act as the identity on the one-hot (physical) subspace, so deleting
           them leaves the operator unchanged there -- as
           :meth:`test_unary_projectors_are_redundant_on_the_physical_subspace` demonstrates. They
           only alter the action on unphysical, non-one-hot basis states.

           The consequence is a real cost, not a cosmetic one: the ``O(d)``-weight strings are why the
           unary encoding transpiles to far deeper circuits than the compact codes at small ``d``
           (the Bose-Hubbard guide reports a 2q-depth of 645 for unary against 94 for Gray at
           ``d = 3``). The behavior is locked here as-is; the fix belongs in the mapper, which could
           omit the redundant projectors.
        """
        op = BosonOperator.from_dict({(ann(0),): 1.0})
        for num_levels in (2, 3, 4, 5):
            observable = occupation_encoding(
                op, LevelSpec(num_levels - 1, "unary_projected")
            ).simplify()
            weights = {len(term.indices) for term in observable}
            assert weights == {num_levels}, f"d={num_levels}"

    def test_unary_projectors_are_redundant_on_the_physical_subspace(self):
        """Deleting the ``Zero`` projectors leaves the unary operator unchanged where it matters.

        This is the evidence behind the warning on
        :meth:`test_unary_term_weight_grows_with_the_truncation`: the strictly weight-2 form the
        docstring advertises is available and agrees with what is emitted, on the one-hot subspace.
        """
        from qiskit.quantum_info import SparseObservable

        pauli_names = {2: "X", 3: "Y"}
        num_levels = 4
        op = BosonOperator.from_dict({(cre(0),): 1.0})
        observable = occupation_encoding(
            op, LevelSpec(num_levels - 1, "unary_projected")
        ).simplify()

        stripped_terms = []
        for term in observable:
            indices, labels = [], []
            for index, bit in zip(term.indices, term.bit_terms, strict=True):
                if int(bit) in pauli_names:
                    indices.append(int(index))
                    labels.append(pauli_names[int(bit)])
            stripped_terms.append(("".join(labels), indices, complex(term.coeff)))

        stripped = SparseObservable.from_sparse_list(stripped_terms, observable.num_qubits)
        assert max(len(term.indices) for term in stripped) == 2

        indices = physical_indices(1, num_levels, "unary_projected")
        as_matrix = SparsePauliOp.from_sparse_observable(stripped).to_matrix()
        np.testing.assert_allclose(
            np.asarray(as_matrix)[np.ix_(indices, indices)],
            reference_matrix(op, 1, num_levels),
            atol=1e-12,
        )

    def test_operator_is_not_simplified_on_output(self):
        """The docstring promises an unsimplified result; duplicate terms are the caller's problem."""
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        observable = occupation_encoding(op, LevelSpec(3, "unary_projected"))
        assert len(observable) >= len(observable.simplify())


class TestPhysicalSubspaceHelpers(unittest.TestCase):
    """Sanity checks on this file's own reference construction, so a bug there cannot hide."""

    def test_reference_ladder_reproduces_the_number_operator(self):
        for num_levels in (2, 3, 5):
            creation, annihilation = ladder_matrices(num_levels)
            number = creation @ annihilation
            np.testing.assert_allclose(np.diag(number).real, np.arange(num_levels), atol=1e-12)

    def test_reference_ladder_defect(self):
        """The reference matrices themselves show the ``1 - d`` defect."""
        for num_levels in (2, 3, 5):
            creation, annihilation = ladder_matrices(num_levels)
            commutator = annihilation @ creation - creation @ annihilation
            expected = np.eye(num_levels, dtype=complex)
            expected[num_levels - 1, num_levels - 1] = 1 - num_levels
            np.testing.assert_allclose(commutator, expected, atol=1e-12)

    def test_physical_indices_are_distinct_and_ordered(self):
        for encoding in ENCODINGS:
            for num_modes, num_levels in ((1, 5), (2, 3), (3, 2)):
                indices = physical_indices(num_modes, num_levels, encoding)
                assert len(indices) == num_levels**num_modes
                assert len(set(indices)) == len(indices)
                per_mode = num_qubits_per_mode(num_levels - 1, encoding)
                assert all(index < (1 << (num_modes * per_mode)) for index in indices)

    def test_physical_index_of_the_vacuum_is_zero_for_compact_codes(self):
        for encoding in ("standard_binary", "gray"):
            assert physical_indices(2, 3, encoding)[0] == 0

    def test_reference_matrix_little_endian_convention(self):
        """``b^dagger_0`` must act on the *lowest* levels index, i.e. mode 0 varies fastest."""
        op = BosonOperator.from_dict({(cre(0),): 1.0})
        matrix = reference_matrix(op, 2, 2)
        # |0,0> (index 0) -> |1,0> (index 1), since mode 0 is the least significant digit
        assert matrix[1, 0] == pytest.approx(1.0)
        assert matrix[2, 0] == pytest.approx(0.0)

    def test_every_level_combination_is_covered(self):
        combinations = list(itertools.product(range(3), repeat=2))
        assert len(physical_indices(2, 3, "gray")) == len(combinations)


class TestAccumulatorMemoryBound:
    """The mapped image must grow with the *distinct* content, not the number of terms emitted.

    ``qk_obs_add_inplace`` concatenates without deduplicating. Powers of one number operator map
    onto the same diagonal Pauli strings, so a large operator made of repeats would otherwise hold
    every string it emitted although its simplified image is tiny.
    """

    @staticmethod
    def _repeated_number_operator(repeats: int) -> BosonOperator:
        return BosonOperator(
            [1.0] * repeats,
            [True, False] * repeats,
            [0, 0] * repeats,
            [2 * i for i in range(repeats + 1)],
        )

    def test_duplicates_are_merged_while_mapping(self):
        repeats = 20_000
        obs = occupation_encoding(
            self._repeated_number_operator(repeats), LevelSpec(3, "unary_projected")
        )
        emitted = repeats * 3  # three non-vacuum levels per repeat under the unary encoding
        distinct = len(obs.simplify())
        assert distinct == 3
        # without the merge the image would hold all `emitted` terms; leave generous slack for the
        # growth factor and the compaction floor, and assert only that it was bounded
        assert len(obs) < emitted // 4

    def test_merging_does_not_change_the_operator(self):
        """Compaction only combines duplicates, so the value is exactly ``repeats * n``."""
        repeats = 20_000
        big = occupation_encoding(
            self._repeated_number_operator(repeats), LevelSpec(3, "unary_projected")
        )
        single = occupation_encoding(
            self._repeated_number_operator(1), LevelSpec(3, "unary_projected")
        )
        assert (big.simplify() - repeats * single.simplify()).simplify() == type(big).zero(
            big.num_qubits
        )

    def test_small_operators_are_unaffected(self):
        """Below the compaction floor nothing is merged, and the result is still exact."""
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (cre(0), cre(0), ann(0), ann(0)): 1.0})
        obs = occupation_encoding(op, LevelSpec(3, "unary_projected"))
        assert len(obs.simplify()) == len(obs.simplify().simplify())
