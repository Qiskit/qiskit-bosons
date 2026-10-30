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

"""InitializeModes tests.

The API-level fingerprint of bosons having no exclusion principle is that this gate takes a list of
occupation *numbers* rather than one flag per mode. These tests pin that, plus the validation and
value semantics around it.
"""

from __future__ import annotations

import unittest

import pytest
from qiskit_bosons.circuit import BosonicGate
from qiskit_bosons.circuit.library import InitializeModes


class TestConstruction(unittest.TestCase):
    def test_occupation_is_integers_not_flags(self):
        """A mode may hold several bosons, so the entries are counts."""
        gate = InitializeModes([2, 0, 3])
        assert gate.occupation == (2, 0, 3)
        assert all(isinstance(level, int) for level in gate.occupation)

    def test_booleans_are_coerced_to_integers(self):
        """A fermionic-style flag list still works, but becomes 0/1 integers."""
        gate = InitializeModes([True, False, True])
        assert gate.occupation == (1, 0, 1)
        assert all(type(level) is int for level in gate.occupation)

    def test_floats_are_truncated_to_integers(self):
        gate = InitializeModes([2.0, 0.0])
        assert gate.occupation == (2, 0)

    def test_occupation_is_an_immutable_tuple(self):
        gate = InitializeModes([1, 2])
        assert isinstance(gate.occupation, tuple)
        with pytest.raises(TypeError):
            gate.occupation[0] = 5  # type: ignore[index]

    def test_num_modes_is_the_occupation_length(self):
        gate = InitializeModes([0, 1, 2, 3])
        assert gate.num_modes == 4

    def test_is_a_bosonic_gate(self):
        assert isinstance(InitializeModes([1]), BosonicGate)
        assert InitializeModes([1]).name == "InitializeModes"

    def test_empty_occupation(self):
        gate = InitializeModes([])
        assert gate.occupation == ()
        assert gate.num_modes == 0
        assert gate.max_occupation == 0

    def test_high_occupation_is_allowed(self):
        """No cutoff lives on the gate, so an arbitrarily large occupation constructs fine."""
        gate = InitializeModes([100])
        assert gate.max_occupation == 100

    def test_negative_occupation_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            InitializeModes([1, -1])

    def test_negative_occupation_raises_even_alone(self):
        with pytest.raises(ValueError, match="non-negative"):
            InitializeModes([-3])

    def test_error_message_shows_the_offending_list(self):
        with pytest.raises(ValueError, match=r"\[0, -2\]"):
            InitializeModes([0, -2])


class TestMaxOccupation(unittest.TestCase):
    """The truncation level this gate demands of a downstream encoding."""

    def test_is_the_largest_level(self):
        assert InitializeModes([1, 0, 4, 2]).max_occupation == 4

    def test_vacuum_demands_nothing(self):
        assert InitializeModes([0, 0, 0]).max_occupation == 0

    def test_empty_defaults_to_zero(self):
        assert InitializeModes([]).max_occupation == 0

    def test_is_not_the_total(self):
        """It is a per-mode bound, not the total boson number."""
        gate = InitializeModes([1, 1, 1])
        assert gate.max_occupation == 1
        assert sum(gate.occupation) == 3


class TestVacuum(unittest.TestCase):
    def test_all_zero(self):
        gate = InitializeModes.vacuum(4)
        assert gate.occupation == (0, 0, 0, 0)
        assert gate.num_modes == 4

    def test_equals_the_explicit_construction(self):
        assert InitializeModes.vacuum(3) == InitializeModes([0, 0, 0])

    def test_zero_modes(self):
        assert InitializeModes.vacuum(0).occupation == ()


class TestValueSemantics(unittest.TestCase):
    def test_equal_occupations_compare_equal(self):
        assert InitializeModes([1, 0, 2]) == InitializeModes([1, 0, 2])

    def test_different_occupations_compare_unequal(self):
        assert InitializeModes([1, 0]) != InitializeModes([0, 1])
        assert InitializeModes([1]) != InitializeModes([2])

    def test_different_lengths_compare_unequal(self):
        assert InitializeModes([1, 0]) != InitializeModes([1])

    def test_boolean_and_integer_forms_compare_equal(self):
        """Because the flags are coerced, these are the same gate."""
        assert InitializeModes([True, False]) == InitializeModes([1, 0])

    def test_comparison_against_another_type_is_not_implemented(self):
        assert InitializeModes([1]).__eq__("not a gate") is NotImplemented
        assert InitializeModes([1]) != "not a gate"

    def test_hash_matches_equality(self):
        assert hash(InitializeModes([1, 0, 2])) == hash(InitializeModes([1, 0, 2]))
        assert hash(InitializeModes([1, 0])) != hash(InitializeModes([0, 1]))

    def test_usable_as_a_dict_key(self):
        counts = {InitializeModes([1, 0]): "a"}
        counts[InitializeModes([1, 0])] = "b"
        assert len(counts) == 1
        assert counts[InitializeModes([1, 0])] == "b"

    def test_hash_is_type_qualified(self):
        """The type name enters the hash, so a same-shaped tuple does not collide by construction."""
        assert hash(InitializeModes([1, 0])) != hash((1, 0))
