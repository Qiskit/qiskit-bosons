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

"""The shared operator contract, mirroring ``qiskit-fermions``.

Upstream exercises this mixin against four operator types. Here there is a single one,
:class:`.BosonOperator`, but the mixin is kept in the same shape so that a future second operator
type inherits the contract rather than copying it, and so that the two packages stay diffable.
"""

from __future__ import annotations

import copy
from typing import Protocol

import pytest


class _OperatorClassProvider(Protocol):
    """The single hook a test class provides: which operator type to exercise."""

    @staticmethod
    def get_class() -> type: ...


class OperatorContractTests:
    """Tests of the arithmetic and iteration contract shared by every operator type.

    Every test here is expressed purely in terms of ``zero``, ``one`` and the *identity* term
    ``{(): coeff}``, so none of them needs to know how a particular operator family spells its
    generators. Anything specific to the bosonic algebra (repeated modes, contractions, truncation)
    stays in ``test_boson_operator.py``.

    Subclasses supply :meth:`get_class` and inherit these.
    """

    # The subclass provides this; declared for the type checker's benefit.
    get_class: staticmethod

    def test_zero(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.zero()
        assert op == cls.from_dict({})

    def test_one(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.one()
        assert op == cls.from_dict({(): 1})

    def test_add(self: _OperatorClassProvider):
        cls = self.get_class()
        one = cls.one()
        two = cls.from_dict({(): 2})
        three = one + two
        assert three.equiv(cls.from_dict({(): 3}))

    def test_iadd(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.one()
        two = cls.from_dict({(): 2})
        op += two
        assert op.equiv(cls.from_dict({(): 3}))

    def test_sub(self: _OperatorClassProvider):
        cls = self.get_class()
        one = cls.one()
        two = cls.from_dict({(): 2})
        new_one = two - one
        assert new_one.equiv(one)

    def test_isub(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.from_dict({(): 2})
        one = cls.one()
        op -= one
        assert op.equiv(one)

    def test_mul(self: _OperatorClassProvider):
        cls = self.get_class()
        one = cls.one()
        three = one * 3
        assert three.equiv(cls.from_dict({(): 3}))

    def test_rmul(self: _OperatorClassProvider):
        cls = self.get_class()
        one = cls.one()
        three = 3 * one
        assert three.equiv(cls.from_dict({(): 3}))

    def test_imul(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.one()
        op *= 3
        assert op.equiv(cls.from_dict({(): 3}))

    def test_div(self: _OperatorClassProvider):
        cls = self.get_class()
        three = cls.from_dict({(): 3})
        one_half = three / 2.0
        assert one_half.equiv(cls.from_dict({(): 1.5}))

    def test_idiv(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.from_dict({(): 3})
        op /= 2.0
        assert op.equiv(cls.from_dict({(): 1.5}))

    def test_neg(self: _OperatorClassProvider):
        cls = self.get_class()
        one = cls.one()
        assert (-one).equiv(cls.from_dict({(): -1}))

    def test_iter(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.one()
        assert list(op.iter_terms()) == [([], 1)]

    def test_iter_with_groups(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.one()
        op.groups = [0]
        assert list(op.iter_terms_with_groups()) == [([], 1, 0)]

    def test_equiv(self: _OperatorClassProvider):
        cls = self.get_class()
        op = cls.from_dict({(): 1e-7})
        zero = cls.zero()
        assert not op.equiv(zero)
        assert op.equiv(zero, 1e-6)
        assert not op.equiv(zero, 1e-8)

    def test_matmul(self: _OperatorClassProvider):
        """``@`` composes two operators."""
        cls = self.get_class()
        two = cls.from_dict({(): 2})
        three = cls.from_dict({(): 3})
        assert (two @ three).equiv(cls.from_dict({(): 6}))

    def test_imatmul(self: _OperatorClassProvider):
        """``@=`` composes in place, matching ``@``."""
        cls = self.get_class()
        op = cls.from_dict({(): 2})
        op @= cls.from_dict({(): 3})
        assert op.equiv(cls.from_dict({(): 6}))

    def test_matmul_leaves_operands_untouched(self: _OperatorClassProvider):
        """The out-of-place form must not mutate either operand."""
        cls = self.get_class()
        two = cls.from_dict({(): 2})
        three = cls.from_dict({(): 3})
        _ = two @ three
        assert two.equiv(cls.from_dict({(): 2}))
        assert three.equiv(cls.from_dict({(): 3}))

    def test_deepcopy_is_independent(self: _OperatorClassProvider):
        """``deepcopy`` must copy the group assignment, not share it."""
        cls = self.get_class()
        op = cls.from_dict({(): 1})
        op.groups = [0]

        clone = copy.deepcopy(op)
        clone.groups = [5]
        clone += cls.from_dict({(): 1})

        assert op.groups == [0], "mutating the copy's groups changed the original"
        assert op.equiv(cls.from_dict({(): 1})), "mutating the copy changed the original's terms"

    def test_pow_modulo_is_accepted(self: _OperatorClassProvider):
        """``pow(op, exponent, None)`` must behave like the two-argument form."""
        cls = self.get_class()
        two = cls.from_dict({(): 2})
        assert two.__pow__(2, None).equiv(cls.from_dict({(): 4}))
        assert pow(two, 2, None).equiv(cls.from_dict({(): 4}))

    def test_set_groups_err(self: _OperatorClassProvider):
        """The ``groups`` setter rejects an array that does not hold one index per term."""
        cls = self.get_class()
        op = cls.one() + cls.one()
        assert len(op.get_coeffs()) == 2

        with pytest.raises(ValueError, match="expected one group index per term"):
            op.groups = [0]

        with pytest.raises(ValueError, match="expected one group index per term"):
            op.groups = [0, 0, 1]

        assert not op.has_groups(), "a rejected assignment changed the operator"

        op.groups = [0, 1]
        op.groups = None
        assert not op.has_groups(), "clearing the group indices must always be allowed"

    def test_setstate_rejects_a_length_mismatch(self: _OperatorClassProvider):
        """Unpickling validates the group indices rather than trusting them."""
        cls = self.get_class()
        op = cls.one()

        with pytest.raises(ValueError, match="expected one group index per term"):
            op.__setstate__([0, 1])

        assert not op.has_groups(), "a rejected payload changed the operator"
