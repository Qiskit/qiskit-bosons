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

"""The dimensionless quadrature operators."""

from __future__ import annotations

from qiskit_bosons.operators.boson_action import ann, cre

from .. import BosonOperator


def position(mode: int) -> BosonOperator:
    r"""Builds the dimensionless position (quadrature) operator of a single mode.

    .. math::

        \hat{q} = \frac{1}{\sqrt{2}}\left(b^\dagger + b\right) \, .

    .. doctest::

        >>> from qiskit_bosons.operators.library import position
        >>> position(0).is_hermitian()
        True

    Args:
        mode: the mode index.

    Returns:
        The position operator as a :class:`~qiskit_bosons.operators.BosonOperator`.
    """
    scale = 2.0**-0.5
    return BosonOperator.from_dict({(cre(mode),): scale, (ann(mode),): scale})


def momentum(mode: int) -> BosonOperator:
    r"""Builds the dimensionless momentum (quadrature) operator of a single mode.

    .. math::

        \hat{p} = \frac{i}{\sqrt{2}}\left(b^\dagger - b\right) \, .

    .. doctest::

        >>> from qiskit_bosons.operators.library import momentum
        >>> momentum(0).is_hermitian()
        True

    Args:
        mode: the mode index.

    Returns:
        The momentum operator as a :class:`~qiskit_bosons.operators.BosonOperator`.
    """
    scale = 2.0**-0.5
    return BosonOperator.from_dict({(cre(mode),): 1j * scale, (ann(mode),): -1j * scale})
