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

"""Operator terms grouping.

No automatic grouping strategy ships today. Group indices are assigned directly via
:attr:`~qiskit_bosons.operators.BosonOperator.groups` -- which is how a
:class:`~qiskit_bosons.circuit.library.Evolution` gate is split into one sub-evolution per group.

What this module does provide is the *analysis* of an existing grouping. Group indices carry no
intrinsic meaning, so none of these functions reports whether a grouping is "correct"; each answers
one narrow question about it, so that an assumption a downstream consumer makes can be checked up
front rather than being paid for on every call.
"""

from qiskit_bosons._lib.operators.operators_terms.grouping.analysis import (
    group_coeff_means,
    groups_are_hermitian,
    groups_have_uniform_coeffs,
)

__all__ = [
    "group_coeff_means",
    "groups_are_hermitian",
    "groups_have_uniform_coeffs",
]
