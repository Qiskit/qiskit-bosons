// This code is a Qiskit project.
//
// (C) Copyright IBM 2026.
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at https://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use crate::operators::boson_operator::PyBosonOperator;
use pyo3::prelude::*;
use pyo3_stub_gen::derive::*;
use qiskit_bosons_core::operators::terms::filtering::diagonal::filter_diagonal_terms;

/// Filters out the terms of an operator that are diagonal in the occupation-number basis.
///
/// This removes every term from the provided
/// :class:`~qiskit_bosons.operators.BosonOperator` that is a product of number operators
/// (:math:`b^\dagger_i b_i`). This includes the constant term (a product of zero number operators),
/// single number operators, as well as higher-order products such as
/// :math:`n_i n_j = b^\dagger_i b^\dagger_j b_j b_i`. A mode may repeat, so products such as
/// :math:`b^\dagger_i b^\dagger_i b_i b_i = n_i (n_i - 1)` are diagonal, too.
///
/// .. hint::
///    A diagonal term only multiplies each Fock state by a phase. Filtering such terms is therefore
///    only appropriate when that phase is irrelevant to what you compute, for example when only
///    the sampled occupations of a time-evolved state matter. It is *not* a general simplification:
///    on a superposition of Fock states, the phases differ between its components.
///
/// .. caution::
///    The provided operator *must* be normal-ordered! This is an underlying assumption of the
///    implementation that is *not* being verified. See
///    :meth:`~qiskit_bosons.operators.BosonOperator.normal_ordered` for how to get an operator
///    of that form.
///
/// .. note::
///    The operator is modified *in place*. If it tracks group indices (see
///    :attr:`~qiskit_bosons.operators.BosonOperator.groups`), surviving terms retain their
///    relative grouping but the group indices are reassigned to a contiguous range starting from 0.
///    You must therefore *not* rely on the specific group index of any term being preserved across
///    a call to this function.
///
/// .. doctest::
///
///     >>> from qiskit_bosons.operators import BosonOperator
///     >>> from qiskit_bosons.operators.terms.filtering import filter_diagonal_terms
///     >>> op = BosonOperator.from_dict(
///     ...     {
///     ...         (): 1.0,  # constant (diagonal)
///     ...         ((True, 0), (False, 0)): 2.0,  # n_0 (diagonal)
///     ...         ((True, 0), (True, 1), (False, 1), (False, 0)): 3.0,  # n_0 n_1 (diagonal)
///     ...         ((True, 0), (False, 1)): 4.0,  # b†_0 b_1 (off-diagonal)
///     ...     }
///     ... )
///     >>> filter_diagonal_terms(op)
///     >>> list(op.iter_terms())
///     [([(True, 0), (False, 1)], (4+0j))]
///
/// Args:
///     op: the normal-ordered operator whose diagonal terms to filter out.
#[gen_stub_pyfunction(module = "qiskit_bosons._lib.operators.operators_terms.filtering.diagonal")]
#[pyfunction(name = "filter_diagonal_terms")]
pub fn py_filter_diagonal_terms(op: &Bound<PyBosonOperator>) {
    filter_diagonal_terms(&mut op.borrow_mut().inner);
}

#[pymodule]
pub mod diagonal {
    #[pymodule_export]
    use super::py_filter_diagonal_terms;
}
