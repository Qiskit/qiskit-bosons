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

use std::collections::HashSet;

use num_complex::Complex64;
use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::{class::basic::CompareOp, exceptions::PyNotImplementedError};
use pyo3_stub_gen::derive::*;
use std::collections::HashMap;

use qiskit_bosons_core::operators::boson_operator::BosonOperator;
use qiskit_bosons_core::operators::{OperatorMacro, OperatorTrait};

pub type PyBosonAction = (bool, u32);

crate::declare_operator_iters!(
    BosonOperatorDataIter,
    BosonOperatorDataGroupIter,
    "qiskit_bosons.operators.boson_operator",
    "qiskit_bosons._lib.operators.boson_operator",
    "BosonOperatorDataIter",
    "BosonOperatorDataGroupIter",
    PyBosonAction,
    "tuple[list[tuple[bool, int]], complex]",
    "tuple[list[tuple[bool, int]], complex, int]"
);

/// A bosonic operator.
///
/// Definition
/// ==========
///
/// This operator is defined by a linear combination of products of bosonic creation and
/// annihilation operators. The individual factors fulfill the canonical **commutation**
/// relations: [1]_
///
/// .. math::
///
///     \left[b^\dagger_i, b^\dagger_j\right] =
///     \left[b_i, b_j\right] = 0,~~\text{and}~~
///     \left[b_i, b^\dagger_j\right] = \delta_{ij} \, .
///
/// This makes the definition of the entire operator the following:
///
/// .. math::
///
///    \text{\texttt{BosonOperator}} = \sum_i c_i \bigotimes_j \hat{B_j} \, ,
///
/// where :math:`\hat{B_j} \in \{ b_j, b^\dagger_j \}` and :math:`c_i` is the (complex) coefficient
/// making up the linear combination of products. The index :math:`j` can take any value between 0
/// and the number of bosonic modes acted upon by the operator minus 1.
///
/// .. rubric:: Algebra
///
/// Two consequences of the commutation relations above shape this entire class:
///
/// * **A mode index may repeat.** :math:`b^\dagger_i b^\dagger_i \neq 0`: a mode can hold
///   arbitrarily many bosons, so a mode index may legitimately repeat within a single term. Terms
///   with repeated indices are never silently dropped.
/// * **Reordering never introduces a sign.** Swapping two factors is free unless they are a
///   creation/annihilation pair on the *same* mode, in which case
///   :math:`b_i b^\dagger_i = b^\dagger_i b_i + 1`. See :meth:`.normal_ordered`.
///
/// .. rubric:: Truncation
///
/// A bosonic mode has an **infinite-dimensional** Fock space. This operator is exact and carries no
/// occupation cutoff; truncation is introduced only when mapping to qubits, where it is a required,
/// explicit parameter. See :mod:`~qiskit_bosons.mappers` and :meth:`.max_occupation`.
///
/// .. _BosonOperator-implementation:
///
/// Implementation
/// ==============
///
/// This class stores the terms and coefficients in multiple sparse vectors, akin to the
/// `compressed sparse row format
/// <https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)>`_
/// commonly used for sparse matrices. More concretely, a single operator contains 4 arrays:
///
/// .. table::
///
///    ============== =================================================================================
///    ``coeffs``     A vector of complex coefficients consisting of two 64-bit floating point numbers.
///    ``actions``    A vector of booleans storing the nature of the second-quantization actions.
///    ``modes``      A vector of 32-bit integers storing the bosonic mode indices acted upon.
///    ``boundaries`` A vector of integers indicating the boundaries in ``actions`` and ``modes``.
///    ============== =================================================================================
///
/// Entries in ``actions`` indicate creation (annihilation) operators by ``True`` (``False``).
///
/// .. note::
///    The same mode index may appear multiple times within one term (e.g. ``+0 +0`` for
///    :math:`b^\dagger_0 b^\dagger_0`). Such repetitions are meaningful and are preserved verbatim.
///
/// .. note::
///    You may access **read-only copies** of these internal arrays via their respective methods:
///    :meth:`.get_coeffs`, :meth:`.get_actions`, :meth:`.get_modes`, and :meth:`.get_boundaries`.
///
/// This data structure allows for very efficient construction and manipulation of operators.
/// However, it implies that duplicate terms may be contained in an operator at any moment.
/// These must be resolved manually through the use of :meth:`.simplify`.
///
/// Construction
/// ------------
///
/// An operator can be constructed directly by providing the arrays outlined above:
///
/// .. doctest::
///
///     >>> from qiskit_bosons.operators import BosonOperator
///     >>> coeffs = [1.0, 2.0, -3.0, 4.0j, -0.5j]
///     >>> actions = [True, False, False, True, True, True, False, False]
///     >>> modes = [0, 0, 0, 1, 0, 1, 2, 3]
///     >>> boundaries = [0, 0, 1, 2, 4, 8]
///     >>> op = BosonOperator(coeffs, actions, modes, boundaries)
///     >>> print(format(op))
///       1.000000e0 +0.000000e0j * ()
///      -3.000000e0 +0.000000e0j * (-0)
///       0.000000e0 +4.000000e0j * (-0 +1)
///       2.000000e0 +0.000000e0j * (+0)
///      -0.000000e0-5.000000e-1j * (+0 +1 -2 -3)
///
/// For convenience, it is possible to construct an operator from a Python dictionary like so:
///
/// .. doctest::
///
///     >>> from qiskit_bosons.operators import cre, ann
///     >>> op = BosonOperator.from_dict(
///     ...     {
///     ...         (): 1.0,
///     ...         (cre(0),): 2.0,
///     ...         (ann(0),): -3.0,
///     ...         (ann(0), cre(1)): 4.0j,
///     ...         (cre(0), cre(1), ann(2), ann(3)): -0.5j,
///     ...     }
///     ... )
///     >>> print(format(op))
///       1.000000e0 +0.000000e0j * ()
///      -3.000000e0 +0.000000e0j * (-0)
///       0.000000e0 +4.000000e0j * (-0 +1)
///       2.000000e0 +0.000000e0j * (+0)
///      -0.000000e0-5.000000e-1j * (+0 +1 -2 -3)
///
/// In this example, we have leveraged :func:`.cre` and :func:`.ann` for creating the creation and
/// annihilation operators at the specified modes.
///
/// In addition, the following construction and quick helper methods are available:
///
/// .. autosummary::
///
///    zero
///    one
///    from_terms
///    from_terms_with_groups
///
/// Formatting
/// ----------
///
/// In the examples above, the constructed operators have been printed using the output from
/// :py:func:`format`, which results in a human-readable form of the operator.
///
/// .. doctest::
///
///     >>> print(format(op))
///       1.000000e0 +0.000000e0j * ()
///      -3.000000e0 +0.000000e0j * (-0)
///       0.000000e0 +4.000000e0j * (-0 +1)
///       2.000000e0 +0.000000e0j * (+0)
///      -0.000000e0-5.000000e-1j * (+0 +1 -2 -3)
///
/// .. note::
///    The printing order of ``format(op)`` gets explicitly sorted before printing. As such, it
///    does not reflect the order of the terms inside the operator.
///
/// An alternative form can be obtained from the :py:func:`repr` function, which results in a
/// Python-interpretable representation. In other words, this output can readily be copied and
/// pasted into a Python shell:
///
/// .. doctest::
///
///     >>> print(repr(op))
///     BosonOperator.from_dict({...})
///
/// Finally, for large operators both of these outputs may be very long and undesirable. Then, a
/// very simple form with minimal information can be obtained from the :py:func:`str` function:
///
/// .. doctest::
///
///     >>> print(str(op))
///     <BosonOperator with 5 terms>
///
/// Iteration
/// ---------
///
/// Since the underlying data structure is implemented in Rust and has a non-trivial layout, it
/// cannot be iterated over directly:
///
/// .. doctest::
///
///     >>> list(iter(op))
///     Traceback (most recent call last):
///       ...
///     TypeError: 'qiskit_bosons.operators.boson_operator.BosonOperator' object is not iterable
///
/// Instead, this class provides custom iterators to fulfill this purpose:
///
/// .. doctest::
///
///     >>> list(sorted(op.iter_terms()))
///     [([], (1+0j)), ([(False, 0)], (-3+0j)), ([(False, 0), (True, 1)], 4j), ([(True, 0)], (2+0j)), ([(True, 0), (True, 1), (False, 2), (False, 3)], (-0-0.5j))]
///
/// See also:
///     :meth:`iter_terms`
///         For more relevant implementation details.
///
/// The table below lists all available iterators:
///
/// .. autosummary::
///
///    iter_terms
///    iter_terms_with_groups
///
/// Arithmetics
/// -----------
///
/// The following arithmetic operations are supported:
///
/// Addition/Subtraction
/// ^^^^^^^^^^^^^^^^^^^^
///
/// .. doctest::
///
///     >>> op = BosonOperator.one()
///     >>> (op + op).simplify()
///     BosonOperator.from_dict({(): 2+0j})
///     >>> (op - op).simplify()
///     BosonOperator.from_dict({})
///     >>> op += op
///     >>> op.simplify()
///     BosonOperator.from_dict({(): 2+0j})
///     >>> op -= op
///     >>> op.simplify()
///     BosonOperator.from_dict({})
///
/// Scalar Multiplication/Divison
/// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
///
/// .. doctest::
///
///     >>> op = BosonOperator.one()
///     >>> (2 * op).simplify()
///     BosonOperator.from_dict({(): 2+0j})
///     >>> (op / 2).simplify()
///     BosonOperator.from_dict({(): 0.5+0j})
///     >>> op *= 2
///     >>> op.simplify()
///     BosonOperator.from_dict({(): 2+0j})
///     >>> op /= 2
///     >>> op.simplify()
///     BosonOperator.from_dict({(): 1+0j})
///
/// Operator Composition
/// ^^^^^^^^^^^^^^^^^^^^
///
/// .. note::
///    Operator composition corresponds to left-multiplication: ``c = a & b`` corresponds to
///    :math:`C = B A`. In other words, the composition of two operators returns a resulting
///    operator that performs "first ``a`` and then ``b``".
///
/// .. doctest::
///
///     >>> op1 = BosonOperator.from_dict({(): 2.0, (cre(0),): 3.0})
///     >>> op2 = BosonOperator.from_dict({(): 1.5, (ann(1),): 4.0})
///     >>> comp = (op1 & op2).simplify()
///     >>> print(format(comp))
///       3.000000e0 +0.000000e0j * ()
///       8.000000e0 +0.000000e0j * (-1)
///       1.200000e1 +0.000000e0j * (-1 +0)
///       4.500000e0 +0.000000e0j * (+0)
///     >>> op2 &= op1
///     >>> print(format(op2.simplify()))
///       3.000000e0 +0.000000e0j * ()
///       8.000000e0 +0.000000e0j * (-1)
///       4.500000e0 +0.000000e0j * (+0)
///       1.200000e1 +0.000000e0j * (+0 -1)
///     >>> squared = (op1 ** 2).simplify()
///     >>> print(format(squared))
///       4.000000e0 +0.000000e0j * ()
///       1.200000e1 +0.000000e0j * (+0)
///       9.000000e0 +0.000000e0j * (+0 +0)
///
/// .. note::
///    For convenience, the right-multiplication is implemented by ``c = a @ b`` (resulting in
///    :math:`C = A B`).
///
/// .. doctest::
///
///     >>> (op1 @ op2).equiv(op2 & op1)
///     True
///
/// Other Operations
/// ^^^^^^^^^^^^^^^^
///
/// In addition to the magic methods that correspond to the arithmetic operations outlined above,
/// the following methods are available:
///
/// .. autosummary::
///
///    adjoint
///    ichop
///    simplify
///    normal_ordered
///    relabel_modes
///
/// Properties
/// ^^^^^^^^^^
///
/// Finally, various methods exist to check certain properties of an operator:
///
/// .. autosummary::
///
///    is_hermitian
///    get_support
///    min_num_modes
///    max_rank
///    max_occupation
///    mode_occupations
///    conserves_particle_number
///    conserves_sector
///
/// .. [1] https://en.wikipedia.org/wiki/Second_quantization#Boson_creation_and_annihilation_operators
#[gen_stub_pyclass]
#[gen_stub(module = "qiskit_bosons._lib.operators.boson_operator")]
#[pyclass(
    module = "qiskit_bosons.operators.boson_operator",
    name = "BosonOperator"
)]
#[derive(Clone)]
pub struct PyBosonOperator {
    pub inner: BosonOperator,
}

impl From<BosonOperator> for PyBosonOperator {
    fn from(inner: BosonOperator) -> Self {
        Self { inner }
    }
}

crate::impl_operator_magic_methods!(PyBosonOperator, "BosonOperator");

#[gen_stub_pymethods]
#[pymethods]
impl PyBosonOperator {
    #[new]
    fn new(
        coeffs: Vec<Complex64>,
        actions: Vec<bool>,
        modes: Vec<u32>,
        boundaries: Vec<usize>,
    ) -> Self {
        Self {
            inner: BosonOperator {
                coeffs,
                actions,
                modes,
                boundaries,
                groups: None,
            },
        }
    }

    /// Constructs a new operator from a dictionary.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict(
    ///     ...     {
    ///     ...         (): 1.0-1.0j,
    ///     ...         ((True, 0), (False, 1)): 2.0,
    ///     ...     }
    ///     ... )
    ///     >>> print(format(op))
    ///       1.000000e0 -1.000000e0j * ()
    ///       2.000000e0 +0.000000e0j * (+0 -1)
    ///
    /// .. note::
    ///    The order of the terms inside the resulting operator is not guaranteed, and it can differ
    ///    between runs. The printing order of ``format(op)`` is sorted explicitly, so it does not
    ///    reveal it. Use :func:`.canonical_order` if you need a fixed term order.
    ///
    /// Args:
    ///     data: a dictionary mapping tuples of terms to complex coefficients. Each key is a tuple
    ///         of ``(bool, int)`` pairs. You may use :func:`.cre` and :func:`.ann` to simplify
    ///         their construction.
    ///
    /// Returns:
    ///     A new operator.
    #[classmethod]
    fn from_dict(_cls: &Bound<'_, PyType>, data: HashMap<Vec<(bool, u32)>, Complex64>) -> Self {
        let mut coeffs = vec![];
        let mut actions = vec![];
        let mut modes = vec![];
        let mut boundaries = vec![0];

        data.iter().for_each(|(terms, coeff)| {
            coeffs.push(*coeff);
            terms.iter().for_each(|(action, idx)| {
                actions.push(*action);
                modes.push(*idx);
            });
            boundaries.push(modes.len());
        });

        Self {
            inner: BosonOperator {
                coeffs,
                actions,
                modes,
                boundaries,
                groups: None,
            },
        }
    }

    /// Returns a read-only list of the operator's coefficients.
    ///
    /// .. note::
    ///    This method returns a **copy** of the internal data.
    ///
    /// .. seealso::
    ///    The explanation of the internal data structure,
    ///    :ref:`here <BosonOperator-implementation>`.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.one()
    ///     >>> op += -1j * BosonOperator.one()
    ///     >>> op.get_coeffs()
    ///     [(1+0j), -1j]
    ///
    /// Returns:
    ///     A list of the operator's coefficients.
    fn get_coeffs(&self) -> Vec<Complex64> {
        self.inner.coeffs().to_vec()
    }

    /// Returns a read-only list of the operator's actions.
    ///
    /// .. note::
    ///    This method returns a **copy** of the internal data.
    ///
    /// .. seealso::
    ///    The explanation of the internal data structure,
    ///    :ref:`here <BosonOperator-implementation>`.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.one()
    ///     >>> op += BosonOperator.from_dict({((True, 0), (False, 1)): 1.0})
    ///     >>> op.get_actions()
    ///     [True, False]
    ///
    /// Returns:
    ///     A list of the operator's actions.
    fn get_actions(&self) -> Vec<bool> {
        self.inner.actions().to_vec()
    }

    /// Returns a read-only list of the operator's acted-upon mode indices.
    ///
    /// .. note::
    ///    This method returns a **copy** of the internal data.
    ///
    /// .. seealso::
    ///    The explanation of the internal data structure,
    ///    :ref:`here <BosonOperator-implementation>`.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.one()
    ///     >>> op += BosonOperator.from_dict({((True, 0), (False, 1)): 1.0})
    ///     >>> op.get_modes()
    ///     [0, 1]
    ///
    /// Returns:
    ///     A list of the operator's modes.
    fn get_modes(&self) -> Vec<u32> {
        self.inner.modes().to_vec()
    }

    /// Returns a read-only list of the indices indicating the boundaries between operator terms.
    ///
    /// .. note::
    ///    This method returns a **copy** of the internal data.
    ///
    /// .. seealso::
    ///    The explanation of the internal data structure,
    ///    :ref:`here <BosonOperator-implementation>`.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.one()
    ///     >>> op += BosonOperator.from_dict({((True, 0), (False, 1)): 1.0})
    ///     >>> op.get_boundaries()
    ///     [0, 0, 2]
    ///
    /// Returns:
    ///     A list of the operator's terms boundaries.
    fn get_boundaries(&self) -> Vec<usize> {
        self.inner.boundaries().to_vec()
    }

    /// Returns the set of mode indices which this operator acts upon.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict(
    ///     ...     {
    ///     ...         ((True, 0), (False, 4)): 1,
    ///     ...         ((True, 1), (True, 3), (False, 4), (False, 7)): 1,
    ///     ...     }
    ///     ... )
    ///     >>> assert op.get_support() == {0, 1, 3, 4, 7}
    ///
    /// Returns:
    ///     The set of mode indices which this operator acts upon.
    fn get_support(&self) -> HashSet<u32> {
        self.inner.get_support()
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp, _py: Python<'_>) -> PyResult<bool> {
        let eq = self.inner.coeffs == other.inner.coeffs
            && self.inner.actions == other.inner.actions
            && self.inner.modes == other.inner.modes
            && self.inner.boundaries == other.inner.boundaries;
        match op {
            CompareOp::Eq => Ok(eq),
            CompareOp::Ne => Ok(!eq),
            _ => Err(PyErr::new::<PyNotImplementedError, _>("")),
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        let mut items_str = Vec::new();
        for term in self.inner.iter() {
            let key_parts: Vec<String> = term
                .iter()
                .map(|(action, orb)| {
                    format!("({}, {})", if *action { "True" } else { "False" }, orb)
                })
                .collect();
            let key_str = format!("({})", key_parts.join(", "));
            let val_str = format!("{}{:+}j", term.coeff.re, term.coeff.im);
            items_str.push(format!("{key_str}: {val_str}"));
        }
        Ok(format!(
            "BosonOperator.from_dict({{{}}})",
            items_str.join(", ")
        ))
    }

    fn __format__(&self, _format_spec: &str) -> PyResult<String> {
        let mut sorted: Vec<_> = self.inner.iter().collect();
        sorted.sort_by_key(|&term| term.into_vec());
        let mut items_str = Vec::new();
        for term in sorted {
            let key_parts: Vec<String> = term
                .iter()
                .map(|(action, orb)| format!("{}{}", if *action { "+" } else { "-" }, orb))
                .collect();
            let key_str = format!("({})", key_parts.join(" "));
            let val_str = format!("{:12.6e}{:+12.6e}j", term.coeff.re, term.coeff.im);
            items_str.push(format!("{val_str} * {key_str}"));
        }
        Ok(items_str.join("\n").to_string())
    }

    /// Constructs the additive identity operator.
    ///
    /// Adding the operator that is constructed by this method to another one has no effect.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({(): 2.0})
    ///     >>> zero = BosonOperator.zero()
    ///     >>> op + zero == op
    ///     True
    ///
    /// ..
    #[classmethod]
    fn zero(_cls: &Bound<'_, PyType>) -> Self {
        BosonOperator::zero().into()
    }

    /// Constructs the multiplicative identity operator.
    ///
    /// Composing the operator that is constructed by this method with another one has no effect.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({(): 2.0})
    ///     >>> one = BosonOperator.one()
    ///     >>> op & one == op
    ///     True
    ///
    /// ..
    #[classmethod]
    fn one(_cls: &Bound<'_, PyType>) -> Self {
        BosonOperator::one().into()
    }

    /// Returns an equivalent but simplified operator.
    ///
    /// The simplification process first sums all coefficients that belong to equal terms and then
    /// only retains those whose total coefficient exceeds the specified tolerance (just like
    /// :meth:`.ichop`).
    ///
    /// When an operator has been arithmetically manipulated or constructed in a way that does not
    /// guarantee unique terms, this method should be called before applying any method that
    /// filters numerically small coefficients to avoid loss of information. See the example below
    /// which showcases how :meth:`.ichop` can truncate terms that sum to a total coefficient
    /// magnitude which should not be truncated:
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> coeffs = [1e-5] * int(1e5)
    ///     >>> boundaries = [0] + [0] * int(1e5)
    ///     >>> op = BosonOperator(coeffs, [], [], boundaries)
    ///     >>> canon = op.simplify(1e-4)
    ///     >>> assert canon.equiv(op.one(), 1e-6)
    ///     >>> op.ichop(1e-4)
    ///     >>> assert op.equiv(op.zero(), 1e-6)
    ///
    /// Args:
    ///     atol: the absolute tolerance for the cutoff. This value defaults to ``1e-8``.
    ///
    /// Returns:
    ///     An equivalent but simplified operator.
    #[pyo3(signature = (atol=1e-8))]
    fn simplify(&self, atol: f64) -> Self {
        self.inner.simplify(atol).into()
    }

    /// Removes terms whose coefficient magnitude lies below the provided threshold.
    ///
    /// This method modifies the operator *in place* and returns ``None``.
    ///
    /// .. caution::
    ///    This method truncates coefficients greedily! If the acted upon operator may contain
    ///    separate coefficients for duplicate terms consider calling :meth:`.simplify` instead!
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({(): 1e-4, ((True, 0),): 1e-6, ((False, 0),): 1e-10})
    ///     >>> print(format(op))
    ///       1.000000e-4 +0.000000e0j * ()
    ///      1.000000e-10 +0.000000e0j * (-0)
    ///       1.000000e-6 +0.000000e0j * (+0)
    ///     >>> op.ichop()
    ///     >>> print(format(op))
    ///       1.000000e-4 +0.000000e0j * ()
    ///       1.000000e-6 +0.000000e0j * (+0)
    ///     >>> op.ichop(1e-5)
    ///     >>> print(format(op))
    ///       1.000000e-4 +0.000000e0j * ()
    ///
    /// Args:
    ///     atol: the absolute tolerance for the cutoff. This value defaults to ``1e-8``.
    #[pyo3(signature = (atol=1e-8))]
    fn ichop(&mut self, atol: f64) {
        self.inner.ichop(atol);
    }

    /// An iterator over the operator's terms.
    ///
    /// .. warning::
    ///    Mutating the iteration items does **not** affect the underlying operator data.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({(): 2.0, ((True, 0),): 1.0, ((False, 1),): -1.0j})
    ///     >>> list(sorted(op.iter_terms()))
    ///     [([], (2+0j)), ([(False, 1)], (-0-1j)), ([(True, 0)], (1+0j))]
    ///
    /// ..
    fn iter_terms(slf: PyRef<'_, Self>) -> PyResult<Py<BosonOperatorDataIter>> {
        let vectorized: Vec<(Vec<PyBosonAction>, Complex64)> = slf
            .inner
            .iter()
            .map(|term| (term.into_vec(), term.coeff))
            .collect();
        let iter = BosonOperatorDataIter {
            inner: vectorized.into_iter(),
        };
        Py::new(slf.py(), iter)
    }

    /// Constructs a new operator from an iterator of terms (see also :meth:`.iter_terms`).
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({(): 2.0, ((True, 0),): 1.0, ((False, 1),): -1.0j})
    ///     >>> op.equiv(BosonOperator.from_terms(op.iter_terms()))
    ///     True
    ///
    /// Args:
    ///     terms: an iterator of terms as produced by :meth:`.iter_terms`.
    ///
    /// Returns:
    ///     A new operator.
    #[classmethod]
    fn from_terms(_cls: &Bound<'_, PyType>, terms: &Bound<'_, PyAny>) -> PyResult<Self> {
        let mut inner = BosonOperator::zero();
        // We build the operator term-by-term via `_append_term` rather than routing through the
        // core `OperatorTrait::from_terms`. The latter takes an iterator of borrowed term *views*
        // (`TermView<'a>`) all tied to a single lifetime `'a`, whereas here each term arrives as an
        // owned Python tuple. Adapting to `from_terms` would require materializing every term into
        // owned buffers up front and keeping them all alive for the duration of the iterator (so the
        // views can borrow from them) — a full extra copy of the operator's data. Appending directly
        // keeps construction streaming, one term at a time, with no such buffer.
        terms.try_iter()?.try_for_each(|item| -> PyResult<()> {
            let (term, coeff) = item?.extract::<(Vec<PyBosonAction>, Complex64)>()?;
            let (actions, modes): (Vec<bool>, Vec<u32>) = term.into_iter().unzip();
            inner._append_term(coeff, &actions, &modes);
            Ok(())
        })?;
        Ok(inner.into())
    }

    /// An iterator over the operator's terms with their associated group index.
    ///
    /// .. warning::
    ///    Mutating the iteration items does **not** affect the underlying operator data.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator(
    ///     ...     [2.0, 1.0, -1.0],
    ///     ...     [True, False, True, False],
    ///     ...     [0, 1, 1, 0],
    ///     ...     [0, 0, 2, 4],
    ///     ... )
    ///     >>> op.groups = [0, 1, 1]
    ///     >>> list(sorted(op.iter_terms_with_groups()))
    ///     [([], (2+0j), 0), ([(True, 0), (False, 1)], (1+0j), 1), ([(True, 1), (False, 0)], (-1+0j), 1)]
    ///
    /// ..
    fn iter_terms_with_groups(slf: PyRef<'_, Self>) -> PyResult<Py<BosonOperatorDataGroupIter>> {
        let vectorized: Vec<(Vec<PyBosonAction>, Complex64, u32)> = slf
            .inner
            .iter_with_groups()
            .map(|term| (term.into_vec(), term.coeff, term.group))
            .collect();
        let iter = BosonOperatorDataGroupIter {
            inner: vectorized.into_iter(),
        };
        Py::new(slf.py(), iter)
    }

    /// Constructs a new operator from an iterator of terms with groups (see also
    /// :meth:`.iter_terms_with_groups`).
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator(
    ///     ...     [2.0, 1.0, -1.0],
    ///     ...     [True, False, True, False],
    ///     ...     [0, 1, 1, 0],
    ///     ...     [0, 0, 2, 4],
    ///     ... )
    ///     >>> op.groups = [0, 1, 1]
    ///     >>> reconstructed = BosonOperator.from_terms_with_groups(op.iter_terms_with_groups())
    ///     >>> op.equiv(reconstructed) and op.groups == reconstructed.groups
    ///     True
    ///
    /// Args:
    ///     terms: an iterator of terms as produced by :meth:`.iter_terms_with_groups`.
    ///
    /// Returns:
    ///     A new operator.
    #[classmethod]
    fn from_terms_with_groups(
        _cls: &Bound<'_, PyType>,
        terms: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let mut inner = BosonOperator::zero();
        let mut groups = vec![];
        // See `from_terms` for why we append directly instead of using the core `from_terms*`.
        terms.try_iter()?.try_for_each(|item| -> PyResult<()> {
            let (term, coeff, group) = item?.extract::<(Vec<PyBosonAction>, Complex64, u32)>()?;
            let (actions, modes): (Vec<bool>, Vec<u32>) = term.into_iter().unzip();
            inner._append_term(coeff, &actions, &modes);
            groups.push(group);
            Ok(())
        })?;
        // Built in lockstep with the terms above, so the length check of `set_groups` is redundant.
        inner._write_groups(Some(groups));
        Ok(inner.into())
    }

    /// An optional vector of `group indices` for each term.
    ///
    /// For more information refer to the :ref:`grouping_explanation` guide.
    #[getter]
    pub fn get_groups(&self) -> Option<Vec<u32>> {
        self.inner.groups.clone()
    }

    /// Sets the :attr:`groups` attribute.
    ///
    /// Raises:
    ///     ValueError: if the number of group indices differs from the number of terms. Assigning
    ///         ``None`` to clear the group indices is always allowed. This is validated here because
    ///         nothing downstream re-checks it: too few indices would silently drop the trailing
    ///         terms wherever terms are iterated together with their groups, and too many would make
    ///         :meth:`num_groups` report groups that no term carries.
    #[setter]
    pub fn set_groups(&mut self, groups: Option<Vec<u32>>) -> PyResult<()> {
        self.inner.set_groups(groups).map_err(crate::value_err)
    }

    /// Returns whether this operator tracks group indices.
    ///
    /// This is equivalent to (but cheaper than) checking ``op.groups is not None``, because it does
    /// not copy the group indices out of the operator in order to inspect them.
    ///
    /// .. note::
    ///    This returns ``True`` even when :attr:`groups` is an empty list, which is the state of a
    ///    grouped operator that holds no terms.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator(
    ///     ...     [1.0, 2.0, -1.0, -2.0],
    ///     ...     [True, False, True, False, True, False, True, False],
    ///     ...     [0, 1, 2, 3, 1, 0, 3, 2],
    ///     ...     [0, 2, 4, 6, 8],
    ///     ... )
    ///     >>> op.has_groups()
    ///     False
    ///     >>> op.groups = [0, 1, 0, 1]
    ///     >>> op.has_groups()
    ///     True
    ///
    /// Returns:
    ///     Whether :attr:`groups` is set on this operator.
    pub fn has_groups(&self) -> bool {
        self.inner.has_groups()
    }

    /// Returns the number of groups.
    ///
    /// If :attr:`groups` is ``None``, this function also returns ``None``. Otherwise, it will
    /// return the number of groups which is defined to be the largest occurring group index plus
    /// 1 (which may therefore be used as the index for the next group).
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator(
    ///     ...     [1.0, 2.0, -1.0, -2.0],
    ///     ...     [True, False, True, False, True, False, True, False],
    ///     ...     [0, 1, 2, 3, 1, 0, 3, 2],
    ///     ...     [0, 2, 4, 6, 8],
    ///     ... )
    ///     >>> op.groups = [0, 1, 0, 1]
    ///     >>> op.num_groups()
    ///     2
    ///
    /// Returns:
    ///     The largest group index in :attr:`groups` plus 1.
    pub fn num_groups(&self) -> Option<u32> {
        self.inner.num_groups()
    }

    /// Splits this operator into an optional list of new operators based on :attr:`groups`.
    ///
    /// If :attr:`groups` is ``None``, this function also returns ``None``. Otherwise, if
    /// ``group_indices`` is ``None`` (the default), it returns a list of one new operator for
    /// every group index in :attr:`groups`, in index order. If ``group_indices`` is given, only
    /// the requested indices are built, in the given order: this avoids the cost of constructing
    /// operators for groups that are never used, which is especially beneficial when only a small
    /// number of groups out of a much larger total are needed, e.g. when subsampling groups for a
    /// randomized product formula. A duplicate index in ``group_indices`` is returned once per
    /// occurrence.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator(
    ///     ...     [1.0, 2.0, -1.0, -2.0],
    ///     ...     [True, False, True, False, True, False, True, False],
    ///     ...     [0, 1, 2, 3, 1, 0, 3, 2],
    ///     ...     [0, 2, 4, 6, 8],
    ///     ... )
    ///     >>> print(op.split_out_groups())
    ///     None
    ///     >>> op.groups = [0, 1, 0, 1]
    ///     >>> groups = op.split_out_groups()
    ///     >>> for g in groups:
    ///     ...     print(list(sorted(g.iter_terms())))
    ///     [([(True, 0), (False, 1)], (1+0j)), ([(True, 1), (False, 0)], (-1+0j))]
    ///     [([(True, 2), (False, 3)], (2+0j)), ([(True, 3), (False, 2)], (-2+0j))]
    ///     >>> groups = op.split_out_groups(group_indices=[1])
    ///     >>> for g in groups:
    ///     ...     print(list(sorted(g.iter_terms())))
    ///     [([(True, 2), (False, 3)], (2+0j)), ([(True, 3), (False, 2)], (-2+0j))]
    ///
    /// Args:
    ///     group_indices: the group indices for which to build operators, in the desired output
    ///         order. When omitted, every group is built, in index order.
    ///
    /// Returns:
    ///     An optional vector of one new operator for each requested group index.
    #[pyo3(signature = (group_indices=None))]
    fn split_out_groups(
        slf: PyRef<'_, Self>,
        group_indices: Option<Vec<u32>>,
    ) -> Option<Vec<Self>> {
        let groups = slf.inner.split_out_groups(group_indices.as_deref());
        match groups {
            None => None,
            Some(g) => {
                let mut out = Vec::with_capacity(g.len());
                g.into_iter().for_each(|group_op| out.push(group_op.into()));
                Some(out)
            }
        }
    }

    /// Returns the Hermitian conjugate (or adjoint) of this operator.
    ///
    /// This affects the terms and coefficients as follows:
    ///
    /// - the actions in each term reverse their order and flip between creation and annihilation
    /// - the coefficients are complex conjugated
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({(): -1.0j, ((True, 0), (False, 1)): 1.0})
    ///     >>> adj = op.adjoint()
    ///     >>> print(format(adj))
    ///      -0.000000e0 +1.000000e0j * ()
    ///       1.000000e0 -0.000000e0j * (+1 -0)
    ///
    /// ..
    fn adjoint(&self) -> Self {
        self.inner.adjoint().into()
    }

    /// Checks this operator for equivalence with another operator.
    ///
    /// Equivalence in this context means approximate equality up to the specified absolute
    /// tolerance. To be more precise, this method returns ``True``, when all the absolute values
    /// of the coefficients in the difference ``other - self`` are below the specified threshold
    /// ``atol``.
    ///
    /// .. note::
    ///    This is the mathematical comparison you almost always want. It differs from the ``==``
    ///    operator, which tests exact equality of the *stored* terms (their coefficients, actions,
    ///    modes, and internal term boundaries) with no tolerance and no simplification. Two
    ///    mathematically equal operators can therefore compare unequal under ``==`` if they are
    ///    stored differently -- for example an unsimplified ``a + a`` versus ``2 * a``, or terms
    ///    held in a different order. Use ``equiv`` to compare operators up to numerical tolerance.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({(): 1e-7})
    ///     >>> zero = BosonOperator.zero()
    ///     >>> op.equiv(zero)
    ///     False
    ///     >>> op.equiv(zero, 1e-6)
    ///     True
    ///     >>> op.equiv(zero, 1e-9)
    ///     False
    ///
    /// Args:
    ///     other: the other operator to compare with.
    ///     atol: the absolute tolerance for the comparison. This value defaults to ``1e-8``.
    #[pyo3(signature = (other, atol=1e-8))]
    fn equiv(&self, other: &Self, atol: f64) -> bool {
        self.inner.equiv(&other.inner, atol)
    }

    /// Returns an equivalent operator with normal ordered terms.
    ///
    /// The normal order of an operator term is defined such that all creation actions appear
    /// before all annihilation actions.
    /// Within each group, the acted-upon modes are ordered lexicographically. Whether their order
    /// is ascending or descending depends upon the value of the ``sandwich`` argument:
    ///
    /// - ``None`` (default): both groups are ordered lexicographically descending (e.g.
    ///   ``+1 +0 -1 -0``)
    /// - ``True``: larger indices appear towards the middle, i.e. creation actions are
    ///   lexicographically ascending while annihilation ones are descending (e.g.
    ///   ``+0 +1 -1 -0``)
    /// - ``False``: smaller indices appear towards the middle, i.e. creation actions are
    ///   lexicographically descending while annihilation ones are ascending (e.g.
    ///   ``+1 +0 -0 -1``)
    ///
    /// .. note::
    ///    When a term is being reordered, the canonical **commutation** relations have to be taken
    ///    into account, :math:`b_i b^\dagger_j = \delta_{ij} + b^\dagger_j b_i`, implying that the
    ///    number of terms may change. Note the sign: the contracted term is *added*. Reordering
    ///    itself never introduces a sign, and no term ever vanishes -- :math:`b^\dagger_i
    ///    b^\dagger_i` is a perfectly good operator.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({((False, 1), (True, 1), (False, 0), (True, 0)): 1})
    ///     >>> print(format(op.normal_ordered().simplify()))
    ///       1.000000e0 +0.000000e0j * ()
    ///       1.000000e0 +0.000000e0j * (+0 -0)
    ///       1.000000e0 +0.000000e0j * (+1 -1)
    ///       1.000000e0 +0.000000e0j * (+1 +0 -1 -0)
    ///     >>> print(format(op.normal_ordered(sandwich=True).simplify()))
    ///       1.000000e0 +0.000000e0j * ()
    ///       1.000000e0 +0.000000e0j * (+0 -0)
    ///       1.000000e0 +0.000000e0j * (+0 +1 -1 -0)
    ///       1.000000e0 +0.000000e0j * (+1 -1)
    ///     >>> print(format(op.normal_ordered(sandwich=False).simplify()))
    ///       1.000000e0 +0.000000e0j * ()
    ///       1.000000e0 +0.000000e0j * (+0 -0)
    ///       1.000000e0 +0.000000e0j * (+1 -1)
    ///       1.000000e0 +0.000000e0j * (+1 +0 -0 -1)
    ///
    /// Returns:
    ///     An equivalent but normal-ordered operator.
    #[pyo3(signature = (sandwich=None))]
    fn normal_ordered(&self, sandwich: Option<bool>) -> Self {
        self.inner.normal_ordered(sandwich).into()
    }

    /// Returns whether this operator is Hermitian.
    ///
    /// .. note::
    ///    This check is implemented using :meth:`.equiv` on the :meth:`.normal_ordered` difference
    ///    of ``self`` and its :meth:`.adjoint` and :meth:`.zero`.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({
    ///     ...     ((True, 0), (False, 1)): 1.00001j,
    ///     ...     ((True, 1), (False, 0)): -1j,
    ///     ... })
    ///     >>> op.is_hermitian()
    ///     False
    ///     >>> op.is_hermitian(1e-4)
    ///     True
    ///
    /// Args:
    ///     atol: The numerical accuracy upto which coefficients are considered equal. This value
    ///         defaults to ``1e-8``.
    ///
    /// Returns:
    ///     Whether this operator is Hermitian.
    #[pyo3(signature = (atol=1e-8))]
    fn is_hermitian(&self, atol: f64) -> bool {
        self.inner.is_hermitian(atol)
    }

    /// Returns the maximum rank of the terms in this operator.
    ///
    /// .. note::
    ///    The length of the longest term can depend on the operator's form which means that (for
    ///    example) operator simplification or normal-ordering can result in a different maximum
    ///    rank.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({
    ///     ...     ((True, 0), (False, 1), (True, 2), (False, 3)): 1,
    ///     ... })
    ///     >>> op.max_rank()
    ///     4
    ///
    /// Returns:
    ///     The maximum rank of this operator.
    fn max_rank(&self) -> u32 {
        self.inner.max_rank()
    }

    /// Returns the highest occupation any single mode is driven to by any one term.
    ///
    /// This is the largest number of creation operators acting on one and the same mode within a
    /// single term, which is the maximum of :meth:`mode_occupations` over all modes. It measures how
    /// far one application of a term can raise a mode, so it is a **lower bound** on the truncation
    /// level a boson-to-qubit encoding needs in order to represent this operator without clipping
    /// only when the state starts in the vacuum: applying :math:`b^\dagger_i b^\dagger_i` to the
    /// vacuum populates occupation 2, which a cutoff of ``max_occupation < 2`` cannot hold. For a
    /// state with :math:`n` bosons in a mode, the mode can reach :math:`n` plus this value.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator, ann, cre
    ///     >>> hopping = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
    ///     >>> hopping.max_occupation()
    ///     1
    ///     >>> interaction = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
    ///     >>> interaction.max_occupation()
    ///     2
    ///
    /// .. warning::
    ///    This is only a lower bound, **not** a sufficient cutoff. The occupation reachable by
    ///    *repeated* application of the operator -- as happens over a time evolution -- is
    ///    unbounded, and a mode rotation moves amplitude between occupation sectors, so acting on
    ///    :math:`\lvert 1, 1 \rangle` can populate :math:`\lvert 2, 0 \rangle`. Choosing the cutoff
    ///    remains a physical modeling decision; this only rules out cutoffs that are certainly
    ///    too small.
    ///
    /// Returns:
    ///     The largest per-mode creation count over all terms.
    fn max_occupation(&self) -> u32 {
        self.inner.max_occupation()
    }

    /// Returns, for each mode, the highest occupation it is driven to by any one term.
    ///
    /// The occupation of a mode in a term is the number of creation operators acting on that mode
    /// within the term. Annihilation operators are not subtracted, so the entry is how far one
    /// application of a term can **raise** the mode, whatever the term does to it afterwards.
    /// Because it measures a rise, it is a lower bound on the cutoff of the mode only when the mode
    /// starts empty. A mode that starts in occupation :math:`n` can reach :math:`n` plus the entry,
    /// so add the initial occupation of the state when you pick a cutoff. Like
    /// :meth:`max_occupation`, the result is not a sufficient cutoff.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator, ann, cre
    ///     >>> op = BosonOperator.from_dict({
    ///     ...     (cre(0), cre(0), ann(3)): 1.0,
    ///     ...     (cre(3), cre(3), cre(3)): 1.0,
    ///     ... })
    ///     >>> op.mode_occupations()
    ///     [2, 0, 0, 3]
    ///     >>> op.mode_occupations([3, 1])
    ///     [3, 0]
    ///
    /// Args:
    ///     modes: the modes to report, in the order to report them. If omitted (``None``), every
    ///         mode from ``0`` to the largest one the operator acts on is reported, with zero for
    ///         the modes it does not touch. A requested mode the operator does not touch,
    ///         including one beyond :meth:`min_num_modes`, reports zero.
    ///
    /// Returns:
    ///     One bound per reported mode.
    #[pyo3(signature = (modes=None))]
    fn mode_occupations(&self, modes: Option<Vec<u32>>) -> Vec<u32> {
        self.inner.mode_occupations(modes.as_deref())
    }

    /// Returns the smallest number of modes that can host this operator.
    ///
    /// This is one more than the largest mode index the operator acts on, or ``0`` if it acts on
    /// no mode.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator, ann, cre
    ///     >>> op = BosonOperator.from_dict({(cre(2), ann(5)): 1.0})
    ///     >>> op.min_num_modes()
    ///     6
    ///     >>> BosonOperator.zero().min_num_modes()
    ///     0
    ///
    /// Returns:
    ///     The minimum number of modes.
    fn min_num_modes(&self) -> u32 {
        self.inner.min_num_modes()
    }

    /// Returns whether this operator is particle-number conserving.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({((True, 0), (False, 1)): 1})
    ///     >>> op.conserves_particle_number()
    ///     True
    ///     >>> op = BosonOperator.from_dict({((True, 0),): 1})
    ///     >>> op.conserves_particle_number()
    ///     False
    ///
    /// Returns:
    ///     Whether this operator is particle-number conserving.
    fn conserves_particle_number(&self) -> bool {
        self.inner.conserves_particle_number()
    }

    /// Returns whether every term conserves particle number within each mode block.
    ///
    /// ``block_sizes`` partitions the mode range into consecutive, non-overlapping blocks: block
    /// ``b`` spans modes ``[start_b, start_b + block_sizes[b])`` where ``start_b`` is the sum of the
    /// preceding block sizes. A term conserves the sector if and only if, in *every* block, its
    /// number of creation operators equals its number of annihilation operators. A term acting on a
    /// mode beyond the final block does not conserve the sector.
    ///
    /// An empty ``block_sizes`` treats all modes as a single block, making this equivalent to
    /// :meth:`conserves_particle_number`.
    ///
    /// Blocks are useful for structured bosonic systems. A vibrational Hamiltonian in the
    /// Christiansen :math:`n`-mode form, for instance, assigns one block of `modals` to each physical
    /// vibrational mode, and a physical state carries exactly one excitation per block; checking
    /// per-block conservation verifies that an operator respects that structure.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({((True, 0), (False, 2)): 1})
    ///     >>> op.conserves_sector([4])  # one block spanning all 4 modes
    ///     True
    ///     >>> op.conserves_sector([2, 2])  # moves a boson from the first block to the second
    ///     False
    ///
    /// Args:
    ///     block_sizes: the sizes of the consecutive mode blocks that each must be conserved.
    ///
    /// Returns:
    ///     Whether every term conserves particle number within each mode block.
    fn conserves_sector(&self, block_sizes: Vec<u32>) -> bool {
        self.inner.conserves_sector(&block_sizes)
    }

    /// Returns a new operator with relabeled modes.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.operators import BosonOperator
    ///     >>> op = BosonOperator.from_dict({
    ///     ...     ((True, 0), (False, 1)): 1,
    ///     ...     ((True, 0), (False, 1), (True, 2), (False, 3)): 1,
    ///     ... })
    ///     >>> permutation = [5, 6, 4, 3]
    ///     >>> relabeled = op.relabel_modes(permutation)
    ///     >>> print(format(relabeled))
    ///       1.000000e0 +0.000000e0j * (+5 -6)
    ///       1.000000e0 +0.000000e0j * (+5 -6 +4 -3)
    ///
    /// Args:
    ///     permutation: the index permutation list. Mode ``i`` is relabeled to ``permutation[i]``,
    ///         so the list must contain no duplicate entries and must be long enough to index every
    ///         mode the operator acts upon (its length must exceed the operator's largest mode
    ///         index).
    ///
    /// Returns:
    ///     A new operator with its modes relabeled.
    ///
    /// Raises:
    ///     ValueError: if ``permutation`` contains duplicate entries, or is too short to relabel
    ///         some mode the operator acts upon.
    fn relabel_modes(&self, permutation: Vec<u32>) -> PyResult<Self> {
        self.inner
            .relabel_modes(permutation)
            .map(Into::into)
            .map_err(crate::value_err)
    }

    /// Returns the constructor arguments needed to pickle this operator.
    ///
    /// Together with :meth:`__getstate__`/:meth:`__setstate__` (which round-trip
    /// :attr:`groups`), this makes instances of this class picklable.
    fn __getnewargs__(&self) -> (Vec<Complex64>, Vec<bool>, Vec<u32>, Vec<usize>) {
        (
            self.inner.coeffs.clone(),
            self.inner.actions.clone(),
            self.inner.modes.clone(),
            self.inner.boundaries.clone(),
        )
    }

    /// Returns the pickled state of this operator (its :attr:`groups`).
    fn __getstate__(&self) -> Option<Vec<u32>> {
        self.inner.groups.clone()
    }

    /// Restores this operator's :attr:`groups` from its pickled state.
    ///
    /// Raises:
    ///     ValueError: if the pickled group indices do not number one per term. This is the one
    ///         route by which externally supplied group indices reach an operator without passing
    ///         through the :attr:`groups` setter, so it is validated rather than trusted: a payload
    ///         from an incompatible version, or a hand-edited one, would otherwise construct an
    ///         operator whose group indices nothing downstream re-checks.
    fn __setstate__(&mut self, state: Option<Vec<u32>>) -> PyResult<()> {
        self.inner.set_groups(state).map_err(crate::value_err)
    }
}

#[pymodule]
pub mod boson_operator {
    #[pymodule_export]
    use super::PyBosonOperator;
}
