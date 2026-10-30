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

use super::encoding::{EncodingArg, PyLevelEncoding};
use crate::operators::boson_operator::PyBosonOperator;
use pyo3::prelude::*;
use pyo3_stub_gen::derive::*;
use qiskit_bosons_core::mappers::levels::spec::{LevelSpec, ModeLayout};

/// The cutoffs of a spec as given from Python: one integer shared by every mode, or one per mode.
///
/// A sequence is never broadcast, whatever its length: only a bare integer is.
#[derive(FromPyObject)]
pub enum CutoffArg {
    Shared(u32),
    PerMode(Vec<u32>),
}

/// The truncation of a register of bosonic modes and the encoding that lays it out on qubits.
///
/// A spec says how many levels each mode keeps and how those levels are placed on qubits, which in
/// turn fixes how many qubits each mode occupies and where its qubits start. Everything that has to
/// agree about that (the mapper, the transpiler layout, the state preparation) is given the same
/// spec, so it cannot disagree.
///
/// The cutoff of a mode is its ``max_occupation``, the largest boson number it can hold, so the mode
/// keeps :math:`d = \texttt{max\_occupation} + 1` levels. There is deliberately **no default**:
/// truncation is a physical modeling decision, and the cutoff of every mode must bound the
/// occupation `reachable` in that mode, not merely the occupation of the initial state.
///
/// The encoding has no default either. No :class:`.LevelEncoding` is best in general, and the choice
/// decides how many qubits every mode occupies, so it is made explicitly.
///
/// A single integer applies to every mode and leaves the number of modes open, to be fixed by
/// whatever the spec is applied to. A sequence gives one cutoff per mode and so also fixes the
/// number of modes; a sequence is never broadcast, so one whose length differs from the number of
/// modes it is used with is an error, including a sequence of length one.
///
/// .. doctest::
///
///     >>> from qiskit_bosons.mappers import LevelSpec
///     >>> spec = LevelSpec([4, 2, 2], "gray")
///     >>> spec.widths
///     [3, 2, 2]
///     >>> spec.offsets
///     [0, 3, 5]
///     >>> spec.total_width
///     7
///     >>> LevelSpec(2, "unary_projected").is_uniform
///     True
///
/// .. note::
///    Under a compact encoding (:attr:`.LevelEncoding.STANDARD_BINARY` or :attr:`.LevelEncoding.GRAY`)
///    a mode with a single level (``max_occupation = 0``) occupies **zero** qubits, because its
///    only level needs no qubit to be told apart from any other. Such a mode is allowed and does
///    not move the qubits of the modes after it, but a spec in which every mode would occupy zero
///    qubits is rejected.
///
/// .. note::
///    The unary encodings hold at most 128 levels per mode, that is ``max_occupation`` up to 127,
///    because they identify a level by one bit of a 128-bit mask. A larger cutoff under
///    :attr:`.LevelEncoding.UNARY_PROJECTED` or :attr:`.LevelEncoding.UNARY_UNPROJECTED` raises a
///    :class:`ValueError`, which asks for a compact encoding instead. This is a fixed limit of the
///    implementation and cannot be changed at runtime. The compact encodings are not affected by it.
///
/// Args:
///     max_occupation: the cutoff of every mode, either one integer shared by all modes or one per
///         mode.
///     encoding: the level encoding, either a :class:`.LevelEncoding` or its string name.
///
/// Raises:
///     ValueError: if no mode would occupy a qubit, if a cutoff is too large (for a unary
///         encoding, if it exceeds 127), or if ``encoding`` is a string that is not a recognized
///         name.
///     TypeError: if ``max_occupation`` is neither an integer nor a sequence of integers, or if
///         ``encoding`` is neither a :class:`.LevelEncoding` nor a string.
#[gen_stub_pyclass]
#[gen_stub(module = "qiskit_bosons._lib.mappers.levels")]
#[pyclass(
    module = "qiskit_bosons.mappers.levels",
    name = "LevelSpec",
    frozen,
    eq,
    hash
)]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PyLevelSpec {
    pub inner: LevelSpec,
}

#[gen_stub_pymethods]
#[pymethods]
impl PyLevelSpec {
    #[new]
    #[pyo3(signature = (max_occupation, encoding))]
    fn new(
        #[gen_stub(override_type(type_repr = "int | collections.abc.Sequence[int]", imports = ("collections.abc")))]
        max_occupation: CutoffArg,
        // `LevelEncoding` lives in the same stub module, so naming it needs no import.
        #[gen_stub(override_type(type_repr = "LevelEncoding | str"))] encoding: EncodingArg,
    ) -> PyResult<Self> {
        let encoding = encoding.resolve()?;
        let inner = match max_occupation {
            CutoffArg::Shared(max_occupation) => LevelSpec::uniform(max_occupation, encoding),
            CutoffArg::PerMode(max_occupations) => LevelSpec::per_mode(max_occupations, encoding),
        }
        .map_err(crate::value_err)?;
        Ok(Self { inner })
    }

    /// The :class:`.LevelEncoding` of this spec.
    #[getter]
    fn encoding(&self) -> PyLevelEncoding {
        self.inner.encoding().into()
    }

    /// The shared cutoff of every mode if this spec was built from a single integer, and the
    /// per-mode list of cutoffs otherwise.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "int | list[int]"))]
    fn max_occupation(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match (
            self.inner.uniform_max_occupation(),
            self.inner.max_occupations(),
        ) {
            (Some(shared), _) => Ok(shared.into_pyobject(py)?.into_any().unbind()),
            (None, Some(per_mode)) => Ok(per_mode.to_vec().into_pyobject(py)?.into_any().unbind()),
            (None, None) => unreachable!("a spec holds either a shared or a per-mode cutoff"),
        }
    }

    /// Whether every mode shares one cutoff.
    ///
    /// A spec built from a single integer leaves the number of modes open, so its :attr:`widths`,
    /// :attr:`offsets` and :attr:`total_width` are only defined once :meth:`resolve` has fixed it.
    #[getter]
    fn is_uniform(&self) -> bool {
        self.inner.is_uniform()
    }

    /// The number of modes this spec covers, or ``None`` for a spec built from a single integer,
    /// which leaves the number of modes open.
    #[getter]
    fn num_modes(&self) -> Option<u32> {
        self.inner.num_modes()
    }

    /// The number of qubits each mode occupies.
    ///
    /// This is :meth:`.LevelEncoding.num_qubits` of each mode's level count. It is only defined for
    /// a spec that covers a definite number of modes.
    ///
    /// Raises:
    ///     ValueError: if the spec was built from a single integer, so that its number of modes is
    ///         open. Use :meth:`resolve` to fix it first.
    #[getter]
    fn widths(&self) -> PyResult<Vec<u32>> {
        Ok(self.fixed_layout()?.widths)
    }

    /// The index of the first qubit of each mode.
    ///
    /// This is the running sum of :attr:`widths`, so a mode that occupies no qubit does not move the
    /// next one. It is only defined for a spec that covers a definite number of modes.
    ///
    /// Raises:
    ///     ValueError: if the spec was built from a single integer. Use :meth:`resolve` first.
    #[getter]
    fn offsets(&self) -> PyResult<Vec<u32>> {
        Ok(self.fixed_layout()?.offsets)
    }

    /// The number of qubits all modes together occupy.
    ///
    /// Raises:
    ///     ValueError: if the spec was built from a single integer. Use :meth:`resolve` first.
    #[getter]
    fn total_width(&self) -> PyResult<u32> {
        Ok(self.fixed_layout()?.total_width)
    }

    /// Returns this spec applied to exactly ``num_modes`` modes.
    ///
    /// A spec built from a single integer is expanded to ``num_modes`` modes. A spec with one cutoff
    /// per mode must already cover exactly that many.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.mappers import LevelSpec
    ///     >>> LevelSpec(2, "unary_projected").resolve(3).widths
    ///     [3, 3, 3]
    ///
    /// Args:
    ///     num_modes: the number of modes.
    ///
    /// Returns:
    ///     A spec with one cutoff per mode.
    ///
    /// Raises:
    ///     ValueError: if the spec covers a different number of modes.
    fn resolve(&self, num_modes: u32) -> PyResult<Self> {
        let layout = self
            .inner
            .mode_layout(num_modes)
            .map_err(crate::value_err)?;
        let inner = LevelSpec::per_mode(layout.max_occupations, self.inner.encoding())
            .map_err(crate::value_err)?;
        Ok(Self { inner })
    }

    /// Returns by how much each mode's cutoff falls short of an operator's occupation bounds.
    ///
    /// For every mode this is the amount by which the cutoff is below the bound that
    /// :meth:`.BosonOperator.mode_occupations` reports for ``operator``, or zero if it is not. So an
    /// all-zero result means that no bound is violated, and a non-zero sum means that something is
    /// truncated.
    ///
    /// That bound is how far one application of a term can raise a mode, and it takes no account of
    /// the occupation the mode starts in. A mode that starts with :math:`n` bosons can reach
    /// :math:`n` plus the bound, so the cutoff of a mode needs to be at least that. For example, a
    /// beamsplitter generator raises each of its two modes by one, so this method reports no
    /// shortfall at a cutoff of one, although the beamsplitter acting on
    /// :math:`\lvert 1, 1 \rangle` populates :math:`\lvert 2, 0 \rangle` and needs a cutoff of two.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.mappers import LevelSpec
    ///     >>> from qiskit_bosons.operators import BosonOperator, ann, cre
    ///     >>> op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0, (cre(1), ann(1)): 1.0})
    ///     >>> LevelSpec([1, 3], "unary_projected").truncation(op)
    ///     [1, 0]
    ///     >>> sum(LevelSpec([2, 2], "unary_projected").truncation(op)) == 0
    ///     True
    ///
    /// .. warning::
    ///    This measures the lower-bound heuristic of :meth:`.BosonOperator.mode_occupations` and
    ///    nothing else. It does not say whether the mapping is accurate: the bound is not a
    ///    sufficient cutoff, so an all-zero result is no guarantee. Conversely a model can report a
    ///    shortfall although its mapping is exact. For example a hard-core model at
    ///    ``max_occupation = 1`` contains the term :math:`b^\dagger b^\dagger b b`, which reports a
    ///    bound of 2 but is exactly zero on the two retained levels.
    ///
    /// Args:
    ///     operator: the operator to check.
    ///
    /// Returns:
    ///     A list with one shortfall per mode. For a spec with one cutoff per mode, it covers every
    ///     mode of the spec. For a spec built from a single integer, it covers the modes up to the
    ///     largest one the operator acts on.
    ///
    /// Raises:
    ///     ValueError: if the operator acts on more modes than the spec covers.
    fn truncation(&self, operator: PyRef<'_, PyBosonOperator>) -> PyResult<Vec<u32>> {
        let bounds = operator.inner.mode_occupations(None);
        self.inner.shortfalls(&bounds).map_err(crate::value_err)
    }

    fn __repr__(&self) -> String {
        let max_occupation = match (
            self.inner.uniform_max_occupation(),
            self.inner.max_occupations(),
        ) {
            (Some(shared), _) => shared.to_string(),
            (None, Some(per_mode)) => format!("{per_mode:?}"),
            (None, None) => unreachable!("a spec holds either a shared or a per-mode cutoff"),
        };
        format!(
            "LevelSpec({max_occupation}, {})",
            self.encoding().__repr__()
        )
    }

    /// Returns the callable and arguments that rebuild this spec when unpickled.
    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, SpecArgs)> {
        let cls = py.get_type::<Self>().into_any();
        Ok((cls, (self.max_occupation(py)?, self.encoding().name())))
    }
}

/// The constructor arguments of a spec: its cutoffs and the string name of its encoding.
type SpecArgs = (Py<PyAny>, &'static str);

impl PyLevelSpec {
    fn fixed_layout(&self) -> PyResult<ModeLayout> {
        match self.inner.num_modes() {
            Some(num_modes) => self.inner.mode_layout(num_modes).map_err(crate::value_err),
            None => Err(crate::value_err(
                "this spec was built from a single cutoff, so its number of modes is open; \
                 call `resolve(num_modes)` first",
            )),
        }
    }
}
