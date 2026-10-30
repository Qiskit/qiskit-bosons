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

use pyo3::prelude::*;
use pyo3_stub_gen::derive::*;
use qiskit_bosons_core::mappers::levels::encoding::{LevelEncoding, MAX_UNARY_LEVELS};

/// The layout of the levels of a single system on qubits.
///
/// A level encoding decides which computational basis state of a block of qubits represents each of
/// the :math:`d` levels of one system. It is independent of what the levels are, so mappers in
/// different bases (for example :func:`~qiskit_bosons.mappers.library.occupation_encoding`, whose
/// levels are boson occupation numbers) share the same choice.
///
/// .. list-table::
///    :header-rows: 1
///
///    * - Encoding
///      - Qubits per system
///      - Level :math:`n` is stored as
///    * - :attr:`UNARY_PROJECTED`
///      - :math:`d`
///      - the one-hot pattern with qubit :math:`n` set, with a projector on every idle qubit
///    * - :attr:`UNARY_UNPROJECTED`
///      - :math:`d`
///      - the same one-hot pattern, without the projectors
///    * - :attr:`STANDARD_BINARY`
///      - :math:`\lceil \log_2 d \rceil`
///      - the binary digits of :math:`n`
///    * - :attr:`GRAY`
///      - :math:`\lceil \log_2 d \rceil`
///      - the Gray code :math:`n \oplus (n \gg 1)`
///
/// There is no universally best choice, so all of them are offered on an equal footing. [1]_ The
/// two unary encodings share their bit patterns and give ladder terms with only two non-diagonal
/// factors, at a cost of :math:`d` qubits. They differ in the operator built from those patterns.
/// :attr:`UNARY_UNPROJECTED` leaves the idle qubits alone, which is the weight-2 form whose gate
/// count grows only linearly in :math:`d`. :attr:`UNARY_PROJECTED` also puts a projector
/// :math:`\lvert 0 \rangle\langle 0 \rvert` on each idle qubit, so its terms touch all :math:`d`
/// qubits and are costlier to synthesize. The two are equal on the one-hot subspace, and the
/// unprojected operator is undefined outside it, so the choice is a modeling decision and not a
/// pure optimization. The two compact encodings need only :math:`\lceil \log_2 d \rceil` qubits but
/// produce terms with more non-diagonal factors, and their projectors cannot be dropped, because
/// there the zero bits are part of the level number.
///
/// Sawaya et al. find that the compact codes are frequently better in `both` qubit and gate count
/// -- at :math:`d = 4` the unary encoding needs roughly three times `more` two-qubit gates,
/// whereas by :math:`d = 10` it needs about half. They further find the Gray code preferable for
/// tridiagonal (hopping) terms and standard binary for diagonal (number-operator) terms.
///
/// Note that :attr:`STANDARD_BINARY` and :attr:`GRAY` emit the same `number` of terms; the Gray
/// code's advantage lies entirely in the weight `distribution`, since consecutive levels differ in
/// exactly one bit.
///
/// Every function that takes a level encoding also accepts its string name, which is what
/// :func:`str` returns: ``"unary_projected"``, ``"unary_unprojected"``, ``"standard_binary"`` or
/// ``"gray"``. :meth:`from_name` is the inverse, and rejects any other name, including the former
/// ``"unary"``:
///
/// .. doctest::
///
///     >>> from qiskit_bosons.mappers import LevelEncoding
///     >>> str(LevelEncoding.STANDARD_BINARY)
///     'standard_binary'
///     >>> LevelEncoding.from_name("gray")
///     LevelEncoding.GRAY
///     >>> LevelEncoding.from_name(str(LevelEncoding.UNARY_UNPROJECTED)) == LevelEncoding.UNARY_UNPROJECTED
///     True
///     >>> LevelEncoding.from_name("unary")
///     Traceback (most recent call last):
///     ...
///     ValueError: unknown encoding "unary"; expected one of 'unary_projected', 'unary_unprojected', 'standard_binary', 'gray'
///
/// .. [1] N. P. D. Sawaya et al., Resource-efficient digital quantum simulation of d-level systems
///        for photonic, vibrational, and spin-s Hamiltonians, npj Quantum Information 6, 49 (2020),
///        `arXiv:1909.12847 <https://arxiv.org/abs/1909.12847>`_.
#[gen_stub_pyclass_enum(module = "qiskit_bosons._lib.mappers.levels")]
#[pyclass(
    module = "qiskit_bosons.mappers.levels",
    name = "LevelEncoding",
    frozen,
    eq,
    eq_int,
    hash
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PyLevelEncoding {
    /// One qubit per level; level :math:`n` sets qubit :math:`n`. Every idle qubit carries a
    /// projector.
    ///
    /// A system with :math:`d` levels costs :math:`d` qubits, and only :math:`d` of the
    /// :math:`2^d` basis states of those qubits are ever populated. An operator that moves the level
    /// up or down by one has two non-diagonal factors, and this variant also puts a projector
    /// :math:`\lvert 0 \rangle\langle 0 \rvert` on each of the other qubits of the system. Its
    /// terms therefore touch all :math:`d` qubits, vanish outside the one-hot subspace, and are
    /// costlier to synthesize than those of :attr:`UNARY_UNPROJECTED`, by a margin that grows
    /// quickly with :math:`d`.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.mappers import LevelEncoding
    ///     >>> [LevelEncoding.UNARY_PROJECTED.num_qubits(d) for d in (1, 2, 3, 4, 5, 9)]
    ///     [1, 2, 3, 4, 5, 9]
    ///     >>> [format(LevelEncoding.UNARY_PROJECTED.level_to_bits(n), "05b") for n in range(5)]
    ///     ['00001', '00010', '00100', '01000', '10000']
    ///
    /// The bitmask of a level has qubit :math:`q` in bit :math:`q`, so qubit 0 is the rightmost digit.
    #[pyo3(name = "UNARY_PROJECTED")]
    UnaryProjected,
    /// One qubit per level, as for :attr:`UNARY_PROJECTED`, but without the projectors on the idle
    /// qubits.
    ///
    /// The bit patterns, the qubit count and the matrix elements are those of
    /// :attr:`UNARY_PROJECTED`. Only the projectors are left out, so an operator that moves the
    /// level up or down by one keeps the weight-2 form with two non-diagonal factors and acts as the
    /// identity on every other qubit, however large :math:`d` is. Its gate count then grows only
    /// linearly in :math:`d`.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.mappers import LevelEncoding
    ///     >>> [LevelEncoding.UNARY_UNPROJECTED.num_qubits(d) for d in (1, 2, 3, 4, 5, 9)]
    ///     [1, 2, 3, 4, 5, 9]
    ///     >>> [format(LevelEncoding.UNARY_UNPROJECTED.level_to_bits(n), "05b") for n in range(5)]
    ///     ['00001', '00010', '00100', '01000', '10000']
    ///
    /// The two unary variants are equal on the one-hot subspace. This one is **undefined outside
    /// it**, since nothing suppresses the evolution of an unphysical state, so choosing it is a
    /// modeling decision and not a pure optimization. It is only meaningful for the unary code.
    #[pyo3(name = "UNARY_UNPROJECTED")]
    UnaryUnprojected,
    /// The binary digits of the level.
    ///
    /// A system with :math:`d` levels costs :math:`\lceil \log_2 d \rceil` qubits. A system with a
    /// single level costs none, because its only level needs no qubit to be told apart from any
    /// other. When :math:`d` is not a power of two, the surplus basis states are never populated.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.mappers import LevelEncoding
    ///     >>> [LevelEncoding.STANDARD_BINARY.num_qubits(d) for d in (1, 2, 3, 4, 5, 9)]
    ///     [0, 1, 2, 2, 3, 4]
    ///     >>> [format(LevelEncoding.STANDARD_BINARY.level_to_bits(n), "03b") for n in range(5)]
    ///     ['000', '001', '010', '011', '100']
    ///
    /// Adjacent levels can differ in many bits (``011`` and ``100`` above differ in all three),
    /// which is what makes the operators that connect them higher-weight than under
    /// :attr:`GRAY`.
    #[pyo3(name = "STANDARD_BINARY")]
    StandardBinary,
    /// The binary-reflected Gray code of the level.
    ///
    /// It costs the same :math:`\lceil \log_2 d \rceil` qubits as :attr:`STANDARD_BINARY`, and the
    /// same zero qubits for a single level, but consecutive levels differ in exactly one bit.
    ///
    /// .. doctest::
    ///
    ///     >>> from qiskit_bosons.mappers import LevelEncoding
    ///     >>> [LevelEncoding.GRAY.num_qubits(d) for d in (1, 2, 3, 4, 5, 9)]
    ///     [0, 1, 2, 2, 3, 4]
    ///     >>> [format(LevelEncoding.GRAY.level_to_bits(n), "03b") for n in range(5)]
    ///     ['000', '001', '011', '010', '110']
    ///
    /// Each code in the list above differs from its neighbor in a single digit, which keeps the
    /// operators that connect adjacent levels low-weight.
    #[pyo3(name = "GRAY")]
    Gray,
}

impl PyLevelEncoding {
    /// The string names, in declaration order, as accepted by [`Self::from_name`].
    const NAMES: [&'static str; 4] = [
        "unary_projected",
        "unary_unprojected",
        "standard_binary",
        "gray",
    ];

    /// Returns the native encoding this variant stands for.
    ///
    /// The inverse of the `From<LevelEncoding>` conversion below; the two are the only places that
    /// pair the variants up.
    pub fn inner(self) -> LevelEncoding {
        match self {
            PyLevelEncoding::UnaryProjected => LevelEncoding::UnaryProjected,
            PyLevelEncoding::UnaryUnprojected => LevelEncoding::UnaryUnprojected,
            PyLevelEncoding::StandardBinary => LevelEncoding::StandardBinary,
            PyLevelEncoding::Gray => LevelEncoding::Gray,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            PyLevelEncoding::UnaryProjected => Self::NAMES[0],
            PyLevelEncoding::UnaryUnprojected => Self::NAMES[1],
            PyLevelEncoding::StandardBinary => Self::NAMES[2],
            PyLevelEncoding::Gray => Self::NAMES[3],
        }
    }

    fn parse(name: &str) -> PyResult<Self> {
        match name {
            "unary_projected" => Ok(PyLevelEncoding::UnaryProjected),
            "unary_unprojected" => Ok(PyLevelEncoding::UnaryUnprojected),
            "standard_binary" => Ok(PyLevelEncoding::StandardBinary),
            "gray" => Ok(PyLevelEncoding::Gray),
            other => Err(crate::value_err(format!(
                "unknown encoding {other:?}; expected one of 'unary_projected', \
                 'unary_unprojected', 'standard_binary', 'gray'"
            ))),
        }
    }
}

impl From<LevelEncoding> for PyLevelEncoding {
    /// Returns the variant that stands for a native encoding.
    fn from(encoding: LevelEncoding) -> Self {
        match encoding {
            LevelEncoding::UnaryProjected => PyLevelEncoding::UnaryProjected,
            LevelEncoding::UnaryUnprojected => PyLevelEncoding::UnaryUnprojected,
            LevelEncoding::StandardBinary => PyLevelEncoding::StandardBinary,
            LevelEncoding::Gray => PyLevelEncoding::Gray,
        }
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyLevelEncoding {
    /// Returns the number of qubits needed to represent a system with ``num_levels`` levels.
    ///
    /// Args:
    ///     num_levels: the number of levels :math:`d` of the system.
    ///
    /// Returns:
    ///     :math:`d` for the two unary encodings and :math:`\lceil \log_2 d \rceil` for the two
    ///     compact encodings. The compact encodings can address more levels than :math:`d` whenever
    ///     it is not a power of two; those surplus levels are never populated by the operators built
    ///     from the encoding.
    pub fn num_qubits(&self, num_levels: u32) -> u32 {
        self.inner().num_qubits(num_levels)
    }

    /// Returns the qubit basis state representing level ``level``.
    ///
    /// Args:
    ///     level: the level, counted from zero.
    ///
    /// Returns:
    ///     A bitmask whose bit :math:`q` is the value of qubit :math:`q` of the system.
    ///
    /// Raises:
    ///     ValueError: for a unary encoding if ``level`` is 128 or more, since the unary encodings
    ///         hold at most 128 levels.
    pub fn level_to_bits(&self, level: u32) -> PyResult<u128> {
        if self.inner().is_unary() && level >= MAX_UNARY_LEVELS {
            return Err(crate::value_err(format!(
                "the unary encoding holds at most {MAX_UNARY_LEVELS} levels, but level {level} \
                 was requested"
            )));
        }
        Ok(self.inner().level_to_bits(level))
    }

    /// Returns the encoding with the given string name.
    ///
    /// Args:
    ///     name: one of ``"unary_projected"``, ``"unary_unprojected"``, ``"standard_binary"`` or
    ///         ``"gray"``.
    ///
    /// Raises:
    ///     ValueError: if ``name`` is not a recognized name.
    #[staticmethod]
    pub fn from_name(name: &str) -> PyResult<Self> {
        Self::parse(name)
    }

    fn __str__(&self) -> &'static str {
        self.name()
    }

    /// Returns the callable and arguments that rebuild this encoding when unpickled.
    ///
    /// The encoding is rebuilt from its string name, so a pickle stays valid however the variants
    /// are numbered internally.
    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, (&'static str,))> {
        let from_name = py.get_type::<Self>().getattr("from_name")?;
        Ok((from_name, (self.name(),)))
    }

    pub(crate) fn __repr__(&self) -> String {
        let variant = match self {
            PyLevelEncoding::UnaryProjected => "UNARY_PROJECTED",
            PyLevelEncoding::UnaryUnprojected => "UNARY_UNPROJECTED",
            PyLevelEncoding::StandardBinary => "STANDARD_BINARY",
            PyLevelEncoding::Gray => "GRAY",
        };
        format!("LevelEncoding.{variant}")
    }
}

/// A level encoding given either as a [`PyLevelEncoding`] or as its string name.
///
/// This is the single point at which the Python-facing argument becomes the native enum
/// ([`Self::resolve`]), so every binding that takes an encoding accepts exactly the same spellings
/// and raises the same errors.
#[derive(FromPyObject)]
pub enum EncodingArg {
    Level(PyLevelEncoding),
    Name(String),
}

impl EncodingArg {
    /// Resolves the argument to the native encoding.
    ///
    /// There is deliberately no default encoding: the choice of layout is the caller's to make.
    ///
    /// Raises a `ValueError` for a string that is not a recognized name.
    pub fn resolve(self) -> PyResult<LevelEncoding> {
        match self {
            EncodingArg::Level(level_encoding) => Ok(level_encoding.inner()),
            EncodingArg::Name(name) => Ok(PyLevelEncoding::parse(&name)?.inner()),
        }
    }
}
