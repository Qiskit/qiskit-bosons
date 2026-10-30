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

// `commutators` attaches staticmethods onto the operator pyclasses themselves (see the
// `impl_commutators!` macro), so it exports no module of its own -- hence the empty `#[pymodule]`
// below. It must still be declared for those `#[pymethods]` blocks to be compiled in.
pub mod commutators;

// NOTE: the operator library exposes no native *module* members today. The module stays registered
// so that the public `qiskit_bosons.operators.library` namespace exists and native additions have a
// home.
#[pymodule]
pub mod operators_library {}
