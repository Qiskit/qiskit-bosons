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

pub mod encoding;
pub mod spec;

// The physical module holding the framework-level level types.
//
// The classes declare their public path through `#[pyclass(module = ...)]`, and
// `python/qiskit_bosons/mappers/__init__.py` re-exports them, so a user sees
// `qiskit_bosons.mappers.LevelEncoding` either way.
#[pymodule]
pub mod levels {
    #[pymodule_export]
    use super::encoding::PyLevelEncoding;

    #[pymodule_export]
    use super::spec::PyLevelSpec;
}
