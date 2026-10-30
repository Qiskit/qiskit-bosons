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

//! The native core of `qiskit-bosons`.
//!
//! This crate holds all native logic and **must not depend on `pyo3`**: the bindings live in
//! sibling crates (`pyext` today, a `cext` in the future) so that the core stays independent of any
//! single binding technology.
//!
//! # FFI backends
//!
//! Reaching Qiskit's `SparseObservable` from native code requires an FFI backend, which is
//! feature-gated behind a single `ffi` alias:
//!
//! - `pyext` -> `qiskit-pyo3-ffi` (the only backend wired up today)
//! - `cext` -> `qiskit-sys` (not currently present; see below)
//!
//! Only the `pyext` backend exists at the moment, but the *seam* is deliberately preserved: every
//! FFI-touching module selects its backend through `cfg_select!` on these features rather than
//! importing a binding crate directly. Reintroducing a C API therefore means adding a `cext`
//! feature plus its `cfg_select!` arm -- not untangling Python from the core.
//!
//! The backends are **mutually exclusive**: they bind the same `ffi` alias to different crates, so
//! enabling both would let Cargo's feature unification compile against the wrong one. The guard
//! below fails the build with a clear message rather than a confusing downstream type mismatch. It
//! is written against `cext` even though no such feature is declared yet, so that the invariant is
//! enforced from the moment that feature reappears.
#[cfg(all(feature = "cext", feature = "pyext"))]
compile_error!("features `cext` and `pyext` are mutually exclusive; enable exactly one");

pub mod linalg;
pub mod mappers;
pub mod operators;
pub mod random;
pub mod testing;
