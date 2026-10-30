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

//! Systems of `d` levels on qubits.
//!
//! Everything in this module is independent of what the levels *are*: the occupation numbers of a
//! truncated bosonic mode, vibrational modals or the points of a position grid all need the same
//! machinery. A mapper in [`library`](crate::mappers::library) decides which operator acts on the
//! levels and builds its matrix elements; the pieces here do the rest.
//!
//! - [`encoding`] holds [`LevelEncoding`](encoding::LevelEncoding), the layout of the levels of one
//!   system on qubits.
//! - [`spec`] holds [`LevelSpec`](spec::LevelSpec), which gives every mode of a register its own
//!   number of levels and derives from them how many qubits each mode occupies and where its qubits
//!   start.
//! - [`matrix`] turns a [`MatrixElement`](matrix::MatrixElement) of an operator on one system into
//!   qubit operators, given the bit patterns an encoding assigns to its levels.
//!
//! Only [`matrix`] reaches the Qiskit FFI backend. The other two are plain Rust.

pub mod encoding;
pub mod matrix;
pub mod spec;
