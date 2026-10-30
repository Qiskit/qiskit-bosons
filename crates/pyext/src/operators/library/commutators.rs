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
use qiskit_bosons_core::operators::library::commutators::{
    anti_commutator, commutator, double_commutator,
};

/// Attaches the three commutator staticmethods that back the `SupportsCommutators` protocol.
///
/// Only one operator type exists today, but the macro is kept so that additional operator types get
/// the protocol by adding a single line -- and so the methods stay verbatim-identical across types.
macro_rules! impl_commutators {
    ($name:ty) => {
        #[gen_stub_pymethods]
        #[pymethods]
        impl $name {
            #[staticmethod]
            pub fn _commutator_(op_a: Self, op_b: Self) -> Self {
                commutator(&op_a.inner, &op_b.inner).into()
            }

            #[staticmethod]
            pub fn _anti_commutator_(op_a: Self, op_b: Self) -> Self {
                anti_commutator(&op_a.inner, &op_b.inner).into()
            }

            #[staticmethod]
            pub fn _double_commutator_(op_a: Self, op_b: Self, op_c: Self, sign: bool) -> Self {
                double_commutator(&op_a.inner, &op_b.inner, &op_c.inner, sign).into()
            }
        }
    };
}

impl_commutators!(PyBosonOperator);
