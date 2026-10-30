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

//! Commutator generators.
//!
//! These are generic over [`OperatorMacro`], so they are algebra-agnostic: the bosonic commutation
//! relations enter through [`BosonOperator`](crate::operators::boson_operator::BosonOperator)'s own
//! composition and normal ordering, not here. The bosonic algebra is *defined* by
//! `[b_i, b_j^dag] = delta_ij`, so the commutator of two operators is the natural way to interrogate
//! it.
//!
//! # Product order
//!
//! These functions build their products with [`OperatorMacro::__matmul__`], which is the
//! mathematical product: `a.__matmul__(b)` is `A B`. So `commutator(op_a, op_b)` is `[A, B]`, as the
//! names suggest.
//!
//! This is deliberately **not** [`OperatorMacro::__and__`]. Following Qiskit, `__and__` is
//! `compose`, i.e. *left* multiplication: `A & B == B @ A`. Composing the brackets out of it would
//! silently return `[B, A] = -[A, B]`, negating every commutator.

use crate::operators::{OperatorMacro, OperatorTrait, ScaledTerm};
use num_complex::Complex64;

pub fn commutator<T>(op_a: &T, op_b: &T) -> T
where
    T: OperatorMacro,
{
    op_a.__matmul__(op_b).__sub__(&op_b.__matmul__(op_a))
}

pub fn anti_commutator<T>(op_a: &T, op_b: &T) -> T
where
    T: OperatorMacro,
{
    op_a.__matmul__(op_b).__add__(&op_b.__matmul__(op_a))
}

pub fn double_commutator<T>(op_a: &T, op_b: &T, op_c: &T, sign: bool) -> T
where
    T: OperatorMacro + OperatorTrait,
{
    let sign_num = if sign {
        Complex64::new(1.0, 0.0)
    } else {
        Complex64::new(-1.0, 0.0)
    };

    let op_ab = op_a.__matmul__(op_b);
    let op_ba = op_b.__matmul__(op_a);
    let op_ac = op_a.__matmul__(op_c);
    let op_ca = op_c.__matmul__(op_a);

    let op_abc = op_ab.__matmul__(op_c);
    let op_cab = op_c.__matmul__(&op_ab);
    let op_bac = op_ba.__matmul__(op_c);
    let op_cba = op_c.__matmul__(&op_ba);
    let op_acb = op_ac.__matmul__(op_b);
    let op_bca = op_b.__matmul__(&op_ca);

    // Rescales every term of `op` by `weight`, allocating nothing: a term view owns its coefficient
    // and only borrows the index data. A local `fn` rather than a closure because the borrow of
    // `op` outlives the call, and a closure's return type cannot be generic over that lifetime.
    fn weighted<'a, T: OperatorTrait>(
        op: &'a T,
        weight: Complex64,
    ) -> impl Iterator<Item = T::TermView<'a>> {
        op.iter().map(move |term| term.scaled(weight))
    }

    let half = Complex64::new(0.5, 0.0);

    // The result is the six triple products, each carrying a fixed weight:
    //
    //     ABC - s CBA - BAC/2 + s CAB/2 - ACB/2 + s BCA/2,    s = sign_num = +-1
    //
    // which reproduces both documented forms: `sign = false` gives
    // `(2 ABC + 2 CBA - BAC - CAB - ACB - BCA)/2` and `sign = true` gives
    // `(2 ABC - 2 CBA - BAC + CAB - ACB + BCA)/2`. Mind the `+s BCA/2`: the four half-weighted
    // products form a single sum, so grouping them into two halves in order to subtract one from
    // the other would flip that sign. Folding each weight into its own block removes that hazard
    // rather than merely warning about it.
    //
    // Appending all six blocks in a single pass replaces a chain of nine
    // `__mul__`/`__add__`/`__sub__` calls, each of which copied its left operand - an operand that
    // grew with every step, so early terms were copied over and over. Here every term is copied
    // exactly once. All weights are +-1 or +-1/2, so folding them into the coefficients is exact.
    T::from_terms(
        weighted(&op_abc, Complex64::new(1.0, 0.0))
            .chain(weighted(&op_cba, -sign_num))
            .chain(weighted(&op_bac, -half))
            .chain(weighted(&op_cab, sign_num * half))
            .chain(weighted(&op_acb, -half))
            .chain(weighted(&op_bca, sign_num * half)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::operators::OperatorTrait;
    use crate::operators::boson_operator::BosonOperator;

    fn op(terms: &[(Complex64, &[(bool, u32)])]) -> BosonOperator {
        let mut out = BosonOperator::zero();
        for (coeff, factors) in terms {
            let actions: Vec<bool> = factors.iter().map(|f| f.0).collect();
            let modes: Vec<u32> = factors.iter().map(|f| f.1).collect();
            out._append_term(*coeff, &actions, &modes);
        }
        out
    }

    fn c(re: f64) -> Complex64 {
        Complex64::new(re, 0.0)
    }

    const CRE: bool = true;
    const ANN: bool = false;

    /// The defining bosonic relation: `[b_0, b_0^dag] = 1`.
    ///
    /// This is the whole point of having a commutator generator for a bosonic package -- it lets the
    /// canonical commutation relation be *checked* rather than assumed.
    #[test]
    fn test_canonical_commutation_relation() {
        let b = op(&[(c(1.0), &[(ANN, 0)])]);
        let bdag = op(&[(c(1.0), &[(CRE, 0)])]);

        let comm = commutator(&b, &bdag).normal_ordered(None);
        assert!(comm.equiv(&BosonOperator::one(), 1e-12));
        // antisymmetry: the reversed order is `[b^dag, b] = -1`
        let comm = commutator(&bdag, &b).normal_ordered(None);
        assert!(comm.equiv(&BosonOperator::one().__neg__(), 1e-12));
    }

    /// Ladder operators on distinct modes commute: `[b_0, b_1^dag] = 0`.
    #[test]
    fn test_distinct_modes_commute() {
        let b0 = op(&[(c(1.0), &[(ANN, 0)])]);
        let b1dag = op(&[(c(1.0), &[(CRE, 1)])]);

        let comm = commutator(&b0, &b1dag).normal_ordered(None);
        assert!(comm.equiv(&BosonOperator::zero(), 1e-12));

        // ...and so do two creations, with no sign flip
        let b0dag = op(&[(c(1.0), &[(CRE, 0)])]);
        assert!(
            commutator(&b0dag, &b1dag)
                .normal_ordered(None)
                .equiv(&BosonOperator::zero(), 1e-12)
        );
    }

    /// `[n_0, b_0^dag] = b_0^dag`: the number operator raises by one.
    #[test]
    fn test_number_operator_commutator() {
        let n0 = op(&[(c(1.0), &[(CRE, 0), (ANN, 0)])]);
        let bdag = op(&[(c(1.0), &[(CRE, 0)])]);

        let comm = commutator(&n0, &bdag).normal_ordered(None);
        assert!(comm.equiv(&bdag, 1e-12));

        // and lowers by one for the annihilation operator: `[n_0, b_0] = -b_0`
        let b = op(&[(c(1.0), &[(ANN, 0)])]);
        let comm = commutator(&n0, &b).normal_ordered(None);
        assert!(comm.equiv(&op(&[(c(-1.0), &[(ANN, 0)])]), 1e-12));
    }

    /// The anti-commutator is *not* the natural bosonic bracket, but it is still well defined.
    ///
    /// `{b_0, b_0^dag} = 2 b_0^dag b_0 + 1`.
    #[test]
    fn test_anti_commutator() {
        let b = op(&[(c(1.0), &[(ANN, 0)])]);
        let bdag = op(&[(c(1.0), &[(CRE, 0)])]);

        let anti = anti_commutator(&b, &bdag).normal_ordered(None);
        let expected = op(&[(c(2.0), &[(CRE, 0), (ANN, 0)]), (c(1.0), &[])]);
        assert!(anti.equiv(&expected, 1e-12));
    }

    /// An operator commutes with itself, and with any scalar multiple of itself.
    #[test]
    fn test_self_commutator_vanishes() {
        let o = op(&[
            (c(1.0), &[(CRE, 0), (ANN, 1)]),
            (c(2.0), &[(CRE, 1), (ANN, 0)]),
        ]);
        assert!(
            commutator(&o, &o)
                .normal_ordered(None)
                .equiv(&BosonOperator::zero(), 1e-12)
        );
        let scaled = o.__mul__(c(3.0));
        assert!(
            commutator(&o, &scaled)
                .normal_ordered(None)
                .equiv(&BosonOperator::zero(), 1e-12)
        );
    }

    /// Antisymmetry: `[A, B] = -[B, A]`.
    #[test]
    fn test_commutator_antisymmetry() {
        let a = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        let b = op(&[(c(1.0), &[(CRE, 1), (ANN, 0)])]);

        let ab = commutator(&a, &b).normal_ordered(None);
        let ba = commutator(&b, &a).normal_ordered(None);
        assert!(ab.equiv(&ba.__neg__(), 1e-12));
    }

    /// The total number operator commutes with any number-conserving term.
    #[test]
    fn test_number_conservation_via_commutator() {
        let total_n = op(&[
            (c(1.0), &[(CRE, 0), (ANN, 0)]),
            (c(1.0), &[(CRE, 1), (ANN, 1)]),
        ]);
        // hopping conserves the total number
        let hop = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        assert!(
            commutator(&total_n, &hop)
                .normal_ordered(None)
                .equiv(&BosonOperator::zero(), 1e-12)
        );
        // a bare creation does not
        let bare = op(&[(c(1.0), &[(CRE, 0)])]);
        assert!(
            !commutator(&total_n, &bare)
                .normal_ordered(None)
                .equiv(&BosonOperator::zero(), 1e-12)
        );
    }

    /// Pins the six weights of the double commutator, independently of term order.
    ///
    /// Builds `ABC - s CBA - BAC/2 + s CAB/2 - ACB/2 + s BCA/2` explicitly from the ten products
    /// and compares with `equiv`. This is the check that catches a wrong sign or a wrong weight --
    /// notably the `+s BCA/2` term, whose sign a two-halves grouping of the half-weighted sum would
    /// flip.
    fn assert_double_commutator_weights<T>(op_a: &T, op_b: &T, op_c: &T)
    where
        T: OperatorMacro + OperatorTrait + std::fmt::Debug,
    {
        let half = Complex64::new(0.5, 0.0);
        for sign in [false, true] {
            let sign_num = if sign {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(-1.0, 0.0)
            };

            let op_ab = op_a.__matmul__(op_b);
            let op_ba = op_b.__matmul__(op_a);
            let op_ac = op_a.__matmul__(op_c);
            let op_ca = op_c.__matmul__(op_a);

            let expected = op_ab
                .__matmul__(op_c)
                .__sub__(&op_c.__matmul__(&op_ba).__mul__(sign_num))
                .__sub__(&op_ba.__matmul__(op_c).__mul__(half))
                .__add__(&op_c.__matmul__(&op_ab).__mul__(sign_num * half))
                .__sub__(&op_ac.__matmul__(op_b).__mul__(half))
                .__add__(&op_b.__matmul__(&op_ca).__mul__(sign_num * half));

            let actual = double_commutator(op_a, op_b, op_c, sign);
            assert!(
                actual.equiv(&expected, 1e-12),
                "weights disagree for sign={sign}"
            );
            // Guard against a slip to `from_terms_with_groups`: every product comes from `matmul`,
            // which tracks no groups, so neither may the sum.
            assert!(actual.groups().is_none());
        }
    }

    #[test]
    fn test_double_commutator_weights() {
        let op_a = op(&[
            (Complex64::new(1.0, 0.25), &[(CRE, 0), (ANN, 1)]),
            (Complex64::new(-0.4, 0.0), &[(CRE, 2)]),
        ]);
        let op_b = op(&[(Complex64::new(0.5, -0.75), &[(CRE, 1), (CRE, 2)])]);
        let op_c = op(&[
            (Complex64::new(-0.3, 1.5), &[(ANN, 2)]),
            (Complex64::new(0.0, 0.6), &[(ANN, 0), (CRE, 1)]),
        ]);
        assert_double_commutator_weights(&op_a, &op_b, &op_c);
    }

    /// The double commutator is linear in each argument.
    ///
    /// Rather than assert a closed form (the formula is McWeeny's Eq. (13.6.18) and does *not*
    /// collapse to a plain nested commutator -- with `C = 1` and `sign = false` it reduces to the
    /// bare product `A B`, not to zero), this pins the structural property that must hold for any
    /// correct implementation: scaling one argument scales the result.
    #[test]
    fn test_double_commutator_is_linear() {
        let a = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        let b = op(&[(c(1.0), &[(CRE, 1), (ANN, 0)])]);
        let cc = op(&[(c(1.0), &[(CRE, 0), (ANN, 0)])]);

        for sign in [true, false] {
            let base = double_commutator(&a, &b, &cc, sign).normal_ordered(None);
            let scaled = double_commutator(&a.__mul__(c(3.0)), &b, &cc, sign).normal_ordered(None);
            assert!(
                scaled.equiv(&base.__mul__(c(3.0)), 1e-12),
                "double_commutator must be linear in its first argument (sign={sign})"
            );
        }
    }

    /// The two `sign` variants genuinely differ, and both are non-trivial.
    ///
    /// A guard against an implementation that ignores `sign`, or that collapses to zero.
    #[test]
    fn test_double_commutator_sign_matters() {
        let a = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        let b = op(&[(c(1.0), &[(CRE, 1), (ANN, 0)])]);
        let cc = op(&[(c(1.0), &[(CRE, 0), (ANN, 0)])]);

        let plus = double_commutator(&a, &b, &cc, true).normal_ordered(None);
        let minus = double_commutator(&a, &b, &cc, false).normal_ordered(None);

        assert!(!plus.equiv(&BosonOperator::zero(), 1e-12));
        assert!(!minus.equiv(&BosonOperator::zero(), 1e-12));
        assert!(
            !plus.equiv(&minus, 1e-12),
            "the `sign` flag must change the result"
        );
    }
}
