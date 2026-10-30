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

use crate::operators::boson_operator::{BosonOperator, BosonOperatorTermView};
use crate::operators::terms::filtering::retain_terms;
use std::iter::zip;

/// Returns whether a *normal-ordered* term is diagonal in the occupation-number basis.
///
/// A term is diagonal exactly when it is a product of number operators, i.e. every acted-upon mode
/// appears equally many times as a creation and as an annihilation. For a normal-ordered term (all
/// creations preceding all annihilations) this is equivalent to the multiset of creation modes
/// being equal to the multiset of annihilation modes.
fn _is_diagonal(term: BosonOperatorTermView) -> bool {
    let mut creations: Vec<u32> = zip(term.actions, term.modes)
        .filter_map(|(&action, &mode)| if action { Some(mode) } else { None })
        .collect();
    let mut annihilations: Vec<u32> = zip(term.actions, term.modes)
        .filter_map(|(&action, &mode)| if !action { Some(mode) } else { None })
        .collect();
    creations.sort_unstable();
    annihilations.sort_unstable();
    creations == annihilations
}

/// Filters out the terms of an operator that are diagonal in the occupation-number basis.
///
/// Every term that is a product of number operators (`b^\dagger_i b_i`) is dropped from `op`. This
/// includes the constant term (a product of zero number operators), single number operators, as
/// well as higher-order products such as `n_i n_j = b^\dagger_i b^\dagger_j b_j b_i`. Such terms
/// are diagonal in the occupation-number basis.
///
/// # Assumptions
///
/// The terms of `op` are assumed to be *normal-ordered*! This is *not* being verified. See
/// [`BosonOperator::normal_ordered`] for how to obtain an operator of that form.
///
/// This is mutating `op` in place. If `op` tracks group indices (see
/// [`BosonOperator::groups`]), surviving terms retain their relative grouping but the group
/// indices are reassigned to a contiguous range starting from 0; callers must therefore *not* rely
/// on the specific group index of any term being preserved across a call to this function.
pub fn filter_diagonal_terms(op: &mut BosonOperator) {
    retain_terms(op, |term| !_is_diagonal(term));
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::*;

    use crate::operators::OperatorTrait;

    #[test]
    fn test_filter_drops_constant_number_and_products() {
        // Terms (in order):
        //   - constant                              (diagonal -> drop)
        //   - n_0 = b†_0 b_0                         (diagonal -> drop)
        //   - n_0 n_1 = b†_0 b†_1 b_1 b_0            (diagonal -> drop)
        //   - b†_0 b_1                               (off-diagonal -> keep)
        let mut op = BosonOperator {
            coeffs: vec![
                Complex64::new(1.0, 0.0),
                Complex64::new(2.0, 0.0),
                Complex64::new(3.0, 0.0),
                Complex64::new(4.0, 0.0),
            ],
            actions: vec![true, false, true, true, false, false, true, false],
            modes: vec![0, 0, 0, 1, 1, 0, 0, 1],
            boundaries: vec![0, 0, 2, 6, 8],
            groups: None,
        };

        filter_diagonal_terms(&mut op);

        let expected = BosonOperator {
            coeffs: vec![Complex64::new(4.0, 0.0)],
            actions: vec![true, false],
            modes: vec![0, 1],
            boundaries: vec![0, 2],
            groups: None,
        };
        assert_eq!(op, expected);
    }

    #[test]
    fn test_filter_keeps_off_diagonal_hopping() {
        // A single off-diagonal hopping term must be kept untouched.
        let mut op = BosonOperator {
            coeffs: vec![Complex64::new(1.0, 0.0)],
            actions: vec![true, false],
            modes: vec![0, 1],
            boundaries: vec![0, 2],
            groups: None,
        };
        let expected = op.clone();
        filter_diagonal_terms(&mut op);
        assert_eq!(op, expected);
    }

    #[test]
    fn test_filter_none_groups_stay_none() {
        let mut op = BosonOperator {
            coeffs: vec![Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)],
            actions: vec![true, false, true, false],
            modes: vec![0, 0, 0, 1],
            boundaries: vec![0, 2, 4],
            groups: None,
        };
        filter_diagonal_terms(&mut op);
        assert!(op.groups.is_none());
    }

    #[test]
    fn test_filter_reindexes_groups() {
        // Two off-diagonal hopping terms sharing group 5, with a diagonal number operator in
        // group 2 in between. The number operator is dropped; the survivors must remain grouped
        // together and be re-indexed to a contiguous index (0).
        let mut op = BosonOperator {
            coeffs: vec![
                Complex64::new(1.0, 0.0),
                Complex64::new(2.0, 0.0),
                Complex64::new(3.0, 0.0),
            ],
            actions: vec![true, false, true, false, true, false],
            modes: vec![0, 1, 0, 0, 1, 0],
            boundaries: vec![0, 2, 4, 6],
            groups: Some(vec![5, 2, 5]),
        };

        filter_diagonal_terms(&mut op);

        assert_eq!(op.groups, Some(vec![0, 0]));
        assert_eq!(op.num_groups(), Some(1));
    }

    /// A repeated mode index is diagonal only if it appears the same number of times on both sides.
    ///
    /// This is the case a fermionic implementation never has to handle: `b†_0 b†_0 b_0 b_0` (the
    /// `n_0 (n_0 - 1)` interaction of a Bose-Hubbard model) is diagonal, while `b†_0 b†_0 b_0` is
    /// not. Multiset -- not set -- equality is what distinguishes them.
    #[test]
    fn test_filter_handles_repeated_modes() {
        let mut op = BosonOperator::zero();
        // b†_0 b†_0 b_0 b_0  -> diagonal (drop)
        op._append_term(
            Complex64::new(1.0, 0.0),
            &[true, true, false, false],
            &[0, 0, 0, 0],
        );
        // b†_0 b†_0 b_0      -> NOT diagonal (keep): unequal multiplicities
        op._append_term(Complex64::new(2.0, 0.0), &[true, true, false], &[0, 0, 0]);
        // b†_0 b†_1 b_1 b_0  -> diagonal (drop): n_0 n_1
        op._append_term(
            Complex64::new(3.0, 0.0),
            &[true, true, false, false],
            &[0, 1, 1, 0],
        );

        filter_diagonal_terms(&mut op);

        assert_eq!(op.coeffs, vec![Complex64::new(2.0, 0.0)]);
        for term in op.iter() {
            assert!(!_is_diagonal(term), "a diagonal term survived filtering");
        }
    }

    /// Filtering a Bose-Hubbard Hamiltonian must leave exactly the hopping terms.
    #[test]
    fn test_filter_bose_hubbard() {
        // H = -t sum_<ij> (b†_i b_j + b†_j b_i) + U/2 sum_i n_i (n_i - 1)
        // on two sites, normal-ordered.
        let t = Complex64::new(-1.0, 0.0);
        let u = Complex64::new(2.0, 0.0);
        let mut op = BosonOperator::zero();
        // hopping (off-diagonal)
        op._append_term(t, &[true, false], &[0, 1]);
        op._append_term(t, &[true, false], &[1, 0]);
        // on-site interaction (diagonal)
        op._append_term(u, &[true, true, false, false], &[0, 0, 0, 0]);
        op._append_term(u, &[true, true, false, false], &[1, 1, 1, 1]);
        // chemical-potential-like number terms (diagonal)
        op._append_term(Complex64::new(0.5, 0.0), &[true, false], &[0, 0]);
        op._append_term(Complex64::new(0.5, 0.0), &[true, false], &[1, 1]);

        let num_terms_before = op.coeffs.len();
        filter_diagonal_terms(&mut op);

        assert_eq!(
            op.coeffs.len(),
            2,
            "only the two hopping terms should survive"
        );
        assert!(op.coeffs.len() < num_terms_before);
        for term in op.iter() {
            assert!(!_is_diagonal(term), "a diagonal term survived filtering");
        }
    }
}
