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

//! The occupation-basis boson-to-qubit mapping.
//!
//! This mapper works in the Fock basis: level `n` of a mode is the state with `n` bosons. How those
//! levels are laid out on qubits is the choice of a [`LevelEncoding`]; this module supplies the
//! bosonic part, namely the truncated matrices of the ladder operators.
//!
//! # Encodings map states, not algebras
//!
//! A bosonic mode has an **infinite-dimensional** Fock space, so any qubit encoding must first
//! truncate it at some maximum occupation `n_max` (giving `d = n_max + 1` levels). Truncation
//! breaks the algebra. Writing `b̄` for the truncated annihilation operator, Somma et al. give
//! (arXiv:quant-ph/0304063, Sec. 4, Eq. 22)
//!
//! ```text
//! [b̄_i, b̄_j] = 0,   [b̄_i, b̄†_j] = δ_ij [ 1 - ((N_P+1)/N_P!) (b̄†_i)^{N_P} (b̄_i)^{N_P} ]
//! ```
//!
//! with `N_P = d - 1`. Equivalently, and more usefully for testing, the defect is a single rank-one
//! projector onto the top level:
//!
//! ```text
//! [b̄, b̄†] = 1 - d |d-1><d-1|
//! ```
//!
//! so the canonical relation holds on every level *except* the highest, where it is off by `-d`.
//! The practical consequences are worth stating plainly:
//!
//! 1. The span of `b̄†, b̄` is **not closed under the commutator**. Somma et al. conclude from this
//!    that "a direct mapping between the bosonic algebra and the spin-1/2 algebra (such as [the]
//!    Jordan-Wigner [transformation]) is not possible". A bosonic encoding is therefore a mapping
//!    of *states* (and of operators acting within the truncated subspace), not an algebra
//!    isomorphism.
//! 2. **Mapping does not commute with normal ordering.** Encoding an operator and encoding its
//!    normal-ordered form give different qubit operators whenever the reordering requires a
//!    contraction, differing by exactly the projector defect above: truncated `b b^dag` has top
//!    diagonal entry `0` (`b^dag` sends the top level out of the retained space, so there is nothing
//!    left for `b` to act on) while `b^dag b + 1` has `d`. The two agree on every other level, and
//!    differ by `d` on the top one. Normal ordering is exact only in the untruncated space, so
//!    truncating before versus after applying it is not the same operation.
//! 3. `d` must be chosen to bound the total occupation *reachable* by the operator being simulated,
//!    not merely the occupation of the initial state. A two-mode beamsplitter acting on `|1,1>` can
//!    populate `|2,0>`, so `d = 2` is not sufficient for that state even though no mode starts above
//!    occupation 1.
//!
//! Because truncation is a modeling choice rather than an implementation detail, the cutoffs of a
//! [`LevelSpec`] are required, explicit parameters of every mapping in this module -- never a silent
//! default. Each mode has its own cutoff, so the reachable occupation must be bounded mode by mode.
//!
//! # Choosing an encoding
//!
//! There is no universally best [`LevelEncoding`]; see its documentation for the trade-offs. For
//! bosonic modes the unary encodings give every ladder term only two non-diagonal factors at a cost
//! of `d` qubits per mode, while the compact encodings need only `ceil(log2 d)`. The unary encoding
//! comes in two variants. [`LevelEncoding::UnaryProjected`] keeps a projector on the idle qubits of
//! the mode and is the costlier of the two to synthesize. [`LevelEncoding::UnaryUnprojected`] drops
//! them, giving the weight-2 form, and is equal to the projected one on the one-hot subspace but
//! undefined outside it.
//!
//! # How a term is mapped
//!
//! A term of a [`BosonOperator`] is a product of ladder operators spread over several modes. Because
//! operators acting on *different* modes commute, the term factorizes: the ladder operators are
//! bucketed by mode, each per-mode product is contracted into a single truncated `d x d` matrix (via
//! [`single_mode_matrix`], which multiplies *before* encoding as Sawaya et al. recommend), and the
//! per-mode matrices are then encoded and composed.
//!
//! The bucketing is valid precisely because distinct modes commute. Within one mode the original
//! operator order is preserved exactly, which is what makes `b b†` and `b† b` map to different (and
//! correctly different) qubit operators.
//!
//! Each matrix element `<row| O |col>` becomes Pauli and projector terms through
//! [`encode_matrix_element`], which is shared by every mapper in a level basis; see
//! [`crate::mappers::levels::matrix`] for how a single element is expanded. For the unary encodings
//! a ladder term has exactly two non-diagonal factors (`X`/`Y`) per mode, whatever the truncation
//! level. Under [`LevelEncoding::UnaryProjected`] it also carries a projector `|0><0|` on each of the
//! mode's other qubits, so it touches all `d` of them; the projectors are stored as single bit
//! terms, but expanding them into Paulis costs two terms each. Under
//! [`LevelEncoding::UnaryUnprojected`] those projectors are left out, so the term is weight 2, which
//! is the form behind the `O(d)` gate scaling of the unary code.

use crate::mappers::levels::encoding::LevelEncoding;
use crate::mappers::levels::matrix::{EncodedTerm, MatrixElement, encode_matrix_element};
use crate::mappers::levels::spec::{LevelSpec, LevelSpecError, ModeLayout};
use crate::mappers::qk_obs::{QkComplex64, ffi_backend as ffi, map_operator, pauli_sum, qk_coeff};
use crate::operators::OperatorTrait;
use crate::operators::boson_operator::BosonOperator;
use num_complex::Complex64;
use std::collections::BTreeMap;

/// Builds the dense `d x d` matrix of a product of ladder operators acting on a single mode.
///
/// `actions` lists the ladder operators applied to this mode, **in operator order** (i.e. rightmost
/// element of the mathematical product first in application order); `true` denotes a creation
/// operator. The result is truncated to `d = max_occupation + 1` levels: any amplitude that would
/// leave the retained range is dropped, which is exactly the truncation whose algebraic cost is
/// described in the module documentation.
///
/// # Why multiply here rather than after encoding
///
/// Sawaya et al. note that for non-unary encodings one must "multiply operators in the matrix
/// representation before encoding": encoding each factor separately and composing the resulting
/// Pauli sums yields algebraically correct but needlessly long expressions (their example is
/// `n^2` at `d = 3`). Composing in the `d x d` matrix representation first and encoding the product
/// once keeps the output minimal. Doing so also makes the truncation unambiguous: the product is
/// truncated once, as a whole.
pub fn single_mode_matrix(actions: &[bool], max_occupation: u32) -> Vec<MatrixElement> {
    let dim = max_occupation + 1;
    // Track, for each starting level, the level it has been carried to and the accumulated
    // amplitude. `None` marks a column whose amplitude has left the truncated range.
    let mut carried: Vec<Option<(u32, f64)>> = (0..dim).map(|level| Some((level, 1.0))).collect();

    // Ladder operators are listed in operator order, so the *last* entry acts first.
    for &is_creation in actions.iter().rev() {
        for slot in carried.iter_mut() {
            let Some((level, amplitude)) = *slot else {
                continue;
            };
            if is_creation {
                // b† |n> = sqrt(n+1) |n+1>, dropped if it would exceed the truncation
                let next = level + 1;
                if next >= dim {
                    *slot = None;
                } else {
                    *slot = Some((next, amplitude * f64::from(next).sqrt()));
                }
            } else if level == 0 {
                // b |0> = 0 -- an exact algebraic zero, not a truncation artifact
                *slot = None;
            } else {
                // b |n> = sqrt(n) |n-1>
                *slot = Some((level - 1, amplitude * f64::from(level).sqrt()));
            }
        }
    }

    carried
        .into_iter()
        .enumerate()
        .filter_map(|(col, carried)| {
            carried.map(|(row, amplitude)| MatrixElement {
                row,
                col: col as u32,
                value: Complex64::new(amplitude, 0.0),
            })
        })
        .collect()
}

/// Maps a [`BosonOperator`] onto a Qiskit `SparseObservable` under an occupation encoding.
///
/// Each mode is truncated to the levels its cutoff in `spec` keeps and laid out on the qubits the
/// spec assigns it: mode `m` starts at the spec's offset for `m`, so modes with different cutoffs
/// may take a different number of qubits. A uniform spec is applied to the `op.min_num_modes()`
/// modes the operator needs; a per-mode spec keeps all of its own modes and must cover at least
/// that many.
///
/// The register has `num_qubits` qubits if given, which may exceed what the modes need (the
/// surplus qubits are idle and sit at the end), and is exactly as wide as the modes otherwise.
///
/// # Truncation
///
/// The cutoffs are required, with no default: truncation is a physical modeling decision, and each
/// must bound the occupation *reachable* in its mode by the operator, not merely that of the initial
/// state. See the module documentation for the exact algebraic cost of truncating.
///
/// # Errors
///
/// Returns a [`LevelSpecError`] if `op` acts on more modes than a per-mode `spec` covers, or if
/// `num_qubits` is smaller than the register the modes need.
pub fn boson_occupation_encoding(
    op: &BosonOperator,
    spec: &LevelSpec,
    num_qubits: Option<u32>,
) -> Result<*mut ffi::QkObs, LevelSpecError> {
    // This also guards against out-of-range indices, which would otherwise reach `qk_obs_*` and
    // abort the process with a non-unwinding panic rather than raising.
    let (layout, num_qubits) = spec.mapping_layout(op.min_num_modes(), num_qubits)?;
    let encoding = spec.encoding();

    // The terms are mapped in parallel and merged under a memory bound; see `map_operator`, which is
    // shared with every other mapper emitting a `SparseObservable`. Each term is handed over with a
    // unit coefficient and scaled in the one `qk_obs_scaled_add_inplace` there, instead of costing an
    // extra observable per term.
    Ok(map_operator(op.iter(), num_qubits, |term, num_qubits| {
        (
            map_term(term.actions, term.modes, num_qubits, &layout, encoding),
            qk_coeff(term.coeff),
        )
    }))
}

/// Maps one term, whose coefficient is applied by the caller, onto a single observable.
///
/// All of the term's Pauli strings are built into one observable, so that `map_operator` can merge
/// the term's own duplicates before it is added to an accumulator.
fn map_term(
    actions: &[bool],
    modes: &[u32],
    num_qubits: u32,
    layout: &ModeLayout,
    encoding: LevelEncoding,
) -> *mut ffi::QkObs {
    // Bucket this term's ladder operators by mode, preserving their relative order *within* each
    // mode. Valid because operators on distinct modes commute; order within a mode is what
    // distinguishes `b b†` from `b† b`, so it must survive. A `BTreeMap` keeps the iteration
    // order deterministic (by mode index), which keeps the emitted term order reproducible.
    let mut per_mode_actions: BTreeMap<u32, Vec<bool>> = BTreeMap::new();
    for (&action, &mode) in actions.iter().zip(modes) {
        per_mode_actions.entry(mode).or_default().push(action);
    }

    // Start from a unit coefficient and compose one mode's contribution at a time.
    let mut acc: Vec<EncodedTerm> = vec![(Vec::new(), Complex64::new(1.0, 0.0))];

    for (&mode, mode_actions) in &per_mode_actions {
        // Each mode is truncated at its own cutoff and sits at its own offset with its own width.
        let mode = mode as usize;
        let elements = single_mode_matrix(mode_actions, layout.max_occupations[mode]);
        if elements.is_empty() {
            // The truncated product annihilates every level: the whole term is zero.
            acc.clear();
            break;
        }

        let offset = layout.offsets[mode];
        let width = layout.widths[mode];
        let mut next: Vec<EncodedTerm> = Vec::new();
        for element in &elements {
            let row_bits = encoding.level_to_bits(element.row);
            let col_bits = encoding.level_to_bits(element.col);
            let encoded = encode_matrix_element(
                row_bits,
                col_bits,
                offset,
                width,
                element.value,
                encoding.drops_idle_projectors(),
            );
            for (bits, coeff) in &encoded {
                for (acc_bits, acc_coeff) in &acc {
                    let mut merged = acc_bits.clone();
                    merged.extend_from_slice(bits);
                    next.push((merged, acc_coeff * coeff));
                }
            }
        }
        acc = next;
    }

    // Emit the accumulated strings as one multi-term observable. `SparseObservable` requires each
    // term's indices to be sorted ascending; modes were visited in ascending order and qubits within
    // a mode likewise, so the list is already sorted. Sort defensively anyway -- an unsorted index
    // list is rejected by `qk_obs_new` at runtime, which would be a much harder failure to diagnose
    // than a redundant sort.
    let mut coeffs: Vec<QkComplex64> = Vec::with_capacity(acc.len());
    let mut bit_terms: Vec<ffi::QkBitTerm> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut boundaries: Vec<usize> = vec![0];
    for (mut bits, coeff) in acc {
        if coeff == Complex64::new(0.0, 0.0) {
            continue;
        }
        bits.sort_by_key(|(qubit, _)| *qubit);
        coeffs.push(qk_coeff(coeff));
        bit_terms.extend(bits.iter().map(|(_, bit_term)| *bit_term));
        indices.extend(bits.iter().map(|(qubit, _)| *qubit));
        boundaries.push(bit_terms.len());
    }

    pauli_sum(num_qubits, &coeffs, &bit_terms, &indices, &boundaries)
}

#[cfg(test)]
mod tests {
    use super::*;

    const CRE: bool = true;
    const ANN: bool = false;

    /// Densifies a matrix-element list for easy comparison.
    fn dense(elements: &[MatrixElement], dim: usize) -> Vec<Vec<f64>> {
        let mut out = vec![vec![0.0; dim]; dim];
        for e in elements {
            out[e.row as usize][e.col as usize] = e.value.re;
        }
        out
    }

    /// `b†` has the sqrt(n+1) subdiagonal, and truncation drops the top creation.
    #[test]
    fn test_single_mode_matrix_creation() {
        // d = 4
        let m = dense(&single_mode_matrix(&[CRE], 3), 4);
        // b†|0> = |1>, b†|1> = sqrt2|2>, b†|2> = sqrt3|3>, b†|3> -> dropped
        assert!((m[1][0] - 1.0).abs() < 1e-12);
        assert!((m[2][1] - 2.0_f64.sqrt()).abs() < 1e-12);
        assert!((m[3][2] - 3.0_f64.sqrt()).abs() < 1e-12);
        // nothing in the last column: the amplitude left the truncated space
        assert!(m.iter().all(|row| row[3].abs() < 1e-12));
    }

    /// `b` has the sqrt(n) superdiagonal, and `b|0> = 0` exactly.
    #[test]
    fn test_single_mode_matrix_annihilation() {
        let m = dense(&single_mode_matrix(&[ANN], 3), 4);
        assert!((m[0][1] - 1.0).abs() < 1e-12);
        assert!((m[1][2] - 2.0_f64.sqrt()).abs() < 1e-12);
        assert!((m[2][3] - 3.0_f64.sqrt()).abs() < 1e-12);
        assert!(m.iter().all(|row| row[0].abs() < 1e-12), "b|0> must vanish");
    }

    /// The number operator `n = b† b` is diagonal with entries `0, 1, ..., d-1`.
    ///
    /// Note the operator order: `[CRE, ANN]` is the product `b† b`, whose rightmost factor `b`
    /// acts first.
    #[test]
    fn test_single_mode_matrix_number_operator() {
        let dim = 5usize;
        let m = dense(&single_mode_matrix(&[CRE, ANN], dim as u32 - 1), dim);
        for (row, values) in m.iter().enumerate() {
            for (col, &value) in values.iter().enumerate() {
                let expected = if row == col { row as f64 } else { 0.0 };
                assert!(
                    (value - expected).abs() < 1e-12,
                    "n[{row}][{col}] = {value} != {expected}"
                );
            }
        }
    }

    /// `b b†` is diagonal with entries `n+1` -- and, crucially, its top entry is *clipped* by the
    /// truncation, which is precisely the `[b, b†] = 1 - d|d-1><d-1|` defect.
    #[test]
    fn test_truncated_commutator_defect() {
        let dim = 4usize;
        let max_occ = dim as u32 - 1;
        // b b† : rightmost factor b† acts first
        let bbd = dense(&single_mode_matrix(&[ANN, CRE], max_occ), dim);
        // b† b
        let dbb = dense(&single_mode_matrix(&[CRE, ANN], max_occ), dim);

        // [b, b†] = b b† - b† b must equal diag(1, 1, ..., 1, 1 - d)
        for level in 0..dim {
            let commutator = bbd[level][level] - dbb[level][level];
            let expected = if level == dim - 1 {
                1.0 - dim as f64
            } else {
                1.0
            };
            assert!(
                (commutator - expected).abs() < 1e-12,
                "[b,b†][{level}][{level}] = {commutator} != {expected}"
            );
        }
    }

    /// A repeated creation `b† b†` populates level 2.
    #[test]
    fn test_single_mode_matrix_double_creation() {
        let m = dense(&single_mode_matrix(&[CRE, CRE], 3), 4);
        // b† b† |0> = sqrt(1)*sqrt(2) |2>
        assert!((m[2][0] - 2.0_f64.sqrt()).abs() < 1e-12);
        // b† b† |1> = sqrt(2)*sqrt(3) |3>
        assert!((m[3][1] - 6.0_f64.sqrt()).abs() < 1e-12);
        assert!(m[0][0].abs() < 1e-12);
    }

    /// An empty action list is the identity.
    #[test]
    fn test_single_mode_matrix_identity() {
        let m = dense(&single_mode_matrix(&[], 2), 3);
        for (row, values) in m.iter().enumerate() {
            for (col, &value) in values.iter().enumerate() {
                let expected = if row == col { 1.0 } else { 0.0 };
                assert!((value - expected).abs() < 1e-12);
            }
        }
    }

    /// Truncating at `max_occupation = 0` leaves only the vacuum, where every ladder operator dies.
    #[test]
    fn test_single_mode_matrix_trivial_truncation() {
        assert!(single_mode_matrix(&[CRE], 0).is_empty());
        assert!(single_mode_matrix(&[ANN], 0).is_empty());
        // ...but the identity survives
        assert_eq!(single_mode_matrix(&[], 0).len(), 1);
    }

    /// Composing in the matrix representation matches composing the factors' matrices by hand.
    ///
    /// This is what licenses the "multiply before encoding" strategy: the fused product is the same
    /// operator as the product of the factors, so nothing is lost by never materializing the
    /// factors separately.
    #[test]
    fn test_matrix_composition_matches_manual_product() {
        let dim = 5usize;
        let max_occ = dim as u32 - 1;

        // n^2 = b† b b† b, built in one shot
        let fused = dense(&single_mode_matrix(&[CRE, ANN, CRE, ANN], max_occ), dim);

        // ...versus the square of the number-operator matrix
        let n = dense(&single_mode_matrix(&[CRE, ANN], max_occ), dim);
        let mut manual = vec![vec![0.0; dim]; dim];
        for (i, out_row) in manual.iter_mut().enumerate() {
            for (j, out) in out_row.iter_mut().enumerate() {
                *out = n[i].iter().zip(n.iter()).map(|(a, row)| a * row[j]).sum();
            }
        }

        for (i, (fused_row, manual_row)) in fused.iter().zip(manual.iter()).enumerate() {
            for (j, (&f, &mm)) in fused_row.iter().zip(manual_row.iter()).enumerate() {
                assert!(
                    (f - mm).abs() < 1e-12,
                    "n^2[{i}][{j}]: fused {f} != manual {mm}"
                );
            }
        }
    }
}
