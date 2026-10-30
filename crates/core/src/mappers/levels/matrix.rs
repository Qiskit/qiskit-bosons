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

//! Encoding the matrix elements of an operator on a `d`-level system as qubit operators.
//!
//! A mapper decides which operator acts on the levels of a system and expresses it as a list of
//! [`MatrixElement`]s. This module turns one such element into Pauli and projector terms, given the
//! bit patterns that a [`LevelEncoding`](crate::mappers::levels::encoding::LevelEncoding) assigns to
//! its row and column. Nothing here depends on what the levels are, so a mapper for any level basis
//! can reuse it as it is.
//!
//! Each matrix element `<row| O |col>` becomes one term. In the qubit basis, moving amplitude from
//! the bit pattern of `col` to that of `row` is a product of single-qubit transfer operators, one per
//! qubit of the system:
//!
//! | `row` bit | `col` bit | single-qubit factor                    |
//! |-----------|-----------|----------------------------------------|
//! | 0         | 0         | `|0><0|` (Qiskit's `Zero` bit term)    |
//! | 1         | 1         | `|1><1|` (Qiskit's `One` bit term)     |
//! | 1         | 0         | `|1><0| = (X - iY)/2`                  |
//! | 0         | 1         | `|0><1| = (X + iY)/2`                  |
//!
//! Qiskit's `SparseObservable` carries the projector bit terms natively, so the diagonal cases need
//! no Pauli expansion at all; only the off-diagonal ones are expanded into their two Pauli terms.

use crate::mappers::levels::encoding::LevelBits;
use crate::mappers::qk_obs::{
    QK_BIT_TERM_ONE, QK_BIT_TERM_X, QK_BIT_TERM_Y, QK_BIT_TERM_ZERO, ffi_backend as ffi,
};
use num_complex::Complex64;

/// A single matrix element `<row| O |col>` of an operator on a `d`-level system.
///
/// The encodings map *matrix elements*, so this is the currency in which single-system operators are
/// expressed before being encoded. It carries no assumption about what the levels are, so a mapper
/// for any level basis can produce these.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MatrixElement {
    /// The row index (the level the element maps *to*).
    pub row: u32,
    /// The column index (the level the element maps *from*).
    pub col: u32,
    /// The value of the matrix element.
    pub value: Complex64,
}

/// A single-qubit factor of an encoded matrix element, before Pauli expansion.
///
/// The two diagonal cases are single `SparseObservable` bit terms; the two off-diagonal ones each
/// expand into a sum of two Pauli terms, which is what forces the cartesian product in
/// [`encode_matrix_element`].
#[derive(Clone, Copy, PartialEq)]
enum QubitFactor {
    /// `|0><0|`
    Zero,
    /// `|1><1|`
    One,
    /// `|1><0|` -- sets the qubit, i.e. raises this digit of the level.
    Raise,
    /// `|0><1|` -- clears the qubit, i.e. lowers this digit of the level.
    Lower,
}

impl QubitFactor {
    /// Returns the Pauli expansion of this factor as `(bit_term, coefficient)` pairs.
    ///
    /// The diagonal projectors are returned as a single native bit term with coefficient 1, so no
    /// expansion happens for them.
    ///
    /// # Sign convention
    ///
    /// In Qiskit's convention `(X + iY) / 2 = |0><1|` and `(X - iY) / 2 = |1><0|`. Note that these
    /// are the *opposite* way round from the `sigma^± = (X ± iY) / 2` shorthand commonly written in
    /// the bosonic-encoding literature, where `sigma^+` denotes the raising operator `|1><0|`.
    /// Getting this backwards silently transposes the whole mapped operator -- it maps `b` to `b†` --
    /// which is why the assignment is spelled out here and asserted directly in the tests.
    fn expansion(self) -> Vec<(ffi::QkBitTerm, Complex64)> {
        match self {
            QubitFactor::Zero => vec![(QK_BIT_TERM_ZERO, Complex64::new(1.0, 0.0))],
            QubitFactor::One => vec![(QK_BIT_TERM_ONE, Complex64::new(1.0, 0.0))],
            // |1><0| = (X - iY)/2
            QubitFactor::Raise => vec![
                (QK_BIT_TERM_X, Complex64::new(0.5, 0.0)),
                (QK_BIT_TERM_Y, Complex64::new(0.0, -0.5)),
            ],
            // |0><1| = (X + iY)/2
            QubitFactor::Lower => vec![
                (QK_BIT_TERM_X, Complex64::new(0.5, 0.0)),
                (QK_BIT_TERM_Y, Complex64::new(0.0, 0.5)),
            ],
        }
    }
}

/// Returns the single-qubit factor moving bit `col_bit` to bit `row_bit`.
fn qubit_factor(row_bit: bool, col_bit: bool) -> QubitFactor {
    match (row_bit, col_bit) {
        (false, false) => QubitFactor::Zero,
        (true, true) => QubitFactor::One,
        (true, false) => QubitFactor::Raise,
        (false, true) => QubitFactor::Lower,
    }
}

/// One fully-expanded Pauli/projector term: `(qubit, bit_term)` pairs plus a coefficient.
pub type EncodedTerm = (Vec<(u32, ffi::QkBitTerm)>, Complex64);

/// Expands one matrix element of one system into `SparseObservable` terms.
///
/// `qubit_offset` is the index of this system's first qubit in the global register, and `num_qubits`
/// the number of qubits this system occupies. `row_bits`/`col_bits` are the encoded bit patterns.
///
/// Nothing here depends on what the levels are: it consumes the bit patterns a
/// [`LevelEncoding`](crate::mappers::levels::encoding::LevelEncoding) assigns and the matrix
/// element, so a mapper in any level basis can reuse it as is.
///
/// The result is the cartesian product over the system's qubits of each factor's Pauli expansion. A
/// factor that is a native projector contributes a single branch, so an all-diagonal element yields
/// exactly one term, while an element flipping `k` qubits yields `2^k`.
///
/// With `drop_zero_projectors` the `|0><0|` factor of every qubit that is `0` in both patterns is
/// left out, so that qubit carries the identity instead. The `|1><1|` factors and the factors of
/// the qubits that flip are always kept. This is only meaningful when the zero bits are not part of
/// the level number, as in the unary encodings: see
/// [`LevelEncoding::drops_idle_projectors`](crate::mappers::levels::encoding::LevelEncoding::drops_idle_projectors).
/// The terms with and without it agree on every state in which each idle qubit is `0`.
///
/// # Panics
///
/// Panics if `num_qubits` exceeds the width of a [`LevelBits`] pattern: the patterns have no bit for
/// the qubits beyond it. The patterns a
/// [`LevelEncoding`](crate::mappers::levels::encoding::LevelEncoding) assigns never need more.
pub fn encode_matrix_element(
    row_bits: LevelBits,
    col_bits: LevelBits,
    qubit_offset: u32,
    num_qubits: u32,
    value: Complex64,
    drop_zero_projectors: bool,
) -> Vec<EncodedTerm> {
    assert!(
        num_qubits <= LevelBits::BITS,
        "a bit pattern holds at most {} qubits, but {num_qubits} were requested",
        LevelBits::BITS
    );
    let mut terms: Vec<EncodedTerm> = vec![(Vec::new(), value)];

    for local in 0..num_qubits {
        let row_bit = (row_bits >> local) & 1 == 1;
        let col_bit = (col_bits >> local) & 1 == 1;
        let factor = qubit_factor(row_bit, col_bit);
        // Only the `|0><0|` factor can be left out: it is the one that says "this qubit is idle".
        if drop_zero_projectors && factor == QubitFactor::Zero {
            continue;
        }
        let expansion = factor.expansion();

        let mut next: Vec<EncodedTerm> = Vec::with_capacity(terms.len() * expansion.len());
        for (bits, coeff) in &terms {
            for (bit_term, factor_coeff) in &expansion {
                let mut extended = bits.clone();
                extended.push((qubit_offset + local, *bit_term));
                next.push((extended, coeff * factor_coeff));
            }
        }
        terms = next;
    }

    terms
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Densifies a single-qubit factor's Pauli expansion into its 2x2 matrix.
    ///
    /// Uses Qiskit's Pauli conventions: `X = [[0,1],[1,0]]`, `Y = [[0,-i],[i,0]]`, and the
    /// projectors `Zero = |0><0|`, `One = |1><1|`. Index `[r][c]` is `<r| M |c>`.
    fn densify(factor: QubitFactor) -> [[Complex64; 2]; 2] {
        let zero = Complex64::new(0.0, 0.0);
        let mut out = [[zero; 2], [zero; 2]];
        for (bit_term, coeff) in factor.expansion() {
            let m = if bit_term == QK_BIT_TERM_X {
                [
                    [zero, Complex64::new(1.0, 0.0)],
                    [Complex64::new(1.0, 0.0), zero],
                ]
            } else if bit_term == QK_BIT_TERM_Y {
                [
                    [zero, Complex64::new(0.0, -1.0)],
                    [Complex64::new(0.0, 1.0), zero],
                ]
            } else if bit_term == QK_BIT_TERM_ZERO {
                [[Complex64::new(1.0, 0.0), zero], [zero, zero]]
            } else if bit_term == QK_BIT_TERM_ONE {
                [[zero, zero], [zero, Complex64::new(1.0, 0.0)]]
            } else {
                panic!("unexpected bit term in expansion")
            };
            for r in 0..2 {
                for c in 0..2 {
                    out[r][c] += coeff * m[r][c];
                }
            }
        }
        out
    }

    fn assert_matrix(actual: [[Complex64; 2]; 2], expected: [[f64; 2]; 2], label: &str) {
        for r in 0..2 {
            for c in 0..2 {
                let got = actual[r][c];
                assert!(
                    (got.re - expected[r][c]).abs() < 1e-12 && got.im.abs() < 1e-12,
                    "{label}[{r}][{c}] = {got} != {}",
                    expected[r][c]
                );
            }
        }
    }

    /// Pins the raise/lower Pauli convention.
    ///
    /// This guards a bug that is invisible in any Hermitian test case: swapping the two expansions
    /// transposes every mapped operator, silently turning `b` into `b†`. The trap is that the
    /// bosonic-encoding literature writes `sigma^+ = (X + iY)/2` for the *raising* operator, whereas
    /// in Qiskit's matrix convention `(X + iY)/2` is `|0><1|`, which *lowers*. The number operator
    /// `n = b† b` is symmetric and passes either way, so only an asymmetric operator catches it.
    #[test]
    fn test_raise_lower_pauli_convention() {
        // |1><0| : maps |0> to |1>, so the nonzero entry is row 1, column 0.
        assert_matrix(
            densify(QubitFactor::Raise),
            [[0.0, 0.0], [1.0, 0.0]],
            "Raise",
        );
        // |0><1| : maps |1> to |0>, so the nonzero entry is row 0, column 1.
        assert_matrix(
            densify(QubitFactor::Lower),
            [[0.0, 1.0], [0.0, 0.0]],
            "Lower",
        );
        // and the diagonal projectors
        assert_matrix(densify(QubitFactor::Zero), [[1.0, 0.0], [0.0, 0.0]], "Zero");
        assert_matrix(densify(QubitFactor::One), [[0.0, 0.0], [0.0, 1.0]], "One");
    }

    /// `qubit_factor` selects the factor that carries `col_bit` to `row_bit`.
    #[test]
    fn test_qubit_factor_selection() {
        assert!(qubit_factor(false, false) == QubitFactor::Zero);
        assert!(qubit_factor(true, true) == QubitFactor::One);
        // row=1, col=0 -> must RAISE the qubit
        assert!(qubit_factor(true, false) == QubitFactor::Raise);
        // row=0, col=1 -> must LOWER it
        assert!(qubit_factor(false, true) == QubitFactor::Lower);
    }

    /// A diagonal element expands to exactly one term; each flipped qubit doubles the count.
    #[test]
    fn test_encode_matrix_element_term_counts() {
        let one = Complex64::new(1.0, 0.0);
        // 2 qubits, no flips (row == col): a single projector-only term
        assert_eq!(encode_matrix_element(0b10, 0b10, 0, 2, one, false).len(), 1);
        // one flipped qubit -> 2 Pauli terms
        assert_eq!(encode_matrix_element(0b11, 0b10, 0, 2, one, false).len(), 2);
        // two flipped qubits -> 4
        assert_eq!(encode_matrix_element(0b01, 0b10, 0, 2, one, false).len(), 4);
        // three flipped qubits -> 8
        assert_eq!(
            encode_matrix_element(0b111, 0b000, 0, 3, one, false).len(),
            8
        );
    }

    /// The qubit offset places a system's terms on the right global qubits.
    #[test]
    fn test_encode_matrix_element_respects_offset() {
        let one = Complex64::new(1.0, 0.0);
        let terms = encode_matrix_element(0b01, 0b01, 5, 2, one, false);
        assert_eq!(terms.len(), 1);
        let indices: Vec<u32> = terms[0].0.iter().map(|(q, _)| *q).collect();
        assert_eq!(indices, vec![5, 6]);
    }

    /// Collects the bit terms of every encoded term, to look at which factors are present.
    fn bit_terms(terms: &[EncodedTerm]) -> Vec<Vec<(u32, ffi::QkBitTerm)>> {
        terms.iter().map(|(bits, _)| bits.clone()).collect()
    }

    /// Dropping the zero projectors leaves out exactly the `|0><0|` factors and nothing else: no
    /// `|0><0|` is left, every other factor is kept, and the element has strictly fewer factors.
    #[test]
    fn test_dropping_zero_projectors_only_removes_zero_factors() {
        let one = Complex64::new(1.0, 0.0);
        // a unary-style ladder element on 4 qubits: level 1 -> level 2 flips qubits 1 and 2
        // (patterns 0b0010 -> 0b0100), qubits 0 and 3 are idle and read 0 in both
        let (row, col) = (0b0100, 0b0010);
        let kept = bit_terms(&encode_matrix_element(row, col, 0, 4, one, false));
        let dropped = bit_terms(&encode_matrix_element(row, col, 0, 4, one, true));

        assert_eq!(
            kept.len(),
            dropped.len(),
            "the Pauli branching is unchanged"
        );
        for (full, slim) in kept.iter().zip(&dropped) {
            assert!(slim.len() < full.len(), "strictly fewer factors");
            assert_eq!(full.len(), 4);
            assert_eq!(slim.len(), 2, "only the two flipped qubits remain");
            assert!(slim.iter().all(|(_, term)| *term != QK_BIT_TERM_ZERO));
            // the surviving factors are exactly the full term's non-zero-projector factors
            let expected: Vec<_> = full
                .iter()
                .filter(|(_, term)| *term != QK_BIT_TERM_ZERO)
                .cloned()
                .collect();
            assert_eq!(slim, &expected);
        }
        // and the full term really did carry zero projectors on the idle qubits 0 and 3
        assert!(kept.iter().all(|t| {
            t.iter()
                .filter(|(_, term)| *term == QK_BIT_TERM_ZERO)
                .count()
                == 2
        }));
    }

    /// The `|1><1|` projector is never dropped, even when everything else is.
    #[test]
    fn test_dropping_zero_projectors_never_drops_one_projectors() {
        let one = Complex64::new(1.0, 0.0);
        // a diagonal element on level 2 of 4 qubits (pattern 0b0100): one `|1><1|`, three `|0><0|`
        let kept = bit_terms(&encode_matrix_element(0b0100, 0b0100, 0, 4, one, false));
        let dropped = bit_terms(&encode_matrix_element(0b0100, 0b0100, 0, 4, one, true));
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].len(), 4);
        assert_eq!(dropped, vec![vec![(2, QK_BIT_TERM_ONE)]]);
    }

    /// The vacuum projector `|0..0><0..0|` becomes the identity, which is an empty term.
    #[test]
    fn test_dropping_zero_projectors_of_the_vacuum_leaves_the_identity() {
        let one = Complex64::new(1.0, 0.0);
        let dropped = encode_matrix_element(0, 0, 0, 3, one, true);
        assert_eq!(dropped.len(), 1);
        assert!(dropped[0].0.is_empty());
        assert_eq!(dropped[0].1, one);
    }

    /// Without the flag nothing changes, and the flag has no effect on an element without idle zeros.
    #[test]
    fn test_the_flag_is_a_no_op_without_idle_zero_qubits() {
        let one = Complex64::new(1.0, 0.0);
        // every qubit flips: there is no `|0><0|` factor to drop
        let kept = bit_terms(&encode_matrix_element(0b11, 0b00, 0, 2, one, false));
        let dropped = bit_terms(&encode_matrix_element(0b11, 0b00, 0, 2, one, true));
        assert_eq!(kept, dropped);
    }
}
