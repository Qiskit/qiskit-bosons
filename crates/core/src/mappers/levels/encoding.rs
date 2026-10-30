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

//! Layouts of a `d`-level system on qubits.
//!
//! A [`LevelEncoding`] decides which computational basis state of a block of qubits represents each
//! of the `d` levels of a single system, and nothing else. It is deliberately independent of what the
//! levels *are*: occupation numbers of a truncated bosonic mode, vibrational modals, points of a
//! position grid and the `2s + 1` states of a spin-`s` particle all need the same choice to be made.
//! The mappers in [`library`](crate::mappers::library) decide which operator acts on those levels and
//! build the matrix elements; [`matrix`](crate::mappers::levels::matrix) then turns each matrix
//! element into qubit operators using the bit patterns the encoding assigns.
//!
//! Sawaya et al. (npj Quantum Information 6, 49 (2020), arXiv:1909.12847) compare exactly the
//! options implemented here for `d`-level systems. Their conclusion is explicitly *not* that one
//! encoding dominates:
//!
//! > although the unary code tends to require more qubits but fewer operations, it is often the case
//! > that the SB or Gray code is more efficient both in terms of qubit counts and operation counts.
//!
//! Concretely: at `d = 4` the unary code needs ~3x *more* two-qubit gates than the compact codes,
//! while by `d = 10` it needs about half as many. They further find that the Gray code suits
//! tridiagonal (hopping) terms while standard binary suits diagonal (number-operator) terms -- so
//! the optimum can even differ between terms of a single Hamiltonian.
//!
//! This module therefore offers them all as first-class choices rather than blessing one. The unary
//! code comes in two variants that share the same bit patterns and differ in the operator that is
//! built from them. See [`LevelEncoding`] for the per-variant trade-offs.

/// The bitmask that identifies a level: bit `q` is the value of qubit `q` of the system.
///
/// This is the one place the width of a bit pattern is chosen. The widths that follow from it,
/// [`MAX_UNARY_LEVELS`] and the number of qubits `encode_matrix_element` accepts, are derived from
/// [`LevelBits::BITS`], so changing the type changes them with it.
///
/// [`LevelBits::BITS`]: u128::BITS
pub type LevelBits = u128;

/// The largest number of levels the unary encoding can hold.
///
/// The unary encoding needs one bit of a [`LevelBits`] per level. Everything that has to agree about
/// this limit takes it from here: [`LevelEncoding::level_to_bits`] asserts it, and a
/// [`LevelSpec`](crate::mappers::levels::spec::LevelSpec) refuses to be built beyond it, so that a
/// level count which does not fit is reported as an error rather than reaching a shift that
/// overflows. The compact encodings need `ceil(log2 d)` bits and are not limited by it.
pub const MAX_UNARY_LEVELS: u32 = LevelBits::BITS;

/// The layout of the `d` levels of a single system on qubits.
///
/// Every variant represents `d` levels, numbered `0..d`. They differ in how the levels are laid out
/// across qubits, which trades qubit count against Pauli weight (and hence two-qubit gate count).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum LevelEncoding {
    /// One qubit per level (aka "unary", "one-hot", or "direct"), with a projector on every idle
    /// qubit.
    ///
    /// Level `n` is the computational basis state with qubit `n` set and all others clear, so a
    /// system costs `d` qubits. A ladder-like operator moving the level up or down by one is, as a
    /// sum of terms with two non-diagonal factors,
    ///
    /// ```text
    /// sum_{n=0}^{d-2} c_n |1><0|_{n+1} (x) |0><1|_n ,
    /// ```
    ///
    /// each moving the single excitation up by one level. This variant keeps the one-hot constraint
    /// in the operator: the matrix-element encoder in [`levels::matrix`](crate::mappers::levels::matrix)
    /// puts a projector `|0><0|` on each of the other qubits of the system, so the terms it emits
    /// touch all `d` qubits and vanish outside the one-hot subspace. That makes them costlier to
    /// synthesize than the weight-2 form of [`LevelEncoding::UnaryUnprojected`], by a margin that
    /// grows quickly with `d`.
    ///
    /// The two unary variants share their bit patterns and agree on the one-hot subspace. The
    /// trade-off of the encoding is the qubit cost and the fact that only `d` of the `2^d` basis
    /// states are physical: there is a large unphysical subspace, which a simulation must avoid
    /// leaking into (see [`LevelEncoding::num_qubits`]).
    UnaryProjected,
    /// One qubit per level, without the projectors on the idle qubits.
    ///
    /// The bit patterns are those of [`LevelEncoding::UnaryProjected`], and so is the matrix element
    /// each operator term is built from. This variant drops the `|0><0|` projector the encoder would
    /// put on each idle qubit, so the ladder terms keep the weight-2 form
    ///
    /// ```text
    /// sum_{n=0}^{d-2} c_n |1><0|_{n+1} (x) |0><1|_n ,
    /// ```
    ///
    /// which acts as the identity on every other qubit. Because the Pauli weight of this form is 2
    /// regardless of `d`, the gate count grows only linearly in `d` (Sawaya et al. bound a
    /// tridiagonal operator by `O(d)` CNOTs), which is what makes the unary code attractive at large
    /// `d`.
    ///
    /// The operator is equal to that of [`LevelEncoding::UnaryProjected`] on the one-hot subspace
    /// and is **undefined outside it**: it does not suppress the evolution of an unphysical state.
    /// Choosing it is a modeling decision and not a pure optimization. It is only meaningful for
    /// the unary code. Dropping the projectors of a compact encoding changes the operator on the
    /// physical states, because there the zero bits are part of the level number.
    UnaryUnprojected,
    /// Standard binary: level `n` is stored as the binary representation of `n`.
    ///
    /// A system costs only `ceil(log2 d)` qubits, an exponential saving over the unary encodings, at
    /// the price of higher-weight Pauli strings: adjacent levels can
    /// differ in up to `log2 d` bits, and every differing bit position is forced to carry a
    /// non-identity Pauli. Sawaya et al. bound a tridiagonal operator by `O(d^2 log d)` CNOTs here.
    ///
    /// They find standard binary preferable for *diagonal* (number-operator) terms, where the
    /// bit-pattern distance never enters.
    StandardBinary,
    /// Gray code: level `n` is stored as `n ^ (n >> 1)`.
    ///
    /// Costs the same `ceil(log2 d)` qubits as [`LevelEncoding::StandardBinary`] and (perhaps
    /// surprisingly) yields exactly the same *number* of Pauli terms. The advantage is entirely in
    /// the weight *distribution*: consecutive levels differ in exactly one bit (Hamming distance 1)
    /// by construction, so the hopping terms that dominate a tridiagonal Hamiltonian produce
    /// low-weight strings where standard binary produces maximal-weight ones.
    ///
    /// Sawaya et al. accordingly find the Gray code preferable for *tridiagonal* (hopping) terms.
    /// Note that this advantage is specific to couplings between adjacent levels; it does not
    /// extend to e.g. a pentadiagonal structure.
    Gray,
}

impl LevelEncoding {
    /// Returns whether this is one of the unary encodings, which use one qubit per level.
    ///
    /// This is the one place that decides what is subject to the unary limit
    /// [`MAX_UNARY_LEVELS`], and which width formula applies.
    pub fn is_unary(&self) -> bool {
        matches!(
            self,
            LevelEncoding::UnaryProjected | LevelEncoding::UnaryUnprojected
        )
    }

    /// Returns whether the operator built from this encoding omits the `|0><0|` projector on the
    /// idle qubits of a system.
    ///
    /// Only [`LevelEncoding::UnaryUnprojected`] does. The compact encodings never can, because
    /// their zero bits are part of the level number.
    pub fn drops_idle_projectors(&self) -> bool {
        matches!(self, LevelEncoding::UnaryUnprojected)
    }

    /// Returns the number of qubits needed to represent a system with `num_levels` levels.
    ///
    /// Note that the compact encodings can address `2^ceil(log2 num_levels)` levels, which exceeds
    /// `num_levels` whenever it is not a power of two. Those surplus levels are unphysical and are
    /// never populated by the operators built from the encoding.
    pub fn num_qubits(&self, num_levels: u32) -> u32 {
        if self.is_unary() {
            num_levels
        } else {
            num_levels.next_power_of_two().trailing_zeros()
        }
    }

    /// Returns the qubit basis state (as a bitmask) representing level `level`.
    ///
    /// The bit at position `q` of the returned mask is the value of the system's qubit `q`.
    ///
    /// # Panics
    ///
    /// Panics for the unary encodings if `level` is [`MAX_UNARY_LEVELS`] or more, as the mask then
    /// has no bit for it. A [`LevelSpec`](crate::mappers::levels::spec::LevelSpec) never lets such
    /// a level reach this method. The compact encodings cannot overflow the mask.
    pub fn level_to_bits(&self, level: u32) -> LevelBits {
        match self {
            LevelEncoding::UnaryProjected | LevelEncoding::UnaryUnprojected => {
                assert!(
                    level < MAX_UNARY_LEVELS,
                    "the unary encoding holds at most {MAX_UNARY_LEVELS} levels, but level {level} \
                     was requested"
                );
                LevelBits::from(1u8) << level
            }
            LevelEncoding::StandardBinary => LevelBits::from(level),
            // The standard binary-reflected Gray code.
            LevelEncoding::Gray => LevelBits::from(level ^ (level >> 1)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNARY: [LevelEncoding; 2] = [
        LevelEncoding::UnaryProjected,
        LevelEncoding::UnaryUnprojected,
    ];

    #[test]
    fn test_helpers_single_out_the_unary_variants() {
        for enc in UNARY {
            assert!(enc.is_unary());
        }
        assert!(!LevelEncoding::StandardBinary.is_unary());
        assert!(!LevelEncoding::Gray.is_unary());

        // only the unprojected unary form drops the idle projectors
        assert!(LevelEncoding::UnaryUnprojected.drops_idle_projectors());
        assert!(!LevelEncoding::UnaryProjected.drops_idle_projectors());
        assert!(!LevelEncoding::StandardBinary.drops_idle_projectors());
        assert!(!LevelEncoding::Gray.drops_idle_projectors());
    }

    #[test]
    fn test_num_qubits() {
        // unary: one qubit per level, for both variants
        for enc in UNARY {
            assert_eq!(enc.num_qubits(2), 2);
            assert_eq!(enc.num_qubits(4), 4);
            assert_eq!(enc.num_qubits(8), 8);
        }

        // compact: ceil(log2 d)
        for enc in [LevelEncoding::StandardBinary, LevelEncoding::Gray] {
            assert_eq!(enc.num_qubits(2), 1); // d=2
            assert_eq!(enc.num_qubits(3), 2); // d=3 -> 2 qubits
            assert_eq!(enc.num_qubits(4), 2); // d=4
            assert_eq!(enc.num_qubits(5), 3); // d=5 -> 3 qubits
            assert_eq!(enc.num_qubits(8), 3); // d=8
            assert_eq!(enc.num_qubits(9), 4); // d=9 -> 4 qubits
        }
    }

    #[test]
    fn test_level_to_bits_unary_is_one_hot() {
        for enc in UNARY {
            for level in 0..5 {
                let bits = enc.level_to_bits(level);
                assert_eq!(bits.count_ones(), 1, "unary must set exactly one bit");
                assert_eq!(bits.trailing_zeros(), level);
            }
        }
    }

    /// Both unary variants assign the same bit patterns; they differ only in the operator.
    #[test]
    fn test_the_unary_variants_share_their_bit_patterns() {
        for level in 0..MAX_UNARY_LEVELS {
            assert_eq!(
                LevelEncoding::UnaryProjected.level_to_bits(level),
                LevelEncoding::UnaryUnprojected.level_to_bits(level)
            );
        }
    }

    /// Unary holds every level up to the limit, including the ones a `u32` mask could not.
    #[test]
    fn test_level_to_bits_unary_reaches_the_limit() {
        for enc in UNARY {
            for level in [31, 32, 33, 64, MAX_UNARY_LEVELS - 1] {
                let bits = enc.level_to_bits(level);
                assert_eq!(bits, LevelBits::from(1u8) << level);
                assert_eq!(bits.trailing_zeros(), level);
            }
        }
    }

    #[test]
    #[should_panic(expected = "the unary encoding holds at most")]
    fn test_level_to_bits_projected_rejects_beyond_the_limit() {
        LevelEncoding::UnaryProjected.level_to_bits(MAX_UNARY_LEVELS);
    }

    #[test]
    #[should_panic(expected = "the unary encoding holds at most")]
    fn test_level_to_bits_unprojected_rejects_beyond_the_limit() {
        LevelEncoding::UnaryUnprojected.level_to_bits(MAX_UNARY_LEVELS);
    }

    #[test]
    fn test_level_to_bits_standard_binary_is_identity() {
        for level in 0..16 {
            assert_eq!(
                LevelEncoding::StandardBinary.level_to_bits(level),
                LevelBits::from(level)
            );
        }
    }

    /// The compact encodings are not limited by the width of the mask, even at the largest level.
    #[test]
    fn test_level_to_bits_compact_handles_the_largest_level() {
        let level = (1u32 << 31) - 1;
        assert_eq!(
            LevelEncoding::StandardBinary.level_to_bits(level),
            LevelBits::from(level)
        );
        assert_eq!(
            LevelEncoding::Gray.level_to_bits(level),
            LevelBits::from(level ^ (level >> 1))
        );
    }

    /// The defining property of the Gray code: adjacent levels differ in exactly one bit.
    ///
    /// This is the whole reason the Gray code helps for hopping terms, so it is asserted directly
    /// rather than by spot-checking a table.
    #[test]
    fn test_gray_code_adjacent_hamming_distance_is_one() {
        for level in 0..31u32 {
            let a = LevelEncoding::Gray.level_to_bits(level);
            let b = LevelEncoding::Gray.level_to_bits(level + 1);
            assert_eq!(
                (a ^ b).count_ones(),
                1,
                "Gray levels {level} and {} must differ in exactly 1 bit",
                level + 1
            );
        }
        // ...and it is a genuine permutation of the levels (injective).
        let codes: Vec<LevelBits> = (0..16)
            .map(|l| LevelEncoding::Gray.level_to_bits(l))
            .collect();
        let mut sorted = codes.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), codes.len());
    }

    /// Standard binary reaches the maximal Hamming distance that the Gray code avoids.
    #[test]
    fn test_standard_binary_has_large_adjacent_hamming_distance() {
        // 3 -> 4 flips 3 bits in standard binary (011 -> 100) but only 1 in Gray.
        let sb_a = LevelEncoding::StandardBinary.level_to_bits(3);
        let sb_b = LevelEncoding::StandardBinary.level_to_bits(4);
        assert_eq!((sb_a ^ sb_b).count_ones(), 3);

        let g_a = LevelEncoding::Gray.level_to_bits(3);
        let g_b = LevelEncoding::Gray.level_to_bits(4);
        assert_eq!((g_a ^ g_b).count_ones(), 1);
    }
}
