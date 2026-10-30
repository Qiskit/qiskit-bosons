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

//! The truncation of a register of modes, together with the encoding that lays it out on qubits.
//!
//! A [`LevelSpec`] is the one place that knows how many levels each mode keeps and how many qubits
//! each mode therefore occupies. Every consumer (the mapper, the transpiler layout, the state
//! preparation) derives its qubit offsets from it, so they cannot disagree about where a mode
//! lives.
//!
//! The cutoff of a mode is its `max_occupation`, so the mode keeps `d = max_occupation + 1` levels.
//! A spec either holds a single cutoff shared by every mode, in which case the number of modes is
//! left open and fixed later by whatever the spec is applied to, or one cutoff per mode, which also
//! fixes the number of modes.
//!
//! # Mode widths
//!
//! The width of a mode is [`LevelEncoding::num_qubits`] of its level count. For the compact
//! encodings a mode with a single level (`max_occupation = 0`) has width zero: its only level needs
//! no qubit to be told apart. Such a mode is allowed and simply does not advance the qubit offset,
//! but a spec whose modes all have width zero would span no qubits at all and is rejected.

use crate::mappers::levels::encoding::{LevelEncoding, MAX_UNARY_LEVELS};
use thiserror::Error;

/// The largest cutoff a spec accepts, so that the level count `max_occupation + 1` stays within the
/// range for which [`LevelEncoding::num_qubits`] is defined.
///
/// The unary encodings are limited further, to [`MAX_UNARY_LEVELS`] levels; see [`check_cutoff`].
const MAX_CUTOFF: u32 = (1 << 31) - 1;

/// The ways in which building or applying a [`LevelSpec`] can fail.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum LevelSpecError {
    #[error(
        "the spec would span no qubits: under a compact encoding a mode with a single level \
         (max_occupation = 0) occupies zero qubits, so at least one mode needs max_occupation >= 1 \
         (or use a unary encoding)"
    )]
    NoQubits,
    #[error("max_occupation = {0} is too large; it must not exceed {MAX_CUTOFF}")]
    CutoffTooLarge(u32),
    #[error(
        "max_occupation = {0} is too large for the unary encodings, which hold at most \
         {MAX_UNARY_LEVELS} levels per mode (max_occupation = {max}); use a compact encoding \
         for more levels",
        max = MAX_UNARY_LEVELS - 1
    )]
    UnaryCutoffTooLarge(u32),
    #[error("the spec covers {expected} modes, but {got} were expected")]
    ModeCountMismatch { expected: usize, got: usize },
    #[error("the operator acts on {needed} modes, but the spec covers only {available}")]
    TooFewModes { needed: u32, available: u32 },
    #[error("{num_qubits} qubits are too few: the spec needs at least {needed}")]
    NumQubitsTooSmall { num_qubits: u32, needed: u32 },
    #[error("the register would need more than u32::MAX qubits")]
    TooManyQubits,
}

/// The cutoffs of a register: one shared value, or one value per mode.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Cutoffs {
    Uniform(u32),
    PerMode(Vec<u32>),
}

/// The truncation of every mode of a register and the encoding that lays the levels out on qubits.
///
/// See the [module documentation](self) for what a spec is.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct LevelSpec {
    cutoffs: Cutoffs,
    encoding: LevelEncoding,
}

/// A [`LevelSpec`] applied to a definite number of modes.
///
/// All vectors hold one entry per mode. Mode `m` occupies the qubits `offsets[m]..offsets[m] +
/// widths[m]`, and `offsets` is the running sum of `widths`, so a mode of width zero does not move
/// the next one.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModeLayout {
    /// The cutoff of each mode.
    pub max_occupations: Vec<u32>,
    /// The number of qubits of each mode.
    pub widths: Vec<u32>,
    /// The index of the first qubit of each mode.
    pub offsets: Vec<u32>,
    /// The number of qubits all modes together occupy.
    pub total_width: u32,
}

impl LevelSpec {
    /// Builds a spec in which every mode has the cutoff `max_occupation`.
    ///
    /// The number of modes is left open. Fails if the cutoff is too large, or if the encoding would
    /// give every mode width zero.
    pub fn uniform(max_occupation: u32, encoding: LevelEncoding) -> Result<Self, LevelSpecError> {
        check_cutoff(max_occupation, encoding)?;
        if width(encoding, max_occupation) == 0 {
            return Err(LevelSpecError::NoQubits);
        }
        Ok(Self {
            cutoffs: Cutoffs::Uniform(max_occupation),
            encoding,
        })
    }

    /// Builds a spec with one cutoff per mode, which also fixes the number of modes.
    ///
    /// Fails if a cutoff is too large, or if no mode would occupy a qubit (which includes an empty
    /// list of cutoffs).
    pub fn per_mode(
        max_occupations: Vec<u32>,
        encoding: LevelEncoding,
    ) -> Result<Self, LevelSpecError> {
        for &max_occupation in &max_occupations {
            check_cutoff(max_occupation, encoding)?;
        }
        if max_occupations
            .iter()
            .all(|&max_occupation| width(encoding, max_occupation) == 0)
        {
            return Err(LevelSpecError::NoQubits);
        }
        Ok(Self {
            cutoffs: Cutoffs::PerMode(max_occupations),
            encoding,
        })
    }

    /// Returns the encoding that lays the levels out on qubits.
    pub fn encoding(&self) -> LevelEncoding {
        self.encoding
    }

    /// Returns whether every mode shares one cutoff.
    pub fn is_uniform(&self) -> bool {
        matches!(self.cutoffs, Cutoffs::Uniform(_))
    }

    /// Returns the number of modes the spec covers, or `None` if it leaves that open.
    pub fn num_modes(&self) -> Option<u32> {
        match &self.cutoffs {
            Cutoffs::Uniform(_) => None,
            Cutoffs::PerMode(cutoffs) => Some(cutoffs.len() as u32),
        }
    }

    /// Returns the shared cutoff of a uniform spec, or `None` for a per-mode one.
    pub fn uniform_max_occupation(&self) -> Option<u32> {
        match &self.cutoffs {
            Cutoffs::Uniform(max_occupation) => Some(*max_occupation),
            Cutoffs::PerMode(_) => None,
        }
    }

    /// Returns the cutoff of every mode of a per-mode spec, or `None` for a uniform one.
    pub fn max_occupations(&self) -> Option<&[u32]> {
        match &self.cutoffs {
            Cutoffs::Uniform(_) => None,
            Cutoffs::PerMode(cutoffs) => Some(cutoffs),
        }
    }

    /// Applies the spec to exactly `num_modes` modes.
    ///
    /// A uniform spec is expanded to that many modes. A per-mode spec must cover exactly that many.
    pub fn mode_layout(&self, num_modes: u32) -> Result<ModeLayout, LevelSpecError> {
        let max_occupations = match &self.cutoffs {
            Cutoffs::Uniform(max_occupation) => vec![*max_occupation; num_modes as usize],
            Cutoffs::PerMode(cutoffs) => {
                if cutoffs.len() != num_modes as usize {
                    return Err(LevelSpecError::ModeCountMismatch {
                        expected: cutoffs.len(),
                        got: num_modes as usize,
                    });
                }
                cutoffs.clone()
            }
        };
        self.build_layout(max_occupations)
    }

    /// Applies the spec to an operator acting on `min_num_modes` modes, as a mapper needs it.
    ///
    /// A per-mode spec must cover at least that many modes and keeps all of its own. A uniform spec
    /// is expanded to exactly `min_num_modes` modes. With `num_qubits` given, the register may be
    /// wider than the modes need, in which case the surplus qubits are idle and sit at the end;
    /// without it the register is exactly as wide as the modes. Returns the layout and the
    /// register width.
    pub fn mapping_layout(
        &self,
        min_num_modes: u32,
        num_qubits: Option<u32>,
    ) -> Result<(ModeLayout, u32), LevelSpecError> {
        let layout = match &self.cutoffs {
            Cutoffs::Uniform(_) => self.mode_layout(min_num_modes)?,
            Cutoffs::PerMode(cutoffs) => {
                let available = cutoffs.len() as u32;
                if min_num_modes > available {
                    return Err(LevelSpecError::TooFewModes {
                        needed: min_num_modes,
                        available,
                    });
                }
                self.mode_layout(available)?
            }
        };
        let num_qubits = match num_qubits {
            None => layout.total_width,
            Some(num_qubits) if num_qubits >= layout.total_width => num_qubits,
            Some(num_qubits) => {
                return Err(LevelSpecError::NumQubitsTooSmall {
                    num_qubits,
                    needed: layout.total_width,
                });
            }
        };
        Ok((layout, num_qubits))
    }

    /// Returns by how much each mode's cutoff falls short of the given occupation lower bounds.
    ///
    /// `bounds[m]` is a lower bound on the cutoff mode `m` needs, for example from
    /// `BosonOperator::mode_occupations`. The result holds, per mode, `bound - max_occupation`
    /// clipped at zero, so an all-zero result means no bound is violated. It has one entry per mode
    /// of a per-mode spec, with a missing bound counting as zero, and one entry per bound for a
    /// uniform one.
    ///
    /// This measures the lower-bound heuristic only. It does not say whether the mapping is
    /// accurate, and a model whose top-level terms are exactly zero can report a shortfall although
    /// the mapping is exact.
    pub fn shortfalls(&self, bounds: &[u32]) -> Result<Vec<u32>, LevelSpecError> {
        match &self.cutoffs {
            Cutoffs::Uniform(max_occupation) => Ok(bounds
                .iter()
                .map(|bound| bound.saturating_sub(*max_occupation))
                .collect()),
            Cutoffs::PerMode(cutoffs) => {
                if bounds.len() > cutoffs.len() {
                    return Err(LevelSpecError::TooFewModes {
                        needed: bounds.len() as u32,
                        available: cutoffs.len() as u32,
                    });
                }
                Ok(cutoffs
                    .iter()
                    .enumerate()
                    .map(|(mode, max_occupation)| {
                        bounds
                            .get(mode)
                            .copied()
                            .unwrap_or(0)
                            .saturating_sub(*max_occupation)
                    })
                    .collect())
            }
        }
    }

    fn build_layout(&self, max_occupations: Vec<u32>) -> Result<ModeLayout, LevelSpecError> {
        let mut widths = Vec::with_capacity(max_occupations.len());
        let mut offsets = Vec::with_capacity(max_occupations.len());
        // The running total is the offset of the next mode and, once every mode is placed, the
        // width of the whole register.
        let mut total_width = 0u32;
        for &max_occupation in &max_occupations {
            let width = width(self.encoding, max_occupation);
            widths.push(width);
            offsets.push(total_width);
            total_width = total_width
                .checked_add(width)
                .ok_or(LevelSpecError::TooManyQubits)?;
        }
        Ok(ModeLayout {
            max_occupations,
            widths,
            offsets,
            total_width,
        })
    }
}

/// Checks that a cutoff can be laid out by `encoding`.
///
/// Every encoding is limited by [`MAX_CUTOFF`]. The unary encodings, projected and unprojected
/// alike, are limited further, because they identify a level by a single bit of a mask that holds
/// [`MAX_UNARY_LEVELS`] of them. This is the one place that limit is enforced for a spec, so nothing
/// built from a spec can reach the shift in [`LevelEncoding::level_to_bits`] with a level beyond it.
fn check_cutoff(max_occupation: u32, encoding: LevelEncoding) -> Result<(), LevelSpecError> {
    if max_occupation > MAX_CUTOFF {
        return Err(LevelSpecError::CutoffTooLarge(max_occupation));
    }
    if encoding.is_unary() && max_occupation >= MAX_UNARY_LEVELS {
        return Err(LevelSpecError::UnaryCutoffTooLarge(max_occupation));
    }
    Ok(())
}

fn width(encoding: LevelEncoding, max_occupation: u32) -> u32 {
    encoding.num_qubits(max_occupation + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ENCODINGS: [LevelEncoding; 4] = [
        LevelEncoding::UnaryProjected,
        LevelEncoding::UnaryUnprojected,
        LevelEncoding::StandardBinary,
        LevelEncoding::Gray,
    ];

    const UNARY: [LevelEncoding; 2] = [
        LevelEncoding::UnaryProjected,
        LevelEncoding::UnaryUnprojected,
    ];

    #[test]
    fn test_widths_and_offsets() {
        let spec = LevelSpec::per_mode(vec![4, 2, 2], LevelEncoding::Gray).unwrap();
        let layout = spec.mode_layout(3).unwrap();
        // d = 5, 3, 3 levels under a compact encoding
        assert_eq!(layout.widths, vec![3, 2, 2]);
        assert_eq!(layout.offsets, vec![0, 3, 5]);
        assert_eq!(layout.total_width, 7);

        let unary = LevelSpec::per_mode(vec![4, 2, 2], LevelEncoding::UnaryProjected).unwrap();
        let layout = unary.mode_layout(3).unwrap();
        assert_eq!(layout.widths, vec![5, 3, 3]);
        assert_eq!(layout.offsets, vec![0, 5, 8]);
        assert_eq!(layout.total_width, 11);
    }

    /// A uniform spec and a per-mode spec repeating its cutoff describe the same register.
    #[test]
    fn test_uniform_matches_repeated_per_mode() {
        for encoding in ENCODINGS {
            for max_occupation in 0..6 {
                let uniform = LevelSpec::uniform(max_occupation, encoding);
                let repeated = LevelSpec::per_mode(vec![max_occupation; 4], encoding);
                // the two are rejected together
                assert_eq!(uniform.is_err(), repeated.is_err());
                let (Ok(uniform), Ok(repeated)) = (uniform, repeated) else {
                    continue;
                };
                assert_eq!(uniform.mode_layout(4), repeated.mode_layout(4));
                assert_eq!(
                    uniform.mapping_layout(4, Some(40)),
                    repeated.mapping_layout(4, Some(40))
                );
            }
        }
    }

    /// A single-level mode takes no qubit under a compact encoding and does not move the next mode.
    #[test]
    fn test_width_zero_mode_does_not_advance_the_offset() {
        let spec = LevelSpec::per_mode(vec![3, 0, 3], LevelEncoding::StandardBinary).unwrap();
        let layout = spec.mode_layout(3).unwrap();
        assert_eq!(layout.widths, vec![2, 0, 2]);
        assert_eq!(layout.offsets, vec![0, 2, 2]);
        assert_eq!(layout.total_width, 4);

        // the same mode costs a qubit under the unary encoding
        let unary = LevelSpec::per_mode(vec![3, 0, 3], LevelEncoding::UnaryProjected).unwrap();
        assert_eq!(unary.mode_layout(3).unwrap().widths, vec![4, 1, 4]);
    }

    #[test]
    fn test_empty_register_is_rejected() {
        for encoding in [LevelEncoding::StandardBinary, LevelEncoding::Gray] {
            assert_eq!(
                LevelSpec::uniform(0, encoding).unwrap_err(),
                LevelSpecError::NoQubits
            );
            assert_eq!(
                LevelSpec::per_mode(vec![0, 0], encoding).unwrap_err(),
                LevelSpecError::NoQubits
            );
            assert_eq!(
                LevelSpec::per_mode(vec![], encoding).unwrap_err(),
                LevelSpecError::NoQubits
            );
        }
        // the unary encoding gives even a single-level mode one qubit
        assert!(LevelSpec::uniform(0, LevelEncoding::UnaryProjected).is_ok());
        assert!(LevelSpec::uniform(0, LevelEncoding::UnaryUnprojected).is_ok());
        assert_eq!(
            LevelSpec::per_mode(vec![], LevelEncoding::UnaryProjected).unwrap_err(),
            LevelSpecError::NoQubits
        );
    }

    #[test]
    fn test_cutoff_too_large_is_rejected() {
        assert_eq!(
            LevelSpec::uniform(u32::MAX, LevelEncoding::Gray).unwrap_err(),
            LevelSpecError::CutoffTooLarge(u32::MAX)
        );
        assert_eq!(
            LevelSpec::per_mode(vec![1, u32::MAX], LevelEncoding::Gray).unwrap_err(),
            LevelSpecError::CutoffTooLarge(u32::MAX)
        );
    }

    /// The unary encodings, projected and unprojected alike, hold exactly `MAX_UNARY_LEVELS` levels,
    /// no more, and the limit is the same for a shared cutoff and for one cutoff per mode.
    #[test]
    fn test_unary_cutoff_limit() {
        let last = MAX_UNARY_LEVELS - 1;

        for encoding in UNARY {
            // the largest cutoff fits and gets one qubit per level
            let uniform = LevelSpec::uniform(last, encoding).unwrap();
            assert_eq!(
                uniform.mode_layout(2).unwrap().widths,
                vec![MAX_UNARY_LEVELS; 2]
            );
            let per_mode = LevelSpec::per_mode(vec![last, 3], encoding).unwrap();
            assert_eq!(
                per_mode.mode_layout(2).unwrap().offsets,
                vec![0, MAX_UNARY_LEVELS]
            );

            // one level more does not
            assert_eq!(
                LevelSpec::uniform(MAX_UNARY_LEVELS, encoding).unwrap_err(),
                LevelSpecError::UnaryCutoffTooLarge(MAX_UNARY_LEVELS)
            );
            assert_eq!(
                LevelSpec::per_mode(vec![1, MAX_UNARY_LEVELS], encoding).unwrap_err(),
                LevelSpecError::UnaryCutoffTooLarge(MAX_UNARY_LEVELS)
            );
        }
    }

    /// The compact encodings are not limited by the width of the mask, only by `MAX_CUTOFF`.
    #[test]
    fn test_compact_encodings_are_not_limited_to_the_unary_levels() {
        // Both cutoffs exceed `MAX_UNARY_LEVELS`: 512 levels need 9 qubits, and the 513 levels of
        // `max_occupation = 512` need one more.
        for encoding in [LevelEncoding::StandardBinary, LevelEncoding::Gray] {
            let spec = LevelSpec::uniform(511, encoding).unwrap();
            assert_eq!(spec.mode_layout(1).unwrap().widths, vec![9]);
            let spec = LevelSpec::uniform(512, encoding).unwrap();
            assert_eq!(spec.mode_layout(1).unwrap().widths, vec![10]);
        }
    }

    #[test]
    fn test_mode_layout_checks_the_mode_count() {
        let spec = LevelSpec::per_mode(vec![1, 2], LevelEncoding::UnaryProjected).unwrap();
        assert_eq!(
            spec.mode_layout(3).unwrap_err(),
            LevelSpecError::ModeCountMismatch {
                expected: 2,
                got: 3
            }
        );
        // a uniform spec fits any mode count
        let uniform = LevelSpec::uniform(1, LevelEncoding::UnaryProjected).unwrap();
        assert_eq!(uniform.mode_layout(5).unwrap().total_width, 10);
        assert_eq!(uniform.num_modes(), None);
        assert_eq!(spec.num_modes(), Some(2));
    }

    #[test]
    fn test_mapping_layout_pads_and_checks() {
        let uniform = LevelSpec::uniform(2, LevelEncoding::UnaryProjected).unwrap();
        // no register width: exactly as wide as the modes
        let (layout, num_qubits) = uniform.mapping_layout(2, None).unwrap();
        assert_eq!((layout.total_width, num_qubits), (6, 6));
        // a wider register keeps the layout and idles the rest
        let (layout, num_qubits) = uniform.mapping_layout(2, Some(9)).unwrap();
        assert_eq!((layout.total_width, num_qubits), (6, 9));
        assert_eq!(
            uniform.mapping_layout(2, Some(5)).unwrap_err(),
            LevelSpecError::NumQubitsTooSmall {
                num_qubits: 5,
                needed: 6
            }
        );

        // a per-mode spec keeps all of its modes even if the operator touches fewer
        let per_mode = LevelSpec::per_mode(vec![2, 2, 1], LevelEncoding::UnaryProjected).unwrap();
        let (layout, num_qubits) = per_mode.mapping_layout(1, None).unwrap();
        assert_eq!((layout.widths.len(), num_qubits), (3, 8));
        assert_eq!(
            per_mode.mapping_layout(4, None).unwrap_err(),
            LevelSpecError::TooFewModes {
                needed: 4,
                available: 3
            }
        );
    }

    #[test]
    fn test_shortfalls() {
        let uniform = LevelSpec::uniform(1, LevelEncoding::UnaryProjected).unwrap();
        assert_eq!(uniform.shortfalls(&[2, 0, 1]).unwrap(), vec![1, 0, 0]);
        assert_eq!(uniform.shortfalls(&[]).unwrap(), Vec::<u32>::new());

        let per_mode =
            LevelSpec::per_mode(vec![3, 1, 0, 5], LevelEncoding::UnaryProjected).unwrap();
        // untouched trailing modes count as zero and the result covers every mode of the spec
        assert_eq!(per_mode.shortfalls(&[2, 4, 1]).unwrap(), vec![0, 3, 1, 0]);
        assert_eq!(
            per_mode.shortfalls(&[0; 5]).unwrap_err(),
            LevelSpecError::TooFewModes {
                needed: 5,
                available: 4
            }
        );
    }
}
