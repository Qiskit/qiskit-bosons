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

//! The bosonic second-quantization operator.
//!
//! # Data layout and algebra
//!
//! Terms are stored as a CSR-style flattening of `(action, mode)` pairs. Three facts about the
//! bosonic algebra shape every method below:
//!
//! 1. **Commutation relations**: `[b_i, b_j^dag] = delta_ij` with `[b_i, b_j] =
//!    [b_i^dag, b_j^dag] = 0`. Reordering two ladder operators never introduces a sign.
//! 2. **No exclusion principle**: `b_i^dag b_i^dag` is a perfectly good operator (it creates two
//!    bosons in mode `i`). Repeated `(action, mode)` pairs must be preserved rather than
//!    annihilating the term.
//! 3. **Unbounded occupation**: a bosonic mode has an infinite-dimensional Fock space. This
//!    operator is exact and truncation-free; an occupation cutoff enters only at the
//!    boson-to-qubit mapping stage (see [`crate::mappers`]).
//!
//! Points 1 and 2 together determine the normal-ordering routine: contracting
//! `b_i b_i^dag -> b_i^dag b_i + 1` adds a term with an unflipped coefficient, and no term is ever
//! dropped as identically zero.

use crate::operators::terms::grouping::group_term_range;
use crate::operators::{
    CoherenceError, GroupedTerm, OperatorMacro, OperatorTrait, ScaledTerm, TermSortKey,
};
use num_complex::{Complex64, ComplexFloat};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::iter::zip;
use std::ops::{
    Add, AddAssign, BitAnd, BitAndAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// A single bosonic ladder action: `(is_creation, mode_index)`.
pub type BosonAction<'a> = (&'a bool, &'a u32);

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BosonOperatorTermView<'a> {
    pub coeff: Complex64,
    pub actions: &'a [bool],
    pub modes: &'a [u32],
}

impl BosonOperatorTermView<'_> {
    pub fn iter(&'_ self) -> impl ExactSizeIterator<Item = BosonAction<'_>> + '_ {
        zip(self.actions, self.modes)
    }

    pub fn to_vec(&'_ self) -> Vec<BosonAction<'_>> {
        zip(self.actions, self.modes).collect()
    }

    pub fn into_vec(&'_ self) -> Vec<(bool, u32)> {
        zip(self.actions.iter().copied(), self.modes.iter().copied()).collect()
    }
}

impl TermSortKey for BosonOperatorTermView<'_> {
    fn sort_key(&self) -> impl Ord {
        // Compare the operator string position-by-position, each factor being an (action, mode)
        // pair. This matches `into_vec` (and hence the sorted display order), so the canonical
        // order agrees with how terms are printed.
        self.into_vec()
    }
}

impl ScaledTerm for BosonOperatorTermView<'_> {
    fn scaled(mut self, factor: Complex64) -> Self {
        // Only the coefficient changes; the index data stays borrowed from the operator this view
        // came from, so nothing is copied.
        self.coeff *= factor;
        self
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BosonOperatorGroupTermView<'a> {
    pub coeff: Complex64,
    pub actions: &'a [bool],
    pub modes: &'a [u32],
    pub group: u32,
}

impl BosonOperatorGroupTermView<'_> {
    pub fn iter(&'_ self) -> impl ExactSizeIterator<Item = BosonAction<'_>> + '_ {
        zip(self.actions, self.modes)
    }

    pub fn to_vec(&'_ self) -> Vec<BosonAction<'_>> {
        zip(self.actions, self.modes).collect()
    }

    pub fn into_vec(&'_ self) -> Vec<(bool, u32)> {
        zip(self.actions.iter().copied(), self.modes.iter().copied()).collect()
    }
}

impl TermSortKey for BosonOperatorGroupTermView<'_> {
    fn sort_key(&self) -> impl Ord {
        // Match the ungrouped `TermView` key exactly (ignoring `group`), so ordering a grouped
        // operator agrees with ordering the same terms ungrouped.
        self.into_vec()
    }
}

impl ScaledTerm for BosonOperatorGroupTermView<'_> {
    fn scaled(mut self, factor: Complex64) -> Self {
        self.coeff *= factor;
        self
    }
}

impl GroupedTerm for BosonOperatorGroupTermView<'_> {
    fn group(&self) -> u32 {
        self.group
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BosonOperator {
    pub coeffs: Vec<Complex64>,
    pub actions: Vec<bool>,
    pub modes: Vec<u32>,
    pub boundaries: Vec<usize>,
    pub groups: Option<Vec<u32>>,
}

crate::impl_operator_macro!(BosonOperator);

impl BosonOperator {
    #[inline]
    pub fn actions(&self) -> &[bool] {
        &self.actions
    }

    #[inline]
    pub fn modes(&self) -> &[u32] {
        &self.modes
    }

    #[inline]
    pub fn boundaries(&self) -> &[usize] {
        &self.boundaries
    }

    /// Appends a single term to the operator. Low-level building block for construction.
    ///
    /// # Warning
    ///
    /// This does **not** maintain `groups`: it pushes a coefficient and a boundary but no group
    /// index. Only call it while `self.groups` is `None` (e.g. on a fresh [`Self::zero`]); calling
    /// it on an operator that tracks groups leaves `coeffs.len()` out of sync with `groups.len()`,
    /// after which [`iter_with_groups`](Self::iter_with_groups) silently drops the trailing terms.
    pub fn _append_term(&mut self, coeff: Complex64, actions: &[bool], modes: &[u32]) {
        self.coeffs.push(coeff);
        self.actions.extend_from_slice(actions);
        self.modes.extend_from_slice(modes);
        self.boundaries.push(self.modes.len());
    }

    /// Splits this operator into new operators based on its [`groups`](Self::groups).
    ///
    /// If `group_indices` is `None`, every group is materialized once, in index order (`0` to
    /// [`num_groups`](Self::num_groups) `- 1`). If `group_indices` is `Some`, only the requested
    /// indices are materialized, in the given order; a duplicate index is materialized once
    /// internally and cloned once per occurrence in the output. This is significantly cheaper than
    /// requesting all indices when only a small number of groups out of a much larger total are
    /// needed, e.g. when subsampling groups for a randomized product formula, since terms
    /// belonging to a group that is not requested are skipped rather than appended anywhere.
    ///
    /// Returns `None` if `self.groups` is `None`, regardless of `group_indices`.
    pub fn split_out_groups(&self, group_indices: Option<&[u32]>) -> Option<Vec<Self>> {
        match group_indices {
            None => {
                let mut groups = vec![Self::zero(); self.num_groups()? as usize];
                for term in self.iter_with_groups() {
                    groups[term.group as usize]._append_term(term.coeff, term.actions, term.modes);
                }
                Some(groups)
            }
            Some(group_indices) => {
                let groups = self.groups.as_ref()?;
                if groups.is_sorted() {
                    // Each group is one contiguous run of terms, so its range can be binary-searched
                    // and its terms gathered directly. That makes a lookup scale with the number of
                    // groups requested rather than the number of terms held, which is what makes
                    // repeatedly sampling a few groups out of a large operator affordable (see
                    // `group_order`, which establishes this layout). The check itself is a
                    // short-circuiting scan over a contiguous `&[u32]`, orders of magnitude cheaper
                    // than the term walk below.
                    return Some(
                        group_indices
                            .iter()
                            .map(|&idx| {
                                let mut acc = Self::zero();
                                // A requested index that no term carries yields a zero operator,
                                // matching the fallback below.
                                if let Some(range) = group_term_range(groups, idx) {
                                    for i in range {
                                        let (start, end) =
                                            (self.boundaries[i], self.boundaries[i + 1]);
                                        acc._append_term(
                                            self.coeffs[i],
                                            &self.actions[start..end],
                                            &self.modes[start..end],
                                        );
                                    }
                                }
                                acc
                            })
                            .collect(),
                    );
                }
                let mut wanted: HashMap<u32, Self> = group_indices
                    .iter()
                    .map(|&idx| (idx, Self::zero()))
                    .collect();
                for term in self.iter_with_groups() {
                    if let Some(acc) = wanted.get_mut(&term.group) {
                        acc._append_term(term.coeff, term.actions, term.modes);
                    }
                }
                Some(
                    group_indices
                        .iter()
                        .map(|idx| {
                            wanted
                                .get(idx)
                                .expect("every requested index was inserted above")
                                .clone()
                        })
                        .collect(),
                )
            }
        }
    }

    /// Returns the normal-ordered form of this operator.
    ///
    /// Normal order places all creation operators to the left of all annihilation operators, with
    /// each of the two runs sorted by mode index. No sign tracking is involved: bosonic ladder
    /// operators commute, up to the contraction term generated by `[b_i, b_i^dag] = 1`.
    ///
    /// The `sandwich` argument selects the ordering: `None` requests plain normal order, while
    /// `Some(reverse)` requests the "sandwiched" order used by expectation-value expressions, in
    /// which `reverse` selects the sort direction of each of the two runs.
    ///
    /// See [`_normal_ordered_term`] for the per-term algorithm.
    pub fn normal_ordered(&self, sandwich: Option<bool>) -> Self {
        let mut result = Self::zero();
        self.iter()
            .for_each(|term| result.__iadd__(&_normal_ordered_term(term, sandwich)));
        result
    }

    pub fn max_rank(&self) -> u32 {
        self.boundaries
            .windows(2)
            .map(|p| p[1] - p[0])
            .max()
            .unwrap_or(0) as u32
    }

    /// Returns whether every term conserves the total boson number.
    ///
    /// A term conserves boson number iff it contains as many creation as annihilation operators.
    /// A mode may appear more than once in a term (there is no exclusion principle), so this counts
    /// actions rather than distinct modes.
    pub fn conserves_particle_number(&self) -> bool {
        // Particle-number conservation is sector conservation with all modes in one block.
        self.conserves_sector(&[])
    }

    /// Returns whether every term conserves boson number within each mode block.
    ///
    /// `block_sizes` partitions the mode range into consecutive, non-overlapping blocks: block `b`
    /// spans modes `[start_b, start_b + block_sizes[b])` where `start_b = sum(block_sizes[..b])`. A
    /// term conserves the sector iff, in *every* block, its number of creation operators equals its
    /// number of annihilation operators. A term touching a mode beyond the last block does not
    /// conserve the sector (there is no block whose count it could balance).
    ///
    /// An empty `block_sizes` treats all modes as a single unbounded block, making this equivalent
    /// to [`Self::conserves_particle_number`].
    ///
    /// Blocks matter for bosonic systems too: a vibrational Hamiltonian in the Christiansen form
    /// groups the modals belonging to one vibrational mode into a block, and physical states carry
    /// exactly one excitation per block.
    pub fn conserves_sector(&self, block_sizes: &[u32]) -> bool {
        // The exclusive upper bound of each block; the last entry is the total mode count. Empty
        // `block_sizes` yields no bounds, which the block lookup below treats as one unbounded block.
        let bounds: Vec<u32> = block_sizes
            .iter()
            .scan(0, |acc, size| {
                *acc += size;
                Some(*acc)
            })
            .collect();
        // Net (creations - annihilations) per block, reused across terms.
        let mut net = vec![0i64; block_sizes.len().max(1)];
        for term in self.iter() {
            net.iter_mut().for_each(|n| *n = 0);
            let mut in_range = true;
            for (&action, &mode) in term.iter() {
                // The block containing `mode`: the first bound strictly greater than it. With no
                // bounds every mode falls into the single block 0.
                let block = match bounds.iter().position(|&b| mode < b) {
                    Some(b) => b,
                    None if bounds.is_empty() => 0,
                    None => {
                        in_range = false;
                        break;
                    }
                };
                net[block] += if action { 1 } else { -1 };
            }
            if !in_range || net.iter().any(|&n| n != 0) {
                return false;
            }
        }
        true
    }

    /// Returns, for each mode, the highest occupation it is driven to by any one term.
    ///
    /// The occupation of a mode in a term is the number of creation operators acting on that mode
    /// within the term. It is a *lower bound* on the truncation level a boson-to-qubit mapping needs
    /// for that mode in order to represent this operator without clipping: `b_i^dag b_i^dag` applied
    /// to the vacuum populates occupation 2, which a cutoff of `max_occupation < 2` cannot represent.
    ///
    /// It is only a lower bound, not a sufficient cutoff: the occupation reachable by *repeated*
    /// application (e.g. over a time evolution) is unbounded, so choosing the cutoff remains a
    /// physical modeling decision.
    ///
    /// If `modes` is `None`, the result has one entry per mode, `0..min_num_modes()`, with a zero for
    /// every mode the operator does not touch. Otherwise it holds one entry per requested mode, in the
    /// order requested, and a mode the operator does not touch (including one beyond
    /// `min_num_modes()`) reports zero.
    pub fn mode_occupations(&self, modes: Option<&[u32]>) -> Vec<u32> {
        let mut highest = vec![0u32; self.min_num_modes() as usize];
        let mut counts: HashMap<u32, u32> = HashMap::new();
        for term in self.iter() {
            counts.clear();
            for (&action, &mode) in term.iter() {
                if action {
                    let c = counts.entry(mode).or_insert(0);
                    *c += 1;
                    let slot = &mut highest[mode as usize];
                    *slot = (*slot).max(*c);
                }
            }
        }
        match modes {
            None => highest,
            Some(modes) => modes
                .iter()
                .map(|&mode| highest.get(mode as usize).copied().unwrap_or(0))
                .collect(),
        }
    }

    /// Returns the highest occupation any single mode is driven to by any one term.
    ///
    /// This is the maximum of [`Self::mode_occupations`] over all modes, and carries the same caveat:
    /// it is a lower bound on the truncation level, not a sufficient cutoff.
    pub fn max_occupation(&self) -> u32 {
        self.mode_occupations(None).into_iter().max().unwrap_or(0)
    }
}

/// Normal-orders a single term, returning the (possibly multi-term) result.
///
/// The algorithm is an insertion sort over the term's ladder operators, driven by a worklist so
/// that the extra terms produced by a contraction are themselves normal-ordered. Three properties
/// characterize it:
///
/// - **No sign tracking at all.** `[b_i, b_j] = 0` and `[b_i^dag, b_j^dag] = 0`, so swapping two
///   like operators is free, and swapping `b_i b_j^dag -> b_j^dag b_i` for `i != j` is free too.
/// - **The contraction adds.** For `i == j`, `b_i b_i^dag = b_i^dag b_i + 1`, so pushing an
///   annihilation past a creation on the same mode spawns a new term with that pair deleted and the
///   coefficient **unchanged in sign**.
/// - **Nothing vanishes.** Equal adjacent like-operators are kept, since `b_i^dag b_i^dag != 0`, so
///   there is no zero-term short circuit.
fn _normal_ordered_term(term_view: BosonOperatorTermView, sandwich: Option<bool>) -> BosonOperator {
    let mut coeffs = vec![];
    let mut actions = vec![];
    let mut modes = vec![];
    let mut boundaries = vec![0];

    let mut stack = vec![(term_view.into_vec(), term_view.coeff)];
    while let Some((mut term, coeff)) = stack.pop() {
        for i in 1..term.len() {
            // shift the operator at index i to the left until it's in the correct location
            for j in (1..=i).rev() {
                let (action_right, index_right) = term[j];
                let (action_left, index_left) = term[j - 1];
                if action_right == action_left {
                    // Both create or both destroy: these always commute, including when the mode
                    // indices are equal. Sort by mode index; `sandwich` selects the direction for
                    // each run.
                    match (index_right.cmp(&index_left), sandwich, action_left) {
                        (Ordering::Equal, _, _) => {
                            // identical operators: already ordered, and emphatically not zero
                        }
                        (Ordering::Greater, None, _)
                        | (Ordering::Less, Some(true), true)
                        | (Ordering::Less, Some(false), false)
                        | (Ordering::Greater, Some(true), false)
                        | (Ordering::Greater, Some(false), true) => {
                            term.swap(j - 1, j);
                        }
                        (Ordering::Less, None, _)
                        | (Ordering::Greater, Some(true), true)
                        | (Ordering::Greater, Some(false), false)
                        | (Ordering::Less, Some(true), false)
                        | (Ordering::Less, Some(false), true) => {}
                    }
                } else if action_right && !action_left {
                    // Annihilation on the left, creation on the right: swap so creations migrate
                    // leftwards. On the same mode this is `b_i b_i^dag = b_i^dag b_i + 1`, whose
                    // `+1` spawns a contracted term.
                    if index_right == index_left {
                        let mut new_term: Vec<(bool, u32)> = Vec::new();
                        new_term.extend_from_slice(&term[..j - 1]);
                        new_term.extend_from_slice(&term[j + 1..]);
                        // NOTE: `+coeff`, not `-coeff` -- this is the whole bosonic difference.
                        stack.push((new_term, coeff));
                    }
                    term.swap(j - 1, j);
                }
            }
        }
        coeffs.push(coeff);
        term.iter().for_each(|&(a, i)| {
            actions.push(a);
            modes.push(i);
        });
        boundaries.push(modes.len())
    }
    BosonOperator {
        coeffs,
        actions,
        modes,
        boundaries,
        groups: None,
    }
}

fn _compose(
    a: &BosonOperator,
    b: &BosonOperator,
) -> (Vec<Complex64>, Vec<bool>, Vec<u32>, Vec<usize>) {
    // The output size is known exactly: one term per (left, right) pair, each holding the factors
    // of both inputs. Reserving up front turns the inner loop's repeated `extend_from_slice` into
    // plain memcpy without the reallocation-and-copy rounds an unreserved vector would incur.
    let num_terms = a.coeffs.len() * b.coeffs.len();
    let num_factors = a.modes.len() * b.coeffs.len() + b.modes.len() * a.coeffs.len();
    let mut coeffs = Vec::with_capacity(num_terms);
    let mut actions = Vec::with_capacity(num_factors);
    let mut modes = Vec::with_capacity(num_factors);
    let mut boundaries = Vec::with_capacity(num_terms + 1);
    boundaries.push(0);

    for left in a.iter() {
        for right in b.iter() {
            coeffs.push(left.coeff * right.coeff);
            actions.extend_from_slice(right.actions);
            actions.extend_from_slice(left.actions);
            modes.extend_from_slice(right.modes);
            modes.extend_from_slice(left.modes);
            boundaries.push(modes.len());
        }
    }
    (coeffs, actions, modes, boundaries)
}

impl OperatorTrait for BosonOperator {
    type TermView<'a> = BosonOperatorTermView<'a>;
    type GroupTermView<'a> = BosonOperatorGroupTermView<'a>;

    fn zero() -> Self {
        Self {
            coeffs: vec![],
            actions: vec![],
            modes: vec![],
            boundaries: vec![0],
            groups: None,
        }
    }

    fn one() -> Self {
        Self {
            coeffs: vec![Complex64::new(1.0, 0.0)],
            actions: vec![],
            modes: vec![],
            boundaries: vec![0, 0],
            groups: None,
        }
    }

    /// Returns whether two operators are algebraically equal, within `atol`.
    ///
    /// The difference is brought to normal order before being compared. Two bosonic operator
    /// strings can be algebraically equal while being structurally different (`b b^dag` versus
    /// `b^dag b + 1`), and the purely structural [`Self::simplify`] would not identify them. Normal
    /// order is the canonical form that does.
    fn equiv(&self, other: &Self, atol: f64) -> bool {
        let mut diff = self.__sub__(other).normal_ordered(None);
        diff = diff.simplify(atol);
        for c in diff.coeffs {
            if c.abs() > atol {
                return false;
            }
        }
        true
    }

    /// Returns whether the operator equals its own adjoint, within `atol`.
    ///
    /// Compared in normal order, for the same reason as [`Self::equiv`]: a Hermitian operator such
    /// as `b b^dag + b^dag b` is only recognizably equal to its adjoint once the commutation
    /// relations have been applied.
    fn is_hermitian(&self, atol: f64) -> bool {
        let mut diff = (self.__sub__(&self.adjoint())).normal_ordered(None);
        diff.ichop(atol);
        diff.equiv(&Self::zero(), atol)
    }

    /// Collects structurally identical terms and drops those below `atol`.
    ///
    /// This is a purely *structural* simplification: it combines terms whose `(actions, modes)`
    /// strings match exactly. It does **not** apply the commutation relations, so `b b^dag` and
    /// `b^dag b + 1` are left distinct. Use [`Self::normal_ordered`] first for an algebraic
    /// canonical form (which is what [`Self::equiv`] does).
    fn simplify(&self, atol: f64) -> Self {
        let mut terms = HashMap::new();
        for term in self.iter() {
            terms
                .entry((term.modes, term.actions))
                .and_modify(|c| *c += term.coeff)
                .or_insert(term.coeff);
        }
        let mut out = Self::zero();
        terms
            .iter()
            .filter(|((_, _), coeff)| coeff.abs() > atol)
            .for_each(|((modes, actions), coeff)| out._append_term(*coeff, actions, modes));
        out
    }

    fn adjoint(&self) -> Self {
        let mut coeffs = vec![];
        let mut actions = vec![];
        let mut modes = vec![];

        self.iter().for_each(|term| {
            coeffs.push(term.coeff.conj());
            actions.extend(term.actions.iter().rev().map(|a| !a));
            modes.extend(term.modes.iter().rev());
        });

        Self {
            coeffs,
            actions,
            modes,
            boundaries: self.boundaries.clone(),
            groups: self.groups.clone(),
        }
    }

    fn __iadd__(&mut self, other: &Self) {
        self.coeffs.extend_from_slice(&other.coeffs);
        self.actions.extend_from_slice(&other.actions);
        self.modes.extend_from_slice(&other.modes);
        let offset = self.boundaries[self.boundaries.len() - 1];
        self.boundaries
            .extend(other.boundaries[1..].iter().map(|b| b + offset));
        self.clear_groups();
    }

    fn __imul__(&mut self, other: Complex64) {
        self.coeffs.iter_mut().for_each(|c| *c *= other);
    }

    fn __iscaled_add__(&mut self, other: &Self, factor: Complex64) {
        self.coeffs.extend(other.coeffs.iter().map(|c| c * factor));
        self.actions.extend_from_slice(&other.actions);
        self.modes.extend_from_slice(&other.modes);
        let offset = self.boundaries[self.boundaries.len() - 1];
        self.boundaries
            .extend(other.boundaries[1..].iter().map(|b| b + offset));
        self.clear_groups();
    }

    fn composed(&self, other: &Self) -> Self {
        let (coeffs, actions, modes, boundaries) = _compose(self, other);
        Self {
            coeffs,
            actions,
            modes,
            boundaries,
            groups: None,
        }
    }

    fn matmul(&self, other: &Self) -> Self {
        let (coeffs, actions, modes, boundaries) = _compose(other, self);
        Self {
            coeffs,
            actions,
            modes,
            boundaries,
            groups: None,
        }
    }

    fn __iand__(&mut self, other: &Self) {
        *self = self.composed(other);
    }

    fn __imatmul__(&mut self, other: &Self) {
        *self = self.matmul(other);
    }

    fn ichop(&mut self, atol: f64) {
        let mut coeffs = vec![];
        let mut actions = vec![];
        let mut modes = vec![];
        let mut boundaries = vec![0];

        self.iter()
            .filter(|term| term.coeff.abs() > atol)
            .for_each(|term| {
                coeffs.push(term.coeff);
                actions.extend_from_slice(term.actions);
                modes.extend_from_slice(term.modes);
                boundaries.push(modes.len());
            });

        self.coeffs = coeffs;
        self.actions = actions;
        self.modes = modes;
        self.boundaries = boundaries;
        self.clear_groups();
    }

    fn iter(&self) -> impl ExactSizeIterator<Item = Self::TermView<'_>> {
        self.coeffs.iter().enumerate().map(|(i, coeff)| {
            let start = self.boundaries[i];
            let end = self.boundaries[i + 1];
            BosonOperatorTermView {
                coeff: *coeff,
                actions: &self.actions[start..end],
                modes: &self.modes[start..end],
            }
        })
    }

    #[inline]
    fn coeffs(&self) -> &[Complex64] {
        &self.coeffs
    }

    fn groups(&self) -> Option<&[u32]> {
        self.groups.as_deref()
    }

    #[inline]
    fn _write_groups(&mut self, groups: Option<Vec<u32>>) {
        self.groups = groups;
    }

    fn iter_with_groups(&self) -> impl ExactSizeIterator<Item = Self::GroupTermView<'_>> {
        if let Some(groups) = &self.groups {
            return self
                .coeffs
                .iter()
                .zip(groups)
                .enumerate()
                .map(|(i, (coeff, gidx))| {
                    let start = self.boundaries[i];
                    let end = self.boundaries[i + 1];
                    BosonOperatorGroupTermView {
                        coeff: *coeff,
                        actions: &self.actions[start..end],
                        modes: &self.modes[start..end],
                        group: *gidx,
                    }
                });
        }
        panic!("This method can only be called when groups are present!");
    }

    fn from_terms<'a, I>(terms: I) -> Self
    where
        Self: 'a,
        I: IntoIterator<Item = Self::TermView<'a>>,
    {
        let mut out = Self::zero();
        for term in terms {
            out._append_term(term.coeff, term.actions, term.modes);
        }
        out
    }

    fn from_terms_with_groups<'a, I>(terms: I) -> Self
    where
        Self: 'a,
        I: IntoIterator<Item = Self::GroupTermView<'a>>,
    {
        let mut out = Self::zero();
        let mut groups = Vec::new();
        for term in terms {
            out._append_term(term.coeff, term.actions, term.modes);
            groups.push(term.group);
        }
        // Built in lockstep with the terms above, so the length check of `set_groups` is redundant.
        out._write_groups(Some(groups));
        out
    }

    fn get_support(&self) -> HashSet<u32> {
        self.modes.iter().copied().collect()
    }

    fn max_mode(&self) -> Option<u32> {
        self.modes.iter().copied().max()
    }

    fn relabel_modes(&self, permutation: Vec<u32>) -> Result<Self, CoherenceError> {
        if permutation.iter().collect::<HashSet<_>>().len() != permutation.len() {
            return Err(CoherenceError::DuplicateIndices);
        }
        // The relabeled modes are computed before anything is copied, so that the error path does
        // not first clone the entire operator only to discard it.
        let modes: Vec<u32> = self
            .modes
            .iter()
            .map(|&idx| {
                permutation
                    .get(idx as usize)
                    .copied()
                    .ok_or(CoherenceError::IndexMapTooSmall)
            })
            .collect::<Result<_, _>>()?;
        Ok(Self {
            coeffs: self.coeffs.clone(),
            actions: self.actions.clone(),
            modes,
            boundaries: self.boundaries.clone(),
            // Relabeling permutes mode indices without reordering, splitting or merging terms, so
            // any grouping of those terms carries over unchanged. This is the one operation that
            // preserves `groups` rather than dropping it.
            groups: self.groups.clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds an operator from a list of `(coeff, [(action, mode), ...])` terms.
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

    #[test]
    fn test_zero() {
        let zero = BosonOperator::zero();
        assert_eq!(
            zero,
            BosonOperator {
                coeffs: vec![],
                actions: vec![],
                modes: vec![],
                boundaries: vec![0],
                groups: None,
            }
        );
    }

    #[test]
    fn test_one() {
        let one = BosonOperator::one();
        assert_eq!(
            one,
            BosonOperator {
                coeffs: vec![Complex64::new(1.0, 0.0)],
                actions: vec![],
                modes: vec![],
                boundaries: vec![0, 0],
                groups: None,
            }
        );
    }

    #[test]
    fn test_scaled_add() {
        let one = BosonOperator::one();
        let two = op(&[(c(2.0), &[])]);

        // The appended coefficients are scaled, the existing ones are untouched.
        let result = two.__scaled_add__(&one, Complex64::new(0.0, 3.0));
        assert_eq!(
            result,
            BosonOperator {
                coeffs: vec![Complex64::new(2.0, 0.0), Complex64::new(0.0, 3.0)],
                actions: vec![],
                modes: vec![],
                boundaries: vec![0, 0, 0],
                groups: None,
            }
        );
    }

    /// `__isub__` is defined in terms of `__iscaled_add__`, so the two must agree exactly.
    #[test]
    fn test_scaled_add_subsumes_sub() {
        // A repeated mode on each side: no exclusion principle, so nothing may be dropped.
        let hop = op(&[(Complex64::new(1.5, -0.5), &[(CRE, 0), (CRE, 0), (ANN, 1)])]);
        let num = op(&[(Complex64::new(0.25, 2.0), &[(CRE, 1), (ANN, 1)])]);

        let mut subtracted = hop.clone();
        subtracted.__isub__(&num);

        let mut scaled = hop.clone();
        scaled.__iscaled_add__(&num, Complex64::new(-1.0, 0.0));

        assert_eq!(subtracted, scaled);

        // A unit factor is plain addition.
        let mut added = hop.clone();
        added.__iadd__(&num);
        let mut scaled_one = hop.clone();
        scaled_one.__iscaled_add__(&num, Complex64::new(1.0, 0.0));
        assert_eq!(added, scaled_one);
    }

    /// Appending terms invalidates the one-index-per-term invariant, so the grouping must go.
    #[test]
    fn test_scaled_add_clears_groups() {
        let mut grouped = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        grouped.set_groups(Some(vec![0])).unwrap();
        assert!(grouped.has_groups());

        grouped.__iscaled_add__(&BosonOperator::one(), Complex64::new(2.0, 0.0));
        assert!(!grouped.has_groups());
    }

    #[test]
    fn test_sub_matches_add_of_negation() {
        let op1 = op(&[
            (Complex64::new(1.5, -0.5), &[(CRE, 0), (ANN, 1)]),
            (
                Complex64::new(-2.0, 1.0),
                &[(CRE, 0), (CRE, 0), (ANN, 0), (ANN, 0)],
            ),
        ]);
        let op2 = op(&[(Complex64::new(0.25, 2.0), &[(CRE, 1), (ANN, 0)])]);

        // The formulation `__sub__` used before it was fused, spelled out here so that the fused
        // version is pinned to it. Complex coefficients matter: a sign error in the real part alone
        // would survive real-only operands.
        let mut expected = op1.clone();
        expected.__iadd__(&op2.__neg__());

        assert_eq!(op1.__sub__(&op2), expected);

        let mut in_place = op1.clone();
        in_place.__isub__(&op2);
        assert_eq!(in_place, expected);
    }

    #[test]
    fn test_add_sub_mul_div_neg() {
        let one = BosonOperator::one();
        let two = op(&[(c(2.0), &[])]);

        assert!((one.clone() + two.clone()).equiv(&op(&[(c(3.0), &[])]), 1e-12));
        assert!((two.clone() - one.clone()).equiv(&one, 1e-12));
        assert!((one.clone() * c(3.0)).equiv(&op(&[(c(3.0), &[])]), 1e-12));
        assert!((two.clone() / c(2.0)).equiv(&one, 1e-12));
        assert!((-one.clone()).equiv(&op(&[(c(-1.0), &[])]), 1e-12));

        let mut acc = BosonOperator::one();
        acc += two.clone();
        assert!(acc.equiv(&op(&[(c(3.0), &[])]), 1e-12));
        acc -= two;
        assert!(acc.equiv(&one, 1e-12));
    }

    /// `b_0^dag b_0^dag` is a legitimate, nonzero bosonic operator -- the single most important
    /// departure from the fermionic case, where `a^dag a^dag = 0`.
    #[test]
    fn test_repeated_creation_is_not_zero() {
        let double = op(&[(c(1.0), &[(CRE, 0), (CRE, 0)])]);

        // it must survive normal ordering unchanged...
        let ordered = double.normal_ordered(None);
        assert!(ordered.equiv(&double, 1e-12));
        // ...and must not be equivalent to zero
        assert!(!ordered.equiv(&BosonOperator::zero(), 1e-12));

        // the same holds for the square of a bare creation operator
        let squared = op(&[(c(1.0), &[(CRE, 0)])]).__pow__(2);
        assert!(squared.equiv(&double, 1e-12));
    }

    /// The canonical commutator `[b, b^dag] = 1`, i.e. `b_0 b_0^dag - b_0^dag b_0 = 1`.
    #[test]
    fn test_canonical_commutator() {
        let bbd = op(&[(c(1.0), &[(ANN, 0), (CRE, 0)])]);
        let dbb = op(&[(c(1.0), &[(CRE, 0), (ANN, 0)])]);

        let commutator = (bbd - dbb).normal_ordered(None);
        assert!(commutator.equiv(&BosonOperator::one(), 1e-12));
    }

    /// `b_0 b_0^dag` normal-orders to `b_0^dag b_0 + 1` -- note the **plus**, where the fermionic
    /// analogue produces `-a^dag a + 1`.
    #[test]
    fn test_normal_ordered_contraction_sign() {
        let bbd = op(&[(c(1.0), &[(ANN, 0), (CRE, 0)])]);
        let expected = op(&[(c(1.0), &[(CRE, 0), (ANN, 0)]), (c(1.0), &[])]);

        let ordered = bbd.normal_ordered(None);
        assert!(
            ordered.equiv(&expected, 1e-12),
            "b b^dag must normal-order to +b^dag b + 1"
        );

        // Explicitly assert the fermionic answer is NOT what we produce.
        let fermionic = op(&[(c(-1.0), &[(CRE, 0), (ANN, 0)]), (c(1.0), &[])]);
        assert!(!ordered.equiv(&fermionic, 1e-12));
    }

    /// Distinct modes commute freely and produce no contraction term.
    #[test]
    fn test_distinct_modes_commute() {
        let b0_b1d = op(&[(c(1.0), &[(ANN, 0), (CRE, 1)])]);
        let b1d_b0 = op(&[(c(1.0), &[(CRE, 1), (ANN, 0)])]);

        // [b_0, b_1^dag] = 0
        assert!(
            (b0_b1d.clone() - b1d_b0)
                .normal_ordered(None)
                .equiv(&BosonOperator::zero(), 1e-12)
        );

        // normal ordering moves the creation left without spawning an extra term
        let ordered = b0_b1d.normal_ordered(None).simplify(1e-12);
        assert_eq!(ordered.coeffs.len(), 1);
    }

    /// Two creations on distinct modes commute: `b_0^dag b_1^dag = b_1^dag b_0^dag`, no sign.
    #[test]
    fn test_creations_commute_without_sign() {
        let a = op(&[(c(1.0), &[(CRE, 0), (CRE, 1)])]);
        let b = op(&[(c(1.0), &[(CRE, 1), (CRE, 0)])]);
        assert!(
            (a - b)
                .normal_ordered(None)
                .equiv(&BosonOperator::zero(), 1e-12),
            "bosonic creation operators must commute with NO sign flip"
        );
    }

    /// `b b b^dag b^dag` on one mode expands via two contractions.
    ///
    /// `b b b^dag b^dag = b^dag b^dag b b + 4 b^dag b + 2`. The constant is pinned by acting on the
    /// vacuum: `b b b^dag b^dag |0> = 2 |0>`.
    #[test]
    fn test_double_contraction() {
        let lhs = op(&[(c(1.0), &[(ANN, 0), (ANN, 0), (CRE, 0), (CRE, 0)])]);
        let expected = op(&[
            (c(1.0), &[(CRE, 0), (CRE, 0), (ANN, 0), (ANN, 0)]),
            (c(4.0), &[(CRE, 0), (ANN, 0)]),
            (c(2.0), &[]),
        ]);
        assert!(lhs.normal_ordered(None).equiv(&expected, 1e-12));
    }

    #[test]
    fn test_normal_ordered_already_ordered_is_fixed_point() {
        let already = op(&[(c(1.0), &[(CRE, 1), (CRE, 0), (ANN, 1), (ANN, 0)])]);
        let ordered = already.normal_ordered(None);
        assert!(ordered.equiv(&already, 1e-12));
        // idempotence
        assert!(ordered.normal_ordered(None).equiv(&ordered, 1e-12));
    }

    #[test]
    fn test_adjoint() {
        // (b_0^dag b_1)^dag = b_1^dag b_0
        let o = op(&[(Complex64::new(0.0, 2.0), &[(CRE, 0), (ANN, 1)])]);
        let expected = op(&[(Complex64::new(0.0, -2.0), &[(CRE, 1), (ANN, 0)])]);
        assert!(o.adjoint().equiv(&expected, 1e-12));
    }

    #[test]
    fn test_is_hermitian() {
        // n_0 = b_0^dag b_0 is Hermitian
        let n0 = op(&[(c(1.0), &[(CRE, 0), (ANN, 0)])]);
        assert!(n0.is_hermitian(1e-12));

        // b_0^dag b_1 alone is not
        let hop = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        assert!(!hop.is_hermitian(1e-12));

        // but the symmetrized combination is
        let sym = op(&[
            (c(1.0), &[(CRE, 0), (ANN, 1)]),
            (c(1.0), &[(CRE, 1), (ANN, 0)]),
        ]);
        assert!(sym.is_hermitian(1e-12));

        // a bare creation operator is never Hermitian
        assert!(!op(&[(c(1.0), &[(CRE, 0)])]).is_hermitian(1e-12));
    }

    #[test]
    fn test_conserves_particle_number() {
        assert!(op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]).conserves_particle_number());
        assert!(!op(&[(c(1.0), &[(CRE, 0)])]).conserves_particle_number());
        // a two-boson creation does not conserve number
        assert!(!op(&[(c(1.0), &[(CRE, 0), (CRE, 0)])]).conserves_particle_number());
        // ...but a pair-hopping term does
        assert!(
            op(&[(c(1.0), &[(CRE, 0), (CRE, 0), (ANN, 1), (ANN, 1)])]).conserves_particle_number()
        );
    }

    #[test]
    fn test_conserves_sector() {
        let hop = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        assert!(hop.conserves_sector(&[2]));
        assert!(hop.conserves_sector(&[]));
        // crossing a block boundary breaks per-block conservation
        assert!(!hop.conserves_sector(&[1, 1]));

        // a mode beyond the last block cannot conserve
        let out_of_range = op(&[(c(1.0), &[(CRE, 0), (ANN, 4)])]);
        assert!(!out_of_range.conserves_sector(&[2, 2]));
    }

    #[test]
    fn test_max_rank() {
        assert_eq!(BosonOperator::zero().max_rank(), 0);
        assert_eq!(BosonOperator::one().max_rank(), 0);
        assert_eq!(op(&[(c(1.0), &[(CRE, 0)])]).max_rank(), 1);
        assert_eq!(
            op(&[(c(1.0), &[(CRE, 0), (CRE, 0), (ANN, 1), (ANN, 1)])]).max_rank(),
            4
        );
    }

    #[test]
    fn test_max_occupation() {
        assert_eq!(BosonOperator::zero().max_occupation(), 0);
        // annihilations do not raise occupation
        assert_eq!(op(&[(c(1.0), &[(ANN, 0), (ANN, 0)])]).max_occupation(), 0);
        assert_eq!(op(&[(c(1.0), &[(CRE, 0)])]).max_occupation(), 1);
        assert_eq!(
            op(&[(c(1.0), &[(CRE, 0), (CRE, 0), (CRE, 0)])]).max_occupation(),
            3
        );
        // counted per mode, not per term: two modes each raised once is still 1
        assert_eq!(op(&[(c(1.0), &[(CRE, 0), (CRE, 1)])]).max_occupation(), 1);
        // maximum is taken across terms
        assert_eq!(
            op(&[(c(1.0), &[(CRE, 0)]), (c(1.0), &[(CRE, 1), (CRE, 1)])]).max_occupation(),
            2
        );
    }

    #[test]
    fn test_mode_occupations() {
        assert_eq!(
            BosonOperator::zero().mode_occupations(None),
            Vec::<u32>::new()
        );
        assert_eq!(
            BosonOperator::one().mode_occupations(None),
            Vec::<u32>::new()
        );

        // modes 0 and 3 are driven, 1 and 2 are untouched and report zero; the list is as long as
        // `min_num_modes()`
        let o = op(&[
            (c(1.0), &[(CRE, 0), (CRE, 0), (ANN, 3)]),
            (c(1.0), &[(CRE, 3), (CRE, 3), (CRE, 3)]),
        ]);
        assert_eq!(o.min_num_modes(), 4);
        assert_eq!(o.mode_occupations(None), vec![2, 0, 0, 3]);

        // a query returns the requested entries in the requested order, zero beyond the operator
        assert_eq!(o.mode_occupations(Some(&[3, 1, 0])), vec![3, 0, 2]);
        assert_eq!(o.mode_occupations(Some(&[7])), vec![0]);
        assert_eq!(o.mode_occupations(Some(&[])), Vec::<u32>::new());

        // the maximum over modes is `max_occupation`
        assert_eq!(o.max_occupation(), 3);
    }

    #[test]
    fn test_max_mode_and_min_num_modes() {
        assert_eq!(BosonOperator::zero().max_mode(), None);
        assert_eq!(BosonOperator::one().max_mode(), None);
        let o = op(&[(c(1.0), &[(CRE, 2), (ANN, 5)])]);
        assert_eq!(o.max_mode(), Some(5));
        assert_eq!(o.min_num_modes(), 6);
        assert_eq!(o.max_mode(), o.get_support().into_iter().max());
    }

    #[test]
    fn test_compose_and_matmul() {
        let a = op(&[(c(2.0), &[(CRE, 0)])]);
        let b = op(&[(c(3.0), &[(ANN, 1)])]);

        // `a & b` == "first a, then b" == B A
        let and = (a.clone() & b.clone()).simplify(1e-12);
        assert_eq!(and.coeffs.len(), 1);
        assert!(and.equiv(&op(&[(c(6.0), &[(ANN, 1), (CRE, 0)])]), 1e-12));

        // `a @ b` == A B, the reverse
        assert!(a.__matmul__(&b).equiv(&b.__and__(&a), 1e-12));
    }

    #[test]
    fn test_ichop() {
        let mut o = op(&[(c(1e-4), &[(CRE, 0)]), (c(1e-8), &[(ANN, 0)])]);
        o.ichop(1e-6);
        assert_eq!(o.coeffs.len(), 1);
        assert!(o.equiv(&op(&[(c(1e-4), &[(CRE, 0)])]), 1e-12));
    }

    /// `ichop` must not conjugate the coefficients it keeps.
    #[test]
    fn test_ichop_preserves_complex_coeff() {
        let mut o = op(&[(Complex64::new(0.0, 1.0), &[(CRE, 0)])]);
        o.ichop(1e-12);
        assert_eq!(o.coeffs, vec![Complex64::new(0.0, 1.0)]);
    }

    #[test]
    fn test_get_support() {
        let o = op(&[
            (c(1.0), &[(CRE, 0), (ANN, 4)]),
            (c(1.0), &[(CRE, 1), (CRE, 3), (ANN, 4), (ANN, 7)]),
        ]);
        assert_eq!(o.get_support(), HashSet::from([0, 1, 3, 4, 7]));
    }

    #[test]
    fn test_relabel_modes() {
        let o = op(&[(c(1.0), &[(CRE, 0), (ANN, 1)])]);
        let relabeled = o.relabel_modes(vec![4, 2]).unwrap();
        assert!(relabeled.equiv(&op(&[(c(1.0), &[(CRE, 4), (ANN, 2)])]), 1e-12));

        assert!(matches!(
            o.relabel_modes(vec![4, 4]),
            Err(CoherenceError::DuplicateIndices)
        ));
        assert!(matches!(
            o.relabel_modes(vec![4]),
            Err(CoherenceError::IndexMapTooSmall)
        ));
    }

    #[test]
    fn test_has_and_num_groups() {
        let mut zero = BosonOperator::zero();
        assert!(!zero.has_groups());
        assert!(zero.num_groups().is_none());

        zero.groups = Some(vec![]);
        assert!(zero.has_groups());
        assert_eq!(zero.num_groups(), Some(0));

        let mut one = BosonOperator::one();
        one.groups = Some(vec![0]);
        assert!(one.has_groups());
        assert_eq!(one.num_groups(), Some(1));
    }

    #[test]
    fn test_split_out_groups_binary_search_path_matches_scan() {
        // The same four terms with the same group sizes, once interleaved (so the scan runs) and once
        // group-ordered (so the binary search runs). Both must yield identical per-group operators.
        let terms = [
            (c(1.0), &[(CRE, 0), (ANN, 1)][..]),
            (c(2.0), &[(CRE, 1), (ANN, 0)][..]),
            (c(3.0), &[(CRE, 0), (ANN, 0)][..]),
            (c(4.0), &[(CRE, 1), (ANN, 1)][..]),
        ];
        let mut interleaved = op(&terms);
        interleaved.groups = Some(vec![1, 0, 1, 0]);
        assert!(!interleaved.groups.as_ref().unwrap().is_sorted());

        let ordered = crate::operators::terms::ordering::group::group_order(&interleaved);
        assert!(ordered.groups.as_ref().unwrap().is_sorted());

        let request = [1, 0, 7, 1];
        let scanned = interleaved.split_out_groups(Some(&request)).unwrap();
        let searched = ordered.split_out_groups(Some(&request)).unwrap();
        assert_eq!(scanned.len(), searched.len());
        for (a, b) in scanned.iter().zip(&searched) {
            assert!(a.equiv(b, 1e-12));
        }
        // index 7 is carried by no term in either layout and gives the zero operator
        assert!(searched[2].coeffs.is_empty());
        // the full split agrees as well
        let all_scanned = interleaved.split_out_groups(None).unwrap();
        let all_searched = ordered.split_out_groups(None).unwrap();
        for (a, b) in all_scanned.iter().zip(&all_searched) {
            assert!(a.equiv(b, 1e-12));
        }
    }

    #[test]
    fn test_split_out_groups() {
        let mut o = op(&[
            (c(1.0), &[(CRE, 0), (ANN, 1)]),
            (c(1.0), &[(CRE, 1), (ANN, 0)]),
            (c(2.0), &[(CRE, 0), (ANN, 0)]),
        ]);
        o.groups = Some(vec![0, 0, 1]);

        let groups = o.split_out_groups(None).unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].coeffs.len(), 2);
        assert_eq!(groups[1].coeffs.len(), 1);

        // order preservation and partial selection
        let selected = o.split_out_groups(Some(&[1, 0])).unwrap();
        assert!(selected[0].equiv(&groups[1], 1e-12));
        assert!(selected[1].equiv(&groups[0], 1e-12));

        // duplicates are materialized once per occurrence
        assert_eq!(o.split_out_groups(Some(&[0, 0])).unwrap().len(), 2);
        // empty request -> empty (not None) result
        assert_eq!(o.split_out_groups(Some(&[])).unwrap().len(), 0);

        // an ungrouped operator yields None
        let ungrouped = BosonOperator::one();
        assert!(ungrouped.split_out_groups(None).is_none());
        assert!(ungrouped.split_out_groups(Some(&[0])).is_none());
    }

    #[test]
    #[should_panic(expected = "This method can only be called when groups are present!")]
    fn test_iter_with_groups_panics_without_groups() {
        let o = BosonOperator::zero();
        let _ = o.iter_with_groups();
    }

    #[test]
    fn test_iter_from_terms_round_trip() {
        let o = op(&[
            (c(1.0), &[(CRE, 0), (ANN, 1)]),
            (c(2.0), &[(CRE, 0), (CRE, 0), (ANN, 1), (ANN, 1)]),
        ]);
        assert_eq!(BosonOperator::from_terms(o.iter()), o);
    }

    #[test]
    fn test_iter_with_groups_from_terms_with_groups_round_trip() {
        let mut o = op(&[
            (c(1.0), &[(CRE, 0), (ANN, 1)]),
            (c(2.0), &[(CRE, 0), (ANN, 0)]),
        ]);
        o.groups = Some(vec![0, 1]);
        assert_eq!(
            BosonOperator::from_terms_with_groups(o.iter_with_groups()),
            o
        );
    }

    /// `n_0^2 = b^dag b b^dag b` normal-orders to `b^dag b^dag b b + n_0`.
    #[test]
    fn test_number_operator_square() {
        let n0 = op(&[(c(1.0), &[(CRE, 0), (ANN, 0)])]);
        let squared = n0.__pow__(2).normal_ordered(None);
        let expected = op(&[
            (c(1.0), &[(CRE, 0), (CRE, 0), (ANN, 0), (ANN, 0)]),
            (c(1.0), &[(CRE, 0), (ANN, 0)]),
        ]);
        assert!(squared.equiv(&expected, 1e-12));
    }
}
