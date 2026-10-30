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

//! Shared machinery for mappers that emit a Qiskit `SparseObservable`.
//!
//! A mapper's choice of encoding decides how a single term of an operator turns into Pauli strings.
//! Everything wrapped around that choice (collecting the strings of a term into one observable and
//! accumulating the terms in parallel under a memory bound) is identical whatever the encoding, and
//! lives here.

use rayon::prelude::*;
use std::sync::{Arc, Mutex};

cfg_select! {
    feature = "pyext" => {
        pub(crate) use num_complex::Complex64 as QkComplex64;
        extern crate qiskit_pyo3_ffi as ffi;
        pub(crate) use qiskit_pyo3_ffi as ffi_backend;
        pub(crate) use ffi::QkBitTerm::One as QK_BIT_TERM_ONE;
        pub(crate) use ffi::QkBitTerm::X as QK_BIT_TERM_X;
        pub(crate) use ffi::QkBitTerm::Y as QK_BIT_TERM_Y;
        pub(crate) use ffi::QkBitTerm::Zero as QK_BIT_TERM_ZERO;
    }
    feature = "cext" => {
        extern crate qiskit_sys as ffi;
        pub(crate) use ffi::QkComplex64 as QkComplex64;
        pub(crate) use qiskit_sys as ffi_backend;
        pub(crate) use ffi::QkBitTerm_QkBitTerm_One as QK_BIT_TERM_ONE;
        pub(crate) use ffi::QkBitTerm_QkBitTerm_X as QK_BIT_TERM_X;
        pub(crate) use ffi::QkBitTerm_QkBitTerm_Y as QK_BIT_TERM_Y;
        pub(crate) use ffi::QkBitTerm_QkBitTerm_Zero as QK_BIT_TERM_ZERO;
    }
}

/// Builds one observable holding several Pauli/projector terms at once.
///
/// This is the one place a mapper calls `qk_obs_new`. Handing a whole term's strings over in a single
/// call keeps a term down to one FFI object, rather than one per string. The buffers are borrowed:
/// `qk_obs_new` copies them into the observable it builds and never writes through them.
pub(crate) fn pauli_sum(
    num_qubits: u32,
    coeffs: &[QkComplex64],
    bit_terms: &[ffi::QkBitTerm],
    indices: &[u32],
    boundaries: &[usize],
) -> *mut ffi::QkObs {
    debug_assert_eq!(indices.len(), bit_terms.len());
    debug_assert_eq!(boundaries.len(), coeffs.len() + 1);
    debug_assert_eq!(boundaries.last().copied(), Some(bit_terms.len()));
    // SAFETY: `qk_obs_new` only reads its four buffers, copying them into the observable it returns,
    // so casting away their constness never leads to a write through a shared borrow.
    unsafe {
        ffi::qk_obs_new(
            num_qubits,
            coeffs.len().try_into().unwrap(),
            bit_terms.len().try_into().unwrap(),
            coeffs.as_ptr().cast_mut(),
            bit_terms.as_ptr().cast_mut(),
            indices.as_ptr().cast_mut(),
            boundaries.as_ptr().cast_mut(),
        )
    }
}

/// Tolerance used for every `qk_obs_canonicalize` call that happens *while* the operator is built.
///
/// This is deliberately far below any physically meaningful coefficient: canonicalization is used
/// here only to *merge* duplicate Pauli terms, never to discard physics. Raising it would silently
/// change the mapped operator.
const CANONICALIZE_TOL: f64 = 1e-18;

/// Relative tolerance of the final merge, as a fraction of the largest coefficient magnitude.
///
/// DEVELOPER NOTE: the final merge cannot use the absolute [`CANONICALIZE_TOL`]. Do not replace the
/// relative tolerance with a fixed absolute one.
///
/// Cancellations between the Pauli images of a Hermitian operator are exact in an encoding such as
/// Jordan-Wigner: the imaginary parts of `(X - iY)/2` and `(X + iY)/2` are `+-0.5`, so they add to
/// exactly zero, and an absolute tolerance far below any coefficient would suffice to drop them. The
/// occupation encodings do not have that property. A matrix element of the truncated ladder operator
/// is `sqrt(n)`, which is not exactly representable, and for the compact encodings (Gray in
/// particular) the products that must cancel against each other reach the same Pauli string through
/// different orders of floating-point operations. What is left over is rounding noise of about
/// `1e-16` times the size of the coefficients involved.
///
/// That noise matters because it is *imaginary*. A Hermitian group maps to a real operator only once
/// the terms that cancel have been merged away, and `PauliEvolutionGate` rejects an operator that
/// carries any imaginary part at all, however small. Measured on a hopping pair, a merge at the
/// absolute `1e-18` leaves `5.6e-17` behind at unit coefficients, `1.7e-13` at `1e3` and `5.8e-11` at
/// `1e6`: the residue scales with the coefficients, so no fixed absolute tolerance is right at every
/// scale. It has to be below every real coefficient yet above the noise, and that window closes as
/// the coefficients grow.
///
/// Hence a tolerance relative to the largest coefficient magnitude. `1e-15` sits about a factor of
/// ten above the noise that was measured (`~1e-16` relative) and many orders of magnitude below any
/// coefficient that carries physics. It is applied to the *final* result only: merging while the
/// operator is built stays at the absolute tolerance, so no partial sum is ever rounded away.
const RELATIVE_FINAL_TOL: f64 = 1e-15;

/// Weight a single bit term carries in [`Wrapper::cost`]: a `QkBitTerm` plus its `u32` qubit index.
const WEIGHT_PER_BIT_TERM: usize = size_of::<u32>() + 1;

/// Weight a single Pauli term carries on top of its bit terms: a `Complex64` coefficient plus the
/// `usize` boundary delimiting it.
///
/// Counted separately from [`WEIGHT_PER_BIT_TERM`] because it is charged per *term* rather than per
/// bit term, and the two cannot be folded into one figure without assuming an average Pauli weight.
/// It also dominates for low-weight terms -- and an identity term has no bit terms at all, so a
/// measure built on bit terms alone would score it as free. That this weight is *non-zero* is the
/// load-bearing property; the particular figures only set the relative cost of a term against a bit
/// term.
const WEIGHT_PER_TERM: usize = size_of::<f64>() * 2 + size_of::<usize>();

/// Smallest accumulator size worth compacting at all, in the units of [`Wrapper::cost`].
///
/// Below this the duplication cannot amount to enough memory to be worth a merge pass, so small
/// operators are mapped without ever being compacted.
const MIN_COMPACTION_FLOOR: usize = (1 << 16) * WEIGHT_PER_BIT_TERM;

/// How far an accumulator may grow past its last merged size before it is merged again.
///
/// This is what bounds the accumulators: a merge leaves behind exactly the distinct content, so the
/// next merge is triggered at this multiple of it and the peak is `factor * distinct` -- proportional
/// to the operator being built, whatever the duplication happens to be.
///
/// The trigger is recomputed from the *measured* post-merge size every cycle, so it tracks the
/// distinct content rather than compounding against it.
const GROWTH_FACTOR: f64 = 1.5;

/// Canonicalizes `obs` at tolerance `tol`, freeing the original, to merge duplicate Pauli terms.
///
/// Takes ownership of the pointer it is given and returns a new one; the input must not be used
/// afterwards.
///
/// # Safety
///
/// `obs` must be a valid, uniquely-owned `QkObs` pointer.
unsafe fn compact_with(obs: *mut ffi::QkObs, tol: f64) -> *mut ffi::QkObs {
    let compacted = unsafe { ffi::qk_obs_canonicalize(obs, tol) };
    unsafe { ffi::qk_obs_free(obs) };
    compacted
}

/// Returns the largest coefficient magnitude of `obs`, or `0.0` if it holds no terms.
///
/// # Safety
///
/// `obs` must be a valid `QkObs` pointer.
unsafe fn largest_coeff_magnitude(obs: *mut ffi::QkObs) -> f64 {
    let num_terms = unsafe { ffi::qk_obs_num_terms(obs) };
    if num_terms == 0 {
        return 0.0;
    }
    let coeffs = unsafe { ffi::qk_obs_coeffs(obs) };
    // SAFETY: `qk_obs_coeffs` points at `qk_obs_num_terms` contiguous coefficients, which stay valid
    // for as long as `obs` is neither modified nor freed; both are respected for the read below.
    let coeffs = unsafe { std::slice::from_raw_parts(coeffs, num_terms) };
    coeffs
        .iter()
        .map(|coeff| coeff.re.hypot(coeff.im))
        .fold(0.0, f64::max)
}

// NOTE: https://stackoverflow.com/a/50341075
struct Wrapper {
    ptr: *mut ffi::QkObs,
    /// [`Wrapper::cost`] of this accumulator just after it was last compacted, or 0 if it never was.
    ///
    /// Compaction is triggered by growth relative to this rather than by an absolute size, so that
    /// an accumulator which is already mostly distinct terms is not compacted repeatedly to no
    /// effect.
    compacted_cost: usize,
}
unsafe impl Send for Wrapper {}

impl Wrapper {
    /// Creates a new, empty accumulator.
    fn zero(num_qubits: u32) -> Self {
        Self {
            ptr: unsafe { ffi::qk_obs_zero(num_qubits) },
            compacted_cost: 0,
        }
    }

    /// Returns a weighted size for this accumulator's term buffers.
    ///
    /// Both triggers are driven by this rather than by `qk_obs_len` alone. A term carries its
    /// coefficient and boundary on top of its bit terms, and that part is charged per term: an
    /// *identity* term has no bit terms whatsoever, so a measure built on `qk_obs_len` alone scores
    /// it as free and neither trigger ever fires however many of them pile up.
    fn cost(&self) -> usize {
        let num_bit_terms = unsafe { ffi::qk_obs_len(self.ptr) };
        let num_terms = unsafe { ffi::qk_obs_num_terms(self.ptr) };
        num_bit_terms * WEIGHT_PER_BIT_TERM + num_terms * WEIGHT_PER_TERM
    }

    /// Merges duplicate terms once the accumulator has grown enough past its last merged size to be
    /// worth the work.
    fn compact_if_grown(&mut self) {
        let cost = self.cost();

        let grown = (self.compacted_cost as f64 * GROWTH_FACTOR) as usize;
        if cost <= MIN_COMPACTION_FLOOR.max(grown) {
            return;
        }

        self.ptr = unsafe { compact_with(self.ptr, CANONICALIZE_TOL) };
        self.compacted_cost = self.cost();
    }

    /// Merges duplicate terms unconditionally, at a tolerance relative to the largest coefficient.
    ///
    /// Used for the result that is handed back, so that its size does not depend on whether the last
    /// combine happened to trip the growth trigger, and so that it is free of rounding noise (see
    /// [`RELATIVE_FINAL_TOL`]).
    ///
    /// DEVELOPER NOTE: this does *not* return early when nothing was added since the last
    /// compaction. That shortcut is only sound when the final merge runs at the same tolerance as the
    /// compactions before it, because then there is provably nothing left to merge. Here the final
    /// tolerance is deliberately coarser, so an accumulator that was compacted at the absolute
    /// tolerance still holds noise that only this pass removes.
    fn compact_final(&mut self) {
        let largest = unsafe { largest_coeff_magnitude(self.ptr) };
        let tol = (RELATIVE_FINAL_TOL * largest).max(CANONICALIZE_TOL);

        self.ptr = unsafe { compact_with(self.ptr, tol) };
        self.compacted_cost = self.cost();
    }
}

/// Maps every term of an operator onto a `QkObs` in parallel, merging duplicates as it goes.
///
/// This is the single copy of the mapping driver shared by every mapper that emits a `QkObs`: the
/// thread pool, the per-worker accumulators and the compaction schedule live here and nowhere else,
/// so a fix to any of them lands once.
///
/// Each term's own image is merged *before* it is added to an accumulator. Without that a single
/// term whose encoding produces many duplicate strings would be held in full next to the
/// accumulator, and the memory a term can cost would be bounded by nothing but the term itself.
///
/// `map_term` returns the term's *unscaled* image together with the coefficient to apply, so that the
/// scaling stays in the one `qk_obs_scaled_add_inplace` below instead of costing an extra observable
/// per term. The returned pointer is consumed (freed) here.
///
/// # Memory
///
/// Each worker thread accumulates into an observable of its own, and terms are handed to whichever
/// thread is free rather than partitioned by the Pauli strings they produce, so every accumulator
/// ends up holding roughly a full copy of the mapped operator. Peak memory therefore grows with the
/// number of threads. `RAYON_NUM_THREADS` is the control for it.
///
/// # Panics
///
/// `map_term` must not unwind: it runs inside `for_each`, and a panic there would leak every
/// accumulator.
pub(crate) fn map_operator<T, I, F>(terms: I, num_qubits: u32, map_term: F) -> *mut ffi::QkObs
where
    I: Iterator<Item = T> + Send,
    T: Send,
    F: Fn(&T, u32) -> (*mut ffi::QkObs, QkComplex64) + Sync,
{
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(0)
        .build()
        .unwrap();

    let mut qubit_ops = vec![];
    for _ in 0..pool.current_num_threads() {
        qubit_ops.push(Arc::new(Mutex::new(Wrapper::zero(num_qubits))));
    }

    pool.install(|| {
        terms.par_bridge().for_each(|term| {
            let (mapped_term, qk_coeff) = map_term(&term, num_qubits);

            let canon_term = unsafe { compact_with(mapped_term, CANONICALIZE_TOL) };

            // NOTE: `mut` so that the accumulator can be swapped for its compacted self below.
            // This relies on there being exactly one accumulator per worker thread, which is what
            // makes the `lock()` uncontended: `par_bridge` only ever runs this closure on a pool
            // worker, so `current_thread_index()` is both `Some` and unique per concurrent call.
            let mut qubit_op = qubit_ops[pool.current_thread_index().unwrap()]
                // this should never lock because we have one item per thread
                .lock()
                .unwrap();

            unsafe { ffi::qk_obs_scaled_add_inplace(qubit_op.ptr, canon_term, &qk_coeff) };

            unsafe { ffi::qk_obs_free(canon_term) };

            // Merge the duplicates accumulated so far. Without this the accumulator grows with the
            // number of emitted terms instead of the number of distinct ones.
            qubit_op.compact_if_grown();
        });
    });

    let mapped_operator: Wrapper = qubit_ops
        .par_iter()
        .fold(|| Wrapper::zero(num_qubits), {
            |mut op1: Wrapper, op2| {
                let op_locked = op2.lock().unwrap();
                unsafe { ffi::qk_obs_add_inplace(op1.ptr, op_locked.ptr) };
                unsafe { ffi::qk_obs_free(op_locked.ptr) };
                // Compact as we go: `qk_obs_add_inplace` concatenates, so folding the per-thread
                // accumulators would otherwise hold their combined size at once.
                op1.compact_if_grown();
                op1
            }
        })
        .reduce(|| Wrapper::zero(num_qubits), {
            |op1, op2| {
                let num_add_terms1 = unsafe { ffi::qk_obs_num_terms(op1.ptr) };
                let num_add_terms2 = unsafe { ffi::qk_obs_num_terms(op2.ptr) };
                // Add into whichever side is larger so the bigger buffer is never copied, then
                // merge the duplicates the concatenation introduced.
                let mut acc = if num_add_terms1 > num_add_terms2 {
                    unsafe { ffi::qk_obs_add_inplace(op1.ptr, op2.ptr) };
                    unsafe { ffi::qk_obs_free(op2.ptr) };
                    op1
                } else {
                    unsafe { ffi::qk_obs_add_inplace(op2.ptr, op1.ptr) };
                    unsafe { ffi::qk_obs_free(op1.ptr) };
                    op2
                };
                acc.compact_if_grown();
                acc
            }
        });

    // Merge the result unconditionally. The combine above stops as soon as the growth trigger is
    // satisfied, which would otherwise leave the size of the returned operator dependent on how the
    // work happened to be split -- the caller would see very different term counts for the same
    // operator from one run to the next.
    let mut mapped_operator = mapped_operator;
    mapped_operator.compact_final();

    mapped_operator.ptr
}

/// Returns the `QkComplex64` mirror of a complex number.
pub(crate) fn qk_coeff(coeff: num_complex::Complex64) -> QkComplex64 {
    QkComplex64 {
        re: coeff.re,
        im: coeff.im,
    }
}
