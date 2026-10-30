# Rust crates

This document explains the structure and roles of the various Rust crates that are being developed
in this package.

## `core`

This crate implements the core functionality of this entire package in Rust.
Its module structure mirrors that of the public Python API.

Below are some guidelines that should be adhered to when working on this crate:

- This crate shall **NOT** depend on `pyo3`. This enforces a clean separation of the Rust core from
  the Python API and will allow the C API built on top of this core to be fully independent from
  Python (in the future, once the Qiskit Rust core fulfills the same criterion).
- As long as we do not plan a public Rust API, this crate can depend on Qiskit's FFI crates. If we
  desire to publish a Rust crate in the future, this dependency is problematic: `qiskit-pyo3-ffi` is
  vendored in this repository precisely because it is not published yet, and a published crate
  cannot depend on an unpublished one. If that time comes, we may choose to separate those parts
  of this crate that depend on Qiskit out into a separate crate. That should leave us with a core data
  model crate in Rust, whose usage may be useful to other Rust packages (e.g. `ffsim`).

The Qiskit FFI backend is feature-gated, because a `pyext` and a `cext` binding crate require
*different, mutually exclusive* backends:

- The `pyext` feature pulls in `qiskit-pyo3-ffi` (Qiskit's PyO3 FFI crate). **This is the only
  backend wired up today.**
- The `cext` feature would pull in `qiskit-sys` (raw bindings to the compiled `libqiskit`). It is
  declared as an **empty placeholder**: no such dependency exists yet, but keeping the feature
  declared lets the `cfg_select!` arms and the mutual-exclusion guard compile without tripping
  rustc's `unexpected_cfgs` lint.

The default feature set enables neither. Enabling both at once is a hard error — a `compile_error!`
in `src/lib.rs` guards against it, since Cargo's feature unification would otherwise silently
compile the crate against the wrong backend.

Every FFI-touching module selects its backend through `cfg_select!` on these features rather than
importing a binding crate directly. That is the seam a future C API plugs into: reintroducing one
means wiring `cext = ["dep:qiskit-sys"]` and adding the crate, **not** untangling Python from the
core. Keep `pyo3` types out of `core` so it stays that way.

## `pyext`

This crate defines the public Python API on top of the `core` crate.
As such, it obviously depends on the `pyo3` crate.
It depends on `core` with its `pyext` feature and links Qiskit through the `qiskit-pyo3-ffi` crate.
Using `qiskit-pyo3-ffi` (rather than linking against a locally compiled Qiskit acceleration binary)
removes the hard-coded link step and is what unlocks Windows support for the Python package.

Just like the `core` crate, the module structure mirrors that of the public Python API and we
leverage `pyo3`'s modern `#[pymodule]` and `#[pymodule_export]` macros to implement the same Python
module layout. Note, that we try to automatically register these nested Python modules in
`sys.modules` during the loading of `python/qiskit_bosons/__init__.py`.

Additionally, the `pyo3-stub-gen` crate is leveraged to generate
[Python stub files](https://typing.python.org/en/latest/guides/writing_stubs.html) to provide proper
type hinting to the Python package as well as to extract the Python API documentation and doctests
from the `pyext` docstrings.

## `cext` (future work)

There is no `cext` crate yet. When one is added it will define the public C API, mirroring the module
structure for ease of navigation even though C has no modules, depend on `core` with its `cext`
feature, and link Qiskit's compiled core through a `qiskit-sys` crate. Its header would be generated
with `cbindgen` during the build.

The `core` crate is written to accommodate this without rework — see the feature-gating note above.

<!-- vim: set tw=100: -->
