# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`qiskit-bosons` extends Qiskit with tools for bosonic systems. Following Qiskit's design philosophy: the **core is written in Rust**, with **first-party Python bindings**, and APIs designed to feel like Qiskit's.

Scope: a bosonic operator data structure, a framework + library for boson-to-qubit encodings, and a framework + library for quantum-circuit synthesis, integrated into Qiskit's transpiler pipeline.

The older [`qiskit-fermions`](https://github.com/Qiskit/qiskit-fermions) follows the same architecture and release pipeline. We deliberately keep the two packages aligned where it makes sense, so it is often the right reference when adding something new here.

## Repository layout note

This checkout is a **bare repo + worktrees** setup. The project source lives in named worktree directories (`main/`, `bosons/`, …), *not* in the top-level directory which only holds `.bare/` and the worktree dirs. Run commands and edit files inside a worktree.

## Crate architecture

The Rust workspace (`crates/`) has two crates whose **module structure mirrors the public Python API** for easy navigation:

- **`core`** — all native logic (operators, mappers, linalg, random). **Must NOT depend on `pyo3`.** This keeps the Rust core cleanly separable from Python so a future C API can be Python-independent. Its FFI backend is feature-gated.
- **`pyext`** — the public **Python** API built on `core` via `pyo3`. Exposes a `#[pymodule] _lib` whose nested modules are registered into `sys.modules` at import time by `python/qiskit_bosons/__init__.py`. Depends on `core` with its `pyext` feature and links Qiskit through `qiskit-pyo3-ffi` (no locally compiled Qiskit binary required — this is what unlocks Windows support and keeps CI from building Qiskit's C API).

One vendored support crate:
- **`qiskit-pyo3-ffi`** — Qiskit's PyO3 FFI crate, vendored until published. Excluded from lint/typos checks; don't hand-edit.

### The `cext` seam is preserved but unused

There is **no `cext` crate today** (the C API is future work), but the *seam* it plugs into is deliberately intact:

- Every FFI-touching module in `core` selects its backend via `cfg_select!` on the `pyext` / `cext` features rather than importing a binding crate directly.
- `core/Cargo.toml` declares an **empty placeholder `cext` feature**, so the `cfg_select!` arms and the mutual-exclusion `compile_error!` in `crates/core/src/lib.rs` keep compiling without tripping rustc's `unexpected_cfgs` lint.
- The two backends must stay **mutually exclusive** — they bind the same `ffi` alias to different crates.

So reintroducing a C API means wiring `cext = ["dep:qiskit-sys"]` and adding the crate, **not** untangling Python from the core. Keep it that way: don't let `pyo3` types leak into `core`.

### Adding a feature (the intended pattern)

1. Implement + test in the Rust `core` crate (minimal internal docs; Rust API is not public).
2. Expose via `pyext` with tests in `tests/python/` and full API docs (docstrings live in the `pyext` Rust source).

Features may also be implemented in the Python API directly first and refactored down to `core` later.

## Build & test

**Always use the Makefile targets, not raw `cargo`.**

| Task | Command |
|------|---------|
| Build Python extension (release) | `make pyext` |
| Build Python extension (dev) | `make pyext-dev` |
| Regenerate `.pyi` stubs | `make pystubs` / `make pystubs-dev` |
| Run Rust tests | `make testrust` (core crate, `pyext` backend) |
| Run Python tests | `make testpython` |
| Run optional-dependency tests | `make testoptional` (currently selects nothing; no optional deps yet) |
| Run guide doctests | `make doctest` (every `docs/**/*.rst`; in practice the guides are what carry them) |
| Lint (rust + python) | `make lint` |
| Auto-fix formatting | `make style` |
| Build docs | `make docs` (Sphinx runs with `-W`: any warning fails the build) |
| Bump version | `make version-bump VERSION=0.2.0.dev0` |

> **`make testrust` needs a Python interpreter.** `core` is tested under the **`pyext`** backend (the only one wired up), which links `libpython`. A bare `cargo test` on FFI-touching modules will abort with a non-unwinding panic; that is why CI runs `testrust` *after* the package install. There are no `QISKIT_LIB`/`QISKIT_INCLUDE` variables in this project — nothing links Qiskit's C library.

**Running a single Python test:**
```bash
python -m pytest -s -p no:doctest tests/python/operators/test_boson_operator.py::TestName::test_method
```
(`-p no:doctest` because doctests are collected separately via Sybil — see below.)

Editable install: `SETUPTOOLS_RUST_CARGO_PROFILE=release pip install -e .` (the `all` extra exists but is an empty placeholder, and CI installs `.[all]`; there are no real optional-dependency extras yet).

### Important build coupling

`make testpython` is what actually *runs* the Python doctests (the docstring/`.pyi` doctests under `python/qiskit_bosons/`, collected via Sybil during a normal `pytest` run), but it does **not** build anything itself — it is a bare `pytest` invocation with no prerequisites. The stubs carry the API docstrings *and doctests*, so they must be regenerated for those doctests to reflect current source: run **`make pyext-dev` first** (it depends on `pystubs-dev`, so it recompiles `pyext` *and* regenerates the `.pyi` stubs). Skipping it silently tests a stale `_lib*.so` against stale stubs. CI sidesteps this by running `make pystubs testpython` after the package install.

> **Type stubs (`.pyi`) — invariants.** Generated by `pyo3-stub-gen`, rooted at the physical `qiskit_bosons._lib` module, **committed** to git, and packaged into the wheel. Never hand-edit a `.pyi` — change the Rust docstrings/signatures in `crates/pyext/src` and run `make pystubs`. `make lint` runs `pystubs-check` (regenerate + `git diff --exit-code`), so stale or hand-edited stubs fail CI. Pure-Python modules are typed **inline**. Full rationale in [`tests/README.md`](tests/README.md).

> **Doctest gotcha.** A `.. doctest::` block must **not** end with its last output line immediately followed by the closing `"""` at a shallower indent — the doctest parser reads the quotes as expected output and errors with "inconsistent leading whitespace". End the docstring with a line of prose after the block.

## Docs & docstrings — where they live

- **Python API docs** come from `python/qiskit_bosons/` source and, for `pyext`-implemented components, from **docstrings written in the `pyext` Rust source**, surfaced through the `.pyi` stubs. Layout config is in `docs/pydoc/`. Docstring convention: **Google** (`ruff` pydocstyle).
- **Doctests** are collected by **Sybil** (`python/conftest.py`, `docs/conftest.py`) from `*.py`/`*.pyi`/`*.rst`, with `ELLIPSIS | NORMALIZE_WHITESPACE | NUMBER` flags. Docstring and `.pyi` doctests run via `make testpython`; `make doctest` collects every `docs/**/*.rst`, which in practice means the guides.
- **Guides** live in `docs/guides/`. They use `.. plot::` with `:context:` `:nofigs:` `:include-source:` (not `.. doctest::`) so Sphinx and Sybil read the same `>>>` blocks. Conventions worth matching: cast to `bool(...)` so doctests print `True` rather than `np.True_`; a hidden guard block is a `.. plot::` with `:context:` `:nofigs:` and *no* `:include-source:`; close matplotlib figures with `plt.close("all")` at the end or the next guide's figures get corrupted.

## Testing structure

- **Rust unit tests** live inline in the crate source; run via `make testrust`.
- **Python tests** in `tests/python/` (mirrors the API tree: `operators/`, `mappers/`, `circuit/`, `transpiler/`).
- Keep tests **laptop-scale**: 3 modes and `d <= 3` for anything building dense matrices. The guide doctests run in ~1.5s and the docs build in a few minutes; keep it that way.
- **Optional-dependency tests are guarded, and the `test` dependency group deliberately omits those deps.** This package declares none yet, but the machinery is in place: `make testpython` from `--group test` alone runs with them absent, which is what proves the guards hold (a *failure* there means a missing guard). `make testoptional` is its exact complement, selecting every guarded test once the extras are installed. `make doctest` is deliberately different: the guides under `docs/` are documentation, so it always runs *with* the optional deps. `tests/python/conftest.py` infers the `optionaldep` marker from the guards themselves (never apply it by hand). Read [`tests/README.md`](tests/README.md) before changing the group split.

## Transpiler / synthesis plugins

Synthesis comes in two stages, each exposed as **stevedore entry points** (see `pyproject.toml` `[project.entry-points]`) and each driven by a pass with the same `config` scheme (`"Method"`, `("Method", args)`, `("Method", args, kwargs)` or an instance):

- **boson-to-boson** (optional, in bosonic space): group `qiskit_bosons.circuit.synthesis`, named `<Gate>.<Method>`: `Evolution.LieTrotter` / `SuzukiTrotter` / `QDrift`, `ModeRotation.Clements`. A plugin is any object with `synthesize(gate) -> BosonicCircuit` (the `B2BSynthesisPlugin` protocol, no base class), managed by `B2BSynthesisPluginManager` and applied by the `B2BSynthesis` pass (`transpiler/passes/optimization/`). Implementations are in `python/qiskit_bosons/circuit/library/synthesis/`.
- **boson-to-qubit**: group `qiskit_bosons.transpiler.synthesis`, e.g. `Evolution.MapperFn`, `InitializeModes.Occupation`. Pass implementations are in `python/qiskit_bosons/transpiler/passes/synthesis/`.

After editing `pyproject.toml` entry points, re-run the editable install, or `B2BSynthesisPluginManager` will not see them.

One thing to know first:

- **The encoding and truncation enter here, at transpile time**, not in the circuit. That is the package's central design point. A bosonic circuit is encoding-agnostic; `OccupationB2QLayout` + the synthesis plugins are what pin it down.

> **Grouped evolutions and product formulas.** Each group of a grouped operator becomes its own `PauliEvolutionGate`, so a `product_formula` handed to `MapperFnEvolutionSynthesis` subdivides terms *within* one group and never interleaves the groups. Since groups do not generally commute, raising its repetition count does **not** reduce the outer Trotter error. To converge the *outer* error, Trotterize in **bosonic** space, where the grouping still exists: give the `Evolution` gate a `synthesis=` (`BosonicLieTrotter(reps)` / `BosonicSuzukiTrotter(order, reps)`) or apply one to a whole circuit with `B2BSynthesis({"Evolution": ...})`. Also note each group must be **Hermitian** on its own — a lone hopping term `b†_i b_j` is not, and `PauliEvolutionGate` rejects the resulting complex coefficients; `groups_are_hermitian` checks it.

> **Mode rotations: the preset decomposes them, the framework does not require it.** `ModeRotation(synthesis=None)` is left whole (like an atomic `Evolution`). `BosonicClementsDecomposition` (`circuit/library/synthesis/`) splits it into atomic, single-group `Evolution` gates of the Givens mesh, and the preset's `optimization` stage applies it with `B2BSynthesis({"ModeRotation": "Clements"})`, whose method wins over the gate's own, as for an `Evolution`. The library ships **no** `ModeRotation` boson-to-qubit plugin, but any plugin registered under that name works for a gate left whole (the guide shows a ten-line one that maps `to_operator()` once); `B2QSynthesis` names both remedies if neither exists. The mesh is exact only on the untruncated Fock space, and a cutoff of the total boson number keeps every factor within it; the direct path truncates the whole generator once and so equals `occupation_encoding` at any cutoff.

> **Neighboring evolutions of one operator are fused, but not by `equiv`.** `BosonicSuzukiTrotter(merge_adjacent=True)` (the default) fuses the same-factor seams of its own schedule, and the `MergeEvolutions` pass does it across gates. The pass compares operators **without** normal ordering (`(a - b).simplify()` must vanish), not with `BosonOperator.equiv`: the truncation does not commute with normal ordering, so `b b†` and `b†b + 1` are `equiv` but map to different matrices, and merging them changes the circuit. Keep that distinction in any new pass that compares operators.

> **Iterate a DAG with `topological_op_nodes()`, not `op_nodes()`.** `op_nodes()` returns the nodes in *storage* order, which a pass that substitutes a node by a sub-circuit (`B2BSynthesis`) leaves out of sync with the wire dependencies. A pass that rebuilds its output node by node (`B2QSynthesis`, `RelabelModes`) must therefore iterate in circuit order, or the factors of successive substitutions interleave: a well-formed circuit for a different unitary, because gates on shared qubits do not commute. Substituting in place is safe with either; `tests/python/transpiler/passes/test_dag_order.py` pins the rebuilding passes.

> **Atomic evolutions.** The factors a synthesis (`BosonicQDrift` included) emits are `Evolution(..., atomic=True)` and have *no* definition, so `decompose()` leaves them alone. This is load-bearing, not cosmetic: splitting a Hermitian group further, term by term, yields non-Hermitian factors whose exponentials are not unitary. `synthesis` and `atomic` are read-only attributes of every `BosonicGate` because the definition is cached; `with_synthesis(...)` returns a copy.

> **The mapper merges its own result, with a relative tolerance.** `occupation_encoding` is driven by `crates/core/src/mappers/qk_obs.rs` (parallel, each term canonicalized before it is accumulated, a final unconditional merge), and `MapperFnEvolutionSynthesis` passes the mapper's output on untouched so its term order survives. The tolerance of the *final* merge (`RELATIVE_FINAL_TOL`) is relative to the largest coefficient, not absolute. The occupation encodings' cancellations are not exact (Gray in particular leaves rounding noise of ~`1e-16` times the coefficient size), that noise is *imaginary*, and `PauliEvolutionGate` rejects any imaginary part. No fixed absolute tolerance works at every scale. Do not replace it with one without re-reading the developer note there. A custom mapper must itself return real coefficients; `simplify(...)` pins a canonical order and `group_wise(...)` emits each group contiguously.
>
> **Never convert a mapped operator with `SparsePauliOp.from_sparse_observable` outside a tiny test.** It expands every projector into two Paulis, so a `unary_projected` diagonal term over `w` qubits becomes `2**w` terms (this exhausted memory once), and a single `b†` at `d` levels expands to `2**(d-2)` times as many Paulis as it holds terms (measured: 64x at `d = 8`). `unary_unprojected` drops the `|0><0|` projectors and only its `|1><1|` factors expand, so the expansion stays about linear in `d` (measured: 28 terms stored and 28 expanded for `b†` at `d = 8`), but the compact encodings and any operator with many `|1><1|` factors can still blow up. Read `SparseObservable.coeffs` directly.

## Versioning

`VERSION.txt` (committed) is the source of truth for Python; `make version-bump VERSION=...` syncs it with the Rust `Cargo.toml` workspace version and refreshes `Cargo.lock` via `tools/set_version.py`. `version.py` appends a git sha for dev builds and is excluded from coverage. Use PEP 440 forms (`0.2.0.dev0`, `1.0.0rc1`, `1.0.0`).

**Pre-release status:** there is no tagged release yet, so `main` is an **unstable development branch**. Until the first tagged release lands, we don't track release notes — breaking changes can be made freely and don't need `reno` entries.

**Releasing and CI:** the release pipeline (wheels via `cibuildwheel`, `release.yml`, `stable/X.Y` branches with Mergify backports, GitHub Pages docs previews) is described in [`docs/release-process.md`](docs/release-process.md). All test/lint/docs workflows are reusable `workflow_call` units driven by `branch-protection.yml` (whose single `Finalize` job is the only required check) and `on-nightly.yml`. Qiskit is bounded to `qiskit>=2.5.0,<2.7`: the vendored `qiskit-pyo3-ffi` snapshot must come from the **oldest** supported minor (Qiskit only appends C API slots), and C API functions newer than the lower bound must not be used (see the crate's README).

## Conventions

- **Every source file requires a copyright header** — `tools/verify_headers.py` enforces this in `make lint` (run over `crates docs python tests tools`, ignoring vendored `crates/qiskit-pyo3-ffi`).
- Python: `ruff` (line-length 100, target py311, Google docstrings), `mypy` (py311).
- Rust: `cargo fmt`, edition 2024, MSRV 1.95.
- Release notes: managed via `reno` (`releasenotes/`), but **not tracked yet**.
- **Commutator argument order:** `commutator(a, b)` is the mathematical `[A, B]`. It composes with `__matmul__` (right multiplication, `A B`) and **not** `__and__` — following Qiskit, `__and__` is `compose`, i.e. *left* multiplication: `A & B == B @ A`. Building the brackets out of it silently negates every commutator.
