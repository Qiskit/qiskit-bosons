# ==============================================================================
# Recipes for Linting
# ==============================================================================
.PHONY: lint style

lint:
	cargo metadata --format-version=1 --locked >/dev/null
	cargo fmt --check
	# `core` is linted under the `pyext` FFI backend, which is the only backend
	# wired up today.  The backend features remain mutually exclusive (see
	# `crates/core/src/lib.rs`), so when a `cext` crate is reintroduced it must be
	# clippy-checked as a *separate* compilation unit rather than folded into a
	# workspace-wide `cargo clippy --all-targets`.
	cargo clippy -p qiskit-bosons-core --all-targets --features pyext -- -D warnings
	cargo clippy -p qiskit-bosons-pyext --all-targets -- -D warnings
	tox -e lint
	$(MAKE) pystubs-check

style:
	cargo fmt
	tox -e style

# ==============================================================================
# Recipes for Docs
# ==============================================================================
.PHONY: docs docsclean

docs:
	sphinx-build -W -j auto -T -E --keep-going -b html docs/ docs/_build/html

docsclean:
	rm -rf docs/stubs/ docs/_build

# ==============================================================================
# Recipes for Rust
# ==============================================================================
.PHONY: testrust
testrust:
	cargo test -p qiskit-bosons-core --no-default-features --features pyext

.PHONY: rustcoverage
rustcoverage: export RUSTFLAGS:=-Cinstrument-coverage
rustcoverage: testrust

# ==============================================================================
# Recipes for Python Building
# ==============================================================================

# `pystubs` and `pystubs-dev` are conflicting rules - they both attempt to
# generate the Python stub files, but they differ between release and dev mode.
#
# `stub_gen` only ever *writes* stub files; it never deletes them.  So a stub for a
# module that was renamed or removed on the Rust side would linger.  Both rules
# therefore depend on `pystubs-clean` to wipe the tree first, making a plain
# `make pystubs` a complete regenerate-from-scratch that cannot leave orphans behind.
.PHONY: pystubs pystubs-dev pystubs-clean
pystubs: pystubs-clean
	cargo run --release --bin stub_gen -p qiskit-bosons-pyext --no-default-features

pystubs-dev: pystubs-clean
	cargo run --bin stub_gen -p qiskit-bosons-pyext --no-default-features

pystubs-clean:
	find python/ -name '*.pyi' -delete

# Verify the committed `_lib` stubs are exactly what `pystubs` regenerates: fail on any
# drift or hand-edit.  This is what enforces "committed but always generator-fresh".
# (We deliberately do not run `mypy.stubtest`: on this project it only emits benign
# positional-only dunder notes and PyO3 nested-submodule false positives -- no real
# findings -- so the drift guard here plus the downstream consumer type-check are the
# authoritative gates.)
.PHONY: pystubs-check
pystubs-check: pystubs
	git diff --exit-code -- 'python/qiskit_bosons/_lib/**/*.pyi'

# `pyext` and `pyext-dev` are conflicting rules - they both attempt to "install"
# the compiled Rust acceleration library, but they differ between release and
# dev mode.
.PHONY: pyext pyext-dev
pyext: pystubs
	python setup.py build_rust --inplace --release

pyext-dev: pystubs-dev
	python setup.py build_rust --inplace

# ==============================================================================
# Recipes for Python Testing
# ==============================================================================

.PHONY: testpython
testpython:
	python -m pytest -s -p no:doctest

# Runs every test that needs an optional dependency: the `skipif`-marked ones plus the
# `optionaldep`-marked ones, which `tests/python/conftest.py` marks automatically because a
# module-level `pytest.importorskip` leaves nothing for `-m skipif` to match on.  The optional
# dependencies must actually be installed, or this selects tests that then skip; see the "Optional
# dependencies" section of tests/README.md.  This package declares no optional dependency today, so
# the selection is currently empty (pytest's "no tests ran" exit code 5 is therefore accepted).
.PHONY: testoptional
testoptional:
	python -m pytest -s -p no:doctest -m "skipif or optionaldep"; test $$? -eq 0 -o $$? -eq 5

.PHONY: pycoverage
pycoverage:
	python -m pytest -s -p no:doctest --cov=python/qiskit_bosons/

# ==============================================================================
# Recipes for Documentation Testing
# ==============================================================================

.PHONY: doctest
doctest:
	python -m pytest -s -p no:doctest docs/

# ==============================================================================
# Coverage recipes
# ==============================================================================

.PHONY: coveragecombine
coveragecombine:
	grcov . --binary-path target/debug/ --source-dir . --output-type lcov --output-path rust.info --llvm --branch --parallel --keep-only 'crates/*'
	coverage lcov -o python.info
	# `coverage lcov` emits one function record per definition, so a property's getter
	# and setter yield two records with the same name. lcov >= 2.4 rejects that as
	# inconsistent; downgrade it to a warning (the duplicates get merged into one).
	lcov --ignore-errors inconsistent --add-tracefile python.info --add-tracefile rust.info --output-file coveralls.info

.PHONY: coveragereport
coveragereport: coveragecombine
	genhtml -o htmlcov coveralls.info

.PHONY: coverage
coverage: rustcoverage pycoverage coveragereport

.PHONY: coverageclean
coverageclean:
	rm *.profraw python.info rust.info coveralls.info
	rm -rf htmlcov/

# ==============================================================================
# Recipes for Versioning
# ==============================================================================
.PHONY: version-bump

# Set the project version across Python (VERSION.txt) and Rust (Cargo.toml), then refresh
# Cargo.lock.  Pass the new version in PEP 440 form, e.g.
#
#     make version-bump VERSION=0.2.0.dev0
#     make version-bump VERSION=1.0.0rc1
#     make version-bump VERSION=1.0.0
#
version-bump:
ifndef VERSION
	$(error VERSION is required, e.g. `make version-bump VERSION=0.2.0.dev0`)
endif
	python tools/set_version.py "$(VERSION)"
	# Regenerate Cargo.lock so the workspace crate versions match Cargo.toml.
	cargo metadata --format-version=1 >/dev/null
