<!-- SHIELDS -->
<div align="left">

  [![Release](https://img.shields.io/pypi/v/qiskit-bosons?label=Release&logo=pypi&logoColor=white)](https://pypi.org/project/qiskit-bosons/)
  [![License](https://img.shields.io/github/license/Qiskit/qiskit-bosons?label=License)](LICENSE.txt)
  [![Docs](https://img.shields.io/badge/%F0%9F%93%84%20Docs-dev-blue.svg)](https://qiskit.github.io/qiskit-bosons/)
  [![Qiskit](https://img.shields.io/badge/Qiskit%20-%20%3E%3D2.5.0%2C%3C2.7%20-%20%236133BD?logo=Qiskit)](https://github.com/Qiskit/qiskit)
  [![Python](https://img.shields.io/badge/python-3.11%7C3.12%7C3.13%7C3.14-blue.svg)](https://www.python.org/)
  [![rustc](https://img.shields.io/badge/rustc-1.95+-blue.svg)](https://rust-lang.github.io/rfcs/2495-min-rust-version.html)
  ![Platform](https://img.shields.io/badge/%F0%9F%92%BB%20Platform-Linux%20%7C%20macOS%20%7C%20Windows-informational)
  [![Tests](https://github.com/Qiskit/qiskit-bosons/actions/workflows/test_development_versions.yml/badge.svg)](https://github.com/Qiskit/qiskit-bosons/actions/workflows/test_development_versions.yml)
  [![Coverage](https://coveralls.io/repos/github/Qiskit/qiskit-bosons/badge.svg?branch=main)](https://coveralls.io/github/Qiskit/qiskit-bosons?branch=main)
</div>

# Qiskit Bosons

> [!WARNING]
> This package is under active development!
> If you have feedback, please [open an issue](https://github.com/Qiskit/qiskit-bosons/issues/new/choose).

`qiskit-bosons` extends the Qiskit SDK with tools for working on bosonic systems. Within its scope
it provides the following:

- An efficient data structure for the representation and manipulation of bosonic operators
- A framework for implementing boson-to-qubit encodings, built around shared level encodings (unary
  with or without projectors, standard binary, Gray), and a library of common mappers
- A framework and library of gates for expressing bosonic circuits
- A transpilation pipeline integrated with the Qiskit transpiler process to synthesize bosonic
  circuits into qubit-based circuits

----------------------------------------------------------------------------------------------------

### Documentation

[Documentation](https://quantum.cloud.ibm.com/docs/addons/qiskit-bosons) for this package will be
available on IBM Quantum Platform.

Until then, you can browse the docs of the development branch on
[GitHub Pages](https://qiskit.github.io/qiskit-bosons/dev/index.html). Previews of the development
branch and of labeled pull requests are published there by CI; see
[`docs/README.md`](docs/README.md#previews).

----------------------------------------------------------------------------------------------------

### Installation

Install this package via `pip`, when possible:

```bash
pip install 'qiskit-bosons'
```

For more installation information, refer to these [installation instructions](docs/install.rst).

----------------------------------------------------------------------------------------------------

### Get started

This package covers the whole path from a bosonic operator to a qubit circuit. You can work in the
bosonic regime with both, operator and circuit representation, before you choose a truncation and a
boson-to-qubit encoding that plug into the multi-representation transpiler pipeline of Qiskit.

All of these features are explained in [several guides](docs/guides/) to help you get started.
For an end-to-end overview, start with the
[Bose-Hubbard guide](docs/guides/bose_hubbard.rst).

----------------------------------------------------------------------------------------------------

### Technical discussion

#### Design intentions

This package is deliberately designed to align with Qiskit: it builds on a core implemented in Rust
and provides first-party language bindings to Python. Its API intends to draw parallels to Qiskit in
order to seamlessly integrate into the workflows of users with experience in programming Qiskit.
While there exists no C API right now, the architecture is set up to support one if the need for one
ever arises in the future. The architecture of this package is also deliberately analogous to that
of [`qiskit-fermions`](https://github.com/Qiskit/qiskit-fermions).

A core principle to the design of `qiskit-bosons` was the decoupling of its bosonic circuit
representation from its qubitized form. To be more precise: the bosonic circuits are meaningful by
themselves and do not require a mapping to qubit space to be interpretable. Furthermore, the bosonic
circuit representation cannot make any assumptions about the truncation and boson-to-qubit encoding
applied later on. Both enter only at transpile time.

----------------------------------------------------------------------------------------------------

### Contributing

The source code is available [on GitHub](https://github.com/Qiskit/qiskit-bosons).

The developer guide is located at
[CONTRIBUTING.md](https://github.com/Qiskit/qiskit-bosons/blob/main/CONTRIBUTING.md)
in the root of this project's repository.
By participating, you are expected to uphold Qiskit's
[code of conduct](https://github.com/Qiskit/qiskit/blob/main/CODE_OF_CONDUCT.md).

We use [GitHub issues](https://github.com/Qiskit/qiskit-bosons/issues/new/choose)
for tracking requests and bugs.

----------------------------------------------------------------------------------------------------

### Citing this package

If you use this package in your research, use the [CITATION.bib](CITATION.bib) file in this project’s
repository to cite the appropriate reference(s).

----------------------------------------------------------------------------------------------------

### License

[Apache License 2.0](LICENSE.txt)

----------------------------------------------------------------------------------------------------

### Deprecation policy

> [!IMPORTANT]
> This policy does **not** until the first release is tagged.

This package follows [semantic versioning](https://semver.org/) and is guided by the principles in
[Qiskit's deprecation policy](https://github.com/Qiskit/qiskit/blob/main/DEPRECATION.md). Breaking
changes are made only occasionally, to improve the user experience. When possible, old interfaces are
kept and marked as deprecated for as long as they can co-exist with the new ones. Each substantial
improvement, breaking change, or deprecation is documented in the
[release notes](https://qiskit.github.io/qiskit-bosons/dev/release-notes.html).

<!-- vim: set tw=100: -->
