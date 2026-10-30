# Overview

This page summarizes the guides that are available in addition to the extensive API documentation of
this package.

To get started with the breadth of features this package supports, check out this end-to-end guide:

- [Simulate Bose-Hubbard dynamics with a choice of encoding](bose_hubbard.html.md#bose-hubbard-getting-started)

The following sections outline additional guides that dive into more specific topics in detail.

## Operators and their structure

A bosonic operator is a sum of ordered products of creation and annihilation operators, and a mode
can appear more than once in a term. Bosonic operators can also tag their terms with group indices,
which later stages use to split a time evolution along physically meaningful pieces. Mappers turn
an operator into a qubit operator under a truncation and an encoding. Start here if you are
building a Hamiltonian.

- [Design principles of operator representations](operators.html.md#operators-explanation)
- [Group operator terms: use the operator structure](grouping.html.md#grouping-explanation)
- [Map operator representations](mappers.html.md#mappers-explanation)

## Bosonic circuits

You can express circuits directly on bosonic modes, without choosing a truncation or an encoding.

- [Work with bosonic circuits](circuit.html.md#bosonic-circuit-explanation)
- [Rotate bosonic modes](mode_rotation.html.md#mode-rotation-explanation)

## Transpiling bosonic circuits

Guides explaining how a [`BosonicCircuit`](../stubs/qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) becomes a [`QuantumCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.QuantumCircuit),
and how to control the Trotter error along the way.

- [Transpile bosonic circuits](transpilation.html.md#transpilation-explanation)
- [Synthesize an Evolution gate in bosonic space](bosonic_synthesis.html.md#bosonic-synthesis-explanation)
