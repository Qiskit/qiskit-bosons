<a id="module-qiskit_bosons.mappers.library"></a>

# Mapper Library

This module provides the boson-to-qubit mappers.

[`occupation_encoding()`](../stubs/qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding) maps a [`BosonOperator`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator) in the occupation
(Fock) basis, truncating each mode and laying its levels out on qubits as a
[`LevelSpec`](../stubs/qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) describes.

A bosonic encoding is **not** an algebra isomorphism: a bosonic mode’s Fock space is
infinite-dimensional and must be truncated first, which breaks the canonical commutation relation on
the highest retained level. See [`occupation_encoding()`](../stubs/qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding) for the exact statement and its
consequences.
