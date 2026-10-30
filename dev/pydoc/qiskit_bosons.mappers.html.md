<a id="module-qiskit_bosons.mappers"></a>

# Representation Mappers

This module provides a framework for implementing custom boson-to-qubit mapper routines.

## Level encodings

A [`LevelEncoding`](../stubs/qiskit_bosons.mappers.LevelEncoding.html.md#qiskit_bosons.mappers.LevelEncoding) decides how the $d$ levels of a single system are laid out on qubits.
It says nothing about what the levels are, so it is shared by every mapper that works in a basis of
discrete levels: the occupation-number basis of a truncated bosonic mode, vibrational modals, or the
points of a position grid. The mappers in [`library`](qiskit_bosons.mappers.library.html.md#module-qiskit_bosons.mappers.library) are given one as
part of a [`LevelSpec`](../stubs/qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) (see below); for example
[`occupation_encoding()`](../stubs/qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding) maps a
[`BosonOperator`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator) in the occupation basis.

A [`LevelSpec`](../stubs/qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) combines a level encoding with the truncation of every mode of a register. It
holds one cutoff per mode (or one cutoff shared by all modes) and from them derives how many qubits
each mode occupies and where its qubits start, so everything that has to agree about that is handed
the same spec.

## Generator-based mappers

The [`library`](qiskit_bosons.mappers.library.html.md#module-qiskit_bosons.mappers.library) submodule ships the ready-made mappers. The generator
function below is the extension point for writing your own: supply the qubit operator that a single
[`BosonAction`](../stubs/qiskit_bosons.operators.boson_action.BosonAction.html.md#qiskit_bosons.operators.boson_action.BosonAction) maps to, and it lifts that into a mapper
over whole operators.
