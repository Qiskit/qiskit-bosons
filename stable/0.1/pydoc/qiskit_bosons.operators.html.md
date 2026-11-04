<a id="module-qiskit_bosons.operators"></a>

# Operator Representations

This module provides the data structures for representing bosonic operators.

## Boson Operator

This operator represents bosonic operators in terms of the second-quantization creation and
annihilation operators, which obey the canonical commutation relations

$$
\left[b_i, b^\dagger_j\right] = \delta_{ij}, \qquad
\left[b_i, b_j\right] = \left[b^\dagger_i, b^\dagger_j\right] = 0 \, .
$$

#### NOTE
A bosonic mode has an **infinite-dimensional** Fock space. [`BosonOperator`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator) is exact and
carries no occupation cutoff; a truncation is introduced only when mapping to qubits (see
[`qiskit_bosons.mappers`](qiskit_bosons.mappers.html.md#module-qiskit_bosons.mappers)). Two consequences are visible directly in this API: a mode index may
repeat within a single term, and [`max_occupation()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.max_occupation) and
[`mode_occupations()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.mode_occupations) report the cutoff such a term demands.

## Operator Library

A library of the model Hamiltonians that recur in bosonic simulation is provided in
[`qiskit_bosons.operators.library`](qiskit_bosons.operators.library.html.md#module-qiskit_bosons.operators.library).
