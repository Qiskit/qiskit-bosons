<a id="module-qiskit_bosons.operators.terms"></a>

# Operator Terms

This module collects routines that operate on the individual terms of an operator, partitioning
them based on their structure. The submodules below group related functionality.

## Grouping

Terms may be tagged with group indices via
[`groups`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.groups), which an
[`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate then honors by decomposing into one
sub-evolution per group. Grouping terms *before* the boson-to-qubit mapping is what lets a grouping
chosen from the problem’s structure survive into the circuit.

No automatic grouping strategy ships yet; assign the indices directly. What this module does
provide is the *analysis* of an existing grouping: [`group_coeff_means()`](../stubs/qiskit_bosons.operators.terms.group_coeff_means.html.md#qiskit_bosons.operators.terms.group_coeff_means) reduces it to one
sampling weight per group, [`groups_are_hermitian()`](../stubs/qiskit_bosons.operators.terms.groups_are_hermitian.html.md#qiskit_bosons.operators.terms.groups_are_hermitian) reports whether each group is separately
Hermitian (which a randomized product formula, and Qiskit’s `PauliEvolutionGate`, require), and
[`groups_have_uniform_coeffs()`](../stubs/qiskit_bosons.operators.terms.groups_have_uniform_coeffs.html.md#qiskit_bosons.operators.terms.groups_have_uniform_coeffs) whether its coefficients are numerically equal.

## Filtering

This module provides convenience functions for removing terms from an operator that do not
contribute meaningfully to a downstream computation.

### Members

## Ordering

This module provides functions to reorder the terms of an operator. While the
resulting operators are mathematically identical, the term order can have
significant implications on algorithm behavior at runtime, due to the term
iteration order.

### Members

### Group-ordered operators

[`group_order()`](../stubs/qiskit_bosons.operators.terms.group_order.html.md#qiskit_bosons.operators.terms.group_order) makes each group one contiguous run of terms. That layout is what lets
[`split_out_groups()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.split_out_groups) binary-search the requested groups
instead of scanning every term, so ordering once up front pays off before repeatedly transpiling with
[`BosonicQDrift`](../stubs/qiskit_bosons.circuit.library.synthesis.BosonicQDrift.html.md#qiskit_bosons.circuit.library.synthesis.BosonicQDrift), which looks up its sampled groups on
every call.
