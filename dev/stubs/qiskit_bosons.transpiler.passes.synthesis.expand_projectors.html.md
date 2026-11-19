# expand_projectors

### expand_projectors(mapper_fn)

Wraps a mapper function so that it expands the projectors of its result into Paulis.

The occupation encodings write the number operator as a projector,
$\lvert 1 \rangle\langle 1 \rvert$, which a [`SparseObservable`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.quantum_info.SparseObservable)
stores as a single factor. [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) expands it to
$(I - Z)/2$ only when it synthesizes the gate, term by term and without merging what
cancels. The `Z` pieces of different terms then end up as separate rotations, interleaved with
the terms they do not commute with, and the result is not the evolution that the group stands for
even where the group is Hermitian and its Pauli form is a sum of commuting pieces.

This wrapper does that expansion up front and merges the Pauli terms that coincide, so a product
formula sees the operator that was meant. Each Pauli string takes the place of its first
occurrence, so the order that [`group_wise()`](qiskit_bosons.transpiler.passes.synthesis.group_wise.html.md#qiskit_bosons.transpiler.passes.synthesis.group_wise) or the mapper chose survives. Terms whose summed
coefficient is below a small fraction of the largest are dropped.

Use it to wrap any `MapperFunction`, **inside** [`group_wise()`](qiskit_bosons.transpiler.passes.synthesis.group_wise.html.md#qiskit_bosons.transpiler.passes.synthesis.group_wise):

```default
MapperFnEvolutionSynthesis(group_wise(expand_projectors(mapper_fn)))
```

An operator without projectors is returned with its terms merged and nothing else changed.

#### WARNING
The expansion is exponential in the number of qubits that a diagonal term acts on: a
projector product over $w$ qubits becomes $2^w$ Pauli terms. This is harmless for
the compact encodings at a small cutoff, but the unary encodings put every level of a mode on
its own qubit and can blow up. Apply it to the groups that need it, and do not use it on a
large operator.

#### CAUTION
Wrapping [`group_wise()`](qiskit_bosons.transpiler.passes.synthesis.group_wise.html.md#qiskit_bosons.transpiler.passes.synthesis.group_wise) with this function, as in `expand_projectors(group_wise(...))`,
merges terms *across* groups and so discards the grouping. Put it inside.

* **Parameters:**
  **mapper_fn** ([*Callable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable) *[* *[*[*OperatorTrait*](qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait) *,* [*int*](https://docs.python.org/3/builtins/functions.html#int) *]* *,* [*SparseObservable*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.quantum_info.SparseObservable) *]*) – the mapper function to wrap.
* **Returns:**
  A mapper function that returns the Pauli expansion of whatever `mapper_fn` returns. It
  does not alter the imaginary parts: an operator that is not Hermitian is still rejected by
  [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate), exactly as without this wrapper.
* **Return type:**
  [*Callable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable)[[[*OperatorTrait*](qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait), [int](https://docs.python.org/3/builtins/functions.html#int)], [*SparseObservable*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.quantum_info.SparseObservable)]
