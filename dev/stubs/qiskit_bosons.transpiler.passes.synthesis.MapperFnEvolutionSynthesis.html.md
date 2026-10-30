# MapperFnEvolutionSynthesis

### *class* MapperFnEvolutionSynthesis(mapper_fn, product_formula=None)

A [`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) for transpiling [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) under a custom mapping.

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

This plugin maps the bosonic Hamiltonian $H$ of the incoming [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate to
a qubit operator using [`mapper_fn`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.mapper_fn) and emits a
[`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate). It thereby preserves the
$e^{-i t H}$ convention of the [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate, with the same evolution time
$t$.

How that [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) is subsequently decomposed into
basis gates is governed by the product formula passed as [`product_formula`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.product_formula). Leaving it at
its default (`None`) defers to the [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate)’s own
default synthesis (a first-order [`LieTrotter`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.synthesis.LieTrotter) decomposition with a
single repetition). Supplying an explicit
[`EvolutionSynthesis`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.synthesis.EvolutionSynthesis) (for example a higher-order
[`SuzukiTrotter`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.synthesis.SuzukiTrotter) or one with several repetitions) selects a different
Trotter-Suzuki product formula, trading circuit depth for a smaller Trotter error.

#### NOTE
The operator returned by [`mapper_fn`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.mapper_fn) is passed on as-is, in particular **without** being
simplified. A product formula synthesizes the Pauli terms in the order it receives them, so the
term order that [`mapper_fn`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.mapper_fn) produces is part of its output and is preserved here. A
[`mapper_fn`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.mapper_fn) that maps an operator group by group, for instance, emits the terms of each
group together, which lets the terms of one group be scheduled in parallel where their supports
are disjoint. See [`group_wise()`](qiskit_bosons.transpiler.passes.synthesis.group_wise.html.md#qiskit_bosons.transpiler.passes.synthesis.group_wise), which wraps any mapper to do exactly that.

The occupation encodings write a ladder operator as `(X - iY)/2`, so the Pauli terms of a
Hermitian group carry imaginary coefficients that cancel only once the terms are merged.
[`occupation_encoding()`](qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding) merges its result, using a
tolerance relative to the largest coefficient so that the floating-point noise of the
cancellation is removed too; a mapper of your own has to return an operator with real
coefficients, since [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) rejects complex ones.

#### CAUTION
A consequence of preserving that order is that the synthesized circuit is only as reproducible
as [`mapper_fn`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.mapper_fn) is. The operators of the Rust core do **not** preserve the order in which
their terms were added, so a mapper that walks an operator’s terms can emit them in a different
order from one run to the next. The circuits that result are all equally valid (they
approximate the same evolution of the same operator) but they need not be identical, and
metrics such as depth or gate count can vary between them.

Wrap the mapper in [`simplify()`](qiskit_bosons.transpiler.passes.synthesis.simplify.html.md#qiskit_bosons.transpiler.passes.synthesis.simplify) to pin a canonical order where that matters:

```default
MapperFnEvolutionSynthesis(simplify(mapper_fn))
```

#### IMPORTANT
[`product_formula`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.product_formula) only subdivides the terms **within** the operator of a single
[`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate. It cannot interleave *separate* gates, so when a grouped operator
has already been decomposed into one [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) per group (see
[Group operator terms: use the operator structure](../guides/grouping.html.md#grouping-explanation)), raising its repetition count does **not** reduce the outer
Trotter error between the groups – which do not generally commute. The result stays a single
outer Trotter step, and its error does not vanish as the repetition count grows.

To converge the grouped evolution, repeat the whole thing: append `steps`
[`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates at `time / steps` each. This is the flip side of grouping being
decided upstream – the partition survives into the circuit, and so does responsibility for
stepping it.

Initializing this transpiler pass plugin can be done with the arguments listed below.

* **Parameters:**
  * **mapper_fn** (*MapperFunction*) – the boson-to-qubit operator mapping function.
  * **product_formula** (*EvolutionSynthesis* *|* *None*) – the product formula with which to synthesize the emitted
    [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate). If `None` (the default), the
    gate’s own default synthesis is used (a first-order
    [`LieTrotter`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.synthesis.LieTrotter) decomposition with a single repetition).

### Attributes

#### mapper_fn *: MapperFunction*

The boson-to-qubit operator mapping function.

The two input arguments should be the following:

1. the operator to be mapped.
2. the number of qubits that the resulting operator should be defined on.

#### NOTE
It is the user’s responsibility to ensure that this function is in-sync with the global
transpilation [`B2QLayout`](qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout) setting.

#### product_formula *: EvolutionSynthesis | [None](https://docs.python.org/3/builtins/constants.html#None)*

The product formula used to synthesize the emitted
[`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate), or `None` to defer to that gate’s own
default synthesis.

### Methods

#### run(in_node, out_dag, , b2q_layout)

Runs this transpilation plugin.

The bosonic Hamiltonian of the incoming [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate is mapped to a qubit
operator via [`mapper_fn`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.mapper_fn) and appended to `out_dag` as a
[`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) implementing $e^{-i t H}$ with the
original evolution time $t$ and the [`product_formula`](#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.product_formula) synthesis.

* **Parameters:**
  * **in_node** ([*DAGOpNode*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode)) – the input boson-based circuit instruction. When this plugin gets called, the
    `in_node.op` attribute must be of type [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution).
  * **out_dag** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)) – the output qubit-based circuit.
  * **b2q_layout** ([*dict*](https://docs.python.org/3/builtins/stdtypes.html#dict) *[*[*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *,* [*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *]*) – the global transpilation [`B2QLayout`](qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout)
    setting.
* **Return type:**
  None

#### SEE ALSO
The documentation of [`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) for more detailed explanations of the
arguments.

* **Raises:**
  [**NotImplementedError**](https://docs.python.org/3/builtins/exceptions.html#NotImplementedError) – when `in_node` acts on bosonic modes that are spread across
      multiple [`BosonicRegister`](qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister) instances.
* **Parameters:**
  * **in_node** ([*DAGOpNode*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode))
  * **out_dag** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit))
  * **b2q_layout** ([*dict*](https://docs.python.org/3/builtins/stdtypes.html#dict) *[*[*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *,* [*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *]*)
* **Return type:**
  None
