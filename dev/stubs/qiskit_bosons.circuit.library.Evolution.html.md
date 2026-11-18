# Evolution

### *class* Evolution(num_modes, operator, time=1.0, , synthesis=None, atomic=False)

Implements the time evolution of a bosonic operator.

Bases: [`BosonicGate`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate)

Given a bosonic `operator` $H$ and an evolution `time` $t$, this gate
implements the unitary

$$
U = e^{-i t H}.
$$

For $U$ to be unitary, $H$ must be Hermitian. This is the caller’s responsibility
and is not verified; [`is_hermitian()`](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.is_hermitian) can check it.

#### NOTE
Breaking this gate down in bosonic space (into a circuit of smaller [`Evolution`](#qiskit_bosons.circuit.library.Evolution)
gates) is **optional**, and it is what the [`synthesis`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.synthesis) attribute
governs. The gate can equally be handed straight to the boson-to-qubit stage, which maps and
synthesizes it whole, however many terms its `operator` holds. See
[`synthesis`](../pydoc/qiskit_bosons.circuit.library.synthesis.html.md#module-qiskit_bosons.circuit.library.synthesis).

When it is decomposed (by [`decompose()`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit.decompose), for instance)
[`synthesis`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.synthesis) decides how. It defaults to [`BosonicLieTrotter`](qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter.html.md#qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter),
which splits the evolution group-by-group when the `operator` has
[`groups`](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.groups) assigned (see
[Group operator terms: use the operator structure](../guides/grouping.html.md#grouping-explanation)), and term-by-term otherwise (exact only when the individual
factors mutually commute). Another method can be chosen, for a whole circuit,
with the [`B2BSynthesis`](qiskit_bosons.transpiler.passes.B2BSynthesis.html.md#qiskit_bosons.transpiler.passes.B2BSynthesis) transpiler pass; [`BosonicQDrift`](qiskit_bosons.circuit.library.synthesis.BosonicQDrift.html.md#qiskit_bosons.circuit.library.synthesis.BosonicQDrift), for example,
replaces the evolution with a randomized sample of its terms. The gates it emits are [`atomic`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.atomic), so
[`synthesis`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.synthesis) no longer applies to them.

#### NOTE
This gate carries **no occupation cutoff**. It is defined purely in terms of the bosonic
operator, whose Fock space is infinite-dimensional. The truncation is introduced later, by the
boson-to-qubit encoding chosen during transpilation – which is precisely what lets one
ansatz be evaluated at several truncation levels without rebuilding it. See
[`occupation_encoding()`](qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding).

Initializing an instance of this gate can be done with the arguments listed below.

* **Parameters:**
  * **num_modes** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – the number of bosonic modes on which this gate acts.
  * **operator** ([*OperatorTrait*](qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait)) – the Hermitian operator $H$ under which to time evolve the acted-upon
    bosonic modes.
  * **time** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – the evolution time $t$ entering the exponent of $e^{-i t H}$. A
    negative value evolves backwards in time.
  * **synthesis** ([*B2BSynthesisPlugin*](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin) *|* *None*) – the boson-to-boson synthesis method with which to decompose this gate, when
    it gets decomposed at all. If `None` (the default), a
    [`BosonicLieTrotter`](qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter.html.md#qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter) instance is used.
  * **atomic** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – whether this gate is a fully decomposed factor which must not be broken down any
    further in bosonic space. See [`BosonicGate.atomic`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.atomic).

### Attributes

#### operator

The operator under which to time evolve the acted-upon bosonic modes.

### Methods

#### equiv(other, atol=1e-10)

Whether `other` is an evolution of the same operator for the same time, up to `atol`.

The operators are compared with [`BosonOperator.equiv()`](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.equiv), and numeric times up to
`atol` as well. Like the circuit, this is a comparison in the algebra of the bosons:
$b b^\dagger$ and $b^\dagger b + 1$ are the same operator, although a truncation
(which only enters when the circuit is lowered to qubits) maps them differently on the top
level of a mode.

* **Parameters:**
  * **other** ([*object*](https://docs.python.org/3/builtins/functions.html#object)) – the operation to compare with.
  * **atol** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – the absolute tolerance of the comparison.
* **Returns:**
  Whether the two evolutions are the same.
* **Return type:**
  [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### inverse(annotated=False)

Returns the inverse of this gate, $e^{+i t H}$.

* **Parameters:**
  **annotated** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – ignored. The inverse of this gate is another [`Evolution`](#qiskit_bosons.circuit.library.Evolution), so it never
  needs to be an [`AnnotatedOperation`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.AnnotatedOperation).
* **Returns:**
  An [`Evolution`](#qiskit_bosons.circuit.library.Evolution) gate evolving the same operator for the negated time.
* **Return type:**
  [*Evolution*](#qiskit_bosons.circuit.library.Evolution)
