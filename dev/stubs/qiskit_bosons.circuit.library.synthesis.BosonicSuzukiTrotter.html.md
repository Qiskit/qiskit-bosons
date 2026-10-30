# BosonicSuzukiTrotter

### *class* BosonicSuzukiTrotter(order=2, reps=1, , merge_adjacent=True)

The higher-order Suzuki-Trotter product formula, applied in bosonic space.

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

A product formula composes the factors of an operator symmetrically, so that the leading error
terms cancel. For $\texttt{order} = 2$ and factors $H_j$, the approximation is the
palindrome

$$
e^{-i t \sum_j H_j} \approx
\left( \prod_{j<n} e^{-i \frac{t}{2} H_j} \right) e^{-i t H_n}
\left( \prod_{j<n} e^{-i \frac{t}{2} H_j} \right)^{\mathrm{R}},
$$

where $\mathrm{R}$ denotes the reversed order. Higher even orders are built from this by
Suzuki’s recursion, each level composing five sweeps of the level below.

The evolved operator is split by its [`groups`](qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.groups) when
it has them (see [Group operator terms: use the operator structure](../guides/grouping.html.md#grouping-explanation)), and term-by-term otherwise. Every order splits it
the same way: what a higher order changes is the ordering of the factors and their time scales, not
the partition.

Splitting along the groups is something this stage can do and Qiskit’s qubit-side
[`SuzukiTrotter`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.synthesis.SuzukiTrotter) cannot: once the operator has been mapped and handed to
a [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate), the grouping is gone, and the formula can
only subdivide the terms within the one operator it is given.

#### NOTE
A higher order buys accuracy with depth: an order-$k$ formula emits roughly
$5^{(k-2)/2}$ times as many factors as an order-2 one. Increasing [`reps`](#qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter.reps) instead
divides the evolution time into more, shorter steps, which reduces the error at a different
rate for the same kind of cost. Which is the better trade depends on the operator.

#### NOTE
Applying this in bosonic space at all is optional, and it composes with the product formula
chosen for the boson-to-qubit stage. Where both approximate, the accuracy of the result is
governed by the weaker of the two. Bosonic factors are, moreover, exact *before* the
occupation truncation and only approximately so afterwards (the truncated algebra does not
close), so a reference value should be computed at the same `max_occupation`.

#### CAUTION
Each factor must be Hermitian for its exponential to be unitary, which does not follow from
their sum being Hermitian: splitting a Hermitian operator can produce non-Hermitian groups (for
example, separating the hopping term $b^\dagger_0 b_1$ from its adjoint
$b^\dagger_1 b_0$). It is the caller’s responsibility to group accordingly; this is not
verified, but `groups_are_hermitian()` checks it.

An operator with no groups is split term by term, which has the same effect, since an
individual term is generally not Hermitian (the number operator $b^\dagger_i b_i$ and
the like being exceptions). The resulting factors are valid bosonic operators and the gates
are built without complaint, but the boson-to-qubit stage then rejects them:
[`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) requires real coefficients and raises
`ValueError: Operator contains complex coefficients, which are not supported`.

Note that the symmetrization partially cancels the error of a non-Hermitian factor, so an
incorrectly grouped operator can look markedly better at an even order than at first order
while still being wrong.

Initializing an instance of this synthesis method can be done with the arguments below.

* **Parameters:**
  * **order** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – the order of the product formula. Must be `1` (which reduces to
    [`BosonicLieTrotter`](qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter.html.md#qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter)) or a positive even integer, since the Suzuki formulas
    are symmetric.
  * **reps** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – the number of times to repeat the formula, each repetition evolving for
    `time / reps`.
  * **merge_adjacent** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – whether to fuse neighboring gates of the same factor into one, summing
    their evolution times. See [`merge_adjacent`](#qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter.merge_adjacent).
* **Raises:**
  [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if `order` is neither `1` nor a positive even integer, or if `reps` is
      not positive.

### Attributes

#### order

The order of the product formula.

#### reps

The number of times the product formula is repeated.

#### merge_adjacent

Whether neighboring gates of the same factor are fused into one.

The palindromic structure makes the last factor of one repetition (and of each sub-sweep of
a higher order) the first of the next. Their exponentials combine exactly, as
$e^{-i a H} e^{-i b H} = e^{-i (a + b) H}$, so fusing them shortens the circuit without
changing the product formula. Disable this to obtain the textbook gate sequence.

This only sees the seams created by a single synthesis. Merging across gate boundaries is
the job of [`MergeEvolutions`](qiskit_bosons.transpiler.passes.MergeEvolutions.html.md#qiskit_bosons.transpiler.passes.MergeEvolutions).

### Methods

#### synthesize(gate)

Synthesizes the gate into the ordered factors of the product formula.

See the class documentation for the formula this implements.

* **Parameters:**
  **gate** ([*Evolution*](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution)) – the gate to synthesize.
* **Returns:**
  A [`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) holding one narrowed [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate per factor of
  the formula. A factor may appear more than once, at different evolution times.
* **Return type:**
  [BosonicCircuit](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit)
