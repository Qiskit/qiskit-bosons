<a id="release-notes"></a>

# Release Notes

<a id="release-notes-upcoming-release-stable-0-1"></a>

## Upcoming release (`stable/0.1`)

<a id="release-notes-upcoming-release-stable-0-1-new-features"></a>

### New Features

- Added `expand_projectors()`, a mapper function adapter next
  to `simplify()` and
  `group_wise()`. It expands the projectors of a mapped
  operator into Pauli terms and merges the ones that coincide, keeping the order in which the
  terms first appear.

  This matters for groups that combine terms which do not commute with the projectors they
  contain, for instance the hopping and density terms of a spin-1/2 chain written with hard-core
  bosons. [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) expands a projector only while
  synthesizing the gate, and does so term by term, so the pieces that cancel are rotated one by one
  between the terms they do not commute with. A group that is exactly solvable as a sum of
  commuting Paulis then no longer produces its exact evolution, and `simplify` cannot help
  because it merges before the expansion. Wrap the mapper inside `group_wise`:
  ```default
  MapperFnEvolutionSynthesis(group_wise(expand_projectors(mapper_fn)))
  ```

  The expansion is exponential in the number of qubits that a diagonal term acts on, so it is
  meant for the compact encodings or a small cutoff.

- [`Evolution`](stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) has a new `equiv` method, which compares
  two evolutions up to a tolerance on their times and their operators, and its instances are now
  hashable.

<a id="release-notes-upcoming-release-stable-0-1-upgrade-notes"></a>

### Upgrade Notes

- [`MergeEvolutions`](stubs/qiskit_bosons.transpiler.passes.MergeEvolutions.html.md#qiskit_bosons.transpiler.passes.MergeEvolutions) now merges two neighboring
  [`Evolution`](stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates whose operators are
  [`equiv()`](stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.equiv), that is equal in normal order, such as
  $b b^\dagger$ and $b^\dagger b + 1$. Before, it merged only operators that were the
  same products of ladder operators as written. A bosonic circuit lives in the untruncated algebra,
  where the merge is exact. The truncation that lowering to qubits introduces does not commute
  with normal ordering, so the lowered circuit of such a merge can differ from the one of the
  unmerged gates on the highest level of a mode, which is an artifact of the cutoff either way.
  To keep the lowered circuit the same with and without the pass, write the operators of
  evolutions on the same modes in one order, for example in normal order.

<a id="release-notes-upcoming-release-stable-0-1-bug-fixes"></a>

### Bug Fixes

- Two [`Evolution`](stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates for the same time on as many modes
  no longer compare equal when their operators differ. The gate used Qiskit’s equality of
  instructions, which compares the parameters and the definitions but not the operator, which is
  no parameter. It now compares its number of modes, its time and its operator, the latter with
  [`equiv()`](stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.equiv), that is in the algebra of the bosons,
  where the circuit lives: $b b^\dagger$ and $b^\dagger b + 1$ are the same operator,
  and the order of the terms does not matter. The truncation of a mode only enters when the
  circuit is lowered to qubits. The [`synthesis`](stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.synthesis) and
  [`atomic`](stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.atomic) flag are not compared: they decide how the
  gate is broken down, not what it is.

<a id="release-notes-0-1-0a1"></a>

## 0.1.0a1

<a id="release-notes-0-1-0a1-prelude"></a>

### Prelude

This is the first release of Qiskit Bosons, an extension of Qiskit with tools for bosonic
systems.  It provides a bosonic operator data structure, a framework and library for encoding
those operators into qubits (mappers), and a framework and library for bosonic circuits that a
transpilation pipeline synthesizes into qubit circuits.

Following Qiskit’s design philosophy, the core of this package is written in Rust and exposed
through first-party Python bindings.

This is an alpha release: the API may change in response to feedback.  Please [open an issue](https://github.com/Qiskit/qiskit-bosons/issues/new/choose) to share yours.
