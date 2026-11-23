<a id="module-qiskit_bosons.circuit.library.synthesis"></a>

# Circuit Library Synthesis (`B2B`)

This module provides the synthesis methods with which a bosonic gate is broken down into a
[`BosonicCircuit`](../stubs/qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) of smaller bosonic gates. Pass one to the gate’s `synthesis` argument to
choose how it gets decomposed. For an [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate that is its time evolution:

```python
Evolution(num_modes, operator, time=1.0, synthesis=BosonicLieTrotter())
```

or choose them for a whole circuit at once with the [`B2BSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2BSynthesis.html.md#qiskit_bosons.transpiler.passes.B2BSynthesis) transpiler pass. Any class
with a `synthesize(gate)` method can leverage a plugin (see [`B2BSynthesisPlugin`](../stubs/qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin)), and the ones
registered under the `qiskit_bosons.circuit.synthesis` entry point are available by name through
the [`B2BSynthesisPluginManager`](../stubs/qiskit_bosons.circuit.library.synthesis.B2BSynthesisPluginManager.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPluginManager).

Both sides of this rewrite live in bosonic space, so these methods can exploit the
[`groups`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.groups) of the evolved operator (see
[Group operator terms: use the operator structure](../guides/grouping.html.md#grouping-explanation)) – unlike the boson-to-qubit synthesis that follows, which only sees
the mapped operator, where the grouping is no longer available. That is what lets the *outer*
Trotter error between non-commuting groups be reduced by a higher order or more repetitions, which a
product formula handed to `MapperFnEvolutionSynthesis`
cannot do.

A [`ModeRotation`](../stubs/qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation) gate is broken down the same way, into [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates, one for
every two-mode rotation and single-mode phase of its mesh. Its `synthesis` is `None` by default,
which leaves the gate whole:

```python
ModeRotation(unitary, synthesis=BosonicClementsDecomposition())
```

#### IMPORTANT
This boson-to-boson step is **optional**. An [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate can be handed straight to
the boson-to-qubit stage regardless of how many terms its operator holds, and so can a
[`ModeRotation`](../stubs/qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation), given a plugin registered for it. Decomposing in bosonic space first is a
choice, taken to obtain cheaper factors or to expose structure that the later stage can exploit.

<a id="builtin-b2b-synthesis-plugins"></a>

## Builtin plugins

The table below lists the builtin plugins, organized by the [`BosonicGate`](../stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) that they act
upon. These are the entry points of the `qiskit_bosons.circuit.synthesis` group, named
`<gate>.<method>`:

#### An overview of the builtin `qiskit_bosons.circuit.synthesis` entry-points.

| Circuit instruction name                                                                                 | Plugin method name   | Reference                                                                                                                                          |
|----------------------------------------------------------------------------------------------------------|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution)       | `LieTrotter`         | [`BosonicLieTrotter`](../stubs/qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter.html.md#qiskit_bosons.circuit.library.synthesis.BosonicLieTrotter)                       |
| [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution)       | `SuzukiTrotter`      | [`BosonicSuzukiTrotter`](../stubs/qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter.html.md#qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter)                 |
| [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution)       | `QDrift`             | [`BosonicQDrift`](../stubs/qiskit_bosons.circuit.library.synthesis.BosonicQDrift.html.md#qiskit_bosons.circuit.library.synthesis.BosonicQDrift)                               |
| [`ModeRotation`](../stubs/qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation) | `Clements`           | [`BosonicClementsDecomposition`](../stubs/qiskit_bosons.circuit.library.synthesis.BosonicClementsDecomposition.html.md#qiskit_bosons.circuit.library.synthesis.BosonicClementsDecomposition) |
