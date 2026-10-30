# B2BSynthesisPlugin

### *class* B2BSynthesisPlugin(\*args, \*\*kwargs)

The protocol for plugins to the [`B2BSynthesis`](qiskit_bosons.transpiler.passes.B2BSynthesis.html.md#qiskit_bosons.transpiler.passes.B2BSynthesis) transpiler pass.

Bases: [`Protocol`](https://docs.python.org/3/library/typing.html#typing.Protocol)

A plugin rewrites a gate into a [`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) of *smaller* [`BosonicGate`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate)
instances. Both sides of that rewrite live in bosonic space, so neither an occupation cutoff
nor an encoding is involved: how the smaller gates are then mapped to qubits is left to the
boson-to-qubit stage (see [`B2QSynthesis`](qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis)). It is the counterpart of
[`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) one step earlier. Because the bosonic structure of the gate
survives the rewrite, the later stage can still exploit it, in particular the
[`groups`](qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.groups) of an evolved operator (see
[Group operator terms: use the operator structure](../guides/grouping.html.md#grouping-explanation)).

The gates a plugin emits should be marked [`atomic`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate.atomic), and an [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution)
among them should carry an operator with a single group, so that the boson-to-qubit stage maps
every one of them whole.

#### IMPORTANT
What the circuit has to be equivalent to is specified by the *gate* a plugin is registered
for, not by this protocol, and it may differ between gates: a circuit that implements the same
unitary is the common case (up to the error of the method for an [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution), and on
the untruncated Fock space for a [`ModeRotation`](qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation)), but a gate may settle for less.
Consult the documentation of the gate whose method you implement, and document any such
caveat there.

#### IMPORTANT
This boson-to-boson step is **optional**. A gate can be handed to the boson-to-qubit stage
whole, where a plugin of that stage registered for it can map it directly. Where an
approximation is involved at both levels, the accuracy of the result is governed by the weaker
of the two.

#### NOTE
Qiskit’s [`EvolutionSynthesis`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.synthesis.EvolutionSynthesis) is the closest analogue, but it belongs
to the qubit layer: it synthesizes a [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) into
a [`QuantumCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.QuantumCircuit), after the operator has already been mapped and the
occupation truncation has been applied. A bosonic synthesis, by contrast, carries no cutoff.

### Methods

#### synthesize(gate)

Synthesizes the provided gate.

* **Parameters:**
  **gate** ([*BosonicGate*](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate)) – the gate to synthesize, which is of the one type that this plugin is registered for
  (the first part of its entry-point name).
* **Returns:**
  A [`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) on `gate.num_modes` modes implementing the gate, in the
  sense that the documentation of the gate defines.
* **Return type:**
  [BosonicCircuit](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit)
