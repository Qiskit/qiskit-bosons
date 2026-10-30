<a id="module-qiskit_bosons.transpiler"></a>

# Transpiler

## Overview

Transpilation is the process of rewriting a given input circuit to conform with desired criteria
such as topological and operational constraints of the hardware used to execute the final circuit.
We are not going to explain this in more detail here, and instead refer to Qiskit’s documentation of
the [`qiskit.transpiler`](https://quantum.cloud.ibm.com/docs/api/qiskit/transpiler#module-qiskit.transpiler) module.

The focus here lies on explaining how we achieve the transpilation of a [`BosonicCircuit`](../stubs/qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit)
to a [`QuantumCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.QuantumCircuit).

## Stages

Conceptually, we split the transpilation process into several stages:

| Stage                                                       | Description                                |
|-------------------------------------------------------------|--------------------------------------------|
| [Input](#qiskit-bosons-transpiler-stage-input)        | convert to DAG data structure              |
| [Optimization](#qiskit-bosons-transpiler-stage-optimization) | bosonic-level optimization                 |
| [Layout](#qiskit-bosons-transpiler-stage-layout)       | boson-to-qubit layouting                   |
| [Synthesis](#qiskit-bosons-transpiler-stage-synthesis)    | boson-to-qubit synthesis                   |
| [Qubit](#qiskit-bosons-transpiler-stage-qubit)        | continued transpilation on the qubit-level |
| [Output](#qiskit-bosons-transpiler-stage-output)       | convert from DAG data structure            |

<a id="qiskit-bosons-transpiler-stage-input"></a>

### Input

The various transpiler passes are implemented to work with [`BosonicDAGCircuit`](../stubs/qiskit_bosons.circuit.BosonicDAGCircuit.html.md#qiskit_bosons.circuit.BosonicDAGCircuit) as the
underlying data structure to store the circuit operations. Such directed acyclic graphs (DAGs)
provide an efficient data model for the traversal and manipulation of circuits.

However, an end-user is more likely to work with a [`BosonicCircuit`](../stubs/qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) as it provides a more
intuitive interface and data model. As such, this input stage simply runs the
[`BosonicCircuitToDAG`](../stubs/qiskit_bosons.transpiler.BosonicCircuitToDAG.html.md#qiskit_bosons.transpiler.BosonicCircuitToDAG) transpiler pass.

<a id="qiskit-bosons-transpiler-stage-optimization"></a>

### Optimization

This stage of the transpilation pipeline can implement circuit optimizations while preserving the
type of circuit to be an instance of [`BosonicDAGCircuit`](../stubs/qiskit_bosons.circuit.BosonicDAGCircuit.html.md#qiskit_bosons.circuit.BosonicDAGCircuit). As such, no qubit information is
required (or necessarily available) at this point in the transpilation pipeline.

<a id="qiskit-bosons-transpiler-stage-layout"></a>

### Layout

One global configuration setting for the transpilation process is the boson-to-qubit “layout”.
This must be provided by the user and it must match the provided boson-to-qubit mapping (in the
sense that, if a chosen mapping encodes a fixed number of bosons with a different number of
qubits, the configured layout must account for that).

In the general case, bosonic modes are not always encoded with an occupation-basis into the qubit
register. Consequently, we cannot associate a single bosonic mode with a single qubit.
Therefore, the user-provided boson-to-qubit layout ([`B2QLayout`](../stubs/qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout)) associates
[`BosonicRegister`](../stubs/qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister) instances with
[`QuantumRegister`](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) ones.

The way to configure this global setting, is by placing a [`B2QLayout`](../stubs/qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout) instance in the
`b2q_layout` field of the `property_set`. For more details
refer to [Layouting Passes](qiskit_bosons.transpiler.passes.html.md#qiskit-bosons-transpiler-passes-layout).

<a id="qiskit-bosons-transpiler-stage-synthesis"></a>

### Synthesis

At its core, the transpilation is handled by the [`B2QSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis) transpiler pass. It is
conceptually similar to Qiskit’s [`HighLevelSynthesis`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.transpiler.passes.HighLevelSynthesis) pass, which
uses various plugins for transpiling high-level circuit instructions. For more details, refer to
the documentation of [`B2QSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis) directly.

How a given [`BosonicGate`](../stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) can be synthesized in terms of qubit-based operations will
depend on the particular gate type as well as the user-chosen boson-to-qubit mapping. For more
details, refer to [Transpiler Pass Synthesis (B2Q)](qiskit_bosons.transpiler.passes.synthesis.html.md#qiskit-bosons-transpiler-passes-synthesis-plugins).

<a id="qiskit-bosons-transpiler-stage-qubit"></a>

### Qubit

At this point in the transpilation process, we have reached a
[`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit) instance and can continue to use Qiskit’s transpilation
pipeline as one would usually.

#### HINT
Additional transpiler passes for optimizations on the qubit-level that take into account the
knowledge of a circuit originating from a [`BosonicCircuit`](../stubs/qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) may be added in the future!

<a id="qiskit-bosons-transpiler-stage-output"></a>

### Output

This stage implements effectively the reverse of [Input](#qiskit-bosons-transpiler-stage-input), by
calling the [`QuantumDAGToCircuit`](../stubs/qiskit_bosons.transpiler.QuantumDAGToCircuit.html.md#qiskit_bosons.transpiler.QuantumDAGToCircuit) transpiler pass.

## Pass Managers

All of the stages above are the default stages by how the [`presets`](qiskit_bosons.transpiler.presets.html.md#module-qiskit_bosons.transpiler.presets)
orchestrate the various `MultiStagePassManager`.

The individual stages can be either a single pass (when they change from one internal representation
(IR) to another, such as the [Input](#qiskit-bosons-transpiler-stage-input),
[Synthesis](#qiskit-bosons-transpiler-stage-synthesis), and
[Output](#qiskit-bosons-transpiler-stage-output) stages above) or a stage can be a
[`BasePassManager`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.BasePassManager) whose internal passes can be modified.

For the stages operating on bosonic circuits ([Optimization](#qiskit-bosons-transpiler-stage-optimization)
and [Layout](#qiskit-bosons-transpiler-stage-layout)), the type of passmanager to use is
[`BosonicPassManager`](../stubs/qiskit_bosons.transpiler.BosonicPassManager.html.md#qiskit_bosons.transpiler.BosonicPassManager).

### Conversion Passes

Some very basic conversion passes are provided directly by this module:

### Presets

For user convenience, the [`qiskit_bosons.transpiler.presets`](qiskit_bosons.transpiler.presets.html.md#module-qiskit_bosons.transpiler.presets) module provides a number of
functions for quickly building pre-defined transpilation pipelines.
