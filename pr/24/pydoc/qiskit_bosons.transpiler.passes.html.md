<a id="module-qiskit_bosons.transpiler.passes"></a>

# Transpiler Passes

This module provides various transpiler passes for the stages explained in
[`qiskit_bosons.transpiler`](qiskit_bosons.transpiler.html.md#module-qiskit_bosons.transpiler).

<a id="qiskit-bosons-transpiler-passes-optimization"></a>

## Optimization Passes

These passes provide different kinds of optimization of [`BosonicDAGCircuit`](../stubs/qiskit_bosons.circuit.BosonicDAGCircuit.html.md#qiskit_bosons.circuit.BosonicDAGCircuit) instances.

<a id="qiskit-bosons-transpiler-passes-layout"></a>

## Layouting Passes

These passes are designed to help configuring the global
[`B2QLayout`](../stubs/qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout) setting for the transpilation process. This setting
needs to be placed in the `b2q_layout` field of the
[`property_set`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.PassManagerState#property_set), from where it will be read during the
[Synthesis Passes](#qiskit-bosons-transpiler-passes-synthesis).

<a id="qiskit-bosons-transpiler-passes-synthesis"></a>

## Synthesis Passes

The main logic for mapping a [`BosonicDAGCircuit`](../stubs/qiskit_bosons.circuit.BosonicDAGCircuit.html.md#qiskit_bosons.circuit.BosonicDAGCircuit) to a
[`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit) is implemented by a single synthesis pass, namely the
[`B2QSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis). It is conceptually similar to Qiskit’s
[`HighLevelSynthesis`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.transpiler.passes.HighLevelSynthesis) because it simply iterates over the
[`BosonicDAGCircuit`](../stubs/qiskit_bosons.circuit.BosonicDAGCircuit.html.md#qiskit_bosons.circuit.BosonicDAGCircuit) instructions and delegates the mapping to qubit-based instructions to
[Transpiler Pass Synthesis (B2Q)](qiskit_bosons.transpiler.passes.synthesis.html.md#qiskit-bosons-transpiler-passes-synthesis-plugins) for each of the encountered types of
[`BosonicGate`](../stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate).
