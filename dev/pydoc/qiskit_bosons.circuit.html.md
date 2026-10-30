<a id="module-qiskit_bosons.circuit"></a>

# Circuit Representation

With this module you can build [`QuantumCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.QuantumCircuit) objects that are
interpreted as acting on bosonic modes rather than qubits.

The following objects are simple type aliases to make this re-interpretation more clear:

For actually expressing your circuits, you should use the [`BosonicCircuit`](../stubs/qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) class along with
the [`library`](qiskit_bosons.circuit.library.html.md#module-qiskit_bosons.circuit.library) of [`BosonicGate`](../stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) implementations.
