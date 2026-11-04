# B2QSynthesisPlugin

### *class* B2QSynthesisPlugin(\*args, \*\*kwargs)

The protocol for plugins to the [`B2QSynthesis`](qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis) transpiler pass.

Bases: [`Protocol`](https://docs.python.org/3/library/typing.html#typing.Protocol)

### Methods

#### run(in_node, out_dag, , b2q_layout)

Translates the provided boson-based circuit instruction to a qubit-based one.

* **Parameters:**
  * **in_node** ([*DAGOpNode*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode)) – a boson-based circuit instruction stored in a
    [`DAGOpNode`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode). Specifically, this guarantees that
    [`op`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode#op) is of type [`BosonicGate`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate).
  * **out_dag** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)) – the qubit-based [`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit) into which this plugin
    must insert the translated circuit instruction.
  * **b2q_layout** ([*dict*](https://docs.python.org/3/builtins/stdtypes.html#dict) *[*[*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *,* [*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *]*) – the [`B2QLayout`](qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout) setting that is global to
    the transpilation process. It is the plugin’s responsibility to respect this mapping
    of [`BosonicRegister`](qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister) to
    [`QuantumRegister`](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister).
* **Return type:**
  None
