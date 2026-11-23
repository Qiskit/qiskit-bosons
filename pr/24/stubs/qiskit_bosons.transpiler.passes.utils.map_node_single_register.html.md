# map_node_single_register

### map_node_single_register(in_node, b2q_layout)

Map the single bosonic mode register acted upon by the circuit operation.

Given a [`DAGOpNode`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode) and a
[`B2QLayout`](qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout) this function will:

1. ensure that all acted-upon bosonic modes belong to the same register
2. map the node’s indices to their global indices with respect to said bosonic mode register
3. map this bosonic mode register to the corresponding qubit register

#### CAUTION
The returned global bosonic mode register indices are knowingly dependent on the order of
the bosonic modes during the construction of their original
[`BosonicRegister`](qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister). This is unlikely to cause any issues but
makes the usage of [`RelabelModes`](qiskit_bosons.transpiler.passes.RelabelModes.html.md#qiskit_bosons.transpiler.passes.RelabelModes) possibly more generally.

* **Parameters:**
  * **in_node** ([*DAGOpNode*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode)) – the node whose indices to parse.
  * **b2q_layout** ([*dict*](https://docs.python.org/3/builtins/stdtypes.html#dict) *[*[*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *,* [*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *]*) – the mapping of bosonic to qubit registers.
* **Returns:**
  The list of the global bosonic mode indices and mapped qubit register instance.
* **Raises:**
  [**NotImplementedError**](https://docs.python.org/3/builtins/exceptions.html#NotImplementedError) – when the provided node acts on bosonic modes that do not all belong
      to the same [`BosonicRegister`](qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister).
* **Return type:**
  [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[list](https://docs.python.org/3/builtins/stdtypes.html#list)[[int](https://docs.python.org/3/builtins/functions.html#int)], [*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister)]
