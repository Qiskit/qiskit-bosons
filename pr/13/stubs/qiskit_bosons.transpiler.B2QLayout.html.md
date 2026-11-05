# B2QLayout

### B2QLayout

A mapping of bosonic mode registers to quantum registers.

Users must provide a data structure of this type to the circuit transpiler. In a trivial case, every
bosonic mode gets mapped to a single qubit. However, boson-to-qubit mappings generally use several
qubits per mode, which is why we do not associate single bosonic modes with single qubits.
Crucially, the size of the registers on either side of this mapping may differ.

alias of [`dict`](https://docs.python.org/3/builtins/stdtypes.html#dict)[[`QuantumRegister`](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister), [`QuantumRegister`](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister)]
