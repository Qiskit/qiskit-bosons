# OccupationInitializeModesSynthesis

### *class* OccupationInitializeModesSynthesis(spec)

A [`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) for transpiling [`InitializeModes`](qiskit_bosons.circuit.library.InitializeModes.html.md#qiskit_bosons.circuit.library.InitializeModes).

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

Prepares each mode’s Fock state by setting the qubits of the bit pattern that the chosen
occupation encoding assigns to that mode’s occupation level. This is a plain layer of
[`XGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.XGate) on the qubits the pattern sets – no entangling gate is
needed for a product Fock state.

#### WARNING
This plugin makes the following assumptions, none of which it verifies:

- **The qubits start in the all-zero state, and the gate is the first one acting on its
  modes.** The pattern is applied with `X` gates, so it *toggles* the qubits rather than
  setting them. A second [`InitializeModes`](qiskit_bosons.circuit.library.InitializeModes.html.md#qiskit_bosons.circuit.library.InitializeModes) gate on the same modes undoes the first
  one or, depending on the encoding, leaves the register in a state that no level of the
  truncation represents.
- **All-zero is the vacuum only under the compact encodings.** Under both unary encodings the
  pattern of level 0 sets a qubit, so a mode that no [`InitializeModes`](qiskit_bosons.circuit.library.InitializeModes.html.md#qiskit_bosons.circuit.library.InitializeModes) gate acts on
  starts outside the one-hot subspace. Place an [`InitializeModes`](qiskit_bosons.circuit.library.InitializeModes.html.md#qiskit_bosons.circuit.library.InitializeModes) gate on every mode
  of such a circuit, with [`InitializeModes.vacuum()`](qiskit_bosons.circuit.library.InitializeModes.html.md#qiskit_bosons.circuit.library.InitializeModes.vacuum) for the empty ones.
- **The** [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) **is the one used by the layout and the mapper.** The layout
  derives the qubit register from its spec and the mapper maps the operators with its own;
  this plugin places the pattern at the offsets of the spec it is given. If they differ, the
  qubits are set without an error, but not the qubits that the rest of the pipeline reads as
  that occupation. The preset builds all three from one spec.

#### NOTE
This plugin needs to know the **encoding** and the **truncation**, because “which qubits
represent occupation 2” is exactly what an encoding defines. They are given by a
[`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) passed to the constructor rather than being read off the gate: the gate
is encoding-agnostic by design, and deferring this choice to transpile time is what lets one
bosonic circuit be lowered under several encodings. Each mode is checked against, and placed
at the qubit offset of, its own entry in the spec.

Initializing this plugin can be done with the arguments listed below.

* **Parameters:**
  **spec** ([*LevelSpec*](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec)) – the [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) with which the rest of the pipeline truncates and encodes
  its modes. It must be the same spec as the layout and the mapper use.

### Attributes

#### spec

The [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) placing each mode’s levels onto qubits.

### Methods

#### run(in_node, out_dag, , b2q_layout)

Runs this transpilation plugin.

* **Parameters:**
  * **in_node** ([*DAGOpNode*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode)) – the input boson-based circuit instruction. When this plugin gets called, the
    `in_node.op` attribute must be of type [`InitializeModes`](qiskit_bosons.circuit.library.InitializeModes.html.md#qiskit_bosons.circuit.library.InitializeModes).
  * **out_dag** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)) – the output qubit-based circuit.
  * **b2q_layout** ([*dict*](https://docs.python.org/3/builtins/stdtypes.html#dict) *[*[*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *,* [*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *]*) – the global transpilation [`B2QLayout`](qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout)
    setting.
* **Return type:**
  None

#### SEE ALSO
The documentation of [`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) for more detailed explanations of the
arguments.

* **Raises:**
  * [**NotImplementedError**](https://docs.python.org/3/builtins/exceptions.html#NotImplementedError) – when `in_node` acts on bosonic modes that are spread across
        multiple [`BosonicRegister`](qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister) instances.
  * [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if the gate declares an occupation that the truncation of its mode cannot
        represent, or if the register has a different number of modes than [`spec`](#qiskit_bosons.transpiler.passes.synthesis.OccupationInitializeModesSynthesis.spec)
        covers.
* **Parameters:**
  * **in_node** ([*DAGOpNode*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode))
  * **out_dag** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit))
  * **b2q_layout** ([*dict*](https://docs.python.org/3/builtins/stdtypes.html#dict) *[*[*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *,* [*QuantumRegister*](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#qiskit.circuit.QuantumRegister) *]*)
* **Return type:**
  None
