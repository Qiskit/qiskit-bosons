# OccupationB2QLayout

### *class* OccupationB2QLayout(spec)

Maps bosonic modes to the qubits a [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) assigns them.

Bases: [`GenericPass`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.GenericPass)[[`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit), [`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)]

Each bosonic mode is truncated to the levels its cutoff in [`spec`](#qiskit_bosons.transpiler.passes.OccupationB2QLayout.spec) keeps and spans
[`LevelEncoding.num_qubits()`](qiskit_bosons.mappers.LevelEncoding.html.md#qiskit_bosons.mappers.LevelEncoding.num_qubits) of that many levels consecutive qubits, so the register of
$M$ modes is mapped to a qubit register of [`LevelSpec.total_width`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec.total_width) qubits and mode
$m$ starts at qubit [`LevelSpec.offsets`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec.offsets) `[m]`.

A boson-to-qubit encoding generally needs more than one qubit per mode ($d$ of them for the
unary encodings), so the qubit register is wider than the mode register. Unlike
[`CustomB2QLayout`](qiskit_bosons.transpiler.passes.CustomB2QLayout.html.md#qiskit_bosons.transpiler.passes.CustomB2QLayout) this derives the registers from the circuit being transpiled rather
than requiring them up front, which is what lets one pass manager serve circuits of different
widths, as long as [`spec`](#qiskit_bosons.transpiler.passes.OccupationB2QLayout.spec) leaves the number of modes open.

#### NOTE
The truncation and the encoding are given by the [`spec`](#qiskit_bosons.transpiler.passes.OccupationB2QLayout.spec) of this **pass**, not by the
circuit. That is the point: the bosonic circuit upstream of the transpiler is independent of
both, so the same
circuit can be lowered at several truncation levels – or under several encodings – by
swapping this pass alone.

#### WARNING
A spec whose modes have **different** cutoffs fixes the qubit offset of each mode by its
position in the register, so reordering the modes afterwards would silently attach the wrong
cutoff to a mode. This pass therefore refuses a circuit to which [`RelabelModes`](qiskit_bosons.transpiler.passes.RelabelModes.html.md#qiskit_bosons.transpiler.passes.RelabelModes) has
already been applied (it leaves a `"permutation"` entry in the circuit’s metadata) and
raises a [`ValueError`](https://docs.python.org/3/builtins/exceptions.html#ValueError). A spec with one shared cutoff is unaffected. Note that
[`CustomB2QLayout`](qiskit_bosons.transpiler.passes.CustomB2QLayout.html.md#qiskit_bosons.transpiler.passes.CustomB2QLayout) does not perform this check: the user is responsible for the
consistency of a layout they provide themselves.

Initializing this transpiler pass can be done with the arguments listed below.

* **Parameters:**
  **spec** ([*LevelSpec*](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec)) – the [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) with the cutoff of every mode and the encoding that lays
  the occupation levels out on qubits. A spec with one cutoff per mode only supports
  circuits with exactly that many modes.

### Attributes

#### spec

The [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) this layout is derived from.

### Methods

#### run(dag)

Runs this analysis pass.

* **Parameters:**
  **dag** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)) – the DAG circuit representation of a [`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit).
* **Returns:**
  The unchanged input. But the `property_set` of this transpilation pipeline will
  be updated.
* **Raises:**
  [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if a register has a different number of modes than [`spec`](#qiskit_bosons.transpiler.passes.OccupationB2QLayout.spec) covers, or
      if [`spec`](#qiskit_bosons.transpiler.passes.OccupationB2QLayout.spec) gives the modes different cutoffs and [`RelabelModes`](qiskit_bosons.transpiler.passes.RelabelModes.html.md#qiskit_bosons.transpiler.passes.RelabelModes) has
      been applied to the circuit.
* **Return type:**
  [*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)

### Inherited Methods

#### execute(passmanager_ir, state, callback=None)

Execute optimization task for input Qiskit IR.

* **Parameters:**
  * **passmanager_ir** (*IR*) – Qiskit IR to optimize.
  * **state** ([*PassManagerState*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.PassManagerState)) – State associated with workflow execution by the pass manager itself.
  * **callback** ([*Callable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable) *[* *[*[*Task*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.Task) *,* *IR_OUT* *,* [*PropertySet*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.PropertySet) *,* [*float*](https://docs.python.org/3/builtins/functions.html#float) *,* [*int*](https://docs.python.org/3/builtins/functions.html#int) *]* *,* *None* *]*  *|* *None*) – A callback function which is called per execution of optimization task.
* **Returns:**
  Optimized Qiskit IR and state of the workflow.
* **Return type:**
  [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[*IR_OUT*, [*PassManagerState*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.PassManagerState)]

#### name()

Name of the pass.

* **Return type:**
  [str](https://docs.python.org/3/builtins/stdtypes.html#str)

#### update_status(state, run_state)

Update workflow status.

* **Parameters:**
  * **state** ([*PassManagerState*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.PassManagerState)) – Pass manager state to update.
  * **run_state** (*RunState*) – Completion status of current task.
* **Returns:**
  Updated pass manager state.
* **Return type:**
  [*PassManagerState*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.PassManagerState)
