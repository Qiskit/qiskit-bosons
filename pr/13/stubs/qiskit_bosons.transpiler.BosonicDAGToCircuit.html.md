# BosonicDAGToCircuit

### *class* BosonicDAGToCircuit

Converts a [`BosonicDAGCircuit`](qiskit_bosons.circuit.BosonicDAGCircuit.html.md#qiskit_bosons.circuit.BosonicDAGCircuit) to a [`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit).

Bases: [`GenericPass`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.GenericPass)[[`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit), [`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit)]

### Methods

#### run(passmanager_ir)

Runs this conversion pass.

* **Parameters:**
  **passmanager_ir** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)) – the input circuit to convert.
* **Returns:**
  The converted circuit.
* **Return type:**
  [*BosonicCircuit*](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit)

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
