# MergeEvolutions

### *class* MergeEvolutions(atol=1e-12)

A transpilation pass merging consecutive [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates of the same operator.

Bases: [`GenericPass`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.GenericPass)[[`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit), [`DAGCircuit`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)]

Two evolutions under the same operator $H$ satisfy

$$
e^{-i t_2 H} e^{-i t_1 H} = e^{-i (t_1 + t_2) H},
$$

because $H$ commutes with itself. This pass finds neighboring [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates
of the same operator and replaces them by one gate whose evolution time is the sum of theirs.
The rewrite is exact, so it changes neither the simulated state nor the Trotter error of the
circuit; it only removes a gate, and with it the synthesized operations that gate would have
produced.

Such neighbors arise wherever two sequences of factors meet on the same factor. Examples are the
seams between the repetitions of [`BosonicCircuit.repeat()`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit.repeat), or between consecutive samples
of [`BosonicQDrift`](qiskit_bosons.circuit.library.synthesis.BosonicQDrift.html.md#qiskit_bosons.circuit.library.synthesis.BosonicQDrift). [`BosonicSuzukiTrotter`](qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter.html.md#qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter) already fuses the seams within
its own output (see [`BosonicSuzukiTrotter.merge_adjacent`](qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter.html.md#qiskit_bosons.circuit.library.synthesis.BosonicSuzukiTrotter.merge_adjacent)); this pass is the general
mechanism which also acts across gate boundaries.

Two gates merge only when all of the following hold:

- they act on the identical wires, in the identical order, with nothing else in between;
- their operators are of the same type and [`BosonOperator.equiv()`](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.equiv) up to [`atol`](#qiskit_bosons.transpiler.passes.MergeEvolutions.atol), that
  is equal in the algebra of the bosons, see the note below;
- they are either both `Evolution.atomic` or both not, and in the latter case they are
  decomposed by the same kind of `Evolution.synthesis` with identical settings.

The merged gate keeps the operator, and thus the
[`groups`](qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.groups), of the first gate of the run.

#### NOTE
The merge is exact in the algebra of the bosons, where the circuit lives. A truncation, which
only enters when the circuit is lowered to qubits, does not commute with normal ordering:
$b b^\dagger$ and $b^\dagger b + 1$ are the same operator, but they map to
matrices that differ on the highest level of the truncated mode. So when two evolutions whose
operators only agree in normal order are merged, the lowered circuit can differ from the one
of the unmerged gates on that level, which is an artifact of the cutoff either way. Write the
operators in one order (for example normal order) if the lowered circuit has to be the same
with and without this pass.

#### CAUTION
This is an early development prototype. Beware of changes to its interface without warning
during the pre-release development of this package.

#### SEE ALSO
[`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution).

Initializing this transpiler pass can be done with the arguments listed below.

* **Parameters:**
  **atol** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – the absolute tolerance up to which the operators of two gates are considered
  equal.

### Attributes

#### atol

The absolute tolerance up to which two operators are considered equal.

### Methods

#### run(dag)

Runs this transpilation pass.

Collects every maximal run of consecutive [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates on the same wires,
splits it into stretches of mergeable gates and replaces each stretch of length two or more
with a single [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) carrying the summed time. Everything else is left
untouched. The input DAG is modified in place.

A `Barrier` interrupts a run, so the gates on either side of
one are not merged. Place a barrier between two gates to keep them apart.

* **Parameters:**
  **dag** ([*DAGCircuit*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGCircuit)) – the input circuit with boson-based instructions. Only
  [`DAGOpNode`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode) with [`BosonicGate`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) instances as their
  [`op`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.dagcircuit.DAGOpNode#op) are supported, plus
  `Barrier`, which is carried through untouched.
* **Returns:**
  The output circuit which is still acting on a bosonic register.
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
