# generate_preset_occupation_pass_manager

### generate_preset_occupation_pass_manager(spec, \*\*kwargs)

Generates a preset transpiler pipeline for the occupation basis, from a [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec).

The truncation and the encoding are given by the `spec` of this function rather than by the
circuit: a bosonic circuit is encoding-agnostic, and choosing them here is what lets the same
circuit be lowered at several truncation levels without being rebuilt.

The layout, the state preparation and the operator mapping of the pipeline are all derived from
the one `spec`, so they cannot disagree about which qubits a mode occupies.

#### NOTE
A spec built from a single cutoff leaves the number of modes open, so the pipeline serves
circuits of any width. A spec with one cutoff per mode fixes it, and the pipeline then only
supports circuits with exactly that many modes.

* **Parameters:**
  * **spec** ([*LevelSpec*](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec)) – the [`LevelSpec`](qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) with the cutoff of every mode and the encoding that lays the
    occupation levels out on qubits.
  * **kwargs** – any additional keyword arguments are forwarded to
    [`generate_preset_pass_manager()`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.transpiler.generate_preset_pass_manager) whose output is used
    for the `qubit` stage. Without any, that stage uses Qiskit’s defaults.
* **Returns:**
  The preset staged boson-to-qubit transpiler pipeline.
* **Return type:**
  [*MultiStagePassManager*](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.passmanager.MultiStagePassManager)
