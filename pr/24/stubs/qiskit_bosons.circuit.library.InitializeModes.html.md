# InitializeModes

### *class* InitializeModes(occupation)

Prepares a bosonic Fock state.

Bases: [`BosonicGate`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate)

This gate declares the occupation of the modes it is placed on, i.e. the Fock state

$$
\lvert n_0, n_1, \ldots, n_{M-1} \rangle
= \prod_m \frac{(b^\dagger_m)^{n_m}}{\sqrt{n_m!}} \lvert \text{vac} \rangle \, ,
$$

where $n_m$ is `occupation[m]`. Synthesis emits the gates that set each mode’s qubits to
the pattern encoding its occupation level.

#### NOTE
The occupation entries are **non-negative integers**, not booleans: a single mode may hold
arbitrarily many bosons.

#### NOTE
This gate does not itself fix a truncation. The occupation it declares must of course be
representable by whatever cutoff the boson-to-qubit encoding uses (an occupation of 3
requires `max_occupation >= 3`) which the synthesis plugin checks at transpile time, since
that is where the cutoff is known. [`max_occupation`](#qiskit_bosons.circuit.library.InitializeModes.max_occupation) reports the largest level this gate
declares.

#### CAUTION
This is an early development prototype. Beware of changes to its interface without warning
during the pre-release development of this package.

Initializing an instance of this gate can be done with the arguments listed below.

* **Parameters:**
  **occupation** (*Sequence* *[*[*int*](https://docs.python.org/3/builtins/functions.html#int) *]*) – the boson number of each mode in the [`BosonicRegister`](qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister) being
  initialized by this gate.
* **Raises:**
  [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if any occupation entry is negative.

### Attributes

#### max_occupation

The largest occupation level this gate declares.

A boson-to-qubit encoding must retain at least this many levels to represent the declared
Fock state; see [`occupation_encoding()`](qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding).

#### occupation *: [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[int](https://docs.python.org/3/builtins/functions.html#int), ...]*

The boson number of each mode in the [`BosonicRegister`](qiskit_bosons.circuit.BosonicRegister.html.md#qiskit_bosons.circuit.BosonicRegister) initialized by this gate.

### Methods

#### *classmethod* vacuum(num_modes)

Builds the gate for the vacuum state $\lvert 0, \ldots, 0 \rangle$.

* **Parameters:**
  **num_modes** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – the number of bosonic modes to initialize.
* **Returns:**
  An [`InitializeModes`](#qiskit_bosons.circuit.library.InitializeModes) gate declaring zero bosons in every mode.
* **Return type:**
  [*InitializeModes*](#qiskit_bosons.circuit.library.InitializeModes)
