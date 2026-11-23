# BosonicGate

### *class* BosonicGate(name, num_modes, , params=None, , label=None, synthesis=None, atomic=False)

The base class for all bosonic gates.

Bases: [`Gate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.Gate)

To ensure consistency only subclasses of this gate class can be added to instances of
[`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit). As such, this class (mostly) serves as a type (for the time being).

### Breaking a gate down in bosonic space

Every bosonic gate has a [`synthesis`](#qiskit_bosons.circuit.BosonicGate.synthesis) and an [`atomic`](#qiskit_bosons.circuit.BosonicGate.atomic) flag, which the boson-to-boson
synthesis framework relies on (see [`synthesis`](../pydoc/qiskit_bosons.circuit.library.synthesis.html.md#module-qiskit_bosons.circuit.library.synthesis)). The
[`synthesis`](#qiskit_bosons.circuit.BosonicGate.synthesis) is the method that rewrites the gate into a circuit of smaller bosonic gates,
and it is `None` for a gate that does not break down. A gate that is [`atomic`](#qiskit_bosons.circuit.BosonicGate.atomic) is never
broken down, however it is asked to. Both are chosen when the gate is built and are read-only
afterwards.

#### CAUTION
Since this is a subclass of [`Gate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.Gate) the documentation of its methods
may refer to qubits. Those references should be interpreted as referring to bosonic modes
in the context of instances of this subclass.

It may also happen that some of the inherited methods may not always make sense because of
this re-interpretation. You have been warned.

#### NOTE
A bosonic gate carries **no occupation cutoff**: it acts on bosonic modes, whose Fock space is
infinite-dimensional. A gate becomes a finite-dimensional operator only once a boson-to-qubit
encoding fixes a truncation, which happens during transpilation. Simulating a bosonic circuit
therefore means transpiling it to qubits first and then simulating the resulting qubit
circuit with any qubit-level simulator; there is deliberately no mode-space state-vector
protocol on this class.

* **Parameters:**
  * **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – The name of the gate.
  * **num_qubits** – The number of qubits the gate acts on.
  * **params** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list) *|* *None*) – A list of parameters.
  * **label** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*) – An optional label for the gate.
  * **num_modes** ([*int*](https://docs.python.org/3/builtins/functions.html#int))
  * **synthesis** ([*B2BSynthesisPlugin*](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin) *|* *None*)
  * **atomic** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool))

### Attributes

#### atomic

Whether this gate is a factor that must not be broken down any further.

A synthesis marks the factors it emits as atomic, because they are the result of the split it
chose to perform: breaking them down again would discard that choice, and would split a
Hermitian group into individually non-Hermitian terms whose exponentials are not unitary. A
gate that is atomic for the stronger reason that there is nothing left to split is marked the
same way.

An atomic gate has no definition at all, so [`decompose()`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit.decompose) leaves it in
place rather than expanding it, and the [`B2BSynthesis`](qiskit_bosons.transpiler.passes.B2BSynthesis.html.md#qiskit_bosons.transpiler.passes.B2BSynthesis) pass leaves it alone. It still
maps and synthesizes normally at the boson-to-qubit level, which is where a factor is turned
into actual operations.

This is read-only, for the same reason as [`synthesis`](#qiskit_bosons.circuit.BosonicGate.synthesis).

#### num_modes

The number of bosonic modes that this gate acts upon.

#### synthesis

The boson-to-boson synthesis method with which this gate is broken down, if any.

This is the customizable part of breaking a gate down: any object implementing the
[`B2BSynthesisPlugin`](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin) protocol can be given, and the [`B2BSynthesis`](qiskit_bosons.transpiler.passes.B2BSynthesis.html.md#qiskit_bosons.transpiler.passes.B2BSynthesis) transpiler
pass selects one per gate type for a whole circuit. `None` marks a gate that does not break
down in bosonic space.

What the circuit that the method produces is equivalent to is specified by the gate, and
documented there. Where the gate’s decomposition is that circuit,
[`decompose()`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit.decompose) yields it; a gate may deviate from this and document so.

This is read-only: the circuit that [`decompose()`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit.decompose) yields is built once
and then cached, so a synthesis method assigned after the fact could not take effect. Pass it
to the constructor instead, or use [`with_synthesis()`](#qiskit_bosons.circuit.BosonicGate.with_synthesis).

### Inherited Attributes

#### base_class

Get the base class of this instruction.  This is guaranteed to be in the inheritance tree
of `self`.

The “base class” of an instruction is the lowest class in its inheritance tree that the
object should be considered entirely compatible with for \_all_ circuit applications.  This
typically means that the subclass is defined purely to offer some sort of programmer
convenience over the base class, and the base class is the “true” class for a behavioral
perspective.  In particular, you should *not* override [`base_class`](#qiskit_bosons.circuit.BosonicGate.base_class) if you are
defining a custom version of an instruction that will be implemented differently by
hardware, such as an alternative measurement strategy, or a version of a parametrized gate
with a particular set of parameters for the purposes of distinguishing it in a
`Target` from the full parametrized gate.

This is often exactly equivalent to `type(obj)`, except in the case of singleton instances
of standard-library instructions.  These singleton instances are special subclasses of their
base class, and this property will return that base.  For example:

```default
>>> isinstance(XGate(), XGate)
True
>>> type(XGate()) is XGate
False
>>> XGate().base_class is XGate
True
```

In general, you should not rely on the precise class of an instruction; within a given
circuit, it is expected that `Instruction.name` should be a more suitable
discriminator in most situations.

#### decompositions

Get the decompositions of the instruction from the SessionEquivalenceLibrary.

#### definition

Return definition in terms of other basic gates.

#### label

Return instruction label

#### mutable

Is this instance is a mutable unique instance or not.

If this attribute is `False` the gate instance is a shared singleton
and is not mutable.

#### name

Return the name.

#### num_clbits

Return the number of clbits.

#### num_qubits

Return the number of qubits.

#### params

The parameters of this `Instruction`.  Ideally these will be gate angles.

### Methods

#### with_synthesis(synthesis)

Returns a copy of this gate that is broken down with another synthesis method.

The gate itself is left as it is, because a gate instance can be shared between circuits.

* **Parameters:**
  **synthesis** ([*B2BSynthesisPlugin*](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin) *|* *None*) – the boson-to-boson synthesis method of the copy.
* **Returns:**
  A shallow copy of this gate whose [`synthesis`](#qiskit_bosons.circuit.BosonicGate.synthesis) is `synthesis`.
* **Return type:**
  Self

### Inherited Methods

#### add_decomposition(decomposition)

Add a decomposition of the instruction to the SessionEquivalenceLibrary.

#### broadcast_arguments(qargs, cargs)

Validation and handling of the arguments and its relationship.

For example, `cx([q[0],q[1]], q[2])` means `cx(q[0], q[2]); cx(q[1], q[2])`. This
method yields the arguments in the right grouping. In the given example:

```default
in: [[q[0],q[1]], q[2]],[]
outs: [q[0], q[2]], []
      [q[1], q[2]], []
```

The general broadcasting rules are:

> * If len(qargs) == 1:
>   ```default
>   [q[0], q[1]] -> [q[0]],[q[1]]
>   ```
> * If len(qargs) == 2:
>   ```default
>   [[q[0], q[1]], [r[0], r[1]]] -> [q[0], r[0]], [q[1], r[1]]
>   [[q[0]], [r[0], r[1]]]       -> [q[0], r[0]], [q[0], r[1]]
>   [[q[0], q[1]], [r[0]]]       -> [q[0], r[0]], [q[1], r[0]]
>   ```
> * If len(qargs) >= 3:
>   ```default
>   [q[0], q[1]], [r[0], r[1]],  ...] -> [q[0], r[0], ...], [q[1], r[1], ...]
>   ```
* **Parameters:**
  * **qargs** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – List of quantum bit arguments.
  * **cargs** ([*list*](https://docs.python.org/3/builtins/stdtypes.html#list)) – List of classical bit arguments.
* **Returns:**
  A tuple with single arguments.
* **Raises:**
  **CircuitError** – If the input is not valid. For example, the number of
      arguments does not match the gate expectation.
* **Return type:**
  [*Iterable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[list](https://docs.python.org/3/builtins/stdtypes.html#list), [list](https://docs.python.org/3/builtins/stdtypes.html#list)]]

#### control(num_ctrl_qubits=1, label=None, ctrl_state=None, annotated=None)

Return the controlled version of itself.

The controlled gate is implemented as `ControlledGate` when `annotated`
is `False`, and as `AnnotatedOperation` when `annotated` is `True`.

#### Deprecated
Deprecated since version 2.3: `qiskit.circuit.gate.Gate.control()`’s argument `annotated` is deprecated as of Qiskit 2.3. It will be removed in Qiskit 3.0. The method Gate.control() no longer accepts annotated=None. The new default is annotated=True, which represents the controlled gate as an AnnotatedOperation (unless a dedicated controlled-gate class already exists). You can explicitly set annotated=False to preserve the previous behavior. However, using annotated=True is recommended, as it defers construction of the controlled circuit to transpiler, and furthermore enables additional controlled-gate optimizations (typically leading to higher-quality circuits).

* **Parameters:**
  * **num_ctrl_qubits** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – Number of controls to add. Defaults to `1`.
  * **label** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*) – Optional gate label. Defaults to `None`.
    Ignored if the controlled gate is implemented as an annotated operation.
  * **ctrl_state** ([*int*](https://docs.python.org/3/builtins/functions.html#int) *|* [*str*](https://docs.python.org/3/builtins/stdtypes.html#str) *|* *None*) – The control state of the gate, specified either as an integer or a bitstring
    (e.g. `"110"`). If `None`, defaults to the all-ones state `2**num_ctrl_qubits - 1`.
  * **annotated** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool) *|* *None*) – Indicates whether the controlled gate should be implemented as a controlled gate
    or as an annotated operation. If `None`, treated as `False`.
* **Returns:**
  A controlled version of this gate.
* **Raises:**
  **QiskitError** – invalid `num_ctrl_qubits` or `ctrl_state`.

#### copy(name=None)

Copy of the instruction.

* **Parameters:**
  **name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – name to be given to the copied circuit, if `None` then the name stays the same.
* **Returns:**
  a copy of the current instruction, with the name updated if it
  was provided
* **Return type:**
  [qiskit.circuit.Instruction](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.Instruction)

#### inverse(annotated=False)

Invert this instruction.

If annotated is False, the inverse instruction is implemented as
a fresh instruction with the recursively inverted definition.

If annotated is True, the inverse instruction is implemented as
`AnnotatedOperation`, and corresponds to the given instruction
annotated with the “inverse modifier”.

Special instructions inheriting from Instruction can
implement their own inverse (e.g. T and Tdg, Barrier, etc.)
In particular, they can choose how to handle the argument `annotated`
which may include ignoring it and always returning a concrete gate class
if the inverse is defined as a standard gate.

* **Parameters:**
  **annotated** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – if set to True the output inverse gate will be returned
  as `AnnotatedOperation`.
* **Returns:**
  The inverse operation.
* **Raises:**
  **CircuitError** – if the instruction is not composite
      and an inverse has not been implemented for it.

#### is_parameterized()

Return whether the `Instruction` contains [compile-time parameters](https://quantum.cloud.ibm.com/docs/api/qiskit/circuit#circuit-compile-time-parameters).

#### power(exponent, annotated=False)

Raise this gate to the power of `exponent`.

Implemented either as a unitary gate (ref. `UnitaryGate`)
or as an annotated operation (ref. `AnnotatedOperation`). In the case of several standard
gates, such as `RXGate`, when the power of a gate can be expressed in terms of another
standard gate that is returned directly.

* **Parameters:**
  * **exponent** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – the power to raise the gate to
  * **annotated** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – indicates whether the power gate can be implemented
    as an annotated operation. In the case of several standard
    gates, such as `RXGate`, this argument is ignored when
    the power of a gate can be expressed in terms of another
    standard gate.
* **Returns:**
  An operation implementing `gate^exponent`
* **Raises:**
  **CircuitError** – If gate is not unitary

#### repeat(n)

Creates an instruction with `self` repeated $n$ times.

* **Parameters:**
  **n** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – Number of times to repeat the instruction
* **Returns:**
  Containing the definition.
* **Return type:**
  [qiskit.circuit.Instruction](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.Instruction)
* **Raises:**
  **CircuitError** – If n < 1.

#### reverse_ops()

For a composite instruction, reverse the order of sub-instructions.

This is done by recursively reversing all sub-instructions.
It does not invert any gate.

* **Returns:**
  a new instruction with
  : sub-instructions reversed.
* **Return type:**
  [qiskit.circuit.Instruction](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.Instruction)

#### soft_compare(other)

Soft comparison between gates. Their names, number of qubits, and classical
bit numbers must match. The number of parameters must match. Each parameter
is compared. If one is a ParameterExpression then it is not taken into
account.

* **Parameters:**
  **other** (*instruction*) – other instruction.
* **Returns:**
  are self and other equal up to parameter expressions.
* **Return type:**
  [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### to_matrix()

Return a Numpy.array for the gate unitary matrix.

* **Returns:**
  if the Gate subclass has a matrix definition.
* **Return type:**
  np.ndarray
* **Raises:**
  **CircuitError** – If a Gate subclass does not implement this method an
      exception will be raised when this base class method is called.

#### to_mutable()

Return a mutable copy of this gate.

This method will return a new mutable copy of this gate instance.
If a singleton instance is being used this will be a new unique
instance that can be mutated. If the instance is already mutable it
will be a deepcopy of that instance.

#### validate_parameter(parameter)

Gate parameters should be int, float, or ParameterExpression
