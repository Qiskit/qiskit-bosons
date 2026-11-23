# map_boson_action_generators

### map_boson_action_generators(operator, map_action, identity, compose=None)

Map a [`BosonOperator`](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator) to another operator type.

This is a generic function to aid in implementing new mappers for [`BosonOperator`](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator)
instances. At its core, it simply iterates over the terms of the operator, mapping each
encountered [`BosonAction`](qiskit_bosons.operators.boson_action.BosonAction.html.md#qiskit_bosons.operators.boson_action.BosonAction) with the user-provided `map_action` function. In
combination with the user-provided `identity` generator, this allows mapping to arbitrary
output types.

#### NOTE
The output type `T` must support multiplication by a scalar via `__mul__`.
If `compose=None` it must also support composition of two instances via `__and__`.
[`SparseObservable`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.quantum_info.SparseObservable) does not, so pass its
[`compose()`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.quantum_info.SparseObservable#compose) explicitly as `compose`; a type with
an `__and__` (such as [`SparsePauliOp`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.quantum_info.SparsePauliOp)) can rely on the default.

#### NOTE
This function is an extension point for your own mappers. It is a minimal illustration and
not a replacement for [`occupation_encoding()`](qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding), which is
implemented in Rust, is parallelized, and bounds the memory it uses while assembling the
result. Nothing in the library calls this function.

The example maps a mode with two levels (`max_occupation = 1`) to one qubit, where the
annihilation operator is $\lvert 0 \rangle\langle 1 \rvert = (X + iY) / 2$ and the
creation operator is its adjoint. Modes with more levels need the encoding of
[`occupation_encoding()`](qiskit_bosons.mappers.library.occupation_encoding.html.md#qiskit_bosons.mappers.library.occupation_encoding) instead.

```pycon
>>> from qiskit_bosons.mappers import map_boson_action_generators
>>> from qiskit_bosons.operators import BosonAction, BosonOperator, ann, cre
>>> from qiskit.quantum_info import SparsePauliOp
>>>
>>> num_qubits = 2
>>> def map_action(action: BosonAction) -> SparsePauliOp:
...     is_creation, mode = action
...     return SparsePauliOp.from_sparse_list(
...         [("X", [mode], 0.5), ("Y", [mode], -0.5j if is_creation else 0.5j)],
...         num_qubits=num_qubits,
...     )
>>>
>>> def identity() -> SparsePauliOp:
...     return SparsePauliOp.from_sparse_list([("", [], 1)], num_qubits)
>>>
>>> op = BosonOperator.from_dict({(cre(0), ann(1)): 2.0})
>>> qop = map_boson_action_generators(op, map_action, identity).simplify()
>>> print([(label, complex(coeff)) for label, coeff in sorted(qop.label_iter())])
[('XX', (0.5+0j)), ('XY', -0.5j), ('YX', 0.5j), ('YY', (0.5+0j))]
```

* **Parameters:**
  * **operator** ([*BosonOperator*](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator)) – the operator to be mapped.
  * **map_action** ([*Callable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable) *[* *[*[*tuple*](https://docs.python.org/3/builtins/stdtypes.html#tuple) *[*[*bool*](https://docs.python.org/3/builtins/functions.html#bool) *,* [*int*](https://docs.python.org/3/builtins/functions.html#int) *]* *]* *,* *T* *]*) – the function to map a single [`BosonAction`](qiskit_bosons.operators.boson_action.BosonAction.html.md#qiskit_bosons.operators.boson_action.BosonAction) to the desired output type.
  * **identity** ([*Callable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable) *[* *[* *]* *,* *T* *]*) – the function to generate the multiplicative identity instance of the output type.
  * **compose** ([*Callable*](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable) *[* *[**T* *,* *T* *]* *,* *T* *]*  *|* *None*) – an optional function to implement the composition logic of two output type
    instances. If this is not provided, it will default to using [`operator.and_()`](https://docs.python.org/3/library/operator.html#operator.and_).
* **Returns:**
  The mapped operator.
* **Return type:**
  *T*
