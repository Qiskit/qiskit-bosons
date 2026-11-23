# SupportsCommutators

### *class* SupportsCommutators(\*args, \*\*kwargs)

A Protocol indicating support for efficient commutator generation.

Bases: [`Protocol`](https://docs.python.org/3/library/typing.html#typing.Protocol)

Implementation of this protocol requires the three methods below. See [`commutator()`](qiskit_bosons.operators.library.commutator.html.md#qiskit_bosons.operators.library.commutator),
[`anti_commutator()`](qiskit_bosons.operators.library.anti_commutator.html.md#qiskit_bosons.operators.library.anti_commutator), and [`double_commutator()`](qiskit_bosons.operators.library.double_commutator.html.md#qiskit_bosons.operators.library.double_commutator) for the type-agnostic helper functions
dispatching to these methods.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> cre_0 = BosonOperator.from_dict({((True, 0),): 1.0})
>>> ann_0 = BosonOperator.from_dict({((False, 0),): 1.0})
>>> comm = BosonOperator._commutator_(ann_0, cre_0)
>>> print(format(comm.normal_ordered().simplify()))
  1.000000e0 +0.000000e0j * ()
```

### Protocol Methods

#### *static* \_anti_commutator_(op_a, op_b)

Computes the anti-commutator of two operator instances.

See [`anti_commutator()`](qiskit_bosons.operators.library.anti_commutator.html.md#qiskit_bosons.operators.library.anti_commutator) for more details.

* **Parameters:**
  * **op_a** (*T*)
  * **op_b** (*T*)
* **Return type:**
  *T*

#### *static* \_commutator_(op_a, op_b)

Computes the commutator of two operator instances.

See [`commutator()`](qiskit_bosons.operators.library.commutator.html.md#qiskit_bosons.operators.library.commutator) for more details.

* **Parameters:**
  * **op_a** (*T*)
  * **op_b** (*T*)
* **Return type:**
  *T*

#### *static* \_double_commutator_(op_a, op_b, op_c, sign)

Computes the double-commutator of three operator instances.

See [`double_commutator()`](qiskit_bosons.operators.library.double_commutator.html.md#qiskit_bosons.operators.library.double_commutator) for more details.

* **Parameters:**
  * **op_a** (*T*)
  * **op_b** (*T*)
  * **op_c** (*T*)
  * **sign** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool))
* **Return type:**
  *T*
