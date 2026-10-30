# commutator

### commutator(op_a, op_b)

Computes the commutator of two operators.

The commutator is defined as:

$$
[A, B] = AB - BA
$$

Operators which support this method must implement the [`SupportsCommutators`](qiskit_bosons.protocols.SupportsCommutators.html.md#qiskit_bosons.protocols.SupportsCommutators) protocol.

#### NOTE
Both inputs must be of the same operator type. This will also determine the output type.

```pycon
>>> from qiskit_bosons.operators import BosonOperator, ann, cre
>>> from qiskit_bosons.operators.library import commutator
>>> b = BosonOperator.from_dict({(ann(0),): 1.0})
>>> b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
>>> commutator(b, b_dag).equiv(BosonOperator.one())
True
>>> n = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
>>> commutator(n, b_dag).equiv(b_dag)
True
>>> b1_dag = BosonOperator.from_dict({(cre(1),): 1.0})
>>> commutator(b_dag, b1_dag).equiv(BosonOperator.zero())
True
```

* **Parameters:**
  * **op_a** (*T*) – the operator $A$ above.
  * **op_b** (*T*) – the operator $B$ above.
* **Returns:**
  The commutator $[A, B]$.
* **Return type:**
  *T*
