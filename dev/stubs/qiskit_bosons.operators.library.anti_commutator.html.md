# anti_commutator

### anti_commutator(op_a, op_b)

Computes the anti-commutator of two operators.

#### NOTE
The *commutator* is the bracket that defines the bosonic algebra; the anti-commutator is well
defined but carries no comparable structural meaning. It is provided for completeness and for
operators built from products of ladder operators, where either bracket may be the useful one.

The anti-commutator is defined as:

$$
\{A, B\} = AB + BA
$$

Operators which support this method must implement the [`SupportsCommutators`](qiskit_bosons.protocols.SupportsCommutators.html.md#qiskit_bosons.protocols.SupportsCommutators) protocol.

#### NOTE
Both inputs must be of the same operator type. This will also determine the output type.

```pycon
>>> from qiskit_bosons.operators import BosonOperator, ann, cre
>>> from qiskit_bosons.operators.library import anti_commutator
>>> b = BosonOperator.from_dict({(ann(0),): 1.0})
>>> b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
>>> expected = BosonOperator.from_dict({(cre(0), ann(0)): 2.0, (): 1.0})
>>> anti_commutator(b, b_dag).equiv(expected)
True
```

* **Parameters:**
  * **op_a** (*T*) – the operator $A$ above.
  * **op_b** (*T*) – the operator $B$ above.
* **Returns:**
  The anti-commutator $\{A, B\}$.
* **Return type:**
  *T*
