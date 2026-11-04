# double_commutator

### double_commutator(op_a, op_b, op_c, sign)

Computes the double-commutator of three operators.

The double-commutator is defined as follows (see also Chapter 13.6, *Equation of motion
methods*, page 479 of <sup>[1](#id3)</sup>; the `sign=False` case is the symmetric double commutator
$[A, B, C]$ appearing in Equation (2) of <sup>[2](#id4)</sup>):

If `sign` is `False`, it returns

$$
[[A, B], C]/2 + [A, [B, C]]/2 = (2ABC + 2CBA - BAC - CAB - ACB - BCA)/2.

$$

If `sign` is `True`, it returns

$$
\{[A, B], C\}/2 + \{A, [B, C]\}/2 = (2ABC - 2CBA - BAC + CAB - ACB + BCA)/2.

$$

Operators which support this method must implement the [`SupportsCommutators`](qiskit_bosons.protocols.SupportsCommutators.html.md#qiskit_bosons.protocols.SupportsCommutators) protocol.

#### NOTE
All three inputs must be of the same operator type. This will also determine the output type.

```pycon
>>> from qiskit_bosons.operators import BosonOperator, ann, cre
>>> from qiskit_bosons.operators.library import double_commutator
>>> op1 = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
>>> op2 = BosonOperator.from_dict({(cre(1), ann(0)): 1.0})
>>> op3 = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
>>> result = double_commutator(op1, op2, op3, False)
>>> other = double_commutator(op1, op2, op3, True)
>>> result.equiv(other)
False
>>> scaled = double_commutator(3.0 * op1, op2, op3, False)
>>> scaled.equiv(3.0 * result)
True
```

* <a id='id3'>**[1]**</a> R. McWeeny. Methods of Molecular Quantum Mechanics. 2nd Edition, Academic Press, 1992. ISBN 0-12-486552-6.
* <a id='id4'>**[2]**</a> P. J. Ollitrault, A. Miessen, I. Tavernelli. Molecular Quantum Dynamics: A Quantum Computing Perspective. Acc. Chem. Res. 54, 4229-4238 (2021). [https://doi.org/10.1021/acs.accounts.1c00514](https://doi.org/10.1021/acs.accounts.1c00514)
* **Parameters:**
  * **op_a** (*T*) – the operator $A$ above.
  * **op_b** (*T*) – the operator $B$ above.
  * **op_c** (*T*) – the operator $C$ above.
  * **sign** ([*bool*](https://docs.python.org/3/builtins/functions.html#bool)) – the nature of the outer (anti-)commutator as per the definition above.
* **Returns:**
  The double-commutator as per the definition above.
* **Return type:**
  *T*
