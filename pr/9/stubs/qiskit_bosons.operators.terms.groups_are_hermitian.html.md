# groups_are_hermitian

### groups_are_hermitian(op, , atol=1e-08)

Returns, for each group, whether the operator formed by that group’s terms is Hermitian.

Groups prescribe no meaning of their own; this checks one *common* convention, namely that each
group is separately Hermitian. That is the property a randomized product formula relies on when it
samples whole groups, since only a Hermitian group has a unitary time evolution. It is equally
what Qiskit’s [`PauliEvolutionGate`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.PauliEvolutionGate) requires of each group once
mapped: a lone hopping term `b†_i b_j` is not Hermitian, so it has to share a group with its
adjoint. If the operator tracks no groups, this returns `None`.

An empty group counts as Hermitian, because the zero operator is.

#### NOTE
This is *not* implied by (nor does it imply) [`groups_have_uniform_coeffs()`](qiskit_bosons.operators.terms.groups_have_uniform_coeffs.html.md#qiskit_bosons.operators.terms.groups_have_uniform_coeffs). Uniform
coefficients do not make a group Hermitian, and a Hermitian group may mix magnitudes.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> from qiskit_bosons.operators.terms.grouping import groups_are_hermitian
>>> # built from the sparse arrays rather than a dictionary, because only the former
>>> # fixes the term order that the group indices are paired with
>>> op = BosonOperator(
...     [1.0, 1.0, 1.0],
...     [True, False, True, False, True, False],
...     [0, 1, 1, 0, 2, 3],
...     [0, 2, 4, 6],
... )
>>> op.groups = [0, 0, 1]  # the leading hopping pair, then one unpaired term
>>> groups_are_hermitian(op)
[True, False]
```

* **Parameters:**
  * **op** – the operator whose groups to check.
  * **atol** – The numerical accuracy upto which coefficients are considered equal. This value
    defaults to `1e-8`.
* **Returns:**
  One flag per group index, or `None` if the operator tracks no groups.
* **Raises:**
  [**TypeError**](https://docs.python.org/3/builtins/exceptions.html#TypeError) – if `op` is not a supported operator type (see
      [`OperatorTrait`](qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait)).
