<a id="module-qiskit_bosons.operators.library"></a>

# Operator Library

Building on top of the operator representations provided by [`operators`](qiskit_bosons.operators.html.md#module-qiskit_bosons.operators), this
module provides tools to simplify the construction of commonly used operators.

Model Hamiltonians are built directly from your own lattice rather than from helpers here; the
[Bose-Hubbard guide](../guides/bose_hubbard.html.md#bose-hubbard-getting-started) works one through end to end.

## Commutator Generators

The brackets through which the bosonic algebra is expressed. The commutator is the defining one:
$[b_i, b_j^\dagger] = \delta_{ij}$ *is* the bosonic algebra, so these are how that structure
gets interrogated in practice (checking whether two terms commute, whether an operator conserves
boson number, and so on).

## Single-mode building blocks

The dimensionless quadrature operators, which are the pieces most other bosonic operators are
written in terms of and whose normalization is easy to get wrong.

These operators are Hermitian and satisfy the commutation relation:

$$
[\hat{q}_i, \hat{p}_j] = i\delta_{ij}
$$

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> from qiskit_bosons.operators.library import commutator, momentum, position
>>> commutator(position(0), momentum(0)).equiv(1j * BosonOperator.one())
True
>>> commutator(position(0), momentum(1)).equiv(BosonOperator.zero())
True
```
