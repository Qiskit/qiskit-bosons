<a id="operators-explanation"></a>

# Design principles of operator representations

This guide explains how a [`BosonOperator`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator) is organized and how to work with it: how to
build one, how to inspect and combine operators, and how to bring two operators into a form where
they can be compared. Start here when you write down a Hamiltonian.

A [`BosonOperator`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator) is an exact sum of products of bosonic creation and annihilation
operators with complex coefficients. It carries no occupation cutoff. A bosonic mode holds
arbitrarily many bosons, so its Fock space is infinite-dimensional, and a truncation enters only
when the operator is mapped onto qubits (see [Truncate the Fock space](mappers.html.md#truncation-explanation)). Everything in this
guide therefore happens in the exact algebra of the canonical commutation relations,

$$
[b_i, b_j^\dagger] = \delta_{ij}, \qquad [b_i, b_j] = [b_i^\dagger, b_j^\dagger] = 0 \, .
$$

Two consequences shape the interface. A mode index can repeat within a single term, because
$b_i^\dagger b_i^\dagger$ is a valid operator that creates two bosons. Also, reordering two
factors never introduces a sign, although moving an annihilation operator past a creation operator
on the same mode adds a term, $b_i b_i^\dagger = b_i^\dagger b_i + 1$.

## Storage

An operator keeps its terms in four flat arrays, similar to the compressed sparse row format of a
sparse matrix:

| `coeffs`     | One complex coefficient per term.                                                           |
|--------------|---------------------------------------------------------------------------------------------|
| `actions`    | One boolean per factor: `True` for a creation and `False` for an<br/>annihilation operator. |
| `modes`      | One mode index per factor.                                                                  |
| `boundaries` | The positions in `actions` and `modes` where each term starts and ends.                     |

The constructor takes these four arrays directly. The operator below has three terms, namely
$2 b_0^\dagger b_1$, $-i b_1^\dagger b_1^\dagger$, and the constant $0.5$, and you
can read the arrays back with the `get_*` methods:

```pycon
>>> from qiskit_bosons.operators import BosonOperator, ann, cre
>>> op = BosonOperator(
...     coeffs=[2.0, -1.0j, 0.5],
...     actions=[True, False, True, True],
...     modes=[0, 1, 1, 1],
...     boundaries=[0, 2, 4, 4],
... )
>>> op.get_coeffs()
[(2+0j), (-0-1j), (0.5+0j)]
>>> op.get_actions()
[True, False, True, True]
>>> op.get_modes()
[0, 1, 1, 1]
>>> op.get_boundaries()
[0, 2, 4, 4]
>>> len(op)
3
```

The first term spans positions 0 to 2 of `actions` and `modes`, the second spans 2 to 4, and the
third spans 4 to 4, which is how the constant term, a term without factors, is stored. The repeated
mode `1` in the second term is preserved as written; a term is never dropped or merged because a
mode index appears twice.

#### NOTE
An operator built with [`from_dict()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.from_dict) does not store its terms in the order
of the dictionary, so its arrays can list the terms in a different order from run to run. The
operator itself is the same, and nothing below that compares or prints operators depends on the
stored order.

The layout makes construction and arithmetic cheap, with one consequence: an operator can hold the
same term several times at once. Sums and products append terms and do not merge them, so call
[`simplify()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.simplify) when you need unique terms (see [Do arithmetic](#operators-arithmetic)).

## Build operators

The most direct way to build an operator is [`from_dict()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.from_dict), which takes a
dictionary from a tuple of factors to a coefficient. The functions [`cre()`](../stubs/qiskit_bosons.operators.cre.html.md#qiskit_bosons.operators.cre) and [`ann()`](../stubs/qiskit_bosons.operators.ann.html.md#qiskit_bosons.operators.ann)
create the factors, and the empty tuple is the constant term:

```pycon
>>> b = BosonOperator.from_dict({(ann(0),): 1.0})
>>> b_dag = BosonOperator.from_dict({(cre(0),): 1.0})
>>> number = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
>>> print(format(number))
  1.000000e0 +0.000000e0j * (+0 -0)
```

In the printed form, `+0` stands for $b_0^\dagger$ and `-0` for $b_0$. The
constructors [`zero()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.zero) and [`one()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.one) return the additive and
the multiplicative identity, and [`from_terms()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.from_terms) rebuilds an operator from the
iterator of [`iter_terms()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.iter_terms).

The operator library provides the dimensionless position and momentum quadratures,

$$
\hat{q} = \frac{b^\dagger + b}{\sqrt{2}}, \qquad
\hat{p} = \frac{i \left(b^\dagger - b\right)}{\sqrt{2}} \, ,
$$

as [`position()`](../stubs/qiskit_bosons.operators.library.position.html.md#qiskit_bosons.operators.library.position) and
[`momentum()`](../stubs/qiskit_bosons.operators.library.momentum.html.md#qiskit_bosons.operators.library.momentum). Together with the commutator functions of the
same module, you can check the algebra directly:

```pycon
>>> from qiskit_bosons.operators.library import commutator, momentum, position
>>> q, p = position(0), momentum(0)
>>> q.is_hermitian() and p.is_hermitian()
True
>>> commutator(q, p).equiv(1j * BosonOperator.one())
True
>>> commutator(b, b_dag).equiv(BosonOperator.one())
True
```

The module also provides [`anti_commutator()`](../stubs/qiskit_bosons.operators.library.anti_commutator.html.md#qiskit_bosons.operators.library.anti_commutator) and
[`double_commutator()`](../stubs/qiskit_bosons.operators.library.double_commutator.html.md#qiskit_bosons.operators.library.double_commutator). All three take operators of the same
type and return a new operator.

<a id="operators-product-order"></a>

### Mind the order of a product

Two operators can be multiplied in two ways, and they are easy to confuse. The `@` operator is the
mathematical product, so `a @ b` is $AB$. The `&` operator is *composition*, which follows
Qiskit’s convention of reading from left to right: `a & b` first applies `a` and then `b`, so
it is $BA$. The commutator functions use `@`, so `commutator(a, b)` is $[A, B]$.
Building the bracket yourself out of `&` silently negates it:

```pycon
>>> (b @ b_dag).equiv(b_dag @ b + BosonOperator.one())
True
>>> (b & b_dag).equiv(b_dag @ b)
True
>>> ((b @ b_dag) - (b_dag @ b)).equiv(BosonOperator.one())
True
>>> ((b & b_dag) - (b_dag & b)).equiv(-BosonOperator.one())
True
```

## Inspect an operator

An operator is not directly iterable, because its terms live in the flat arrays above. Use
[`iter_terms()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.iter_terms) instead, which yields a list of `(is_creation, mode)` pairs and
a coefficient for each term. The iterator returns copies, so changing them does not change the
operator.

```pycon
>>> hopping = BosonOperator.from_dict(
...     {
...         (cre(0), ann(1)): -1.0,
...         (cre(1), ann(0)): -1.0,
...         (cre(2), cre(2), ann(2), ann(2)): 0.5,
...     }
... )
>>> len(hopping)
3
>>> sorted(hopping.iter_terms())
[([(True, 0), (False, 1)], (-1+0j)), ([(True, 1), (False, 0)], (-1+0j)), ([(True, 2), (True, 2), (False, 2), (False, 2)], (0.5+0j))]
```

Several methods summarize the structure of an operator:

```pycon
>>> hopping.get_support()
{0, 1, 2}
>>> hopping.min_num_modes()
3
>>> hopping.max_rank()
4
```

[`get_support()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.get_support) is the set of modes the operator acts on,
[`min_num_modes()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.min_num_modes) is the smallest register that can hold it,
and [`max_rank()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.max_rank) is the length of its longest term.

#### NOTE
A **mode** is an abstract index that labels one bosonic degree of freedom. It has no meaning of
its own: you determine which site of a lattice, which vibrational coordinate or which optical port a mode
stands for. The [`BosonicCircuit`](../stubs/qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) uses the same indices (see
[Work with bosonic circuits](circuit.html.md#bosonic-circuit-explanation)), so an operator and the gates that carry it agree on what
mode `2` is without either one knowing what it represents. Nothing in the operator fixes how
many levels a mode has, which is also why the register size is only a lower bound;
[`min_num_modes()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.min_num_modes) is the smallest register that holds the operator,
not a property of it.

The occupation methods report how many bosons a single term can place in one mode. This is a lower
bound on the cutoff that a mapping needs. It is not a sufficient cutoff, because repeated
applications of the operator, as in a time evolution, can populate higher levels (see
[Truncate the Fock space](mappers.html.md#truncation-explanation)):

```pycon
>>> hopping.mode_occupations()
[1, 1, 2]
>>> hopping.mode_occupations([2, 0])
[2, 1]
>>> hopping.max_occupation()
2
```

Two further checks look at conserved quantities. An operator conserves the particle number when
every term has as many creation as annihilation operators. Passing a list of block sizes tests the
same condition within each block of consecutive modes:

```pycon
>>> hopping.conserves_particle_number()
True
>>> move = BosonOperator.from_dict({(cre(0), ann(2)): 1.0})
>>> move.conserves_sector([3])
True
>>> move.conserves_sector([2, 1])
False
```

Finally, [`relabel_modes()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.relabel_modes) renames modes: mode `i` becomes `permutation[i]`.
The list must be long enough to index every mode of the operator and must not repeat an entry:

```pycon
>>> print(format(move.relabel_modes([5, 6, 7])))
  1.000000e0 +0.000000e0j * (+5 -7)
```

<a id="operators-arithmetic"></a>

## Do arithmetic

Operators support addition, subtraction, negation, and multiplication and division by a number, in
the usual operator form or in place:

```pycon
>>> total = 2 * b_dag + b_dag - 3j * b
>>> len(total)
3
>>> len(total.simplify())
2
>>> total.simplify().equiv(3 * b_dag - 3j * b)
True
```

The three terms of `total` were never merged, so the operator held `b_dag` twice until
[`simplify()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.simplify) combined them. The method sums the coefficients of identical terms
and drops those whose total magnitude is below a tolerance. It is a structural operation: it
compares the written factors and does *not* apply the commutation relations. The related method
[`ichop()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.ichop) removes small coefficients in place without merging first, so call
[`simplify()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.simplify) before [`ichop()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.ichop) when an operator can hold duplicate terms.

Products use `@` or `&` as described above, and `**` raises an operator to an integer power:

```pycon
>>> x = b_dag + b
>>> print(format((x ** 2).simplify()))
  1.000000e0 +0.000000e0j * (-0 -0)
  1.000000e0 +0.000000e0j * (-0 +0)
  1.000000e0 +0.000000e0j * (+0 -0)
  1.000000e0 +0.000000e0j * (+0 +0)
```

The result lists the four orderings as written. Reducing them to a single form needs the
commutation relations, which the next section applies.

The adjoint reverses the order of the factors in each term, exchanges creation and annihilation
operators, and conjugates the coefficients. The method [`is_hermitian()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.is_hermitian) tests
whether an operator equals its adjoint:

```pycon
>>> hop = BosonOperator.from_dict({(cre(0), ann(1)): 1j})
>>> print(format(hop.adjoint()))
  0.000000e0 -1.000000e0j * (+1 -0)
>>> hop.is_hermitian()
False
>>> (hop + hop.adjoint()).is_hermitian()
True
```

A lone hopping term is not Hermitian, which matters when you group terms for a time evolution: each
group must be Hermitian on its own (see [Group operator terms: use the operator structure](grouping.html.md#grouping-explanation)).

## Normal order and compare operators

Two products of ladder operators can be the same operator and still look different, because
$b\,b^\dagger$ and $b^\dagger b + 1$ are equal in the algebra. Operators that are equal
in this sense need not be stored equally, so comparing them needs the commutation relations. The
method [`normal_ordered()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.normal_ordered) applies them: it moves every creation operator to the
left of every annihilation operator and adds a term whenever a pair on the same mode is exchanged.

```pycon
>>> product = BosonOperator.from_dict({(ann(0), cre(0), ann(0), cre(0)): 1.0})
>>> ordered = product.normal_ordered()
>>> len(ordered)
5
>>> print(format(ordered.simplify()))
  1.000000e0 +0.000000e0j * ()
  3.000000e0 +0.000000e0j * (+0 -0)
  1.000000e0 +0.000000e0j * (+0 +0 -0 -0)
```

The ordered operator is returned as an unmerged sum, with five stored terms that reduce to three
after [`simplify()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.simplify) merges them. Normal ordering never makes a term vanish:
$b^\dagger b^\dagger$ and every other repeated factor stay.

The `sandwich` argument sets how the modes inside the two groups of factors are arranged. By
default both groups are sorted in descending order. With `sandwich=True`, the creation operators
are ascending and the annihilation operators descending, which puts the larger indices towards the
middle:

```pycon
>>> pair = BosonOperator.from_dict({(cre(0), cre(1), ann(0), ann(1)): 1.0})
>>> list(pair.normal_ordered().iter_terms())
[([(True, 1), (True, 0), (False, 1), (False, 0)], (1+0j))]
>>> list(pair.normal_ordered(sandwich=True).iter_terms())
[([(True, 0), (True, 1), (False, 1), (False, 0)], (1+0j))]
```

The `==` operator tests exact equality of the stored terms, with no tolerance and no
simplification, so it can report a difference between two operators that are mathematically equal.
The method [`equiv()`](../stubs/qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator.equiv) brings the difference of the two operators to normal order
first and then tests that every coefficient is below a tolerance. It is the comparison you almost
always want:

```pycon
>>> reordered = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (): 1.0})
>>> contracted = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
>>> contracted == reordered
False
>>> contracted.simplify() == reordered.simplify()
False
>>> contracted.equiv(reordered)
True
```

The method `simplify` alone leaves the two forms apart, because it merges only identical written
terms. Likewise, [`is_hermitian()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.is_hermitian) compares an operator with its adjoint in
normal order, so it recognizes $b\,b^\dagger + b^\dagger b$ as Hermitian without further work.

#### NOTE
Normal ordering is an exact identity in the untruncated algebra, and `equiv` relies on that. It
is *not* an identity after a truncation. A mapped operator depends on the form you map, so the
mapped images of `contracted` and `reordered` above differ, on the highest retained level
only. See [Truncate the Fock space](mappers.html.md#truncation-explanation) for the size of the difference and how to choose a form.

## Order and filter terms

Reordering the *terms* of an operator, as opposed to the factors inside a term, does not change the
operator but can change the behavior of algorithms that read the terms in order, such as a
randomized product formula. The functions of [`qiskit_bosons.operators.terms`](../pydoc/qiskit_bosons.operators.terms.html.md#module-qiskit_bosons.operators.terms) return a new
operator and keep any group indices attached to their terms.

[`canonical_order()`](../stubs/qiskit_bosons.operators.terms.canonical_order.html.md#qiskit_bosons.operators.terms.canonical_order) sorts the terms into a fixed order that
depends only on the written factors, so an operator reads the same no matter how it was assembled.
[`order_terms()`](../stubs/qiskit_bosons.operators.terms.order_terms.html.md#qiskit_bosons.operators.terms.order_terms) sorts by a key that you supply, which receives
each term as the tuple that [`iter_terms()`](../stubs/qiskit_bosons.operators.OperatorTrait.html.md#qiskit_bosons.operators.OperatorTrait.iter_terms) yields:

```pycon
>>> from qiskit_bosons.operators.terms import canonical_order, order_terms
>>> mixed = BosonOperator.from_dict({(cre(1), ann(0)): 1.0, (cre(0), ann(1)): 3.0, (): 2.0})
>>> list(canonical_order(mixed).iter_terms())
[([], (2+0j)), ([(True, 0), (False, 1)], (3+0j)), ([(True, 1), (False, 0)], (1+0j))]
>>> by_weight = order_terms(mixed, key=lambda term: abs(term[1]), reverse=True)
>>> [coeff for _, coeff in by_weight.iter_terms()]
[(3+0j), (2+0j), (1+0j)]
```

If an operator carries group indices (see [Group operator terms: use the operator structure](grouping.html.md#grouping-explanation)), the function
[`group_order()`](../stubs/qiskit_bosons.operators.terms.group_order.html.md#qiskit_bosons.operators.terms.group_order) makes each group one contiguous run of terms:

```pycon
>>> from qiskit_bosons.operators.terms import group_order
>>> grouped = BosonOperator.from_dict({(cre(1), ann(0)): 1.0, (cre(0), ann(1)): 3.0, (): 2.0})
>>> grouped.groups = [1, 0, 1]
>>> group_order(grouped).groups
[0, 1, 1]
```

The last helper removes terms rather than reordering them. A term that is a product of number
operators, $b_i^\dagger b_i$, only multiplies each Fock state by a phase, so
[`filter_diagonal_terms()`](../stubs/qiskit_bosons.operators.terms.filter_diagonal_terms.html.md#qiskit_bosons.operators.terms.filter_diagonal_terms) can drop the constant, the number
operators, and their products when that phase is irrelevant for what you compute. It changes the
operator *in place* and requires a normal-ordered operator, which it does not check:

```pycon
>>> from qiskit_bosons.operators.terms import filter_diagonal_terms
>>> model = BosonOperator.from_dict(
...     {
...         (cre(0), ann(1)): -1.0,
...         (cre(1), ann(0)): -1.0,
...         (cre(0), cre(0), ann(0), ann(0)): 0.5,
...         (cre(0), ann(0)): 0.2,
...         (): 0.1,
...     }
... )
>>> filter_diagonal_terms(model)
>>> len(model)
2
>>> model.is_hermitian()
True
```

Only the two hopping terms remain. The on-site interaction and the number operator are diagonal in
the occupation basis, so they are removed together with the constant. The repeated mode in the
interaction is handled correctly: $b_0^\dagger b_0^\dagger b_0 b_0$ is diagonal, while
$b_0^\dagger b_0^\dagger b_0$ would not be.

## Next steps

- [Group operator terms: use the operator structure](grouping.html.md#grouping-explanation) shows how to tag terms with group indices.
- [Map operator representations](mappers.html.md#mappers-explanation) turns an operator into a qubit operator, including the choice of
  encoding and the effect of a truncation.
- [Work with bosonic circuits](circuit.html.md#bosonic-circuit-explanation) places an operator into a circuit as an
  [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate.
- [Simulate Bose-Hubbard dynamics with a choice of encoding](bose_hubbard.html.md#bose-hubbard-getting-started) builds the Hamiltonian of a Bose-Hubbard chain and simulates
  it end to end.
- The [`qiskit_bosons.operators`](../pydoc/qiskit_bosons.operators.html.md#module-qiskit_bosons.operators) and [`qiskit_bosons.operators.terms`](../pydoc/qiskit_bosons.operators.terms.html.md#module-qiskit_bosons.operators.terms) modules document the
  complete interface.
