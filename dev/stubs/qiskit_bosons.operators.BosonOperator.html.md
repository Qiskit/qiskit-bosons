# BosonOperator

### *class* BosonOperator(coeffs, actions, modes, boundaries)

A bosonic operator.

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

## Definition

This operator is defined by a linear combination of products of bosonic creation and
annihilation operators. The individual factors fulfill the canonical **commutation**
relations: <sup>[1](#id2)</sup>

$$
\left[b^\dagger_i, b^\dagger_j\right] =
\left[b_i, b_j\right] = 0,~~\text{and}~~
\left[b_i, b^\dagger_j\right] = \delta_{ij} \, .
$$

This makes the definition of the entire operator the following:

$$
\text{\texttt{BosonOperator}} = \sum_i c_i \bigotimes_j \hat{B_j} \, ,
$$

where $\hat{B_j} \in \{ b_j, b^\dagger_j \}$ and $c_i$ is the (complex) coefficient
making up the linear combination of products. The index $j$ can take any value between 0
and the number of bosonic modes acted upon by the operator minus 1.

### Algebra

Two consequences of the commutation relations above shape this entire class:

* **A mode index may repeat.** $b^\dagger_i b^\dagger_i \neq 0$: a mode can hold
  arbitrarily many bosons, so a mode index may legitimately repeat within a single term. Terms
  with repeated indices are never silently dropped.
* **Reordering never introduces a sign.** Swapping two factors is free unless they are a
  creation/annihilation pair on the *same* mode, in which case
  $b_i b^\dagger_i = b^\dagger_i b_i + 1$. See [`normal_ordered()`](#qiskit_bosons.operators.BosonOperator.normal_ordered).

### Truncation

A bosonic mode has an **infinite-dimensional** Fock space. This operator is exact and carries no
occupation cutoff; truncation is introduced only when mapping to qubits, where it is a required,
explicit parameter. See [`mappers`](../pydoc/qiskit_bosons.mappers.html.md#module-qiskit_bosons.mappers) and [`max_occupation()`](#qiskit_bosons.operators.BosonOperator.max_occupation).

<a id="bosonoperator-implementation"></a>

## Implementation

This class stores the terms and coefficients in multiple sparse vectors, akin to the
[compressed sparse row format](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format))
commonly used for sparse matrices. More concretely, a single operator contains 4 arrays:

| `coeffs`     | A vector of complex coefficients consisting of two 64-bit floating point numbers.   |
|--------------|-------------------------------------------------------------------------------------|
| `actions`    | A vector of booleans storing the nature of the second-quantization actions.         |
| `modes`      | A vector of 32-bit integers storing the bosonic mode indices acted upon.            |
| `boundaries` | A vector of integers indicating the boundaries in `actions` and `modes`.            |

Entries in `actions` indicate creation (annihilation) operators by `True` (`False`).

#### NOTE
The same mode index may appear multiple times within one term (e.g. `+0 +0` for
$b^\dagger_0 b^\dagger_0$). Such repetitions are meaningful and are preserved verbatim.

#### NOTE
You may access **read-only copies** of these internal arrays via their respective methods:
[`get_coeffs()`](#qiskit_bosons.operators.BosonOperator.get_coeffs), [`get_actions()`](#qiskit_bosons.operators.BosonOperator.get_actions), [`get_modes()`](#qiskit_bosons.operators.BosonOperator.get_modes), and [`get_boundaries()`](#qiskit_bosons.operators.BosonOperator.get_boundaries).

This data structure allows for very efficient construction and manipulation of operators.
However, it implies that duplicate terms may be contained in an operator at any moment.
These must be resolved manually through the use of [`simplify()`](#qiskit_bosons.operators.BosonOperator.simplify).

### Construction

An operator can be constructed directly by providing the arrays outlined above:

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> coeffs = [1.0, 2.0, -3.0, 4.0j, -0.5j]
>>> actions = [True, False, False, True, True, True, False, False]
>>> modes = [0, 0, 0, 1, 0, 1, 2, 3]
>>> boundaries = [0, 0, 1, 2, 4, 8]
>>> op = BosonOperator(coeffs, actions, modes, boundaries)
>>> print(format(op))
  1.000000e0 +0.000000e0j * ()
 -3.000000e0 +0.000000e0j * (-0)
  0.000000e0 +4.000000e0j * (-0 +1)
  2.000000e0 +0.000000e0j * (+0)
 -0.000000e0-5.000000e-1j * (+0 +1 -2 -3)
```

For convenience, it is possible to construct an operator from a Python dictionary like so:

```pycon
>>> from qiskit_bosons.operators import cre, ann
>>> op = BosonOperator.from_dict(
...     {
...         (): 1.0,
...         (cre(0),): 2.0,
...         (ann(0),): -3.0,
...         (ann(0), cre(1)): 4.0j,
...         (cre(0), cre(1), ann(2), ann(3)): -0.5j,
...     }
... )
>>> print(format(op))
  1.000000e0 +0.000000e0j * ()
 -3.000000e0 +0.000000e0j * (-0)
  0.000000e0 +4.000000e0j * (-0 +1)
  2.000000e0 +0.000000e0j * (+0)
 -0.000000e0-5.000000e-1j * (+0 +1 -2 -3)
```

In this example, we have leveraged [`cre()`](qiskit_bosons.operators.cre.html.md#qiskit_bosons.operators.cre) and [`ann()`](qiskit_bosons.operators.ann.html.md#qiskit_bosons.operators.ann) for creating the creation and
annihilation operators at the specified modes.

In addition, the following construction and quick helper methods are available:

### Formatting

In the examples above, the constructed operators have been printed using the output from
[`format()`](https://docs.python.org/3/builtins/functions.html#format), which results in a human-readable form of the operator.

```pycon
>>> print(format(op))
  1.000000e0 +0.000000e0j * ()
 -3.000000e0 +0.000000e0j * (-0)
  0.000000e0 +4.000000e0j * (-0 +1)
  2.000000e0 +0.000000e0j * (+0)
 -0.000000e0-5.000000e-1j * (+0 +1 -2 -3)
```

#### NOTE
The printing order of `format(op)` gets explicitly sorted before printing. As such, it
does not reflect the order of the terms inside the operator.

An alternative form can be obtained from the [`repr()`](https://docs.python.org/3/builtins/functions.html#repr) function, which results in a
Python-interpretable representation. In other words, this output can readily be copied and
pasted into a Python shell:

```pycon
>>> print(repr(op))
BosonOperator.from_dict({...})
```

Finally, for large operators both of these outputs may be very long and undesirable. Then, a
very simple form with minimal information can be obtained from the `str()` function:

```pycon
>>> print(str(op))
<BosonOperator with 5 terms>
```

### Iteration

Since the underlying data structure is implemented in Rust and has a non-trivial layout, it
cannot be iterated over directly:

```pycon
>>> list(iter(op))
Traceback (most recent call last):
  ...
TypeError: 'qiskit_bosons.operators.boson_operator.BosonOperator' object is not iterable
```

Instead, this class provides custom iterators to fulfill this purpose:

```pycon
>>> list(sorted(op.iter_terms()))
[([], (1+0j)), ([(False, 0)], (-3+0j)), ([(False, 0), (True, 1)], 4j), ([(True, 0)], (2+0j)), ([(True, 0), (True, 1), (False, 2), (False, 3)], (-0-0.5j))]
```

#### SEE ALSO
[`iter_terms()`](#qiskit_bosons.operators.BosonOperator.iter_terms)
: For more relevant implementation details.

The table below lists all available iterators:

### Arithmetics

The following arithmetic operations are supported:

#### Addition/Subtraction

```pycon
>>> op = BosonOperator.one()
>>> (op + op).simplify()
BosonOperator.from_dict({(): 2+0j})
>>> (op - op).simplify()
BosonOperator.from_dict({})
>>> op += op
>>> op.simplify()
BosonOperator.from_dict({(): 2+0j})
>>> op -= op
>>> op.simplify()
BosonOperator.from_dict({})
```

#### Scalar Multiplication/Divison

```pycon
>>> op = BosonOperator.one()
>>> (2 * op).simplify()
BosonOperator.from_dict({(): 2+0j})
>>> (op / 2).simplify()
BosonOperator.from_dict({(): 0.5+0j})
>>> op *= 2
>>> op.simplify()
BosonOperator.from_dict({(): 2+0j})
>>> op /= 2
>>> op.simplify()
BosonOperator.from_dict({(): 1+0j})
```

#### Operator Composition

#### NOTE
Operator composition corresponds to left-multiplication: `c = a & b` corresponds to
$C = B A$. In other words, the composition of two operators returns a resulting
operator that performs “first `a` and then `b`”.

```pycon
>>> op1 = BosonOperator.from_dict({(): 2.0, (cre(0),): 3.0})
>>> op2 = BosonOperator.from_dict({(): 1.5, (ann(1),): 4.0})
>>> comp = (op1 & op2).simplify()
>>> print(format(comp))
  3.000000e0 +0.000000e0j * ()
  8.000000e0 +0.000000e0j * (-1)
  1.200000e1 +0.000000e0j * (-1 +0)
  4.500000e0 +0.000000e0j * (+0)
>>> op2 &= op1
>>> print(format(op2.simplify()))
  3.000000e0 +0.000000e0j * ()
  8.000000e0 +0.000000e0j * (-1)
  4.500000e0 +0.000000e0j * (+0)
  1.200000e1 +0.000000e0j * (+0 -1)
>>> squared = (op1 ** 2).simplify()
>>> print(format(squared))
  4.000000e0 +0.000000e0j * ()
  1.200000e1 +0.000000e0j * (+0)
  9.000000e0 +0.000000e0j * (+0 +0)
```

#### NOTE
For convenience, the right-multiplication is implemented by `c = a @ b` (resulting in
$C = A B$).

```pycon
>>> (op1 @ op2).equiv(op2 & op1)
True
```

#### Other Operations

In addition to the magic methods that correspond to the arithmetic operations outlined above,
the following methods are available:

#### Properties

Finally, various methods exist to check certain properties of an operator:

* <a id='id2'>**[1]**</a> [https://en.wikipedia.org/wiki/Second_quantization#Boson_creation_and_annihilation_operators](https://en.wikipedia.org/wiki/Second_quantization#Boson_creation_and_annihilation_operators)

### Attributes

#### groups

An optional vector of group indices for each term.

For more information refer to the [Group operator terms: use the operator structure](../guides/grouping.html.md#grouping-explanation) guide.

### Methods

#### adjoint()

Returns the Hermitian conjugate (or adjoint) of this operator.

This affects the terms and coefficients as follows:

- the actions in each term reverse their order and flip between creation and annihilation
- the coefficients are complex conjugated

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({(): -1.0j, ((True, 0), (False, 1)): 1.0})
>>> adj = op.adjoint()
>>> print(format(adj))
 -0.000000e0 +1.000000e0j * ()
  1.000000e0 -0.000000e0j * (+1 -0)
```

#### conserves_particle_number()

Returns whether this operator is particle-number conserving.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({((True, 0), (False, 1)): 1})
>>> op.conserves_particle_number()
True
>>> op = BosonOperator.from_dict({((True, 0),): 1})
>>> op.conserves_particle_number()
False
```

* **Returns:**
  Whether this operator is particle-number conserving.

#### conserves_sector(block_sizes)

Returns whether every term conserves particle number within each mode block.

`block_sizes` partitions the mode range into consecutive, non-overlapping blocks: block
`b` spans modes `[start_b, start_b + block_sizes[b])` where `start_b` is the sum of the
preceding block sizes. A term conserves the sector if and only if, in *every* block, its
number of creation operators equals its number of annihilation operators. A term acting on a
mode beyond the final block does not conserve the sector.

An empty `block_sizes` treats all modes as a single block, making this equivalent to
[`conserves_particle_number()`](#qiskit_bosons.operators.BosonOperator.conserves_particle_number).

Blocks are useful for structured bosonic systems. A vibrational Hamiltonian in the
Christiansen $n$-mode form, for instance, assigns one block of modals to each physical
vibrational mode, and a physical state carries exactly one excitation per block; checking
per-block conservation verifies that an operator respects that structure.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({((True, 0), (False, 2)): 1})
>>> op.conserves_sector([4])  # one block spanning all 4 modes
True
>>> op.conserves_sector([2, 2])  # moves a boson from the first block to the second
False
```

* **Parameters:**
  **block_sizes** – the sizes of the consecutive mode blocks that each must be conserved.
* **Returns:**
  Whether every term conserves particle number within each mode block.

#### equiv(other, atol=1e-08)

Checks this operator for equivalence with another operator.

Equivalence in this context means approximate equality up to the specified absolute
tolerance. To be more precise, this method returns `True`, when all the absolute values
of the coefficients in the difference `other - self` are below the specified threshold
`atol`.

#### NOTE
This is the mathematical comparison you almost always want. It differs from the `==`
operator, which tests exact equality of the *stored* terms (their coefficients, actions,
modes, and internal term boundaries) with no tolerance and no simplification. Two
mathematically equal operators can therefore compare unequal under `==` if they are
stored differently – for example an unsimplified `a + a` versus `2 * a`, or terms
held in a different order. Use `equiv` to compare operators up to numerical tolerance.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({(): 1e-7})
>>> zero = BosonOperator.zero()
>>> op.equiv(zero)
False
>>> op.equiv(zero, 1e-6)
True
>>> op.equiv(zero, 1e-9)
False
```

* **Parameters:**
  * **other** – the other operator to compare with.
  * **atol** – the absolute tolerance for the comparison. This value defaults to `1e-8`.

#### *classmethod* from_dict(data)

Constructs a new operator from a dictionary.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict(
...     {
...         (): 1.0-1.0j,
...         ((True, 0), (False, 1)): 2.0,
...     }
... )
>>> print(format(op))
  1.000000e0 -1.000000e0j * ()
  2.000000e0 +0.000000e0j * (+0 -1)
```

#### NOTE
The order of the terms inside the resulting operator is not guaranteed, and it can differ
between runs. The printing order of `format(op)` is sorted explicitly, so it does not
reveal it. Use [`canonical_order()`](qiskit_bosons.operators.terms.canonical_order.html.md#qiskit_bosons.operators.terms.canonical_order) if you need a fixed term order.

* **Parameters:**
  **data** – a dictionary mapping tuples of terms to complex coefficients. Each key is a tuple
  of `(bool, int)` pairs. You may use [`cre()`](qiskit_bosons.operators.cre.html.md#qiskit_bosons.operators.cre) and [`ann()`](qiskit_bosons.operators.ann.html.md#qiskit_bosons.operators.ann) to simplify
  their construction.
* **Returns:**
  A new operator.

#### *classmethod* from_terms(terms)

Constructs a new operator from an iterator of terms (see also [`iter_terms()`](#qiskit_bosons.operators.BosonOperator.iter_terms)).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({(): 2.0, ((True, 0),): 1.0, ((False, 1),): -1.0j})
>>> op.equiv(BosonOperator.from_terms(op.iter_terms()))
True
```

* **Parameters:**
  **terms** – an iterator of terms as produced by [`iter_terms()`](#qiskit_bosons.operators.BosonOperator.iter_terms).
* **Returns:**
  A new operator.

#### *classmethod* from_terms_with_groups(terms)

Constructs a new operator from an iterator of terms with groups (see also
[`iter_terms_with_groups()`](#qiskit_bosons.operators.BosonOperator.iter_terms_with_groups)).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator(
...     [2.0, 1.0, -1.0],
...     [True, False, True, False],
...     [0, 1, 1, 0],
...     [0, 0, 2, 4],
... )
>>> op.groups = [0, 1, 1]
>>> reconstructed = BosonOperator.from_terms_with_groups(op.iter_terms_with_groups())
>>> op.equiv(reconstructed) and op.groups == reconstructed.groups
True
```

* **Parameters:**
  **terms** – an iterator of terms as produced by [`iter_terms_with_groups()`](#qiskit_bosons.operators.BosonOperator.iter_terms_with_groups).
* **Returns:**
  A new operator.

#### get_actions()

Returns a read-only list of the operator’s actions.

#### NOTE
This method returns a **copy** of the internal data.

#### SEE ALSO
The explanation of the internal data structure,
[here](#bosonoperator-implementation).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.one()
>>> op += BosonOperator.from_dict({((True, 0), (False, 1)): 1.0})
>>> op.get_actions()
[True, False]
```

* **Returns:**
  A list of the operator’s actions.

#### get_boundaries()

Returns a read-only list of the indices indicating the boundaries between operator terms.

#### NOTE
This method returns a **copy** of the internal data.

#### SEE ALSO
The explanation of the internal data structure,
[here](#bosonoperator-implementation).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.one()
>>> op += BosonOperator.from_dict({((True, 0), (False, 1)): 1.0})
>>> op.get_boundaries()
[0, 0, 2]
```

* **Returns:**
  A list of the operator’s terms boundaries.

#### get_coeffs()

Returns a read-only list of the operator’s coefficients.

#### NOTE
This method returns a **copy** of the internal data.

#### SEE ALSO
The explanation of the internal data structure,
[here](#bosonoperator-implementation).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.one()
>>> op += -1j * BosonOperator.one()
>>> op.get_coeffs()
[(1+0j), -1j]
```

* **Returns:**
  A list of the operator’s coefficients.

#### get_modes()

Returns a read-only list of the operator’s acted-upon mode indices.

#### NOTE
This method returns a **copy** of the internal data.

#### SEE ALSO
The explanation of the internal data structure,
[here](#bosonoperator-implementation).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.one()
>>> op += BosonOperator.from_dict({((True, 0), (False, 1)): 1.0})
>>> op.get_modes()
[0, 1]
```

* **Returns:**
  A list of the operator’s modes.

#### get_support()

Returns the set of mode indices which this operator acts upon.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict(
...     {
...         ((True, 0), (False, 4)): 1,
...         ((True, 1), (True, 3), (False, 4), (False, 7)): 1,
...     }
... )
>>> assert op.get_support() == {0, 1, 3, 4, 7}
```

* **Returns:**
  The set of mode indices which this operator acts upon.

#### has_groups()

Returns whether this operator tracks group indices.

This is equivalent to (but cheaper than) checking `op.groups is not None`, because it does
not copy the group indices out of the operator in order to inspect them.

#### NOTE
This returns `True` even when [`groups`](#qiskit_bosons.operators.BosonOperator.groups) is an empty list, which is the state of a
grouped operator that holds no terms.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator(
...     [1.0, 2.0, -1.0, -2.0],
...     [True, False, True, False, True, False, True, False],
...     [0, 1, 2, 3, 1, 0, 3, 2],
...     [0, 2, 4, 6, 8],
... )
>>> op.has_groups()
False
>>> op.groups = [0, 1, 0, 1]
>>> op.has_groups()
True
```

* **Returns:**
  Whether [`groups`](#qiskit_bosons.operators.BosonOperator.groups) is set on this operator.

#### ichop(atol=1e-08)

Removes terms whose coefficient magnitude lies below the provided threshold.

This method modifies the operator *in place* and returns `None`.

#### CAUTION
This method truncates coefficients greedily! If the acted upon operator may contain
separate coefficients for duplicate terms consider calling [`simplify()`](#qiskit_bosons.operators.BosonOperator.simplify) instead!

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({(): 1e-4, ((True, 0),): 1e-6, ((False, 0),): 1e-10})
>>> print(format(op))
  1.000000e-4 +0.000000e0j * ()
 1.000000e-10 +0.000000e0j * (-0)
  1.000000e-6 +0.000000e0j * (+0)
>>> op.ichop()
>>> print(format(op))
  1.000000e-4 +0.000000e0j * ()
  1.000000e-6 +0.000000e0j * (+0)
>>> op.ichop(1e-5)
>>> print(format(op))
  1.000000e-4 +0.000000e0j * ()
```

* **Parameters:**
  **atol** – the absolute tolerance for the cutoff. This value defaults to `1e-8`.

#### is_hermitian(atol=1e-08)

Returns whether this operator is Hermitian.

#### NOTE
This check is implemented using [`equiv()`](#qiskit_bosons.operators.BosonOperator.equiv) on the [`normal_ordered()`](#qiskit_bosons.operators.BosonOperator.normal_ordered) difference
of `self` and its [`adjoint()`](#qiskit_bosons.operators.BosonOperator.adjoint) and [`zero()`](#qiskit_bosons.operators.BosonOperator.zero).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({
...     ((True, 0), (False, 1)): 1.00001j,
...     ((True, 1), (False, 0)): -1j,
... })
>>> op.is_hermitian()
False
>>> op.is_hermitian(1e-4)
True
```

* **Parameters:**
  **atol** – The numerical accuracy upto which coefficients are considered equal. This value
  defaults to `1e-8`.
* **Returns:**
  Whether this operator is Hermitian.

#### iter_terms()

An iterator over the operator’s terms.

#### WARNING
Mutating the iteration items does **not** affect the underlying operator data.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({(): 2.0, ((True, 0),): 1.0, ((False, 1),): -1.0j})
>>> list(sorted(op.iter_terms()))
[([], (2+0j)), ([(False, 1)], (-0-1j)), ([(True, 0)], (1+0j))]
```

#### iter_terms_with_groups()

An iterator over the operator’s terms with their associated group index.

#### WARNING
Mutating the iteration items does **not** affect the underlying operator data.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator(
...     [2.0, 1.0, -1.0],
...     [True, False, True, False],
...     [0, 1, 1, 0],
...     [0, 0, 2, 4],
... )
>>> op.groups = [0, 1, 1]
>>> list(sorted(op.iter_terms_with_groups()))
[([], (2+0j), 0), ([(True, 0), (False, 1)], (1+0j), 1), ([(True, 1), (False, 0)], (-1+0j), 1)]
```

#### max_occupation()

Returns the highest occupation any single mode is driven to by any one term.

This is the largest number of creation operators acting on one and the same mode within a
single term, which is the maximum of [`mode_occupations()`](#qiskit_bosons.operators.BosonOperator.mode_occupations) over all modes. It measures how
far one application of a term can raise a mode, so it is a **lower bound** on the truncation
level a boson-to-qubit encoding needs in order to represent this operator without clipping
only when the state starts in the vacuum: applying $b^\dagger_i b^\dagger_i$ to the
vacuum populates occupation 2, which a cutoff of `max_occupation < 2` cannot hold. For a
state with $n$ bosons in a mode, the mode can reach $n$ plus this value.

```pycon
>>> from qiskit_bosons.operators import BosonOperator, ann, cre
>>> hopping = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
>>> hopping.max_occupation()
1
>>> interaction = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
>>> interaction.max_occupation()
2
```

#### WARNING
This is only a lower bound, **not** a sufficient cutoff. The occupation reachable by
*repeated* application of the operator – as happens over a time evolution – is
unbounded, and a mode rotation moves amplitude between occupation sectors, so acting on
$\lvert 1, 1 \rangle$ can populate $\lvert 2, 0 \rangle$. Choosing the cutoff
remains a physical modeling decision; this only rules out cutoffs that are certainly
too small.

* **Returns:**
  The largest per-mode creation count over all terms.

#### max_rank()

Returns the maximum rank of the terms in this operator.

#### NOTE
The length of the longest term can depend on the operator’s form which means that (for
example) operator simplification or normal-ordering can result in a different maximum
rank.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({
...     ((True, 0), (False, 1), (True, 2), (False, 3)): 1,
... })
>>> op.max_rank()
4
```

* **Returns:**
  The maximum rank of this operator.

#### min_num_modes()

Returns the smallest number of modes that can host this operator.

This is one more than the largest mode index the operator acts on, or `0` if it acts on
no mode.

```pycon
>>> from qiskit_bosons.operators import BosonOperator, ann, cre
>>> op = BosonOperator.from_dict({(cre(2), ann(5)): 1.0})
>>> op.min_num_modes()
6
>>> BosonOperator.zero().min_num_modes()
0
```

* **Returns:**
  The minimum number of modes.

#### mode_occupations(modes=None)

Returns, for each mode, the highest occupation it is driven to by any one term.

The occupation of a mode in a term is the number of creation operators acting on that mode
within the term. Annihilation operators are not subtracted, so the entry is how far one
application of a term can **raise** the mode, whatever the term does to it afterwards.
Because it measures a rise, it is a lower bound on the cutoff of the mode only when the mode
starts empty. A mode that starts in occupation $n$ can reach $n$ plus the entry,
so add the initial occupation of the state when you pick a cutoff. Like
[`max_occupation()`](#qiskit_bosons.operators.BosonOperator.max_occupation), the result is not a sufficient cutoff.

```pycon
>>> from qiskit_bosons.operators import BosonOperator, ann, cre
>>> op = BosonOperator.from_dict({
...     (cre(0), cre(0), ann(3)): 1.0,
...     (cre(3), cre(3), cre(3)): 1.0,
... })
>>> op.mode_occupations()
[2, 0, 0, 3]
>>> op.mode_occupations([3, 1])
[3, 0]
```

* **Parameters:**
  **modes** – the modes to report, in the order to report them. If omitted (`None`), every
  mode from `0` to the largest one the operator acts on is reported, with zero for
  the modes it does not touch. A requested mode the operator does not touch,
  including one beyond [`min_num_modes()`](#qiskit_bosons.operators.BosonOperator.min_num_modes), reports zero.
* **Returns:**
  One bound per reported mode.

#### normal_ordered(sandwich=None)

Returns an equivalent operator with normal ordered terms.

The normal order of an operator term is defined such that all creation actions appear
before all annihilation actions.
Within each group, the acted-upon modes are ordered lexicographically. Whether their order
is ascending or descending depends upon the value of the `sandwich` argument:

- `None` (default): both groups are ordered lexicographically descending (e.g.
  `+1 +0 -1 -0`)
- `True`: larger indices appear towards the middle, i.e. creation actions are
  lexicographically ascending while annihilation ones are descending (e.g.
  `+0 +1 -1 -0`)
- `False`: smaller indices appear towards the middle, i.e. creation actions are
  lexicographically descending while annihilation ones are ascending (e.g.
  `+1 +0 -0 -1`)

#### NOTE
When a term is being reordered, the canonical **commutation** relations have to be taken
into account, $b_i b^\dagger_j = \delta_{ij} + b^\dagger_j b_i$, implying that the
number of terms may change. Note the sign: the contracted term is *added*. Reordering
itself never introduces a sign, and no term ever vanishes – $b^\dagger_i
b^\dagger_i$ is a perfectly good operator.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({((False, 1), (True, 1), (False, 0), (True, 0)): 1})
>>> print(format(op.normal_ordered().simplify()))
  1.000000e0 +0.000000e0j * ()
  1.000000e0 +0.000000e0j * (+0 -0)
  1.000000e0 +0.000000e0j * (+1 -1)
  1.000000e0 +0.000000e0j * (+1 +0 -1 -0)
>>> print(format(op.normal_ordered(sandwich=True).simplify()))
  1.000000e0 +0.000000e0j * ()
  1.000000e0 +0.000000e0j * (+0 -0)
  1.000000e0 +0.000000e0j * (+0 +1 -1 -0)
  1.000000e0 +0.000000e0j * (+1 -1)
>>> print(format(op.normal_ordered(sandwich=False).simplify()))
  1.000000e0 +0.000000e0j * ()
  1.000000e0 +0.000000e0j * (+0 -0)
  1.000000e0 +0.000000e0j * (+1 -1)
  1.000000e0 +0.000000e0j * (+1 +0 -0 -1)
```

* **Returns:**
  An equivalent but normal-ordered operator.

#### num_groups()

Returns the number of groups.

If [`groups`](#qiskit_bosons.operators.BosonOperator.groups) is `None`, this function also returns `None`. Otherwise, it will
return the number of groups which is defined to be the largest occurring group index plus
1 (which may therefore be used as the index for the next group).

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator(
...     [1.0, 2.0, -1.0, -2.0],
...     [True, False, True, False, True, False, True, False],
...     [0, 1, 2, 3, 1, 0, 3, 2],
...     [0, 2, 4, 6, 8],
... )
>>> op.groups = [0, 1, 0, 1]
>>> op.num_groups()
2
```

* **Returns:**
  The largest group index in [`groups`](#qiskit_bosons.operators.BosonOperator.groups) plus 1.

#### *classmethod* one()

Constructs the multiplicative identity operator.

Composing the operator that is constructed by this method with another one has no effect.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({(): 2.0})
>>> one = BosonOperator.one()
>>> op & one == op
True
```

#### relabel_modes(permutation)

Returns a new operator with relabeled modes.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({
...     ((True, 0), (False, 1)): 1,
...     ((True, 0), (False, 1), (True, 2), (False, 3)): 1,
... })
>>> permutation = [5, 6, 4, 3]
>>> relabeled = op.relabel_modes(permutation)
>>> print(format(relabeled))
  1.000000e0 +0.000000e0j * (+5 -6)
  1.000000e0 +0.000000e0j * (+5 -6 +4 -3)
```

* **Parameters:**
  **permutation** – the index permutation list. Mode `i` is relabeled to `permutation[i]`,
  so the list must contain no duplicate entries and must be long enough to index every
  mode the operator acts upon (its length must exceed the operator’s largest mode
  index).
* **Returns:**
  A new operator with its modes relabeled.
* **Raises:**
  [**ValueError**](https://docs.python.org/3/builtins/exceptions.html#ValueError) – if `permutation` contains duplicate entries, or is too short to relabel
      some mode the operator acts upon.

#### simplify(atol=1e-08)

Returns an equivalent but simplified operator.

The simplification process first sums all coefficients that belong to equal terms and then
only retains those whose total coefficient exceeds the specified tolerance (just like
[`ichop()`](#qiskit_bosons.operators.BosonOperator.ichop)).

When an operator has been arithmetically manipulated or constructed in a way that does not
guarantee unique terms, this method should be called before applying any method that
filters numerically small coefficients to avoid loss of information. See the example below
which showcases how [`ichop()`](#qiskit_bosons.operators.BosonOperator.ichop) can truncate terms that sum to a total coefficient
magnitude which should not be truncated:

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> coeffs = [1e-5] * int(1e5)
>>> boundaries = [0] + [0] * int(1e5)
>>> op = BosonOperator(coeffs, [], [], boundaries)
>>> canon = op.simplify(1e-4)
>>> assert canon.equiv(op.one(), 1e-6)
>>> op.ichop(1e-4)
>>> assert op.equiv(op.zero(), 1e-6)
```

* **Parameters:**
  **atol** – the absolute tolerance for the cutoff. This value defaults to `1e-8`.
* **Returns:**
  An equivalent but simplified operator.

#### split_out_groups(group_indices=None)

Splits this operator into an optional list of new operators based on [`groups`](#qiskit_bosons.operators.BosonOperator.groups).

If [`groups`](#qiskit_bosons.operators.BosonOperator.groups) is `None`, this function also returns `None`. Otherwise, if
`group_indices` is `None` (the default), it returns a list of one new operator for
every group index in [`groups`](#qiskit_bosons.operators.BosonOperator.groups), in index order. If `group_indices` is given, only
the requested indices are built, in the given order: this avoids the cost of constructing
operators for groups that are never used, which is especially beneficial when only a small
number of groups out of a much larger total are needed, e.g. when subsampling groups for a
randomized product formula. A duplicate index in `group_indices` is returned once per
occurrence.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator(
...     [1.0, 2.0, -1.0, -2.0],
...     [True, False, True, False, True, False, True, False],
...     [0, 1, 2, 3, 1, 0, 3, 2],
...     [0, 2, 4, 6, 8],
... )
>>> print(op.split_out_groups())
None
>>> op.groups = [0, 1, 0, 1]
>>> groups = op.split_out_groups()
>>> for g in groups:
...     print(list(sorted(g.iter_terms())))
[([(True, 0), (False, 1)], (1+0j)), ([(True, 1), (False, 0)], (-1+0j))]
[([(True, 2), (False, 3)], (2+0j)), ([(True, 3), (False, 2)], (-2+0j))]
>>> groups = op.split_out_groups(group_indices=[1])
>>> for g in groups:
...     print(list(sorted(g.iter_terms())))
[([(True, 2), (False, 3)], (2+0j)), ([(True, 3), (False, 2)], (-2+0j))]
```

* **Parameters:**
  **group_indices** – the group indices for which to build operators, in the desired output
  order. When omitted, every group is built, in index order.
* **Returns:**
  An optional vector of one new operator for each requested group index.

#### *classmethod* zero()

Constructs the additive identity operator.

Adding the operator that is constructed by this method to another one has no effect.

```pycon
>>> from qiskit_bosons.operators import BosonOperator
>>> op = BosonOperator.from_dict({(): 2.0})
>>> zero = BosonOperator.zero()
>>> op + zero == op
True
```

### Protocol Methods

#### *static* \_anti_commutator_(op_a, op_b)

#### *static* \_commutator_(op_a, op_b)

#### *static* \_double_commutator_(op_a, op_b, op_c, sign)
