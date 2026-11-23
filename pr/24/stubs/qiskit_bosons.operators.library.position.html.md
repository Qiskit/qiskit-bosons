# position

### position(mode)

Builds the dimensionless position (quadrature) operator of a single mode.

$$
\hat{q} = \frac{1}{\sqrt{2}}\left(b^\dagger + b\right) \, .
$$

```pycon
>>> from qiskit_bosons.operators.library import position
>>> position(0).is_hermitian()
True
```

* **Parameters:**
  **mode** ([*int*](https://docs.python.org/3/builtins/functions.html#int)) – the mode index.
* **Returns:**
  The position operator as a [`BosonOperator`](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator).
* **Return type:**
  [*BosonOperator*](qiskit_bosons.operators.BosonOperator.html.md#qiskit_bosons.operators.BosonOperator)
