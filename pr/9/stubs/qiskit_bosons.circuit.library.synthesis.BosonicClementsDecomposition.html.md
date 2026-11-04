# BosonicClementsDecomposition

### *class* BosonicClementsDecomposition(, tol=1e-12)

The decomposition of a mode rotation into a mesh of two-mode rotations and phases.

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

The unitary $U$ of the [`ModeRotation`](qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation) is factorized with
[`givens_decomposition()`](qiskit_bosons.linalg.givens_decomposition.html.md#qiskit_bosons.linalg.givens_decomposition) into at most $M(M-1)/2$ rotations between *adjacent* modes
plus $M$ single-mode phases, the rectangular mesh of Clements et al. <sup>[1](#id2)</sup>. Every factor
$F$ is itself a mode rotation, whose many-body unitary is $e^{-\mathrm{i} H_F}$ for
the Hermitian generator $H_F$ of [`ModeRotation.to_operator()`](qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation.to_operator). It is therefore emitted
as an [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate at unit time on the one or two modes it acts on.

The emitted gates are `atomic`, and the operator of each holds a single group.
The terms of a two-mode generator do not commute, so splitting it up any further would change the
factor: the boson-to-qubit stage maps and synthesizes it whole.

Because $U \mapsto \exp(\sum_{ij} \log(U)_{ij} b^\dagger_i b_j)$ is a group homomorphism,
the product of the factors equals the full rotation exactly on the untruncated Fock space. The
mesh itself therefore adds no error on its own. What it buys is locality: every factor acts on
the modes of its own, instead of one evolution reaching all of them.

#### WARNING
A truncation does not preserve the homomorphism, because the truncated algebra is not closed.
The product of the truncated factors then differs from the truncation of the full rotation
whenever one of the factors takes a mode past the cutoff. A cutoff of at least the total boson
number keeps every factor within it, as each factor conserves that number. A smaller cutoff
can still be enough, but is then a property of the state and the rotation that has to be
converged. Increase the cutoff until the result stops changing.

#### NOTE
The emitted [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates are only exact once the boson-to-qubit stage synthesizes
them well. The terms of a two-mode generator do not commute, so a single-step product formula
leaves a sizeable error in each factor. Converge it with the `product_formula` of the
boson-to-qubit plugin, for example [`LieTrotter`](https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.synthesis.LieTrotter) with more repetitions.

* <a id='id2'>**[1]**</a> W. R. Clements et al., Optimal design for universal multiport interferometers, Optica 3, 1460 (2016), [arXiv:1603.08788](https://arxiv.org/abs/1603.08788).

Initializing an instance of this synthesis method can be done with the arguments below.

* **Parameters:**
  **tol** ([*float*](https://docs.python.org/3/builtins/functions.html#float)) – factors and entries below this magnitude are dropped. This is applied to the
  entries of the generator of a factor, see [`ModeRotation.to_operator()`](qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation.to_operator), and to
  the angle of a single-mode phase.

### Attributes

#### tol

The magnitude below which entries of a generator and angles of a phase are dropped.

### Methods

#### synthesize(gate)

Synthesizes the gate into the mesh of atomic [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gates.

See the class documentation for the decomposition this implements.

* **Parameters:**
  **gate** ([*ModeRotation*](qiskit_bosons.circuit.library.ModeRotation.html.md#qiskit_bosons.circuit.library.ModeRotation)) – the gate to synthesize.
* **Returns:**
  A [`BosonicCircuit`](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit) holding one [`Evolution`](qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution) gate per two-mode rotation
  and per non-trivial phase, in the order in which they are applied. The rotations come
  first, then the phases.
* **Return type:**
  [BosonicCircuit](qiskit_bosons.circuit.BosonicCircuit.html.md#qiskit_bosons.circuit.BosonicCircuit)
