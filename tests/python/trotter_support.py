# This code is a Qiskit project.
#
# (C) Copyright IBM 2026.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at https://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Shared helpers for tests of the bosonic Trotterization framework.

Everything here is laptop-scale: three modes, three occupation levels (``d = 3``).

The numerical checks isolate the *outer* Trotter error, the one between groups that do not commute.
They do so by simulating the transpiled circuit **without** decomposing its
:class:`~qiskit.circuit.library.PauliEvolutionGate`, which :class:`~qiskit.quantum_info.Statevector`
applies as the exact exponential of its group. What is left is exactly the product of the group
exponentials that the bosonic synthesis chose, so its distance from the exact evolution measures that
synthesis and nothing else.
"""

from __future__ import annotations

import numpy as np
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply

from .utils import preset_without_qubit_stage

NUM_SITES = 3
MAX_OCCUPATION = 2  # d = 3 levels per mode
ENCODING = "unary_projected"
OCCUPATION = [1, 0, 1]  # one boson on each end of the chain


def bose_hubbard(
    num_sites: int = NUM_SITES, tunneling: float = 1.0, interaction: float = 2.0
) -> BosonOperator:
    """The normal-ordered Bose-Hubbard chain, grouped so that every group is Hermitian.

    One group per bond (its two hopping terms are each other's adjoint, so the group is Hermitian
    although a lone ``b†_i b_j`` is not) plus one group for the mutually commuting on-site terms.
    """
    edges = [(site, site + 1) for site in range(num_sites - 1)]

    terms: dict = {}
    for left, right in edges:
        terms[(cre(left), ann(right))] = -tunneling
        terms[(cre(right), ann(left))] = -tunneling
    for site in range(num_sites):
        terms[(cre(site), cre(site), ann(site), ann(site))] = 0.5 * interaction

    operator = BosonOperator.from_dict(terms)

    groups = []
    for term, _ in operator.iter_terms():
        modes = {mode for _, mode in term}
        groups.append(len(edges) if len(modes) == 1 else edges.index(tuple(sorted(modes))))
    operator.groups = groups
    return operator


def initial_state(circuit_width: int | None = None) -> Statevector:
    """The unary-encoded product Fock state ``OCCUPATION``, as a qubit state vector."""
    per_mode = MAX_OCCUPATION + 1  # unary: one qubit per level
    index = 0
    for mode, level in enumerate(OCCUPATION):
        index |= (1 << level) << (mode * per_mode)
    width = circuit_width or NUM_SITES * per_mode
    return Statevector.from_int(index, 2**width)


def exact_state(operator: BosonOperator, time: float) -> np.ndarray:
    """``exp(-i t H)`` applied to the initial state, with ``H`` mapped at the same truncation."""
    observable = occupation_encoding(operator, LevelSpec(MAX_OCCUPATION, ENCODING)).simplify()
    matrix = csc_matrix(SparsePauliOp.from_sparse_observable(observable).to_matrix(sparse=True))
    return np.asarray(expm_multiply(-1j * time * matrix, initial_state().data))


def circuit_state(circuit: BosonicCircuit) -> np.ndarray:
    """Transpiles a bosonic circuit under the unary encoding and simulates it.

    The ``PauliEvolutionGate`` of every emitted factor is deliberately left un-decomposed, so each is
    applied as the exact exponential of its group (see the module docstring). The circuit starts from
    the all-zero qubit state because it carries its own ``InitializeModes``.
    """
    lowered = preset_without_qubit_stage(LevelSpec(MAX_OCCUPATION, ENCODING)).run(circuit)
    return np.asarray(Statevector.from_int(0, 2**lowered.num_qubits).evolve(lowered).data)


def evolution_circuit(operator: BosonOperator, time: float, **kwargs) -> BosonicCircuit:
    """A circuit holding the initial occupation followed by one ``Evolution`` gate."""
    circuit = BosonicCircuit(NUM_SITES)
    circuit.append(InitializeModes(OCCUPATION), circuit.modes)
    circuit.append(Evolution(NUM_SITES, operator, time=time, **kwargs), circuit.modes)
    return circuit


def distance(left: np.ndarray, right: np.ndarray) -> float:
    """The largest amplitude difference between two state vectors."""
    return float(np.max(np.abs(left - right)))
