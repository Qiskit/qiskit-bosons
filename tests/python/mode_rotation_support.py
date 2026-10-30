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

"""Mode rotation syntheses written outside the library.

The tests use them to show that the framework does not depend on the one decomposition and the one
pattern that the library ships: a different boson-to-boson decomposition (Reck) and a boson-to-qubit
plugin that maps the generator of a whole :class:`.ModeRotation` directly.
"""

from __future__ import annotations

import numpy as np
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.dagcircuit import DAGCircuit, DAGOpNode
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, ModeRotation
from qiskit_bosons.transpiler import B2QLayout
from qiskit_bosons.transpiler.passes.utils import map_node_single_register


def reck_factors(unitary: np.ndarray) -> list[tuple[tuple[int, ...], np.ndarray]]:
    """Factorizes a unitary into adjacent-mode blocks in the triangular arrangement of Reck et al.

    The columns of ``unitary`` are eliminated one after the other from the bottom up with two-mode
    rotations, which leaves a diagonal. Returns the ``(modes, block)`` pairs in the order in which
    the rotations are to be applied: the diagonal phases first, then the rotations.
    """
    num_modes = len(unitary)
    remaining = unitary.astype(complex).copy()
    steps = []
    for column in range(num_modes - 1):
        for row in range(num_modes - 1, column, -1):
            upper, lower = remaining[row - 1, column], remaining[row, column]
            if abs(lower) < 1e-14:
                continue
            norm = np.hypot(abs(upper), abs(lower))
            elimination = np.array([[np.conj(upper), np.conj(lower)], [-lower, upper]]) / norm
            remaining[[row - 1, row], :] = elimination @ remaining[[row - 1, row], :]
            steps.append(((row - 1, row), elimination))
    assert np.allclose(remaining, np.diag(np.diag(remaining)), atol=1e-10)

    factors: list[tuple[tuple[int, ...], np.ndarray]] = [
        ((mode,), np.array([[remaining[mode, mode]]])) for mode in range(num_modes)
    ]
    factors += [(modes, block.conj().T) for modes, block in reversed(steps)]
    return factors


class BosonicReckDecomposition:
    """The triangular mesh of Reck et al., through the same plugin protocol as the Clements mesh.

    It needs more depth than the rectangular arrangement and is only used to show that the interface
    is not tied to one decomposition.
    """

    def synthesize(self, gate: ModeRotation) -> BosonicCircuit:
        """Synthesizes the gate into atomic, single-group Evolution gates."""
        definition = BosonicCircuit(gate.num_modes)
        for modes, block in reck_factors(gate.rotation_unitary):
            generator = ModeRotation(block).to_operator()
            if len(generator) == 0:
                continue
            generator.groups = [0] * len(generator)
            definition.append(
                Evolution(len(modes), generator, time=1.0, atomic=True),
                [definition.modes[mode] for mode in modes],
            )
        return definition


class DirectModeRotationSynthesis:
    """A boson-to-qubit plugin mapping the generator of the whole rotation, with no decomposition.

    This is the whole plugin. The generator is the Hermitian operator of
    :meth:`.ModeRotation.to_operator`, so the rotation is a single evolution at unit time. The
    truncation is applied once, to the generator of the complete rotation.
    """

    def __init__(self, mapper_fn, product_formula=None) -> None:
        self.mapper_fn = mapper_fn
        self.product_formula = product_formula

    def run(self, in_node: DAGOpNode, out_dag: DAGCircuit, *, b2q_layout: B2QLayout) -> None:
        """Runs this transpilation plugin."""
        modes, qreg = map_node_single_register(in_node, b2q_layout)
        generator = in_node.op.to_operator().relabel_modes(list(modes))
        out_dag.apply_operation_back(
            PauliEvolutionGate(
                self.mapper_fn(generator, len(qreg)), time=1.0, synthesis=self.product_formula
            ),
            tuple(qreg),
        )
