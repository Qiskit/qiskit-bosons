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

"""The Clements decomposition of a mode rotation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from qiskit_bosons.linalg import givens_decomposition

if TYPE_CHECKING:
    from qiskit_bosons.circuit import BosonicCircuit

    from ..mode_rotation import ModeRotation


class BosonicClementsDecomposition:
    r"""The decomposition of a mode rotation into a mesh of two-mode rotations and phases.

    The unitary :math:`U` of the :class:`.ModeRotation` is factorized with
    :func:`.givens_decomposition` into at most :math:`M(M-1)/2` rotations between *adjacent* modes
    plus :math:`M` single-mode phases, the rectangular mesh of Clements et al. [1]_. Every factor
    :math:`F` is itself a mode rotation, whose many-body unitary is :math:`e^{-\mathrm{i} H_F}` for
    the Hermitian generator :math:`H_F` of :meth:`.ModeRotation.to_operator`. It is therefore emitted
    as an :class:`.Evolution` gate at unit time on the one or two modes it acts on.

    The emitted gates are :attr:`~.Evolution.atomic`, and the operator of each holds a single group.
    The terms of a two-mode generator do not commute, so splitting it up any further would change the
    factor: the boson-to-qubit stage maps and synthesizes it whole.

    Because :math:`U \mapsto \exp(\sum_{ij} \log(U)_{ij} b^\dagger_i b_j)` is a group homomorphism,
    the product of the factors equals the full rotation exactly on the untruncated Fock space. The
    mesh itself therefore adds no error on its own. What it buys is locality: every factor acts on
    the modes of its own, instead of one evolution reaching all of them.

    .. warning::
       A truncation does not preserve the homomorphism, because the truncated algebra is not closed.
       The product of the truncated factors then differs from the truncation of the full rotation
       whenever one of the factors takes a mode past the cutoff. A cutoff of at least the total boson
       number keeps every factor within it, as each factor conserves that number. A smaller cutoff
       can still be enough, but is then a property of the state and the rotation that has to be
       converged. Increase the cutoff until the result stops changing.

    .. note::
       The emitted :class:`.Evolution` gates are only exact once the boson-to-qubit stage synthesizes
       them well. The terms of a two-mode generator do not commute, so a single-step product formula
       leaves a sizeable error in each factor. Converge it with the ``product_formula`` of the
       boson-to-qubit plugin, for example :class:`~qiskit.synthesis.LieTrotter` with more repetitions.

    .. [1] W. R. Clements et al., Optimal design for universal multiport interferometers,
       Optica 3, 1460 (2016), `arXiv:1603.08788 <https://arxiv.org/abs/1603.08788>`_.
    """

    def __init__(self, *, tol: float = 1e-12) -> None:
        """Initializing an instance of this synthesis method can be done with the arguments below.

        Args:
            tol: factors and entries below this magnitude are dropped. This is applied to the
                entries of the generator of a factor, see :meth:`.ModeRotation.to_operator`, and to
                the angle of a single-mode phase.
        """
        self.tol = tol
        """The magnitude below which entries of a generator and angles of a phase are dropped."""

    def synthesize(self, gate: ModeRotation) -> BosonicCircuit:
        """Synthesizes the gate into the mesh of atomic :class:`.Evolution` gates.

        See the class documentation for the decomposition this implements.

        Args:
            gate: the gate to synthesize.

        Returns:
            A :class:`.BosonicCircuit` holding one :class:`.Evolution` gate per two-mode rotation
            and per non-trivial phase, in the order in which they are applied. The rotations come
            first, then the phases.
        """
        # deferred imports: `qiskit_bosons.circuit` imports this subpackage's siblings at module load
        from qiskit_bosons.circuit import BosonicCircuit

        from ..evolution import Evolution
        from ..mode_rotation import ModeRotation

        rotations, phases = givens_decomposition(gate.rotation_unitary)
        definition = BosonicCircuit(gate.num_modes)

        def append_factor(block: np.ndarray, modes: tuple[int, ...]) -> None:
            # the factor acts on its own modes only, so its generator is built on local modes 0, 1
            generator = ModeRotation(block).to_operator(tol=self.tol)
            if len(generator) == 0:
                return
            # one group: the factor is a single unit and must not be split up any further
            generator.groups = [0] * len(generator)
            definition.append(
                Evolution(len(modes), generator, time=1.0, atomic=True),
                [definition.modes[mode] for mode in modes],
            )

        # U = D conj(G_N) ... conj(G_1) and U -> exp(sum log(U)_ij b†_i b_j) is a group homomorphism,
        # so the circuit applies conj(G_1) first and the phases last.
        for c, s, i, j in rotations:
            append_factor(np.array([[c, np.conj(s)], [-s, c]], dtype=complex), (i, j))
        for mode, phase in enumerate(phases):
            angle = float(np.angle(phase))
            if abs(angle) > self.tol:
                append_factor(np.array([[np.exp(1j * angle)]]), (mode,))

        return definition
