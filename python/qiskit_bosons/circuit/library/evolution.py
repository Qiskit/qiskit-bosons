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

"""Hamiltonian evolution gate."""

from __future__ import annotations

from typing import TYPE_CHECKING

from qiskit_bosons.operators import OperatorTrait

from .. import BosonicGate

if TYPE_CHECKING:
    from .synthesis import B2BSynthesisPlugin


class Evolution(BosonicGate):
    r"""Implements the time evolution of a bosonic operator.

    Given a bosonic ``operator`` :math:`H` and an evolution ``time`` :math:`t`, this gate
    implements the unitary

    .. math::

        U = e^{-i t H}.

    For :math:`U` to be unitary, :math:`H` must be Hermitian. This is the caller's responsibility
    and is not verified; :meth:`~qiskit_bosons.operators.BosonOperator.is_hermitian` can check it.

    .. note::
       Breaking this gate down in bosonic space (into a circuit of smaller :class:`.Evolution`
       gates) is **optional**, and it is what the :attr:`~.BosonicGate.synthesis` attribute
       governs. The gate can equally be handed straight to the boson-to-qubit stage, which maps and
       synthesizes it whole, however many terms its ``operator`` holds. See
       :mod:`~qiskit_bosons.circuit.library.synthesis`.

       When it `is` decomposed (by :meth:`~.BosonicCircuit.decompose`, for instance)
       :attr:`~.BosonicGate.synthesis` decides how. It defaults to :class:`.BosonicLieTrotter`,
       which splits the evolution group-by-group when the ``operator`` has
       :attr:`~qiskit_bosons.operators.BosonOperator.groups` assigned (see
       :ref:`grouping_explanation`), and term-by-term otherwise (exact only when the individual
       factors mutually commute). Another method can be chosen, for a whole circuit,
       with the :class:`.B2BSynthesis` transpiler pass; :class:`.BosonicQDrift`, for example,
       replaces the evolution with a randomized sample of its terms. The gates it emits are :attr:`~.BosonicGate.atomic`, so
       :attr:`~.BosonicGate.synthesis` no longer applies to them.

    .. note::
       This gate carries **no occupation cutoff**. It is defined purely in terms of the bosonic
       operator, whose Fock space is infinite-dimensional. The truncation is introduced later, by the
       boson-to-qubit encoding chosen during transpilation -- which is precisely what lets one
       ansatz be evaluated at several truncation levels without rebuilding it. See
       :func:`~qiskit_bosons.mappers.library.occupation_encoding`.
    """

    def __init__(
        self,
        num_modes: int,
        operator: OperatorTrait,
        time: float = 1.0,
        *,
        synthesis: B2BSynthesisPlugin | None = None,
        atomic: bool = False,
    ) -> None:
        r"""Initializing an instance of this gate can be done with the arguments listed below.

        Args:
            num_modes: the number of bosonic modes on which this gate acts.
            operator: the Hermitian operator :math:`H` under which to time evolve the acted-upon
                bosonic modes.
            time: the evolution time :math:`t` entering the exponent of :math:`e^{-i t H}`. A
                negative value evolves backwards in time.
            synthesis: the boson-to-boson synthesis method with which to decompose this gate, when
                it gets decomposed at all. If ``None`` (the default), a
                :class:`.BosonicLieTrotter` instance is used.
            atomic: whether this gate is a fully decomposed factor which must not be broken down any
                further in bosonic space. See :attr:`.BosonicGate.atomic`.
        """
        from .synthesis import BosonicLieTrotter

        self.operator = operator
        """The operator under which to time evolve the acted-upon bosonic modes."""

        super().__init__(
            "Evolution",
            num_modes,
            [time],
            label=f"evolve({' '.join(str(operator).split())})",
            synthesis=synthesis or BosonicLieTrotter(),
            atomic=atomic,
        )

    def inverse(self, annotated: bool = False) -> Evolution:
        r"""Returns the inverse of this gate, :math:`e^{+i t H}`.

        Args:
            annotated: ignored. The inverse of this gate is another :class:`.Evolution`, so it never
                needs to be an :class:`~qiskit.circuit.AnnotatedOperation`.

        Returns:
            An :class:`.Evolution` gate evolving the same operator for the negated time.
        """
        # `Instruction.inverse` would recurse through the definition and return a plain `Gate`,
        # losing the operator and this gate's synthesis method; negating the time is both cheaper and
        # exact.
        return Evolution(
            self.num_modes,
            self.operator,
            time=-self.params[0],
            synthesis=self.synthesis,
            atomic=self.atomic,
        )
