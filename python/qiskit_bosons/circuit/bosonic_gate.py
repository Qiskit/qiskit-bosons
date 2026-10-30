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

"""BosonicGate."""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING, Self, cast

from qiskit.circuit import Gate

if TYPE_CHECKING:
    from .library.synthesis import B2BSynthesisPlugin


class BosonicGate(Gate):
    """The base class for all bosonic gates.

    To ensure consistency only subclasses of this gate class can be added to instances of
    :class:`.BosonicCircuit`. As such, this class (mostly) serves as a type (for the time being).

    .. rubric:: Breaking a gate down in bosonic space

    Every bosonic gate has a :attr:`synthesis` and an :attr:`atomic` flag, which the boson-to-boson
    synthesis framework relies on (see :mod:`~qiskit_bosons.circuit.library.synthesis`). The
    :attr:`synthesis` is the method that rewrites the gate into a circuit of smaller bosonic gates,
    and it is ``None`` for a gate that does not break down. A gate that is :attr:`atomic` is never
    broken down, however it is asked to. Both are chosen when the gate is built and are read-only
    afterwards.

    .. caution::
       Since this is a subclass of :class:`~qiskit.circuit.Gate` the documentation of its methods
       may refer to `qubits`. Those references should be interpreted as referring to `bosonic modes`
       in the context of instances of this subclass.

       It may also happen that some of the inherited methods may not always make sense because of
       this `re-interpretation`. You have been warned.

    .. note::
       A bosonic gate carries **no occupation cutoff**: it acts on bosonic modes, whose Fock space is
       infinite-dimensional. A gate becomes a finite-dimensional operator only once a boson-to-qubit
       encoding fixes a truncation, which happens during transpilation. Simulating a bosonic circuit
       therefore means transpiling it to qubits first and then simulating the resulting qubit
       circuit with any qubit-level simulator; there is deliberately no mode-space state-vector
       protocol on this class.
    """

    def __init__(  # noqa: D107
        self,
        name: str,
        num_modes: int,
        /,
        params: list | None = None,
        *,
        label: str | None = None,
        synthesis: B2BSynthesisPlugin | None = None,
        atomic: bool = False,
    ) -> None:
        if params is None:
            params = []

        # both read-only (see their properties): the definition is cached the first time it is built,
        # so a later re-assignment of either would be silently ignored
        self._synthesis = synthesis
        self._atomic = atomic

        super().__init__(name, num_modes, params, label)

    @property
    def synthesis(self) -> B2BSynthesisPlugin | None:
        """The boson-to-boson synthesis method with which this gate is broken down, if any.

        This is the customizable part of breaking a gate down: any object implementing the
        :class:`.B2BSynthesisPlugin` protocol can be given, and the :class:`.B2BSynthesis` transpiler
        pass selects one per gate type for a whole circuit. ``None`` marks a gate that does not break
        down in bosonic space.

        What the circuit that the method produces is equivalent to is specified by the gate, and
        documented there. Where the gate's decomposition is that circuit,
        :meth:`~.BosonicCircuit.decompose` yields it; a gate may deviate from this and document so.

        This is read-only: the circuit that :meth:`~.BosonicCircuit.decompose` yields is built once
        and then cached, so a synthesis method assigned after the fact could not take effect. Pass it
        to the constructor instead, or use :meth:`with_synthesis`.
        """
        return self._synthesis

    @property
    def atomic(self) -> bool:
        r"""Whether this gate is a factor that must not be broken down any further.

        A synthesis marks the factors it emits as atomic, because they are the result of the split it
        chose to perform: breaking them down again would discard that choice, and would split a
        Hermitian group into individually non-Hermitian terms whose exponentials are not unitary. A
        gate that is atomic for the stronger reason that there is nothing left to split is marked the
        same way.

        An atomic gate has no definition at all, so :meth:`~.BosonicCircuit.decompose` leaves it in
        place rather than expanding it, and the :class:`.B2BSynthesis` pass leaves it alone. It still
        maps and synthesizes normally at the boson-to-qubit level, which is where a factor is turned
        into actual operations.

        This is read-only, for the same reason as :attr:`synthesis`.
        """
        return self._atomic

    def with_synthesis(self, synthesis: B2BSynthesisPlugin | None) -> Self:
        """Returns a copy of this gate that is broken down with another synthesis method.

        The gate itself is left as it is, because a gate instance can be shared between circuits.

        Args:
            synthesis: the boson-to-boson synthesis method of the copy.

        Returns:
            A shallow copy of this gate whose :attr:`synthesis` is ``synthesis``.
        """
        gate = copy.copy(self)
        gate._synthesis = synthesis
        # the cached definition belongs to the method it was built with
        gate._definition = None
        return gate

    def _define(self) -> None:
        # A gate that is atomic, or that has no synthesis, has no bosonic definition: it is left
        # untouched by `decompose()` and lowered directly by the boson-to-qubit stage. Note that this
        # must leave `_definition` as `None` rather than assigning an empty circuit -- an empty
        # definition would make `decompose()` *drop* the gate instead of keeping it.
        if self._atomic or self._synthesis is None:
            return

        self._definition = self._synthesis.synthesize(self)._inner

    @property
    def num_modes(self) -> int:
        """The number of bosonic modes that this gate acts upon."""
        return cast(int, self._num_qubits)
