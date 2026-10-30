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

"""BosonicCircuit."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, cast

from qiskit.circuit import Instruction, QuantumCircuit, QuantumRegister

from . import BosonicMode, BosonicSpecifier
from .bosonic_gate import BosonicGate

if TYPE_CHECKING:
    from . import BosonicRegister


class BosonicCircuit:
    """A wrapper around :class:`~qiskit.circuit.QuantumCircuit` for expressing bosonic circuits.

    This class maintains a reduced API compared to the full API of the underlying
    :class:`~qiskit.circuit.QuantumCircuit`. This is done to avoid exposing (among other methods)
    the ability to apply qubit-based gates onto a bosonic circuit, which would not be a
    well-defined operation in the general case.

    :meth:`barrier` is re-exposed because it is an exception to that rule: a barrier acts on modes
    without interpreting them and carries no unitary content, so it is as well-defined on a bosonic
    register as it is on a qubit one.
    """

    def __init__(self, num_modes: int) -> None:
        """Initializing a circuit instance can be done with the arguments listed below.

        Args:
            num_modes: the number of bosonic modes on which this circuit acts.
        """
        self.register: BosonicRegister = QuantumRegister(num_modes, "b")
        """The inner circuit's :type:`~qiskit_bosons.circuit.BosonicRegister`."""
        self._inner = QuantumCircuit(self.register)

    @property
    def metadata(self) -> dict:
        """Re-exposes :external:attr:`~qiskit.circuit.QuantumCircuit.metadata`."""
        return cast(dict, self._inner.metadata)

    @metadata.setter
    def metadata(self, metadata: dict) -> None:
        self._inner.metadata = metadata

    @property
    def modes(self) -> list[BosonicMode]:
        """The bosonic mode ``bits`` that this circuit acts upon."""
        return cast(list[BosonicMode], self._inner.qubits)

    def append(
        self,
        gate: BosonicGate,
        bargs: BosonicSpecifier,
        cargs: None = None,
        *,
        copy: bool = True,
    ) -> None:
        """Appends a :class:`.BosonicGate` to this circuit.

        Args:
            gate: the bosonic gate to apply.
            bargs: the bosonic modes on which this gate acts.
            cargs: the classical bits on which this gate acts.

              .. warning::
                 No gates of this kind are currently supported.

            copy: forwarded to :meth:`~qiskit.circuit.QuantumCircuit.append`.

        Raises:
            ValueError: if the provided ``gate`` is not an instance of :class:`.BosonicGate`.
        """
        if not isinstance(gate, BosonicGate):
            raise ValueError(f"Unsupported instruction type: {type(gate)}")
        self._inner.append(gate, bargs, cargs, copy=copy)

    def barrier(self, *bargs: BosonicSpecifier, label: str | None = None) -> None:
        """Re-exposes :external:meth:`~qiskit.circuit.QuantumCircuit.barrier`.

        A barrier carries no unitary effect. It marks a synchronization point that optimization
        passes will not fuse gates across.

        The boson-to-qubit synthesis stage maps a barrier onto every qubit of each register it
        touched, so it constrains the subsequent qubit-level transpilation too. Because a
        :type:`~qiskit_bosons.transpiler.B2QLayout` relates whole registers, whose two sides differ in
        size (an encoding needs several qubits per mode), a barrier is deliberately not narrowed to
        the qubits encoding the modes it was placed on.

        Args:
            bargs: the bosonic modes on which to place the barrier. When omitted, the barrier spans
                every mode of this circuit.
            label: the string label of the barrier.
        """
        self._inner.barrier(*bargs, label=label)

    def compose(
        self,
        other: BosonicCircuit,
        bargs: BosonicSpecifier | None = None,
        *,
        front: bool = False,
        inplace: bool = False,
    ) -> BosonicCircuit | None:
        """Composes another :class:`.BosonicCircuit` onto this one.

        The instructions of ``other`` are inlined into this circuit, rather than wrapped into a single
        gate. This mirrors :external:meth:`~qiskit.circuit.QuantumCircuit.compose` with its default
        ``wrap=False``. The ``wrap=True`` counterpart is deliberately not offered: a dedicated
        encapsulation mechanism is planned instead, so adding a second way to box up a block of
        instructions here would leave two patterns solving the same problem.

        Args:
            other: the circuit to compose onto this one. It may act on fewer modes than this circuit,
                in which case ``bargs`` selects the modes it is placed on.
            bargs: the bosonic modes of this circuit onto which ``other`` is placed. Defaults to all
                of this circuit's modes, in order.
            front: whether to compose ``other`` before rather than after this circuit's instructions.
            inplace: whether to modify this circuit rather than return a new one.

        Returns:
            The composed circuit, or ``None`` if ``inplace`` is set.

        Raises:
            ValueError: if ``other`` is not an instance of :class:`.BosonicCircuit`. Only bosonic
                circuits are accepted because the plain
                :external:meth:`~qiskit.circuit.QuantumCircuit.compose` would happily inline
                qubit-based instructions, bypassing the type guard that :meth:`append` enforces.
            CircuitError: if ``other`` acts on more modes than this circuit, or ``bargs`` does not
                match its width.
        """
        if not isinstance(other, BosonicCircuit):
            raise ValueError(f"Unsupported circuit type: {type(other)}")

        # Every `BosonicCircuit` is populated through the guarded `append`, so the instructions of
        # `other` are already known to be `BosonicGate`s (or barriers); the operand type check above
        # is all that is needed and keeps this O(1) rather than re-walking the instructions.
        if inplace:
            self._inner.compose(other._inner, bargs, front=front, inplace=True)
            return None

        inner_composed = self._inner.compose(other._inner, bargs, front=front, inplace=False)
        out = BosonicCircuit(len(self.register))
        # `QuantumCircuit.compose` returns a circuit carrying a *new* register object, so re-sync it
        # with ours (as `decompose` does) to keep the mode identities valid.
        out.register = self.register
        out._inner = inner_composed
        return out

    def count_ops(self) -> OrderedDict[str, int]:
        """Re-exposes :external:meth:`~qiskit.circuit.QuantumCircuit.count_ops`."""
        return cast(OrderedDict[str, int], self._inner.count_ops())

    def decompose(
        self,
        gates_to_decompose: (
            str | type[Instruction] | Sequence[str | type[Instruction]] | None
        ) = None,
        reps: int = 1,
    ) -> BosonicCircuit:
        """Re-exposes :external:meth:`~qiskit.circuit.QuantumCircuit.decompose`."""
        inner_decomposed = self._inner.decompose(gates_to_decompose=gates_to_decompose, reps=reps)
        out = BosonicCircuit(len(self.register))
        out.register = self.register
        out._inner = inner_decomposed
        return out

    def draw(self, *args, **kwargs) -> Any:
        """Directly exposes the inner circuit's :meth:`~qiskit.circuit.QuantumCircuit.draw` method."""
        return self._inner.draw(*args, **kwargs)

    def repeat(self, reps: int, *, insert_barriers: bool = False) -> BosonicCircuit:
        """Repeats this circuit ``reps`` times.

        The repetitions are `flattened`: the returned circuit holds this circuit's instructions
        ``reps`` times in sequence, rather than one gate wrapping the circuit as
        :external:meth:`~qiskit.circuit.QuantumCircuit.repeat` produces. As a consequence
        :meth:`count_ops` scales with ``reps``, and :meth:`decompose` decomposes the repeated body
        rather than peeling off a repetition.

        Flattening also leaves the gates that meet at a repetition boundary adjacent, so optimization
        passes can act across it. Set ``insert_barriers`` to keep the repetitions apart instead, which
        places a :meth:`barrier` that those passes do not fuse across.

        .. note::
           This is unrelated to the ``reps`` argument of :class:`.BosonicSuzukiTrotter`, which
           subdivides a `single` :class:`.Evolution` into more, shorter factors. Repeating a circuit
           holding one :class:`.Evolution` merely multiplies its total evolution time. To converge a
           grouped evolution, give the gate a synthesis with several ``reps`` and a ``time`` of
           ``total / steps``, or repeat a *single* step ``steps`` times.

        .. seealso::
           :ref:`bosonic_circuit_explanation` works through a repeated Trotter step, including the
           boundary fusion and how to prevent it.

        Args:
            reps: how often this circuit should be repeated. Repeating it zero times yields an empty
                circuit, which is the identity.
            insert_barriers: whether to place a :meth:`barrier` between consecutive repetitions. No
                barrier is placed after the final repetition, so the returned circuit ends on the same
                instruction it would without this argument.

        Returns:
            A new circuit containing ``reps`` repetitions of this circuit.

        Raises:
            ValueError: if ``reps`` is negative.
        """
        if reps < 0:
            raise ValueError(f"The number of repetitions must be non-negative, but got {reps}.")

        out = BosonicCircuit(len(self.register))
        out.register = self.register
        # `copy_empty_like` keeps the register, mode objects, name and metadata of this circuit, so the
        # repeated circuit differs from it only in the instructions composed in below.
        out._inner = self._inner.copy_empty_like()

        # NOTE: the repetitions share this circuit's gate instances rather than copying them per
        # repetition. `QuantumCircuit.compose` only copies operations carrying a
        # `ParameterExpression`, and no bosonic gate is parameterized in that sense. Sharing is safe
        # here: the gates expose their configuration read-only (e.g. `Evolution.synthesis`), the
        # transpiler substitutes nodes instead of mutating operations in place, and a gate's cached
        # definition is a deterministic function of that configuration.
        for rep in range(reps):
            out.compose(self, inplace=True)
            if insert_barriers and rep != reps - 1:
                out.barrier()

        return out
