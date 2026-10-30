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

"""Mode initialization gate."""

from __future__ import annotations

from collections.abc import Sequence

from .. import BosonicGate


class InitializeModes(BosonicGate):
    r"""Prepares a bosonic Fock state.

    This gate declares the occupation of the modes it is placed on, i.e. the Fock state

    .. math::

        \lvert n_0, n_1, \ldots, n_{M-1} \rangle
        = \prod_m \frac{(b^\dagger_m)^{n_m}}{\sqrt{n_m!}} \lvert \text{vac} \rangle \, ,

    where :math:`n_m` is ``occupation[m]``. Synthesis emits the gates that set each mode's qubits to
    the pattern encoding its occupation level.

    .. note::
       The occupation entries are **non-negative integers**, not booleans: a single mode may hold
       arbitrarily many bosons.

    .. note::
       This gate does not itself fix a truncation. The occupation it declares must of course be
       representable by whatever cutoff the boson-to-qubit encoding uses (an occupation of 3
       requires ``max_occupation >= 3``) which the synthesis plugin checks at transpile time, since
       that is where the cutoff is known. :attr:`max_occupation` reports the largest level this gate
       declares.

    .. caution::
       This is an early development prototype. Beware of changes to its interface without warning
       during the pre-release development of this package.
    """

    def __init__(self, occupation: Sequence[int]) -> None:
        """Initializing an instance of this gate can be done with the arguments listed below.

        Args:
            occupation: the boson number of each mode in the :class:`.BosonicRegister` being
                initialized by this gate.

        Raises:
            ValueError: if any occupation entry is negative.
        """
        occ = [int(n) for n in occupation]
        if any(n < 0 for n in occ):
            raise ValueError(f"occupation must be non-negative, got {occ!r}.")

        self.occupation: tuple[int, ...] = tuple(occ)
        """The boson number of each mode in the :class:`~.BosonicRegister` initialized by this gate."""

        super().__init__("InitializeModes", len(self.occupation), [])

    @property
    def max_occupation(self) -> int:
        """The largest occupation level this gate declares.

        A boson-to-qubit encoding must retain at least this many levels to represent the declared
        Fock state; see :func:`~qiskit_bosons.mappers.library.occupation_encoding`.
        """
        return max(self.occupation, default=0)

    @classmethod
    def vacuum(cls, num_modes: int) -> InitializeModes:
        r"""Builds the gate for the vacuum state :math:`\lvert 0, \ldots, 0 \rangle`.

        Args:
            num_modes: the number of bosonic modes to initialize.

        Returns:
            An :class:`InitializeModes` gate declaring zero bosons in every mode.
        """
        return cls([0] * num_modes)

    def __eq__(self, other: object) -> bool:
        """Compares two gates by their declared occupation."""
        if not isinstance(other, InitializeModes):
            return NotImplemented
        return self.occupation == other.occupation

    def __hash__(self) -> int:
        """Hashes this gate by its declared occupation."""
        return hash((type(self).__name__, self.occupation))
