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

"""An occupation-encoding boson-to-qubit layout."""

from __future__ import annotations

from qiskit.circuit import QuantumRegister

from qiskit_bosons.circuit import BosonicDAGCircuit
from qiskit_bosons.mappers import LevelSpec

from ... import B2QLayout, BosonicDAGCircuitPass


class OccupationB2QLayout(BosonicDAGCircuitPass):
    r"""Maps bosonic modes to the qubits a :class:`.LevelSpec` assigns them.

    Each bosonic mode is truncated to the levels its cutoff in :attr:`spec` keeps and spans
    :meth:`.LevelEncoding.num_qubits` of that many levels consecutive qubits, so the register of
    :math:`M` modes is mapped to a qubit register of :attr:`.LevelSpec.total_width` qubits and mode
    :math:`m` starts at qubit :attr:`.LevelSpec.offsets` ``[m]``.

    A boson-to-qubit encoding generally needs more than one qubit per mode (:math:`d` of them for the
    unary encodings), so the qubit register is wider than the mode register. Unlike
    :class:`.CustomB2QLayout` this derives the registers from the circuit being transpiled rather
    than requiring them up front, which is what lets one pass manager serve circuits of different
    widths, as long as :attr:`spec` leaves the number of modes open.

    .. note::
       The truncation and the encoding are given by the :attr:`spec` of this **pass**, not by the
       circuit. That is the point: the bosonic circuit upstream of the transpiler is independent of
       both, so the same
       circuit can be lowered at several truncation levels -- or under several encodings -- by
       swapping this pass alone.

    .. warning::
       A spec whose modes have **different** cutoffs fixes the qubit offset of each mode by its
       position in the register, so reordering the modes afterwards would silently attach the wrong
       cutoff to a mode. This pass therefore refuses a circuit to which :class:`.RelabelModes` has
       already been applied (it leaves a ``"permutation"`` entry in the circuit's metadata) and
       raises a :class:`ValueError`. A spec with one shared cutoff is unaffected. Note that
       :class:`.CustomB2QLayout` does not perform this check: the user is responsible for the
       consistency of a layout they provide themselves.
    """

    def __init__(self, spec: LevelSpec) -> None:
        """Initializing this transpiler pass can be done with the arguments listed below.

        Args:
            spec: the :class:`.LevelSpec` with the cutoff of every mode and the encoding that lays
                the occupation levels out on qubits. A spec with one cutoff per mode only supports
                circuits with exactly that many modes.
        """
        super().__init__()

        self.spec = spec
        """The :class:`.LevelSpec` this layout is derived from."""

    def run(self, dag: BosonicDAGCircuit) -> BosonicDAGCircuit:
        """Runs this analysis pass.

        Args:
            dag: the DAG circuit representation of a :class:`.BosonicCircuit`.

        Returns:
            The unchanged input. But the :attr:`.property_set` of this transpilation pipeline will
            be updated.

        Raises:
            ValueError: if a register has a different number of modes than :attr:`spec` covers, or
                if :attr:`spec` gives the modes different cutoffs and :class:`.RelabelModes` has
                been applied to the circuit.
        """
        if not self.spec.is_uniform and "permutation" in dag.metadata:
            raise ValueError(
                "the circuit has been relabeled by RelabelModes, but the LevelSpec gives its "
                "modes different cutoffs, which would then be attached to the wrong modes. "
                "Use a spec with a single shared cutoff, or do not relabel the modes."
            )

        layout: B2QLayout = {}
        for breg in dag.qregs.values():
            # resolving a spec with one cutoff per mode checks the number of modes, and one with a
            # shared cutoff is expanded to this register's width
            resolved = self.spec.resolve(len(breg))
            layout[breg] = QuantumRegister(resolved.total_width, "q")

        self.property_set["b2q_layout"] = layout

        return dag
