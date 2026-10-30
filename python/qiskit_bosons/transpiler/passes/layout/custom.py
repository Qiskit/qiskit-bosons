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

"""A custom boson-to-qubit layout."""

from __future__ import annotations

from qiskit_bosons.circuit import BosonicDAGCircuit

from ... import B2QLayout, BosonicDAGCircuitPass


class CustomB2QLayout(BosonicDAGCircuitPass):
    """Sets the user-provided :type:`~qiskit_bosons.transpiler.B2QLayout` in the transpiler.

    This pass simply populates the ``b2q_layout`` field of the
    :attr:`~qiskit.passmanager.PassManagerState.property_set`, with the value of :attr:`layout`.

    .. note::
       It is the user's responsibility to ensure that subsequently used transpilation passes (and
       plugins) are configured to adhere to this :type:`~qiskit_bosons.transpiler.B2QLayout`
       correctly.
    """

    def __init__(self, layout: B2QLayout) -> None:
        """Initializing this transpiler pass can be done with the arguments listed below.

        Args:
            layout: the custom ``b2q_layout`` to use.
        """
        super().__init__()

        self.layout: B2QLayout = layout
        """The user-provided mapping of :type:`~qiskit_bosons.circuit.BosonicRegister` to
        :class:`~qiskit.circuit.QuantumRegister`."""

    def run(self, dag: BosonicDAGCircuit) -> BosonicDAGCircuit:
        """Runs this analysis pass.

        Args:
            dag: the DAG circuit representation of a :class:`.BosonicCircuit`.

        Returns:
            The unchanged input. But the :attr:`.property_set` of this transpilation pipeline will
            be updated.
        """
        self.property_set["b2q_layout"] = self.layout

        return dag
