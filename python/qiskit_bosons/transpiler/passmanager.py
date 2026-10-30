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

"""The passmanager interface for bosonic circuits."""

from qiskit.passmanager import BasePassManager

from qiskit_bosons.circuit import BosonicCircuit, BosonicDAGCircuit

from .converters import BosonicCircuitToDAG, BosonicDAGToCircuit


class BosonicPassManager(BasePassManager[BosonicDAGCircuit]):
    """A transpiler pass manager converting one :class:`.BosonicCircuit` to another.

    The internal representation (IR) of this pass manager is :class:`.BosonicDAGCircuit`, but the
    input and output type are :class:`.BosonicCircuit`.
    """

    def _passmanager_frontend(self, input_program: BosonicCircuit, **kwargs) -> BosonicDAGCircuit:
        return BosonicCircuitToDAG().run(input_program)

    def _passmanager_backend(
        self, passmanager_ir: BosonicDAGCircuit, in_program: BosonicCircuit, **kwargs
    ) -> BosonicCircuit:
        return BosonicDAGToCircuit().run(passmanager_ir)
