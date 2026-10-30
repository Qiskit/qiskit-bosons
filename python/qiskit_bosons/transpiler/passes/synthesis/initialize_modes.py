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

"""Mode initialization synthesis."""

from __future__ import annotations

from typing import cast

from qiskit.circuit.library import XGate
from qiskit.dagcircuit import DAGCircuit, DAGOpNode

from qiskit_bosons.mappers import LevelSpec

from ... import B2QLayout
from ..utils import map_node_single_register


class OccupationInitializeModesSynthesis:
    """A :class:`.B2QSynthesisPlugin` for transpiling :class:`.InitializeModes`.

    Prepares each mode's Fock state by setting the qubits of the bit pattern that the chosen
    occupation encoding assigns to that mode's occupation level. This is a plain layer of
    :class:`~qiskit.circuit.library.XGate` on the qubits the pattern sets -- no entangling gate is
    needed for a product Fock state.

    .. warning::
       This plugin makes the following assumptions, none of which it verifies:

       - **The qubits start in the all-zero state, and the gate is the first one acting on its
         modes.** The pattern is applied with ``X`` gates, so it *toggles* the qubits rather than
         setting them. A second :class:`.InitializeModes` gate on the same modes undoes the first
         one or, depending on the encoding, leaves the register in a state that no level of the
         truncation represents.
       - **All-zero is the vacuum only under the compact encodings.** Under both unary encodings the
         pattern of level 0 sets a qubit, so a mode that no :class:`.InitializeModes` gate acts on
         starts outside the one-hot subspace. Place an :class:`.InitializeModes` gate on every mode
         of such a circuit, with :meth:`.InitializeModes.vacuum` for the empty ones.
       - **The** :class:`.LevelSpec` **is the one used by the layout and the mapper.** The layout
         derives the qubit register from its spec and the mapper maps the operators with its own;
         this plugin places the pattern at the offsets of the spec it is given. If they differ, the
         qubits are set without an error, but not the qubits that the rest of the pipeline reads as
         that occupation. The preset builds all three from one spec.

    .. note::
       This plugin needs to know the **encoding** and the **truncation**, because "which qubits
       represent occupation 2" is exactly what an encoding defines. They are given by a
       :class:`.LevelSpec` passed to the constructor rather than being read off the gate: the gate
       is encoding-agnostic by design, and deferring this choice to transpile time is what lets one
       bosonic circuit be lowered under several encodings. Each mode is checked against, and placed
       at the qubit offset of, its own entry in the spec.
    """

    def __init__(self, spec: LevelSpec) -> None:
        """Initializing this plugin can be done with the arguments listed below.

        Args:
            spec: the :class:`.LevelSpec` with which the rest of the pipeline truncates and encodes
                its modes. It must be the same spec as the layout and the mapper use.
        """
        self.spec = spec
        """The :class:`.LevelSpec` placing each mode's levels onto qubits."""

    def run(self, in_node: DAGOpNode, out_dag: DAGCircuit, *, b2q_layout: B2QLayout) -> None:
        """Runs this transpilation plugin.

        Args:
            in_node: the input boson-based circuit instruction. When this plugin gets called, the
                ``in_node.op`` attribute `must` be of type :class:`.InitializeModes`.
            out_dag: the output qubit-based circuit.
            b2q_layout: the global transpilation :class:`~qiskit_bosons.transpiler.B2QLayout`
                setting.

        .. seealso::
           The documentation of :class:`.B2QSynthesisPlugin` for more detailed explanations of the
           arguments.

        Raises:
            NotImplementedError: when ``in_node`` acts on bosonic modes that are spread across
                multiple :type:`~qiskit_bosons.circuit.BosonicRegister` instances.
            ValueError: if the gate declares an occupation that the truncation of its mode cannot
                represent, or if the register has a different number of modes than :attr:`spec`
                covers.
        """
        breg_indices, qreg = map_node_single_register(in_node, b2q_layout)
        occupation = in_node.op.occupation

        # The register's width in qubits is what the layout derived from the same spec, so resolving
        # the spec against the number of modes that width stands for checks the two agree. A spec
        # with one shared cutoff is expanded to as many modes as the gate reaches.
        spec = self.spec.resolve(max(max(breg_indices, default=-1) + 1, self.spec.num_modes or 0))
        if len(qreg) < spec.total_width:
            raise ValueError(
                f"the qubit register has {len(qreg)} qubits, but the LevelSpec needs "
                f"{spec.total_width}. The layout must be derived from the same spec."
            )
        encoding = spec.encoding
        # a resolved spec always holds one cutoff per mode
        cutoffs = cast(list[int], spec.max_occupation)
        offsets = spec.offsets
        widths = spec.widths

        for local_mode, level in enumerate(occupation):
            mode = breg_indices[local_mode]
            if level > cutoffs[mode]:
                raise ValueError(
                    f"InitializeModes declares an occupation of {level} for mode {mode}, which "
                    f"exceeds the truncation max_occupation={cutoffs[mode]} this pipeline was "
                    "configured with. Raise max_occupation (and rebuild the pass manager) to "
                    "represent this state."
                )
            bits = encoding.level_to_bits(level)
            for offset in range(widths[mode]):
                if (bits >> offset) & 1:
                    out_dag.apply_operation_back(XGate(), (qreg[offsets[mode] + offset],))
