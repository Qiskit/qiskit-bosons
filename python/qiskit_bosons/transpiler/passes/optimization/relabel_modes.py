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

"""An optimization pass to relabel bosonic modes."""

from __future__ import annotations

from qiskit import QuantumRegister

from qiskit_bosons.circuit import BosonicDAGCircuit, BosonicRegister

from ... import BosonicDAGCircuitPass


class RelabelModes(BosonicDAGCircuitPass):
    """A transpilation pass to relabel the bosonic modes.

    This pass reorders the bosonic modes of a circuit. The reordering is described by a
    ``permutation`` list, read as a mapping from *original* to *new* mode index: original mode ``i``
    is placed at new index ``permutation[i]`` in the relabeled circuit. For example,
    ``permutation = [0, 2, 4, 1, 3, 5]`` places original mode ``1`` at new index ``2``, original
    mode ``2`` at new index ``4``, and so on. Being a permutation, each index must appear exactly
    once.

    .. rubric:: Post-processing

    The relabeling reorders the bosonic modes, which in turn influences the boson-to-qubit
    mapping chosen by a later synthesis stage (and thus the achievable circuit depth). Because the
    mode order changes, any bitstring sampled from the final circuit is expressed in the *new* mode
    order and must be mapped back to the *original* order before it can be interpreted. The
    relabeling that was actually applied is recorded in a ``permutation`` field of the returned
    :class:`.BosonicDAGCircuit`'s :attr:`~qiskit.dagcircuit.DAGCircuit.metadata`.

    .. important::
       Prefer reading the relabeling from the circuit metadata rather than from :attr:`permutation`,
       since the metadata is guaranteed to reflect what the pass actually did.

    .. warning::
       This pass cannot be combined with a :class:`.LevelSpec` that gives the modes **different**
       cutoffs. The cutoffs are attached to the modes by their position in the register, so
       reordering the modes would silently attach a cutoff to the wrong mode. :class:`.OccupationB2QLayout`
       raises a :class:`ValueError` for a relabeled circuit in that case. A spec with a single
       shared cutoff is unaffected.

    Conceptually, undoing the relabeling assigns to each original mode ``m`` the value that was
    measured for new mode ``permutation[m]``. In practice this is complicated by the fact that
    :class:`.BosonicRegister` modes and Qiskit's classical bits run in opposite (little-endian)
    order, so the mode-space gather turns into an index negation (``~idx``) followed by a final
    reversal (``[::-1]``) on the counts bitstrings.

    The example below relabels a six-mode system to an interleaved mode order, a common trick to
    reduce the implementation depth, and then undoes the relabeling on the sampled counts. It is
    synthesized with ``max_occupation=1`` under the ``"standard_binary"`` encoding, which places
    exactly one qubit per mode; this keeps the counts bitstrings in one-to-one correspondence with the
    modes, so the post-processing below is a pure mode permutation. Under a truncation or encoding
    that spends several qubits per mode, the gather has to be applied per mode block instead of per
    bit:

    .. doctest::

       >>> from qiskit.passmanager import MultiStagePassManager
       >>> from qiskit.providers.basic_provider import BasicSimulator
       >>> from qiskit_bosons.circuit import BosonicCircuit
       >>> from qiskit_bosons.circuit.library import InitializeModes
       >>> from qiskit_bosons.mappers import LevelSpec
       >>> from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
       >>> from qiskit_bosons.transpiler.passes import (
       ...     B2QSynthesis, B2QSynthesisPluginManager, OccupationB2QLayout, RelabelModes,
       ... )
       >>>
       >>> # modes 0, 1 and 3 each hold a single boson
       >>> circ = BosonicCircuit(6)
       >>> circ.append(InitializeModes([1, 1, 0, 1, 0, 0]), circ.modes)
       >>>
       >>> synth_plugins = B2QSynthesisPluginManager()
       >>> synth = B2QSynthesis()
       >>> synth.methods["InitializeModes"] = synth_plugins.method(
       ...     "InitializeModes", "Occupation"
       ... )(LevelSpec(1, "standard_binary"))
       >>>
       >>> # map the original mode order onto the interleaved one
       >>> relabel = RelabelModes(permutation=[0, 2, 4, 1, 3, 5])
       >>>
       >>> pm = MultiStagePassManager(
       ...     init=BosonicCircuitToDAG(),
       ...     optimization=relabel,
       ...     layout=OccupationB2QLayout(LevelSpec(1, "standard_binary")),
       ...     synthesis=synth,
       ...     output=QuantumDAGToCircuit(),
       ... )
       >>>
       >>> qcirc = pm.run(circ)
       >>> qcirc.measure_all()
       >>>
       >>> bit_permutation = qcirc.metadata.get("permutation")
       >>> print(bit_permutation)
       [0, 2, 4, 1, 3, 5]
       >>>
       >>> res = BasicSimulator().run(qcirc, shots=1).result()
       >>> counts = res.get_counts()
       >>> print(counts)  # measured in the interleaved ordering
       {'000111': 1}
       >>>
       >>> # undo the relabeling to recover the counts in the original mode ordering
       >>> post_processed = {
       ...     "".join(bitstring[~idx] for idx in bit_permutation)[::-1]: count
       ...     for bitstring, count in counts.items()
       ... }
       >>> print(post_processed)  # recovered in the original mode ordering
       {'001011': 1}

    The recovered bitstring matches the occupation the circuit was initialized with, read in the
    little-endian order of the classical register.
    """

    def __init__(self, permutation: list[int]) -> None:
        """Initializing this transpiler pass can be done with the arguments listed below.

        Args:
            permutation: the index permutation used to relabel the bosonic mode indices. See
                :attr:`.permutation` for details.
        """
        super().__init__()

        self.permutation = permutation
        """The index permutation used to relabel the bosonic mode indices.

        A ``list[int]`` mapping original mode index ``i`` to new mode index ``permutation[i]`` (see
        the class docstring for details). Its length has to match the number of bosonic modes of the
        circuit being transpiled, and each index has to appear exactly once.

        .. note::
           The permutation is always explicit, chosen for whatever locality the target hardware or
           encoding rewards. An occupation encoding emits no parity strings, so the span of a term
           does not drive its Pauli weight.
        """

    def run(self, dag: BosonicDAGCircuit) -> BosonicDAGCircuit:
        """Runs this transpilation pass.

        Args:
            dag: the input circuit with boson-based instructions. Only
                :class:`~qiskit.dagcircuit.DAGOpNode` with :class:`.BosonicGate` instances as their
                :attr:`~qiskit.dagcircuit.DAGOpNode.op` are supported, plus
                :class:`~qiskit.circuit.library.Barrier`, which is relabeled along with the modes it
                spans.

        Returns:
            The output circuit which is still acting on a bosonic register.

        Raises:
            NotImplementedError: when the provided input circuit has more than a single register.
        """
        if len(dag.qregs) > 1:
            raise NotImplementedError(
                "Cannot apply the RelabelModes pass to a circuit with more than a single register."
            )

        permutation = self.permutation

        out_dag = dag.copy_empty_like()
        out_dag.metadata["permutation"] = permutation.copy()

        # NOTE: we already know out_dag.qregs to contain only a single entry
        for orig_name, orig_register in out_dag.qregs.items():
            out_dag.remove_qregs(orig_register)
            relabeled_register: BosonicRegister = QuantumRegister(
                name=f"{orig_name}'", bits=[orig_register[idx] for idx in permutation]
            )
            out_dag.add_qreg(relabeled_register)

        # circuit order, not storage order: the latter is out of sync with the wire dependencies
        # after a pass substituted a node by a sub-circuit, and the output is rebuilt node by node
        for node in dag.topological_op_nodes():
            orig_indices = [orig_register.index(mode) for mode in node.qargs]
            out_dag.apply_operation_back(
                node.op, qargs=[relabeled_register[idx] for idx in orig_indices]
            )

        return out_dag
