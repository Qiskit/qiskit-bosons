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

"""Boson-to-qubit circuit synthesis pass."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, TypeAlias, cast

from qiskit.circuit import Barrier
from qiskit.dagcircuit import DAGCircuit
from qiskit.passmanager import GenericPass

from qiskit_bosons.circuit import BosonicDAGCircuit, BosonicGate

from ... import B2QLayout
from .plugin import B2QSynthesisPlugin, B2QSynthesisPluginManager

B2QSynthesisConfig: TypeAlias = dict[
    str,
    str
    | tuple[str, Sequence[Any]]
    | tuple[str, Sequence[Any], Mapping[str, Any]]
    | B2QSynthesisPlugin,
]
"""The dictionary type used to configure the :attr:`.B2QSynthesis.methods`.

The keys of this dictionary must be names of :class:`.BosonicGate` circuit instructions.
The values can be one of four types:

1. ``str``: the simplest scenario simply specifies the name of the plugin method to use for the
   corresponding key in this dictionary. The plugin is extracted from the
   :class:`.B2QSynthesisPluginManager` and its ``__init__`` method may **not** require any
   arguments.
2. ``tuple[str, Sequence[Any]]``: the first ``str`` is the same as the above, and the second
   ``Sequence[Any]`` can be used to provide any positional arguments to the plugin's ``__init__``
   method.
3. ``tuple[str, Sequence[Any], Mapping[str, Any]]``: the first two fields are the same as in 2. and
   the third ``Mapping[str, Any]`` can be used to provide any keyword arguments to the plugin's
   ``__init__`` method.
4. an instance of a :class:`.B2QSynthesisPlugin`, which is used as is. This is how a plugin that is
   not registered through an entry point enters the pass.
"""


class B2QSynthesis(GenericPass[BosonicDAGCircuit, DAGCircuit]):
    """A transpilation pass to map boson-based circuit instructions to qubit-based ones.

    This transpilation pass works similarly to Qiskit's
    :class:`~qiskit.transpiler.passes.HighLevelSynthesis` pass; given an input
    :class:`.BosonicDAGCircuit`, it iterates the contained instructions and delegates the
    translation to qubit-based instructions to matching :attr:`methods`. The insertion of the
    qubit-based circuit instructions into the output :class:`~qiskit.dagcircuit.DAGCircuit` is also
    left to the plugin method. This pass will merely have prepared the
    :class:`~qiskit.circuit.QuantumRegister` according to the global transpilation
    :class:`~qiskit_bosons.transpiler.B2QLayout` setting.

    .. rubric:: Usage

    There are two ways to configure the plugins used by this transpiler pass. The examples below
    show the equivalent configuration to the :func:`.generate_preset_occupation_pass_manager`.

    1. Providing a ``config`` during initialization:

       .. doctest::

          >>> from qiskit_bosons.mappers import LevelSpec
          >>> from qiskit_bosons.mappers.library import occupation_encoding
          >>> from qiskit_bosons.transpiler import passes
          >>>
          >>> # The encoding and the truncation are bound here, at transpile time -- the bosonic
          >>> # circuit upstream of this pass is independent of both.
          >>> spec = LevelSpec(2, "unary_projected")
          >>> def mapper_fn(operator, num_qubits):
          ...     return occupation_encoding(operator, spec, num_qubits)
          >>> config = {
          ...     "Evolution": ("MapperFn", (mapper_fn,)),
          ...     "InitializeModes": ("Occupation", (spec,)),
          ... }
          >>> synth = passes.B2QSynthesis(config)

    2. Manually populating the :attr:`methods` attribute:

       .. doctest::

          >>> from qiskit_bosons.transpiler import passes
          >>>
          >>> synth = passes.B2QSynthesis()
          >>> synth.methods["Evolution"] = passes.MapperFnEvolutionSynthesis(mapper_fn)
          >>> synth.methods["InitializeModes"] = (
          ...     passes.OccupationInitializeModesSynthesis(spec)
          ... )

       Through this manual approach it is also possible to inject custom
       :class:`.B2QSynthesisPlugin` implementations without requiring them to be registered through
       an `entry-point <https://setuptools.pypa.io/en/latest/userguide/entry_point.html>`_.

    """

    def __init__(self, config: B2QSynthesisConfig | None = None) -> None:
        """Initializing this transpiler pass can be done with the arguments listed below.

        Args:
            config: an optional dictionary to pre-populate the :attr:`methods` with plugins provided
                by the :class:`.B2QSynthesisPluginManager`, or with plugin instances. See
                :type:`.B2QSynthesisConfig`.

        Raises:
            KeyError: if ``config`` names an instruction or a method for which no plugin is
                installed.
            TypeError: if ``config`` gives a value that is none of the forms of
                :type:`.B2QSynthesisConfig`.
        """
        super().__init__()

        self.methods: dict[str, B2QSynthesisPlugin] = {}
        """A dictionary of boson-to-qubit circuit instruction transpilation methods.

        For this transpilation pass to have any effect, this dictionary must be populated with
        instances of the :class:`.B2QSynthesisPlugin` protocol. The keys of this dictionary
        correspond to the ``__name__`` of a circuit instruction.

        All available transpilation plugins are managed by the :class:`.B2QSynthesisPluginManager`.

        .. note::
           This dictionary is **empty** by default! You can pre-populate it with plugins provided by
           the :class:`.B2QSynthesisPluginManager` by providing a ``config`` argument during
           initialization of this transpiler pass instance.
        """

        if config is not None:
            pm = B2QSynthesisPluginManager()
            for op_name, method_spec in config.items():
                init_args: Sequence[Any] = ()
                init_kwargs: Mapping[str, Any] = {}
                match method_spec:
                    case str():
                        method_name = method_spec
                    case (str() as name, args):
                        method_name = name
                        init_args = args
                    case (str() as name, args, kwargs):
                        method_name = name
                        init_args = args
                        init_kwargs = kwargs
                    case _ if callable(getattr(method_spec, "run", None)):
                        self.methods[op_name] = method_spec
                        continue
                    case _:
                        raise TypeError(
                            f"The synthesis for '{op_name}' must be a method name, a tuple of a "
                            f"method name and its arguments, or a plugin instance, not {method_spec!r}."
                        )

                if op_name not in pm.plugins:
                    raise KeyError(
                        f"No boson-to-qubit synthesis is installed for '{op_name}'. Instructions "
                        f"that have one: {sorted(pm.op_names())}."
                    )
                if method_name not in pm.method_names(op_name):
                    raise KeyError(
                        f"No synthesis method '{method_name}' is installed for '{op_name}'. "
                        f"Installed: {sorted(pm.method_names(op_name))}."
                    )

                method = pm.method(op_name, method_name)
                self.methods[op_name] = method(*init_args, **init_kwargs)

    def run(self, dag: BosonicDAGCircuit) -> DAGCircuit:
        """Runs this transpilation pass.

        Args:
            dag: the input circuit with boson-based instructions. Only
                :class:`~qiskit.dagcircuit.DAGOpNode` with :class:`.BosonicGate` instances as their
                :attr:`~qiskit.dagcircuit.DAGOpNode.op` are supported, plus
                :class:`~qiskit.circuit.library.Barrier`, which is mapped onto every qubit of each
                register it touched.

        Returns:
            The output circuit with qubit-based instructions.

        Raises:
            ValueError: when a :class:`~qiskit.dagcircuit.DAGOpNode` is encountered whose
                :attr:`~qiskit.dagcircuit.DAGOpNode.op` is neither of type :class:`.BosonicGate` nor
                a :class:`~qiskit.circuit.library.Barrier`.
            TypeError: when a :class:`.BosonicGate` type is encountered for which no translation
                plugin is present in :attr:`methods`.
        """
        b2q_layout = cast(B2QLayout, self.property_set["b2q_layout"])

        # Build the qubit-based output DAG from scratch rather than copying the bosonic input and
        # deleting its mode wires: ``copy_empty_like`` followed by ``remove_qregs``/``remove_qubits``
        # leaves the DAG's internal wire graph in a state that intermittently trips a spurious
        # ``not a DAG`` panic in downstream analysis (e.g. Qiskit's ``Depth`` pass). Starting fresh
        # and only adding the mapped qubit registers keeps the graph consistent.
        out_dag = DAGCircuit()
        out_dag.name = dag.name
        out_dag.global_phase = dag.global_phase
        out_dag.metadata = dag.metadata
        for qreg in b2q_layout.values():
            out_dag.add_qreg(qreg)
        for creg in dag.cregs.values():
            out_dag.add_creg(creg)

        # Iterate in circuit order. `op_nodes` returns the nodes in storage order, which a pass that
        # substitutes a node by a sub-circuit (such as a decomposition) leaves out of sync with the
        # wire dependencies: the factors of two successive substitutions then interleave, and the
        # operations are appended to the output in the wrong order.
        for node in dag.topological_op_nodes():
            # Key on `base_class`, not `type(...)`: a singleton instruction is an instance of a
            # generated `_Singleton...` subclass, whose name no entry point could match.
            # `base_class` is the user-facing class, and is identical to `type(...)` for every
            # non-singleton instruction, so this is the right key in both cases.
            op_type = node.op.base_class

            # A barrier is not a BosonicGate but is supported: it carries no bosonic content, so it
            # needs no synthesis plugin, only a placement on the mapped qubits. Because a B2QLayout
            # relates whole registers (and an encoding needs several qubits per mode, so the two
            # sides differ in size), the modes a bosonic barrier covers have no per-qubit image in
            # general. Span every qubit of each mapped register it touched instead: that is correct
            # under any layout, and a barrier can span several registers, which is why this does not
            # go through `map_node_single_register`.
            if isinstance(node.op, Barrier):
                qargs = [
                    qubit
                    for breg, qreg in b2q_layout.items()
                    if any(mode in breg for mode in node.qargs)
                    for qubit in qreg
                ]
                out_dag.apply_operation_back(Barrier(len(qargs), label=node.op.label), qargs=qargs)
                continue

            if not isinstance(node.op, BosonicGate):
                raise ValueError(f"Encountered an unsupported circuit instruction type: {op_type}")

            plugin = self.methods.get(op_type.__name__, None)
            if plugin is None:
                hint = ""
                if op_type.__name__ == "ModeRotation":
                    hint = (
                        " The library ships no ModeRotation plugin: either break the gate down in "
                        "bosonic space first, for example with the B2BSynthesis pass "
                        "in the optimization stage, or register a plugin that maps it directly."
                    )
                raise TypeError(
                    "No plugin registered for transpiling a circuit instruction of type: "
                    f"{op_type}.{hint}"
                )

            plugin.run(node, out_dag, b2q_layout=b2q_layout)

        return out_dag
