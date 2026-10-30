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

"""Boson-to-boson gate synthesis pass."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any, TypeAlias

from qiskit.converters import circuit_to_dag
from qiskit.dagcircuit import DAGOpNode

from qiskit_bosons.circuit import BosonicDAGCircuit, BosonicGate
from qiskit_bosons.circuit.library.synthesis import B2BSynthesisPlugin, B2BSynthesisPluginManager

from ... import BosonicDAGCircuitPass

B2BSynthesisConfig: TypeAlias = dict[
    str,
    str
    | tuple[str, Sequence[Any]]
    | tuple[str, Sequence[Any], Mapping[str, Any]]
    | B2BSynthesisPlugin,
]
"""The dictionary type used to configure the :attr:`.B2BSynthesis.methods`.

The keys of this dictionary must be names of :class:`.BosonicGate` classes that can be broken down
in bosonic space. The values can be one of four types:

1. ``str``: the simplest scenario simply specifies the name of the plugin method to use for the
   corresponding key in this dictionary. The plugin is extracted from the
   :class:`.B2BSynthesisPluginManager` and its ``__init__`` method may **not** require any
   arguments.
2. ``tuple[str, Sequence[Any]]``: the first ``str`` is the same as the above, and the second
   ``Sequence[Any]`` can be used to provide any positional arguments to the plugin's ``__init__``
   method.
3. ``tuple[str, Sequence[Any], Mapping[str, Any]]``: the first two fields are the same as in 2. and
   the third ``Mapping[str, Any]`` can be used to provide any keyword arguments to the plugin's
   ``__init__`` method.
4. an instance of a :class:`.B2BSynthesisPlugin`, which is used as is. This is how a plugin that is
   not registered through an entry point enters the pass.
"""


class B2BSynthesis(BosonicDAGCircuitPass):
    r"""A transpilation pass breaking gates down into smaller ones in bosonic space.

    Some gates can be rewritten into a circuit of smaller bosonic gates before they reach the
    boson-to-qubit stage, which then maps and synthesizes those. This pass does that for the gates
    named in :attr:`methods`, with the plugin that :attr:`methods` names for each, so the choice is
    made once for a whole transpilation pipeline instead of gate by gate. It is the bosonic
    counterpart of :class:`.B2QSynthesis`, one stage earlier: where that pass maps every instruction
    to qubits, this one rewrites some gates into others and leaves the rest as they are.

    .. code-block:: python

       pm.optimization = BosonicPassManager(
           [
               B2BSynthesis(
                   {
                       "Evolution": ("SuzukiTrotter", (), {"order": 2, "reps": 4}),
                       "ModeRotation": "Clements",
                   }
               )
           ]
       )

    .. rubric:: Usage

    There are two ways to configure the plugins used by this transpiler pass:

    1. Providing a ``config`` during initialization:

       .. doctest::

          >>> from qiskit_bosons.transpiler import passes
          >>>
          >>> config = {"Evolution": ("SuzukiTrotter", (), {"order": 2})}
          >>> synth = passes.B2BSynthesis(config)

    2. Manually populating the :attr:`methods` attribute:

       .. doctest::

          >>> from qiskit_bosons.circuit.library import synthesis
          >>> from qiskit_bosons.transpiler import passes
          >>>
          >>> synth = passes.B2BSynthesis()
          >>> synth.methods["Evolution"] = synthesis.BosonicSuzukiTrotter(order=2)

       Through this manual approach it is also possible to inject custom
       :class:`.B2BSynthesisPlugin` implementations without requiring them to be registered through
       an `entry-point <https://setuptools.pypa.io/en/latest/userguide/entry_point.html>`_.

    It belongs in the ``optimization`` stage, ahead of the layout and the synthesis. The method named
    here takes precedence over the :attr:`~.BosonicGate.synthesis` of an individual gate. To use the
    synthesis of the gates themselves instead, leave this pass out and call
    :meth:`~.BosonicCircuit.decompose`.

    Each selected gate is replaced by the circuit that its plugin produces. A gate whose type is not a
    key of :attr:`methods` is left untouched, as is any gate that is :attr:`~.BosonicGate.atomic`
    (also with ``apply=False``: it is not even tagged with the plugin) and any node rejected by an
    optional :attr:`filter`. A gate that has no synthesis at all, such as an
    :class:`.InitializeModes`, passes through no matter what.

    This is where the grouping of the operator pays off for the *outer* Trotter error. A product
    formula handed to :class:`.MapperFnEvolutionSynthesis` only subdivides the terms within a single
    group and never interleaves the groups, so raising its repetition count does not reduce the error
    between groups that do not commute. A formula applied here splits along the groups themselves.

    Nothing in this step depends on the truncation or the encoding, which only enter later. The same
    decomposition therefore serves every ``LevelSpec``. See
    :class:`~qiskit_bosons.circuit.library.synthesis.BosonicClementsDecomposition` for what it
    implies for the cutoff.

    .. note::
       Running the pass again is harmless: the factors a plugin emits are
       :attr:`~.BosonicGate.atomic`, and an atomic gate is a terminal factor which this pass leaves
       in place rather than splitting further. The same exemption makes the pass a fixed point over
       the output of :class:`.BosonicQDrift`, whose sampled gates are atomic because the random
       draw *is* the Trotterization it performs.

    .. note::
       Set :attr:`apply` to ``False`` to only *select* the method and leave the expansion to
       something else, such as Qiskit's :class:`~qiskit.transpiler.passes.Decompose` or a
       :meth:`~.BosonicCircuit.decompose` call. Be aware that a gate which is never expanded
       reaches the boson-to-qubit stage whole, where it is mapped without its
       :attr:`~.BosonicGate.synthesis` ever being read, and every method then produces identical
       output.

    .. note::
       Decomposing is a choice, not a requirement. A :class:`.ModeRotation` that is left whole is
       lowered by a boson-to-qubit synthesis plugin registered for it, if there is one; see
       :ref:`mode_rotation_direct_explanation`.

    .. caution::
       Not every method suits every operator. An :class:`.Evolution` whose operator groups all
       mutually commute (a Hamiltonian made only of number operators, say) is synthesized
       exactly at any order, so a higher order only adds depth. Use :attr:`filter` to restrict the
       pass to the gates that benefit.

       Expanding an operator that carries no :attr:`~qiskit_bosons.operators.OperatorTrait.groups`
       splits it term by term, and a lone hopping term :math:`b^\dagger_i b_j` is not Hermitian, so
       its exponential is not unitary. Nothing complains until the mapped operator reaches
       :class:`~qiskit.circuit.library.PauliEvolutionGate`, which raises ``ValueError: Operator
       contains complex coefficients, which are not supported``. Assigning conjugate-paired groups
       is what makes each factor Hermitian; see
       :func:`~qiskit_bosons.operators.terms.grouping.groups_are_hermitian`.
    """

    def __init__(
        self,
        config: B2BSynthesisConfig | None = None,
        *,
        apply: bool = True,
        filter: Callable[[DAGOpNode], bool] | None = None,  # noqa: A002
    ) -> None:
        """Initializing this transpiler pass can be done with the arguments listed below.

        Args:
            config: an optional dictionary to pre-populate the :attr:`methods` with plugins provided
                by the :class:`.B2BSynthesisPluginManager`. See :type:`.B2BSynthesisConfig`.
            apply: whether to expand each selected gate into the factors its plugin produces. When
                ``False``, the gate is only tagged with the plugin (unless it is atomic, which is
                never touched) and something else has to expand it later. See :attr:`apply`.
            filter: an optional predicate deciding which nodes to break down. It is called with the
                :class:`~qiskit.dagcircuit.DAGOpNode` of a gate whose type is in :attr:`methods`,
                and the node is left untouched unless it returns ``True``. If ``None`` (the default),
                every such node is selected.

        Raises:
            KeyError: if ``config`` names a gate or a method for which no plugin is installed.
            TypeError: if ``config`` gives a value that is none of the forms of
                :type:`.B2BSynthesisConfig`.
        """
        super().__init__()

        self.methods: dict[str, B2BSynthesisPlugin] = {}
        """A dictionary of the boson-to-boson synthesis methods, by gate.

        The keys of this dictionary correspond to the ``__name__`` of a :class:`.BosonicGate`, and
        the values are instances of the :class:`.B2BSynthesisPlugin` protocol. A gate that is not a
        key is not touched by this pass, so the dictionary is **empty** by default, which makes the
        pass do nothing. Pre-populate it with plugins provided by the :class:`.B2BSynthesisPluginManager`
        through the ``config`` argument.
        """

        if config is not None:
            manager = B2BSynthesisPluginManager()
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
                    case _ if callable(getattr(method_spec, "synthesize", None)):
                        self.methods[op_name] = method_spec
                        continue
                    case _:
                        raise TypeError(
                            f"The synthesis for '{op_name}' must be a method name, a tuple of a "
                            f"method name and its arguments, or a plugin instance, not {method_spec!r}."
                        )

                if op_name not in manager.plugins:
                    raise KeyError(
                        f"No boson-to-boson synthesis is installed for '{op_name}'. Gates that "
                        f"have one: {sorted(manager.op_names())}."
                    )
                if method_name not in manager.method_names(op_name):
                    raise KeyError(
                        f"No synthesis method '{method_name}' is installed for '{op_name}'. "
                        f"Installed: {sorted(manager.method_names(op_name))}."
                    )

                plugin = manager.method(op_name, method_name)
                self.methods[op_name] = plugin(*init_args, **init_kwargs)

        self.apply = apply
        """Whether to expand each selected gate into the factors its plugin produces.

        When this is ``False``, the pass only *selects* the method: it tags each selected gate with
        its plugin (see :meth:`.BosonicGate.with_synthesis`), except an :attr:`~.BosonicGate.atomic`
        gate, which is left as it is, and leaves the expansion to whatever builds the gate's
        definition later, such as Qiskit's :class:`~qiskit.transpiler.passes.Decompose` or a
        :meth:`~.BosonicCircuit.decompose` call.
        Note that a gate which is never expanded reaches the boson-to-qubit stage whole, where it is
        mapped without its :attr:`~.BosonicGate.synthesis` ever being read.
        """

        self.filter = filter
        """The predicate selecting which nodes to break down."""

    def run(self, dag: BosonicDAGCircuit) -> BosonicDAGCircuit:
        """Runs this transpilation pass.

        Every node whose gate type is in :attr:`methods` and that is accepted by :attr:`filter` is
        replaced by the gates that its plugin produces, or (when :attr:`apply` is ``False``) by an
        equivalent gate merely carrying the plugin. All other nodes are left untouched. The input DAG
        is modified in place.

        Args:
            dag: the input circuit with boson-based instructions. Only
                :class:`~qiskit.dagcircuit.DAGOpNode` with :class:`.BosonicGate` instances as their
                :attr:`~qiskit.dagcircuit.DAGOpNode.op` are supported, plus
                :class:`~qiskit.circuit.library.Barrier`, which is carried through untouched.

        Returns:
            The output circuit which is still acting on a bosonic register.
        """
        # Circuit order, not storage order: the latter is out of sync with the wire dependencies
        # after a pass substituted a node by a sub-circuit, and this pass does so itself.
        # Materialized with `list` because the loop substitutes nodes in place. Qiskit documents
        # `topological_op_nodes` as a generator, which must not be consumed while the DAG changes
        # under it; the current implementation happens to hand out an already collected sequence,
        # but that is not a guarantee. Substituting one node leaves the others valid.
        for node in list(dag.topological_op_nodes()):
            if not isinstance(node.op, BosonicGate):
                continue
            # Key on `base_class`, as `B2QSynthesis` does: it is the user-facing class, identical to
            # `type(...)` for every non-singleton gate and robust against generated subclasses.
            plugin = self.methods.get(node.op.base_class.__name__)
            if plugin is None:
                continue

            # An atomic gate is a terminal factor that must not be broken down any further, so it
            # passes through unchanged: neither expanded nor tagged with the plugin, and not even
            # offered to the filter. This is also what keeps the pass a fixed point over its own
            # output and over `BosonicQDrift`'s.
            if node.op.atomic:
                continue

            if self.filter is not None and not self.filter(node):
                continue

            if not self.apply:
                # A *copy* rather than an assignment to `node.op`: a gate instance can be shared
                # between circuits (Qiskit copies gates with a shallow `__dict__` copy), so mutating
                # it would retroactively change the synthesis of every other circuit holding it.
                dag.substitute_node(node, node.op.with_synthesis(plugin), inplace=True)
                continue

            # Synthesize directly rather than going through the gate's definition: the factors are
            # all that is wanted here, and this avoids both the definition cache and the overhead of
            # running a separate expansion pass over the whole DAG.
            expanded = circuit_to_dag(plugin.synthesize(node.op)._inner)
            dag.substitute_node_with_dag(node, expanded)

        return dag
