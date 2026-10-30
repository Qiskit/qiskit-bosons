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

"""The preset transpiler pipeline based on an occupation boson-to-qubit encoding."""

from __future__ import annotations

import functools

from qiskit.passmanager import MultiStagePassManager
from qiskit.transpiler import generate_preset_pass_manager

from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding

from .. import BosonicCircuitToDAG, BosonicPassManager, QuantumDAGToCircuit
from ..passes import (
    B2BSynthesis,
    B2QSynthesis,
    B2QSynthesisConfig,
    OccupationB2QLayout,
)


def generate_preset_occupation_pass_manager(spec: LevelSpec, **kwargs) -> MultiStagePassManager:
    """Generates a preset transpiler pipeline for the occupation basis, from a :class:`.LevelSpec`.

    The truncation and the encoding are given by the ``spec`` of this function rather than by the
    circuit: a bosonic circuit is encoding-agnostic, and choosing them here is what lets the same
    circuit be lowered at several truncation levels without being rebuilt.

    The layout, the state preparation and the operator mapping of the pipeline are all derived from
    the one ``spec``, so they cannot disagree about which qubits a mode occupies.

    .. note::
       A spec built from a single cutoff leaves the number of modes open, so the pipeline serves
       circuits of any width. A spec with one cutoff per mode fixes it, and the pipeline then only
       supports circuits with exactly that many modes.

    Args:
        spec: the :class:`.LevelSpec` with the cutoff of every mode and the encoding that lays the
            occupation levels out on qubits.
        kwargs: any additional keyword arguments are forwarded to
            :external:func:`~qiskit.transpiler.generate_preset_pass_manager` whose output is used
            for the ``qubit`` stage. Without any, that stage uses Qiskit's defaults.

    Returns:
        The preset staged boson-to-qubit transpiler pipeline.
    """
    # A `ModeRotation` is broken down into `Evolution` gates here, in bosonic space, and those are
    # what the synthesis stage maps. Further bosonic-level optimizations can be attached to this
    # stage (for example `pm.optimization.append(...)`) without altering the structure of the
    # pipeline.
    optimization = BosonicPassManager(B2BSynthesis({"ModeRotation": "Clements"}))

    layout = BosonicPassManager(OccupationB2QLayout(spec))

    mapper_fn = functools.partial(_map, spec=spec)

    config: B2QSynthesisConfig = {
        "Evolution": ("MapperFn", (mapper_fn,)),
        "InitializeModes": ("Occupation", (spec,)),
    }

    synth = B2QSynthesis(config)

    # `MultiStagePassManager` runs its stages in the order given, and the `qubit` stage works on the
    # DAG, so it has to come before `output` converts the DAG into a circuit.
    return MultiStagePassManager(
        input=BosonicCircuitToDAG(),
        optimization=optimization,
        layout=layout,
        synthesis=synth,
        qubit=generate_preset_pass_manager(**kwargs),
        output=QuantumDAGToCircuit(),
    )


def _map(operator, num_qubits, *, spec):
    """Adapts :func:`.occupation_encoding` to the ``(operator, num_qubits)`` mapper-fn signature.

    This is deliberately private and bound to the ``spec`` of the preset: the layout, the state
    preparation and the mapper all take that one object, so a mapper supplied from outside could
    never disagree with them. It is a module-level Python function, rather than a lambda or a
    ``partial`` of :func:`.occupation_encoding`, so that it can be pickled.

    ``num_qubits`` is the width of the register the layout derived from the same ``spec``. The
    mapper checks that the spec fits into it, so that a mismatch between the layout and the encoding
    surfaces as a clear error here rather than as a confusing Pauli-index error.
    """
    return occupation_encoding(operator, spec, num_qubits)
