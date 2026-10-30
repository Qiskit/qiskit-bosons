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

# ruff: noqa: D205,D212,D415
"""
==========
Transpiler
==========

.. currentmodule:: qiskit_bosons.transpiler

--------
Overview
--------

Transpilation is the process of rewriting a given input circuit to conform with desired criteria
such as topological and operational constraints of the hardware used to execute the final circuit.
We are not going to explain this in more detail here, and instead refer to Qiskit's documentation of
the :mod:`qiskit.transpiler` module.

The focus here lies on explaining how we achieve the transpilation of a :class:`.BosonicCircuit`
to a :class:`~qiskit.circuit.QuantumCircuit`.

------
Stages
------

Conceptually, we split the transpilation process into several stages:

================================================== ==========================================
Stage                                              Description
================================================== ==========================================
:ref:`qiskit_bosons-transpiler-stage-input`        convert to DAG data structure
:ref:`qiskit_bosons-transpiler-stage-optimization` bosonic-level optimization
:ref:`qiskit_bosons-transpiler-stage-layout`       boson-to-qubit layouting
:ref:`qiskit_bosons-transpiler-stage-synthesis`    boson-to-qubit synthesis
:ref:`qiskit_bosons-transpiler-stage-qubit`        continued transpilation on the qubit-level
:ref:`qiskit_bosons-transpiler-stage-output`       convert from DAG data structure
================================================== ==========================================

.. _qiskit_bosons-transpiler-stage-input:

Input
^^^^^

The various transpiler passes are implemented to work with :class:`.BosonicDAGCircuit` as the
underlying data structure to store the circuit operations. Such `directed acyclic graphs` (DAGs)
provide an efficient data model for the traversal and manipulation of circuits.

However, an end-user is more likely to work with a :class:`.BosonicCircuit` as it provides a more
intuitive interface and data model. As such, this input stage simply runs the
:class:`.BosonicCircuitToDAG` transpiler pass.

.. _qiskit_bosons-transpiler-stage-optimization:

Optimization
^^^^^^^^^^^^

This stage of the transpilation pipeline can implement circuit optimizations while preserving the
type of circuit to be an instance of :class:`.BosonicDAGCircuit`. As such, no qubit information is
required (or necessarily available) at this point in the transpilation pipeline.

.. _qiskit_bosons-transpiler-stage-layout:

Layout
^^^^^^

One global configuration setting for the transpilation process is the boson-to-qubit `"layout"`.
This must be provided by the user and it must match the provided boson-to-qubit mapping (in the
sense that, if a chosen mapping encodes a fixed number of bosons with a different number of
qubits, the configured layout must account for that).

In the general case, bosonic modes are not always encoded with an occupation-basis into the qubit
register. Consequently, we cannot associate a single bosonic mode with a single qubit.
Therefore, the user-provided boson-to-qubit layout (:type:`B2QLayout`) associates
:type:`~qiskit_bosons.circuit.BosonicRegister` instances with
:class:`~qiskit.circuit.QuantumRegister` ones.

.. autosummary::
   :toctree: ../stubs/
   :template: autosummary/typealias.rst

   B2QLayout

The way to configure this global setting, is by placing a :type:`B2QLayout` instance in the
``b2q_layout`` field of the :attr:`~qiskit.transpiler.PassManager.property_set`. For more details
refer to :ref:`qiskit_bosons-transpiler-passes-layout`.

.. _qiskit_bosons-transpiler-stage-synthesis:

Synthesis
^^^^^^^^^

At its core, the transpilation is handled by the :class:`.B2QSynthesis` transpiler pass. It is
conceptually similar to Qiskit's :class:`~qiskit.transpiler.passes.HighLevelSynthesis` pass, which
uses various `plugins` for transpiling high-level circuit instructions. For more details, refer to
the documentation of :class:`.B2QSynthesis` directly.

How a given :class:`.BosonicGate` can be synthesized in terms of qubit-based operations will
depend on the particular gate type as well as the user-chosen boson-to-qubit mapping. For more
details, refer to :ref:`qiskit_bosons-transpiler-passes-synthesis-plugins`.

.. _qiskit_bosons-transpiler-stage-qubit:

Qubit
^^^^^

At this point in the transpilation process, we have reached a
:class:`~qiskit.dagcircuit.DAGCircuit` instance and can continue to use Qiskit's transpilation
pipeline as one would usually.

.. hint::
   Additional transpiler passes for optimizations on the qubit-level that take into account the
   knowledge of a circuit originating from a :class:`.BosonicCircuit` may be added in the future!

.. _qiskit_bosons-transpiler-stage-output:

Output
^^^^^^

This stage implements effectively the reverse of :ref:`qiskit_bosons-transpiler-stage-input`, by
calling the :class:`.QuantumDAGToCircuit` transpiler pass.

-------------
Pass Managers
-------------

All of the stages above are the default stages by how the :mod:`~qiskit_bosons.transpiler.presets`
orchestrate the various :class:`~qiskit.passmanager.BasePassManager`s in their resulting
:class:`~qiskit.passmanager.MultiStagePassManager`.

The individual stages can be either a single pass (when they change from one internal representation
(IR) to another, such as the :ref:`qiskit_bosons-transpiler-stage-input`,
:ref:`qiskit_bosons-transpiler-stage-synthesis`, and
:ref:`qiskit_bosons-transpiler-stage-output` stages above) or a stage can be a
:class:`~qiskit.passmanager.BasePassManager` whose internal passes can be modified.

For the stages operating on bosonic circuits (:ref:`qiskit_bosons-transpiler-stage-optimization`
and :ref:`qiskit_bosons-transpiler-stage-layout`), the type of passmanager to use is
:class:`.BosonicPassManager`.

.. autosummary::
   :toctree: ../stubs/

   BosonicPassManager

Conversion Passes
^^^^^^^^^^^^^^^^^

Some very basic conversion passes are provided directly by this module:

.. autosummary::
   :toctree: ../stubs/

   BosonicCircuitToDAG
   BosonicDAGToCircuit
   QuantumDAGToCircuit

Presets
^^^^^^^

For user convenience, the :mod:`qiskit_bosons.transpiler.presets` module provides a number of
functions for quickly building pre-defined transpilation pipelines.
"""

from __future__ import annotations

from typing import TypeAlias

from qiskit.circuit import QuantumRegister
from qiskit.passmanager import GenericPass

from qiskit_bosons.circuit import BosonicDAGCircuit, BosonicRegister

from .converters import BosonicCircuitToDAG, BosonicDAGToCircuit, QuantumDAGToCircuit
from .passmanager import BosonicPassManager

BosonicDAGCircuitPass: TypeAlias = GenericPass[BosonicDAGCircuit, BosonicDAGCircuit]
"""The type definition of a bosonic-to-bosonic generic transpiler pass."""

B2QLayout: TypeAlias = dict[BosonicRegister, QuantumRegister]
"""A mapping of bosonic mode registers to quantum registers.

Users must provide a data structure of this type to the circuit transpiler. In a trivial case, every
bosonic mode gets mapped to a single qubit. However, boson-to-qubit mappings generally use several
qubits per mode, which is why we do not associate single bosonic modes with single qubits.
Crucially, the size of the registers on either side of this mapping may differ.
"""

__all__ = [
    "B2QLayout",
    "BosonicCircuitToDAG",
    "BosonicDAGToCircuit",
    "BosonicPassManager",
    "QuantumDAGToCircuit",
]
