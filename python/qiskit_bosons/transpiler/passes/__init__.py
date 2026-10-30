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
=================
Transpiler Passes
=================

.. currentmodule:: qiskit_bosons.transpiler.passes

This module provides various transpiler passes for the stages explained in
:mod:`qiskit_bosons.transpiler`.

.. _qiskit_bosons-transpiler-passes-optimization:

Optimization Passes
-------------------

These passes provide different kinds of optimization of :class:`.BosonicDAGCircuit` instances.

.. autosummary::
   :toctree: ../stubs/

   B2BSynthesis
   ~optimization.b2b_synthesis.B2BSynthesisConfig
   MergeEvolutions
   RelabelModes

.. _qiskit_bosons-transpiler-passes-layout:

Layouting Passes
----------------

These passes are designed to help configuring the global
:type:`~qiskit_bosons.transpiler.B2QLayout` setting for the transpilation process. This setting
needs to be placed in the ``b2q_layout`` field of the
:attr:`~qiskit.passmanager.PassManagerState.property_set`, from where it will be read during the
:ref:`qiskit_bosons-transpiler-passes-synthesis`.

.. autosummary::
   :toctree: ../stubs/

   CustomB2QLayout
   OccupationB2QLayout

.. _qiskit_bosons-transpiler-passes-synthesis:

Synthesis Passes
----------------

The main logic for mapping a :class:`.BosonicDAGCircuit` to a
:class:`~qiskit.dagcircuit.DAGCircuit` is implemented by a single synthesis pass, namely the
:class:`.B2QSynthesis`. It is conceptually similar to Qiskit's
:class:`~qiskit.transpiler.passes.HighLevelSynthesis` because it simply iterates over the
:class:`.BosonicDAGCircuit` instructions and delegates the mapping to qubit-based instructions to
:ref:`qiskit_bosons-transpiler-passes-synthesis-plugins` for each of the encountered types of
:class:`.BosonicGate`.

.. autosummary::
   :toctree: ../stubs/

   B2QSynthesis
   ~synthesis.b2q_synthesis.B2QSynthesisConfig
"""

from .layout import CustomB2QLayout, OccupationB2QLayout
from .optimization import (
    B2BSynthesis,
    B2BSynthesisConfig,
    MergeEvolutions,
    RelabelModes,
)
from .synthesis import (
    B2QSynthesis,
    B2QSynthesisConfig,
    B2QSynthesisPlugin,
    B2QSynthesisPluginManager,
    MapperFnEvolutionSynthesis,
    OccupationInitializeModesSynthesis,
    group_wise,
    simplify,
)

__all__ = [
    "B2BSynthesis",
    "B2BSynthesisConfig",
    "B2QSynthesis",
    "B2QSynthesisConfig",
    "B2QSynthesisPlugin",
    "B2QSynthesisPluginManager",
    "CustomB2QLayout",
    "MapperFnEvolutionSynthesis",
    "MergeEvolutions",
    "OccupationB2QLayout",
    "OccupationInitializeModesSynthesis",
    "RelabelModes",
    "group_wise",
    "simplify",
]
