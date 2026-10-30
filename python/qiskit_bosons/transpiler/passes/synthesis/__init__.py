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

# ruff: noqa: D212,D415
"""
.. _qiskit_bosons-transpiler-passes-synthesis-plugins:

===================================
Transpiler Pass Synthesis (``B2Q``)
===================================

.. currentmodule:: qiskit_bosons.transpiler.passes.synthesis

The :class:`.B2QSynthesis` transpiler pass exposes a `plugin interface` (see
:class:`.B2QSynthesisPlugin`) through which implementations for mapping boson-based instructions
to qubit-based ones can be registered. It is the counterpart of the boson-to-boson synthesis that the
:class:`.B2BSynthesis` pass performs one stage earlier (see
:mod:`~qiskit_bosons.circuit.library.synthesis`).

The available plugins are managed by the :class:`.B2QSynthesisPluginManager`. Choosing the right
combination of synthesis plugins is crucial and the :class:`.B2QSynthesis` transpilation pass does
**not** have any default plugins defined! For this, refer to the transpiler
:mod:`~qiskit_bosons.transpiler.presets`.

.. autosummary::
   :toctree: ../stubs/

   B2QSynthesisPlugin
   B2QSynthesisPluginManager

.. _builtin_b2q_synthesis_plugins:

Builtin plugins
---------------

The plugins to use at runtime have to be registered in the :attr:`.B2QSynthesis.methods` attribute.
This can be done easily through the ``config`` argument of :class:`.B2QSynthesis` which specifies
the plugins to use from the :class:`.B2QSynthesisPluginManager` or it can be done manually by
writing into :attr:`.B2QSynthesis.methods` directly.

There is one plugin per :class:`.BosonicGate` type. Each is handed the gate and the global
:type:`~qiskit_bosons.transpiler.B2QLayout`, and inserts the equivalent qubit-based instructions.

The builtin plugins are where the boson-to-qubit **encoding** and the **truncation** enter the
pipeline: the plugins taking a ``mapper_fn`` do so through the mapper they wrap, and the
``InitializeModes`` plugin through the :class:`.LevelSpec` it is given. This is the point of
deferring the encoding to transpile time: the bosonic circuit upstream of this stage is independent
of both choices.

The table below lists the builtin plugins, organized by the :class:`.BosonicGate` that they act
upon. These are the entry points of the ``qiskit_bosons.transpiler.synthesis`` group, named
``<gate>.<method>``:

.. table:: An overview of the builtin ``qiskit_bosons.transpiler.synthesis`` entry-points.

   +---------------------------+--------------------+---------------------------------------------+
   | Circuit instruction name  | Plugin method name | Reference                                   |
   +===========================+====================+=============================================+
   | :class:`.Evolution`       | ``MapperFn``       | :class:`MapperFnEvolutionSynthesis`         |
   +---------------------------+--------------------+---------------------------------------------+
   | :class:`.InitializeModes` | ``Occupation``     | :class:`OccupationInitializeModesSynthesis` |
   +---------------------------+--------------------+---------------------------------------------+

.. autosummary::
   :toctree: ../stubs/
   :class: sd-d-none

   MapperFnEvolutionSynthesis
   OccupationInitializeModesSynthesis


Additional helper functions
---------------------------

The Pauli term order that a mapper function produces is preserved through synthesis, since a product
formula synthesizes the terms in the order it receives them. These adapters wrap a mapper function to
choose that order:

.. autosummary::
   :toctree: ../stubs/

   group_wise
   simplify
"""

from .b2q_synthesis import B2QSynthesis, B2QSynthesisConfig
from .evolution import MapperFnEvolutionSynthesis, MapperFunction, group_wise, simplify
from .initialize_modes import OccupationInitializeModesSynthesis
from .plugin import B2QSynthesisPlugin, B2QSynthesisPluginManager

__all__ = [
    "B2QSynthesis",
    "B2QSynthesisConfig",
    "B2QSynthesisPlugin",
    "B2QSynthesisPluginManager",
    "MapperFnEvolutionSynthesis",
    "MapperFunction",
    "OccupationInitializeModesSynthesis",
    "group_wise",
    "simplify",
]
