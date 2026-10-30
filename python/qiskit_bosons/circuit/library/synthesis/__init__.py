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
===================================
Circuit Library Synthesis (``B2B``)
===================================

.. currentmodule:: qiskit_bosons.circuit.library.synthesis

This module provides the synthesis methods with which a bosonic gate is broken down into a
:class:`.BosonicCircuit` of smaller bosonic gates. Pass one to the gate's ``synthesis`` argument to
choose how it gets decomposed. For an :class:`.Evolution` gate that is its time evolution:

.. code-block:: python

   Evolution(num_modes, operator, time=1.0, synthesis=BosonicLieTrotter())

or choose them for a whole circuit at once with the :class:`.B2BSynthesis` transpiler pass. Any class
with a ``synthesize(gate)`` method can leverage a plugin (see :class:`.B2BSynthesisPlugin`), and the ones
registered under the ``qiskit_bosons.circuit.synthesis`` entry point are available by name through
the :class:`.B2BSynthesisPluginManager`.

Both sides of this rewrite live in `bosonic` space, so these methods can exploit the
:attr:`~qiskit_bosons.operators.OperatorTrait.groups` of the evolved operator (see
:ref:`grouping_explanation`) -- unlike the boson-to-qubit synthesis that follows, which only sees
the mapped operator, where the grouping is no longer available. That is what lets the *outer*
Trotter error between non-commuting groups be reduced by a higher order or more repetitions, which a
product formula handed to :class:`~qiskit_bosons.transpiler.passes.MapperFnEvolutionSynthesis`
cannot do.

A :class:`.ModeRotation` gate is broken down the same way, into :class:`.Evolution` gates, one for
every two-mode rotation and single-mode phase of its mesh. Its ``synthesis`` is ``None`` by default,
which leaves the gate whole:

.. code-block:: python

   ModeRotation(unitary, synthesis=BosonicClementsDecomposition())

.. important::
   This boson-to-boson step is **optional**. An :class:`.Evolution` gate can be handed straight to
   the boson-to-qubit stage regardless of how many terms its operator holds, and so can a
   :class:`.ModeRotation`, given a plugin registered for it. Decomposing in bosonic space first is a
   choice, taken to obtain cheaper factors or to expose structure that the later stage can exploit.

.. autosummary::
   :toctree: ../stubs/

   B2BSynthesisPlugin
   B2BSynthesisPluginManager

.. _builtin_b2b_synthesis_plugins:

Builtin plugins
---------------

The table below lists the builtin plugins, organized by the :class:`.BosonicGate` that they act
upon. These are the entry points of the ``qiskit_bosons.circuit.synthesis`` group, named
``<gate>.<method>``:

.. table:: An overview of the builtin ``qiskit_bosons.circuit.synthesis`` entry-points.

   +--------------------------+--------------------+---------------------------------------+
   | Circuit instruction name | Plugin method name | Reference                             |
   +==========================+====================+=======================================+
   | :class:`.Evolution`      | ``LieTrotter``     | :class:`BosonicLieTrotter`            |
   +--------------------------+--------------------+---------------------------------------+
   | :class:`.Evolution`      | ``SuzukiTrotter``  | :class:`BosonicSuzukiTrotter`         |
   +--------------------------+--------------------+---------------------------------------+
   | :class:`.Evolution`      | ``QDrift``         | :class:`BosonicQDrift`                |
   +--------------------------+--------------------+---------------------------------------+
   | :class:`.ModeRotation`   | ``Clements``       | :class:`BosonicClementsDecomposition` |
   +--------------------------+--------------------+---------------------------------------+

.. autosummary::
   :toctree: ../stubs/
   :class: sd-d-none

   BosonicLieTrotter
   BosonicSuzukiTrotter
   BosonicQDrift
   BosonicClementsDecomposition
"""

from .clements_decomposition import BosonicClementsDecomposition
from .lie_trotter import BosonicLieTrotter
from .plugin import B2BSynthesisPlugin, B2BSynthesisPluginManager
from .qdrift import BosonicQDrift
from .suzuki_trotter import BosonicSuzukiTrotter

__all__ = [
    "B2BSynthesisPlugin",
    "B2BSynthesisPluginManager",
    "BosonicClementsDecomposition",
    "BosonicLieTrotter",
    "BosonicQDrift",
    "BosonicSuzukiTrotter",
]
