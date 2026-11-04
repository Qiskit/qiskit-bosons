########
Overview
########

.. meta::
   :description: Overview of the Qiskit Bosons guides, grouped by topic, with links to each one.

..
    _This page may *NOT* use any `.. toctree::` directives! Instead, the guides'
    toctree lives in the `docs/guides/index.rst` page. Here, all guides should
    be categorized into a section using a :ref: reference.

This page summarizes the guides that are available in addition to the extensive API documentation of
this package.

To get started with the breadth of features this package supports, check out this end-to-end guide:

- :ref:`bose_hubbard_getting_started`

The following sections outline additional guides that dive into more specific topics in detail.


Operators and their structure
=============================

A bosonic operator is a sum of ordered products of creation and annihilation operators, and a mode
can appear more than once in a term. Bosonic operators can also tag their terms with group indices,
which later stages use to split a time evolution along physically meaningful pieces. Mappers turn
an operator into a qubit operator under a truncation and an encoding. Start here if you are
building a Hamiltonian.

- :ref:`operators_explanation`
- :ref:`grouping_explanation`
- :ref:`mappers_explanation`


Bosonic circuits
================

You can express circuits directly on bosonic modes, without choosing a truncation or an encoding.

- :ref:`bosonic_circuit_explanation`
- :ref:`mode_rotation_explanation`


Transpiling bosonic circuits
============================

Guides explaining how a :class:`.BosonicCircuit` becomes a :class:`~qiskit.circuit.QuantumCircuit`,
and how to control the Trotter error along the way.

- :ref:`transpilation_explanation`
- :ref:`bosonic_synthesis_explanation`
