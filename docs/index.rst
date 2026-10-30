###############
Qiskit Bosons
###############

.. meta::
   :description: Overview of Qiskit Bosons, an extension of the Qiskit SDK for bosonic systems:
                 a bosonic operator data structure, boson-to-qubit encodings with a choice of
                 truncation, bosonic circuits, and a transpilation pipeline with boson-to-boson and
                 boson-to-qubit stages that turns them into qubit circuits.

.. warning::
   This package is under active development!
   If you have feedback, please `open an issue <https://github.com/Qiskit/qiskit-bosons/issues/new/choose>`_.

``qiskit-bosons`` extends the Qiskit SDK with tools for working on bosonic
systems. Within its scope it provides the following:

- An efficient data structure for the representation and manipulation of bosonic
  operators
- A framework for implementing boson-to-qubit encodings, built around shared
  level encodings (unary with or without projectors, standard binary, Gray), and
  a library of common mappers
- A framework and library of gates for expressing bosonic circuits
- A transpilation pipeline integrated with the Qiskit transpiler process to
  synthesize bosonic circuits into qubit-based circuits

Get started
-----------

This package covers the whole path from a bosonic operator to a qubit circuit.
You can work in the bosonic regime with both, operator and circuit
representation, before you choose a truncation and a boson-to-qubit encoding
that plug into the multi-representation transpiler pipeline of Qiskit.

All of these features are explained in :doc:`several guides <guides/overview>` to
help you get started. For an end-to-end overview, start with the
:doc:`Bose-Hubbard guide <guides/bose_hubbard>`.

Technical discussion
--------------------

Design intentions
"""""""""""""""""

This package is deliberately designed to align with Qiskit: it builds on a core
implemented in Rust and provides first-party language bindings to Python. Its
API intends to draw parallels to Qiskit in order to seamlessly integrate into the
workflows of users with experience in programming Qiskit. While there exists no
C API right now, the architecture is set up to support one if the need for one
ever arises in the future. The architecture of this package is also deliberately
analogous to that of `qiskit-fermions
<https://github.com/Qiskit/qiskit-fermions>`_.

A core principle to the design of ``qiskit-bosons`` was the decoupling of its
bosonic circuit representation from its qubitized form. To be more precise: the
bosonic circuits are meaningful by themselves and do not require a mapping to
qubit space to be interpretable. Furthermore, the bosonic circuit representation
cannot make any assumptions about the truncation and boson-to-qubit encoding
applied later on. Both enter only at transpile time.

Deprecation policy
------------------

.. important::

   This policy does **not** apply until the first release is tagged.

We follow `semantic versioning <https://semver.org/>`_ and are guided by the principles in
`Qiskit's deprecation policy <https://github.com/Qiskit/qiskit/blob/main/DEPRECATION.md>`_.
We may occasionally make breaking changes in order to improve the user experience.
When possible, we will keep old interfaces and mark them as deprecated, as long as they can co-exist with the
new ones.
Each substantial improvement, breaking change, or deprecation will be documented in the
:doc:`release notes <release-notes>`.


Contributing
------------

The source code is available `on GitHub
<https://github.com/Qiskit/qiskit-bosons>`_.

The developer guide is located at `CONTRIBUTING.md
<https://github.com/Qiskit/qiskit-bosons/blob/main/CONTRIBUTING.md>`_
in the root of this project's repository.
By participating, you are expected to uphold Qiskit's `code of conduct
<https://github.com/Qiskit/qiskit/blob/main/CODE_OF_CONDUCT.md>`_.

We use `GitHub issues <https://github.com/Qiskit/qiskit-bosons/issues/new/choose>`_
for tracking requests and bugs.


Citing this package
-------------------

If you use this package in your research, use the `CITATION.bib
<https://github.com/Qiskit/qiskit-bosons/blob/main/CITATION.bib>`_ file in
this project's repository to cite the appropriate reference(s).


License
-------

`Apache License 2.0 <https://github.com/Qiskit/qiskit-bosons/blob/main/LICENSE.txt>`_


.. toctree::
  :hidden:

   Documentation Home <self>
   Installation Instructions <install>

.. toctree::
   :hidden:
   :caption: Guides

   Overview <guides/overview>
   Simulate Bose-Hubbard dynamics with a choice of encoding <guides/bose_hubbard>
   More guides <guides/index>

.. toctree::
   :hidden:
   :caption: API reference

   Python API Reference <pydoc/index>
   Release Notes <release-notes>
   GitHub <https://github.com/Qiskit/qiskit-bosons>
