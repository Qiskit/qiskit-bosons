Install the Python API
======================

.. meta::
   :description: Install the Qiskit Bosons Python package from PyPI with pip, or from source for
                 development, and how to verify and test the installation.

Prerequisites
^^^^^^^^^^^^^

First, create a minimal environment with only Python installed in it.
Use `Python virtual environments
<https://docs.python.org/3/tutorial/venv.html>`__.

.. code:: sh

    python3 -m venv /path/to/virtual/environment

Activate your new environment.

.. code:: sh

    source /path/to/virtual/environment/bin/activate

There are two primary ways to install this package: from PyPI or from source.
The preferred method is to install from PyPI:

Install from PyPI
^^^^^^^^^^^^^^^^^

.. important::
   ``qiskit-bosons`` has not yet released to PyPI but will soon. Stay tuned!

.. code:: sh

    pip install qiskit-bosons


Install from source
^^^^^^^^^^^^^^^^^^^

If you plan to develop in the repository you should install from source.

First, Qiskit must be installed in the same Python environment that you install
``qiskit-bosons`` in. The simplest way to do this is from PyPI, but you can
install Qiskit from source if you want to.

.. code:: console

   $ pip install qiskit

Next, clone the ``qiskit-bosons`` repository.

.. code:: sh

    git clone git@github.com:Qiskit/qiskit-bosons.git

Now, install the Rust toolchain, upgrade pip, and enter the repository. Refer to
the `Rust documentation <https://www.rust-lang.org/tools/install>`__ for
instructions on installing the toolchain.

.. code:: sh

    ### <INSTALL RUST HERE> ###
    pip install --upgrade pip
    cd qiskit-bosons

Install the remaining ``build`` dependencies.
If you plan on developing in the repository, install the ``dev`` dependencies.

.. code:: sh

    pip install --group build
    pip install --group dev  # optional

The next step is to install ``qiskit-bosons`` to the virtual environment.

.. code:: sh

    pip install .

.. hint::

   You can also perform an editable install while still compiling the
   underlying Rust crate in ``release`` mode:

   .. code:: sh

       SETUPTOOLS_RUST_CARGO_PROFILE=release pip install -e .

.. hint::

   Some features rely on optional dependencies. To install all of them,
   use the ``all`` extra:

   .. code:: sh

       pip install ".[all]"

   Refer to the ``[project.optional-dependencies]`` section of `pyproject.toml
   <https://github.com/Qiskit/qiskit-bosons/blob/main/pyproject.toml>`_
   to see which dependency groups are available.

You can optionally verify that the installation was successful. The simplest
test is to try and import one of the classes provided by the ``qiskit-bosons``
package, as in the following example:

.. code:: sh

    python -c "from qiskit_bosons.circuit import BosonicCircuit"

If this completes successfully, your installation worked.

.. hint::

   You can also run the entire Python test suite:

   .. code:: sh

       pip install --group test
       make testpython

   To additionally run the tests that exercise optional dependencies, install
   the extras and use the ``testoptional`` target:

   .. code:: sh

       pip install -e ".[all]"
       make testoptional
