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

"""Testing utilities."""

import numpy as np
from qiskit.transpiler import PassManager
from qiskit_bosons.mappers import LevelSpec
from qiskit_bosons.transpiler.presets import generate_preset_occupation_pass_manager


def preset_without_qubit_stage(spec: LevelSpec):
    """Builds the preset pipeline with its ``qubit`` stage emptied.

    The preset always ends in Qiskit's qubit-level transpilation, which synthesizes every
    ``PauliEvolutionGate`` with its product formula. Tests that read the gates the bosonic stages emit,
    or that need them as exact exponentials, run the pipeline without that stage instead.
    """
    pass_manager = generate_preset_occupation_pass_manager(spec)
    pass_manager.qubit = PassManager()
    return pass_manager


def random_unitary(dim: int, *, seed=None, dtype=complex) -> np.ndarray:
    """Draws a Haar-random unitary via the QR decomposition of a complex Gaussian matrix."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((dim, dim)).astype(dtype, copy=False)
    z += 1j * rng.standard_normal((dim, dim)).astype(dtype, copy=False)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r)
    return q * (d / np.abs(d))
