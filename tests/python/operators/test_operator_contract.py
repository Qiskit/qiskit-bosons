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

"""Runs the shared operator contract against :class:`.BosonOperator`."""

from __future__ import annotations

from qiskit_bosons.operators import BosonOperator

from .operator_contract_tests import OperatorContractTests


class TestBosonOperatorContract(OperatorContractTests):
    """The arithmetic and iteration contract, as inherited from the mixin."""

    @staticmethod
    def get_class() -> type[BosonOperator]:
        return BosonOperator
