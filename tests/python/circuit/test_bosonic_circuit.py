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

"""BosonicCircuit tests.

This class deliberately exposes a *reduced* API compared to
:class:`~qiskit.circuit.QuantumCircuit`, so that qubit gates cannot be applied to bosonic modes.
These tests pin that gate keeping and the group-aware decomposition of an ``Evolution``.
"""

from __future__ import annotations

import unittest

import numpy as np
import pytest
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import XGate
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes, ModeRotation
from qiskit_bosons.operators import BosonOperator, ann, cre


def bose_hubbard(num_sites: int, tunneling: float = 1.0, interaction: float = 2.0):
    """The normal-ordered Bose-Hubbard chain, grouped one group per bond plus one diagonal group."""
    edges = [(site, site + 1) for site in range(num_sites - 1)]
    terms: dict = {}
    for left, right in edges:
        terms[(cre(left), ann(right))] = -tunneling
        terms[(cre(right), ann(left))] = -tunneling
    for site in range(num_sites):
        terms[(cre(site), cre(site), ann(site), ann(site))] = 0.5 * interaction

    operator = BosonOperator.from_dict(terms)

    groups = []
    for term, _ in operator.iter_terms():
        modes = {mode for _, mode in term}
        if len(modes) == 1:
            groups.append(len(edges))  # the shared diagonal group
        else:
            groups.append(edges.index(tuple(sorted(modes))))
    operator.groups = groups
    return operator


class TestConstruction(unittest.TestCase):
    def test_modes(self):
        circuit = BosonicCircuit(4)
        assert len(circuit.modes) == 4
        assert len(circuit.register) == 4

    def test_zero_modes(self):
        circuit = BosonicCircuit(0)
        assert circuit.modes == []

    def test_starts_empty(self):
        assert dict(BosonicCircuit(3).count_ops()) == {}

    def test_metadata_round_trip(self):
        circuit = BosonicCircuit(2)
        circuit.metadata = {"truncation": "chosen later"}
        assert circuit.metadata == {"truncation": "chosen later"}


class TestAppendGateKeeping(unittest.TestCase):
    """Only :class:`.BosonicGate` subclasses may be appended."""

    def test_rejects_a_qubit_gate(self):
        circuit = BosonicCircuit(2)
        with pytest.raises(ValueError, match="Unsupported instruction type"):
            circuit.append(XGate(), [circuit.modes[0]])

    def test_rejects_a_plain_quantum_circuit_instruction(self):
        circuit = BosonicCircuit(2)
        inner = QuantumCircuit(1)
        inner.h(0)
        with pytest.raises(ValueError, match="Unsupported instruction type"):
            circuit.append(inner.to_gate(), [circuit.modes[0]])

    def test_rejects_a_non_gate(self):
        circuit = BosonicCircuit(2)
        with pytest.raises(ValueError, match="Unsupported instruction type"):
            circuit.append("not a gate", [circuit.modes[0]])  # type: ignore[arg-type]

    def test_accepts_bosonic_gates(self):
        circuit = BosonicCircuit(3)
        circuit.append(InitializeModes([1, 0, 2]), circuit.modes)
        circuit.append(ModeRotation.beamsplitter(0.3), circuit.modes[:2])
        assert dict(circuit.count_ops()) == {"InitializeModes": 1, "ModeRotation": 1}

    def test_accepts_a_gate_on_a_subset_of_modes(self):
        circuit = BosonicCircuit(4)
        circuit.append(InitializeModes([1, 2]), circuit.modes[1:3])
        assert dict(circuit.count_ops()) == {"InitializeModes": 1}

    def test_rejection_leaves_the_circuit_untouched(self):
        circuit = BosonicCircuit(2)
        with pytest.raises(ValueError):
            circuit.append(XGate(), [circuit.modes[0]])
        assert dict(circuit.count_ops()) == {}


class TestCountOps(unittest.TestCase):
    def test_counts_repeated_gates(self):
        circuit = BosonicCircuit(2)
        for _ in range(3):
            circuit.append(ModeRotation.beamsplitter(0.1), circuit.modes)
        assert dict(circuit.count_ops()) == {"ModeRotation": 3}

    def test_counts_distinct_gates(self):
        circuit = BosonicCircuit(2)
        circuit.append(InitializeModes([1, 0]), circuit.modes)
        circuit.append(
            Evolution(2, BosonOperator.from_dict({(cre(0), ann(0)): 1.0})), circuit.modes
        )
        assert dict(circuit.count_ops()) == {"InitializeModes": 1, "Evolution": 1}


class TestDecompose(unittest.TestCase):
    """A grouped ``Evolution`` decomposes into one sub-``Evolution`` per group."""

    def test_grouped_evolution_splits_per_group(self):
        num_sites = 3
        hamiltonian = bose_hubbard(num_sites)
        assert hamiltonian.num_groups() == 3

        circuit = BosonicCircuit(num_sites)
        circuit.append(InitializeModes([1, 0, 1]), circuit.modes)
        circuit.append(Evolution(num_sites, hamiltonian, time=0.4), circuit.modes)
        assert dict(circuit.count_ops()) == {"InitializeModes": 1, "Evolution": 1}

        decomposed = circuit.decompose()
        assert dict(decomposed.count_ops()) == {"InitializeModes": 1, "Evolution": 3}

    def test_each_group_is_hermitian(self):
        """The whole point of grouping: ``PauliEvolutionGate`` needs a Hermitian generator."""
        hamiltonian = bose_hubbard(3)
        for group in hamiltonian.split_out_groups():
            assert group.is_hermitian()

    def test_ungrouped_evolution_splits_per_term(self):
        operator = BosonOperator.from_dict(
            {(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0, (cre(0), ann(0)): 1.0}
        )
        assert not operator.has_groups()

        circuit = BosonicCircuit(2)
        circuit.append(Evolution(2, operator, time=0.2), circuit.modes)
        decomposed = circuit.decompose()
        assert dict(decomposed.count_ops()) == {"Evolution": 3}

    def test_grouping_is_ignored_without_a_decompose(self):
        """Skip the decomposition and the grouping simply never reaches the circuit."""
        circuit = BosonicCircuit(3)
        circuit.append(Evolution(3, bose_hubbard(3), time=0.4), circuit.modes)
        assert dict(circuit.count_ops()) == {"Evolution": 1}

    def test_decompose_preserves_the_register(self):
        circuit = BosonicCircuit(3)
        circuit.append(Evolution(3, bose_hubbard(3), time=0.4), circuit.modes)
        decomposed = circuit.decompose()
        assert decomposed.register is circuit.register
        assert len(decomposed.modes) == 3

    def test_decompose_returns_a_new_circuit(self):
        circuit = BosonicCircuit(3)
        circuit.append(Evolution(3, bose_hubbard(3), time=0.4), circuit.modes)
        decomposed = circuit.decompose()
        assert decomposed is not circuit
        # the original still holds a single Evolution
        assert dict(circuit.count_ops()) == {"Evolution": 1}

    def test_decompose_reps(self):
        """The per-group factors are atomic, so further repetitions reach a fixed point.

        Splitting a Hermitian group further, term by term, would produce individually non-Hermitian
        factors whose exponentials are not unitary; marking the factors atomic prevents that.
        """
        hamiltonian = bose_hubbard(3)
        circuit = BosonicCircuit(3)
        circuit.append(Evolution(3, hamiltonian, time=0.4), circuit.modes)
        once = circuit.decompose()
        twice = circuit.decompose(reps=2)
        assert once.count_ops()["Evolution"] == 3
        assert twice.count_ops() == once.count_ops()
        assert circuit.decompose(reps=5).count_ops() == once.count_ops()

    def test_decompose_targets_a_named_gate(self):
        circuit = BosonicCircuit(3)
        circuit.append(InitializeModes([1, 0, 1]), circuit.modes)
        circuit.append(Evolution(3, bose_hubbard(3), time=0.4), circuit.modes)
        decomposed = circuit.decompose(gates_to_decompose="Evolution")
        assert dict(decomposed.count_ops()) == {"InitializeModes": 1, "Evolution": 3}

    def test_sub_evolutions_act_on_the_active_modes_only(self):
        """Each group is relabeled onto just the modes its terms touch."""
        num_sites = 3
        circuit = BosonicCircuit(num_sites)
        circuit.append(Evolution(num_sites, bose_hubbard(num_sites), time=0.4), circuit.modes)
        decomposed = circuit.decompose()

        widths = sorted(instruction.operation.num_modes for instruction in decomposed._inner.data)
        # two bond groups on 2 modes each, one diagonal group spanning all 3
        assert widths == [2, 2, 3]

    def test_sub_evolutions_carry_the_same_time(self):
        time = 0.37
        circuit = BosonicCircuit(3)
        circuit.append(Evolution(3, bose_hubbard(3), time=time), circuit.modes)
        for instruction in circuit.decompose()._inner.data:
            assert instruction.operation.params[0] == pytest.approx(time)

    def test_decomposing_the_grouped_operator_conserves_the_total(self):
        """The sub-evolutions' operators reassemble into the original Hamiltonian."""
        num_sites = 3
        hamiltonian = bose_hubbard(num_sites)
        circuit = BosonicCircuit(num_sites)
        circuit.append(Evolution(num_sites, hamiltonian, time=0.4), circuit.modes)

        total = BosonOperator.zero()
        for instruction in circuit.decompose()._inner.data:
            local = instruction.operation.operator
            # map the sub-evolution's local modes back onto the global register
            global_modes = [
                circuit.decompose()._inner.find_bit(mode).index for mode in instruction.qubits
            ]
            permutation = list(global_modes) + [
                mode for mode in range(num_sites) if mode not in global_modes
            ]
            total += local.relabel_modes(permutation)

        assert total.simplify().equiv(hamiltonian)


class TestDraw(unittest.TestCase):
    def test_text_drawing_mentions_the_gates(self):
        circuit = BosonicCircuit(2)
        circuit.append(InitializeModes([1, 0]), circuit.modes)
        # Qiskit's text drawer normalizes the capitalization of a gate name
        drawing = str(circuit.draw(output="text")).lower()
        assert "initializemodes" in drawing

    def test_text_drawing_labels_the_mode_register(self):
        drawing = str(BosonicCircuit(2).draw(output="text"))
        assert "b_0" in drawing and "b_1" in drawing


class TestEvolutionGate(unittest.TestCase):
    def test_carries_no_cutoff(self):
        """The gate is defined purely by the operator; the truncation arrives at transpile time."""
        gate = Evolution(2, BosonOperator.from_dict({(cre(0), ann(0)): 1.0}), time=0.5)
        assert not hasattr(gate, "max_occupation")
        assert gate.params[0] == pytest.approx(0.5)

    def test_default_time_is_one(self):
        gate = Evolution(1, BosonOperator.from_dict({(cre(0), ann(0)): 1.0}))
        assert gate.params[0] == pytest.approx(1.0)

    def test_negative_time_evolves_backwards(self):
        gate = Evolution(1, BosonOperator.from_dict({(cre(0), ann(0)): 1.0}), time=-0.5)
        assert gate.params[0] == pytest.approx(-0.5)

    def test_hermiticity_is_not_verified(self):
        """Documented as the caller's responsibility."""
        non_hermitian = BosonOperator.from_dict({(cre(0), ann(1)): 1.0})
        assert not non_hermitian.is_hermitian()
        assert Evolution(2, non_hermitian).num_modes == 2

    def test_operator_is_retained(self):
        operator = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        assert Evolution(1, operator).operator is operator

    def test_num_modes_may_exceed_the_operator_support(self):
        """Idle modes are padded, exactly as the mapper pads them."""
        operator = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        gate = Evolution(4, operator)
        assert gate.num_modes == 4
        assert np.array_equal(sorted(operator.get_support()), [0])
