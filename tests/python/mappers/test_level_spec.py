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

"""LevelSpec tests: per-mode cutoffs, widths, offsets and the checks that guard them.

A spec gives every mode its own cutoff. The mapper, the layout and the state preparation all derive
their qubit offsets from it, so these tests pin the offsets themselves, and then check that a
mapping with different cutoffs per mode is faithful to a dense truncated-Fock reference.
"""

from __future__ import annotations

import itertools
import pickle
import unittest

import numpy as np
import pytest
from qiskit.passmanager import MultiStagePassManager
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution, InitializeModes
from qiskit_bosons.mappers import LevelEncoding, LevelSpec
from qiskit_bosons.mappers.library import occupation_encoding
from qiskit_bosons.operators import BosonOperator, ann, cre
from qiskit_bosons.transpiler import BosonicCircuitToDAG, QuantumDAGToCircuit
from qiskit_bosons.transpiler.passes import RelabelModes
from qiskit_bosons.transpiler.presets import generate_preset_occupation_pass_manager

from ..utils import preset_without_qubit_stage

ENCODINGS = ("unary_projected", "standard_binary", "gray")


def hamiltonian() -> BosonOperator:
    """A three-mode chain with hopping, an interaction and a number term on every mode."""
    terms = {}
    for a, b in ((0, 1), (1, 2)):
        terms[(cre(a), ann(b))] = -1.0
        terms[(cre(b), ann(a))] = -1.0
    for mode in range(3):
        terms[(cre(mode), cre(mode), ann(mode), ann(mode))] = 0.5
        terms[(cre(mode), ann(mode))] = 0.25 * (mode + 1)
    return BosonOperator.from_dict(terms)


def dense_reference(op: BosonOperator, cutoffs: list[int]) -> np.ndarray:
    """The truncated matrix of ``op``, built from the ladder definition with one factor per mode.

    Mode ``m`` has ``cutoffs[m] + 1`` levels. Qiskit is little-endian, so mode 0 is the **last**
    kron factor.
    """
    dims = [cutoff + 1 for cutoff in cutoffs]
    ladders = []
    for dim in dims:
        creation = np.zeros((dim, dim))
        for level in range(dim - 1):
            creation[level + 1, level] = np.sqrt(level + 1)
        ladders.append((creation, creation.T))

    total = np.zeros((int(np.prod(dims)),) * 2, dtype=complex)
    for term, coeff in op.iter_terms():
        factors = [np.eye(dim, dtype=complex) for dim in dims]
        for action, mode in term:
            factors[mode] = factors[mode] @ ladders[mode][0 if action else 1]
        matrix = np.array([[1.0 + 0.0j]])
        for mode in reversed(range(len(dims))):
            matrix = np.kron(matrix, factors[mode])
        total += coeff * matrix
    return total


def physical_indices(spec: LevelSpec, encoding: LevelEncoding) -> list[int]:
    """The qubit basis indices of the in-code states, ordered like :func:`dense_reference`."""
    cutoffs = spec.max_occupation
    indices = []
    for levels in itertools.product(*(range(cutoff + 1) for cutoff in reversed(cutoffs))):
        value = 0
        for mode, level in zip(reversed(range(len(cutoffs))), levels, strict=True):
            value |= encoding.level_to_bits(level) << spec.offsets[mode]
        indices.append(value)
    return indices


class TestWidthsAndOffsets(unittest.TestCase):
    def test_per_mode_cutoffs_give_per_mode_widths(self):
        for encoding, widths in (
            ("unary_projected", [5, 3, 3]),
            ("standard_binary", [3, 2, 2]),
            ("gray", [3, 2, 2]),
        ):
            spec = LevelSpec([4, 2, 2], encoding)
            assert spec.widths == widths, encoding
            assert spec.total_width == sum(widths)

    def test_offsets_are_the_running_sum_of_the_widths(self):
        spec = LevelSpec([4, 2, 2], "gray")
        assert spec.offsets == [0, 3, 5]
        assert LevelSpec([4, 2, 2], "unary_projected").offsets == [0, 5, 8]

    def test_a_mode_without_qubits_does_not_advance_the_offset(self):
        spec = LevelSpec([3, 0, 3], "standard_binary")
        assert spec.widths == [2, 0, 2]
        assert spec.offsets == [0, 2, 2]
        assert spec.total_width == 4

    def test_a_single_integer_leaves_the_number_of_modes_open(self):
        spec = LevelSpec(2, "unary_projected")
        assert spec.is_uniform and spec.num_modes is None
        for attribute in ("widths", "offsets", "total_width"):
            with pytest.raises(ValueError, match="resolve"):
                getattr(spec, attribute)
        assert spec.resolve(3).widths == [3, 3, 3]

    def test_a_sequence_fixes_the_number_of_modes(self):
        spec = LevelSpec([1, 2, 3], "unary_projected")
        assert not spec.is_uniform
        assert spec.num_modes == 3
        assert spec.max_occupation == [1, 2, 3]

    def test_the_encoding_is_required(self):
        """There is deliberately no default encoding, so leaving it out is an error."""
        for bad in (lambda: LevelSpec(2), lambda: LevelSpec([1, 2]), lambda: LevelSpec(2, None)):
            with pytest.raises(TypeError):
                bad()  # type: ignore[call-arg, arg-type]

    def test_the_encoding_is_resolved_once(self):
        assert LevelSpec(2, "unary_projected").encoding == LevelEncoding.UNARY_PROJECTED
        assert LevelSpec(2, "unary_unprojected").encoding == LevelEncoding.UNARY_UNPROJECTED
        assert LevelSpec(2, "gray").encoding == LevelEncoding.GRAY
        assert LevelSpec(2, LevelEncoding.GRAY) == LevelSpec(2, "gray")

    def test_equality_hash_repr_and_pickle(self):
        spec = LevelSpec([1, 2], "gray")
        assert spec == LevelSpec([1, 2], "gray")
        assert spec != LevelSpec([1, 3], "gray")
        assert spec != LevelSpec([1, 2], "unary_projected")
        assert hash(spec) == hash(LevelSpec([1, 2], LevelEncoding.GRAY))
        assert pickle.loads(pickle.dumps(spec)) == spec
        assert pickle.loads(pickle.dumps(LevelSpec(3, "unary_projected"))) == LevelSpec(
            3, "unary_projected"
        )
        assert repr(spec) == "LevelSpec([1, 2], LevelEncoding.GRAY)"


class TestBroadcastingAndLengths(unittest.TestCase):
    def test_a_single_integer_broadcasts_to_any_number_of_modes(self):
        for num_modes in (1, 2, 5):
            assert (
                LevelSpec(2, "unary_projected").resolve(num_modes).max_occupation == [2] * num_modes
            )

    def test_a_sequence_of_the_wrong_length_raises(self):
        with pytest.raises(ValueError, match="covers 2 modes"):
            LevelSpec([1, 2], "unary_projected").resolve(3)

    def test_a_one_element_sequence_does_not_broadcast(self):
        spec = LevelSpec([2], "unary_projected")
        assert spec.num_modes == 1
        with pytest.raises(ValueError, match="covers 1 modes"):
            spec.resolve(3)

    def test_a_negative_or_non_integer_cutoff_is_rejected(self):
        """The cutoff is an unsigned integer, so anything else fails at the Python boundary."""
        for bad in (-1, 1.5, "two", [1, -1], [1.5]):
            with pytest.raises(TypeError):
                LevelSpec(bad, "unary_projected")  # type: ignore[arg-type]

    def test_an_unknown_encoding_raises_when_the_spec_is_built(self):
        with pytest.raises(ValueError, match="unknown encoding"):
            LevelSpec(2, "binary")


class TestRejectedSpecs(unittest.TestCase):
    def test_a_spec_without_qubits_is_rejected(self):
        """Under a compact encoding a single-level mode has width zero, so all-zero is empty."""
        for encoding in ("standard_binary", "gray"):
            for cutoffs in (0, [0], [0, 0, 0], []):
                with pytest.raises(ValueError, match="no qubits"):
                    LevelSpec(cutoffs, encoding)

    def test_the_unary_encoding_always_has_qubits(self):
        assert LevelSpec(0, "unary_projected").resolve(2).total_width == 2
        with pytest.raises(ValueError, match="no qubits"):
            LevelSpec([], "unary_projected")

    def test_a_mode_without_qubits_is_allowed_next_to_wider_ones(self):
        assert LevelSpec([0, 1], "gray").widths == [0, 1]


class TestUniformMatchesScalar(unittest.TestCase):
    """A sequence repeating one cutoff must reproduce the single-integer result exactly."""

    def test_the_mapped_operator_is_identical(self):
        for encoding in ENCODINGS:
            for cutoff in (1, 2, 3):
                uniform = occupation_encoding(hamiltonian(), LevelSpec(cutoff, encoding))
                repeated = occupation_encoding(hamiltonian(), LevelSpec([cutoff] * 3, encoding))
                assert uniform == repeated, (encoding, cutoff)
                assert uniform.num_qubits == repeated.num_qubits

    def test_the_layout_and_the_widths_are_identical(self):
        for encoding in ENCODINGS:
            uniform = LevelSpec(2, encoding).resolve(3)
            repeated = LevelSpec([2, 2, 2], encoding)
            assert (uniform.widths, uniform.offsets) == (repeated.widths, repeated.offsets)


class TestMixedCutoffsAreFaithful(unittest.TestCase):
    """Different cutoffs per mode against a dense truncated-Fock reference."""

    def test_mapped_operator_matches_the_dense_reference(self):
        for encoding in ENCODINGS:
            for cutoffs in ([1, 2, 1], [2, 1, 2], [2, 2, 1]):
                spec = LevelSpec(cutoffs, encoding)
                observable = occupation_encoding(hamiltonian(), spec)
                dense = np.asarray(SparsePauliOp.from_sparse_observable(observable).to_matrix())
                indices = physical_indices(spec, spec.encoding)
                got = dense[np.ix_(indices, indices)]
                want = dense_reference(hamiltonian(), cutoffs)
                np.testing.assert_allclose(got, want, atol=1e-12, err_msg=f"{encoding} {cutoffs}")

    def test_a_mode_without_qubits_is_faithful_too(self):
        """A frozen mode (one level) among two live ones, under a compact encoding."""
        op = BosonOperator.from_dict(
            {
                (cre(0), ann(2)): 1.0,
                (cre(2), ann(0)): 1.0,
                (cre(0), ann(0)): 0.5,
                (cre(1), ann(1)): 3.0,  # annihilates the frozen mode's only level: zero
            }
        )
        cutoffs = [2, 0, 1]
        spec = LevelSpec(cutoffs, "gray")
        observable = occupation_encoding(op, spec)
        assert observable.num_qubits == 3
        dense = np.asarray(SparsePauliOp.from_sparse_observable(observable).to_matrix())
        indices = physical_indices(spec, spec.encoding)
        np.testing.assert_allclose(
            dense[np.ix_(indices, indices)], dense_reference(op, cutoffs), atol=1e-12
        )

    def test_a_wider_register_leaves_the_extra_qubits_idle(self):
        spec = LevelSpec([2, 1, 2], "gray")
        tight = occupation_encoding(hamiltonian(), spec)
        padded = occupation_encoding(hamiltonian(), spec, num_qubits=tight.num_qubits + 3)
        assert padded.num_qubits == tight.num_qubits + 3
        # the same terms, acting on the same qubits: nothing touches the padding
        assert len(padded) == len(tight)
        assert max(int(i) for term in padded for i in term.indices) < tight.num_qubits


class TestMapperChecks(unittest.TestCase):
    def test_too_few_qubits_raises(self):
        spec = LevelSpec(2, "unary_projected")
        needed = spec.resolve(3).total_width
        with pytest.raises(ValueError, match="too few"):
            occupation_encoding(hamiltonian(), spec, num_qubits=needed - 1)
        assert occupation_encoding(hamiltonian(), spec, num_qubits=needed).num_qubits == needed

    def test_an_operator_wider_than_a_per_mode_spec_raises(self):
        with pytest.raises(ValueError, match="acts on 3 modes"):
            occupation_encoding(hamiltonian(), LevelSpec([2, 2], "unary_projected"))

    def test_a_per_mode_spec_keeps_all_of_its_modes(self):
        op = BosonOperator.from_dict({(cre(0), ann(0)): 1.0})
        assert occupation_encoding(op, LevelSpec([2, 2, 2], "unary_projected")).num_qubits == 9

    def test_a_uniform_spec_is_sized_by_the_operator(self):
        op = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
        assert occupation_encoding(op, LevelSpec(2, "unary_projected")).num_qubits == 6
        assert (
            occupation_encoding(BosonOperator.zero(), LevelSpec(2, "unary_projected")).num_qubits
            == 0
        )


class TestTruncation(unittest.TestCase):
    def test_a_cutoff_below_the_bound_reports_the_shortfall(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})  # drives mode 0 to occupation 2
        assert op.mode_occupations() == [2]
        assert LevelSpec([1], "unary_projected").truncation(op) == [1]
        assert LevelSpec([2], "unary_projected").truncation(op) == [0]
        assert LevelSpec([3], "unary_projected").truncation(op) == [0]

    def test_it_is_reported_per_mode(self):
        op = BosonOperator.from_dict(
            {(cre(0), cre(0), cre(0)): 1.0, (cre(1), ann(1)): 1.0, (cre(2), cre(2)): 1.0}
        )
        assert op.mode_occupations() == [3, 1, 2]
        shortfall = LevelSpec([1, 1, 5], "gray").truncation(op)
        assert shortfall == [2, 0, 0]
        assert sum(shortfall) != 0

    def test_nothing_is_truncated_when_every_cutoff_covers_its_bound(self):
        shortfall = LevelSpec([2, 2, 2], "unary_projected").truncation(hamiltonian())
        assert shortfall == [0, 0, 0]
        assert sum(shortfall) == 0

    def test_a_uniform_spec_is_compared_against_every_mode_of_the_operator(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0, (cre(1), cre(1), cre(1)): 1.0})
        assert LevelSpec(1, "unary_projected").truncation(op) == [1, 2]

    def test_untouched_trailing_modes_have_no_shortfall(self):
        op = BosonOperator.from_dict({(cre(0), cre(0)): 1.0})
        assert LevelSpec([1, 0, 0], "unary_projected").truncation(op) == [1, 0, 0]

    def test_an_operator_wider_than_the_spec_raises(self):
        with pytest.raises(ValueError, match="covers only 1"):
            LevelSpec([5], "unary_projected").truncation(hamiltonian())

    def test_a_hard_core_model_can_report_a_shortfall_although_it_is_exact(self):
        """The report measures the lower-bound heuristic, not whether the mapping is accurate."""
        hard_core = BosonOperator.from_dict({(cre(0), cre(0), ann(0), ann(0)): 1.0})
        spec = LevelSpec([1], "unary_projected")
        assert spec.truncation(hard_core) == [1]
        # b†b†bb is exactly zero on the two retained levels, so the truncated mapping is exact
        assert len(occupation_encoding(hard_core, spec).simplify()) == 0


class TestOperatorModeCounts(unittest.TestCase):
    def test_min_num_modes(self):
        assert BosonOperator.zero().min_num_modes() == 0
        assert BosonOperator.one().min_num_modes() == 0
        assert BosonOperator.from_dict({(cre(2), ann(5)): 1.0}).min_num_modes() == 6

    def test_mode_occupations_default_length_and_zeros(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(3)): 1.0, (cre(3),) * 3: 1.0})
        assert op.mode_occupations() == [2, 0, 0, 3]
        assert BosonOperator.zero().mode_occupations() == []

    def test_mode_occupations_query(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(3)): 1.0, (cre(3),) * 3: 1.0})
        assert op.mode_occupations([3, 1, 0]) == [3, 0, 2]
        assert op.mode_occupations([9]) == [0]
        assert op.mode_occupations([]) == []

    def test_max_occupation_is_the_maximum_over_modes(self):
        op = BosonOperator.from_dict({(cre(0), cre(0), ann(3)): 1.0, (cre(3),) * 3: 1.0})
        assert op.max_occupation() == max(op.mode_occupations()) == 3
        assert BosonOperator.zero().max_occupation() == 0


class TestTranspiler(unittest.TestCase):
    def _circuit(self, occupation, num_modes=3):
        circuit = BosonicCircuit(num_modes)
        circuit.append(InitializeModes(occupation), circuit.modes)
        return circuit

    def test_the_layout_register_is_the_total_width(self):
        spec = LevelSpec([4, 2, 1], "gray")
        lowered = generate_preset_occupation_pass_manager(spec).run(self._circuit([1, 1, 0]))
        assert lowered.num_qubits == spec.total_width == 3 + 2 + 1

    def test_initialize_modes_lands_each_mode_at_its_own_offset(self):
        for encoding in ENCODINGS:
            spec = LevelSpec([3, 1, 2], encoding)
            occupation = [3, 1, 2]
            lowered = generate_preset_occupation_pass_manager(spec).run(self._circuit(occupation))
            index = sum(
                spec.encoding.level_to_bits(level) << offset
                for level, offset in zip(occupation, spec.offsets, strict=True)
            )
            amplitudes = np.asarray(Statevector(lowered).data)
            assert abs(amplitudes[index]) == pytest.approx(1.0), encoding

    def test_initialize_modes_checks_each_mode_against_its_own_cutoff(self):
        spec = LevelSpec([3, 1, 3], "unary_projected")
        pass_manager = generate_preset_occupation_pass_manager(spec)
        pass_manager.run(self._circuit([3, 1, 3]))  # exactly at every cutoff
        with pytest.raises(ValueError, match=r"occupation of 2 for mode 1.*max_occupation=1"):
            pass_manager.run(self._circuit([0, 2, 0]))

    def test_a_per_mode_spec_only_supports_that_many_modes(self):
        pass_manager = generate_preset_occupation_pass_manager(LevelSpec([2, 2], "unary_projected"))
        with pytest.raises(ValueError, match="covers 2 modes"):
            pass_manager.run(self._circuit([0, 0, 0]))

    def test_a_single_cutoff_serves_circuits_of_any_width(self):
        pass_manager = generate_preset_occupation_pass_manager(LevelSpec(2, "unary_projected"))
        widths = [pass_manager.run(self._circuit([0] * n, n)).num_qubits for n in (1, 3, 5)]
        assert widths == [3, 9, 15]

    def test_a_width_zero_mode_runs_through_the_whole_pipeline(self):
        """A frozen mode between two live ones: layout, state preparation and evolution."""
        spec = LevelSpec([1, 0, 1], "gray")
        hop = BosonOperator.from_dict({(cre(0), ann(2)): 1.0, (cre(2), ann(0)): 1.0})
        circuit = BosonicCircuit(3)
        circuit.append(InitializeModes([1, 0, 0]), circuit.modes)
        circuit.append(Evolution(3, hop, time=0.3), circuit.modes)

        lowered = preset_without_qubit_stage(spec).run(circuit)
        assert lowered.num_qubits == 2  # the frozen mode takes no qubit
        assert "PauliEvolution" in lowered.count_ops()

        # the evolved state is the exact rotation between the two live modes
        state = np.asarray(Statevector.from_int(0, 4).evolve(lowered.decompose(reps=4)).data)
        # |1,0,0> -> cos(t)|1,0,0> - i sin(t)|0,0,1>; mode 0 is qubit 0 and mode 2 is qubit 1
        assert abs(state[0b01]) == pytest.approx(np.cos(0.3), abs=1e-6)
        assert abs(state[0b10]) == pytest.approx(np.sin(0.3), abs=1e-6)

    def test_an_operator_living_only_on_a_width_zero_mode_maps_to_nothing(self):
        spec = LevelSpec([1, 0, 1], "gray")
        only_frozen = BosonOperator.from_dict({(cre(1), ann(1)): 1.0})
        assert len(occupation_encoding(only_frozen, spec)) == 0

    def test_relabel_modes_with_different_cutoffs_raises(self):
        """The cutoffs belong to positions, so reordering the modes would misassign them."""

        def run(spec):
            return MultiStagePassManager(
                input=BosonicCircuitToDAG(),
                optimization=RelabelModes(permutation=[1, 0, 2]),
                layout=generate_preset_occupation_pass_manager(spec).layout,
                output=QuantumDAGToCircuit(),
            ).run(self._circuit([0, 0, 0]))

        with pytest.raises(ValueError, match="RelabelModes"):
            run(LevelSpec([1, 2, 1], "unary_projected"))
        # a repeated cutoff is still a per-mode spec, so it is refused as well
        with pytest.raises(ValueError, match="RelabelModes"):
            run(LevelSpec([2, 2, 2], "unary_projected"))
        # a single shared cutoff cannot be misassigned, so it is unaffected; the layout stage
        # completes (there is no synthesis stage here, so the register is left as it is)
        assert run(LevelSpec(2, "unary_projected")).num_qubits == 3


UNARY_VARIANTS = ("unary_projected", "unary_unprojected")


def zero_projector_count(observable) -> int:
    """The number of ``|0><0|`` bit terms in a ``SparseObservable``, read from its term labels."""
    return sum(label.count("0") for label, _, _ in observable.to_sparse_list())


class TestUnaryVariants(unittest.TestCase):
    """The two unary encodings share their layout and differ in the operator built from it."""

    def test_the_unprojected_operator_equals_the_projected_one_on_the_one_hot_codewords(self):
        """Away from the one-hot subspace the two may differ, so only the codewords are compared."""
        op = hamiltonian()
        for max_occupation in (1, 2):
            projected = LevelSpec(max_occupation, "unary_projected").resolve(3)
            unprojected = LevelSpec(max_occupation, "unary_unprojected").resolve(3)
            indices = physical_indices(projected, projected.encoding)
            assert indices == physical_indices(unprojected, unprojected.encoding)

            matrices = []
            for spec in (projected, unprojected):
                observable = occupation_encoding(op, spec)
                matrices.append(SparsePauliOp.from_sparse_observable(observable).to_matrix())
            window = np.ix_(indices, indices)
            np.testing.assert_allclose(matrices[0][window], matrices[1][window], atol=1e-12)

    def test_both_unary_variants_reproduce_the_dense_reference(self):
        op = hamiltonian()
        for encoding in UNARY_VARIANTS:
            spec = LevelSpec(1, encoding).resolve(3)
            observable = occupation_encoding(op, spec)
            matrix = SparsePauliOp.from_sparse_observable(observable).to_matrix()
            indices = physical_indices(spec, spec.encoding)
            np.testing.assert_allclose(
                matrix[np.ix_(indices, indices)],
                dense_reference(op, [1, 1, 1]),
                atol=1e-12,
                err_msg=encoding,
            )

    def test_the_unprojected_operator_has_no_zero_projectors(self):
        op = hamiltonian()
        for max_occupation in (1, 2, 3):
            projected = occupation_encoding(op, LevelSpec(max_occupation, "unary_projected"))
            unprojected = occupation_encoding(op, LevelSpec(max_occupation, "unary_unprojected"))
            assert zero_projector_count(projected) > 0
            assert zero_projector_count(unprojected) == 0

    def test_the_unprojected_terms_do_not_grow_with_the_untouched_qubits(self):
        """A ladder term has two flipped qubits, so it touches two qubits whatever ``d`` is."""
        for max_occupation in (2, 3, 5, 8):
            observable = occupation_encoding(
                BosonOperator.from_dict({(cre(0),): 1.0}),
                LevelSpec(max_occupation, "unary_unprojected"),
            )
            widths = {len(indices) for _, indices, _ in observable.to_sparse_list()}
            assert widths == {2}

    def test_the_projected_terms_touch_every_qubit_of_the_mode(self):
        for max_occupation in (2, 3, 5):
            observable = occupation_encoding(
                BosonOperator.from_dict({(cre(0),): 1.0}),
                LevelSpec(max_occupation, "unary_projected"),
            )
            widths = {len(indices) for _, indices, _ in observable.to_sparse_list()}
            assert widths == {max_occupation + 1}

    def test_the_compact_encodings_are_unaffected_by_the_split(self):
        """Their projectors are part of the level number, so no variant drops them."""
        op = hamiltonian()
        for encoding in ("standard_binary", "gray"):
            assert zero_projector_count(occupation_encoding(op, LevelSpec(3, encoding))) > 0

    def test_both_unary_variants_share_their_layout(self):
        for max_occupation in (0, 1, 4):
            projected = LevelSpec(max_occupation, "unary_projected").resolve(3)
            unprojected = LevelSpec(max_occupation, "unary_unprojected").resolve(3)
            assert projected.widths == unprojected.widths
            assert projected.offsets == unprojected.offsets
            assert projected.total_width == unprojected.total_width

    def test_both_unary_variants_are_limited_to_128_levels(self):
        for encoding in UNARY_VARIANTS:
            assert LevelSpec(127, encoding).resolve(1).widths == [128]
            with pytest.raises(ValueError, match="unary encodings"):
                LevelSpec(128, encoding)
            with pytest.raises(ValueError, match="unary encodings"):
                LevelSpec([1, 128], encoding)
            with pytest.raises(ValueError, match="holds at most"):
                getattr(LevelEncoding, encoding.upper()).level_to_bits(128)

    def test_the_former_unary_name_is_rejected(self):
        for call in (
            lambda: LevelSpec(2, "unary"),
            lambda: LevelEncoding.from_name("unary"),
        ):
            with pytest.raises(ValueError, match="unary_projected"):
                call()

    def test_the_unary_variants_round_trip_through_name_pickle_and_repr(self):
        for variant, name in (
            (LevelEncoding.UNARY_PROJECTED, "unary_projected"),
            (LevelEncoding.UNARY_UNPROJECTED, "unary_unprojected"),
        ):
            assert str(variant) == name
            assert LevelEncoding.from_name(name) == variant
            assert repr(variant) == f"LevelEncoding.{name.upper()}"
            assert pickle.loads(pickle.dumps(variant)) == variant
            assert hash(variant) == hash(LevelEncoding.from_name(name))
        assert LevelEncoding.UNARY_PROJECTED != LevelEncoding.UNARY_UNPROJECTED
        assert LevelSpec(2, "unary_projected") != LevelSpec(2, "unary_unprojected")
        spec = LevelSpec([1, 2], "unary_unprojected")
        assert pickle.loads(pickle.dumps(spec)) == spec
        assert repr(spec) == "LevelSpec([1, 2], LevelEncoding.UNARY_UNPROJECTED)"

    def test_both_variants_prepare_the_same_state(self):
        """State preparation only needs the bit patterns, which the variants share."""
        for max_occupation in (1, 2):
            patterns = [
                [getattr(LevelEncoding, name.upper()).level_to_bits(n) for n in range(5)]
                for name in UNARY_VARIANTS
            ]
            assert patterns[0] == patterns[1]
            del max_occupation
