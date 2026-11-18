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

"""The equality and the equivalence of bosonic Evolution gates."""

from qiskit_bosons.circuit import BosonicCircuit
from qiskit_bosons.circuit.library import Evolution
from qiskit_bosons.operators import BosonOperator, ann, cre


def hopping(scale=1.0):
    return BosonOperator.from_dict({(cre(0), ann(1)): scale, (cre(1), ann(0)): scale})


def test_equality_compares_the_operator_the_time_and_the_modes():
    evo = Evolution(2, hopping(), time=0.5)

    # the synthesis and the atomic flag are not part of what the gate is
    assert evo == Evolution(2, hopping(), time=0.5, atomic=True)
    assert hash(evo) == hash(Evolution(2, hopping(), time=0.5))
    assert evo != Evolution(2, hopping(2.0), time=0.5)
    assert evo != Evolution(2, hopping(), time=0.25)
    assert evo != Evolution(3, hopping(), time=0.5)


def test_atomic_evolutions_of_different_operators_differ():
    # an atomic gate has no definition, so Qiskit's equality of instructions saw only the time
    evo = Evolution(2, hopping(), time=0.5, atomic=True)
    assert evo != Evolution(2, hopping(2.0), time=0.5, atomic=True)


def test_equality_does_not_depend_on_the_order_of_the_terms():
    forward = BosonOperator.from_dict({(cre(0), ann(1)): 1.0, (cre(1), ann(0)): 1.0})
    backward = BosonOperator.from_dict({(cre(1), ann(0)): 1.0, (cre(0), ann(1)): 1.0})
    assert Evolution(2, forward, time=0.5) == Evolution(2, backward, time=0.5)


def test_an_appended_evolution_equals_the_one_given():
    evo = Evolution(2, hopping(), time=0.5)
    circuit = BosonicCircuit(2)
    circuit.append(evo, circuit.modes)
    assert circuit._inner.data[0].operation == evo


def test_circuits_of_evolutions_of_different_operators_differ():
    first, second = BosonicCircuit(2), BosonicCircuit(2)
    first.append(Evolution(2, hopping(), time=0.5), first.modes)
    second.append(Evolution(2, hopping(2.0), time=0.5), second.modes)
    assert first._inner != second._inner


def test_equiv_compares_up_to_a_tolerance():
    evo = Evolution(2, hopping(), time=0.5)
    close = BosonOperator.from_dict({(cre(0), ann(1)): 1.0 + 1e-12, (cre(1), ann(0)): 1.0})

    assert evo.equiv(Evolution(2, close, time=0.5 + 1e-12))
    assert not evo.equiv(Evolution(2, close, time=0.5 + 1e-12), atol=1e-14)
    assert not evo.equiv(Evolution(2, hopping(), time=0.6))
    assert not evo.equiv(Evolution(2, hopping(2.0), time=0.5))


def test_the_operators_are_compared_in_the_algebra_of_the_bosons():
    # b b^dagger and b^dagger b + 1 are the same operator: the circuit lives in the untruncated
    # algebra, and a truncation only enters when it is lowered to qubits
    a = BosonOperator.from_dict({(ann(0), cre(0)): 1.0})
    b = BosonOperator.from_dict({(cre(0), ann(0)): 1.0, (): 1.0})
    assert Evolution(1, a, time=0.1).equiv(Evolution(1, b, time=0.1))
    assert Evolution(1, a, time=0.1) == Evolution(1, b, time=0.1)
