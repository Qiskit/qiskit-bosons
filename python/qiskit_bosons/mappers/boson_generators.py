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

"""BosonOperator mapper."""

from collections.abc import Callable
from operator import and_

from qiskit_bosons.mappers._typing import T
from qiskit_bosons.operators import BosonOperator
from qiskit_bosons.operators.boson_action import BosonAction


def map_boson_action_generators(
    operator: BosonOperator,
    map_action: Callable[[BosonAction], T],
    identity: Callable[[], T],
    compose: Callable[[T, T], T] | None = None,
) -> T:
    r"""Map a :class:`.BosonOperator` to another operator type.

    This is a generic function to aid in implementing new mappers for :class:`.BosonOperator`
    instances. At its core, it simply iterates over the terms of the operator, mapping each
    encountered :class:`.BosonAction` with the user-provided ``map_action`` function. In
    combination with the user-provided ``identity`` generator, this allows mapping to arbitrary
    output types.

    .. note::
       The output type ``T`` must support multiplication by a scalar via ``__mul__``.
       If ``compose=None`` it must also support composition of two instances via ``__and__``.
       :class:`~qiskit.quantum_info.SparseObservable` does not, so pass its
       :meth:`~qiskit.quantum_info.SparseObservable.compose` explicitly as ``compose``; a type with
       an ``__and__`` (such as :class:`~qiskit.quantum_info.SparsePauliOp`) can rely on the default.

    .. note::
       This function is an extension point for your own mappers. It is a minimal illustration and
       not a replacement for :func:`~qiskit_bosons.mappers.library.occupation_encoding`, which is
       implemented in Rust, is parallelized, and bounds the memory it uses while assembling the
       result. Nothing in the library calls this function.

    The example maps a mode with two levels (``max_occupation = 1``) to one qubit, where the
    annihilation operator is :math:`\lvert 0 \rangle\langle 1 \rvert = (X + iY) / 2` and the
    creation operator is its adjoint. Modes with more levels need the encoding of
    :func:`~qiskit_bosons.mappers.library.occupation_encoding` instead.

    .. doctest::

        >>> from qiskit_bosons.mappers import map_boson_action_generators
        >>> from qiskit_bosons.operators import BosonAction, BosonOperator, ann, cre
        >>> from qiskit.quantum_info import SparsePauliOp
        >>>
        >>> num_qubits = 2
        >>> def map_action(action: BosonAction) -> SparsePauliOp:
        ...     is_creation, mode = action
        ...     return SparsePauliOp.from_sparse_list(
        ...         [("X", [mode], 0.5), ("Y", [mode], -0.5j if is_creation else 0.5j)],
        ...         num_qubits=num_qubits,
        ...     )
        >>>
        >>> def identity() -> SparsePauliOp:
        ...     return SparsePauliOp.from_sparse_list([("", [], 1)], num_qubits)
        >>>
        >>> op = BosonOperator.from_dict({(cre(0), ann(1)): 2.0})
        >>> qop = map_boson_action_generators(op, map_action, identity).simplify()
        >>> print([(label, complex(coeff)) for label, coeff in sorted(qop.label_iter())])
        [('XX', (0.5+0j)), ('XY', -0.5j), ('YX', 0.5j), ('YY', (0.5+0j))]

    Args:
        operator: the operator to be mapped.
        map_action: the function to map a single :class:`.BosonAction` to the desired output type.
        identity: the function to generate the multiplicative identity instance of the output type.
        compose: an optional function to implement the composition logic of two output type
            instances. If this is not provided, it will default to using :py:func:`operator.and_`.

    Returns:
        The mapped operator.
    """
    if compose is None:
        compose = and_

    mapped_operator: T = 0 * identity()
    for terms, coeff in operator.iter_terms():
        mapped_terms = identity()

        for term in terms:
            mapped_terms = compose(map_action(term), mapped_terms)

        mapped_operator += coeff * mapped_terms

    return mapped_operator
