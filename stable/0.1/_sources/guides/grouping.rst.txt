.. _grouping_explanation:

Group operator terms: use the operator structure
================================================

.. meta::
   :description: How to assign group indices to the terms of a BosonOperator, check that each group is
                 Hermitian or has uniform coefficients, and lay out the operator so that picking a few
                 groups stays cheap.

A :class:`~qiskit_bosons.operators.BosonOperator` can tag each of its terms with a **group
index**. The tags carry no meaning of their own: they only say which terms belong together.
What the tags *mean* is decided by whatever consumes them. In this package the consumers are
the time-evolution gate, which is split into one factor per group, and the randomized qDRIFT
pass, which samples whole groups.

This guide shows how to assign groups, how to check the properties a consumer relies on, and
how to lay the operator out so that picking a few groups out of a large one stays cheap.

.. note::
   No automatic grouping strategy ships with this package. You assign the indices yourself, from
   what you know about your model.

Assign groups
-------------

The running example is a three-site Bose-Hubbard chain. Each bond contributes a hopping pair, and
the on-site interactions form a diagonal piece. The simplest way to build a grouped operator is
:meth:`~qiskit_bosons.operators.BosonOperator.from_terms_with_groups`, which takes each term
together with its group index as an ``(actions, coefficient, group)`` triple. You decide the group
at the moment you write the term down, when you still know why it belongs there:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from qiskit_bosons.operators import BosonOperator, ann, cre
   >>>
   >>> bonds = [(0, 1), (1, 2)]
   >>> terms = []
   >>> for group, (left, right) in enumerate(bonds):  # one group per bond
   ...     terms.append(([cre(left), ann(right)], -1.0, group))
   ...     terms.append(([cre(right), ann(left)], -1.0, group))
   >>> for site in range(3):  # the on-site interactions share one group
   ...     terms.append(([cre(site), cre(site), ann(site), ann(site)], 1.0, len(bonds)))
   >>> hamiltonian = BosonOperator.from_terms_with_groups(terms)
   >>> len(hamiltonian), hamiltonian.num_groups()
   (7, 3)

The terms are stored in the order you give them, so the group indices cannot end up attached to the
wrong terms.

An operator that already exists has to be grouped after the fact, by assigning the
:attr:`~qiskit_bosons.operators.BosonOperator.groups` attribute. That is where care is needed,
because :meth:`~qiskit_bosons.operators.BosonOperator.from_dict` does not specify the order in which
it stores its terms. Writing the list out by position can then attach an index to the wrong term.
Derive each index from the term itself instead:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> existing = BosonOperator.from_dict({tuple(actions): coeff for actions, coeff, _ in terms})
   >>> existing.groups = [
   ...     len(bonds)
   ...     if len({mode for _, mode in term}) == 1
   ...     else bonds.index(tuple(sorted({mode for _, mode in term})))
   ...     for term, _ in existing.iter_terms()
   ... ]
   >>> existing.num_groups(), existing.equiv(hamiltonian)
   (3, True)

The number of indices has to match the number of terms. Anything else is rejected, because a
too-short list would silently drop trailing terms wherever terms are iterated together with
their groups:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> try:
   ...     hamiltonian.groups = [0, 1]
   ... except ValueError as error:
   ...     print(error)
   expected one group index per term, but got 2 for 7 terms

Check what the consumer relies on
---------------------------------

Three functions in :mod:`qiskit_bosons.operators.terms` answer one narrow question each
about an existing grouping. None of them says a grouping is "correct", because correctness is
relative to the consumer.

**Is each group Hermitian?** This is the property that matters most. The time evolution of a
group is only unitary if the group is Hermitian, and a lone hopping term :math:`b^\dagger_i b_j`
is *not*. Only the pair with its adjoint is Hermitian:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from qiskit_bosons.operators.terms import groups_are_hermitian
   >>> groups_are_hermitian(hamiltonian)
   [True, True, True]

A grouping that splits a hopping pair is caught immediately:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> split = BosonOperator.from_terms_with_groups(
   ...     [(actions, coeff, index) for index, (actions, coeff, _) in enumerate(terms)]
   ... )  # one term per group
   >>> [bool(flag) for flag in groups_are_hermitian(split)]
   [False, False, False, False, True, True, True]

The four hopping terms are caught. The three ``True`` entries are the on-site terms
:math:`b^\dagger_i b^\dagger_i b_i b_i`, which are Hermitian by themselves. A mode can
legitimately repeat within a term, and nothing is dropped for doing so. The flags follow the
order of ``terms``, which is why each one can be read against the term it describes.

**Are the coefficients uniform?** A randomized product formula weighs a group by the mean
*magnitude* of its coefficients, which is a faithful description of the group only if the
magnitudes agree:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from qiskit_bosons.operators.terms import group_coeff_means, groups_have_uniform_coeffs
   >>> groups_have_uniform_coeffs(hamiltonian)
   [True, True, True]
   >>> group_coeff_means(hamiltonian)
   [1.0, 1.0, 1.0]

Neither check implies the other: uniform coefficients do not make a group Hermitian, and a
Hermitian group can mix magnitudes. By default the magnitudes are compared, because that is what
:func:`~qiskit_bosons.operators.terms.group_coeff_means` averages. A Hermitian pair with complex
coefficients has equal magnitudes but unequal coefficients, so it satisfies only the default form:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> conjugate = BosonOperator.from_terms_with_groups(
   ...     [([cre(0), ann(1)], 1.0j, 0), ([cre(1), ann(0)], -1.0j, 0)]
   ... )
   >>> groups_have_uniform_coeffs(conjugate)
   [True]
   >>> groups_have_uniform_coeffs(conjugate, abs=False)
   [False]

.. note::
   :func:`~qiskit_bosons.operators.terms.groups_are_hermitian` is built on
   :meth:`~qiskit_bosons.operators.BosonOperator.is_hermitian`, which compares an operator with its
   adjoint after normal ordering. That is an identity of the untruncated algebra, so a group that
   passes is Hermitian as a bosonic operator. See :ref:`operators_explanation` for how products of
   ladder operators are ordered and compared.

Order the terms by group
------------------------

:meth:`~qiskit_bosons.operators.BosonOperator.split_out_groups` returns one operator per
requested group. When the terms of a group are scattered through the operator, it has to scan
every term to find them. When each group forms one contiguous run, it binary-searches the group
boundaries instead, so the cost of a lookup follows the number of groups *requested* rather
than the number of terms *held*.

:func:`~qiskit_bosons.operators.terms.group_order` produces that layout. It is a stable sort by
group index: the value of the operator does not change, and terms within a group keep their
relative order.

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from qiskit_bosons.operators.terms import group_order
   >>> scattered = BosonOperator.from_terms_with_groups(
   ...     [(actions, coeff, group) for (actions, coeff, _), group in zip(terms, [1, 0, 1, 0, 2, 2, 2])]
   ... )
   >>> ordered = group_order(scattered)
   >>> ordered.groups
   [0, 0, 1, 1, 2, 2, 2]
   >>> ordered.equiv(scattered)
   True

Both layouts return the same operators; only the cost differs:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> request = [2, 0]
   >>> all(
   ...     a.equiv(b)
   ...     for a, b in zip(
   ...         scattered.split_out_groups(request), ordered.split_out_groups(request), strict=True
   ...     )
   ... )
   True

Next steps
----------

The grouping is chosen here, at the bosonic level, while the structure of the model is still
visible.

- :ref:`bosonic_synthesis_explanation` shows how to use it to Trotterize an evolution, and why
  that reduces an error which a qubit-level product formula cannot.
- :ref:`operators_explanation` explains the operator data structure that carries the group
  indices, and how its terms are stored and ordered.
- :ref:`bose_hubbard_grouping` applies a grouping to a complete model.
