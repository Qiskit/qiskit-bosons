<a id="module-qiskit_bosons.protocols"></a>

# Protocols

This module provides the protocols used throughout the package.

A [`typing.Protocol`](https://docs.python.org/3/library/typing.html#typing.Protocol) is an interface that a class can fulfill structurally, without explicitly
inheriting from it. For example, any class that implements a `__len__` method satisfies
[`Sized`](https://docs.python.org/3/library/typing.html#typing.Sized), whether or not it says so anywhere in its class hierarchy.

The protocols here follow the same idea, each centered on a small set of dunder-style methods. A
class implementing them thereby satisfies the protocol and becomes usable wherever it is expected,
without any explicit registration.

| Protocol                                                                                                         | Methods                                                    | Helper functions                                                                                                                                                                                                                                                                                                                                     |
|------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`SupportsCommutators`](../stubs/qiskit_bosons.protocols.SupportsCommutators.html.md#qiskit_bosons.protocols.SupportsCommutators) | `_commutator_`, `_anti_commutator_`, `_double_commutator_` | [`commutator()`](../stubs/qiskit_bosons.operators.library.commutator.html.md#qiskit_bosons.operators.library.commutator), [`anti_commutator()`](../stubs/qiskit_bosons.operators.library.anti_commutator.html.md#qiskit_bosons.operators.library.anti_commutator), [`double_commutator()`](../stubs/qiskit_bosons.operators.library.double_commutator.html.md#qiskit_bosons.operators.library.double_commutator) |
