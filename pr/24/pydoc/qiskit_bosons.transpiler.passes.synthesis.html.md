<a id="module-qiskit_bosons.transpiler.passes.synthesis"></a>

<a id="qiskit-bosons-transpiler-passes-synthesis-plugins"></a>

# Transpiler Pass Synthesis (`B2Q`)

The [`B2QSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis) transpiler pass exposes a plugin interface (see
[`B2QSynthesisPlugin`](../stubs/qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin)) through which implementations for mapping boson-based instructions
to qubit-based ones can be registered. It is the counterpart of the boson-to-boson synthesis that the
[`B2BSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2BSynthesis.html.md#qiskit_bosons.transpiler.passes.B2BSynthesis) pass performs one stage earlier (see
[`synthesis`](qiskit_bosons.circuit.library.synthesis.html.md#module-qiskit_bosons.circuit.library.synthesis)).

The available plugins are managed by the [`B2QSynthesisPluginManager`](../stubs/qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPluginManager.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPluginManager). Choosing the right
combination of synthesis plugins is crucial and the [`B2QSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis) transpilation pass does
**not** have any default plugins defined! For this, refer to the transpiler
[`presets`](qiskit_bosons.transpiler.presets.html.md#module-qiskit_bosons.transpiler.presets).

<a id="builtin-b2q-synthesis-plugins"></a>

## Builtin plugins

The plugins to use at runtime have to be registered in the [`B2QSynthesis.methods`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis.methods) attribute.
This can be done easily through the `config` argument of [`B2QSynthesis`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis) which specifies
the plugins to use from the [`B2QSynthesisPluginManager`](../stubs/qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPluginManager.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPluginManager) or it can be done manually by
writing into [`B2QSynthesis.methods`](../stubs/qiskit_bosons.transpiler.passes.B2QSynthesis.html.md#qiskit_bosons.transpiler.passes.B2QSynthesis.methods) directly.

There is one plugin per [`BosonicGate`](../stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) type. Each is handed the gate and the global
[`B2QLayout`](../stubs/qiskit_bosons.transpiler.B2QLayout.html.md#qiskit_bosons.transpiler.B2QLayout), and inserts the equivalent qubit-based instructions.

The builtin plugins are where the boson-to-qubit **encoding** and the **truncation** enter the
pipeline: the plugins taking a `mapper_fn` do so through the mapper they wrap, and the
`InitializeModes` plugin through the [`LevelSpec`](../stubs/qiskit_bosons.mappers.LevelSpec.html.md#qiskit_bosons.mappers.LevelSpec) it is given. This is the point of
deferring the encoding to transpile time: the bosonic circuit upstream of this stage is independent
of both choices.

The table below lists the builtin plugins, organized by the [`BosonicGate`](../stubs/qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) that they act
upon. These are the entry points of the `qiskit_bosons.transpiler.synthesis` group, named
`<gate>.<method>`:

#### An overview of the builtin `qiskit_bosons.transpiler.synthesis` entry-points.

| Circuit instruction name                                                                                       | Plugin method name   | Reference                                                                                                                                                        |
|----------------------------------------------------------------------------------------------------------------|----------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`Evolution`](../stubs/qiskit_bosons.circuit.library.Evolution.html.md#qiskit_bosons.circuit.library.Evolution)             | `MapperFn`           | [`MapperFnEvolutionSynthesis`](../stubs/qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis.html.md#qiskit_bosons.transpiler.passes.synthesis.MapperFnEvolutionSynthesis)                 |
| [`InitializeModes`](../stubs/qiskit_bosons.circuit.library.InitializeModes.html.md#qiskit_bosons.circuit.library.InitializeModes) | `Occupation`         | [`OccupationInitializeModesSynthesis`](../stubs/qiskit_bosons.transpiler.passes.synthesis.OccupationInitializeModesSynthesis.html.md#qiskit_bosons.transpiler.passes.synthesis.OccupationInitializeModesSynthesis) |

## Additional helper functions

The Pauli term order that a mapper function produces is preserved through synthesis, since a product
formula synthesizes the terms in the order it receives them. These adapters wrap a mapper function to
choose that order:
