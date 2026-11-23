# B2QSynthesisPluginManager

### *class* B2QSynthesisPluginManager

A simple manager of the installed boson-to-qubit synthesis plugins.

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

All plugins that are registered under the `qiskit_bosons.transpiler.synthesis`
[entry-point](https://setuptools.pypa.io/en/latest/userguide/entry_point.html) will be
available to this class and get exposed by the [`plugins`](#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPluginManager.plugins) attribute. Below is an example
from the `pyproject.toml` of this package itself:

```toml
[project.entry-points."qiskit_bosons.transpiler.synthesis"]
"Evolution.MapperFn" = "qiskit_bosons.transpiler.passes.synthesis:MapperFnEvolutionSynthesis"
"InitializeModes.Occupation" = "qiskit_bosons.transpiler.passes.synthesis:OccupationInitializeModesSynthesis"

```

An overview of all builtin plugins is given [here](../pydoc/qiskit_bosons.transpiler.passes.synthesis.html.md#builtin-b2q-synthesis-plugins).

### Attributes

#### plugins *: [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[str](https://docs.python.org/3/builtins/stdtypes.html#str), [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]]*

A dictionary mapping circuit instruction names to their synthesis methods.

The keys of this dictionary are names of [`BosonicGate`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) classes. The values are the
list of installed [`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) methods that may be used to synthesize the
instances of that particular gate to its qubit form.

### Methods

#### method(op_name, method_name)

Returns the requested [`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) type.

* **Parameters:**
  * **op_name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – the name of the operation whose plugin to return.
  * **method_name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – the name of the synthesis plugin.
* **Returns:**
  The [`B2QSynthesisPlugin`](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin) class.
* **Return type:**
  [type](https://docs.python.org/3/builtins/functions.html#type)[[*B2QSynthesisPlugin*](qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin.html.md#qiskit_bosons.transpiler.passes.synthesis.B2QSynthesisPlugin)]

#### method_names(op_name)

Returns the names of all installed plugin methods for a given circuit instruction.

* **Parameters:**
  **op_name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – the name of the operation whose synthesis methods to return.
* **Returns:**
  The list of installed plugin method names.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]

#### op_names()

Returns the names of all circuit instructions for which any plugin is installed.

* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]
