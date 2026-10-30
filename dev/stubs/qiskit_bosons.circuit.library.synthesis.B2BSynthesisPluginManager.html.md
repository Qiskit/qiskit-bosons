# B2BSynthesisPluginManager

### *class* B2BSynthesisPluginManager

A simple manager of the installed boson-to-boson synthesis plugins.

Bases: [`object`](https://docs.python.org/3/builtins/functions.html#object)

All plugins that are registered under the `qiskit_bosons.circuit.synthesis`
[entry-point](https://setuptools.pypa.io/en/latest/userguide/entry_point.html) will be
available to this class and get exposed by the [`plugins`](#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPluginManager.plugins) attribute. The name of an entry
point is the name of the gate class followed by the name of the method, separated by a dot, for
example `ModeRotation.Clements`. Below is the table from the `pyproject.toml` of this package
itself:

```toml
[project.entry-points."qiskit_bosons.circuit.synthesis"]
"Evolution.LieTrotter" = "qiskit_bosons.circuit.library.synthesis:BosonicLieTrotter"
"Evolution.SuzukiTrotter" = "qiskit_bosons.circuit.library.synthesis:BosonicSuzukiTrotter"
"Evolution.QDrift" = "qiskit_bosons.circuit.library.synthesis:BosonicQDrift"
"ModeRotation.Clements" = "qiskit_bosons.circuit.library.synthesis:BosonicClementsDecomposition"

```

An overview of all builtin plugins is given [here](../pydoc/qiskit_bosons.circuit.library.synthesis.html.md#builtin-b2b-synthesis-plugins).

### Attributes

#### plugins *: [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[str](https://docs.python.org/3/builtins/stdtypes.html#str), [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]]*

A dictionary mapping gate names to their synthesis methods.

The keys of this dictionary are names of [`BosonicGate`](qiskit_bosons.circuit.BosonicGate.html.md#qiskit_bosons.circuit.BosonicGate) classes. The values are the
list of installed [`B2BSynthesisPlugin`](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin) methods that may be used to synthesize the
instances of that particular gate.

### Methods

#### method(op_name, method_name)

Returns the requested [`B2BSynthesisPlugin`](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin) type.

* **Parameters:**
  * **op_name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – the name of the gate whose plugin to return.
  * **method_name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – the name of the synthesis plugin.
* **Returns:**
  The [`B2BSynthesisPlugin`](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin) class.
* **Raises:**
  [**KeyError**](https://docs.python.org/3/builtins/exceptions.html#KeyError) – if no plugin is registered under that name.
* **Return type:**
  [type](https://docs.python.org/3/builtins/functions.html#type)[[*B2BSynthesisPlugin*](qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin.html.md#qiskit_bosons.circuit.library.synthesis.B2BSynthesisPlugin)]

#### method_names(op_name)

Returns the names of all installed plugin methods for a given gate.

* **Parameters:**
  **op_name** ([*str*](https://docs.python.org/3/builtins/stdtypes.html#str)) – the name of the gate whose synthesis methods to return.
* **Returns:**
  The list of installed plugin method names.
* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]

#### op_names()

Returns the names of all gates for which any plugin is installed.

* **Return type:**
  [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[str](https://docs.python.org/3/builtins/stdtypes.html#str)]
