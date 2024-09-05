# Development guide
The project aims to use modern tools for streamlining code development, improve code quality and align code style between developers.

## General tools

### pre-commit hooks
In order to speed up development it is useful to receive feedback early in the iteration cycle. [`pre-commit`](https://pre-commit.com/) is a framework that allows for automating this process.
```bash
pip install pre-commit
```

`pre-commit` can be run manually on the changed files using
```bash
pre-commit run
```
or on all files in the repo using
```bash
pre-commit run --all-files
```

The commands being invoked are defined in `.pre-commit-config.yaml`. These commands can also be installed as git pre-commit hooks through
```bash
pre-commit install
```
which is highly recommended.


## Python development tools
> [!NOTE]
> Formatting rules are not currently enforced, but might be at an unspecified time in the future.

### Formatting and linting
The traditional way is to use [`black`](https://github.com/psf/black) for autoformatting, coupled with [`isort`](https://pycqa.github.io/isort/) for automatically sorting imports, and then use one or several other tools for linting. A more modern approach is to use [`ruff`](https://github.com/astral-sh/ruff) for both formatting and linting, as it can largely replace previous tools.

`ruff` is configured in `pyproject.toml` under `[tool.ruff]`.

To make `ruff` ignore formatting rules, one can use the pragma commands `# fmt: on`, `# fmt: off`, and `# fmt: skip`. It also respects `isort` commands (`# isort: skip_file`, `# isort: on`, `# isort: off`, etc) and to suppress lint errors, there is `# noqa`.

#### Automation
For increased quality-of-life it is recommended to install automatic formatting in your IDE of choice and turn on format-on-save. _E.g._, by installing the [`ruff` extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff) to VS Code.

In your `settings.json` add the following settings:
```json
"[python]": {
  "editor.formatOnSave": true,
  "editor.codeActionsOnSave": {
    "source.fixAll": "explicit",
    "source.organizeImports": "explicit"
  },
  "editor.defaultFormatter": "charliermarsh.ruff"
}
```

### Static type checking
Python is a dynamically typed language, but has since Python 3.5 introduced [type-hinting](https://peps.python.org/pep-0484/). Type-hints can be verified using the static type checker [`mypy`](https://github.com/python/mypy).

`mypy` is configured in `pyproject.toml` under `[tool.mypy]`.

To suppress a `mypy` error, one can use the inline pragma `# type: ignore[assignment]`, or `# mypy: disable-error-code=` at the top of the file.

### Unit testing
The two dominant testing frameworks are the built-in [`unittest`](https://docs.python.org/3/library/unittest.html) and the third-party [`pytest`](https://docs.pytest.org/), with the latter currently appearing more popular.

### Code coverage
One of the most popular way of measuring code coverage in python is by using [`Coverage.py`](https://coverage.readthedocs.io/en/latest/), installed through `pip install coverage`. Another tool often cited, `pytest-cov`, internally uses `Coverage.py` as a dependency.

Depending on which testing framework, `coverage run <test command>` runs the tests and saves the collected data in a binary file `.coverage`. A report can then be generated with `coverage report`, or `coverage html` for a nicer presentation.

### Profiling
There are a plethora of ways to profile python code, but a neat tool to use is [`pyinstrument`](https://github.com/joerick/pyinstrument), installed through `pip install pyinstrument`.

This profiler can be used in multiple ways, _e.g._ it lets you quickly profile code with:
```python
from pyinstrument import Profiler

profiler = Profiler()
profiler.start()

# code you want to profile

profiler.stop()
profiler.print()
```

## C++ development tools

### Formatting
Formatting of C++ code is done using [clang-format](https://clang.llvm.org/docs/ClangFormat.html) and the formatting rules are defined in `.clang-format`.

To suppress `clang-format` for a certain section of code, you can simply encapsulate it with
```cpp
// clang-format off
...
// clang-format on
```

#### Automation
As with python, automatic formatting in your IDE is recommended. In VS Code one can install the [`clang-format` extension](https://marketplace.visualstudio.com/items?itemName=xaver.clang-format)

In your `settings.json` simply add the following settings:
```json
"[cpp]": {
  "editor.formatOnSave": true,
  "editor.defaultFormatter": "xaver.clang-format"
}
```

### Linting
Linting in C++ is done using [clang-tidy](https://clang.llvm.org/extra/clang-tidy/), which is configured in `.clang-tidy`.

In order for `clang-tidy` to be able to interpret the project source code, it needs to know how it is built. For simple projects this information can be supplied with a `compile_flags.txt`, listing flags to be used for all files. More complicated projects instead uses `compile_commands.json`, providing a database over compilation commands.

To suppress `clang-tidy` warnings, one can use
`// NOLINT`, `// NOLINTNEXTLINE`, and `// NOLINTBEGIN`+`//NOLINTEND`.

### Static analysis
[`clangd`](https://clangd.llvm.org/) is a _language server_ that can provide smart features such as code completion, compile errors, go-to-definition etc. to your IDE. It also embeds `clang-tidy` (it also uses the `clang-format` engine for formatting), making it preferred over running `clang-tidy` separately. As with `clang-tidy`, `clangd` requires `compile_commands.json` to interpret the source code.

`clangd` recognizes `.clang-tidy` for `clang-tidy`-settings and can itself be configured with a `.clangd` file.

In VS Code one can install the
[`clangd` extension](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd). The extension itself doesn't come with `clangd`, which requires separate installation.

#### Compile commands
The process for generating `compile_commands.json` depends on which build system is being used, and can often be done using the build tool itself. See [official documentation](https://clang.llvm.org/docs/JSONCompilationDatabase.html).

### Code coverage
TODO

## Resources
- [Scientific Python Library Development Guide](https://learn.scientific-python.org/development/guides/style/)
