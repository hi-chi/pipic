# Development guide

The project aims to use modern tools for streamlining code development, improve code quality and align code style between developers.


## Python development tools
> [!NOTE]
> Formatting rules are not currently enforced, but might be at an unspecified time in the future.


### pre-commit hooks
In order to speed up development it is useful to have early feedback in the iteration cycle. [`pre-commit`](https://pre-commit.com/) is a framework that allows for automating this process.
```
pip install pre-commit
```

`pre-commit` can be run manually on the changed files using
```
pre-commit run
```
or on all files in the repo using
```
pre-commit run --all-files
```

The commands to run are defined in `.pre-commit-config.yaml`. These commands can also be installed as git pre-commit hooks through
```
pre-commit install
```
which is highly recommended.


### Formatting and linting

The traditional way was to use [`black`](https://github.com/psf/black) for autoformatting, coupled with [`isort`](https://pycqa.github.io/isort/) for automatically sorting imports, and then use another tool for linting. A more modern approach is to use [`ruff`](https://github.com/astral-sh/ruff) for both formatting and linting, as it can largely replace previous tools.

`ruff` is configured in `pyproject.toml` under `[tool.ruff]`.


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
[`mypy`](https://github.com/python/mypy)

`mypy` is configured in `pyproject.toml` under `[tool.mypy]`.


## Resources
- [Scientific Python Library Development Guide](https://learn.scientific-python.org/development/guides/style/)
