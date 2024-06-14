# Development guide

The project aims to use modern tools for streamlining code development, improve code quality and align code style between developers.

## Python development tools

### Formatting
> [!NOTE]
> Formatting rules are not currently enforced, but might be at an unspecified time in the future.

For autoformatting python-code the project uses [`black`](https://github.com/psf/black), _'the uncompromising Python code formatter'_, which can be installed with
```
pip install black
```

The formatter can be run manually on the python directory using
```
black --config pyproject.toml pipic/
```
Config settings are stored in `pyproject.toml` under `[tool.black]`.

Upon occasion, manual formatting is desired. This can be achived by temporarily turning formatting off: 
```
# fmt off
Manually formatted python code goes here
# fmt on
```

#### Automation
For increased quality-of-life it is recommended to use black for automatic formatting in your IDE of choice. _E.g._, the black formatter extension can be installed on VS Code. Simply turn on format-on-save, and add `--line-length=100` as argument.


### Linting
[`ruff`](https://github.com/astral-sh/ruff)

### Static type checking
[`mypy`](https://github.com/python/mypy)


## C++ development tools
clang-tidy