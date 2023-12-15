# Unit test guide
The following document provides a guide for unit testing $\pi$-PIC in order to maintain correctness during development.

## How to run 
In order to be able to run the tests, the project must first be built (because the binaries must be compiled). To do this, simply start by installing the local project
```
pip install .
```
after which unit tests can be run by writing
```
python -m unittest
```
while standing in the project root folder. This utilises the [unittest](https://docs.python.org/3/library/unittest.html) package, which is part of Python's standard library.

The above command will automatically run python tests from reachable subpackages (any directory with an `__init__.py`), such as `tests/`, with files matching the pattern `test*.py`.

## Creating unit tests
Tests are made by creating a class that inherits from `unittest.TestCase`. Each member function that matches the pattern `test*()` will be run, and is considered to be _one_ test. For more detailed guides, see links below.

## Links
For guides, see _e.g._:
- https://www.browserstack.com/guide/unit-testing-python