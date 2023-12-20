# Making extensions

$\pi$-PIC offers a possibility to develop extensions in Python, C/C++, Fortran and other programming languages that can be used to produce a callable function for Python. Extensions can be developed, compiled and used independently, as well as contributed to the toolbox delivered with $\pi$-PIC (see [Python extensions](https://github.com/hi-chi/pipic/blob/main/pipic/extensions) and [C++ extensions](https://github.com/hi-chi/pipic/blob/main/src/extensions)). Apart from modifying the field state via `fieldLoop()`, the extensions can modify, add and remove particles based on the local field state. For example, this can be used to account for ionization, radiation reaction, QED particle and photon generation, etc.

Developing a new extension implies creating an `<extension>.so` (C/C++ extensions) or `<extension>.py` (Python extensions) file, which we examplify below. After placing this file in a folder with the file of simulation, one can add the handler of an extension to be called by a given $\pi$-PIC container (replace `<*>` with respective entries):
```
import <extension> 
<simulation>.add_cell_handler(name=<extension>.name, subject=’<particleType>’, \
handler=<extension>.handler(<handler parameters>)
```
Here `<particleType>` is a list of particle types to be processed by the extension; `all_types` can be used to process all types, `cells` to call the handler for each cell independently of particles there (can be used to add new particles). In addition, if needed one can pass addressed of the extension's data:
```
data_double=<extension>.data_double(), data_int=<extension>.data_int()
```

Python extensions
--
Creating an extension within Python is based on developing a handler in the form of a [Numba C callback](https://numba.readthedocs.io/en/stable/user/cfunc.html) that is called for each cell and offers a possibility to add new particles and/or modify states of all particles in the cell sequentially for all requested types. When called the handler receives data pointers to be passed to a `CellInterface` structure that offers basic information and options for modifying existing and/or adding new particles within the processed cell (see complete list of interfaces in [`cellinterface.py`](https://github.com/hi-chi/pipic/blob/main/pipic/interfaces/cellinterface.py)). The state of the field can be modified whenever needed using `field_loop()`, see [User interfaces](https://github.com/hi-chi/pipic/blob/main/docs/guides/INTERFACES.md). Despite some limitations on what can be used within Numba C callbacks they can provide sufficiently high performance and flexibility for many cases of interest.

The layout of a Python extension includes the following elements:
1.  Set the name of extension for further reference (should be the same as the name of the extension file)
2.  If needed, allocate data of real and integer type to be accessed from the handler (the access must be thread-safe)
3.  Develop the handler, which should start from initiating the `CellInterface`:
    ```
    @cfunc(handler_callback)
    def Handler(CI_I, CI_D, CI_F, CI_P, CI_NP, data_double, data_int):
        C = CellInterface(CI_I, CI_D, CI_F, CI_P, CI_NP) # unpacking data
        ...
    ```
4.  Develop a function that initiates the exension (i.e. set the internal data based on any parameters needed) and returns the address to the `Handler`:
    ```
    def handler(<parameters>):
        ...
        return Handler.address
    ```
5.  Develop a function that returns the address of the extension's data (if needed), e.g.:
    ```
    def data_double():
        return addressof(DataDouble) 
    ```


As an example see an extension [`x_reflector_py.py`](https://github.com/hi-chi/pipic/blob/main/pipic/extensions/x_reflector_py.py) that can be used to reflect particles from an $x$-limited region and its use in [`x_reflector_py_test.py`](https://github.com/hi-chi/pipic/blob/main/examples/x_reflector_py_test.py) (note that this extension is added to $\pi$-PIC as an example).

C/C++ extensions
--
Creating an extension in C/C++ provides full flexibility and performance, while also permitting the creation of small extensions compiled locally. The layout of development is similar to the one for Python developments:
1.	From `src` folder copy the following interface files in a local folder:
    - `primitives.h`
    - `interfaces.h`
    - `CMakeLists.txt`
2.	In the file `CMakeLists.txt` replace “pipic” with the extension name.
3.	Clone `pybind11` to the local folder"
    ```
    git clone https://github.com/pybind/pybind11
    ```
4.	Create an `<extension>.cpp` file and develop the code of extension following the layout of examples in [`extensions`](https://github.com/hi-chi/pipic/blob/main/src/extensions) (for more controls see methods of `struct cellInterface` in [`interfaces.h`](https://github.com/hi-chi/pipic/blob/main/src/interfaces.h)).
5.	Compile your extension by sequentially running:
    ```
    cmake .
    make
    ```

As an example see an extension [`x_converter_c.cpp`](https://github.com/hi-chi/pipic/blob/main/src/extensions/x_converter_c/x_converter_c.cpp) that can be used to convert type of particles in an $x$-limited region. The use of this extension is exemplified in [`x_reflector_py_test.py`](https://github.com/hi-chi/pipic/blob/main/examples/x_converter_c_test.py) (note that this extension is added to $\pi$-PIC as an example).

