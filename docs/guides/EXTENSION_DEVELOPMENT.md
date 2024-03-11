# Making extensions

$\pi$-PIC offers a possibility to develop extensions in Python, C/C++, Fortran and other programming languages that can be used to produce a callable function for Python. Extensions can be developed, compiled and used independently, as well as contributed to the toolbox delivered with $\pi$-PIC (see [Python extensions](https://github.com/hi-chi/pipic/blob/main/pipic/extensions) and [C++ extensions](https://github.com/hi-chi/pipic/blob/main/src/extensions)). Apart from modifying the field state via `fieldLoop()`, the extensions can modify, add and remove particles based on the local field state. For example, this can be used to account for ionization, radiation reaction, QED particle and photon generation, etc.

The connection of an extension to $\pi$-PIC is established via Python, but utilises direct data access without loss of performance. This is arranged via an interface introduced by the struct cellInerface given in `/pipic/interfaces/cellInterface.py` for Python and `src/interfaces.h` for C/C++ (the interface can be extended or ported to other programming languages by request). The extensions must provide an address to a static function called *Handler* that takes raw data and passes it to `cellIterface` constructor to create an object that offers the access to data in a convenient form. When adding an extension, one can specify which particle types are to be affected in the parameter `subject`. The *Handler* is called by $\pi$-PIC during `advance()` for each subset of particles in each cell (non-empty for the specified types) for each specified type. The cellIterface provides access to the data of all particles in the subset, and to the electromagnetic field at the cell corners. In addition, the interface provides functionality for adding new particles within the same cell (new particles are kept unseen until the next loop). If the subject includes `cells`, the Handler is called for each cell. This can be used, e.g., to make an extension to account for ionisation. 

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

In case of contributing to the toolbox of $\pi$-PIC, the name of an extension should be short (<32 characters) but sufficiently descriptive. One practical choice is to use one word/acronym describing the physics in question and then after an underscore a reference (first author and year) to the article that described the method in use, e.g. `qed_volokitin2023`.

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
Creating an extension in C/C++ provides full flexibility and performance, while also permitting the creation of small extensions compiled locally. The layout of development is similar to the one for Python developments. In case of contributing to the toolbox of $\pi$-PIC, the following routine is recommended (if the extension concerns a modification of $\pi$-PIC core/interface, consider using the procedure described in the [next section](#core-development)): 


1. Create a local folder for your extension, name it as the extension:  
	`mkdir <name>`  
	`cd <name>`  
2. Clone the project from there:  
	`git clone https://github.com/hi-chi/pipic.git`  
	`cd pipic`  
3. Create a branch:   
	`git branch <name>`   
	`git checkout <name>`  
4. (optional) reinstall pipic:  
	`pip uninstall pipic`  
    `python3 -m pip install .`  
    (optional) check that it works:  
	    `python3 examples/basic_example.py`  
5. Go to the src folder and create a folder named by the name of your extension:  
	`cd src/extensions`  
	`mkdir <name>`  
	`cd <name>`  
6. Clone pybind:  
	`git clone https://github.com/pybind/pybind11`  
7. Bring files to local development:  
	`cp ../../primitives.h .`  
	`cp ../../interfaces.h .`  
	`cp ../../CMakeLists.txt .`  
8. Amend CMakeLists.txt:  
	`set(<name>`  
	`<name>.cpp)`  
	`pybind11_add_module(<name>${<name>})`  
	(optional) remove four lines that concern FFTW, and `#include <fftw3.h>` from `primitives.h` 
9. Copy files for testing, e.g.:  
	`cp ../x_converter_c/x_converter_c.cpp .`  
	`cp ../../../examples/x_converter_c_test.py .`  
11. Rename the cpp file to match the name of the extension: `<name>.cpp`  
12. Amend as follows a few lines in `<name>.cpp` for testing:  
	`const string name = "<name>";`  
	`cout << “hi from <name>” << endl;` add this line to the body of `handler()`  
	`PYBIND11_MODULE(<name>, object) {`  
13. Rename the py file to `<name>_test.py` and amend it as follows:  
	`import <name>`   
	`extension_handler = <name>.handler(location=-L/4-L/32, thickness=L/16,`  
	`typeTo=sim.get_type_index('electron'))`  
	`sim.add_handler(name=<name>.name, subject='electron1',`   
	`handler=extension_handler)`  
14. Generate so-file:  
    `cmake .`  
    `make`  
15. Test by running:  
	`python3 <name>.py`  
16. Develop the extension; to recompile use:  
	`make`  
17. For testing and debugging purposes it can be useful to disable the use threads by setting `use_omp = False` among parameters of `advance()`. When commiting changes to git do not add copied/generated files.  

### When finished with the developments and tests, update and upload the branch:  

17. Add a licence statement to `<name>.cpp` and to all needed .h files.  
    In `<name>.cpp` add and an underscore in front of the name in the following line:  
	`PYBIND11_MODULE(_<name>, object) {`  
18. Add the example of use named `<name>_test.py` to the folder examples; amend the import:
	`from pipic.extensions import <name>`  
19. Register the extension in `pipic/extensions/__init__.py`  
20. Test everything by reinstalling from scratch:  
	`pip uninstall pipic`  
	`python3 -m pip install .`  
	`python3 examples/<name>_test.py`  
22. Add the following files to the branch (do not add copied/generated files):  
    `git add src/extensions/<name>/<name>.cpp` and other needed .h files  
	`git add pipic/extensions/__init__.py`  
    `git add examples/<name>_test.py`  
23. Add a note about the user interfaces, assumptions and everything needed to the [documentation](../EXTENSIONS.md).  
24. Upload the branch, make a pull request  


As an example see an extension [`x_converter_c.cpp`](https://github.com/hi-chi/pipic/blob/main/src/extensions/x_converter_c/x_converter_c.cpp) that can be used to convert type of particles in an $x$-limited region. The use of this extension is exemplified in [`x_reflector_py_test.py`](https://github.com/hi-chi/pipic/blob/main/examples/x_converter_c_test.py) (note that this extension is added to $\pi$-PIC as an example).

# Core development

- Create a local folder for your extension:  
    `mkdir <pipic_name>`  
    `cd <pipic_name>`
- Clone the project from there:  
    `git clone https://github.com/hi-chi/pipic.git`  
    `cd pipic`  
- Create a branch:  
    `git branch <name>`   
    `git checkout <name>`  
- Uninstall pipic (if it has been installed previously, see pip list):  
    `pip uninstall pipic`  
- Install from the current location:  
    `python3 -m pip install .`  
- (optional) Check that it works:  
    `python3 examples/basic_example.py`  
- Go to the src folder:  
    `cd src`  
- Clone pybind:  
    `git clone https://github.com/pybind/pybind11`  
- Change source files, e.g. in pipic.h add a line to pipic::pipic():  
    `cout << “running test version” << endl;`  
- In pipic.h rename the struct pipic to test_pipic, three occasions: name (line 28), constructor (line 37) and destructor (line 51). [this is needed to avoid name conflict]  
- In pipic.cpp add a line in the beginning (e.g. after line 41):  
    `typedef test_pipic pipic;`  
- Generate so-file:  
    `cmake .`  
    `make`  
- Copy a test py file into local folder  
- Along with pipic import `_pipic` and call `_pipic.init` instead of `pipic.init`  
- Test it by running a test file  
- To recompile:   
    `make`  
