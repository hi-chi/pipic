import pipic.consts as consts
import pipic.ctypes as ctypes
import pipic.extensions as extensions
import pipic.interfaces as interfaces
import pipic.types as types

__version__: str

electron_charge: float
electron_mass: float
light_velocity: float
proton_mass: float

class init:
    def __init__(
        self,
        solver: str,
        nx: int,
        xmin: float,
        xmax: float,
        ny: int = ...,
        ymin: float = ...,
        ymax: float = ...,
        nz: int = ...,
        zmin: float = ...,
        zmax: float = ...,
    ) -> None: ...
    def add_handler(
        self,
        name: str,
        subject: str,
        handler: int,
        data_double: int = ...,
        data_int: int = ...,
    ) -> None: ...
    def add_particles(
        self,
        name: str,
        number: int,
        charge: float,
        mass: float,
        temperature: float,
        density: int,
        data_double: int = ...,
        data_int: int = ...,
    ) -> None: ...
    def advance(self, time_step: float, number_of_iterations: int = ...) -> None: ...
    def custom_field_loop(
        self,
        number_of_iterations: int,
        it2r: int,
        field2data: int,
        data_double: int = ...,
        data_int: int = ...,
    ) -> None: ...
    def field_loop(
        self,
        handler: int,
        data_double: int = ...,
        data_int: int = ...,
        use_omp: bool = ...,
    ) -> None: ...
    def fourier_solver_settings(
        self, divergence_cleaning: int = ..., sin2_kfilter: int = ...
    ) -> None: ...
    def get_number_of_particles(self) -> int: ...
    def get_type_index(self, type_name: str) -> int: ...
    def log_policy(self, log_to_file: bool = ..., log_to_screen: bool = ...) -> None: ...
    def particle_loop(
        self,
        name: str,
        handler: int,
        data_double: int = ...,
        data_int: int = ...,
    ) -> None: ...
    def set_rng_seed(self, seed: int) -> None: ...
    @property
    def nx(self) -> int: ...
    @property
    def ny(self) -> int: ...
    @property
    def nz(self) -> int: ...
    @property
    def xmax(self) -> float: ...
    @property
    def xmin(self) -> float: ...
    @property
    def ymax(self) -> float: ...
    @property
    def ymin(self) -> float: ...
    @property
    def zmax(self) -> float: ...
    @property
    def zmin(self) -> float: ...
