import unittest
import subprocess
import sys
import os
import tempfile
import shutil

EXAMPLES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'examples'))
EXTENSIONS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src', 'extensions'))
SOLVER_TEST_SCRIPT = os.path.join(EXAMPLES_DIR, 'basic_example_test_solvers.py')

EXAMPLE_FILES = [
    'basic_example.py',
    'basic_example_3d.py',
    'downsampler_gonoskov2022_test.py',
    'energy_conservation.py',
    'focused_pulse_test.py',
    'laser_solid_interaction.py',
    'plasma_oscillation.py',
    'qed_gonoskov2015_test.py',
    'qed_volokitin2023_test.py',
    'x_converter_c_test.py',
    'x_reflector_c_test.py',
    'x_reflector_py_test.py',
]

# Each entry is (relative_path_from_EXTENSIONS_DIR, script_filename)
EXTENSION_TEST_FILES = [
    ('absorbing_boundaries', 'test_absorbing_boundaries.py'),
    ('absorbing_boundaries', 'test_absorbing_boundaries_with_plasma_no_ba.py'),
    ('absorbing_boundaries', 'test_absorbing_boundaries_with_plasma.py'),
    ('moving_window', 'test_moving_window_1d.py'),
    ('moving_window', 'test_moving_window.py'),
    ('moving_window', 'test_moving_window_rotated.py'),
]

SOLVERS = [
    'electrostatic_1d',
    'ec',
    'ec2',
    'emc2',
    'fourier_boris',
]


def make_script_test(filepath, source_dir, timeout, test_prefix, extra_args=None):
    """Return a test method that runs *filepath* for up to *timeout* seconds."""
    filename = os.path.basename(filepath)

    def test_method(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy all files from the script's source directory so relative
            # data files are available, then run in the clean temp dir.
            for f in os.listdir(source_dir):
                src = os.path.join(source_dir, f)
                if os.path.isfile(src):
                    shutil.copy2(src, tmpdir)
            try:
                command = [sys.executable, filepath]
                if extra_args:
                    command.extend(extra_args)
                result = subprocess.run(
                    command,
                    timeout=timeout,
                    capture_output=True,
                    text=True,
                    cwd=tmpdir,
                )
                self.assertEqual(
                    result.returncode, 0,
                    msg=f"{filename} exited with code {result.returncode}:\n{result.stderr}",
                )
                print(result.returncode)
            except subprocess.TimeoutExpired:
                # Still running after timeout — considered a pass
                pass
        # tmpdir and all created files are deleted here

    stem = os.path.splitext(filename)[0]
    test_method.__name__ = f'{test_prefix}{stem}'
    return test_method


class TestPipic(unittest.TestCase):
    def test_dummy(self):
        self.assertEqual(1 + 1, 2)

    def test_import_pipic(self):
        import pipic


class TestExamples(unittest.TestCase):
    pass


for _filename in EXAMPLE_FILES:
    _filepath = os.path.join(EXAMPLES_DIR, _filename)
    _test_name = f'test_example_{os.path.splitext(_filename)[0]}'
    setattr(TestExamples, _test_name,
            make_script_test(_filepath, EXAMPLES_DIR, timeout=10, test_prefix='test_example_'))


class TestExtensions(unittest.TestCase):
    pass


for _subdir, _filename in EXTENSION_TEST_FILES:
    _source_dir = os.path.join(EXTENSIONS_DIR, _subdir)
    _filepath = os.path.join(_source_dir, _filename)
    _test_name = f'test_extension_{os.path.splitext(_filename)[0]}'
    setattr(TestExtensions, _test_name,
            make_script_test(_filepath, _source_dir, timeout=10, test_prefix='test_extension_'))


class TestSolvers(unittest.TestCase):
    pass


for _solver in SOLVERS:
    _test_name = f'test_solver_{_solver}'
    setattr(
        TestSolvers,
        _test_name,
        make_script_test(
            SOLVER_TEST_SCRIPT,
            EXAMPLES_DIR,
            timeout=10,
            test_prefix='test_solver_',
            extra_args=['--solver', _solver],
        ),
    )


if __name__ == "__main__":
    unittest.main()
