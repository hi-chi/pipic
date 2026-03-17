import unittest
import subprocess
import sys
import os
import tempfile
import shutil

TEST_SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test_scripts'))
EXAMPLES_DIR = os.path.join(TEST_SCRIPTS_DIR, 'examples')
EXTENSIONS_DIR = os.path.join(TEST_SCRIPTS_DIR, 'test_extensions')
SOLVER_TEST_SCRIPT = os.path.join(TEST_SCRIPTS_DIR, 'test_solvers', 'basic_example_test_solvers.py')

EXAMPLE_FILES = [
    'basic_example.py',
    'basic_example_3d.py',
    'energy_conservation.py',
    'laser_solid_interaction.py',
    'plasma_oscillation.py',
]

def collect_extension_test_files(root_dir):
    """Return extension test scripts from tests/test_scripts/test_extensions."""
    items = []
    for entry in sorted(os.scandir(root_dir), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        for filename in sorted(os.listdir(entry.path)):
            if filename.endswith('_test.py') or filename.startswith('test_'):
                items.append((entry.name, filename))
    return items


EXTENSION_TEST_FILES = collect_extension_test_files(EXTENSIONS_DIR)

SOLVERS = [
    'electrostatic_1d',
    'ec',
    'ec2',
    'emc2',
    'fourier_boris',
]

PIC_ADVANCEMENTS_PER_TEST = 1


def make_script_test(filepath, source_dir, test_prefix, extra_args=None):
    """Return a test method that runs *filepath* and requires clean completion."""
    filename = os.path.basename(filepath)

    def test_method(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy all files from the script's source directory so relative
            # data files are available, then run in the clean temp dir.
            for f in os.listdir(source_dir):
                src = os.path.join(source_dir, f)
                if os.path.isfile(src):
                    shutil.copy2(src, tmpdir)
            command = [sys.executable, filepath]
            command.extend(["--steps", str(PIC_ADVANCEMENTS_PER_TEST)])
            if extra_args:
                command.extend(extra_args)
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                cwd=tmpdir,
            )
            self.assertEqual(
                result.returncode, 0,
                msg=f"{filename} exited with code {result.returncode}:\n{result.stderr}",
            )
            print(result.returncode)
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
            make_script_test(_filepath, EXAMPLES_DIR, test_prefix='test_example_'))


class TestExtensions(unittest.TestCase):
    pass


for _subdir, _filename in EXTENSION_TEST_FILES:
    _source_dir = os.path.join(EXTENSIONS_DIR, _subdir)
    _filepath = os.path.join(_source_dir, _filename)
    _test_name = f'test_extension_{os.path.splitext(_filename)[0]}'
    setattr(TestExtensions, _test_name,
            make_script_test(_filepath, _source_dir, test_prefix='test_extension_'))


class TestSolvers(unittest.TestCase):
    pass


for _solver in SOLVERS:
    _test_name = f'test_solver_{_solver}'
    setattr(
        TestSolvers,
        _test_name,
        make_script_test(
            SOLVER_TEST_SCRIPT,
            os.path.dirname(SOLVER_TEST_SCRIPT),
            test_prefix='test_solver_',
            extra_args=['--solver', _solver],
        ),
    )


if __name__ == "__main__":
    unittest.main()
