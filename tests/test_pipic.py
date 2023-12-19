import unittest


class TestPipic(unittest.TestCase):

    def test_dummy(self):
        self.assertEqual(1+1, 2)

    def test_import_pipic(self):
        import pipic


if __name__ == '__main__':
    unittest.main()
