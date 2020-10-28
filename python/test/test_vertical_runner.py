import unittest
from python.PCA import vertical_pca_benchmark as benchmark


class TestBenchmark(unittest.TestCase):
    def test_assure_consecutive(self):
        self.assertEqual(benchmark.assure_consecutive([0, 1, 2, 3, 4, 5, 6]), 6)
        self.assertEqual(benchmark.assure_consecutive([1, 2, 3, 4, 5, 6, 8]), 5)
        self.assertEqual(benchmark.assure_consecutive([1, 4, 3, 4, 5, 6, 8]), 0)
        self.assertEqual(benchmark.assure_consecutive([]), -1)
if __name__ == '__main__':
    unittest.main()