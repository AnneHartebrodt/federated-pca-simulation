import unittest
from python.import_export import gwas_import


class TestGWASScaling(unittest.TestCase):
    def test_get_nr_snps(self):
        infile = '/home/anne/Documents/featurecloud/gwas/data/hapmap/thin.rec.raw.T'
        self.assertEqual(gwas_import.get_nr_snps(infile), 9696)

if __name__ == '__main__':
    unittest.main()
