import unittest
import sys
sys.path.append('/Users/carlomazzaferro/Documents/CCBB/neoantigen/VarP-master/')

import os
import glob
from VarP import utils
from shutil import copyfile
import pandas as pd

class TestHandleReaders(unittest.TestCase):

    def setUp(self):
        self.test_directory = os.path.dirname(os.path.realpath('__file__'))
        self.reader_names = ['reader_1', 'reader_2', 'reader_3', 'reader_4']
        self.handler = utils.HandleReaders(self.reader_names)
        self.test_df = pd.DataFrame()

    def test_create_collection_from_reader(self):
        num_vcf_files = 4
        self.create_empty_vcf_files(num_vcf_files-1)

        collection = self.handler.create_collection_from_readers(self.test_directory)
        self.assertEqual(len(collection),num_vcf_files)

    def test_coding_effect_creator(self):
        num_vcf_files = 4
        self.create_empty_vcf_files(num_vcf_files-1)

        var_coll = self.handler.create_collection_from_readers(self.test_directory)[0] #select first
        self.assertEqual(str(type(var_coll)), "<class 'varcode.variant_collection.VariantCollection'>")

    def test_return_df(self):
        num_vcf_files = 4
        self.create_empty_vcf_files(num_vcf_files-1)

        var_coll = self.handler.create_collection_from_readers(self.test_directory) #select firstc
        coding_effects = self.handler.return_list_coding_effects(var_coll)

        protein_list = self.handler.return_protein_list(coding_effects)

        self.assertEqual(len(protein_list), len(coding_effects))

        df = self.handler.return_dataframe(protein_list, coding_effects)

        self.assertEqual(type(self.test_df), type(df))

    @staticmethod
    def create_empty_vcf_files(num_files):
        for i in range(0, num_files):
            copyfile('test_mini.vcf', str(i) + 'test_mini.vcf')
    def tearDown(self):

        #Clean Up
        filelist = glob.glob("*.vcf")
        for f in filelist:
            if f.startswith("t"):
                pass
            else:
                os.remove(f)

if __name__ == '__main__':
    unittest.main()
