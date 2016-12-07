import unittest
import csv
import os
import sys
import glob
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/')

from mhc_parser import utilities, methods, models


class TestModels(unittest.TestCase):

    def setUp(self):

        self.net_mhc_files_test = glob.glob('/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser/' \
                                 'tests/test_files/mhc_preds_test_fasta_reduced/*.xls')
        self.reduced_fasta = '/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser/tests/' \
                             'test_files/test_fasta_reduced.fasta'

        self.model_methods = models.PredictionCollection(self.net_mhc_files_test, self.reduced_fasta)
        self.model_methods.digest_multiple()
        self.df_list = self.model_methods.return_protein_df_list()


    def test_return_df(self):
        print(self.df_list[0])


    if __name__ == '__main__':
        unittest.main()