import unittest
import sys
#quick and dirty way of importing functions
from variantannotation import csv_to_df

sys.path.append('/Users/carlomazzaferro/Documents/Code/variant-annotation/variantannotation')
file_name = "/Users/carlomazzaferro/Desktop/CSV to be tested/Tumor_targeted_processed.csv"
sample_list = csv_to_df.open_and_parse("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/test data/UnitTestData/test_MiniCsv.csv")
sample_df = csv_to_df.parse_to_df(sample_list)


class MongoDBQueryTest(unittest.TestCase):

    def setUp(self):
        self.data = file_name