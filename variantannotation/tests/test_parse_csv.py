import unittest
import sys
#quick and dirty way of importing functions
sys.path.append("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/variantAnnotation")
from csv_to_df import open_and_parse
from csv_to_df import parse_to_df
from csv_to_df import open_and_parse_chunks

file_name = "/Users/carlomazzaferro/Desktop/CSV to be tested/Tumor_targeted_processed.csv"
sample_list = open_and_parse("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/test data/UnitTestData/test_MiniCsv.csv")
sample_df = parse_to_df(sample_list)


class OpenParseTest(unittest.TestCase):

    def setUp(self):
        self.data = file_name

    def test_csv_read_data_headers(self):
        list2 = open_and_parse(self.data)
        self.assertEqual(
            ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene',
             'ExonicFunc.knownGene', 'AAChange.knownGene', 'tfbsConsSites', 'cytoBand', 'targetScanS',
             'genomicSuperDups', 'gwasCatalog', 'esp6500siv2_all', '1000g2015aug_all', 'PopFreqMax', '1000G2012APR_ALL',
             '1000G2012APR_AFR', '1000G2012APR_AMR', '1000G2012APR_ASN', '1000G2012APR_EUR', 'ESP6500si_ALL',
             'ESP6500si_AA', 'ESP6500si_EA', 'CG46', 'clinvar_20140929', 'cosmic70', 'nci60', 'Otherinfo'],

            list2[0]
            )

    def test_parse_to_df(self):
        self.assertEqual(list(sample_df.columns), list(parse_to_df(sample_list).columns))

    def test_parse_csv_chunks(self):
        list3 = open_and_parse_chunks(self.data, 1000, 0)
        list2 = open_and_parse(self.data)
        self.assertEqual(len(list3), len(list2))


"""
    def test_csv_read_data_points(self):
        self.assertEqual(read_data(self.data)[1][7], '87')

    def test_get_min_score_difference(self):
        self.assertEqual(get_min_score_difference(self.parsed_data), 1)

    def test_get_team(self):
        index_value = get_min_score_difference(self.parsed_data)
        self.assertEqual(get_team(index_value), 'Liverpool')
"""

if __name__ == '__main__':
    unittest.main()