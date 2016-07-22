import unittest
import sys
sys.path.append("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes")
from variantAnnotation import utilities
from variantAnnotation import csv_to_df


file_name = "/Users/carlomazzaferro/Desktop/CSV to be tested/Tumor_targeted_processed.csv"
sample_list = csv_to_df.open_and_parse("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/test data/UnitTestData/test_MiniCsv.csv")
sample_df = csv_to_df.parse_to_df(sample_list)

larger_list = csv_to_df.open_and_parse(file_name)
larger_df = csv_to_df.parse_to_df(larger_list)


class TestUtilities(unittest.TestCase):

    def setUp(self):
        self.data = file_name
        self.cytoBand = ['1p36.33', '16p11.1', '16q11.2', '16q21', 'Xp22.32', 'Xp22.2', 'Xp22.11',
                         'Xq12', 'Xq13.1']

        self.cytoBand_split = [
                                ['1', 'p', '36', '.', '33'],
                                ['16', 'p', '11', '.', '1'],
                                ['16', 'q', '11', '.', '2'],
                                ['16', 'q', '21'],
                                ['X', 'p', '22', '.', '32'],
                                ['X', 'p', '22', '.', '2'],
                                ['X', 'p', '22', '.', '11'],
                                ['X', 'q', '12'],
                                ['X', 'q', '13', '.', '1']
                              ]

        self.cytoBand_dict = [
                                {'Arm': 'p', 'Band': 6, 'Chromosome': '1', 'Region': 3, 'Sub_Band': 33},
                                {'Arm': 'p', 'Band': 1, 'Chromosome': '16', 'Region': 1, 'Sub_Band': 1},
                                {'Arm': 'q', 'Band': 1, 'Chromosome': '16', 'Region': 1, 'Sub_Band': 2},
                                {'Arm': 'q', 'Band': 1, 'Chromosome': '16', 'Region': 2},
                                {'Arm': 'p', 'Band': 2, 'Chromosome': 'X', 'Region': 2, 'Sub_Band': 32},
                                {'Arm': 'p', 'Band': 2, 'Chromosome': 'X', 'Region': 2, 'Sub_Band': 2},
                                {'Arm': 'p', 'Band': 2, 'Chromosome': 'X', 'Region': 2, 'Sub_Band': 11},
                                {'Arm': 'q', 'Band': 2, 'Chromosome': 'X', 'Region': 1},
                                {'Arm': 'q', 'Band': 3, 'Chromosome': 'X', 'Region': 1, 'Sub_Band': 1}
                             ]

    def test_split_cytoBand(self):

        cytolist = []
        for i in range(0, len(self.cytoBand)):
            cytolist.append(utilities.split_cytoband(self.cytoBand[i]))

        self.assertEqual(cytolist, self.cytoBand_split)

    def test_cytoBand_to_dict(self):

        cytodict = []
        for i in self.cytoBand:
            print cytodict
            cytodict.append(utilities.lists_to_dict(utilities.split_cytoband(i)))

        self.assertEqual(cytodict, self.cytoBand_dict)


    def

#HGVS_id creation not tested since it comes straight from myvariant.info's implementation: tested before.






    #def test_split_string(self):
    #    list1 = larger_df["ExonicFunc.knownGene"].dropna()


       # self.assertEqual()ls




if __name__ == '__main__':
    unittest.main()
