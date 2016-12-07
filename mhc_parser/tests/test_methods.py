import unittest
import glob
import sys
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/')
from mhc_parser import utilities, methods, models
import pandas


class TestModels(unittest.TestCase):

    def setUp(self):

        self.test_data = [['', '', '', 'HLA-B1501', '', '', '', ''],
                          ['Pos', 'Peptide', 'ID','nM', 'Rank','Core', 'H_Avg_Ranks','N_binders'],
                          [0, 'MDKKYSIGL', 'Test1', 9694.8, 14, 'MDKKYSIGL', 14, 0],
                          [1, 'DKKYSIGLD', 'Test1', 42058.3, 95, 'DKKYSIGLD', 95, 0],
                          [2, 'KKYSIGLDI', 'Test1', 22207.5, 31, 'KKYSIGLDI', 31, 0],
                          [3, 'KYSIGLDIG', 'Test1', 29965.4, 48, 'KYSIGLDIG', 48, 0],
                          [4, 'YSIGLDIGT', 'Test2', 14.2, 12, 'YSIGLDIGT', 2, 1]]

        self.filtered_test_data = [['', '', '', 'HLA-B1501', '', '', '', ''],
                                   ['Pos', 'Peptide', 'ID','nM', 'Rank','Core', 'H_Avg_Ranks','N_binders'],
                                   [4, 'YSIGLDIGT', 'Test2', 14.2, 12, 'YSIGLDIGT', 2, 1]]


        self.net_mhc_file_test = glob.glob('/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser/' \
                                 'tests/test_files/mhc_preds_test_fasta_reduced/*B1501_9.xls')
        self.reduced_fasta = '/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser/tests/' \
                             'test_files/test_fasta_reduced.fasta'

        self.model_methods = models.PredictionCollection(self.net_mhc_file_test, self.reduced_fasta)
        self.col_data = self.test_data[1]
        self.allele = [i for i in self.test_data[0] if i != ''][0]
        self.df = pandas.DataFrame(self.test_data[2:], columns=self.col_data)
        self.df['Allele'] = [self.allele]*len(self.df)
        self.df['Nmer'] = [len(self.test_data[2][1])]*len(self.df)

    def test_filtering(self):


        threshold = 50
        protein_id = 'Test1'

        #parse filtered and unfiltered data
        #methods.filter_low_affinity requires a prediction collection object (dictionary)
        collection = {}
        collection[protein_id] =  {}
        collection[protein_id]['Predictions'] = self.parsing_function(self.test_data)

        #parse filtered data
        parse_filtered = self.parsing_function(self.filtered_test_data)
        #filter using method to be tested
        methods_filtered = methods.filter_low_affinity(collection, protein_id, threshold)

        self.assertEqual(len(parse_filtered), len(methods_filtered))
        self.assertEqual(parse_filtered[0].Allele, methods_filtered[0].Allele)
        self.assertEqual(parse_filtered[0].Nmer, methods_filtered[0].Nmer)
        self.assertEqual(parse_filtered[0].nM, methods_filtered[0].nM)
        self.assertEqual(parse_filtered[0].Peptide, methods_filtered[0].Peptide)
        self.assertEqual(parse_filtered[0].ID, methods_filtered[0].ID)
        self.assertEqual(parse_filtered[0].lists_to_print, methods_filtered[0].lists_to_print)

    def parsing_function(self, data):
        # already tested and working
        iterator = iter(data)
        allele = self.model_methods.get_allele(iterator)
        _ = next(iterator)
        single_preds = []

        for line in iterator:
            pred = self.model_methods.pass_to_pred_class(line, allele)
            single_preds.append(pred)

        return single_preds

    def test_return_df(self):
        protein_id = 'Test1'
        collection = {}
        collection[protein_id] = {}
        collection[protein_id]['Predictions'] = self.parsing_function(self.test_data)

        df_from_metods = methods.to_df(collection, self.model_methods.cols, protein_id)
        self.assertEqual(list(df_from_metods.columns), list(self.df.columns))
        self.assertEqual(df_from_metods.values.tolist(), self.df.values.tolist())

if __name__ == '__main__':
    unittest.main()
