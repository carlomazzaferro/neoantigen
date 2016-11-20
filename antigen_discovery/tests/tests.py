import unittest
import os
import sys
import pandas
import glob


from nepitope import scoring_utils
from nepitope.scoring_utils import Score



file_name = os.path.dirname(os.path.realpath('__file__')) + '/tests/test_scoring.fasta'
scoring = scoring_utils.Score(file_name, [9, 10, 11])

l = [0.86281,  0.65371, 0.84304,  0.78496, 0.82407,  0.77538,  0.83054,  0.68654, 0.77331]
mean_l = reduce(lambda x, y: x + y, l) / float(len(l))


class TestScoringUtils(unittest.TestCase):

    def setUp(self):
        self.data = file_name

        self.test_list = [[0.7815955555555556, 'MDKKYSIGL', 9],
                          [0.7633500000000001, 'DKKYSIGLD', 9],
                          [0.8262666666666667, 'IGTNSVGWA', 9],
                          [0.7322544444444445, 'LEESFLVEE', 9],
                          [0.7556699999999998, 'EESFLVEED', 9],
                          [0.7949977777777778, 'ESFLVEEDK', 9],
                          [0.7851955555555556, 'SFLVEEDKK', 9],
                          [0.7558333333333334, 'FLVEEDKKH', 9],
                          [0.7114655555555555, 'LVEEDKKHE', 9],
                          [0.7133622222222222, 'VEEDKKHER', 9]]

        self.df_cols = ['Align_Col_Number', 'Score', 'Column', 'Peptide']

        self.score_pep_val = [[mean_l, 'MDKKYSIGL', len('MDKKYSIGL')]]

    def test_scoring_utils(self):
        self.assertTrue(self, os.path.exists(os.getcwd()+'/score_conservation.py'))

    def test_conservation_scipt(self):
        curr_dir = os.path.dirname(os.path.realpath('__file__')) + '/tests/'
        scoring.create_score_file(file_name, curr_dir)
        self.assertTrue(os.stat(curr_dir + "stderr.txt").st_size == 0)

    def test_get_dfs(self):
        curr_dir = os.path.dirname(os.path.realpath('__file__')) + '/tests/'
        file_name_ = curr_dir + "processed_testing.txt"

        headers = ['Align_Col_Number', 'Score', 'Column']
        df_ = pandas.read_csv(file_name_, sep='\t', comment="#", names=headers)
        df_ = Score.add_pep_col(df_)

        lst_df = [df_]
        summary = Score.calculate_avg(lst_df)

        self.assertEqual(list(df_.columns), self.df_cols)
        self.assertEqual(summary, self.score_pep_val)

    def test_run_scoring(self):
        curr_dir = os.path.dirname(os.path.realpath('__file__')) + '/tests/'
        nmers = [8,9,10,11]

        list_dfs = []
        for i in nmers:
            for files in glob.glob(curr_dir + "/processed_nmerized_%i*" % i):
                list_dfs.append(scoring_utils.Score.get_dfs(files))
        self.assertEqual(len(list_dfs), len(glob.glob(curr_dir + "/processed_nmerized_*")))

    def test_create_list_from_fasta(self):
        file_name = self.data
        num_lines = sum(1 for line in open(file_name))
        lines_list = scoring_utils.Score.create_lists(file_name)
        lines = [line.rstrip('\n') for line in open(file_name)]

        self.assertEqual(lines[0], lines_list[0])
        self.assertEqual(num_lines, len(lines_list))

    def test_create_lists_from_fasta(self):

        file_name = self.data
        num_lines = sum(1 for line in open(file_name))
        line_list_1, line_list_2 = scoring_utils.Score.create_separate_lists(file_name)
        self.assertEqual(num_lines, len(line_list_1) + len(line_list_2))







if __name__ == '__main__':
    unittest.main()