import os
import shlex
import subprocess
import pandas
import glob


class Score(object):

    score_script_path = "python " + os.path.dirname(os.path.realpath('__file__')) + "/score_conservation.py "

    def __init__(self, fasta_input, nmers):
        """

        :param fasta_input:
        :param nmers:
        """
        self.input = fasta_input
        self.nmers = nmers

    def make_windows(self,  out_nmers_path):
        lines = self.create_lists(self.input)
        length = len(lines[1])
        for j in self.nmers:
            for i in range(0, (length / j) + 1):
                with open(out_nmers_path + 'nmerized_%i_%i' % (j, i), 'w') as outfile:
                    for line in lines:
                        if '>' in line:
                            outfile.write(line + '\n')
                        else:
                            outfile.write(line[i:i + j] + '\n')

    @staticmethod
    def create_lists(fasta_file):

        with open(fasta_file) as infile:
            all_list = []
            peptide = ""
            lines = infile.readlines()
            for i in range(0, len(lines)):
                if lines[i].startswith('>'):
                    all_list.append(lines[i].rstrip())
                else:
                    peptide += lines[i].rstrip()
                try:
                    if lines[i + 1].startswith('>'):
                        all_list.append(peptide)
                        peptide = ""
                except:
                    all_list.append(peptide)
        return all_list

    def run_scoring(self, out_nmers_path):

        for i in self.nmers:
            for files in glob.glob(out_nmers_path + "/nmerized_%i*" %i):
                self.create_score_file(files, out_nmers_path)

        list_dfs = []
        for i in range(0, len(self.nmers)):
            for files in glob.glob(out_nmers_path + "processed_nmerized_%i*" % self.nmers[i]):
                print files
                list_dfs.append(self.get_dfs(files))

        scores = self.calculate_avg(list_dfs)

        return scores


    @staticmethod
    def create_score_file(files, out_nmers_path):

        os.chdir(os.path.dirname(os.path.realpath('__file__')))

        args_str = Score.score_script_path + files
        args = shlex.split(args_str)

        with open(out_nmers_path + "processed_" + files.split('/')[-1] + ".txt", "wb") as out, open(out_nmers_path + "stderr.txt", "wb") as err:
            subprocess.Popen(args, stdout=out, stderr=err)

    @staticmethod
    def get_dfs(files):

        headers = ['Align_Col_Number', 'Score', 'Column']
        df_ = pandas.read_csv(files, sep='\t', comment="#", names=headers)
        return Score.add_pep_col(df_)

    @staticmethod
    def most_common(lst):
        return max(set(lst), key=lst.count)

    @staticmethod
    def add_pep_col(df):

        aa_list = list(df["Column"])
        peptide = []

        for i in aa_list:
            peptide.append(Score.most_common(i))
        pep = ''.join(peptide)
        df['Peptide'] = pep

        return df


    @staticmethod
    def calculate_avg(list_dfs):
        list_values = []
        for i in list_dfs:
            score = i["Score"].mean()
            peptide = list(i["Peptide"])[0]
            nmer = len(peptide)
            list_values.append([score, peptide, nmer])

        return list_values






