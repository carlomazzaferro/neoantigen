import os
import shlex
import subprocess
import pandas
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging


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
        return "All files written to %s" % out_nmers_path

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
                list_dfs.append(self.get_dfs(files))

        scores = self.calculate_avg(list_dfs)

        return scores

    def create_large_fasta(self, out_nmers_path):

        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(message)s',
                            datefmt='%a, %d %b %Y %H:%M:%S')

        for i in range(0, len(self.nmers)):
            list_files = glob.glob(out_nmers_path + "/nmerized_%i*" % self.nmers[i])

            with open(out_nmers_path + 'consolidated_fasta_%i.fasta' % self.nmers[i], 'w') as w_file:
                for filen in list_files:
                    with open(filen, 'rU') as o_file:
                        seq_records = SeqIO.parse(o_file, 'fasta')
                        SeqIO.write(seq_records, w_file, 'fasta')

        consolidated = glob.glob(out_nmers_path + "/consolidated*")

        for i in range(0, len(consolidated)):
            df = Score.return_df(consolidated[i])
            filtered_name = out_nmers_path + 'filtered_' + consolidated[i].split('/')[-1]
            Score.dataframe_to_fasta(df, outfile=filtered_name)
            len_file = Score.get_length(filtered_name)
            print 'File filtered_consolidated_fasta_%i.fasta written to %s' % (self.nmers[i], out_nmers_path)

            if len_file > 1000:
                logging.info('File length > 1000. Will have to split file in smaller chunks. Use split_fasta before '
                                'running netMHCcons prediction.')

    @staticmethod
    def get_length(file_name):
        num_lines = sum(1 for line in open(file_name))
        return num_lines

    @staticmethod
    def return_df(fasta_file):
        recs = SeqIO.parse(fasta_file, 'fasta')
        keys = ['locus_tag', 'translation', 'description']
        data = [(r.name, str(r.seq), str(r.description)) for r in recs]
        df = pandas.DataFrame(data, columns=(keys))
        df = Score.clean_df(df)
        return df.drop_duplicates(subset='translation')

    @staticmethod
    def clean_df(df):
        df = df[df['translation'].str.contains('-') == False]
        return df

    @staticmethod
    def dataframe_to_fasta(df, seqkey='translation', idkey='locus_tag',
                           descrkey='description',
                           outfile='out.faa'):

        seqs = []
        for i, row in df.iterrows():
            rec = SeqRecord(Seq(row[seqkey]), id=row[idkey],
                            description=row[descrkey])
            seqs.append(rec)
        SeqIO.write(seqs, outfile, "fasta")
        return outfile

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


def split_fasta_file(out_nmers_path, nmer):
    """
    Split fasta file so that it contains at each file contains at most ~300 peptides

    :param out_nmers_path: path of filtered_consolidated fasta files
    :return: split files
    """

    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        datefmt='%a, %d %b %Y %H:%M:%S')

    filtered = glob.glob(out_nmers_path + "/filtered_consolidated_fasta_%i*" % nmer)
    for i in range(0, len(filtered)):

        len_file = Score.get_length(filtered[i])
        num_file_splits = (len_file / 300) + 1
        logging.info('Files will be split into %i files' % num_file_splits)
        list_fasta = []

        with open(filtered[i], 'r') as infile:
            for line in infile:
                list_fasta.append(line.rstrip())

        for j in range(0, num_file_splits):
            new_file_names = out_nmers_path + "/split_%i_" % j + filtered[i].split('/')[-1]

            if (j + 1) * 300 < len_file:
                list_slice = list_fasta[j * 300:(j + 1) * 300]
            else:
                list_slice = list_fasta[j * 300::]
            with open(new_file_names, 'w') as outfile:
                for item in list_slice:
                    outfile.write("%s\n" % item)

            print 'File %s written to %s' % (new_file_names.split('/')[-1], out_nmers_path)


class FileConsolidation(object):

    stable_cols = ['Pos', 'Peptide', 'ID']

    def __init__(self, filepath, file_pattern):
        """

        :param fasta_input:
        :param nmers:
        """
        self.filepath = filepath
        self.file_pattern = file_pattern
        self.files = glob.glob(self.filepath + self.file_pattern)
        self.allele_list = self.get_allele_list(self.files)

    def load_batch(self):

        os.chdir(self.filepath)
        list_dfs = []
        list_summaries = []

        for files in self.files:

            df = pandas.read_csv(files, sep='\t', skiprows=1)
            processed, summary = self.concat_sliced(df)
            list_dfs.append(processed)
            list_summaries.append(summary)

        return list_dfs, list_summaries

    def concat_sliced(self, df):

        sliced_cols = self.slice_over_df(df)
        list_dfs = []

        for i in sliced_cols:
            list_dfs.append(self.rename_cols(df[i]))

        major, summary = self.return_concat(self.add_allele_name(list_dfs, self.allele_list))

        return major, summary

    @staticmethod
    def concat_batches(list_dfs, list_summaries):
        conc1 = pandas.concat(list_dfs)
        conc2 = pandas.concat(list_summaries)
        return conc1, conc2

    @staticmethod
    def label_affinity(row):
        if row['nM'] < 50.0:
            return 'High'
        if 50.0 < row['nM'] < 500.0:
            return 'Intermediate'
        if 500.0 < row['nM'] < 5000.0:
            return 'Low'
        if row['nM'] > 5000.0:
            return 'No'

    def aggregate_inf0(self, df1, df2):
        merged = pandas.merge(df1, df2, on=FileConsolidation.stable_cols)
        merged['Affinity Level'] = merged.apply(lambda row: self.label_affinity(row), axis=1)
        return merged

    @staticmethod
    def get_allele_list(files):
        unique_alleles = []
        for filename in files:
            df1 = pandas.read_csv(filename, sep='\t')
            cols = list(df1.columns)

            for item in cols:
                if item.startswith('HL'):
                    unique_alleles.append(item)

        return unique_alleles

    @staticmethod
    def slice_over_df(df):

        all_cols = list(df.columns)
        list_cols = []

        for i in xrange(3, len(all_cols), 3):

            sliced = all_cols[i:(i + 3)]
            sliced.extend(FileConsolidation.stable_cols)
            list_cols.append(sliced)

        return list_cols

    @staticmethod
    def rename_cols(df):

        col_names = list(df.columns)

        for i in col_names:
            if '.' in i:
                df = df.rename(columns={i: i.split('.')[0]})

        return df

    @staticmethod
    def add_allele_name(list_dfs, allele_list):

        for i in range(0, len(list_dfs[:-1])):
            list_dfs[i]['Allele'] = allele_list[i]

        return list_dfs

    @staticmethod
    def return_concat(list_dfs):

        major_df = pandas.concat(list_dfs[:-1])
        summary = list_dfs[-1]

        return major_df, summary

