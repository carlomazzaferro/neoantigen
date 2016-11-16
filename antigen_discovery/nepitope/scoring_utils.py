import os
import shlex
import subprocess
import pandas
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tabulate import tabulate



class Score(object):

    """
    Class that implements score calculations along an MSA for every n-mer specified. The algorithm
    is based on the Jensen-Shannon divergence, and can be found in the file socre_conservation.py.
    This class is aimed at created the necessary inputs for the calculation over every window of an
    alignment. Must be initialized with the fasta file containing the alignment as well as a list of
    n-mers to be used for the scoring.
    """

    def __init__(self, fasta_input, nmers):
        """

        :param fasta_input:fasta input containing MSA
        :param nmers: list of n-mers over which scorng will be calculatared
        """
        self.input = fasta_input
        self.nmers = nmers


    @staticmethod
    def create_separate_lists(fasta_file):
        """
        Creates 2 lists from a fasta file
        :param fasta_file: file
        :return: one list for the IDs in the file and one list for the proteins/peptides in it
        """
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
            j = []
            k = []
            for i in all_list:
                if i.startswith('>'):
                    j.append(i)
                else:
                    k.append(i)
        return j, k

    @staticmethod
    def assign_score_mhc_results(mhc_df, conserv_df):

        alleles = list(mhc_df['Allele'].unique())
        nmers = list(mhc_df['n-mer'].unique())

        list_dfs = []

        for allele in alleles:
            for nmer in nmers:

                sliced = Score._get_df_allele_nmer(mhc_df, allele, nmer)
                score_df = Score._get_nmer_score_df(conserv_df, nmer)
                list_scores = Score._get_score_list(score_df)
                final_df = Score._add_list_as_column(sliced, list_scores)
                list_dfs.append(final_df)

        resulting_df = pandas.concat(list_dfs)
        #resulting_df = resulting_df.loc[resulting_df['Score'] > 0]
        resulting_df['Peptide'] = resulting_df['Peptide'].str.replace('X', '-')

        return resulting_df

    @staticmethod
    def _get_score_list(df):
        """

        :param df:
        :return:
        """
        return list(df['Score'].values)

    @staticmethod
    def _add_list_as_column(df_1, list_1):
        if len(list_1) < len(df_1):
            list_2 = [0] * (len(df_1) - len(list_1))
            list_1 = list_1 + list_2
        df_1['Score'] = list_1
        return df_1

    @staticmethod
    def _get_df_allele_nmer(df, allele, nmer):
        return df.loc[(df['Allele'] == allele) & (df['n-mer'] == nmer)]

    @staticmethod
    def _get_nmer_score_df(df, nmer):
        return df.loc[df['n-mer'] == nmer]

    @staticmethod
    def add_conserv_score_to_df_list(list_mhc_dfs, conserv_df):

        agg_list = []
        for i in list_mhc_dfs:
            agg_list.append(Score.assign_score_mhc_results(i, conserv_df))

        return agg_list

    @staticmethod
    def get_length(file_name):
        """
        Get number of lines in a file
        :param file_name: file name
        :return: number of lines
        """
        num_lines = sum(1 for line in open(file_name))
        return num_lines

    @staticmethod
    def return_df(fasta_file):
        """
        Returns dataframe from fasta file without duplicates
        :param fasta_file: file
        :return: df
        """
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

        with open(out_nmers_path + "stderr.txt") as error_file:
            lines = error_file.readlines()
            if lines:
                print ("Errors occured. Check %s for details" % (out_nmers_path + "stderr.txt"))

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


def get_summary_data(list_dfs):
    """
    Pretty print summary data about each protein
    :param list_dfs: list of dataframes containing info about each protein
    :return: prints to output a nicely formatted tables
    """
    for i in list_dfs:
        prot_name = i.ID.unique()[0]
        align_title = i["Alignment Title"].unique()[0]
        align_hits = i["Hits"].unique()[0]
        prot_len = i.Length.unique()[0]
        num_high_affinity = i.loc[i["Affinity Level"] == 'High']
        num_med_affinity = i.loc[i["Affinity Level"] == 'Intermediate']
        num_low_affinity = i.loc[i["Affinity Level"] == 'Low']
        num_no_affinity = i.loc[i["Affinity Level"] == 'No']
        HA_per_AA = float(i['Length'].unique()[0]) / (len(num_high_affinity))

        lists_to_print = [["Accession Number", prot_name],
                          ["Alignment Title", align_title],
                          ["Protein length", prot_len],
                          ["Alignment Hits", align_hits],
                          ["High affinity peptides", len(num_high_affinity)],
                          ["Medium affinity peptides", len(num_med_affinity)],
                          ["Low affinity peptides", len(num_low_affinity)],
                          ["No affinity peptides", len(num_no_affinity)],
                          ["High affinity per amino acid", HA_per_AA]]

        print (tabulate(lists_to_print))


def add_blast_extra_data(list_dfs, extra_data_file):
    list_ = open(extra_data_file, 'r').readlines()
    list_ = [list_[i:i + 4] for i in range(0, len(list_), 4)]

    df_extra_data = create_df_from_list_(list_)

    list_2 = []
    for i in list_dfs:
        list_2.append(pandas.merge(i, df_extra_data, on='ID'))
    return list_2


def f(x):
    return int(x[0]) / float(x[1])


def create_df_from_list_(list_):
    df = pandas.DataFrame(list_, columns=['ID', 'Alignment Title', 'Length', 'Hits'])
    for column in df:
        df[column] = df[column].str.strip()
    df['ID'] = df['ID'].str[1:]
    df['Identity Percentage'] = df[['Hits', 'Length']].apply(f, axis=1)

    return df