import os
import shlex
import subprocess
import pandas
import glob
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging
from tabulate import tabulate
from skbio import TabularMSA, Protein


class Score(object):

    """
    Class that implements score calculations along an MSA for every n-mer specified. The algorithm
    is based on the Jensen-Shannon divergence, and can be found in the file socre_conservation.py.
    This class is aimed at created the necessary inputs for the calculation over every window of an
    alignment. Must be initialized with the fasta file containing the alignment as well as a list of
    n-mers to be used for the scoring.
    """

    score_script_path = "python " + os.path.dirname(os.path.realpath('__file__')) + "/score_conservation.py "

    def __init__(self, fasta_input, nmers):
        """

        :param fasta_input:fasta input containing MSA
        :param nmers: list of n-mers over which scorng will be calculatared
        """
        self.input = fasta_input
        self.nmers = nmers

    def make_windows(self,  out_nmers_path):
        """
        deprecated method
        :param out_nmers_path: file to the windowized n-mer files
        :return: creates files over which score conservation will be calculated
        """
        lines = self.create_lists(self.input)
        length = len(lines[1])
        print (length)
        for j in self.nmers:
            for i in range(0, length - j +1):
                with open(out_nmers_path + 'nmerized_%i_%i' % (j, i), 'w') as outfile:
                    for line in lines:
                        if '>' in line:
                            outfile.write(line + '\n')
                        else:
                            outfile.write(line[i:i + j] + '\n')
        return "All files written to %s" % out_nmers_path

    @staticmethod
    def create_lists(fasta_file):
        """
        Create list from a fasta file
        :param fasta_file: file
        :return: list with every row as entry in a list
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
        return all_list

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
        """
        Joins files containing windowized n-mers
        :param out_nmers_path: filepath to windows
        :return: large fasta file containing all the files in the path that have pepties of the same length
        """

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

            print ('File filtered_consolidated_fasta_%i.fasta written to %s' % (self.nmers[i], out_nmers_path))

            if len_file > 1000:
                logging.info('File length > 1000. Will have to split file in smaller chunks. Use split_fasta before '
                                'running netMHCcons prediction.')

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

            print('File %s written to %s' % (new_file_names.split('/')[-1], out_nmers_path))

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


class Alignment(object):

    def __init__(self, msa_file, ref_protein_file):
        self.msa_file = msa_file
        self.msa = self.read_msa_file()
        self.positional_conservation = self.get_positional_conservation()
        self.reference_protein = ref_protein_file
        self.reference_protein_string = self.open_fasta_return_single_seq()

    def read_msa_file(self):

        msa = TabularMSA.read(self.msa_file, constructor=Protein)
        msa.reassign_index(minter='id')
        return msa

    def get_positional_conservation(self):

        positional_conservation = self.msa.conservation(metric='inverse_shannon_uncertainty',
                                                        degenerate_mode='nan', gap_mode='include')
        return np.nan_to_num(positional_conservation)

    def create_score_df_from_scikit_bio(self, nmers):

        ls_df = []

        for i in nmers:
            scores = []
            peptides = []

            for j in range(0, len(self.reference_protein_string) - i):
                scores.append(np.mean(self.positional_conservation[j:j + i]))
                peptides.append(self.reference_protein_string[j:j + i])

            df = pandas.DataFrame([scores, peptides], index=['Score', 'Peptide'])
            df = df.T
            df['n-mer'] = i
            ls_df.append(df)

        return pandas.concat(ls_df)

    def open_fasta_return_single_seq(self):

        with open(self.reference_protein) as inf:

            prot = []
            next(inf)
            for line in inf:
                prot.append(line)

            return prot[0]

class AddData (object):

    def __init__(self, msa_file_input, msa_file_output, scores_df, positional_conservation,
                 all_alleles=True, list_alleles=None, pos_cons_treshold=None):

        if pos_cons_treshold is None:
            self.pos_cons_treshold = 0.1
        else:
            self.pos_cons_treshold = pos_cons_treshold
        self.msa_file_input = msa_file_input
        self.msa_file_output = msa_file_output
        self.scores_df = scores_df
        self.all_alleles = all_alleles
        self.list_alleles = list_alleles
        self.positional_conservation = positional_conservation
        self.alleles = self._check_return_alleles(self.scores_df, self.all_alleles, self.list_alleles)
        self.nmers = self._get_nmers_from_affinity_df(self.scores_df)
        self.high_aa_low_cons_df = self._high_aff_low_cons_to_df(self.return_high_affinity_and_not_conserved())

    def open_files(self):

        with open(self.msa_file_input) as inf, open(self.msa_file_output, 'w') as out:
            self.write_conservation_scores(inf, out)
            self.write_affinity_scores(out)

    def write_conservation_scores(self, inf, out):

        for line in inf:
            line = line.replace('X', '-')
            out.write(line)
        out.write('>CONSERVATION_INFO\n')

        for i in self.positional_conservation:
            if i > self.pos_cons_treshold:
                out.write('O')
            else:
                out.write('-')

    def write_affinity_scores(self, out):

        for nmer in self.nmers:
            for allele in self.alleles:

                to_print = self._slice_df(nmer, allele, self.scores_df)
                peps = self._get_peptides(to_print)

                for idx in range(0, len(peps)):

                    if idx > 3250:
                        continue

                    if '--' in peps[idx]:
                        continue

                    if not self._get_affinity_per_peptide(peps[idx], to_print):
                        continue
                    else:
                        self._write_out(nmer, allele, idx, out, peps)

    def high_affinity_low_cons_df(self):
        selected_df = self.scores_df.loc[(self.scores_df['Affinity Level'] == 'High') & (self.scores_df['Score'] < 0.1)]
        selected_df = selected_df.loc[(selected_df['Pos'] < 3250) & (selected_df['Peptide'].str.contains('--') == False)]
        return selected_df

    def return_high_affinity_and_not_conserved(self):

        high_aff_not_cons = []

        for nmer in self.nmers:
            for allele in self.alleles:

                to_print = self._slice_df(nmer, allele, self.scores_df)
                peps = self._get_peptides(to_print)

                for idx in range(0, len(peps)):

                    mean_cons = self._get_mean_pos_cons_per_pep(nmer, idx)
                    if self._get_affinity_per_peptide(peps[idx], to_print):
                        if mean_cons < self.pos_cons_treshold:
                            print (mean_cons)
                            high_aff_not_cons.append([idx, peps[idx]])

        return high_aff_not_cons

    @staticmethod
    def _high_aff_low_cons_to_df(list_of_lists):
        return pandas.DataFrame(list_of_lists, columns=['Peptide Position', 'Peptide'])

    def _get_mean_pos_cons_per_pep(self, nmer, index):

        initial_aminoa_acid = index*nmer
        endind_amino_acid = (index+1)*nmer

        return np.mean(self.positional_conservation[initial_aminoa_acid:endind_amino_acid])

    @staticmethod
    def _write_out(nmer, allele, idx, out, peps):

        out.write('\n>High_Affinity_Loc|n-mer=%i|allele=%s\n' % (nmer, allele))
        out.write('-' * idx)
        out.write(peps[idx])
        out.write('-' * (len(peps) - idx - 1))

    @staticmethod
    def _get_affinity_per_peptide(pep, df):
        aff_per_pep = df.loc[df['Peptide'] == pep]
        if len(aff_per_pep) > 1:
            return False

        if list(aff_per_pep['Affinity Level'].values)[0] == 'High':
            return True
        else:
            return False

    @staticmethod
    def _slice_df(nmer, allele, df):

        to_print = df.loc[(df['n-mer'] == nmer) & (df['Allele'] == allele)]
        to_print['Peptide'] = to_print['Peptide'].str.replace('X', '-')
        return to_print

    @staticmethod
    def _get_peptides(df):
        return list(df['Peptide'].values)

    @staticmethod
    def _check_return_alleles(scores_df, all_alleles, list_alleles):

        if all_alleles:
            alls = list(scores_df.Allele.unique())
        else:
            alls = list_alleles

        if (all_alleles is False) & (list_alleles is None):
            raise ValueError('No allele provided')

        return alls

    @staticmethod
    def _get_nmers_from_affinity_df(scores_df):
        return list(scores_df['n-mer'].unique())
