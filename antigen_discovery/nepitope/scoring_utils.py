import os
import shlex
import subprocess
import pandas
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging
from tabulate import tabulate


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
        print length
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

            print 'File filtered_consolidated_fasta_%i.fasta written to %s' % (self.nmers[i], out_nmers_path)

            if len_file > 1000:
                logging.info('File length > 1000. Will have to split file in smaller chunks. Use split_fasta before '
                                'running netMHCcons prediction.')

    @staticmethod
    def assign_score_mhc_results(mhc_df, conserv_df):
        """
        Method to ensure the data is being assigned in the correct locations
        :param mhc_df:
        :param conserv_df:
        :return:
        """
        alleles = mhc_df['Allele'].unique()
        nmers = mhc_df['n-mer'].unique()

        list_dfs = []

        for i in alleles:
            for j in nmers:

                sliced = mhc_df.loc[(mhc_df['Allele'] == i) & (mhc_df['n-mer'] == j)]
                score_df = conserv_df.loc[conserv_df['n-mer'] == j]
                sliced['Score'] = score_df['Score']
                list_dfs.append(sliced)

        return pandas.concat(list_dfs)

    @staticmethod
    def add_conserv_score_to_df_list(list_mhc_dfs, conserv_df):
        """
        Method for joining conservation score and binding affinity.

        1-log50k  -	 nM  -	Rank   -  Pos   -   Peptide  -  ID  -  Allele  -  Affinity  -  Level  -  n-mer  -  Score
        :param list_mhc_dfs: list of dataframes containing data regarding mhc binding affinity for each protein
        :param conserv_df: list of dataframes containing data regarding conservation scores for each protein
        :return: concatenated dataframes in a list
        """
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
                print "Errors occured. Check %s for details" % (out_nmers_path + "stderr.txt")

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
    """
    Class to ease up the analysis of multiple files/proteins. It will take as inputs a file containing the results
    from netMHCcons and will provide methods to output the data in nicely formatted pandas dataframes that contain
    a variety of accessory information regarding the proteins in question.
    3 different constructors provided, each one well suited for a specific input.
    """
    stable_cols = ['Pos', 'Peptide', 'ID']  #Columns that are always present


    @classmethod
    def load_full_file(cls, file_name):
        """
        For the situation in which a file containing all possible predictions (combinations of alleles and nmers)
        for a single protein.
        :param file_name: ex: file output from netMHCcons for 1 protein
        :return:class with attributes: protein list (should be a list with 1 object), allele list
        """
        obj = cls()
        obj.files = [file_name]
        obj.allele_list = obj.get_allele_list(obj.files)
        obj.protein_list = obj.get_prot_list(obj.files)

        return obj

    @classmethod
    def load_full_multiprot_file(cls, file_name):
        """
        For the situation in which a file containing all possible predictions (combinations of alleles and nmers)
        for multiple proteins.
        :param file_name: file output from netMHCcons for more than 1 protein
        :return:  protein list, allele list
        """
        obj = cls()
        obj.files = [file_name]
        obj.allele_list = obj.get_allele_list(obj.files)
        obj.protein_list = obj.get_prot_list(obj.files)

        return obj

    @classmethod
    def load_batches(cls, filepath, file_names):
        """
        When netMHC or any other prediction method is run locally, then it might be the case that you'll have multiple
        files, each containing prediction for different alleles/n-mers. This method takes care of loading them all
        and consolidating them for later processing.
        :param filepath: path to files
        :param file_names: ame of files
        :return:
        """
        obj = cls()
        obj.filepath = filepath
        files = []
        for i in file_names:
            files.append(filepath + i)
        obj.files = files
        alleles = []
        for i in file_names:
            alleles.append(os.path.splitext(i)[0].split('_')[-1])
        obj.allele_list = alleles

        return obj

    @classmethod
    def load_batch(cls, filepath, file_pattern):
        """
        Same as above but it will load all files give a name pattern, for instance:
        load_batch('/data/predictions', 'netMHC_QQZWC_protein_')
        :param filepath:
        :param file_pattern:
        :return:
        """
        obj = cls()

        obj.filepath = filepath
        obj.file_pattern = file_pattern
        os.chdir(obj.filepath)
        obj.files = glob.glob(filepath + file_pattern)
        obj.allele_list = obj.get_allele_list(obj.files)

        return obj

    def return_df_list(self):
        """
        Returning a list of dataframes from class.files attribute
        :return:
        """
        list_dfs = []
        if len(self.files) > 1:
            for files in self.files:

                df = pandas.read_csv(files, sep='\t', skiprows=1)
                if 'split' in files.split('/')[-1]:

                    processed = self.concat_sliced(df)
                    list_dfs.append(processed)

                else:
                    df = df[['Pos', 'Peptide', 'nM', 'Rank', 'ID']]
                    df["Allele"] = os.path.splitext(files)[0].split('_')[-1]
                    list_dfs.append(df)

            return list_dfs
        else:
            df = pandas.read_csv(self.files[0], sep='\t', skiprows=1)
            processed = self.concat_sliced(df)
            return processed

    def concat_sliced(self, df):
        """
        Concatenate dataframes in list after slicing them: a netMHC prediction comes in wide format, having 3 columns
        for each allele. Here we slice them and concatenate them in order to get a pandas dataframe in long format

        :param df: wide dataframe from netMHC
        :return: long dataframe
        """

        sliced_cols = self.slice_over_df(df)
        list_dfs = []

        for i in sliced_cols:
            list_dfs.append(self.rename_cols(df[i]))

        major = self.return_concat(self.add_allele_name(list_dfs, self.allele_list))

        return major

    def list_df_by_prot(self, df):

        """
        Slice long dataframe per protein and place them in a list
        :param df: long dataframe
        :return: list of dataframes
        """
        list_dfs = []
        for i in self.protein_list:
            df_of_prot = df.loc[df['ID'] == i]
            df_of_prot = self.aggregate_info(df_of_prot)
            list_dfs.append(df_of_prot)

        return list_dfs

    @staticmethod
    def replace_X_with_underscore(df):
        """
        Necessary for proper netMHC processing.
        :param df: df
        :return: df
        """
        df['Peptide'] = df['Peptide'].str.replace('X', '-')
        return df

    @staticmethod
    def concat_batches(list_dfs):
        """
        Concat pandas in a list
        :param list_dfs: list of dfs
        :return: concat'd df
        """
        conc1 = pandas.concat(list_dfs)
        return conc1

    @staticmethod
    def label_affinity(row):
        """
        Function ot be applied to a dataframe to introduce a column with a label for the binding affinity
        :param row: row of df
        :return: binding affinity level
        """
        if row['nM'] < 50.0:
            return 'High'
        if 50.0 < row['nM'] < 500.0:
            return 'Intermediate'
        if 500.0 < row['nM'] < 5000.0:
            return 'Low'
        if row['nM'] > 5000.0:
            return 'No'

    def aggregate_info(self, df1):
        """
        Add extra data to dataframe: affinity level label, and n-mer
        :param df1: df
        :return: df with extra data
        """

        df1['Affinity Level'] = df1.apply(lambda row: self.label_affinity(row), axis=1)
        df1['n-mer'] = df1['Peptide'].str.len()

        return df1

    @staticmethod
    def get_allele_list(files):
        """
        Retrieve alleles from netMHC file
        :param files: netMHC predictions
        :return: list of alleles
        """
        unique_alleles = []

        for i in files:

            df1 = pandas.read_csv(i, sep='\t')
            cols = list(df1.columns)

            for item in cols:
                if item.startswith('HL'):
                    unique_alleles.append(item)

        return unique_alleles

    @staticmethod
    def get_prot_list(files):
        """
        Retrieve list of proteins from netMHC file
        :param files: netMHC predictions
        :return: list of proteins
        """
        unique_prots = []

        for i in files:
            df1 = pandas.read_csv(i, sep='\t', skiprows=1)
            prot_IDs = list(df1['ID'].unique())
            unique_prots.append(prot_IDs)

        unique_prots = [item for sublist in unique_prots for item in sublist]
        return unique_prots

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
        """
        Add allele as a column to dataframes in a list of dataframes
        :param list_dfs: list of dfs
        :param allele_list: allele list retrieved from netMHC predictions
        :return: list of dataframes
        """

        for i in range(0, len(list_dfs[:-1])):
            list_dfs[i]['Allele'] = allele_list[i]

        return list_dfs

    @staticmethod
    def return_concat(list_dfs):

        major_df = pandas.concat(list_dfs[:-1])

        return major_df


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

        print tabulate(lists_to_print)


def add_blast_extra_data(list_dfs, extra_data_file):
    """
    Retrieve data from a blast alignment xml file and place it into a dataframe
    :param list_dfs: list of dataframes
    :param extra_data_file: blast alignment file
    :return: list of dataframes with extra data
    """
    list_ = open(extra_data_file, 'r').readlines()
    list_ = [list_[i:i + 4] for i in xrange(0, len(list_), 4)]
    df_extra_data = pandas.DataFrame(list_, columns=['ID', 'Alignment Title', 'Length', 'Hits'])
    df_extra_data['Identity Percentage'] = df_extra_data.apply(lambda row: (float(row['Hits'])/float(row['length'])))

    list_2 = []
    for i in list_dfs:
        list_2.append(pandas.merge(i, df_extra_data, on='ID'))

    return list_2


