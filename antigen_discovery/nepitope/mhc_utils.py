import pandas
import os
import glob


class FileConsolidation(object):
    """
    Class to ease up the analysis of multiple files/proteins. It will take as inputs a file containing the results
    from netMHCcons and will provide methods to output the data in nicely formatted pandas dataframes that contain
    a variety of accessory information regarding the proteins in question.
    4 different constructors provided, each one well suited for a specific input.
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
        :return: a dataframe if only one file is provided, or a list of dataframes if multiple files are provided
        """
        list_dfs = []
        if len(self.files) > 1:
            for idx, files in enumerate(self.files):

                df = pandas.read_csv(files, sep='\t', skiprows=1)
                if 'split' in files.split('/')[-1]:

                    processed = self.concat_sliced(df)
                    list_dfs.append(processed)

                else:
                    df = df[['Pos', 'Peptide', 'nM', 'Rank', 'ID']]
                    df["Allele"] = self.allele_list[idx]
                    df = self.aggregate_info(df)
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
                if item.startswith('H'):
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

        for i in range(3, len(all_cols), 3):

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
