import pandas
from pymongo import MongoClient
import warnings


class DataEnrichment(object):
    #TODO write tests for pretty much everything
    """
    Class designed to integrate gene expression data and to easily parse/filter/visualize it.
    It is initialized with collection name, database name, pvac_seq file name and gene_expression file name.
    It assumes the existance of a MongoDB instance with stored data
    """
    def __init__(self, db_name, collection_name, pvac_seq_file, gene_file):
        self.collection_name = collection_name
        self.db_name = db_name
        self.pvac_seq_file = pvac_seq_file
        self.gene_file = gene_file

    def retrieve_data_from_db(self, chr=None, start=None, end=None):
        """
        Retrieves an entry at a time from the MongoDB.
        :param chr: chromosome
        :param start: mutation Ref
        :param end: mutation Alt
        :return: list with a single document retrieved from MongoDB
        """
        if chr is not None and start is not None and end is not None:

            client = MongoClient()
            db = getattr(client, self.db_name)
            collection = getattr(db, self.collection_name)
            return list(collection.find({"$and": [{"Chr": chr}, {"Start": start+1}, {"End": end}]}))

        else:
            raise ValueError("Location values not provided")

    def retrieve_multiple(self):
        """
        Bt inheriting the file names from the class, this method retrieves chromosome, ref and alt data from
        dataframes generated from the input files and queries from MongoDB the data associated with those
        variants.
        :return: list of documents with variant data
        """

        list_queries = self.retrieve_chr_start_end()
        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)

        list_out = []
        for i in list_queries:
            query = collection.find({"$and": [{"Chr": i[0]}, {"Start": i[1]+1}, {"End": i[2]}]})
            list_out.append(list(query))

        return list_out

    def retrieve_read_counts(self):
        """
        From a list of documents retreived from MongoDB, retrive the read counts per allele

        :return: read_counts for each variant allele
        """
        list_queries = self.retrieve_multiple()
        read_counts = []
        for i in list_queries:
            info = list([i[0]["Genotype_Call"]["AD"][0],
                         i[0]["Genotype_Call"]["AD"][1],
                         i[0]["Genotype_Call"]["DP"],
                         i[0]["hgvs_key"]])

            read_counts.append(info)

        return read_counts

    def integrate_data(self, type_data=None):
        """
        Integrates data coming from mongoDB with the existent dataframe. type_data
        refers to either of the following options: Normal, Tumor DNA ,Tumor RNA.
        :param type_data: Normal, Tumor DNA ,Tumor RNA
        :return: dataframe with integrated data
        """
        df1, df2 = self.get_dfs()

        # Initialize empty columns to be filled
        df1['Gene Exp FPKM'] = None
        df1['AD_REF'] = None
        df1['AD_ALT'] = None
        df1['DP'] = None
        df1['hgvs_key'] = None

        to_query = self.retrieve_chr_start_end()
        read_and_ids = self.retrieve_read_counts()

        if len(to_query) != len(read_and_ids):
            raise ValueError("Mismatch between number of queries and queried data")

        new_cols = ['AD_REF', 'AD_ALT', 'DP', 'hgvs_key']

        for i in range(0, len(to_query)):
            df1.loc[(df1.Chromosome == to_query[i][0]) & (df1.Start == to_query[i][1]) &
                   (df1.Stop == to_query[i][2]), new_cols] = read_and_ids[i]
            df1 = self.add_gene_data(df1, df2)

        return self.rename_cols(df1, type_data)

    def get_dfs(self):
        """
        interal method to retrieve dfs from files
        :return: a dataframe for each file
        """
        pvac_seq_df = pandas.read_csv(self.pvac_seq_file, sep='\t')
        gene_data_df = pandas.read_csv(self.gene_file, sep='\t')
        gene_data_df = gene_data_df.rename(columns={'FPKM': 'Gen Exp'})
        return pvac_seq_df, gene_data_df

    def retrieve_chr_start_end(self):
        pvac_df = self.get_dfs()[0]
        to_query = pandas.unique(pvac_df[['Chromosome', 'Start', 'Stop']].values)
        return list(to_query)

    @staticmethod
    def rename_cols(df, type_data):
        """
        Internal Method to rename columns and changing counts to int
        :param df: pvaq-seq-data-associated dataframe
        :param type_data: "Normal", "Tumor DNA" or "Tumor RNA"'
        :return: df with renamed columns
        """
        df['AD_ALT'] = df['AD_ALT'].apply(lambda x: int(x))
        df['AD_REF'] = df['AD_REF'].apply(lambda x: int(x))

        if type_data == 'Normal':
            df = df.rename(columns={'AD_REF': 'Normal Ref Count', 'AD_ALT': 'Normal Var Count'})
        elif type_data == 'Tumor DNA':
            df = df.rename(columns={'AD_REF': 'Tumor DNA Ref Count', 'AD_ALT': 'Tumor DNA Var Count'})
        elif type_data == 'Tumor RNA':
            df = df.rename(columns={'AD_REF': 'Tumor RNA Ref Count', 'AD_ALT': 'Tumor RNA Var Count'})
        else:
            raise ValueError('Invalid data type provided. Choose between "Normal", "Tumor DNA" and "Tumor RNA"')

        return df

    @staticmethod
    def add_gene_data(df, gene_data):
        """
        :param df:
        :param gene_data:
        :return: Normal, Tumor DNA ,Tumor RNA
        """
        possible_genes = df['Ensembl Gene ID'].unique()
        for i in possible_genes:
            try:
                df.loc[df['Ensembl Gene ID'] == i, 'Gene Exp FPKM'] = gene_data.loc[gene_data['gene_id'] == i].values[0][-1]
            except:
                warnings.warn("Data for gene {} is not available, replacing it with -1".format(i))
                df.loc[df['Ensembl Gene ID'] == i, 'Gene Exp FPKM'] = -1
        return df


class GeneInfo(object):
    """
    Class that enables the user to retrieve some important data about the genes in the dataset
    """

    def __init__(self, pvac_seq_file, gene_file):
        self.pvac_seq_file = pvac_seq_file
        self.gene_file = gene_file

    def visualize_df(self, pvac_seq=None, gene_data=None):
        """
        Self-explanatory
        :param pvac_seq: True of False, will display if true
        :param gene_data: True of False, will display if true
        :return: dataframe view
        """

        if (pvac_seq is None) & (gene_data is None):
                raise ValueError("No file specification provided")
        else:
            if pvac_seq is True:
                return self.get_dfs()[0]

            if gene_data is True:
                return self.get_dfs()[1]

    def check_valid_gene_ids(self):
        """
        This method will output a summary of the gene data
        :return: None
        """
        df1, df2 = self.get_dfs()
        valid = list(df1['Ensembl Gene ID'].isin(df2['gene_id']))
        valid_num = sum(valid).astype(float)/len(valid)

        possible_genes = df1['Ensembl Gene ID'].unique()
        print "N possible genes: %i \n" % len(possible_genes)
        print "Percentage of found genes: {}%\n".format(valid_num)

        for i in possible_genes:
            dat = (df2.loc[df2['gene_id'] == i].values)
            if len(dat) == 0:
                print 'No gene data found for gene {}'.format(i)
            else:
                print "Gene {}: FPKM: {}".format(dat[0][0], dat[0][-1])

    def get_dfs(self):
        """
        interal method to retrieve dfs from files
        :return: a dataframe for each file
        """
        pvac_seq_df = pandas.read_csv(self.pvac_seq_file, sep='\t')
        gene_data_df = pandas.read_csv(self.gene_file, sep='\t')
        gene_data_df = gene_data_df.rename(columns={'FPKM': 'Gen Exp'})
        return pvac_seq_df, gene_data_df


def write_tsv_out(filepath, df):
    df.to_csv(filepath, sep='\t')




