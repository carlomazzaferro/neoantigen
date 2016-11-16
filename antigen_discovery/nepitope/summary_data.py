import pandas
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy
from nepitope import pep_utils



class SummaryData(object):

    @classmethod
    def summarize_data_for_each_mhc_pred(cls, list_container, original_peptides_data, results_of_interest):
        """
        Only one currently used
        :param list_container:
        :param original_peptides_data:
        :param results_of_interest:
        :return:
        """
        obj = cls()
        obj.original_data = original_peptides_data
        obj.results_of_interest = results_of_interest
        obj.container = list_container
        obj.nmer_list = obj.get_nmer_list()
        obj.alleles = obj.get_allele_list()
        obj.list_high_affinity_peps = obj.get_list_num_high_affinity()
        obj.med_affinity_peps = obj.get_list_num_med_affinity()
        obj.low_affinity_peps = obj.get_list_num_low_affinity()
        obj.no_affinity_peps = obj.get_list_num_no_affinity()
        obj.original_peps = [pep[0] for pep in original_peptides_data]
        obj.original_pos = [int(pep[1]) for pep in original_peptides_data]
        obj.top_20_percent_peptides_per_prediction = obj.get_top_scoring_peptides()
        obj.summary_df = obj.return_df_from_lol()

        return obj

    @classmethod
    def summarize_without_blast_dat(cls, list_container):
        obj = cls()
        obj.container = list_container
        obj.proteins = obj.get_proteins(obj.container)
        obj.my_df = obj.get_df_no_prot_input()
        obj.alleles = obj.get_alleles()
        obj.high_affinity_peps = obj.get_num_high_affinity(obj.my_df)
        obj.med_affinity_peps = obj.get_num_med_affinity(obj.my_df)
        obj.low_affinity_peps = obj.get_num_low_affinity(obj.my_df)
        obj.no_affinity_peps = obj.get_num_no_affinity(obj.my_df)
        obj.high_cons_score = obj.get_high_cons_score_peps()
        obj.num_high_cons_peps = len(obj.high_cons_score)
        obj.high_affinity_and_conserved = obj.get_high_affinity_and_conserved()
        obj.num_high_affinity_and_conserved = len(obj.high_affinity_and_conserved)
        obj.low_cons_score = obj.get_low_cons_score_peps()
        obj.high_affinity_and_not_conserved = obj.get_low_affinity_and_not_conserved()
        obj.num_high_affinity_and_not_conserved = len(obj.high_affinity_and_not_conserved)

        obj.table = [['Accession ID', obj.proteins[0]], ['Num High Affinity Peps', obj.high_affinity_peps],
                     ['Num Med Affinity Peps', obj.med_affinity_peps], ['Num Low Affinity Peps', obj.low_affinity_peps],
                     ['Num No Affinity Peps', obj.no_affinity_peps], ['Num High Affinity and Conserved',
                                                                      obj.num_high_affinity_and_conserved],
                     ['Num High Affinity Not Conserved', obj.num_high_affinity_and_not_conserved, 'AFSFSFD']]

        return obj

    @classmethod
    def summarize_all_data(cls, list_container, fasta_file):
        obj = cls()
        obj.container = list_container
        obj.fasta_file = fasta_file
        obj.proteins = list(pandas.concat(obj.container).ID.unique())
        obj.prot_lengths = obj.get_list_lengths()
        obj.num_prots = len(obj.proteins)
        obj.list_high_affinity_peps = obj.get_list_affinity(level='high')
        obj.list_med_affinity_peps = obj.get_list_affinity(level='med')
        obj.list_low_affinity_peps = obj.get_list_affinity(level='low')
        obj.list_no_affinity_peps = obj.get_list_affinity(level='no')

        obj.data_list = [obj.proteins, obj.list_high_affinity_peps, obj.list_med_affinity_peps,
                         obj.list_low_affinity_peps, obj.list_no_affinity_peps, obj.prot_lengths]

        obj.indexes = ['Protein', 'Num High Affinity Peps', 'Num Med Affinity Peps', 'Num Low Affinity Peps',
                       'Num No Affinity Peps', 'Protein Length']

        return obj

    @classmethod
    def summarize_protein_data(cls, list_container, protein_of_interest):
        """
        For the situation in which a file containing all possible predictions (combinations of alleles and nmers)
        for a single protein.
        :param file_name: ex: file output from netMHCcons for 1 protein
        :return:class with attributes: protein list (should be a list with 1 object), allele list
        """
        obj = cls()
        obj.container = list_container
        obj.proteins = obj.get_proteins(obj.container)
        obj.my_protein = protein_of_interest
        obj.my_df = obj.get_df_from_prot()
        obj.title = obj.my_df["Alignment Title"].unique()[0]
        obj.high_affinity_peps = obj.get_num_high_affinity(obj.my_df)
        obj.med_affinity_peps = obj.get_num_med_affinity(obj.my_df)
        obj.low_affinity_peps = obj.get_num_low_affinity(obj.my_df)
        obj.no_affinity_peps = obj.get_num_no_affinity(obj.my_df)
        obj.high_affinity_per_aa = obj.get_high_affinity_per_aa(obj.my_df)
        obj.length = obj.get_length(obj.my_df)
        obj.single_hits = obj.get_hits(obj.my_df)

        return obj

#--------#---FUNCTIONS USED WITH 1ST CLASSMETHOD SOLELY---#--------#--------#--------

    def write_peptides_to_fasta(self, fasta_dir):
        peps = self.summary_df['top scoring peptides'].values.tolist()
        other_data = self.summary_df[['original peptide', 'allele', 'nmer', 'original pos']].values.tolist()
        with open(fasta_dir, 'w') as out:
            for i in range(len(peps)):
                for pep in peps[i]:
                    out.write(
                        ">Peptide:" + other_data[i][0] + '_' + other_data[i][1] + '_' + str(
                            other_data[i][-1]) + '\n')
                    out.write(pep + '\n')

    def return_df_from_lol(self):
        headers = ['nmer', 'allele', 'num high affinity peps', 'num med affinity peps', 'num low affinity peps',
                   'num no affinity peps', 'original peptide', 'original pos', 'top scoring peptides']

        deep_list = [self.nmer_list, self.alleles, self.list_high_affinity_peps, self.med_affinity_peps,
                     self.low_affinity_peps, self.no_affinity_peps, self.original_peps, self.original_pos,
                     self.top_20_percent_peptides_per_prediction]

        return pandas.DataFrame(deep_list, index=headers).T

    def get_top_scoring_peptides(self):

        peptides_in_list = []
        for i in self.container:
            high_scoring = i.sort_values(by='nM', ascending=False).head(int(len(i) * 0.2))
            high_scoring['Peptide'] = high_scoring['Peptide'].str.replace('X', '-')
            peptides_in_list.append(high_scoring['Peptide'].values.tolist())

        return peptides_in_list

    def get_allele_list(self):
        allele_list = []
        for i in self.container:
            allele_list.append(i['Allele'].unique()[0])
        return allele_list

    def get_nmer_list(self):
        nmer_list = []
        for i in self.container:
            nmer_list.append(i['n-mer'].unique()[0])
        return nmer_list

    def get_list_num_high_affinity(self):
        high_aff = []
        for i in self.container:
            high_aff.append(len(i.loc[i['Affinity Level'] == 'High']))
        return high_aff

    def get_list_num_med_affinity(self):
        med_aff = []
        for i in self.container:
            med_aff.append(len(i.loc[i['Affinity Level'] == 'Intermediate']))
        return med_aff

    def get_list_num_low_affinity(self):
        low_aff = []
        for i in self.container:
            low_aff.append(len(i.loc[i['Affinity Level'] == 'Low']))
        return low_aff

    def get_list_num_no_affinity(self):
        no_aff = []
        for i in self.container:
            no_aff.append(len(i.loc[i['Affinity Level'] == 'No']))
        return no_aff

# --------#--------#--------#--------#--------#--------#--------#--------#--------#--------#--------#

    def return_dataframe(self, num_display=None, rank_by='High Affinity Peptides Per AA'):

        df = pandas.DataFrame(self.data_list, index=self.indexes)
        df = df.T
        if num_display:
            return df.sort_values(by=rank_by).head(num_display)

        else:
            return df.sort_values(by=rank_by)

    def get_df_no_prot_input(self):
        filtered = self.container[0].loc[(self.container[0]['Score'] > 0) &
                                         (self.container[0]['Peptide'].str.contains('--') == False)]
        return filtered

    def get_df_from_prot(self):
        df = pandas.DataFrame()
        for i in self.container:
            if i["ID"].unique()[0] == self.my_protein:
                df = i

        if df.empty:
            return "Protein not found"

        else:
            return df

    def get_list_affinity(self, level=None):
        list_affinities = []

        if level == 'high':
            for i in self.container:
                list_affinities.append(self.get_num_high_affinity(i))

        if level == 'med':
            for i in self.container:
                list_affinities.append(self.get_num_med_affinity(i))

        if level == 'low':
            for i in self.container:
                list_affinities.append(self.get_num_low_affinity(i))

        if level == 'no':
            for i in self.container:
                list_affinities.append(self.get_num_no_affinity(i))

        return list_affinities

    def get_high_affinity_per_aa(self, df):
        return float(df['Length'].unique()[0]) / self.high_affinity_peps

    def get_list_high_affinity_per_aa(self):
        my_list = []
        for i in range(0, len(self.container)):
            my_list.append(float(self.container[i]['Length'].unique()[0]) / self.list_high_affinity_peps[i])
        return my_list

    def get_list_lengths(self):

        lengths = []
        idx, seq = pep_utils.create_separate_lists(self.fasta_file)
        idx = [id_.strip('>') for id_ in idx]
        joint_idx_seq = list(zip(idx, seq))

        for prot in self.proteins:
            for pair in joint_idx_seq:
                if prot[0:15] == pair[0][0:15]:
                    lengths.append(len(pair[1]))
        return lengths

    def get_list_hits(self):

        hits_list = []
        for i in self.container:
            hits_list.append(self.get_hits(i))
        return hits_list

    def get_high_cons_score_peps(self):
        high_scores = self.my_df.loc[self.my_df['Score'] > 0.3]
        return high_scores

    def get_low_cons_score_peps(self):
        low_scores = self.my_df.loc[self.my_df['Score'] < 0.3]
        return low_scores

    def get_high_affinity_and_conserved(self):
        return self.high_cons_score.loc[self.high_cons_score['Affinity Level'] == 'High']

    def get_low_affinity_and_not_conserved(self):
        return self.low_cons_score.loc[self.low_cons_score['Affinity Level'] == 'High']

    def get_alleles(self):
        return list(self.my_df['Allele'].unique())

    @staticmethod
    def get_hits(df):
        return df.Hits.unique()[0]

    @staticmethod
    def get_length(df):
        return df.Length.unique()[0]

    @staticmethod
    def get_num_high_affinity(df):
        return len(df.loc[df["Affinity Level"] == 'High'])

    @staticmethod
    def get_num_med_affinity(df):
        return len(df.loc[df["Affinity Level"] == 'Intermediate'])

    @staticmethod
    def get_num_low_affinity(df):
        return len(df.loc[df["Affinity Level"] == 'Low'])

    @staticmethod
    def get_num_no_affinity(df):
        return len(df.loc[df["Affinity Level"] == 'No'])

    @staticmethod
    def get_proteins(list_dfs):
        """
        Retrieve list of proteins from df list
        :param list_dfs: df list
        :return: list of proteins
        """
        unique_prots = []

        for i in list_dfs:
            prot_IDs = i['ID'].unique()
            unique_prots.append(prot_IDs)

        unique_prots = [item for sublist in unique_prots for item in sublist]
        return unique_prots

    def print_table(self):
        print(tabulate(self.table))

    def get_title(self):

        titles = []
        for i in self.container:
            titles.append(i["Alignment Title"].unique()[0])

        return titles

    """

    def plot_affinity_versus_conservation_score(self):  # , nmer=9):
        out_df = self.my_df
        out_df = out_df.reset_index(drop=True)
        out_df = out_df.loc[out_df['Peptide'].str.contains('--') == False]
        out_df['Indx'] = out_df.index
        # print(out_df['n-mer'].unique())
        numpy.random.seed(0)

        _Affinity = ['High', 'Low', 'Intermediate']
        col = sns.color_palette("hls", n_colors=3)
        sns.set(font_scale=1.8)
        fg = sns.FacetGrid(data=out_df, palette=col, hue='Affinity Level',
                           hue_order=_Affinity, aspect=1.4, size=15)

        fg.map(plt.scatter, 'Indx', 'Score', s=250, linewidth=0.1, edgecolor="white").set_axis_labels("Window",
                                                                                                      "Consevation")
        fg.map(plt.axhline, y=out_df['Score'].mean(), xmin=0, xmax=1,
               color='b', alpha=0.5, ls='dotted', lw=1.3,
               label='Mean Conservation Score').set_axis_labels("Peptide Window", "Consevation").add_legend()
    """