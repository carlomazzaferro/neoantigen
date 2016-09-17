import pandas
from tabulate import tabulate
import matplotlib.pyplot as plt
import seaborn as sns
import numpy

sns.set(style='ticks')


class SummaryData(object):
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
    def summarize_all_data(cls, list_container, show_names=False):
        obj = cls()
        obj.container = list_container
        obj.proteins = obj.get_proteins(obj.container)
        obj.titles = obj.get_title()
        obj.num_prots = len(obj.proteins)
        obj.list_high_affinity_peps = obj.get_list_affinity(level='high')
        obj.list_med_affinity_peps = obj.get_list_affinity(level='med')
        obj.list_low_affinity_peps = obj.get_list_affinity(level='low')
        obj.list_no_affinity_peps = obj.get_list_affinity(level='no')
        obj.lengths = obj.get_list_lengths()
        obj.hits = obj.get_list_hits()
        obj.list_high_affinity_per_aa = obj.get_list_high_affinity_per_aa()

        if show_names:
            obj.display = obj.display_proteins()

        obj.data_list = [obj.proteins, obj.titles, obj.list_high_affinity_per_aa,
                         obj.list_high_affinity_peps, obj.list_med_affinity_peps,
                         obj.list_low_affinity_peps, obj.list_no_affinity_peps,
                         obj.lengths, obj.hits]

        obj.indexes = ['Accession ID', 'Title', 'High Affinity Peptides Per AA',
                       'Num High Affinity Peps', 'Num Med Affinity Peps', 'Num Low Affinity Peps',
                       'Num No Affinity Peps', 'Protein Length', 'Alignment Hits']

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
        obj.my_df = obj.get_
        obj.title = obj.my_df["Alignment Title"].unique()[0]
        obj.high_affinity_peps = obj.get_num_high_affinity(obj.my_df)
        obj.med_affinity_peps = obj.get_num_med_affinity(obj.my_df)
        obj.low_affinity_peps = obj.get_num_low_affinity(obj.my_df)
        obj.no_affinity_peps = obj.get_num_no_affinity(obj.my_df)
        obj.high_affinity_per_aa = obj.get_high_affinity_per_aa(obj.my_df)
        obj.length = obj.get_length(obj.my_df)
        obj.single_hits = obj.get_hits(obj.my_df)

        return obj

    def return_dataframe(self, num_display=None, rank_by='High Affinity Peptides Per AA'):

        df = pandas.DataFrame(self.data_list, index=self.indexes)
        df = df.T
        if num_display:
            return df.sort_values(by=rank_by).head(num_display)

        else:
            return df.sort_values(by=rank_by)

    def get_df_no_prot_input(self):
        return self.container[0]

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
        for i in self.container:
            lengths.append(self.get_length(i))
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

    def display_proteins(self):
        for i in range(0, len(self.proteins)):
            print("Protein Accession Number: %s" % self.proteins[i])
            print("Associated Alignment Title: %s \n" % self.titles[i])

    def plot_affinity_versus_conservation_score(self):  # , nmer=9):
        out_df = self.my_df
        out_df = out_df.reset_index(drop=True)
        out_df['Indx'] = out_df.index
        print(out_df['n-mer'].unique())
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
