import pandas as pd
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class SimilartiyScore(object):
    def __init__(self, dataframe_list, reference_protein_id=None):

        self.df_list = dataframe_list
        self.filtered_df = self._first_pass_filter()
        self.pos_list = self._retrieve_positions_list()
        self.ranges_list = self._retrieve_ranges_list()
        self.prot_list = self._get_protein_list()
        self.referece = self._get_ref_prot(reference_protein_id=reference_protein_id)
        self.matches = self._find_matches()
        self.num_matches = [len(i) for i in self.matches]
        self.total_matches = self._find_total_matches()
        self.num_total_matches = [len(i) for i in self.total_matches]
        self.columns = ['Proteins', 'Matches Positions', 'Total Num Matches', 'Ranges']
        self.rank_by_matches = self._get_df_ranked_by_matches()
        # self.ordered_list =
        self.iterative_ranking = self.create_ranking_by_affiniy_location()

    ##############################################################################################
    # FUNCTIONS FOR ITERATIVE RANKING: These functions are aiming at creating the attribute      #
    # iterative_ranking and ordered_list. It works as follows:                                   #
    #  1. Select a protein that will be used 1st as a reference protein.                         #
    #  2. Calculate the number of 'Matches' (interpreted as number of High Affinity              #
    #     AA that share the same position) between reference protein and every other protein     #
    #  3. Append the one with lowest number of matches to a list (priority queue)                #
    #  4. The protein appended becomes the reference protein                                     #
    #  5. Repeat process with the appended protein, which is compared to all the other           #
    #     proteins (excluding the previous one).  Process goes on until every protein is ranked. #
    ##############################################################################################

    def create_ranking_by_affiniy_location(self):
        """
        :param reference_protein_id: to check the proper names of the protein IDs,
        check the attribute 'prot_list' by running SimilairtyScore.prot_list for
        a complete list.
        :return: list of of peptides ordered as defined above.
        """
        list_ordered = []
        base_df = self.rank_by_matches
        first_entry = self._get_first_entry(base_df, self.referece)
        list_ordered.append(first_entry)
        base_df = base_df[base_df.Proteins != self.referece]

        for i in range(len(base_df)):
            if i == 0:
                df = base_df
            else:
                df = remaining_data

            new_ref, remaining_ranges, remaining_prots = self._get_remaining_data(
                df.sort_values(by='Total Num Matches'))
            new_matches_list = self._find_total_matches_by_ref(new_ref[-1], remaining_ranges)
            num_new_matches = self._get_num_new_matches(new_matches_list)

            remaining_data = self._get_iterative_df_with_matches(remaining_prots, new_matches_list,
                                                                 num_new_matches, remaining_ranges)
            list_ordered.append(new_ref)

            df = self._get_ordered_list_into_pandas(list_ordered)
            df['Percentage Match To Self'] = pd.Series(self._find_percentage_match(df))
            df['Percentage Match To Prev'] = pd.Series(self._find_percentage_match_to_prev(df))

        return df.drop_duplicates(subset='Proteins')

    def _get_ordered_list_into_pandas(self, ordered_list):
        df = pd.DataFrame(ordered_list, columns=self.columns)
        return df

    def _get_remaining_data(self, df):
        return self._get_new_ref(df), self._get_remaining_ranges(df), self._get_remaining_proteins(df)

    @staticmethod
    def _get_first_entry(base_df, reference_protein_id):
        first_item = base_df.loc[base_df['Proteins'] == reference_protein_id].values.tolist()[0]
        first_item[2] = 'All'
        return first_item

    @staticmethod
    def _get_num_new_matches(matches_list):
        return [len(i) for i in matches_list]

    @staticmethod
    def _get_new_ref(df):
        return df.head(1).values.tolist()[0]

    @staticmethod
    def _get_remaining_ranges(df):
        return df['Ranges'][1::].values.tolist()

    @staticmethod
    def _get_remaining_proteins(df):
        return df['Proteins'][1::].values.tolist()

    def _get_iterative_df_with_matches(self, remaining_prots, matches_positions, num_matches, ranges):

        data_ = [remaining_prots, matches_positions, num_matches, ranges]
        df = pd.DataFrame(data_, index=self.columns)

        return df.T

    @staticmethod
    def _find_total_matches_by_ref(ref, remaining_ranges):

        ref_ranges = [item for sublist in ref for item in sublist]
        total_matches = []

        for range_list in remaining_ranges:
            comparison_ranges = [item for sublist in range_list for item in sublist]
            intersec_list = list(set(ref_ranges).intersection(set(comparison_ranges)))
            total_matches.append(intersec_list)

        return total_matches

    ################################################################################
    # FUNCTIONS FOR PLOTTING: Data to plotted are the peptide ranges per protein   #
    # methods available for affinity levels of selected protein against all others #
    # or for pairs of proteins.                                                    #
    ################################################################################

    def plot_affinity_regions(self):
        plot_data = self._retrieve_plot_data(self.ranges_list)
        protein_list = self.prot_list
        self.plotting_function(plot_data, protein_list)

    def plot_iterative_ranking_plot(self):
        to_plot = self._retrieve_plot_data(self.iterative_ranking.Ranges.values.tolist())
        protein_list = self.iterative_ranking.Proteins.values.tolist()
        color_data = self.iterative_ranking['Percentage Match To Prev'][1::].values.tolist()
        color_data.insert(0, 0.0)
        self.plotting_function(to_plot, protein_list, color_data)

    @staticmethod
    def _retrieve_plot_data(range_list):

        concatd_ranges_list = []
        for range_ in range_list:
            concatd_ranges_list.append(list(set([item for sublist in range_ for item in sublist])))

        to_plot = [[ranges, [i + 1] * len(ranges)] for i, ranges in enumerate(concatd_ranges_list)]

        return to_plot

    @staticmethod
    def _get_color_map(color_data, plot_data):

        color_list = []
        for i, color in enumerate(color_data):
            interpolated_col = np.interp(color, [0, max(color_data)], [-1, 1])
            full_list = np.full(len(plot_data[i][0]), interpolated_col)
            color_list.append(full_list)

        return color_list

    def plotting_function(self, plot_data, prot_list, color_data):
        plt.figure(figsize=(15, 25))
        ax = plt.gca()
        ax.set_ylim([-1, len(prot_list) + 1])
        ax.invert_yaxis()
        plt.title('High Affinity Score Locations')
        colors = self._get_color_map(color_data, plot_data)

        for i, coordinate in enumerate(plot_data):
            plt.scatter(coordinate[0], coordinate[1], c=colors[i], linewidth='0', s=50, vmin=-1, vmax=1,
                        label=prot_list[i], cmap='rainbow_r')
            # plt.scatter(coordinate[0], coordinate[1], color=colors, label=prot_list[i])

        plt.legend(bbox_to_anchor=(1.2, 1), loc=1, ncol=1)
        plt.show()

    ##############################################################################
    # FUNCTIONS TO RETRIEVE CLASS ATTRIBUTES: Mostly helper functions to get the #
    # shared data within the class.                                              #
    ##############################################################################

    def _find_percentage_match(self, df):

        high_aa_count_list = self._high_aa_count(df)
        list_matches = df['Total Num Matches'][1::].values.tolist()
        percentage_list = [a / float(b) for a, b in zip(list_matches, high_aa_count_list)]
        percentage_list.insert(0, 100)

        return percentage_list

    def _find_percentage_match_to_prev(self, df):

        list_matches = df['Total Num Matches'][1::].values.tolist()
        ref = df['Len High AA'][:-1].values.tolist()
        percentage_to_ref_list = [a / float(b) for a, b in zip(list_matches, ref)]
        percentage_to_ref_list.insert(0, 100)

        return percentage_to_ref_list

    def _high_aa_count(self, df):

        df['Len High AA'] = df['Ranges'].apply(self._get_high_aa_count)
        return df['Len High AA'][1::].values.tolist()

    @staticmethod
    def _get_high_aa_count(x):
        ref_ranges = [item for sublist in x for item in sublist]  # flatten
        return len(ref_ranges)

    def _get_ref_prot(self, reference_protein_id=None):
        if reference_protein_id:
            return reference_protein_id
        else:
            reference = [self.pos_list[0][j][0] for j in range(len(self.pos_list[0]))]
            return reference

    def _find_matches(self):
        matches = []
        ref = [self.pos_list[0][j][0] for j in range(len(self.pos_list[0]))]

        for i in range(1, len(self.pos_list)):
            yy = [self.pos_list[i][j][0] for j in range(len(self.pos_list[i]))]
            intersec_list = list(set(ref).intersection(set(yy)))


            matches.append(intersec_list)
            matches.append(list(set(ref).intersection(set(yy))))

        return matches

    def _find_total_matches(self):

        ref_ranges = [item for sublist in self.ranges_list[0] for item in sublist]
        total_matches = []

        for range_list in self.ranges_list:
            comparison_ranges = [item for sublist in range_list for item in sublist]
            intersec_list = list(set(ref_ranges).intersection(set(comparison_ranges)))
            total_matches.append(intersec_list)

        return total_matches

    def _get_df_ranked_by_matches(self):

        data_ = [self.prot_list, self.total_matches, self.num_total_matches, self.ranges_list]
        df = pd.DataFrame(data_, index=['Proteins', 'Matches Positions', 'Total Num Matches', 'Ranges'])

        return df.T

    def _retrieve_ranges_list(self):

        ranges_list = []
        for single_list in self.pos_list:
            ranges_list.append([list(range(pos[0], pos[1] + 1)) for pos in single_list])

        return ranges_list

    def _first_pass_filter(self):
        filtered_dfs = []
        for i in self.df_list:
            filt_ = i.loc[i['Affinity Level'] == 'High']
            clean = self._return_clean(filt_)
            filtered_dfs.append(clean)

        return filtered_dfs

    def _get_protein_list(self):
        prots = []
        for i in self.filtered_df:
            prots.append(i.ID.unique()[0])

        if len(prots) == len(self.df_list):
            return prots

        else:
            prots = []
            for i in self.df_list:
                prots.append(i.ID.unique()[0])

            return prots

    @staticmethod
    def _return_clean(df):
        df['Peptide'] = df['Peptide'].str.replace('X', '-')
        df = df.loc[df['Peptide'].str.contains('--') == False]
        return df

    def _retrieve_positions_list(self):
        _pos_ID = []
        for i in self.filtered_df:
            i['Final Pos'] = i['Pos'] + i['n-mer']
            i = i.sort_values(by='Pos')
            _pos_ID.append(i[['Pos', 'Final Pos']].values.tolist())

        return _pos_ID

    # def _plot_in_euclidean_space(some_lists)

    @staticmethod
    def _calc_min_dist(arr1, arr2):
        return min(product(arr1, arr2), key=lambda t: abs(t[0] - t[1]))

    @staticmethod
    def _pad_smaller_list(list_):
        threshold = 5
        for i in range(len(ls_2)):
            distance = ls_2[i] - ls_1[i]
            print(distance)
            if distance > threshold:
                ls_2.insert(i, 0)