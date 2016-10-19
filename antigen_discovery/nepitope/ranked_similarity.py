import pandas
import numpy as np
import matplotlib.pyplot as plt
from difflib import SequenceMatcher


class SimilartiyScore(object):

    def __init__(self, dataframe_list, fasta_file, reference_protein_id):

        self.df_list = dataframe_list
        self.fasta_file = fasta_file
        self.reference_pep = reference_protein_id
        self.filtered_df = self._first_pass_filter_and_add_ranges()
        self.possible_align_ranges = self._get_possible_alignment_ranges()[0]
        self.prot_list = self._get_possible_alignment_ranges()[1]
        self.df_sliced_by_range = self._slice_df_by_range()

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

    def _get_possible_alignment_ranges(self):

        possible_nmers = []
        max_val = []
        proteins = []

        for i in self.filtered_df:
            max_val.append(i.Pos.max())
            possible_nmers = list(i['n-mer'].unique())
            proteins.extend(list(i.ID.unique()))

        true_max = max(max_val)
        # possible_nmers = list(set(possible_nmers))

        ranges = []

        for i, nmer in enumerate(sorted(possible_nmers)):
            ranges.append(list(self._get_ranges(true_max - nmer + 2, nmer)))

        return ranges[-1], list(set(proteins))

    @staticmethod
    def _get_ranges(max_val, nmer):

        my_range = []
        for i in range(max_val):
            my_range.append([list(range(i, i + nmer))])

        return my_range

    def _slice_df_by_range(self):

        dfs_in_range = []
        dff = pandas.concat(self.filtered_df)

        for idx, specific_range in enumerate(self.possible_align_ranges):

            mask = dff.Range.apply(self.is_subset, args=(specific_range,))
            df_to_append = dff[mask]

            if len(df_to_append != 0):
                df_to_append['Group'] = idx
                dfs_in_range.append(df_to_append)

        return dfs_in_range

    @staticmethod
    def is_subset(col_set, range_set):
        return set(col_set[0]).issubset(set(range_set[0]))

    def run_pipe(self):

        ordered_protein_list = []
        ordered_scores_list = []
        ordered_protein_list.append(self.reference_pep)
        ordered_scores_list.append(100)

        for i in range(0, len(self.prot_list)):
            reference_pep = self._update_ref_pep(ordered_protein_list)
            dfs_of_interest = self._get_dfs_of_interest(reference_pep, ordered_protein_list)
            concatd_filtered = self._concat_filter(dfs_of_interest, ordered_protein_list)
            groups = self._get_group_numbers(concatd_filtered)
            subset_by_group = self._subset_by_group(groups, concatd_filtered, reference_pep)
            prot_and_score = self._return_protein_and_associated_score(subset_by_group)
            sorted_scores_series = self._make_scoring(prot_and_score, ordered_protein_list)
            ordered_protein_list.append(self._select_top_peptide(sorted_scores_series)[0])
            ordered_scores_list.append(self._select_top_peptide(sorted_scores_series)[1])

        return ordered_protein_list, ordered_scores_list

    @staticmethod
    def _select_top_peptide(my_series):
        prot_name = my_series.index[0]
        prot_score = my_series.values[0]
        return prot_name, prot_score

    @staticmethod
    def _make_scoring(prot_and_score, ordered_protein_list):

        to_df = pandas.DataFrame(prot_and_score[0], columns=prot_and_score[0][0])
        to_df = to_df.T.sort_values(by=0, ascending=False).T[1::]

        iterprots = iter(prot_and_score)
        next(iterprots)

        for i in iterprots:
            len_1_df = pandas.DataFrame(i, columns=i[0]).T.sort_values(by=0, ascending=False).T[1::]
            to_df = to_df.append(len_1_df)
        sorted_series = to_df.sum().sort_values()

        if len(ordered_protein_list) > 1:
            sorted_series = sorted_series[~sorted_series.index.isin(ordered_protein_list[:-1])]

        return sorted_series

    def _return_protein_and_associated_score(self, subset_by_group):

        prot_and_score = []
        for i in subset_by_group:
            i = i.drop_duplicates(subset='ID').sort_values(by=['Similarity To Ref', 'n-mer', 'ID'], ascending=False)
            present_proteins = set(list(i.ID.values))
            other_proteins = [x for x in self.prot_list if x not in present_proteins]
            all_prots = list(i.ID.values) + other_proteins
            scores = list(i['Similarity To Ref'].values)
            scores.extend([0] * len(other_proteins))
            prot_and_score.append([all_prots, scores])

        return prot_and_score

    def _subset_by_group(self, groups, concatd_filtered, reference_pep):

        subset_by_group = []
        for group in groups:
            ref_pep = concatd_filtered.loc[(concatd_filtered['Group'] == group) &
                                           (concatd_filtered['ID'] == reference_pep)].Peptide.values

            if len(ref_pep) != 0:
                subset_df = concatd_filtered.loc[concatd_filtered['Group'] == group]
                subset_df['Similarity To Ref'] = subset_df['Peptide'].apply(self.similar, args=(ref_pep[0],))
                subset_by_group.append(subset_df)

        return subset_by_group

    @staticmethod
    def similar(a, b):
        return SequenceMatcher(None, a, b).ratio()

    @staticmethod
    def _get_group_numbers(concatd_filtered):
        groups = list(concatd_filtered.Group.unique())
        return groups

    @staticmethod
    def _concat_filter(list_df, ordered_protein_list):
        conctttt = pandas.concat(list_df).drop_duplicates(subset=['1-log50k', 'nM', 'Rank', 'Pos', 'Peptide',
                                                                  'ID', 'Allele', 'Affinity Level', 'n-mer',
                                                                  'Final Pos'])

        if len(ordered_protein_list) < 1:
            return conctttt
        else:
            return conctttt[~conctttt.ID.isin(ordered_protein_list[:-1])]

    def _get_dfs_of_interest(self, reference_pep, ordered_protein_list):

        dfs_containing_reference_pep = []

        for df in self.df_sliced_by_range:
            if reference_pep not in df.ID.unique():
                continue
            else:
                if len(ordered_protein_list) > 1:
                    df = df[~df['ID'].isin(ordered_protein_list[:-1])]
                dfs_containing_reference_pep.append(df)

        return dfs_containing_reference_pep

    @staticmethod
    def _update_ref_pep(ordered_protein_list):
        return ordered_protein_list[-1]

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
        # colors = cm.rainbow(np.linspace(0, 1, len(plot_data)))
        colors = self._get_color_map(color_data, plot_data)
        for i, coordinate in enumerate(plot_data):
            plt.scatter(coordinate[0], coordinate[1], c=colors[i], linewidth='0', s=50, vmin=-1, vmax=1,
                        label=prot_list[i], cmap='rainbow_r')
            # plt.scatter(coordinate[0], coordinate[1], color=colors, label=prot_list[i])

        plt.legend(bbox_to_anchor=(1.5, 1), loc=1, ncol=1)
        plt.show()

    ##############################################################################
    # FUNCTIONS TO RETRIEVE CLASS ATTRIBUTES: Mostly helper functions to get the #
    # shared data within the class.                                              #
    ##############################################################################

    ##### FILTER, ADD RANGE ####
    def _first_pass_filter_and_add_ranges(self):

        filtered_dfs = []
        for i in self.df_list:
            filt_ = i.loc[i['Affinity Level'] == 'High']
            clean = self._return_clean(filt_)
            filtered_dfs.append(clean)

        filtered_dfs_range_added = self._add_df_range(filtered_dfs)

        return filtered_dfs_range_added

    @staticmethod
    def _add_df_range(filtered_df):

        range_added = []

        for i in filtered_df:
            i['Final Pos'] = i['Pos'] + i['n-mer']
            i['Range'] = i[['Pos', 'Final Pos']].apply(lambda x: [list(range(x[0], x[1]))], axis=1)
            i = i.sort_values(by='Pos')
            range_added.append(i)

        return range_added

    @staticmethod
    def _return_clean(df):
        df['Peptide'] = df['Peptide'].str.replace('X', '-')
        df = df.loc[df['Peptide'].str.contains('--') == False]
        return df

