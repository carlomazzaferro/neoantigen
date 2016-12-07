import pandas
import re
import sys
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser')
import utilities


class PairwiseComp(object):

    threshold_levels = {'High': [0, 50],
                        'Intermediate': [50, 500],
                        'Low': [500, 5000],
                        'No': [5000, 100000]}

    def __init__(self, pred_collection, min_nmer, original_fasta, threshold=None):

        """
        :param list_dfs: list of dataframes coming from previous step
        :param min_nmer: minimum length nmer for which matches will be found
        :param original_fasta: fasta file containing original protein sequences
        :param threshold: #######Just High for now - need to be properly set######
        """

        self.pred_col = pred_collection
        self.threshold = self._return_threshold_level(threshold)
        self.filt_dfs = self._first_pass_filter_and_add_ranges()
        self.min_nmer = min_nmer
        self.fasta = original_fasta
        self.peps_and_prots = self._get_id_protseq_tuple()
        self.peps_prots_ranges = self._add_high_aa_location_data()
        self.list_of_scores_dicts = self._get_protein_dict(self.filt_dfs)
        self.original_proteins_df = pandas.DataFrame(self.peps_prots_ranges, columns=['Protein', 'ProtSeq', 'Ranges'])
        self.original_proteins = self.original_proteins_df.Protein.values.tolist()

    def pipe_run(self):
        """
        The actual pipeline for pairwise comparisons.
        :return:
        """

        list_dict_scores = []

        for ref_df in self.filt_dfs:
            """
            Run only for high affinity peptides: filt.dfs are dfs from original proteins that have been filtered by a
            user-set threshold level.
            """
            list_dict_scores.append([self.run_pairwise_comp(ref_df), ref_df.Protein.unique()])

        total_score_dict = self._add_total_score(list_dict_scores)
        return self._return_df_from_dict_of_dicts(total_score_dict)

    def run_pairwise_comp(self, ref_df):
        """
        Pairwise comparison ran between a ref df and every over df. Every df will be the ref df at some point.
        Match scores are saved to a dictionary, which is updated at every round.
        :param ref_df: reference dataframe used to make comparisons with every other protein sequence
        :return: a dictionary containing score data for a single pairwise comparison.
        """

        #List of lists
        ref_df_peps = self._get_all_peptides_from_df(ref_df)                                 #Extract high affinity peptides
        score_dict_per_len = self._get_protein_dict_per_len(self.filt_dfs, ref_df_peps)      #Create scores dictionary

        for prot_name in self.original_proteins:

            prot_seq = self.original_proteins_df.ProtSeq[self.original_proteins_df.Protein == prot_name].values[0]
            ranges = self.original_proteins_df.Ranges[self.original_proteins_df.Protein == prot_name].values[0]
            #Ranges: index data about the location of high affinity peptides in protein being used for comparison
            #Ranges_2: make shallow list from deep list of lists
            ranges_2 = [item for sublist in [i[0] for i in ranges] for item in sublist]

            matches_range = []

            for list_pep in ref_df_peps:
                for single_pep in list_pep:

                    high_aa_count = 0
                    pep_len = len(single_pep)
                    count = prot_seq.count(single_pep)    #Number of times a single pep occurs in the entire prot seq

                    if count > 0:                         #Find locations where matches occur
                        it = re.finditer(single_pep, prot_seq)

                        for i in it:
                            present_range = list(range(i.start(), i.end()))
                            if set(present_range).issubset(set(ranges_2)):
                                high_aa_count += 1
                                matches_range.append(present_range)      #Retain match location data

                    self._update_dict_values_per_len(score_dict_per_len, prot_name, count,
                                                     pep_len, high_aa_count, matches_range)

        return score_dict_per_len

    def _return_df_from_dict_of_dicts(self, total_score_dict):

        """
        Rearrange pandas df as a MultiIndex df:
         ______________ _______________ _______ _______ _______ ________ _________________ ___________
        |              |               |   4   |   5   |  ...  |   11   |   Matches Loc   |   Total   |
        | Outer Prot 1 | Inner Prot 1  |   #   |   #   |  ...  |   #    |       []        |     #     |
        | Outer Prot 1 | Inner Prot 2  |   #   |   #   |  ...  |   #    |       []        |     #     |
        |  ...
        | Outer Prot 1 | Inner Prot n  |   #   |   #   |  ...  |   #    |       []        |     #     |
        | Outer Prot 2 | Inner Prot 1  |   #   |   #   |  ...  |   #    |       []        |     #     |
        | Outer Prot 2 | Inner Prot 2  |   #   |   #   |  ...  |   #    |       []        |     #     |
        |  ...
        | Outer Prot 2 | Inner Prot n  |   #   |   #   |  ...  |   #    |       []        |     #     |
        |  ...
        | Outer Prot n | Inner Prot n  |   #   |   #   |  ...  |   #    |       []        |     #     |
         -------------- --------------- ------- ------- ------- -------- ----------------- -----------
        """

        lsls = []
        for i in total_score_dict:
            reform = {(outerKey, str(innerKey)): values for outerKey, innerDict in i[0].items() for innerKey, values in
                      innerDict.items()}
            lsls.append(pandas.Series(reform, name=i[-1][0]))
        df = pandas.concat(lsls, axis=1).sort_index(axis=1)
        dc = self._rank_by_tuple(df)

        return pandas.DataFrame(dc).T

    @staticmethod
    def _rank_by_tuple(df):

        running_dict = {}
        cols = list(df.columns.values)

        for col_idx in cols:
            for col_col in cols:
                vals = df[col_col].loc[col_idx]
                running_dict[(col_idx, col_col)] = vals    #MultiIndex df

        return running_dict

    def write_csv_filtered_dfs(self, path):

        for i in self.filt_dfs:
            i.to_csv(path + i.Protein.unique()[0] + '.csv', sep=',')

    def _add_high_aa_location_data(self):

        new_list_ranges_added = []
        for tpl in self.peps_and_prots:
            for df in self.filt_dfs:
                if df.Protein.unique() == tpl[0]:
                    new_list_ranges_added.append(tpl + (list(df.Range.values), ))

        return new_list_ranges_added

    def _first_pass_filter_and_add_ranges(self):
        """
        Retrieve only high affinity locations, and add the range of the peptide
        :return: list of dfs
        """

        filtered_dfs = []
        for i in self.list_dfs:

            filt_ = i.loc[(i['Affinity_level'].apply(lambda x: float(x)) > self.threshold[0]) &
                          (i['Affinity_level'].apply(lambda x: float(x)) < self.threshold[1])]

            clean = self._return_clean(filt_)
            filtered_dfs.append(clean)

        filtered_dfs_range_added = self._add_df_range(filtered_dfs)

        return filtered_dfs_range_added

    def _get_id_protseq_tuple(self):

        peps, prots = utilities.create_separate_lists(self.fasta)
        peps_short_named = [pep.strip('>').replace('.', '_')[0:15] for pep in peps]

        return list(zip(peps_short_named, prots))

    @staticmethod
    def _add_df_range(filtered_df):

        range_added = []

        for i in filtered_df:
            i['Final Pos'] = i['Pos'].apply(lambda x: int(x)) + i['nmer'].apply(lambda x: int(x)) -1
            i['Range'] = i[['Pos', 'Final Pos']].apply(lambda x: [list(range(x[0], x[1] + 1))], axis=1)
            i = i.sort_values(by='Pos')
            range_added.append(i)

        return range_added

    @staticmethod
    def _return_clean(df):
        df['Peptide'] = df['Peptide'].str.replace('X', '-')
        df = df.loc[df['Peptide'].str.contains('--') == False]
        return df

    def _add_total_score(self, score_dict_len):

        for dict_ in score_dict_len:
            for key in dict_[0]:
                dict_[0][key]['Total'] = sum(self.excl_list(dict_[0][key].values())) - dict_[0][key]['Num High AA']

        return score_dict_len

    @staticmethod
    def excl_list(it):
        return [i for i in it if type(i) is int]


    @staticmethod
    def _get_protein_dict_per_len(lsss_1, ref_df_peps):

        all_peps_set_list = [item for sublist in ref_df_peps for item in sublist]
        max_len = len(max(all_peps_set_list, key=len))
        possible_pep_len = list(range(5, max_len + 1))

        dict_of_dicts = {}

        for df in lsss_1:
            protein = df.Protein.unique()[0]
            dict_of_dicts[protein] = {num: 0 for num in possible_pep_len}
            dict_of_dicts[protein]['Num High AA'] = 0
            dict_of_dicts[protein]['Matches Loc'] = []

        return dict_of_dicts


    @staticmethod
    def _get_protein_dict(lsss_1):

        proteins = []
        for i in lsss_1:
            proteins.append(i.Protein.unique()[0])

        return {prot: 0 for prot in proteins}

    def _get_all_peptides_from_df(self, df):

        peptides_in_df = df.Peptide.values.tolist()
        all_possible_peps = []

        for i in peptides_in_df:
            all_possible_peps.extend(self._get_possible_peptides(i))

        return all_possible_peps

    def _get_possible_peptides(self, stri):
        # TODO: Clean this up asap. Too hacky, not very readable

        nmer_list = list(range(self.min_nmer, len(stri) + 1))
        length_pep = len(stri)
        num_possible_sub_peps = [length_pep - nmer + 1 for nmer in nmer_list if length_pep - nmer + 1 > 0]

        subs_all = []
        for i in num_possible_sub_peps:
            for nmer, j in zip(nmer_list, range(i)):
                subs = [stri[j:j + nmer] for i, nmer in enumerate(nmer_list)]
                subs_all.extend(subs)

        subs_all = list(set(subs_all))
        lst = [[w for w in subs_all if len(w) == num] for num in set(len(i) for i in subs_all)]

        return lst

    @staticmethod
    def _update_dict_values_per_len(score_dict_per_len, prot_name, num_matches, pep_len, high_aa_matches, matches_range):

        score_dict_per_len[prot_name][pep_len] += num_matches
        score_dict_per_len[prot_name]['Num High AA'] += high_aa_matches
        if matches_range not in score_dict_per_len[prot_name]['Matches Loc']:    #To eliminate duplicates
            score_dict_per_len[prot_name]['Matches Loc'].append(matches_range)

    def _return_threshold_level(self, threshold):
        if not threshold:
            return [0, 50]
        if isinstance(threshold, str):
            return self.threshold_levels[threshold]
        if isinstance(threshold, int):
            return [0, threshold]
