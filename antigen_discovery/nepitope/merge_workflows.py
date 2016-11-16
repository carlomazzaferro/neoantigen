from collections import Counter
import pandas

class MergeSwapsAndComp(object):

    reference = 'StreptococcusPyogenes_reference'

    def __init__(self, pw_comp_df, swaps_df, num_peptide_swaps):

        self.n = num_peptide_swaps
        self.original_df = pw_comp_df
        self.ref_df = self._get_ref_df()
        self.swaps_df = self._get_swaps_df(swaps_df)
        self.loc_matches = self._get_loc_matches()
        self.top_loc_matches = self._Most_Common(self.loc_matches)
        self.priority_df = self.get_priority_df()
        self.swap_list = self.get_priority_swap_list()
        self.top_swap_list = self.select_top(self.swap_list)

    def get_priority_df(self):
        df = pandas.concat(self.get_priority_swap_list())
        return df.groupby(df.index).first()

    def get_priority_swap_list(self):
        list_dfs = []
        for i in self.top_loc_matches:
            list_dfs.append(self.filter_non_most_commons(i[0]))
        return list_dfs

    def filter_non_most_commons(self, idx):
        self.swaps_df['Most Comm'] = self.swaps_df.Range.apply(lambda x: idx in x)
        return self.swaps_df[self.swaps_df['Most Comm'] == True]

    def select_top(self, list_dfs):
        maximal_list = []
        for df in list_dfs:
            ranges = df.Range.values
            df['issubset'] = df['original pos'].apply(lambda x: self.count_occurences(x, ranges, ))
            maximal_list.append(df.ix[df['issubset'].idxmax()])
        df = pandas.DataFrame(maximal_list)

        return df.groupby(df.index).first()

    @staticmethod
    def count_occurences(pos, ranges):
        count = 0
        for i in ranges:
            if pos in i:
                count += 1
        return count

    def _get_ref_df(self):
        ref_df = self.original_df.loc[self.reference]
        return ref_df

    def _get_swaps_df(self, swaps_df):
        swaps_df['Range'] = swaps_df.apply(lambda x: self._add_ranges_list(x['original pos'],
                                                                           x['original pos'] + x['nmer']), axis=1)
        return swaps_df

    @staticmethod
    def _add_ranges_list(a, b):
        return list(range(a, b))

    def _get_loc_matches(self):
        return self._flatten(self.ref_df['Matches Loc'].values.tolist())

    def _flatten(self, x):
        result = []
        for el in x:
            if hasattr(el, "__iter__") and not isinstance(el, str):
                result.extend(self._flatten(el))
            else:
                result.append(el)
        return result

    def _Most_Common(self, lst):
        data = Counter(lst)
        return data.most_common(self.n)



