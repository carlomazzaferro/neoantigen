from collections import Counter
import pandas
import re
from nepitope import pep_utils

class MergeSwapsAndComp(object):

    reference = 'StreptococcusPyogenes_reference'
    num_common_peps = 50
    no_return_cols = ['Most Comm',	'issubset' , 'Final Filter']

    def __init__(self, pw_comp_df, swaps_df, num_peptide_swaps):

        self.n = num_peptide_swaps
        self.original_df = pw_comp_df
        self.ref_df = self._get_ref_df()
        self.swaps_df = self._get_swaps_df(swaps_df)
        self.loc_matches = self._get_loc_matches()
        self.top_loc_matches = self._Most_Common(self.loc_matches)
        self.priority_df = self.get_priority_df()
        self.swap_list = self.get_priority_swap_list()
        self.top_swap_df = self.select_top()


    ##### METHODS AND THEIR HELPERS ######
    def get_modified_fasta(self, fasta, fasta_out_dir):

        idx, seq = pep_utils.create_separate_lists(fasta)
        zipped = list(zip(idx, seq))
        zipped = [list(x) for x in zipped]
        swapped_fasta = []

        peps_to_swap = self._get_swap_peps()
        orig_peps = self._get_orig_peps()
        orig_and_swaps = dict(zip(orig_peps, peps_to_swap))

        for i in zipped:
            if i[0][1::] == self.reference:
                i[1] = self.replace_all(i[1], orig_and_swaps)
            swapped_fasta.append(i)

        with open(fasta_out_dir, 'w') as out:
            for i in swapped_fasta:
                out.write(i[0] + '\n')
                out.write(i[1] + '\n')

    @staticmethod
    def replace_all(text, dic):
        for i, j in dic.items():
            text = text.replace(i, j)
        return text

    def get_modified_df_list(self, list_df):

        peps_to_swap = self._get_swap_peps()
        orig_peps = self._get_orig_peps()
        orig_and_swaps = dict(zip(orig_peps, peps_to_swap))
        print(orig_and_swaps)

        for idx, df in enumerate(list_df):
            if df.ID.unique()[0] == self.reference:
                new_df = df
                new_df.Peptide = new_df.Peptide.replace(orig_and_swaps)
                ix = new_df.loc[new_df.Peptide.isin(peps_to_swap)].index
                print(ix)
                new_df.set_value(ix, 'Affinity Level', 'No')
                list_df[idx] = new_df

        return list_df

    def _get_swap_peps(self):

        swaps = []
        list_swaps = list(self.top_swap_df['top scoring peptides'][0:self.n].values)
        for i in list_swaps:
            swap = list(filter(None, re.sub('[^A-Za-z0-9]+', ',', i).split(',')))[-3]
            swaps.append(swap)

        return swaps

    def _get_orig_peps(self):
        return self.top_swap_df['original peptide'].values.tolist()[0:self.n]

    def final_sort_and_filtering(self, df):
        unique_priority_list = []

        for top_loc_match in self.top_loc_matches:
            sel, df = self._select_top(df, top_loc_match[0])
            unique_priority_list.append(sel)
        final_df = pandas.concat(unique_priority_list)[0:self.n]   #Number of wanted peptides
        final_df = final_df.drop(self.no_return_cols, axis=1)      #Unwanted, used as placeholders

        return final_df

    @staticmethod
    def _select_top(top_swap, top_loc_match):
        top_swap['Final Filter'] = [top_loc_match in i for i in top_swap.Range]
        return top_swap[top_swap['Final Filter'] == True], top_swap[top_swap['Final Filter'] == False]

    #######-------#######------#######-------#######

    def get_priority_df(self):                                  ## ATTRIBUTE METHODS
        df = pandas.concat(self.get_priority_swap_list())
        return df.groupby(df.index).first()

    def get_priority_swap_list(self):                           ## ATTRIBUTE METHODS
        list_dfs = []
        for i in self.top_loc_matches:
            list_dfs.append(self.filter_non_most_commons(i[0]))
        return list_dfs

    def filter_non_most_commons(self, idx):
        self.swaps_df['Most Comm'] = self.swaps_df.Range.apply(lambda x: idx in x)
        return self.swaps_df[self.swaps_df['Most Comm'] == True]

    def select_top(self):                             ## ATTRIBUTE METHODS

        maximal_list = []
        for df in self.swap_list:
            ranges = df.Range.values
            df['issubset'] = df['original pos'].apply(lambda x: self.count_occurences(x, ranges, ))
            maximal_list.append(df.ix[df['issubset'].idxmax()])
        df = pandas.DataFrame(maximal_list)

        return df.groupby(df.index).first()

    @staticmethod
    def count_occurences(pos, ranges):                               ## ATTRIBUTE METHODS
        count = 0
        for i in ranges:
            if pos in i:
                count += 1
        return count

    def _get_ref_df(self):                                   ## ATTRIBUTE METHOD
        ref_df = self.original_df.loc[self.reference]
        return ref_df

    def _get_swaps_df(self, swaps_df):                       ## ATTRIBUTE METHODS
        swaps_df['Range'] = swaps_df.apply(lambda x: self._add_ranges_list(x['original pos'],
                                                                           x['original pos'] + x['nmer']), axis=1)
        return swaps_df

    @staticmethod
    def _add_ranges_list(a, b):                                ## ATTRIBUTE METHODS
        return list(range(a, b))

    def _get_loc_matches(self):                                  ## ATTRIBUTE METHODS
        return self._flatten(self.ref_df['Matches Loc'].values.tolist())

    def _flatten(self, x):                                        ## ATTRIBUTE METHODS
        result = []
        for el in x:
            if hasattr(el, "__iter__") and not isinstance(el, str):
                result.extend(self._flatten(el))
            else:
                result.append(el)
        return result

    def _Most_Common(self, lst):                                    ## ATTRIBUTE METHODS
        data = Counter(lst)
        return data.most_common(self.num_common_peps)
