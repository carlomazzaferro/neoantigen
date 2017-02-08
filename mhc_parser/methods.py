import sys
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/')
from mhc_parser import net_mhc_func
import glob
import os
import re
import numpy as np
from shutil import move, rmtree
import pandas


####### METHODS FOR RUNNING NETMHC ##########

#TODO: test netMHC interface methods

def run_mhc(pred, fasta_location, mhc_path):
    nmer = [pred.Nmer]
    allele = [pred.Allele]
    net_mhc = net_mhc_func.netMHCComand(mhc_path, fasta_location, nmers=nmer, alleles=allele)
    net_mhc.create_text_command(write_to_txt=True)
    net_mhc.run_netMHC()


def retrive_xls(pred, tmp_dir):

    dirs = glob.glob(tmp_dir + '/mhc_preds*')
    list_pred_dat = [pred.Allele, str(pred.Nmer), str(pred.Pos), pred.Peptide, pred.ID.replace('_', '-')]

    for mhc_pred_dir in dirs:
        if all(dat in mhc_pred_dir for dat in list_pred_dat):
            pred_file = glob.glob(mhc_pred_dir + '/*.xls')[0]
            return pred_file


def write_to_fasta(pred, protein_id, tmp_dir):
    print(protein_id, tmp_dir)
    prot_dir = tmp_dir + protein_id
    file_name = "_".join([prot_dir, 'swap', pred.Peptide, 'Pos', str(pred.Pos), 'ID',
                          pred.ID.replace('_', '-'), 'Allele', pred.Allele, 'Nmer', str(pred.Nmer)])

    with open(file_name, "w") as fasta:
        for swap in pred.Swap.swaps:
            fasta.write("".join([">", pred.ID, "_", pred.Peptide, "_", swap, "\n"]))
            fasta.write("".join([swap, '\n']))

    return file_name


def reorg_files(tmp_dir):

    dirs = glob.glob(tmp_dir + '/mhc_preds*')
    final_dest = tmp_dir + '/preds_per_swap'
    try:
        os.mkdir(final_dest)
    except:
        pass
    for pred_dir in dirs:
        file_source = glob.glob(pred_dir + '/*.xls')
        move(file_source[0], final_dest)
    print('Swap predictions regrouped to %s' % final_dest)
    for i in dirs:
        rmtree(i)

######## OTHER METHODS ########


def filter_low_affinity(collection, protein_id, threshold):

    prediction_collection = collection[protein_id]['Predictions']
    below_threshold = []

    for prediction in prediction_collection:
        if prediction.nM < threshold:
            below_threshold.append(prediction)
    return below_threshold


def to_df(dictionary_collection, cols, protein_id):
    protein_data = dictionary_collection[protein_id]
    pred_data = []

    for dat in protein_data['Predictions']:
        elemnts = []
        for elem in dat.data_list:
            elemnts.append(elem)
        elemnts.append(dat.Nmer)
        pred_data.append(elemnts)

    df = pandas.DataFrame(pred_data, columns=cols)
    df.Pos = df.Pos.apply(lambda x: int(x))

    return df


def filter_all(collection, threshold):

    new_collec = {prot: {} for prot in list(collection.keys())}

    for key_ in collection.keys():
        new_collec[key_]['High Affinity Ranges'] = collection[key_]['High Affinity Ranges']
        new_collec[key_]['ProtSeq'] = collection[key_]['ProtSeq']
        new_collec[key_]['Predictions'] = filter_low_affinity(collection, key_, threshold)

    return new_collec

########## Pairwise comparison stuff ############
class PairWiseComp(object):

    def __init__(self, pred_col, threshold, min_nmer):

        self.pred_col = pred_col.dictionary_collection
        self.filtered_col = filter_all(self.pred_col, threshold)
        self.min = min_nmer
        self.proteins = list(self.pred_col.keys())

    def pipe_run(self):
        """
        The actual pipeline for pairwise comparisons.
        :return:
        """
        list_df_scores = []

        for protein in self.proteins:

            dic = self.run_pairwise_comp(protein)
            self._add_total_score(dic)

            list_df_scores.append(pandas.DataFrame(dic).T)

        return list_df_scores

    def get_comp_dic(self, protein):
        tpl = [np.array([protein] * len(self.proteins)), np.array(self.proteins)]
        tpls = list(zip(*tpl))

        dics = {}

        for pair in tpls:
            dics[pair] = {num: 0 for num in list(range(self.min, 12))}
            dics[pair]['Num High AA'] = 0
            dics[pair]['Matches Loc'] = []

        return dics

    def run_pairwise_comp(self, protein):
        """
        Pairwise comparison ran between a ref df and every over df. Every df will be the ref df at some point.
        Match scores are saved to a dictionary, which is updated at every round.
        :param protein: reference dataframe used to make comparisons with every other protein sequence
        :return: a dictionary containing score data for a single pairwise comparison.
        """

        # List of lists
        skel_dic = self.get_comp_dic(protein)
        ref_pred_collection = self.filtered_col[protein]['Predictions']

        for prot_name in self.proteins:

            prot_seq = self.filtered_col[prot_name]['ProtSeq']
            ranges = self.filtered_col[prot_name]['High Affinity Ranges']

            # Ranges: index data about the location of high affinity peptides in protein being used for comparison
            # Ranges_2: make shallow list from deep list of lists
            ranges_2 = [item for sublist in [i[0] for i in ranges] for item in sublist]

            matches_range = []

            for pred in ref_pred_collection:

                pep = pred.Peptide
                high_aa_count = 0
                pep_len = len(pep)
                count = prot_seq.count(pep)  # Number of times a single pep occurs in the entire prot seq

                if count > 0:  # Find locations where matches occur
                    it = re.finditer(pep, prot_seq)

                    for i in it:
                        present_range = list(range(i.start(), i.end()))
                        if set(present_range).issubset(set(ranges_2)):
                            high_aa_count += 1
                            matches_range.append(present_range)  # Retain match location data

                self._update_dict_values_per_len(skel_dic, protein, prot_name, count,
                                                 pep_len, high_aa_count, matches_range)

        return skel_dic

    def get_peps(self, protein):
        peps = []
        collec = self.filtered_col[protein]['Predictions']

        for item in collec:
            peps.append(item.Peptide)

        return peps

    @staticmethod
    def _update_dict_values_per_len(score_dict_per_len, ref_prot, prot_name, num_matches,
                                    pep_len, high_aa_matches, matches_range):

        tpl_key = (ref_prot, prot_name)
        score_dict_per_len[tpl_key][pep_len] += num_matches
        score_dict_per_len[tpl_key]['Num High AA'] += high_aa_matches
        if matches_range not in score_dict_per_len[tpl_key]['Matches Loc']:  # To eliminate duplicates
            score_dict_per_len[tpl_key]['Matches Loc'].append(matches_range)

    def _add_total_score(self, score_dict_len):
        for key in score_dict_len:
            score_dict_len[key]['Total'] = sum(self.excl_list(score_dict_len[key].values())) - \
                                           score_dict_len[key]['Num High AA']

        return score_dict_len

    @staticmethod
    def excl_list(it):
        return [i for i in it if type(i) is int]