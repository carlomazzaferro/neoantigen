import sys
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/')
from mhc_parser import net_mhc_func
import glob
import os
import re
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

########## Pairwise comparison stuff ############

def pipe_run(prediction_collection):
    """
    The actual pipeline for pairwise comparisons.
    :return:
    """
    proteins = list(prediction_collection.dictionary_collection.keys())
    #for i in

    #prediction_per_protein =

    list_dict_scores = []

    for ref_df in a:
        """
        Run only for high affinity peptides: filt.dfs are dfs from original proteins that have been filtered by a
        user-set threshold level.
        """
        list_dict_scores.append([run_pairwise_comp(ref_df), ref_df.Protein.unique()])

    total_score_dict = self._add_total_score(list_dict_scores)
    return self._return_df_from_dict_of_dicts(total_score_dict)


def run_pairwise_comp(self, ref_df):
    """
    Pairwise comparison ran between a ref df and every over df. Every df will be the ref df at some point.
    Match scores are saved to a dictionary, which is updated at every round.
    :param ref_df: reference dataframe used to make comparisons with every other protein sequence
    :return: a dictionary containing score data for a single pairwise comparison.
    """

    # List of lists
    ref_df_peps = self._get_all_peptides_from_df(ref_df)  # Extract high affinity peptides
    score_dict_per_len = self._get_protein_dict_per_len(self.filt_dfs, ref_df_peps)  # Create scores dictionary

    for prot_name in self.original_proteins:

        prot_seq = self.original_proteins_df.ProtSeq[self.original_proteins_df.Protein == prot_name].values[0]
        ranges = self.original_proteins_df.Ranges[self.original_proteins_df.Protein == prot_name].values[0]
        # Ranges: index data about the location of high affinity peptides in protein being used for comparison
        # Ranges_2: make shallow list from deep list of lists
        ranges_2 = [item for sublist in [i[0] for i in ranges] for item in sublist]

        matches_range = []

        for list_pep in ref_df_peps:
            for single_pep in list_pep:

                high_aa_count = 0
                pep_len = len(single_pep)
                count = prot_seq.count(single_pep)  # Number of times a single pep occurs in the entire prot seq

                if count > 0:  # Find locations where matches occur
                    it = re.finditer(single_pep, prot_seq)

                    for i in it:
                        present_range = list(range(i.start(), i.end()))
                        if set(present_range).issubset(set(ranges_2)):
                            high_aa_count += 1
                            matches_range.append(present_range)  # Retain match location data

                self._update_dict_values_per_len(score_dict_per_len, prot_name, count,
                                                 pep_len, high_aa_count, matches_range)

    return score_dict_per_len






