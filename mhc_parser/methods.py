import sys
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/')
from mhc_parser import net_mhc_func
import glob
import os
from shutil import move, rmtree


def run_mhc(pred, fasta_location, mhc_path):
    nmer = [pred.nmer]
    allele = [pred.allele]
    net_mhc = net_mhc_func.netMHCComand(mhc_path, fasta_location, nmers=nmer, alleles=allele)
    net_mhc.create_text_command(write_to_txt=True)
    net_mhc.run_netMHC()


def retrive_xls(pred, tmp_dir):

    dirs = glob.glob(tmp_dir + '/mhc_preds*')
    list_pred_dat = [pred.allele, str(pred.nmer), str(pred.original_position), pred.peptide, pred.protein.replace('_', '-')]

    for mhc_pred_dir in dirs:
        if all(dat in mhc_pred_dir for dat in list_pred_dat):
            pred_file = glob.glob(mhc_pred_dir + '/*.xls')[0]
            return pred_file


def write_to_fasta(pred, protein_id, tmp_dir):

    prot_dir = tmp_dir + protein_id
    file_name = "_".join([prot_dir, 'swap', pred.peptide, 'Pos', pred.original_position, 'ID',
                          pred.protein.replace('_', '-'), 'Allele', pred.allele, 'nmer', str(pred.nmer)])

    with open(file_name, "w") as fasta:
        for swap in pred.Swap.swaps:
            fasta.write("".join([">", pred.protein, "_", pred.peptide, "_", swap, "\n"]))
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


def filter_low_affinity(collection, protein_id, threshold):

    prediction_collection = collection[protein_id]
    below_threshold = []

    for prediction in prediction_collection:
        if prediction.affinity_level < threshold:
            below_threshold.append(prediction)
    return below_threshold


######

