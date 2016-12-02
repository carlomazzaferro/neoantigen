import glob
import sys
import models

if __name__ == '__main__':
    ref = 'S__pyogenes_Cas'

    files = glob.glob('/Users/carlomazzaferro/Desktop/New_General/MICE/mhc_preds_fasta_base_new_prots/*.xls')
    fasta_file = '/Users/carlomazzaferro/Desktop/New_General/MICE/fasta_base_new_prots.fasta'
    pred_col = models.PredictionCollection(files, fasta_file)
    parsed = pred_col.digest_multiple()
    """
    letters = ['L', 'K', 'M', 'P', 'S']
    pep = 'ABCDEFGH'
    ls_pep = list(pep)
    print(ls_pep)

    print(pep)
    for i, let in enumerate(ls_pep):
        print [[letter,   ]
    """

