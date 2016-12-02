import numpy as np
import pandas
from nepitope import pep_utils
from skbio import TabularMSA, Protein
import os
import webbrowser
from urllib import request
from subprocess import Popen, PIPE
import os


class Alignment(object):

    def __init__(self, fasta_file, ref_protein_id):
        """
        Maniputlation of alignment data. Works in conjunction with scikit-bio's TabulatMSA and Protein modules.
        :param msa_file: multiple sequence alignment file in fasta format (Clutal Omega recommended)
        :param ref_protein_file: Fasta file of reference protein
        """
        self.fasta = fasta_file
        self.project_dir = os.path.dirname(self.fasta)
        self.msa_file = self._create_and_viz_alignment()
        self.msa = self.read_msa_file()
        self.reference_protein_id = ref_protein_id
        self.reference_protein_string = self._get_ref_prot_from_id()
        self.positional_conservation = self._get_positional_conservation()

    def read_msa_file(self):

        msa = TabularMSA.read(self.msa_file, constructor=Protein)
        msa.reassign_index(minter='id')

        return msa

    def create_score_df_from_scikit_bio(self, nmers):
        """
        Function to generate a pandas dataframe containing information about how conserved each peptide in the
        reference protein is. Conservation scores are calculated for each nmer passed.
        is within a
        :param nmers: list of nmers of interest
        :return: dataframe with columns 'Score', 'Peptide', 'n-mer', that is the conservation score for each peptide
        identified
        """

        list_of_dfs = []

        for nmer in nmers:
            scores = []
            peptides = []

            for j in range(0, len(self.reference_protein_string) - nmer):

                scores.append(np.mean(self.positional_conservation[j:j + nmer]))   #Mean score for peptide
                peptides.append(self.reference_protein_string[j:j + nmer])

            df = pandas.DataFrame([scores, peptides], index=['Score', 'Peptide'])
            df = df.T
            df['n-mer'] = nmer
            list_of_dfs.append(df)

        return pandas.concat(list_of_dfs)

    def visualize_alignemnt(self):
        url_ = 'file:{}'.format(request.pathname2url(os.path.abspath(self.project_dir + '/alignment/' + 'MSA_easy_viewing.html')))
        webbrowser.open(url_)

    def _get_positional_conservation(self):
        """
        Apply metric to compute conservation for all alignment positions
        :return: conservation at each position, nan's replaced by zeros.
        """

        positional_conservation = self.msa.conservation(metric='inverse_shannon_uncertainty',
                                                        degenerate_mode='nan', gap_mode='include')
        return np.nan_to_num(positional_conservation)

    def _get_ref_prot_from_id(self):
        """
        Returns ref protein string from fasta
        """
        prot_id, prot_seqs = pep_utils.create_separate_lists(self.msa_file)
        prot_id = [prot.strip('>') for prot in prot_id]

        as_tuple = list(zip(prot_id, prot_seqs))

        ref_seq = None

        for tpl in as_tuple:
            if tpl[0] == self.reference_protein_id:
                ref_seq = tpl[1]
        if not ref_seq:
            raise ValueError('Protein ID provided not found in fasta file')
        else:
            return ref_seq

    def _create_and_viz_alignment(self):
        out_dir = os.path.dirname(self.fasta)

        if not os.path.isdir(out_dir + '/alignment'):
            os.mkdir(out_dir + '/alignment')

        out_align = out_dir + '/alignment' + '/MSA.fasta'
        if os.path.isfile(out_align):
            raise FileExistsError('Alignemnt already exists. Delete it or select other project location')

        self._create_fasta_and_html(out_align)

        return out_align

    def _create_fasta_and_html(self, out_align):

        process = Popen(['clustalo', '-i', self.fasta, '--residuenumber', '-o', out_align, '--outfmt=fasta'],
                        stdout=PIPE, stderr=PIPE)

        stdout, stderr = process.communicate()
        if not stderr:
            print('MSA in fasta created to %s' % out_align)
            self._create_html(out_align)
        else:
            print(stderr)

    @staticmethod
    def _create_html(out_dir):

        html_dir = os.path.dirname(out_dir) + '/MSA_easy_viewing.html'

        process = Popen(" ".join(['mview', '-in', 'fasta', '-ruler', 'on', '-html', 'head', '-coloring', 'any',
                                  out_dir, '>', html_dir]), shell=True)

        stdout, stderr = process.communicate()
        if not stderr:
            print('MSA in html created to %s' % html_dir)
        else:
            print(stderr)


class AddData (object):

    def __init__(self, msa_file_input, msa_file_output, scores_df, positional_conservation,
                 all_alleles=True, list_alleles=None, pos_cons_treshold=None):
        """

        :param msa_file_input:
        :param msa_file_output:
        :param scores_df:
        :param positional_conservation:
        :param all_alleles:
        :param list_alleles:
        :param pos_cons_treshold:
        """

        if pos_cons_treshold is None:
            self.pos_cons_treshold = 0.1
        else:
            self.pos_cons_treshold = pos_cons_treshold
        self.msa_file_input = msa_file_input
        self.msa_file_output = msa_file_output
        self.scores_df = scores_df
        self.all_alleles = all_alleles
        self.list_alleles = list_alleles
        self.positional_conservation = positional_conservation
        self.alleles = self._check_return_alleles(self.scores_df, self.all_alleles, self.list_alleles)
        self.nmers = self._get_nmers_from_affinity_df(self.scores_df)
        self.high_aa_low_cons_df = self._high_aff_low_cons_to_df(self.return_high_affinity_and_not_conserved())

    def _create_html(self):

        html_dir = os.path.dirname(self.msa_file_output) + '/MSA_easy_viewing.html'
        process = Popen(" ".join(['mview', '-in', 'fasta', '-ruler', 'on', '-html', 'head', '-coloring', 'any',
                                  self.msa_file_output, '>', html_dir]), shell=True)

        stdout, stderr = process.communicate()
        if not stderr:
            print('MSA in html created to %s' % html_dir)
        else:
            print(stderr)
        return html_dir

    def visualize_alignemnt(self):

        html_dir = self._create_html()
        url_ = 'file:{}'.format(request.pathname2url(html_dir))
        webbrowser.open(url_)

    def open_files(self):

        with open(self.msa_file_input) as inf, open(self.msa_file_output, 'w') as out:
            self.write_conservation_scores(inf, out)
            self.write_affinity_scores(out)

    def write_conservation_scores(self, inf, out):

        for line in inf:
            line = line.replace('X', '-')
            out.write(line)
        out.write('>CONSERVATION_INFO\n')

        for i in self.positional_conservation:
            if i > self.pos_cons_treshold:
                out.write('O')
            else:
                out.write('-')

    def write_affinity_scores(self, out):

        for nmer in self.nmers:
            for allele in self.alleles:

                to_print = self._slice_df(nmer, allele, self.scores_df)
                peps = self._get_peptides(to_print)

                for idx in range(0, len(peps)):

                    if idx > 3250:
                        continue

                    if '--' in peps[idx]:
                        continue

                    if not self._get_affinity_per_peptide(peps[idx], to_print):
                        continue
                    else:
                        self._write_out(nmer, allele, idx, out, peps)

    def high_affinity_low_cons_df(self):
        selected_df = self.scores_df.loc[(self.scores_df['Affinity Level'] == 'High') &
                                         (self.scores_df['Score'] < self.pos_cons_treshold)]

        selected_df = selected_df.loc[(selected_df['Pos'] < 3250) & (selected_df['Peptide'].str.contains('--') == False)]
        return selected_df

    def return_high_affinity_and_not_conserved(self):

        high_aff_not_cons = []

        for nmer in self.nmers:
            for allele in self.alleles:

                to_print = self._slice_df(nmer, allele, self.scores_df)
                peps = self._get_peptides(to_print)

                for idx in range(0, len(peps)):

                    mean_cons = self._get_mean_pos_cons_per_pep(nmer, idx)
                    if self._get_affinity_per_peptide(peps[idx], to_print):
                        if mean_cons < self.pos_cons_treshold:
                            print (mean_cons)
                            high_aff_not_cons.append([idx, peps[idx]])

        return high_aff_not_cons

    @staticmethod
    def _high_aff_low_cons_to_df(list_of_lists):
        return pandas.DataFrame(list_of_lists, columns=['Peptide Position', 'Peptide'])

    def _get_mean_pos_cons_per_pep(self, nmer, index):

        initial_aminoa_acid = index*nmer
        endind_amino_acid = (index+1)*nmer

        return np.mean(self.positional_conservation[initial_aminoa_acid:endind_amino_acid])

    @staticmethod
    def _write_out(nmer, allele, idx, out, peps):

        out.write('\n>High_Affinity_Loc|n-mer=%i|allele=%s\n' % (nmer, allele))
        out.write('-' * idx)
        out.write(peps[idx])
        out.write('-' * (len(peps) - idx - 1))

    @staticmethod
    def _get_affinity_per_peptide(pep, df):
        aff_per_pep = df.loc[df['Peptide'] == pep]
        if len(aff_per_pep) > 1:
            return False

        if list(aff_per_pep['Affinity Level'].values)[0] == 'High':
            return True
        else:
            return False

    @staticmethod
    def _slice_df(nmer, allele, df):

        to_print = df.loc[(df['n-mer'] == nmer) & (df['Allele'] == allele)]
        to_print['Peptide'] = to_print['Peptide'].str.replace('X', '-')
        return to_print

    @staticmethod
    def _get_peptides(df):
        return list(df['Peptide'].values)

    @staticmethod
    def _check_return_alleles(scores_df, all_alleles, list_alleles):

        if all_alleles:
            alls = list(scores_df.Allele.unique())
        else:
            alls = list_alleles

        if (all_alleles is False) & (list_alleles is None):
            raise ValueError('No allele provided')

        return alls

    @staticmethod
    def _get_nmers_from_affinity_df(scores_df):
        return list(scores_df['n-mer'].unique())

"""
class PyhloTree(object):

    def __init__(self):

        self.msa_file
"""
