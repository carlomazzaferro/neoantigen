import os
import glob
from shutil import move, rmtree
from nepitope import net_MHC_interface
import importlib
importlib.reload(net_MHC_interface)


class Swaps(object):

    list_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def __init__(self, high_affinity_df, fasta_file_dir, net_mhc_path, proteins=None):

        self.df = high_affinity_df
        self.protein_ids = self._get_prot_ids(proteins)
        #self.protein_input_df = self.df[self.df.ID.isin(self.protein_ids)]
        self.fasta_dir = fasta_file_dir
        self.mhc_path = net_mhc_path

    def find_swaps_write_to_fasta(self):

        nmers = self._get_nmers(self.df)
        alleles = self._get_alleles(self.df)

        mhc_commands = []
        for prot_id in self.protein_ids:
            try:
                os.mkdir(self.fasta_dir + '%s/' % prot_id)
            except:
                pass
            for nmer in nmers:
                for allele in alleles:
                    sliced = self._slice_df(nmer, allele, prot_id)
                    if self.check_size(sliced):
                        list_of_lists = sliced[['n-mer', 'Allele', 'ID', 'Pos', 'Peptide']].values.tolist()

                        for item in list_of_lists:
                            swaps = self._create_swaps(item[-1])
                            fasta_file = self._open_write_fasta(item, swaps, prot_id)
                            self._create_mhc_command(item, fasta_file)

            self.reorg_files(prot_id)

        return mhc_commands

    def reorg_files(self, prot_id):

        prot_dir = self.fasta_dir + '%s' % prot_id
        dirs = glob.glob(prot_dir + '/mhc_preds*')
        final_dest = prot_dir + '/preds_per_swap'
        try:
            os.mkdir(final_dest)
        except:
            pass
        for i in dirs:
            file_source = glob.glob(i + '/*.xls')
            move(file_source[0], final_dest)
        print('Swap predictions regrouped to %s' % final_dest)
        for i in dirs:
            rmtree(i)

    def _create_mhc_command(self, item, fasta_location):
        nmer = [item[0]]
        allele = [item[1]]
        net_mhc = net_MHC_interface.netMHCComand(self.mhc_path, fasta_location, nmers=nmer, alleles=allele)
        net_mhc.create_text_command(write_to_txt=True)
        net_mhc.run_netMHC()

    def _open_write_fasta(self, data, swaps, prot_id):
        file_name = "_".join([self.fasta_dir + '%s/' % prot_id, 'swap', data[-1], 'Pos', str(data[-2]), 'ID', str(data[-3]).replace('_', '-'),
                              'Allele', str(data[-4]), 'nmer', str(data[-5])])

        with open(file_name + '.fasta', 'w') as inf:
            for swap in swaps:
                inf.write("".join(['>', prot_id, '_', swap, '\n']))
                inf.write(swap + '\n')

        return file_name + '.fasta'

    def _create_swaps(self, peptide):

        list_peps = []
        for i in range(len(peptide)):
            for k in range(len(self.list_AA)):
                list_peps.append(self._insert_aa(peptide, i, self.list_AA[k]))

        return list_peps

    def _slice_df(self, nmer, allele, prot_id):
        return self.df.loc[(self.df['n-mer'] == nmer) & (self.df['Allele'] == allele) & (self.df['ID'] == prot_id)]

    @staticmethod
    def _insert_aa(string, index, aa):
        hash_string = list(string)
        del hash_string[index]
        hash_string.insert(index, aa)
        return "".join(hash_string)

    def _get_prot_ids(self, proteins):
        if proteins == 'All':
            return list(self.df['ID'].unique())
        if isinstance(proteins, list):
            return self.check_existence(proteins)

    def check_existence(self, proteins):
        for protein in proteins:
            if protein not in self.df.ID.unique():
                raise ValueError('Input protein %s not found in csv files' % protein)
        return proteins

    @staticmethod
    def _get_nmers(pepdata):
        return pepdata['n-mer'].unique()

    @staticmethod
    def _get_alleles(pepdata):
        return pepdata['Allele'].unique()

    @staticmethod
    def check_size(sliced):
        if len(sliced) == 0:
            return False
        else:
            return True
