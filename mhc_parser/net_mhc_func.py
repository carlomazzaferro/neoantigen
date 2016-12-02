import os
from subprocess import Popen, PIPE


class netMHCComand(object):

    all_supertypes = ['HLA-A0101', 'HLA-A0201', 'HLA-A0301', 'HLA-A2402',
                      'HLA-A2601', 'HLA-B0702', 'HLA-B0801', 'HLA-B2705',
                      'HLA-B3901', 'HLA-B4001', 'HLA-B5801', 'HLA-B1501']

    default_nmers = [8, 9, 10, 11]

    def __init__(self, netMHC_path, fasta_file, nmers=None, alleles=None):

        self.fasta = fasta_file
        self.basename = os.path.basename(os.path.splitext(self.fasta)[0])
        self.base_dir = os.path.dirname(self.fasta)
        self.mhc_pred_dir = self.base_dir + '/mhc_preds_' + self.basename
        self.nmers = self._get_nmers(nmers)
        self.alleles = self._get_alleles(alleles)
        self.netMHC = netMHC_path

    def _get_alleles(self, alleles):
        """
        Determine if allele input is provided. If None is passed, return all available supertypes.
        :param alleles: input
        :return: alleles
        """

        if alleles == None:
            self.alleles = self.all_supertypes
        else:
            self.alleles = alleles
        return self.alleles

    def _get_nmers(self, nmers):
        """
        Determine if nmer input is provided. If None is passed, return all default alleles.
        :param alleles: input
        :return: alleles
        """

        if nmers == None:
            self.nmers = self.default_nmers
        else:
            self.nmers = nmers
        return self.nmers

    def _det_if_multi_sequence(self):
        """
        Determine if netMHC will be run in a single sequence or in multiple sequences.
        :return: True (multi seq) or False (single seq) in fasta.
        """
        ids, prot_seqs = utilities.create_separate_lists(self.fasta)
        return len(ids) > 1

    def create_text_command(self, write_to_txt=True):
        try:
            os.mkdir(self.mhc_pred_dir)
        except:
            pass

        command_str = []

        for i, allele in enumerate(self.alleles):
            for j, nmer in enumerate(self.nmers):
                file_name_out = self.mhc_pred_dir + '/' + self.basename + '_' + allele + '_' + str(nmer) + '.xls'
                command_str.append('%s -a %s -f %s -l %i -xls -xlsfile %s' % (self.netMHC, allele, self.fasta,
                                                                              nmer, file_name_out))
        if write_to_txt:
            self._write_text_file(self.mhc_pred_dir, command_str)
            return 'netMHC commands written to run_netMHC.txt located at %s' %self.mhc_pred_dir
        else:
            return command_str

    def run_netMHC(self):
        process = Popen(" ".join(['bash', self.mhc_pred_dir + '/run_netMHC.txt']), shell=True)
        stdout, stderr = process.communicate()

        if not stderr:
            print('Predictions being saved to %s' % self.mhc_pred_dir)
        else:
            print(stderr)


    @staticmethod
    def _write_text_file(out_location, command_list):

        with open(out_location + '/run_netMHC.txt', 'w') as out:
            for i in command_list:
                out.write(i + '\n')

