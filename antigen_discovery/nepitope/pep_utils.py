import pandas
from nepitope import scoring_utils
import os
#############################################################################
# UTILITIES TO RETRIEVE ORTHOLOGUES HAVING HIGH AFFINITY FROM MSA ALIGNMENT #
#############################################################################


def retrieve_orthologue_peptides(fasta_file_msa, high_aa_low_cons_df):
    list_of_pos_and_nmers = _retrieve_relevant_positions(high_aa_low_cons_df)
    print (list_of_pos_and_nmers[0])
    pepdata = _iterate_over_msa_file(list_of_pos_and_nmers, fasta_file_msa)
    return _deep_list_to_df(pepdata)


def _iterate_over_msa_file(list_of_pos_and_nmers, fasta_file_msa):
    peps_pos_allele = []
    list_ids, list_proteins = scoring_utils.Score.create_separate_lists(fasta_file_msa)

    for i in range(len(list_proteins) - 1):
        for position in list_of_pos_and_nmers:
            allele = position[2]
            nmer = position[1]
            initial_amino_acid = position[0]
            final_amino_acid = position[0] + nmer
            peps_pos_allele.append([list_proteins[i][initial_amino_acid:final_amino_acid],
                                       nmer, list_ids[i], initial_amino_acid, allele])

    return peps_pos_allele


def write_fasta_targets_out(pepdata, fasta_file_output):
    pepdata_df = _return_ids_peps_df(pepdata)
    list_names = list(pepdata_df['Joined'].values)
    list_peps = list(pepdata_df['Peptide'].values)

    with open(fasta_file_output, 'w') as fasta:
        for i in range(len(list_names)):
            fasta.write(list_names[i] + '\n')
            fasta.write(list_peps[i] + '\n')


def _return_ids_peps_df(pepdata):
    df = pepdata.drop_duplicates(subset='Peptide')
    df = df.loc[df['Peptide'].str.contains('-') == False]
    id_initial_nmer = df['Protein ID'].str.cat([df['Initial AA'].apply(str),
                                                df['n-mer'].apply(str)], sep='_')
    df['Joined'] = id_initial_nmer
    return df[['Joined', 'Peptide']]


def _retrieve_relevant_positions(df):
    return df[['Pos', 'n-mer', 'Allele']].values.tolist()


def _deep_list_to_df(list_):
    headers = ['Peptide', 'n-mer', 'Protein ID', 'Initial AA', 'Allele']
    df = pandas.DataFrame(list_, columns=headers)
    df['Peptide'] = df['Peptide'].str.replace('X', '-')
    df = df.loc[df['Peptide'].str.contains('--') == False]
    return df


def fasta_per_allele_per_nmer(pepdata, filepath):
    nmers = _get_nmers(pepdata)
    alleles = _get_alleles(pepdata)
    prot_ids = _get_prot_ids(pepdata)
    for prot_id in prot_ids:
        for nmer in nmers:
            for allele in alleles:
                sliced = _slice_df(pepdata, nmer, allele)
                if check_size(sliced):
                    new_name = sliced.apply(_apply_transf_to_protein_id, axis=1)
                    sliced['transf'] = pandas.Series(new_name, index=sliced.index)
                    peptide_data = _return_list_of_lists(sliced)
                    _write_fasta(filepath, nmer, allele, peptide_data, prot_id)


def _apply_transf_to_protein_id(x):
    name = x['Protein ID'].split('_')[1]
    new_name = '>' + str(x['Initial AA']) + '_' + name
    return new_name


def _slice_df(pepdata, nmer, allele, prot_id):
    return pepdata.loc[(pepdata['n-mer'] == nmer) & (pepdata['Allele'] == allele) & (pepdata['ID'] == prot_id)]


def check_size(sliced):
    if len(sliced) == 0:
        return False
    else:
        return True


def _return_list_of_lists(sliced):
    return sliced[['Peptide', 'transf']].values.tolist()


def _write_fasta(filepath, nmer, allele, peptide_data, prot_id):
    allele = allele.replace(':', '-')
    with open(filepath + 'fasta_%s_%s_%s.fasta' % (nmer, allele.replace(':', '_'), prot_id), 'w') as out:
        for i in peptide_data:
            out.write(i[1] + '\n')
            out.write(i[0].replace('-', 'X') + '\n')


#####################################################
# UTILITIES TO CREATE PEPTIDE SWAPS IN FASTA FORMAT #
#####################################################

def find_swaps_write_to_fasta(high_priority_df, fasta_files_dir):
    nmers = _get_nmers(high_priority_df)
    alleles = _get_alleles(high_priority_df)
    prot_ids = _get_prot_ids(high_priority_df)

    for prot_id in prot_ids:
        os.mkdir(fasta_files_dir + '%s/' % prot_id)
        for nmer in nmers:
            for allele in alleles:
                sliced = _slice_df(high_priority_df, nmer, allele, prot_id)
                if check_size(sliced):
                    list_of_lists = sliced[['n-mer', 'Allele', 'ID', 'Pos', 'Peptide']].values.tolist()
                    for item in list_of_lists:
                        swaps = _create_swaps(item[-1])
                        _open_write_fasta(fasta_files_dir, item, swaps, prot_id)


def _open_write_fasta(fasta_files_dir, data, swaps, prot_id):
    file_name = "_".join([fasta_files_dir + '%s/' % prot_id, 'swap', data[-1], 'Pos', str(data[-2]), 'ID', str(data[-3]).replace('_', '-'),
                          'Allele', str(data[-4]), 'nmer', str(data[-5])])

    with open(file_name + '.fasta', 'w') as inf:
        for swap in swaps:
            inf.write("".join(['>', prot_id, '_', swap, '\n']))
            inf.write(swap + '\n')


def _create_swaps(peptide):
    list_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    list_peps = []
    for i in range(len(peptide)):
        for k in range(len(list_AA)):
            list_peps.append(_insert_aa(peptide, i, list_AA[k]))

    return list_peps


def _insert_aa(string, index, aa):
    hash_string = list(string)
    del hash_string[index]
    hash_string.insert(index, aa)
    return "".join(hash_string)

def _get_prot_ids(pepdata):
    return pepdata['ID'].unique()

def _get_nmers(pepdata):
    return pepdata['n-mer'].unique()


def _get_alleles(pepdata):
    return pepdata['Allele'].unique()


def check_size(sliced):
    if len(sliced) == 0:
        return False
    else:
        return True

def create_separate_lists(fasta_file):
    """
    Creates 2 lists from a fasta file
    :param fasta_file: file
    :return: one list for the IDs in the file and one list for the proteins/peptides in it
    """
    with open(fasta_file) as infile:
        all_list = []
        peptide = ""
        lines = infile.readlines()
        for i in range(0, len(lines)):
            if lines[i].startswith('>'):
                all_list.append(lines[i].rstrip())
            else:
                peptide += lines[i].rstrip()
            try:
                if lines[i + 1].startswith('>'):
                    all_list.append(peptide)
                    peptide = ""
            except:
                all_list.append(peptide)
        j = []
        k = []
        for i in all_list:
            if i.startswith('>'):
                j.append(i.strip('>'))
            else:
                k.append(i)
    return j, k
