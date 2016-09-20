import pandas
from nepitope import scoring_utils


def retrieve_orthologue_peptides(fasta_file_msa, high_aa_low_cons_df):
    list_of_pos_and_nmers = _retrieve_relevant_positions(high_aa_low_cons_df)
    pepdata = iterate_over_msa_file(list_of_pos_and_nmers, fasta_file_msa)
    return _deep_list_to_df(pepdata)


def iterate_over_msa_file(list_of_pos_and_nmers, fasta_file_msa):
    peps_and_positions = []
    list_ids, list_proteins = scoring_utils.Score.create_separate_lists(fasta_file_msa)

    for i in range(len(list_proteins) - 1):
        for position in list_of_pos_and_nmers:
            nmer = position[1]
            initial_amino_acid = position[0]
            final_amino_acid = position[0] + nmer
            peps_and_positions.append([list_proteins[i][initial_amino_acid:final_amino_acid],
                                       nmer, list_ids[i], initial_amino_acid])

    return peps_and_positions


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
    return df[['Pos', 'n-mer']].values.tolist()


def _deep_list_to_df(list_):
    headers = ['Peptide', 'n-mer', 'Protein ID', 'Initial AA']
    df = pandas.DataFrame(list_, columns=headers)
    df['Peptide'] = df['Peptide'].str.replace('X', '-')
    df = df.loc[df['Peptide'].str.contains('--') == False]
    return df
