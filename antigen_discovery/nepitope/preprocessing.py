import pandas
from nepitope import scoring_utils


def join_files(files):
    list_dfs = []

    for i in files:
        df = load_to_df(i)
        if len(df) != 89:
            print('Wrong Length')
        list_dfs.append(df)
    return list_dfs


def load_to_df(file_):

    with open(file_) as inf:
        names = []
        nums = []
        proteins = []
        name_numbered = []

        for line in inf:
            line = line.split()
            proteins.append(line[-1])
            num = line[0]
            nums.append(num)
            name = line[1:-1]
            name = "".join(name)
            name_numbered.append(num + '_' + name)
            names.append(name)

    df = pandas.DataFrame([name_numbered, nums, proteins], index=['Name', 'Num', 'Sequence']).T

    return df


def test_consistency(list_dfs):
    list_names = []
    for i in list_dfs:
        if len(i.Name.unique()) != 89:
            raise ValueError('Wrong length, can\'t concatenate')

    for i in list_dfs:
        list_names.append(list(i.Name.unique()))

    for i in range(0, len(list_names)):
        for j in range(0, len(list_names)):
            if list_names[i] != list_names[j]:
                raise ValueError('Inconsistent names')


def test_series_consistency(df):
    seq = df.Sequence
    for i in range(0, len(seq)):
        for j in range(0, len(seq)):
            if len(seq[i]) != len(seq[j]):
                raise ValueError('Inconsistent lengths')


def concat_dfs(list_dfs):
    base_df = list_dfs[0]
    del list_dfs[0]

    for i in range(0, len(list_dfs)):
        base_df.Sequence = base_df.Sequence.str.cat(list_dfs[i].Sequence)

    return base_df


def load_motifs_and_conservation_locs(file_):

    with open(file_) as inf:
        data_motifs = []
        data_inf_pos = []

        for line in inf:
            if line.startswith('Motifs'):
                line = line.split()
                data_motifs.append(line[-1])
            else:
                line = line.split()
                data_inf_pos.append(line[-1])

    joined_motifs = "".join(data_motifs)
    joined_inf_pos = "".join(data_inf_pos)

    return joined_motifs, joined_inf_pos


def dealign_and_add_reference_protein(original_fasta, new_fasta, ref_prot, ref_prot_name):

    names, proteins = scoring_utils.Score.create_separate_lists(original_fasta)

    with open(new_fasta, 'w') as fasta:
        for i in range(0, len(names)):
            fasta.write('>' + names[i] + '\n')
            fasta.write(proteins[i].replace('-', "") + '\n')

        fasta.write(ref_prot_name + '\n')
        fasta.write(ref_prot)


