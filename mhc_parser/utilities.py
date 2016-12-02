

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






