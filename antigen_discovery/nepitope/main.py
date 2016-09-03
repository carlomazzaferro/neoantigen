import sys
sys.path.append("/Users/carlomazzaferro/Documents/Code/neantigen/antigen_discovery")

from nepitope import alignment_tool

filepath = '/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/test.fasta'
ID = 'Q99ZW2'
website = 'uniprot'

alignment_tool.fetch_protein(filepath, website=website, ID=ID, fetch_fasta=True)

#run blast
#blastp -out /Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/BLAST_out.xml -outfmt 5 -query /Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/test.txt -db nr -remote -evalue 0.001


xml_filepath = '/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/BLAST_out.xml'

alignment_tool.viz_blast(xml_filepath)

treshold = 70
fasta_out = '/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/fasta_out_from_blast.fasta'

alignment_tool.write_blast_filtered_output(xml_filepath, fasta_out, treshold)

#Run clustal omega
#!clustalo -i /Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/fasta_out_from_blast.fasta --residuenumber -o /Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/msa_aligned.clustal --outfmt=clustal
#viz tool:
#!mview -in clustal -ruler on -html head -coloring any /Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/msa_aligned.clustal > /Users/carlomazzaferro/Documents/Code/clustal.html

import urllib
f = urllib.urlopen("/Users/carlomazzaferro/clustal.html").read()
print f

import os
import webbrowser
try:
    from urllib import pathname2url         # Python 2.x
except:
    from urllib.request import pathname2url # Python 3.x

url = 'file:{}'.format(pathname2url(os.path.abspath("/Users/carlomazzaferro/clustal.html")))
webbrowser.open(url)


msa_alig = '/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/msa_aligned.fasta'
alignment_tool.read_clustal_alignment(msa_alig)

os.chdir(os.path.dirname(os.path.realpath('__file__'))


from nepitope import scoring_utils


msa_alig = '/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/msa_aligned.fasta'
cls = scoring_utils.Score(msa_alig, [9,10,11])
import os

cls.create_score_file('/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/nmers_out/nmerized_9_3', '/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/nmers_out/')



ls = cls.run_scoring('/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/nmers_out/')

