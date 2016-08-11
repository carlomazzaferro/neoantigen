
import os
import numpy as np
import pandas as pd
pd.set_option('display.width', 100)
pd.set_option('max_colwidth', 80)
import matplotlib as mpl
import seaborn as sns
sns.set_context("notebook", font_scale=1.4)
from IPython.display import display, HTML
from epitopepredict import base, sequtils, analysis


fastafile = '/Volumes/Seagate Backup Plus Drive/vcf_files/varcode_to_test/mhcpanii_outs/Peptides_RNA.txt'

#get data in fasta format
zaire = sequtils.fasta2Dataframe(fastafile)


P = base.getPredictor('netmhciipan')
savepath1 = 'netmhciipan'
#run prediction for several alleles and save results to savepath
alleles = ["HLA-DRB1*0101", "HLA-DRB1*0108", "HLA-DRB1*0305", "HLA-DRB1*0401",
           "HLA-DRB1*0404", "HLA-DRB3*0101", "HLA-DRB4*0104"]

df = P.predictProteins(zaire, length=11, alleles=alleles, save=True, path=savepath1)