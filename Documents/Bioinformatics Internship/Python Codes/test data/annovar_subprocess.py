import os
import subprocess

os.chdir("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/variantAnnotation/")
with open("PATHS.txt") as f:
    content = f.readlines()
    print content


ANNOVAR_PATH = content[0].split('=')[1].strip('\n').strip(' ')
INPUT_VCF_PATH = content[1].split('=')[1].strip('\n').strip(' ')
OUTPUT_CSV_PATH = content[2].split('=')[1].strip('\n').strip(' ')






def run_annovar(annovar_path, input_path, output_path):
