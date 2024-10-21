import os
import subprocess
from Bio import SeqIO
from src.CLEAN.utils import *
from src.CLEAN.infer import *

# these are fastas
#input_list = ["data/med4_uniprotkb_proteome_UP000001026_2024_01_13"] 
# input_list = ["data/hm2_uniprotkb_proteome_UP000006538_2024_01_13", "data/ssp7_uniprotkb_proteome_UP000258925_2024_01_13"]
# input_list = ["data/GCA_000011465.1_ASM1146v1_genomic_MED4_prot", 
# 	"data/GCA_000858745.1_ViralProj15134_genomic_SSP_7_prot",
# 	"data/GCA_000892475.1_ViralProj64705_genomic_HM2_prot"
# ]
import sys

# print(sys.argv)

input_list = sys.argv[1:]

# print(input_list)

for i in input_list:
	subprocess.run(['echo', 'ensure the fasta only has IDs in the headers: '])
	fasta_name = i+".fasta"
	output_fasta_file = os.path.splitext(fasta_name)[0] + '_only_ids_in_headers.fasta'
	with open(fasta_name) as input_handle, open(output_fasta_file, "w") as output_handle:
		for record in SeqIO.parse(input_handle, "fasta"):
			# write only the ID as the header to the output fasta
			output_handle.write(f">{record.id}\n{record.seq}\n")
	subprocess.run(['echo', 'about to run retrive_esm1b_embedding() for: ', os.path.basename(output_fasta_file)])
	retrive_esm1b_embedding(os.path.splitext(os.path.basename(output_fasta_file))[0])
	# then create csv since clean's version doesn't work
	# this code keeps only the id from fasta header in the Entry col, not the additional annotation info
	subprocess.run(['echo', 'about to create the .csv from this .fasta: ', i])
	with open('./data/' + os.path.basename(i) + '.csv', 'w', newline='') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t')
		csvwriter.writerow(['Entry', 'EC number', 'Sequence'])
		for record in SeqIO.parse('./data/' + os.path.basename(i) + '.fasta', "fasta"):
			csvwriter.writerow([record.id, "", str(record.seq)])
	subprocess.run(['echo', 'about to run infer_maxsep() for: ', i])
	infer_maxsep('split100', os.path.basename(i), report_metrics=False, pretrained=True, gmm = './data/pretrained/gmm_ensumble.pkl')
	subprocess.run(['echo', 'about to run infer_pvalue() for: ', i])
	infer_pvalue('split100', os.path.basename(i), p_value = 1e-5, nk_random = 20, report_metrics = False, pretrained=True)
