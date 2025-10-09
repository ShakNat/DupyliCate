#code used to obtain pep file without alternate transcripts

import re, os, sys, subprocess, gzip, glob
import pandas as pd
import csv
import shutil
from pathlib import Path

#Removing alternate transcripts from the peptide FASTA file
def pepclean(gff3_input_file, output_file, no_trans_pep):
	if gff3_input_file[-2:].lower() != 'gz':#uncompressed gff file
    	transcripts_per_gene = {}
    	with open(gff3_input_file, "r") as f:
        	line = f.readline()
        	while line:
            	if line[0] != "#":
                	parts = line.strip().split('\t')
                	if len(parts) > 2:
                    	if (parts[2] == "transcript") or (parts[2] == "mRNA"):
                        	partsnew = parts[-1].strip().split(';')
                        	for each in partsnew:
                            	pattern_par = r'^Parent=.*$'
                            	if re.match(pattern_par, each):
                                	partsnew1 = str(each).replace("Parent=", "")
                        	for every in partsnew:
                            	pattern_ID = r'^ID=.*$'
                            	if re.match(pattern_ID, every):
                                	partsnew0 = str(every).replace("ID=", "")
                        	try:
                            	transcripts_per_gene[partsnew1].append(partsnew0)
                        	except KeyError:
                            	transcripts_per_gene.update({partsnew1: [partsnew0]})
            	line = f.readline()
	else:#Dealing with compressed GFF3 file
    	transcripts_per_gene = {}
    	with gzip.open(gff3_input_file, "rt") as f:
        	line = f.readline()
        	while line:
            	if line[0] != "#":
                	parts = line.strip().split('\t')
                	if len(parts) > 2:
                    	if (parts[2] == "transcript") or (parts[2] == "mRNA"):
                        	partsnew = parts[-1].strip().split(';')
                        	for each in partsnew:
                            	pattern_par=r'^Parent=.*$'
                            	if re.match(pattern_par, each):
                                	partsnew1 = str(each).replace("Parent=", "")
                        	for every in partsnew:
                            	pattern_ID = r'^ID=.*$'
                            	if re.match(pattern_ID, every):
                                	partsnew0 = str(every).replace("ID=", "")
                        	try:
                            	transcripts_per_gene[partsnew1].append(partsnew0)
                        	except KeyError:
                            	transcripts_per_gene.update({partsnew1: [partsnew0]})
            	line = f.readline()
	gene_names = list(transcripts_per_gene.keys())
#dealing with uncompressed pep file
	if output_file[-2:].lower() != 'gz':
    	with open(output_file, "r") as f:
        	with open(no_trans_pep, "w") as out:
            	line = f.readlines()
            	pep_dict = {}
            	for each in line:
                	if '>' in each:
                    	ind = line.index(each)
                    	pep_dict.update({str(each).replace('>', '').replace('\n', ''): line[ind + 1]})
            	transcripts = list(pep_dict.keys())
            	for gene in gene_names:
                	trans_length = []
                	for trans in transcripts_per_gene[gene]:
                    	for each in transcripts:
                        	if trans == each:
                            	if len(transcripts_per_gene[gene]) < 2:
                                	out.write('>' + str(each) + '\n' + str(pep_dict[each]))
                            	else:
                                	trans_length.append((each, pep_dict[each]))
                	if len(trans_length) == 0:
                    	pass
                	else:
                    	seq = max(trans_length, key=len)
                    	best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
                    	out.write('>' + str(best_trans) + "\n" + str(seq))


	else:#dealing with compressed pep file
    	with gzip.open(output_file, "rt") as f:
        	with open(no_trans_pep, "w") as out:
            	line = f.readlines()
            	pep_dict = {}
            	for each in line:
                	if '>' in each:
                    	ind = line.index(each)
                    	pep_dict.update({str(each).replace('>', '').replace('\n', ''): line[ind + 1]})
            	transcripts = list(pep_dict.keys())
            	for gene in gene_names:
                	trans_length = []
                	for trans in transcripts_per_gene[gene]:
                    	for each in transcripts:
                        	if trans == each:
                            	if len(transcripts_per_gene[gene]) < 2:
                                	out.write('>' + str(gene) + '\n' + str(pep_dict[each]))
                            	else:
                                	trans_length.append(pep_dict[each])
                	if len(trans_length) == 0:
                    	pass
                	else:
                    	seq = max(trans_length, key=len)
                    	out.write('>' + str(gene) + "\n" + str(seq))
pepclean(GFF_file, PEP_file,path_to_write_pep_file_without_alternate_transcripts )

#code to process GFF files for DupGen_finder

def extract_headers_from_fasta(fasta_file):
	headers = set()
	with open(fasta_file) as f:
    	for line in f:
        	if line.startswith(">"):
            	header = line[1:].strip().split()[0]  # remove ">" and any description
            	headers.add(header)
	return headers


def filter_gff_by_headers(gff_file, headers, output_file):
	with open(gff_file) as gff, open(output_file, "w") as out:
    	for line in gff:
        	if line.startswith("#") or "\tmRNA\t" not in line or "\ttranscript\t" not in line:
            	continue
        	fields = line.strip().split("\t")
        	attributes = fields[8]
        	attr_dict = dict(attr.split("=", 1) for attr in attributes.split(";") if "=" in attr)
        	gene_id = attr_dict.get("ID", "")
        	if gene_id in headers:
            	contig = fields[0]
            	start = fields[3]
            	end = fields[4]
            	out.write(f"{contig}\t{gene_id}\t{start}\t{end}\n")

headers = extract_headers_from_fasta(fasta_file)
filter_gff_by_headers(gff_file, headers, output_file)
