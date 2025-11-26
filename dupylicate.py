### v1.2 ###
### Shakunthala Natarajan ###
### Boas Pucker ###

__usage__ = """
			python3 dupylicate.py

			MANDATORY:
			--gff <full path to folder containing GFF files>
			--fasta <full path to folder containing WGS FASTA files> | --cds <full path to folder containing CDS FASTA files> | --pep <full path to folder containing PEP files>
			--out <full path to output folder>

			OPTIONAL:
			--ref <name of the organism to be taken as reference> Default is NA and the script runs in reference-free mode
			--prokaryote <use the flag if the organisms for analysis are prokaryotes>
			--pseudos_gff <yes | no for optional inclusion or exclusion of pseudogenes with coding features in the GFF file> Default is no 
			--gff_config <full path to txt file containing the gff parameters to be used separated by spaces or tabs;
						  Needs 4 columns in the order 
							(i) base file name - same as the base name you use for the gff file | all in case all the gff 
								have the same attribute pattern
							(ii) child_attribute: attribute field of the mRNA or transcript feature in the file like ID
							(iii) child_parent_linker: attribute field of the mRNA or transcript, CDS, exon features that link them with their 
								  respective parent feature like Parent - Note: base assumption by the tool is that all child levels
								  have the same child-parent linker attribute fields. For eg., if Parent is the child-parent linker in the mRNA feature line,
								  then Parent will be the child-parent linker for all other child-level feature lines in the GFF
							(iv) parent_attribute: attribute field of the gene feature like ID 
			--mode <overlap | strict modes; in the overlap mode genes repeat among the different duplicates classification; 
				   in the strict mode there is no gene repetition and you have a new classification group called mixed 
				   duplicates containing the related connected components from the other three duplicate classes> DEFAULT is overlap
			--to_annotate <full path to file containing names of queries and the corresponding reference organism for GeMoMa annotation separated by comma - Query,Reference -> one pair per line>
			--seq_aligner <choose one among blast | diamond | mmseqs2 > DEFAULT is DIAMOND
			--blast <full path to BLAST if not already in your PATH environment variable>
			--diamond <full path to diamond if not already in your PATH environment variable>
			--mmseqs <full path to mmseqs2 if not already in your PATH environment variable>
			--mafft <full path to MAFFT if not already in your PATH environment variable>
			--tree <yes | no to include or not include phylogeny in ortholog detection> < default is yes>
			--evalue <evalue for alignment> DEFAULT is 1e-5
			--gemoma <full path to GeMoMa if not already in your PATH environment variable>
			--qc < yes | no for quality check with BUSCO> Default is no
			--busco  <full path to BUSCO or choose 'busco_docker' for docker-based BUSCO installation> Default is 'busco'
			--busco_version <BUSCO version in the format vx.x.x - needed only if you have docker-based BUSCO installation> Default is v6.0.0
			--busco_lineage <full path to tab separated config file where the first column is sample organism name and the second column is busco lineage>
			--container_version <docker container version of BUSCO> Default is cv1
			--docker_host_path < full host folder path - needed for docker-based BUSCO installation>
			--docker_container_path < full mount path in the docker container - needed for docker-based BUSCO installation>
			--score <auto for automatic threshold finding and <float number between 0 and 1> for manual threshold finding for segregating singletons and duplicates> DEFAULT is auto
			--self_simcut <float similarity percentage to remove self alignment hits with low similarity percentage> Default is 50.0
			--hits <number of top hits to be considered for finding suitable orthologs in the reference organism> Default is 10
			--ortho_candidates <user specified integer value for listing the potential ortholog candidates for low confidence ortholog cases> Default is 3
			--occupancy <float number cutoff for MAFFT aligned file trimming> Default is 0.1
			--scoreratio < float ratio of forward alignment bit score and self alignment bit score of query to assess if the forward hit is valid> Default is 0.3
			--fwd_simcut < float similarity percentage to remove forward alignment hits with low similarity percentage> Default is 40.0
			--cores <number of cores needed to run Dupylicate analysis> DEFAULT is 4
			--proximity <integer value for the number of intervening genes to detect proximal duplications> DEFAULT is 10
			--synteny_score <float value which is used as a cut-off or threshold for synteny analysis> DEFAULT is 0.5
			--flank <integer value specifying the number of flanking genes to be considered to determine the synteny window size in synteny analysis>
			--side <integer value for synteny support from either side of a flanking region of a synteny window>
			--ka_ks <yes | no to calculate ka, ks values> DEFAULT is 'no'
			--ka_ks_method < MYN | NG methods for Ka/Ks ratio calculation> Default is NG
			--duplicates_analysis <yes | no for further statistical analysis of identified gene duplicates> DEFAULT is no
			--specific_duplicates_analysis <yes | no for further statistical analysis of specified ref genes' gene duplicates> DEFAULT is no
			--ref_free_specific_duplicates_analysis <full path to txt file containing the names of genes whose expression analysis is to be done in the absence of a reference>
			--dpi <resolution value desired for plots low | moderate | high | very high > DEFAULT is moderate
			--analyse_disperse_duplicates <yes | no for statistical analysis of dispersed gene duplicates> DEFAULT is no
			--exp <full path to the folder with expression counts files> 
			--avg_exp_cut < float value cutoff for average expression value across samples for classifying a gene as pseudogene based on gene expression> Default is 1
			--genes <integer value referring to the number of genes to be taken for random pairing to determine threshold 
					value for assessing functional divergence of gene duplicates based on gene expression values> DEFAULT is 10000
			--specific_genes <user given list of specific genes from the reference with one transcript name per line, for analysing gene duplication outputs>
			--clean_up <yes | no cleans up intermediate files/ folders> DEFAULT is YES

			ALLOWED FILE EXTENSIONS:

			<'.gff', '.gff.gz', '.gff3', '.gff3.gz', '.fa', '.fa.gz', '.fna', '.fna.gz','.fasta', '.fasta.gz',
			'.genome.fasta', '.genome.fa', '.genome.fasta.gz', '.genome.fa.gz', '.genome.fna',
			'.genome.fna.gz','.cds.fa', '.cds.fasta','.cds.fa.gz', '.cds.fasta.gz', '.cds.fna',
			'.cds.fna.gz','.pep.fa','.pep.fa.gz','.pep.fasta','.pep.fasta.gz', '.pep.fna',
			'.pep.fna.gz','.tsv','.txt','.txt.gz','.tpms.txt','.tpms.txt.gz'>

			"""

######### imports #########
import re, os, sys, subprocess, gzip, glob
import pickle
import hashlib
import json
import io
import pandas as pd
import csv
import shutil
from pathlib import Path
import fileinput
import numpy as np
import fnmatch
import itertools
from numpy import mean
from scipy import stats
from operator import itemgetter
from datetime import datetime
from collections import Counter
from collections import defaultdict
import logging
try:
	import dendropy
	from dendropy import Tree, PhylogeneticDistanceMatrix, TaxonNamespace
except ImportError:
	pass
# Set DendroPy's logger to only show errors
logging.getLogger('dendropy').setLevel(logging.ERROR)
import traceback
import copy
import seaborn as sns
import math
from math import ceil
import tempfile
import random
import time
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress, gaussian_kde, skew, kurtosis, fisher_exact
from scipy.signal import savgol_filter
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import concurrent.futures
import multiprocessing
from functools import reduce
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
try:
	from tqdm import tqdm
	tqdm_available = True
except ImportError:
	tqdm_available = False
##### end of imports #####
#global definition of dictionary with problematic characters for dendropy processing, busco processing and their respective placeholders
replacements = {
		"'": "§quo",
		":": "§col",
		",": "§com",
		"(": "§lbr",
		")": "§rbr",
		";": "§sco",
		"[": "§lsbr",
		"]": "§rsbr",
		"\\": "§bsl",
		"\"": "§dquo",
		"/": "§fsl"
		}
# Set numpy to use float64 everywhere
np.set_printoptions(precision=15)
#Defining a global dictionary of defaults for json file parameter recording for Dupylicate run
Defaults = {
	'ref': 'NA',
	'org_type': 'eukaryote',
	'pseudos_gff': 'no',
	'gff_config': 'default',
	'out': 'Output_dir',
	'gff': 'GFF file(s)',
	'fasta': 'NA',
	'cds': 'NA',
	'pep': 'NA',
	'qc': 'no',
	'busco': 'NA',
	'busco_version': 'NA',
	'busco_lineage':'auto',
	'container_version':'NA',
	'docker_host_path': 'NA',
	'docker_container_path': 'NA',
	'seq_aligner': 'NA',
	'blast': 'NA',
	'diamond': 'NA',
	'mmseqs': 'NA',
	'tree': 'yes',
	'mafft': 'NA',
	'occupancy': 0.1,
	'fasttree': 'NA',
	'to_annotate': 'NA',
	'gemoma': 'NA',
	'cores': 4,
	'evalue': '1e-5',
	'hits': 10,
	'ortho_candidates': 3,
	'score': 'auto',
	'scoreratio': 0.3,
	'fwd_simcut': 40.0,
	'self_simcut': 50.0,
	'proximity': 10,
	'synteny_score': 0.5,
	'flank': 5,
	'side': 1,
	'ka_ks': 'no',
	'ka_ks_method': 'NA',
	'duplicates_analysis': 'no',
	'specific_duplicates_analysis': 'no',
	'ref_free_specific_duplicates_analysis':'NA',
	'dpi': 'moderate',
	'analyse_disperse_duplicates': 'no',
	'exp': 'NA',
	'avg_exp_cut': 1,
	'genes': 10000,
	'specific_genes': 'NA',
	'mode': 'overlap',
	'clean_up': 'yes',
	'dpi': 'moderate'
}

#function to record the Dupylicate run parameters in a json file
def parse_arguments(arguments):
	parameters = Defaults.copy()
	if '--prokaryote' in arguments:
		parameters['org_type'] = 'prokaryote'
	else:
		parameters['org_type'] = 'eukaryote'
	for i, arg in enumerate(arguments):
		if arg.startswith("--"):
			key = arg[2:]
			if key in Defaults:
				if i + 1 >= len(arguments):
					raise ValueError(f"Missing value for {arg}")
				value = arguments[i + 1]
				default_type = type(Defaults[key])

				# Handle yes/no conversion
				if Defaults[key] in ("yes", "no"):
					val = value.lower()
					if val in {"true", "1", "yes", "y"}:
						parameters[key] = "yes"
					elif val in {"false", "0", "no", "n"}:
						parameters[key] = "no"
					else:
						raise ValueError(f"Invalid yes/no value for --{key}: {value}")
				else:
					parameters[key] = default_type(value)
	return parameters

#function to calculate md5sums of input files
def calculate_md5sum (file_path, block_size = 65536):
	hasher = hashlib.md5()
	if file_path.endswith('.gz'):
		with gzip.open(file_path, 'rb') as f:
			while chunk := f.read(block_size):
				hasher.update(chunk)
	else:
		with open(file_path, 'rb') as f:
			while chunk := f.read(block_size):
				hasher.update(chunk)
	return hasher.hexdigest()

#function to distribute cores for a single level parallelism approach
def distribute_cores_single_level(total_cores, num_tasks):
	# If fewer tasks than cores, give each task multiple cores
	if num_tasks <= total_cores:
		tasks_in_parallel = num_tasks
		cores_per_task = total_cores // num_tasks
		return tasks_in_parallel

	# If more tasks than cores, run as many as possible with 1 core each
	else:
		tasks_in_parallel = total_cores
		cores_per_task = 1
		return tasks_in_parallel

#function to distribute cores for a nested parallelism approach
def distribute_cores(total_cores, num_organisms):
	# Determine optimal distribution based on workload and cores available

	# If only one organism, use all cores for inner processes
	if num_organisms == 1:
		return 1, total_cores

	# Calculate square root as a starting point for balancing
	# This creates a reasonable balance between parallelism levels
	sqrt_cores = int(total_cores ** 0.5)

	# Initial allocation based on square root
	outer_cores = min(num_organisms, max(1, sqrt_cores))

	# Allocate remaining cores to inner processes
	inner_cores = max(1, total_cores // outer_cores)

	# Adjust if cores are wasted
	while (outer_cores * inner_cores) < total_cores and outer_cores < num_organisms:
		outer_cores += 1
		inner_cores = max(1, total_cores // outer_cores)

	# Step back in case balance is off
	if (outer_cores * inner_cores) > total_cores:
		outer_cores -= 1
		inner_cores = max(1, total_cores // outer_cores)

	return outer_cores, inner_cores

#Getting base names of organisms from the folder/ file paths without extensions
def get_basename(folder_path):
	extensions = ['.gff', '.gff.gz', '.gff3', '.gff3.gz','.fa', '.fa.gz', '.fna', '.fna.gz','.fasta', '.fasta.gz',
				  '.cds.fa', '.cds.fasta','.cds.fa.gz', '.cds.fasta.gz', '.cds.fna','.cds.fna.gz',
				  '.pep.fa','.pep.fa.gz','.pep.fasta','.pep.fasta.gz', '.pep.fna','.pep.fna.gz',
				  '.tsv','.txt','.txt.gz','.genome.fasta','.genome.fa','.genome.fasta.gz',
				  '.genome.fa.gz','.genome.fna','.genome.fna.gz','.tpms.txt','.tpms.txt.gz','.pkl','.json']
	extensions_sorted = sorted(extensions, key = len, reverse = True)
	if os.path.isfile(folder_path):
		files = [folder_path]
	elif os.path.isdir(folder_path):
		files = []
		for ele in extensions:
			files.extend(glob.glob(os.path.join(folder_path, ele)))
	input_file_names_wo_ext = []
	for every in files:
		base_name = os.path.basename(str(every))  # getting the file basename
		matched_extension = None
		for ext in extensions_sorted:
			if base_name.endswith(ext):
				matched_extension = ext
				break
		if matched_extension:
			base_name = base_name[:-len(matched_extension)]
			input_file_names_wo_ext.append(base_name)
	return input_file_names_wo_ext

#Function to clean the input FASTA file
def clean_fasta(fasta_file_path, output_fasta, fasta_type,logger):
	orgname = get_basename(fasta_file_path)[0]
	def process_header(line):
		"""Clean header line by removing spaces and tabs from identifier"""
		if " " in line:
			tmp = line.split(' ')[0]
			if '\t' in tmp:
				tmp = tmp.split('\t')[0]
				return tmp + "\n"
			else:
				return tmp + "\n"
		elif "\t" in line:
			tmp = line.split('\t')[0]
			return tmp + "\n"
		else:
			return line.strip() + "\n"

	if fasta_file_path[-2:].lower() != 'gz':  # dealing with uncompressed fasta file
		with open(output_fasta, "w") as out:
			with open(fasta_file_path, "r") as f:
				current_sequence = ""

				for line in f:
					line = line.strip()
					if not line:  # Skip empty lines
						continue

					if line[0] == ">":  # Header line
						# Write previous sequence if exists
						if current_sequence:
							if fasta_type=='pep':
								if 'U' in current_sequence or 'u' in current_sequence:#replace selenocysteine with X
									current_sequence=current_sequence.replace('U','X')
									current_sequence=current_sequence.replace('u', 'x')
									logger.info(f"Replaced selenocysteine with X in a sequence from {orgname}")
								if 'O' in current_sequence or 'o' in current_sequence:#replace pyrroline with X
									current_sequence=current_sequence.replace('O','X')
									current_sequence=current_sequence.replace('o', 'x')
									logger.info(f"Replaced pyrroline with X in a sequence from {orgname}")
							out.write(current_sequence + "\n")
							current_sequence = ""
						# Process and write header
						cleaned_header = process_header(line)
						out.write(cleaned_header)
					else:  # Sequence line
						current_sequence += line

				# Write the last sequence
				if current_sequence:
					if fasta_type == 'pep':
						if 'U' in current_sequence or 'u' in current_sequence:#replace selenocysteine with X
							current_sequence=current_sequence.replace('U', 'X')
							current_sequence=current_sequence.replace('u', 'x')
							logger.info(f"Replaced selenocysteine with X in a sequence from {orgname}")
						if 'O' in current_sequence or 'o' in current_sequence:#replace pyrroline with X
							current_sequence=current_sequence.replace('O', 'X')
							current_sequence=current_sequence.replace('o', 'x')
							logger.info(f"Replaced pyrroline with X in a sequence from {orgname}")
					out.write(current_sequence + "\n")

	else:  # dealing with compressed fasta file
		with gzip.open(output_fasta, "wt") as out:
			with gzip.open(fasta_file_path, "rt") as f:
				current_sequence = ""

				for line in f:
					line = line.strip()
					if not line:  # Skip empty lines
						continue

					if line[0] == ">":  # Header line
						# Write previous sequence if exists
						if current_sequence:
							if fasta_type=='pep':
								if 'U' in current_sequence or 'u' in current_sequence:#replace selenocysteine with X
									current_sequence=current_sequence.replace('U','X')
									current_sequence=current_sequence.replace('u', 'x')
									logger.info(f"Replaced selenocysteine with X in a sequence from {orgname}")
								if 'O' in current_sequence or 'o' in current_sequence:#replace pyrroline with X
									current_sequence=current_sequence.replace('O','X')
									current_sequence=current_sequence.replace('o', 'x')
									logger.info(f"Replaced pyrroline with X in a sequence from {orgname}")
							out.write(current_sequence + "\n")
							current_sequence = ""
						# Process and write header
						cleaned_header = process_header(line)
						out.write(cleaned_header)
					else:  # Sequence line
						current_sequence += line

				# Write the last sequence
				if current_sequence:
					if fasta_type == 'pep':
						if 'U' in current_sequence or 'u' in current_sequence:#replace selenocysteine with X
							current_sequence=current_sequence.replace('U', 'X')
							current_sequence=current_sequence.replace('u', 'x')
							logger.info(f"Replaced selenocysteine with X in a sequence from {orgname}")
						if 'O' in current_sequence or 'o' in current_sequence:#replace pyrroline with X
							current_sequence=current_sequence.replace('O', 'X')
							current_sequence=current_sequence.replace('o', 'x')
							logger.info(f"Replaced pyrroline with X in a sequence from {orgname}")
					out.write(current_sequence + "\n")
	return output_fasta

#Function to parse GFF3 and FASTA assembly files
def parse_fasta(fasta_file):
	fasta_lengths = {}
	fasta_invalid_chars = {}
	if fasta_file[-2:].lower() != 'gz':  # dealing with uncompressed fasta file
		with open(fasta_file, "r") as f:
			header = None
			sequence = []
			for line_num, line in enumerate(f, 1):
				line = line.strip()
				if line.startswith(">"):
					if header:
						seq = "".join(sequence)
						fasta_lengths[header] = len(seq)
						invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
						if invalid_chars:
							fasta_invalid_chars[header] = invalid_chars
					header = line[1:]  # Exclude '>'
					sequence = []
				else:
					sequence.append(line)
			if header:
				seq = "".join(sequence)
				fasta_lengths[header] = len(seq)
				invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
				if invalid_chars:
					fasta_invalid_chars[header] = invalid_chars
	else: # compressed FASTA file
		with gzip.open(fasta_file,'rt')as f:
			header = None
			sequence = []
			for line_num, line in enumerate(f, 1):
				line = line.strip()
				if line.startswith(">"):
					if header:
						seq = "".join(sequence)
						fasta_lengths[header] = len(seq)
						invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
						if invalid_chars:
							fasta_invalid_chars[header] = invalid_chars
					header = line[1:]  # Exclude '>'
					sequence = []
				else:
					sequence.append(line)
			if header:
				seq = "".join(sequence)
				fasta_lengths[header] = len(seq)
				invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
				if invalid_chars:
					fasta_invalid_chars[header] = invalid_chars
	return fasta_lengths, fasta_invalid_chars

#validate GFF3 and cleaned FASTA files
def validate_gff3_and_fasta(gff3_file, fasta_file, error_file, warning_file,child_attribute,parent_attribute,gemoma,parent_attribute_correct):
	id_counter = 0
	parent_id_counter=0
	parent_counter=0
	seen_exact = set()
	seen_features = set()

	exact_duplicates = []
	feature_duplicates = []
	errors = []
	warnings = []
	message = []

	# Parse the FASTA file
	fasta_lengths, fasta_invalid_chars = parse_fasta(fasta_file)

	#making a set of the FASTA headers
	headers = set()
	if fasta_file[-2:].lower() != 'gz':
		with open(fasta_file,'r') as f:
			for line in f:
				if line.startswith('>'):
					headers.add(line[1:].strip())
	else:
		with gzip.open(fasta_file,'rt') as f:
			for line in f:
				if line.startswith('>'):
					headers.add(line[1:].strip())

	#making a set of the GFF contig names
	contigs = set()
	if gff3_file[-2:].lower() != 'gz':
		with open(gff3_file,'r') as f:
			for line in f:
				if not line.startswith('#') and line.strip():
					contigs.add(line.split('\t')[0].strip())
	else:
		with gzip.open(gff3_file,'rt') as f:
			for line in f:
				if not line.startswith('#') and line.strip():
					contigs.add(line.split('\t')[0].strip())

	if len(headers)>len(contigs):
		if contigs<headers:
			pass
		else:
			errors.append(f"FASTA file headers and the contig names in the GFF names do not match")

	if len(contigs) == len(headers):
		if headers == contigs:
			pass
		else:
			errors.append(f"FASTA file headers and the contig names in the GFF names do not match")

	if gff3_file[-2:].lower() != 'gz':
		with open(gff3_file,'r') as f:
			gff_lines = f.readlines()
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	else:
		with gzip.open(gff3_file,'rt') as f:
			gff_lines = f.readlines()
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_gene == False:
		errors.append(f"Level 1 feature 'gene' missing. Please process your GFF file with the AGAT wrapper script first.")
	if has_mrna == False and has_transcript == False:
		errors.append(f"Level 2 feature 'mRNA' or 'transcript' missing. Please process your GFF file with the AGAT wrapper script first.")
	if has_cds == False:
		errors.append(f"Level 3 feature 'CDS' missing. Please process your GFF file with the AGAT wrapper script first.")

	# Parse the GFF3 file
	if gff3_file[-2:].lower() != 'gz':  # dealing with uncompressed gff file
		with open(gff3_file, "r") as gff:
			for line_num, line in enumerate(gff, 1):
				id = None
				parent = None
				if line.startswith("#") or not line.strip():
					continue

				parts = line.strip().split("\t")
				attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
				id = attributes.get(child_attribute, "")
				parent_id=attributes.get(parent_attribute_correct,"")
				parent = attributes.get(parent_attribute, "")
				# Extract fields
				seqname, source, feature, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
				# check for missing attributes
				if has_mrna == True:
					if feature.upper() in ('MRNA') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('MRNA', 'CDS') and parent == "":
						#debug
						print(line)
						parent_counter += 1
				elif has_mrna == False and has_transcript == True:
					if feature.upper() in ('TRANSCRIPT') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('TRANSCRIPT', 'CDS') and parent == "":
						parent_counter += 1
				# Check exact duplicates
				full_line = line.strip()
				if full_line in seen_exact:
					exact_duplicates.append((line_num, full_line))
				else:
					seen_exact.add(full_line)
					# Check feature duplicates
					# Reset feature_key for each iteration
					feature_key = None
					if feature.upper() == 'GENE':
						feature_key = (seqname, feature, start, end, parent_id)
					if has_mrna==True:
						if feature.upper() == 'MRNA':
							feature_key = (seqname, feature, start, end, id, parent)
					elif has_transcript==True:
						if feature.upper() == 'TRANSCRIPT':
							feature_key = (seqname, feature, start, end, id, parent)
				# Only check for duplicates if this feature is being tracked
				if feature_key is not None:
					if feature_key in seen_features:
						feature_duplicates.append((line_num, feature_key))
					else:
						seen_features.add(feature_key)
				if len(parts) < 9:
					errors.append(f"Malformed GFF3 line at line {line_num}: {line.strip()}")
					continue

				contig, _, feature_type, start, end, *_ = parts
				start, end = int(start), int(end)

				# Check if contig name exists in FASTA
				if contig not in fasta_lengths:
					errors.append(f"Line {line_num}: Contig '{contig}' in GFF3 not found in FASTA.")
					continue

				# Check if gene's end position is valid
				if feature_type.lower() == "gene":
					if end > fasta_lengths[contig]:
						errors.append(f"Line {line_num}: Gene end position {end} exceeds length of contig '{contig}' ({fasta_lengths[contig]}).")

					# Check start and end positions
					if start > end or start <= 0 or end <= 0:
						errors.append(f"Line {line_num}: Invalid start ({start}) or end ({end}) position for contig '{contig}'. Start must be < end and both > 0.")
		if len(exact_duplicates)!=0:
			errors.append(f"Found exact duplicated lines in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if len(feature_duplicates)!=0:
			errors.append(f"Found duplicate gene/ transcript level features in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if id_counter!=0:
			errors.append(f"Missing specified child level attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")
		if parent_counter!=0:
			errors.append(f"Missing specified child parent linker in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")
		if parent_id_counter!=0:
			errors.append(f"Missing specified parent attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")

		# Collect FASTA invalid character errors
		for header, invalid_list in fasta_invalid_chars.items():
			for pos, char in invalid_list:
				warnings.append(f"FASTA header '{header}': Invalid character '{char}' at position {pos}.")

		# Write warnings to the output file
		if warnings:
			with open(warning_file, "w") as warning_out:
				message.append('There are some warnings to look at')
				warning_out.write("\n".join(warnings))

		# Write errors to the output file
		if errors:
			with open(error_file, "w") as error_out:
				error_out.write("\n".join(errors))

	else:  # dealing with compressed gff file
		with gzip.open(gff3_file, "rt") as gff:
			for line_num, line in enumerate(gff, 1):
				id = None
				parent = None
				if line.startswith("#") or not line.strip():
					continue

				parts = line.strip().split("\t")
				attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
				id = attributes.get(child_attribute, "")
				parent_id = attributes.get(parent_attribute_correct, "")
				parent = attributes.get(parent_attribute, "")
				# Extract fields
				seqname, source, feature, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
				# check for missing attributes
				if has_mrna == True:
					if feature.upper() in ('MRNA') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('MRNA', 'CDS') and parent == "":
						parent_counter += 1
				elif has_mrna == False and has_transcript == True:
					if feature.upper() in ('TRANSCRIPT') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('TRANSCRIPT', 'CDS') and parent == "":
						parent_counter += 1
				# Check exact duplicates
				full_line = line.strip()
				if full_line in seen_exact:
					exact_duplicates.append((line_num, full_line))
				else:
					seen_exact.add(full_line)
					# Check feature duplicates
					# Reset feature_key for each iteration
					feature_key = None
					if feature.upper() == 'GENE':
						feature_key = (seqname, feature, start, end, parent_id)
					if has_mrna == True:
						if feature.upper() == 'MRNA':
							feature_key = (seqname, feature, start, end, id, parent)
					elif has_transcript == True:
						if feature.upper() == 'TRANSCRIPT':
							feature_key = (seqname, feature, start, end, id, parent)
				# Only check for duplicates if this feature is being tracked
				if feature_key is not None:
					if feature_key in seen_features:
						feature_duplicates.append((line_num, feature_key))
					else:
						seen_features.add(feature_key)
				if len(parts) < 9:
					errors.append(f"Malformed GFF3 line at line {line_num}: {line.strip()}")
					continue

				contig, _, feature_type, start, end, *_ = parts
				start, end = int(start), int(end)

				# Check if contig name exists in FASTA
				if contig not in fasta_lengths:
					errors.append(f"Line {line_num}: Contig '{contig}' in GFF3 not found in FASTA.")
					continue

				# Check if gene's end position is valid
				if feature_type.lower() == "gene":
					if end > fasta_lengths[contig]:
						errors.append(f"Line {line_num}: Gene end position {end} exceeds length of contig '{contig}' ({fasta_lengths[contig]}).")

					# Step 3: Check start and end positions
					if start > end or start <= 0 or end <= 0:
						errors.append(f"Line {line_num}: Invalid start ({start}) or end ({end}) position for contig '{contig}'. Start must be < end and both > 0.")
		if len(exact_duplicates)!=0:
			errors.append(f"Found exact duplicated lines in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if len(feature_duplicates)!=0:
			errors.append(f"Found duplicate gene/ transcript level features in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if id_counter!=0:
			errors.append(f"Missing specified child level attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")
		if parent_counter!=0:
			errors.append(f"Missing specified child parent linker in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")
		if parent_id_counter != 0:
			errors.append(f"Missing specified parent attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")

		# Collect FASTA invalid character errors
		for header, invalid_list in fasta_invalid_chars.items():
			for pos, char in invalid_list:
				warnings.append(f"FASTA header '{header}': Invalid character '{char}' at position {pos}.")

		# Write warnings to the output file
		if warnings:
			with open(warning_file, "w") as warning_out:
				message.append('There are some warnings to look at')
				warning_out.write("\n".join(warnings))

		# Write errors to the output file
		if errors:
			with open(error_file, "w") as error_out:
				error_out.write("\n".join(errors))
	return errors, message

#Parsing CDS FASTA file
def parse_cds_fasta(fasta_file):
	fasta_lengths = {}
	fasta_invalid_chars = {}
	if fasta_file[-2:].lower() != 'gz':#dealing with uncompressed FASTA file
		with open(fasta_file, "r") as f:
			header = None
			sequence = []
			for line_num, line in enumerate(f, 1):
				line = line.strip()
				if line.startswith(">"):
					if header:
						seq = "".join(sequence)
						fasta_lengths[header] = len(seq)
						invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
						if invalid_chars:
							fasta_invalid_chars[header] = invalid_chars
					header = line[1:]  # Exclude '>'
					sequence = []
				else:
					sequence.append(line)
			if header:
				seq = "".join(sequence)
				fasta_lengths[header] = len(seq)
				invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
				if invalid_chars:
					fasta_invalid_chars[header] = invalid_chars
	else:#dealing with compressed FASTA file
		with gzip.open(fasta_file, "rt") as f:
			header = None
			sequence = []
			for line_num, line in enumerate(f, 1):
				line = line.strip()
				if line.startswith(">"):
					if header:
						seq = "".join(sequence)
						fasta_lengths[header] = len(seq)
						invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
						if invalid_chars:
							fasta_invalid_chars[header] = invalid_chars
					header = line[1:]  # Exclude '>'
					sequence = []
				else:
					sequence.append(line)
			if header:
				seq = "".join(sequence)
				fasta_lengths[header] = len(seq)
				invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACGTN"]
				if invalid_chars:
					fasta_invalid_chars[header] = invalid_chars
	return fasta_invalid_chars

#Parsing PEP FASTA file
def parse_pep_fasta(fasta_file):
	fasta_lengths = {}
	fasta_invalid_chars = {}
	if fasta_file[-2:].lower() != 'gz':#dealing with uncompressed FASTA file
		with open(fasta_file, "r") as f:
			header = None
			sequence = []
			for line_num, line in enumerate(f, 1):
				line = line.strip()
				if line.startswith(">"):
					if header:
						seq = "".join(sequence)
						fasta_lengths[header] = len(seq)
						invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACDEFGHIKLMNPQRSTVWY*"]
						if invalid_chars:
							fasta_invalid_chars[header] = invalid_chars
					header = line[1:]  # Exclude '>'
					sequence = []
				else:
					sequence.append(line)
			if header:
				seq = "".join(sequence)
				fasta_lengths[header] = len(seq)
				invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACDEFGHIKLMNPQRSTVWY*"]
				if invalid_chars:
					fasta_invalid_chars[header] = invalid_chars
	else:#dealing with compressed FASTA file
		with gzip.open(fasta_file, "rt") as f:
			header = None
			sequence = []
			for line_num, line in enumerate(f, 1):
				line = line.strip()
				if line.startswith(">"):
					if header:
						seq = "".join(sequence)
						fasta_lengths[header] = len(seq)
						invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACDEFGHIKLMNPQRSTVWY*"]
						if invalid_chars:
							fasta_invalid_chars[header] = invalid_chars
					header = line[1:]  # Exclude '>'
					sequence = []
				else:
					sequence.append(line)
			if header:
				seq = "".join(sequence)
				fasta_lengths[header] = len(seq)
				invalid_chars = [(i + 1, char) for i, char in enumerate(seq) if char.upper() not in "ACDEFGHIKLMNPQRSTVWY*"]
				if invalid_chars:
					fasta_invalid_chars[header] = invalid_chars
	return fasta_invalid_chars

#function to validate GFF3 and cleaned CDS FASTA files
def validate_gff3_and_cds(gff3_file, fasta_file, error_file, warning_file,child_attribute,parent_attribute,gemoma,parent_attribute_correct):
	id_counter = 0
	parent_id_counter = 0
	parent_counter=0
	seen_exact = set()
	seen_features = set()

	exact_duplicates = []
	feature_duplicates = []
	errors = []
	warnings = []
	message = []

	# Parse the FASTA file
	fasta_invalid_chars = parse_cds_fasta(fasta_file)
	headers = set()
	if fasta_file[-2:].lower() != 'gz':  # dealing with uncompressed gff3 file
		with open(fasta_file,'r') as f:
			for line in f:
				if line.startswith('>'):
					headers.add(line[1:].strip())
	else:
		with gzip.open(fasta_file,'rt') as f:
			for line in f:
				if line.startswith('>'):
					headers.add(line[1:].strip())

	if gff3_file[-2:].lower() != 'gz':
		with open(gff3_file,'r') as f:
			gff_lines = f.readlines()
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	else:
		with gzip.open(gff3_file,'rt') as f:
			gff_lines = f.readlines()
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_gene == False:
		errors.append(f"Level 1 feature 'gene' missing. Please process your GFF file with the AGAT wrapper script first.")
	if has_mrna == False and has_transcript == False:
		errors.append(f"Level 2 feature 'mRNA' or 'transcript' missing. Please process your GFF file with the AGAT wrapper script first.")
	if has_cds == False:
		errors.append(f"Level 3 feature 'CDS' missing. Please process your GFF file with the AGAT wrapper script first.")

	name_id_set = set()
	if has_gene:
		name_set = set()
		if gff3_file[-2:].lower() != 'gz':
			with open(gff3_file,'r') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if child_attribute in attributes:
								if attr.strip().startswith(str(child_attribute)+'='):
									name_value = attr.strip().split('=', 1)[1]
									name_set.add(name_value.strip())
		else:
			with gzip.open(gff3_file,'rt') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if child_attribute in attributes:
								if attr.strip().startswith(str(child_attribute)+'='):
									name_value = attr.strip().split('=', 1)[1]
									name_set.add(name_value.strip())
	else:
		name_set = set()
		if gff3_file[-2:].lower() != 'gz':
			with open(gff3_file,'r') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if attr.strip().startswith(str(child_attribute)+'='):
								name_value = attr.strip().split('=', 1)[1]
								name_set.add(name_value.strip())
		else:
			with gzip.open(gff3_file, 'rt') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if attr.strip().startswith(str(child_attribute)+'='):
								name_value = attr.strip().split('=', 1)[1]
								name_set.add(name_value.strip())

	if headers<=name_set:
		pass
	else:
		errors.append('FASTA file headers do not match with the specified child level GFF attributes fields.')
	# Parse the GFF3 file
	if gff3_file[-2:].lower() != 'gz':  # dealing with uncompressed gff file
		with open(gff3_file, "r") as gff:
			for line_num, line in enumerate(gff, 1):
				id = None
				parent = None
				if line.startswith("#") or not line.strip():
					continue

				parts = line.strip().split("\t")
				attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
				id = attributes.get(child_attribute, "")
				parent_id = attributes.get(parent_attribute_correct, "")
				parent = attributes.get(parent_attribute, "")
				# Extract fields
				seqname, source, feature, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
				# check for missing attributes
				if has_mrna == True:
					if feature.upper() in ('MRNA') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('MRNA', 'CDS') and parent == "":
						parent_counter += 1
				elif has_mrna == False and has_transcript == True:
					if feature.upper() in ('TRANSCRIPT') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('TRANSCRIPT', 'CDS') and parent == "":
						parent_counter += 1
				# Check exact duplicates
				full_line = line.strip()
				if full_line in seen_exact:
					exact_duplicates.append((line_num, full_line))
				else:
					seen_exact.add(full_line)
					# Check feature duplicates
					# Reset feature_key for each iteration
					feature_key = None
					if feature.upper() == 'GENE':
						feature_key = (seqname, feature, start, end, parent_id)
					if has_mrna == True:
						if feature.upper() == 'MRNA':
							feature_key = (seqname, feature, start, end, id, parent)
					elif has_transcript == True:
						if feature.upper() == 'TRANSCRIPT':
							feature_key = (seqname, feature, start, end, id, parent)
				# Only check for duplicates if this feature type is being tracked
				if feature_key is not None:
					if feature_key in seen_features:
						feature_duplicates.append((line_num, feature_key))
					else:
						seen_features.add(feature_key)
				if len(parts) < 9:
					errors.append(f"Malformed GFF3 line at line {line_num}: {line.strip()}")
					continue

				contig, _, feature_type, start, end, *_ = parts
				start, end = int(start), int(end)

				# Check if gene's end position is valid
				if feature_type.lower() == "gene":
					# Check start and end positions
					if start > end or start <= 0 or end <= 0:
						errors.append(f"Line {line_num}: Invalid start ({start}) or end ({end}) position for contig '{contig}'. Start must be < end and both > 0.")
		if len(exact_duplicates)!=0:
			errors.append(f"Found exact duplicated lines in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if len(feature_duplicates)!=0:
			errors.append(f"Found duplicate gene/ transcript level features in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if id_counter!=0:
			errors.append(f"Missing specified child attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
		if parent_counter!=0:
			errors.append(f"Missing specified child parent linker in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
		if parent_id_counter != 0:
			errors.append(f"Missing specified parent attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")

		# Collect FASTA invalid character errors
		for header, invalid_list in fasta_invalid_chars.items():
			for pos, char in invalid_list:
				warnings.append(f"FASTA header '{header}': Invalid character '{char}' at position {pos}.")

		# Write warnings to the output file
		if warnings:
			with open(warning_file, "w") as warning_out:
				message.append('There are some warnings to look at')
				warning_out.write("\n".join(warnings))

		# Write errors to the output file
		if errors:
			with open(error_file, "w") as error_out:
				error_out.write("\n".join(errors))

	else:
		with gzip.open(gff3_file, "rt") as gff:
			for line_num, line in enumerate(gff, 1):
				id = None
				parent = None
				if line.startswith("#") or not line.strip():
					continue

				parts = line.strip().split("\t")
				attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
				id = attributes.get(child_attribute, "")
				parent_id=attributes.get(parent_attribute_correct,"")
				parent = attributes.get(parent_attribute, "")
				# Extract fields
				seqname, source, feature, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
				# check for missing attributes
				if has_mrna == True:
					if feature.upper() in ('MRNA') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('MRNA', 'CDS') and parent == "":
						parent_counter += 1
				elif has_mrna == False and has_transcript == True:
					if feature.upper() in ('TRANSCRIPT') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('TRANSCRIPT', 'CDS') and parent == "":
						parent_counter += 1
				# Check exact duplicates
				full_line = line.strip()
				if full_line in seen_exact:
					exact_duplicates.append((line_num, full_line))
				else:
					seen_exact.add(full_line)
					# Check feature duplicates
					# Reset feature_key for each iteration
					feature_key = None
					if feature.upper() == 'GENE':
						feature_key = (seqname, feature, start, end, parent_id)
					if has_mrna == True:
						if feature.upper() == 'MRNA':
							feature_key = (seqname, feature, start, end, id, parent)
					elif has_transcript == True:
						if feature.upper() == 'TRANSCRIPT':
							feature_key = (seqname, feature, start, end, id, parent)
				# Only check for duplicates if this feature type is being tracked
				if feature_key is not None:
					if feature_key in seen_features:
						feature_duplicates.append((line_num, feature_key))
					else:
						seen_features.add(feature_key)
				if len(parts) < 9:
					errors.append(f"Malformed GFF3 line at line {line_num}: {line.strip()}")
					continue

				contig, _, feature_type, start, end, *_ = parts
				start, end = int(start), int(end)

				# Check if gene's end position is valid
				if feature_type.lower() == "gene":
					# Check start and end positions
					if start > end or start <= 0 or end <= 0:
						errors.append(f"Line {line_num}: Invalid start ({start}) or end ({end}) position for contig '{contig}'. Start must be < end and both > 0.")
		if len(exact_duplicates) != 0:
			errors.append(
				f"Found exact duplicated lines in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if len(feature_duplicates) != 0:
			errors.append(f"Found duplicate gene/ transcript level features in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if id_counter != 0:
			errors.append(f"Missing specified child attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
		if parent_counter != 0:
			errors.append(f"Missing specified child parent linker in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
		if parent_id_counter != 0:
			errors.append(f"Missing specified parent attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")

		# Collect FASTA invalid character errors
		for header, invalid_list in fasta_invalid_chars.items():
			for pos, char in invalid_list:
				warnings.append(f"FASTA header '{header}': Invalid character '{char}' at position {pos}.")

		# Write warnings to the output file
		if warnings:
			with open(warning_file, "w") as warning_out:
				message.append('There are some warnings to look at')
				warning_out.write("\n".join(warnings))

		# Write errors to the output file
		if errors:
			with open(error_file, "w") as error_out:
				error_out.write("\n".join(errors))

	return errors, message

#functions to validate GFF3 and PEP FASTA files
def validate_gff3_and_pep(gff3_file, fasta_file, error_file, warning_file,child_attribute,parent_attribute,gemoma,parent_attribute_correct):
	id_counter = 0
	parent_id_counter = 0
	parent_counter = 0
	seen_exact = set()
	seen_features = set()

	exact_duplicates = []
	feature_duplicates = []
	errors = []
	warnings = []
	message = []

	# Parse the FASTA file
	fasta_invalid_chars = parse_pep_fasta(fasta_file)
	headers = set()
	if fasta_file[-2:].lower() != 'gz':  # dealing with uncompressed gff3 file
		with open(fasta_file,'r') as f:
			for line in f:
				if line.startswith('>'):
					headers.add(line[1:].strip())
	else:
		with gzip.open(fasta_file,'rt') as f:
			for line in f:
				if line.startswith('>'):
					headers.add(line[1:].strip())

	if gff3_file[-2:].lower() != 'gz':
		with open(gff3_file,'r') as f:
			gff_lines = f.readlines()
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	else:
		with gzip.open(gff3_file,'rt') as f:
			gff_lines = f.readlines()
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_gene == False:
		errors.append(f"Level 1 feature 'gene' missing. Please process your GFF file with the AGAT wrapper script first.")
	if has_mrna == False and has_transcript == False:
		errors.append(f"Level 2 feature 'mRNA' or 'transcript' missing. Please process your GFF file with the AGAT wrapper script first.")
	if has_cds == False:
		errors.append(f"Level 3 feature 'CDS' missing. Please process your GFF file with the AGAT wrapper script first.")
	if has_gene:
		name_set = set()
		if gff3_file[-2:].lower() != 'gz':
			with open(gff3_file,'r') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if child_attribute in attributes:
								if attr.strip().startswith(str(child_attribute)+'='):
									name_value = attr.strip().split('=', 1)[1]
									name_set.add(name_value.strip())

		else:
			with gzip.open(gff3_file,'rt') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if child_attribute in attributes:
								if attr.strip().startswith(str(child_attribute)+'='):
									name_value = attr.strip().split('=', 1)[1]
									name_set.add(name_value.strip())

	else:
		name_set = set()
		if gff3_file[-2:].lower() != 'gz':
			with open(gff3_file,'r') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if attr.strip().startswith(str(child_attribute)+'='):
								name_value = attr.strip().split('=', 1)[1]
								name_set.add(name_value.strip())
		else:
			with gzip.open(gff3_file,'rt') as f:
				for line in f:
					if line.startswith('#') or not line.strip():
						continue
					parts = line.strip().split('\t')
					if len(parts) < 9:
						continue
					feature = parts[2].upper()
					if feature in ['MRNA', 'TRANSCRIPT']:
						attributes = parts[8]
						for attr in attributes.strip().split(';'):
							if attr.strip().startswith(str(child_attribute)+'='):
								name_value = attr.strip().split('=', 1)[1]
								name_set.add(name_value.strip())

	if headers<=name_set:
		pass
	else:
		errors.append('FASTA file headers do not match with the specified child level GFF attributes fields.')

	# Parse the GFF3 file
	if gff3_file[-2:].lower() != 'gz':#dealing with uncompressed gff3 file
		with open(gff3_file, "r") as gff:
			for line_num, line in enumerate(gff, 1):
				id = None
				parent = None
				if line.startswith("#") or not line.strip():
					continue

				parts = line.strip().split("\t")
				attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
				id = attributes.get(child_attribute, "")
				parent_id=attributes.get(parent_attribute_correct,"")
				parent = attributes.get(parent_attribute, "")
				# Extract fields
				seqname, source, feature, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
				# check for missing attributes
				if has_mrna == True:
					if feature.upper() in ('MRNA') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('MRNA', 'CDS') and parent == "":
						parent_counter += 1
				elif has_mrna == False and has_transcript == True:
					if feature.upper() in ('TRANSCRIPT') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('TRANSCRIPT', 'CDS') and parent == "":
						parent_counter += 1
				# Check exact duplicates
				full_line = line.strip()
				if full_line in seen_exact:
					exact_duplicates.append((line_num, full_line))
				else:
					seen_exact.add(full_line)
					# Check feature duplicates
					# Reset feature_key for each iteration
					feature_key = None
					if feature.upper() == 'GENE':
						feature_key = (seqname, feature, start, end, parent_id)
					if has_mrna == True:
						if feature.upper() == 'MRNA':
							feature_key = (seqname, feature, start, end, id, parent)
					elif has_transcript == True:
						if feature.upper() == 'TRANSCRIPT':
							feature_key = (seqname, feature, start, end, id, parent)
				# Only check for duplicates if we're tracking this feature type
				if feature_key is not None:
					if feature_key in seen_features:
						feature_duplicates.append((line_num, feature_key))
					else:
						seen_features.add(feature_key)
				if len(parts) < 9:
					errors.append(f"Malformed GFF3 line at line {line_num}: {line.strip()}")
					continue

				contig, _, feature_type, start, end, *_ = parts
				start, end = int(start), int(end)
				# Check if gene's end position is valid
				if feature_type.lower() == "gene":
					# Step 3: Check start and end positions
					if start > end or start <= 0 or end <= 0:
						errors.append(f"Line {line_num}: Invalid start ({start}) or end ({end}) position for contig '{contig}'. Start must be < end and both > 0.")
			if len(exact_duplicates) != 0:
				errors.append(f"Found exact duplicated lines in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
			if len(feature_duplicates) != 0:
				errors.append(f"Found duplicate gene/ transcript level features in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
			if id_counter != 0:
				errors.append(f"Missing specified child attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
			if parent_counter != 0:
				errors.append(f"Missing specified child parent linker in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
			if parent_id_counter != 0:
				errors.append(f"Missing specified parent attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")
		# Collect FASTA invalid character errors
		for header, invalid_list in fasta_invalid_chars.items():
			for pos, char in invalid_list:
				warnings.append(f"FASTA header '{header}': Invalid character '{char}' at position {pos}.")

		# Write warnings to the output file
		if warnings:
			with open(warning_file, "w") as warning_out:
				message.append('There are some warnings to look at')
				warning_out.write("\n".join(warnings))

		# Write errors to the output file
		if errors:
			with open(error_file, "w") as error_out:
				error_out.write("\n".join(errors))

	else:#dealing with compressed GFF3 file
		with gzip.open(gff3_file, "rt") as gff:
			for line_num, line in enumerate(gff, 1):
				id = None
				parent = None
				if line.startswith("#") or not line.strip():
					continue

				parts = line.strip().split("\t")
				attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
				id = attributes.get(child_attribute, "")
				parent_id=attributes.get(parent_attribute_correct,"")
				parent = attributes.get(parent_attribute, "")
				# Extract fields
				seqname, source, feature, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
				# check for missing attributes
				if has_mrna == True:
					if feature.upper() in ('MRNA') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('MRNA', 'CDS') and parent == "":
						parent_counter += 1
				elif has_mrna == False and has_transcript == True:
					if feature.upper() in ('TRANSCRIPT') and id == "":
						id_counter += 1
					if feature.upper() in ('GENE') and parent_id == "":
						parent_id_counter += 1
					if feature.upper() in ('TRANSCRIPT', 'CDS') and parent == "":
						parent_counter += 1
				# Check exact duplicates
				full_line = line.strip()
				if full_line in seen_exact:
					exact_duplicates.append((line_num, full_line))
				else:
					seen_exact.add(full_line)
					# Check feature duplicates
					# Reset feature_key for each iteration
					feature_key = None
					if feature.upper()=='GENE':
						feature_key = (seqname, feature, start, end, parent_id)
					if has_mrna == True:
						if feature.upper() == 'MRNA':
							feature_key = (seqname, feature, start, end, id, parent)
					elif has_transcript == True:
						if feature.upper() == 'TRANSCRIPT':
							feature_key = (seqname, feature, start, end, id, parent)
				# Only check for duplicates if we're tracking this feature type
				if feature_key is not None:
					if feature_key in seen_features:
						feature_duplicates.append((line_num, feature_key))
					else:
						seen_features.add(feature_key)
				if len(parts) < 9:
					errors.append(f"Malformed GFF3 line at line {line_num}: {line.strip()}")
					continue

				contig, _, feature_type, start, end, *_ = parts
				start, end = int(start), int(end)

				# Check if gene's end position is valid
				if feature_type.lower() == "gene":
					# Step 3: Check start and end positions
					if start > end or start <= 0 or end <= 0:
						errors.append(f"Line {line_num}: Invalid start ({start}) or end ({end}) position for contig '{contig}'. Start must be < end and both > 0.")
		if len(exact_duplicates)!=0:
			errors.append(f"Found exact duplicated lines in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if len(feature_duplicates)!=0:
			errors.append(f"Found duplicate gene/ transcript level features in the GFF files. Please process your GFF file with the AGAT wrapper script first.")
		if id_counter!=0:
			errors.append(f"Missing specified child attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
		if parent_counter!=0:
			errors.append(f"Missing specified child parent linker in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first.")
		if parent_id_counter!=0:
			errors.append(f"Missing specified parent attribute in the last column of the GFF file. Please process your GFF file with the AGAT wrapper script first to get a standardized attribute field name.")

		# Collect FASTA invalid character errors
		for header, invalid_list in fasta_invalid_chars.items():
			for pos, char in invalid_list:
				warnings.append(f"FASTA header '{header}': Invalid character '{char}' at position {pos}.")

		# Write warnings to the output file
		if warnings:
			with open(warning_file, "w") as warning_out:
				message.append('There are some warnings to look at')
				warning_out.write("\n".join(warnings))
		# Write errors to the output file
		if errors:
			with open(error_file, "w") as error_out:
				error_out.write("\n".join(errors))
	return errors, message

#functions to process validation of fasta and gff3 files simultaneously
def process_fasta_gff3_pair(each_fasta, gff3_input_file, validated_fasta_files, error_files, warning_files, fasta_type, logger,gff_config_parameters,gemoma):
	name_errors = []
	msg_dummy = []

	fasta_file_path = Path(each_fasta)
	# Local import safe for multiprocessing
	fasta_name = get_basename(str(each_fasta))

	# Determine output fasta path
	if str(each_fasta)[-2:].lower() != 'gz':
		output_fasta = os.path.join(validated_fasta_files, fasta_name[0] + '.fa')
	else:
		output_fasta = os.path.join(validated_fasta_files, fasta_name[0] + '.fa.gz')

	# Clean FASTA file
	cleaned_fasta_file = clean_fasta(str(each_fasta), output_fasta, fasta_type,logger)

	# Match the corresponding GFF3 file
	matched_gff3 = None
	for every in gff3_input_file:
		gff_name = get_basename(str(every))
		if fasta_name == gff_name:
			matched_gff3 = every
			gff_file_path = Path(every)
			break
	if 'all' in gff_config_parameters:
		child_attribute=gff_config_parameters['all']['child_attribute']
		parent_attribute=gff_config_parameters['all']['child_parent_linker'] #parent_attribute here means child parent linker - not changing the variable name becuase I have already used this a number of times in the valdiate functions
		parent_attribute_correct=gff_config_parameters['all']['parent_attribute']
	elif 'default' in gff_config_parameters:
		child_attribute = gff_config_parameters['default']['child_attribute']
		parent_attribute = gff_config_parameters['default']['child_parent_linker'] #parent_attribute here means child parent linker - not changing the variable name becuase I have already used this a number of times in the valdiate functions
		parent_attribute_correct = gff_config_parameters['default']['parent_attribute']
	else:
		child_attribute = gff_config_parameters[str(fasta_name[0])]['child_attribute']
		parent_attribute = gff_config_parameters[str(fasta_name[0])]['child_parent_linker']
		parent_attribute_correct = gff_config_parameters[str(fasta_name[0])]['parent_attribute']
	# Validate GFF3 file against cleaned FASTA
	if matched_gff3:
		error_output = os.path.join(error_files, fasta_name[0] + '_errors.txt')
		warning_output = os.path.join(warning_files, fasta_name[0] + '_warnings.txt')
		if fasta_type == 'fasta':
			errors, message = validate_gff3_and_fasta(str(matched_gff3), cleaned_fasta_file, error_output, warning_output,child_attribute,parent_attribute,gemoma,parent_attribute_correct)
		elif fasta_type == 'cds':
			errors, message = validate_gff3_and_cds(str(matched_gff3), cleaned_fasta_file, error_output, warning_output,child_attribute,parent_attribute,gemoma,parent_attribute_correct)
		elif fasta_type == 'pep':
			errors, message = validate_gff3_and_pep(str(matched_gff3), cleaned_fasta_file, error_output, warning_output,child_attribute,parent_attribute,gemoma,parent_attribute_correct)
		return errors, message  # Return errors collected
	else:
		name_error_output = os.path.join(error_files, fasta_name[0] + '_file_name_errors.txt')
		name_errors.append(f"FASTA GFF file name mismatch for '{each_fasta}'")
		with open(name_error_output,'w') as out:
				out.write(str(name_errors)+'\n')
		return name_errors, msg_dummy  # No matching GFF3 found

#master function to parallelize cleaning and validation of input GFF3 and FASTA files
def parallel_clean_and_validate(fasta_input_file, gff3_input_file, validated_fasta_files, error_files, warning_files, fasta_type,cores,logger,gff_config_parameters,gemoma):
	total_errors = []
	total_messages=[]
	total_messages_flat=[]
	i=0
	j=0
	#checking for file extensions
	allowed_patterns = [
		'*.gff', '*.gff.gz', '*.gff3', '*.gff3.gz','*.fa', '*.fa.gz', '*.fna', '*.fna.gz','*.fasta', '*.fasta.gz',
		'*.genome.fasta', '*.genome.fa', '*.genome.fasta.gz', '*.genome.fa.gz', '*.genome.fna',
		'*.genome.fna.gz','*.cds.fa', '*.cds.fasta','*.cds.fa.gz', '*.cds.fasta.gz', '*.cds.fna',
		'*.cds.fna.gz','*.pep.fa','*.pep.fa.gz','*.pep.fasta','*.pep.fasta.gz', '*.pep.fna',
		'*.pep.fna.gz'
	]
	for each in fasta_input_file:
		fasta_path = Path(each)
		fasta_ext_match = [pattern for pattern in allowed_patterns if fnmatch.fnmatch(fasta_path.name, pattern)]
		if not fasta_ext_match:
			logger.error(f"Extension of input FASTA file {each} does not match allowed extensions.")
			i+=1

	for each in gff3_input_file:
		gff_path = Path(each)
		gff_ext_match = [pattern for pattern in allowed_patterns if fnmatch.fnmatch(gff_path.name, pattern)]
		if not gff_ext_match:
			logger.error(f"Extension of input GFF file {each} does not match allowed extensions.")
			j+=1

	if i!=0 or j!=0:
		logger.error('Please check the extensions of your input FASTA and GFF files. They do not match the allowed extensions. Dupylicate analysis exiting...')
		sys.exit()

	# Dynamically detect number of CPU cores
	total_cores = cores
	num_fasta_files = len(fasta_input_file)

	if num_fasta_files == 0:
		raise ValueError("No FASTA input files provided.")
	max_workers = distribute_cores_single_level(total_cores, num_fasta_files)

	logger.info(f"Detected {total_cores} CPU cores.")
	logger.info(f"Parallelizing cleaning and validation using {max_workers} processes...")

	with ProcessPoolExecutor(max_workers=max_workers) as executor:
		futures = []
		for each_fasta in fasta_input_file:
			futures.append(executor.submit(process_fasta_gff3_pair,each_fasta,gff3_input_file,validated_fasta_files,error_files,warning_files,fasta_type,logger,gff_config_parameters,gemoma))
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures), desc="Organisms processed in the FASTA and GFF3 validation-clean step")
		for future in iterator:
			try:
				errors, messages = future.result()
				if errors:
					total_errors.append(errors)
				if messages:
					total_messages.append(messages)
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise
	if total_messages:
		for lists in total_messages:
			for each in lists:
				total_messages_flat.append(each)
		total_messages_flat=set(total_messages_flat)
		for msg in total_messages_flat:
			logger.warning(msg)
	logger.info(f"Completed validation of {len(fasta_input_file)} organism files.")
	return total_errors

#Loading gene IDs from FASTA file
def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	if multiple_fasta_file[-2:].lower() != 'gz':#dealing with uncompressed fasta file
		sequences = {}
		with open(multiple_fasta_file) as f:
			header = f.readline()[1:].strip()
			if " " in header:
				header = header.split(' ')[0]
				if "\t" in header:
					header = header.split('\t')[0]
			elif "\t" in header:
				header = header.split('\t')[0]
			seq = []
			line = f.readline()
			while line:
				if line[0] == '>':
					sequences.update({header: "".join(seq).upper()})
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
						if "\t" in header:
							header = header.split('\t')[0]
					elif "\t" in header:
						header = header.split('\t')[0]
					seq = []
				else:
					seq.append(line.strip())
				line = f.readline()
			sequences.update({header: "".join(seq).upper()})
		return sequences

	else:#dealing with compressed FASTA file(s)
		sequences = {}
		with gzip.open (multiple_fasta_file, "rt") as f:
			header = f.readline()[1:].strip()
			if " " in header:
				header = header.split(' ')[0]
				if "\t" in header:
					header = header.split('\t')[0]
			elif "\t" in header:
				header = header.split('\t')[0]
			seq = []
			line = f.readline()
			while line:
				if line[0] == '>':
					sequences.update({header: "".join(seq).upper()})
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
						if "\t" in header:
							header = header.split('\t')[0]
					elif "\t" in header:
						header = header.split('\t')[0]
					seq = []
				else:
					seq.append(line.strip())
				line = f.readline()
			sequences.update({header: "".join(seq).upper()})
		return sequences

#function for loading transcript information from GFF3 file
def load_transcript_information_from_gff3( gff3_input_file,process_pseudos,child_attribute,child_parent_linker):
	"""! @brief load all transcript information from gff3 file """
	# --- load all data from file --- #
	message = []
	gff_pseudos = set()

	if gff3_input_file[-2:].lower() != 'gz':  # dealing with uncompressed gff file
		information = []
		mrna_dict = {}  # To store mRNA ID -> details mapping
		exon_dict = {}  # To store exon parent -> list of exons mapping
		#code to check for pseudogenes
		with open(gff3_input_file, "r") as f:
			gff_lines = f.readlines()
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
					  if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
					  if not line.startswith('#') and len(line.split('\t')) >= 3)
			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)
		if process_pseudos == 'no':
			if has_mrna:  # checking for pseudogenes when mrna feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith((child_attribute)+'='):
										gene_id = attr[(len(child_attribute)+1):]
										gff_pseudos_genes.add(gene_id)
										break
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'MRNA':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[(len(child_attribute)+1):]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[(len(child_parent_linker)+1):]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			elif has_transcript:  # checking for pseudogenes when transcript feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[3:]
										gff_pseudos_genes.add(gene_id)
										break

					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'TRANSCRIPT':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[(len(child_attribute)+1):]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[(len(child_parent_linker)+1):]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[
										8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			else:  # checking for pseudogenes when no mrna, transcript feature is present and only cds feature is present
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE' or 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[(len(child_attribute)+1):]
										gff_pseudos.add(gene_id)
										break

		with open(gff3_input_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) > 5:
						if has_cds:
							if parts[2].upper() == 'CDS':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None  # Changed from False to None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_parent_linker)+'=') in subp:
												parent = subp.replace((str(child_parent_linker)+'='), "")

										# Check if parent is pseudogene AFTER finding parent
										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue  # Skip this CDS

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if (str(child_parent_linker)+'=') in parts[-1]:
											parent = str(parts[-1]).replace((str(child_parent_linker)+'='), "")

										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
						else:
							# Handle MRNA/TRANSCRIPT and EXON features
							if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
								if len(parts[-1]) > len((str(child_attribute)+'=')):
									if ";" in parts[-1]:
										mrna_id = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_attribute)+'=') in subp:
												mrna_id = subp.replace((str(child_attribute)+'='), "")

										if mrna_id and mrna_id not in gff_pseudos:
											mrna_dict[mrna_id] = {
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
											}

							elif parts[2].upper() == 'EXON':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_parent_linker)+'=') in subp:
												parent = subp.replace((str(child_parent_linker)+'='), "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if (str(child_parent_linker)+'=') in parts[-1]:
											parent = str(parts[-1]).replace((str(child_parent_linker)+'='), "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
				line = f.readline()

		if has_cds:
			# --- sort data by parent --- #
			sorted_data = {}
			for each in information:
				try:
					sorted_data[each['parent']].append(each)
				except KeyError:
					sorted_data.update({each['parent']: [each]})

			final_data = []
			for key in sorted_data.keys():
				if sorted_data[key][0]['orientation'] == '+':
					final_data.append(sorted(sorted_data[key], key=itemgetter('start')))
				else:
					final_data.append(sorted(sorted_data[key], key=itemgetter('start'))[::-1])
		else:
			# Process exon/transcript based annotation
			final_data = []
			for mrna_id, mrna_info in mrna_dict.items():
				if mrna_id in exon_dict:
					exons = exon_dict[mrna_id]

					# Sort the exons based on orientation
					if mrna_info['orientation'] == '+':
						sorted_exons = sorted(exons, key=itemgetter('start'))
					else:
						sorted_exons = sorted(exons, key=itemgetter('start'), reverse=True)
					final_data.append(sorted_exons)
		return final_data, message

	else:#dealing with compressed GFF3 file(s)
		information = []
		mrna_dict = {}  # To store mRNA ID -> details mapping
		exon_dict = {}  # To store exon parent -> list of exons mapping


		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
				   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
						 if not line.startswith('#') and len(line.split('\t')) >= 3)
			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
				  if not line.startswith('#') and len(line.split('\t')) >= 3)

		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			# collecting pseudogenes to skip cds of pseudogenes
			for line in gff_lines:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
						for attr in parts[8].split(';'):
							if attr.startswith(str(child_attribute)+'='):
								gene_id = attr[len(str(child_attribute)+'=')+1:]
								gff_pseudos.add(gene_id)
								break

			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)

		if process_pseudos == 'no':
			if has_mrna:  # checking for pseudogenes when mrna feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[len(str(child_attribute)+'=')+1:]
										gff_pseudos_genes.add(gene_id)
										break
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'MRNA':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[len(str(child_attribute)+'=')+1:]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[len(str(child_parent_linker)+'=')+1:]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			elif has_transcript:  # checking for pseudogenes when transcript feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[len(str(child_attribute)+'=')+1:]
										gff_pseudos_genes.add(gene_id)
										break

					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'TRANSCRIPT':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[len(str(child_attribute)+'=')+1:]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[len(str(child_parent_linker)+'=')+1:]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[
										8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			else:  # checking for pseudogenes when no mrna, transcript feature is present and only cds feature is present
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE' or 'PSEUDO=TRUE' in parts[
								8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										gene_id = attr[len(str(child_attribute)+'=')+1:]
										gff_pseudos.add(gene_id)
										break

		with gzip.open(gff3_input_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) > 5:
						if has_cds:
							if parts[2].upper() == 'CDS':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None  # Changed from False to None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_parent_linker)+'=' in subp:
												parent = subp.replace(str(child_parent_linker)+'=', "")

										# Check if parent is pseudogene AFTER finding parent
										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue  # Skip this CDS

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if str(child_parent_linker)+'=' in parts[-1]:
											parent = str(parts[-1]).replace(str(child_parent_linker)+'=', "")

										if parent:
											if parent in gff_pseudos:
												line = f.readline()
												continue

											information.append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
						else:
							# Handle MRNA/TRANSCRIPT and EXON features
							if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
								if len(parts[-1]) > len(str(child_attribute)+'='):
									if ";" in parts[-1]:
										mrna_id = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_attribute)+'=' in subp:
												mrna_id = subp.replace(str(child_attribute)+'=', "")

										if mrna_id and mrna_id not in gff_pseudos:
											mrna_dict[mrna_id] = {
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
											}

							elif parts[2].upper() == 'EXON':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_parent_linker)+'=' in subp:
												parent = subp.replace(str(child_parent_linker)+'=', "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if str(child_parent_linker)+'=' in parts[-1]:
											parent = str(parts[-1]).replace(str(child_parent_linker)+'=', "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
				line = f.readline()

		if has_cds:
			# --- sort data by parent --- #
			sorted_data = {}
			for each in information:
				try:
					sorted_data[each['parent']].append(each)
				except KeyError:
					sorted_data.update({each['parent']: [each]})

			final_data = []
			for key in sorted_data.keys():
				if sorted_data[key][0]['orientation'] == '+':
					final_data.append(sorted(sorted_data[key], key=itemgetter('start')))
				else:
					final_data.append(sorted(sorted_data[key], key=itemgetter('start'))[::-1])
		else:
			# Process exon/transcript based annotation
			final_data = []
			for mrna_id, mrna_info in mrna_dict.items():
				if mrna_id in exon_dict:
					exons = exon_dict[mrna_id]
					# Sort the exons based on orientation
					if mrna_info['orientation'] == '+':
						sorted_exons = sorted(exons, key=itemgetter('start'))
					else:
						sorted_exons = sorted(exons, key=itemgetter('start'), reverse=True)
					final_data.append(sorted_exons)
		return final_data, message

#function of constructing the CDS FASTA file
def construct_CDS_file( transcript_info, CDS_file, assembly, child_parent_linker):
	"""! @brief construct file with all sequences for translation """
	with open( CDS_file, "w" ) as out:
		for transcript in transcript_info:
			seq = []
			revcomp_status = False
			if transcript[0]['orientation'] == '-':
				revcomp_status = True
			for part in transcript:
				if revcomp_status:
					seq.append( revcomp( assembly[ part['chr'] ][ part['start']-1:part['end'] ] ) )
				else:
					seq.append( assembly[ part['chr'] ][ part['start']-1:part['end'] ] )
				# Get parent ID, handling both formats of CDS feature or mRNA/ exon features' presence
				parent_id = transcript[0]['parent']
				if str(child_parent_linker)+"=" in parent_id:
					parent_id = parent_id.replace(str(child_parent_linker)+"=", "")
			out.write( '>' + str(parent_id) + '\n' + "".join( seq ) + '\n' )

#function for constructing the reverse complement
def revcomp( seq ):
	"""! @brief constructs revcomp """
	new_seq = []
	dictionary = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N','a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n' }
	for nt in seq:
		try:
			new_seq.append( dictionary[ nt ] )
		except KeyError:
			new_seq.append( "N")
	return ''.join( new_seq[::-1] )

#function for translation to peptide sequences - part1
def translate(seq, genetic_code):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	seq = seq.upper()
	peptide = []
	for i in range(int(len(seq) / 3.0)):
		codon = seq[i * 3:i * 3 + 3]
		try:
			peptide.append(genetic_code[codon])
		except:
			peptide.append("*")
	peptide = "".join(peptide)

	if sum([peptide[0] != "M", "*" in peptide[:-1]]) > 0:
		peptide2 = []
		for i in range(int((len(seq) - 1) / 3.0)):
			codon = seq[1 + i * 3:1 + i * 3 + 3]
			try:
				peptide2.append(genetic_code[codon])
			except:
				peptide2.append("*")
		peptide2 = "".join(peptide2)

		if sum([peptide2[0] != "M", "*" in peptide2[:-1]]) > 0:
			peptide3 = []
			for i in range(int((len(seq) - 2) / 3.0)):
				codon = seq[2 + i * 3:2 + i * 3 + 3]
				try:
					peptide3.append(genetic_code[codon])
				except:
					peptide3.append("*")
			peptide3 = "".join(peptide3)

			pep_options = []
			if '*' in peptide:
				final_pep = peptide.split('*')[0] + '*'
				pep_options.append({'seq': final_pep,'cds': seq[0:len(final_pep) * 3], 'stopps': peptide.count('*'),'len': len(final_pep)})
			else:
				pep_options.append({'seq': peptide,'cds': seq[0:len(peptide) * 3],'stopps': peptide.count('*'),'len': len(peptide)})

			if "*" in peptide2:
				final_pep2 = peptide2.split('*')[0] + '*'
				pep_options.append({'seq': final_pep2,'cds': seq[1:1 + len(final_pep2) * 3], 'stopps': peptide2.count('*'),'len': len(final_pep2)})
			else:
				pep_options.append({'seq': peptide2,'cds': seq[1:1 + len(peptide2) * 3],'stopps': peptide2.count('*'),'len': len(peptide2)})

			if "*" in peptide3:
				final_pep3 = peptide3.split('*')[0] + '*'
				pep_options.append({'seq': final_pep3,'cds': seq[2:2 + len(final_pep3) * 3], 'stopps': peptide3.count('*'),'len': len(final_pep3)})
			else:
				pep_options.append({'seq': peptide3,'cds': seq[2:2 + len(peptide3) * 3],'stopps': peptide3.count('*'),'len': len(peptide3)})

			winner = sorted(pep_options, key=lambda x: x['len'])[-1]
			return {'peptide': winner['seq'], 'cds': winner['cds']}

		else:
			pep_options = []
			if '*' in peptide:
				final_pep = peptide.split('*')[0] + '*'
				pep_options.append({'seq': final_pep,'cds': seq[0:len(final_pep) * 3],  'stopps': peptide.count('*'),'len': len(final_pep)})
			else:
				pep_options.append({'seq': peptide,'cds': seq[0:len(peptide) * 3],'stopps': peptide.count('*'),'len': len(peptide)})

			if "*" in peptide2:
				final_pep2 = peptide2.split('*')[0] + '*'
				pep_options.append({'seq': final_pep2,'cds': seq[1:1 + len(final_pep2) * 3],  'stopps': peptide2.count('*'),'len': len(final_pep2)})
			else:
				pep_options.append({'seq': peptide2,'cds': seq[1:1 + len(peptide2) * 3],'stopps': peptide2.count('*'),'len': len(peptide2)})

			winner = sorted(pep_options, key=lambda x: x['len'])[-1]
			return {'peptide': winner['seq'], 'cds': winner['cds']}
	else:
		# Simple case: starts with M, no internal stops
		# Keep the stop codon if it exists
		final_pep = peptide  # Don't strip the stop codon
		return {'peptide': final_pep, 'cds': seq[0:len(final_pep) * 3]}

#function for translation to peptide sequences - part2
def transeq( input_file, output_file, strict_start, strict_end, peptide_to_cds_map ):
	"""! @brief run translation of coding sequences """
	message2 = []
	genetic_code = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q', 'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'}
	sequences = load_sequences( input_file )
	with open( output_file, "w" ) as out:
		keys = sequences.keys()
		for key in keys:
			seq = sequences[ key ]
			if len( seq ) > 9:
				result = translate(seq, genetic_code)
				peptide, cds = result['peptide'], result['cds']
				out.write( '>' + key + '\n' + peptide + '\n' )
				peptide_to_cds_map[key]=str(cds)
			else:
				message2.append (key + " - too short!")
	return message2, peptide_to_cds_map

#function for removing alternate transcripts from the peptide FASTA file
def pepclean(gff3_input_file, output_file, no_trans_pep, no_trans_cds, fasta_type, peptide_to_cds_map,child_attribute,child_parent_linker):
	no_gene_no_parent = False
	has_gene = False
	has_parent =False
	if gff3_input_file[-2:].lower() != 'gz':
		with open(gff3_input_file, "r") as f:
			gff_lines = f.readlines()
			# checking if gene feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)

	else:
		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			# checking if gene feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_mrna:
		coding_feature = 'MRNA'
	else:
		has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
		if has_transcript:
			coding_feature = 'TRANSCRIPT'
		else:
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			if has_cds:
				coding_feature = 'CDS'
			else:
				coding_feature = 'EXON'

	orgname = get_basename(output_file)[0]
	nogene_noparent_counter = 0
	if gff3_input_file[-2:].lower() != 'gz':  # uncompressed gff file
		transcripts_per_gene = {}
		with open(gff3_input_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					no_gene_no_parent = False  # Reset for each line
					parts = line.strip().split('\t')
					if len(parts) > 2:
						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker) + '=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								for each in partsnew:
									pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
									if re.match(pattern_par, each):
										partsnew1 = str(each).replace(str(child_parent_linker) + '=', "")
							for every in partsnew:
								pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
								if re.match(pattern_ID, every):
									partsnew0 = str(every).replace(str(child_attribute) + '=', "")
							try:
								transcripts_per_gene[partsnew1].append(partsnew0)
							except KeyError:
								transcripts_per_gene.update({partsnew1: [partsnew0]})
				line = f.readline()

	else:#compressed gff file
		transcripts_per_gene = {}
		with gzip.open(gff3_input_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					no_gene_no_parent = False  # Reset for each line
					parts = line.strip().split('\t')
					if len(parts) > 2:

						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker)+'=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								for each in partsnew:
									pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
									if re.match(pattern_par, each):
										partsnew1 = str(each).replace(str(child_parent_linker)+'=', "")
							for every in partsnew:
								pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
								if re.match(pattern_ID, every):
									partsnew0 = str(every).replace(str(child_attribute)+'=', "")
							try:
								transcripts_per_gene[partsnew1].append(partsnew0)
							except KeyError:
								transcripts_per_gene.update({partsnew1: [partsnew0]})
				line = f.readline()

	if nogene_noparent_counter != 0:
		no_gene_no_parent = False
	gene_names = list(transcripts_per_gene.keys())
	#dealing with uncompressed pep file
	if output_file[-2:].lower() != 'gz':
		with open(output_file, "r") as f:
			# Only open CDS output file if peptide_to_cds_map is not empty
			if peptide_to_cds_map:
				with open(no_trans_pep, "w") as out, open(no_trans_cds, "w") as cds_out:
					line = f.readlines()
					pep_dict = {}
					for each in line:
						if '>' in each:
							ind = line.index(each)
							pep_dict.update({str(each).replace('>', '').replace('\n', ''): line[ind + 1]})
					transcripts = list(pep_dict.keys())
					if no_gene_no_parent == False:
						for gene in gene_names:
							trans_length = []
							for trans in transcripts_per_gene[gene]:
								for each in transcripts:
									if trans == each:
										if len(transcripts_per_gene[gene]) < 2:
											out.write('>' + str(gene) + '\n' + str(pep_dict[each]))
											cds_out.write('>' + str(gene) + '\n' + str(peptide_to_cds_map[each]) + '\n')
										else:
											trans_length.append((each, pep_dict[each]))
							if len(trans_length) == 0:
								pass
							else:
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(gene) + "\n" + str(seq))
								cds_out.write('>' + str(gene) + '\n' + str(peptide_to_cds_map[best_trans]) + '\n')
					elif no_gene_no_parent == True:
						for gene in gene_names:
							trans_length = []
							for trans in transcripts_per_gene[gene]:
								for each in transcripts:
									if trans == each:
										if len(transcripts_per_gene[gene]) < 2:
											out.write('>' + str(each) + '\n' + str(pep_dict[each]))
											cds_out.write('>' + str(each) + '\n' + str(peptide_to_cds_map[each]) + '\n')
										else:
											trans_length.append((each, pep_dict[each]))
							if len(trans_length) == 0:
								pass
							else:
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(best_trans) + "\n" + str(seq))
								cds_out.write('>' + str(best_trans) + '\n' + str(peptide_to_cds_map[best_trans]) + '\n')
			else:
				# Only write peptide output when no CDS mapping available
				with open(no_trans_pep, "w") as out:
					line = f.readlines()
					pep_dict = {}
					for each in line:
						if '>' in each:
							ind = line.index(each)
							pep_dict.update({str(each).replace('>', '').replace('\n', ''): line[ind + 1]})
					transcripts = list(pep_dict.keys())
					if no_gene_no_parent == False:
						for gene in gene_names:
							trans_length = []
							for trans in transcripts_per_gene[gene]:
								for each in transcripts:
									if trans == each:
										if len(transcripts_per_gene[gene]) < 2:
											out.write('>' + str(gene) + '\n' + str(pep_dict[each]))
										else:
											trans_length.append((each, pep_dict[each]))
							if len(trans_length) == 0:
								pass
							else:
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(gene) + "\n" + str(seq))
					elif no_gene_no_parent == True:
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
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(best_trans) + "\n" + str(seq))

	else:#dealing with compressed pep file
		with gzip.open(output_file, "rt") as f:
			# Only open CDS output file if peptide_to_cds_map is not empty
			if peptide_to_cds_map:
				with open(no_trans_pep, "w") as out, open(no_trans_cds, "w") as cds_out:
					line = f.readlines()
					pep_dict = {}
					for each in line:
						if '>' in each:
							ind = line.index(each)
							pep_dict.update({str(each).replace('>', '').replace('\n', ''): line[ind + 1]})
					transcripts = list(pep_dict.keys())
					if no_gene_no_parent == False:
						for gene in gene_names:
							trans_length = []
							for trans in transcripts_per_gene[gene]:
								for each in transcripts:
									if trans == each:
										if len(transcripts_per_gene[gene]) < 2:
											out.write('>' + str(gene) + '\n' + str(pep_dict[each]))
											cds_out.write('>' + str(gene) + '\n' + str(peptide_to_cds_map[each]) + '\n')
										else:
											trans_length.append((each, pep_dict[each]))
							if len(trans_length) == 0:
								pass
							else:
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(gene) + "\n" + str(seq))
								cds_out.write('>' + str(gene) + '\n' + str(peptide_to_cds_map[best_trans]) + '\n')
					elif no_gene_no_parent == True:
						for gene in gene_names:
							trans_length = []
							for trans in transcripts_per_gene[gene]:
								for each in transcripts:
									if trans == each:
										if len(transcripts_per_gene[gene]) < 2:
											out.write('>' + str(each) + '\n' + str(pep_dict[each]))
											cds_out.write('>' + str(each) + '\n' + str(peptide_to_cds_map[each]) + '\n')
										else:
											trans_length.append((each, pep_dict[each]))
							if len(trans_length) == 0:
								pass
							else:
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(best_trans) + "\n" + str(seq))
								cds_out.write('>' + str(best_trans) + '\n' + str(peptide_to_cds_map[best_trans]) + '\n')
			else:
				# Only write peptide output when no CDS mapping available
				with open(no_trans_pep, "w") as out:
					line = f.readlines()
					pep_dict = {}
					for each in line:
						if '>' in each:
							ind = line.index(each)
							pep_dict.update({str(each).replace('>', '').replace('\n', ''): line[ind + 1]})
					transcripts = list(pep_dict.keys())
					if no_gene_no_parent == False:
						for gene in gene_names:
							trans_length = []
							for trans in transcripts_per_gene[gene]:
								for each in transcripts:
									if trans == each:
										if len(transcripts_per_gene[gene]) < 2:
											out.write('>' + str(gene) + '\n' + str(pep_dict[each]))
										else:
											trans_length.append((each, pep_dict[each]))
							if len(trans_length) == 0:
								pass
							else:
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(gene) + "\n" + str(seq))
					elif no_gene_no_parent == True:
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
								best_trans, seq = max(trans_length, key=lambda x: len(x[1]))
								out.write('>' + str(best_trans) + "\n" + str(seq))

#function for making databases for local alignment
def makedb(no_trans_pep, database, tmp_dir, tool, makeblastdb, blastp, diamond, mmseqs2,logger):
	if tool=='blast':
		cmd = makeblastdb + ' -in ' + no_trans_pep + ' -dbtype prot -parse_seqids' + ' -out ' + database
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
	#for diamond
	elif tool == 'diamond':
		cmd = diamond + ' makedb --in ' + no_trans_pep + ' --db ' + database + ' --quiet '
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
	#for mmseqs2
	elif tool == 'mmseqs2':
		cmd = mmseqs2 + ' createdb ' + no_trans_pep + ' ' + database
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()

#functions to trim an alignment file
def load_alignment(aln_file):
	#load alignment from input file
	sequences, names = {}, []
	with open(aln_file) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				sequences.update({header: "".join(seq)})
				names.append(header)
				header = line.strip()[1:]
				seq = []
			else:
				seq.append(line.strip())
			line = f.readline()
		sequences.update({header: "".join(seq)})
		names.append(header)
	return sequences, names

#function to check the multiple sequence alignment
def aln_len_check(alignment):
	#brief check if all sequences in alignment have the same length
	lengths = []
	for key in list(alignment.keys()):
		lengths.append(len(alignment[key]))
	if len(list(set(lengths))) != 1:  # not all sequences have same length
		for key in list(alignment.keys()):
			sys.stdout.write(key + " - len: " + str(len(alignment[key])) + "\n")
		sys.stdout.flush()
		sys.exit("ERROR: sequences do not have the same length.")

#function for cleaning of multiple seuqence alignment file
def alignment_trimming(aln_file, cln_aln_file, occupancy, sorting):
	#brief remove all alignment columns with insufficient occupancy
	alignment, names = load_alignment(aln_file)

	if sorting:  # alphanumerical sorting of names if activated
		names = sorted(names)

	# --- if there is an alignment (expected case) --- #
	if len(names) > 0:
		aln_len_check(alignment)  # perform alignment check (same sequence lengths?)

		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate(list(alignment.values())[0]):
			counter = 0
			for key in names:
				if alignment[key][idx] != "-":
					counter += 1
			if counter / float(len(list(alignment.keys()))) >= occupancy:
				valid_index.append(idx)

		# --- generate new sequences --- #
		with open(cln_aln_file, "w") as out:
			for key in names:
				seq = alignment[key]
				new_seq = []
				for idx in valid_index:
					new_seq.append(seq[idx])
				new_seq = "".join(new_seq)
				if new_seq.count('-') == len(new_seq):  # exclude sequences that contain only gaps after trimming
					sys.stdout.write("WARNING: only gaps remaining in sequence - " + key + " (sequence not included in output)\n")
					sys.stdout.flush()
				else:
					out.write(">" + key + '\n' + new_seq + '\n')

#function to write FASTA headers and sequences from a dictionary into a file
def write_fasta(query_sequences, ref_sequences, ids, output_path,query_gene,outgroup):
	with open(output_path, "w") as out:
		for i, gene_id in enumerate(ids):
			if i == 0:
				if gene_id in query_sequences:
					# Apply all replacements at once
					safe_gene_id = gene_id
					for char, replacement in replacements.items():
						safe_gene_id = safe_gene_id.replace(char, replacement)
					out.write(f">'{safe_gene_id}'\n")
					query_safe_gene_id = safe_gene_id
					out.write(query_sequences[gene_id] + "\n")
			else:
				if gene_id in ref_sequences:
					safe_gene_id = gene_id
					for char, replacement in replacements.items():
						safe_gene_id = safe_gene_id.replace(char, replacement)
					out.write(f">'{safe_gene_id}_ref'\n")
					if gene_id == outgroup:
						outgroup_safe_gene_id = str(safe_gene_id)+'_ref'
					out.write(ref_sequences[gene_id] + "\n")
	return query_safe_gene_id, outgroup_safe_gene_id

#function to calculate pairwise topological distance from the phylogenetic tree file
def topological_distance(tree, node1_name, node2_name):
	node1 = tree & node1_name
	node2 = tree & node2_name
	lca = tree.get_common_ancestor(node1, node2)

	dist1 = 0
	current = node1
	while current != lca:
		current = current.up
		dist1 += 1

	dist2 = 0
	current = node2
	while current != lca:
		current = current.up
		dist2 += 1

	return dist1 + dist2

#v13 worker function to do MAFFT alignment, cleaning, tree building and phylogenetic distance(s) calculations
def do_mafft_clean_tree(query_gene,test_group,query_sequences, ref_sequences,tmp_dir,mafft,occupancy,fasttree,outgroup, possible_orthologs,logger,top_bit_score,outgroup_bit_score):
	if (len(test_group)-1)>possible_orthologs:
		candidate_orthologs = test_group[1:possible_orthologs+1]
	elif (len(test_group)-1)<=possible_orthologs:
		candidate_orthologs = test_group[1:]

	#code block to calculate score reduction factor for ortholog confidence determination
	score_reduction_factor = 1
	percent_bit_score_diff = ((top_bit_score-outgroup_bit_score)/top_bit_score)*100
	if percent_bit_score_diff <= 10:
		score_reduction_factor = float(score_reduction_factor/2)
	ref_len = (len(test_group) - 1) # total number of ref genes in the test group
	top_rank_hits = math.floor((ref_len/3)+0.5) # calculating the number of top ranking hits in the test group for finding the placement of the ortholog identified later

	try:
		# Write FASTA
		with tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False, suffix=".fasta", mode="w") as tmp_fasta:
			fasta_path = tmp_fasta.name
		modified_query_gene_id, modified_outgroup_gene_id = write_fasta(query_sequences, ref_sequences, test_group,fasta_path, query_gene, outgroup)

		# MAFFT alignment
		aligned_path = fasta_path + "_aligned.fasta"
		cmd = f"{mafft} --quiet --globalpair --maxiterate 1000 --amino {fasta_path} > {aligned_path}"
		subprocess.run(cmd, shell=True)

		# Trimming
		cleaned_aligned_path = fasta_path + '_cleaned_aligned.fasta'
		alignment_trimming(aligned_path, cleaned_aligned_path, occupancy, sorting=True)

		# Build tree
		tree_path = fasta_path + '.nwk'
		cmd = f"{fasttree} -quiet -nopr -wag -boot 10000 {cleaned_aligned_path}"
		with open(tree_path, 'w') as tree_file:
			result = subprocess.run(cmd.split(), stdout=tree_file, stderr=subprocess.PIPE, text=True)

			# Filter stderr to remove only the "Ignored unknown character" messages
			if result.stderr:
				filtered_lines = []
				for line in result.stderr.split('\n'):
					# Skip lines that match the exact FastTree character warning pattern
					if not (line.startswith('Ignored unknown character') and '(seen' in line and 'times)' in line):
						filtered_lines.append(line)

				filtered_stderr = '\n'.join(filtered_lines).strip() # to suppress the ignoreed unknown character messages
				if filtered_stderr:  # Only print if there are other error messages
					print(filtered_stderr, file=sys.stderr)

		# Use DendroPy to compute distances
		if os.path.exists(tree_path) and os.path.getsize(tree_path) != 0:
			taxon_namespace = TaxonNamespace()
			tree = Tree.get_from_path(tree_path, schema="newick", taxon_namespace=taxon_namespace)

			outseq = tree.find_node_with_taxon_label(modified_outgroup_gene_id)
			if outseq is None:
				return query_gene, None, candidate_orthologs, score_reduction_factor
			else:
				tree.to_outgroup_position(outseq)
				tree.is_rooted = True
				pdm = PhylogeneticDistanceMatrix.from_tree(tree)

			try:
				query_taxon = tree.find_node_with_taxon_label(modified_query_gene_id).taxon
			except AttributeError:
				sys.stdout.write(f"Query gene {modified_query_gene_id} not found in tree.\n")
				sys.stdout.flush()
				return query_gene, None, candidate_orthologs, score_reduction_factor

			# Collect edge and patristic distances
			data = []
			for taxon in tree.taxon_namespace:
				if taxon.label != modified_query_gene_id and taxon.label != modified_outgroup_gene_id: # the query gene and outgroup gene are avoided here to skip self pairing of query gene and to avoid inclusion of outgroup in the phylogenetic distance calculation
					try:
						edge_distance = pdm.path_edge_count(query_taxon, taxon)
						patr_distance = pdm.patristic_distance(query_taxon, taxon)
						data.append({
							'gene': taxon.label,
							'edge': edge_distance,
							'patr': patr_distance
						})
					except Exception as e:
						sys.stdout.write(f"Distance calc error: {query_gene} ↔ {taxon.label}: {e}\n")
						sys.stdout.flush()
			# Reverse mapping for decoding to get back original gene names
			reverse_replacements = {v: k for k, v in replacements.items()}
			if data:
				# Sort by edge, then patristic distance
				closest_hit = sorted(data, key=itemgetter('edge', 'patr'))[0]['gene']
				ortholog_hit_with_placeholder = closest_hit.strip("'").split('_ref')[0]
				for replacement, char in reverse_replacements.items():
					ortholog_hit = ortholog_hit_with_placeholder.replace(replacement,char)  # replacing placeholders back with problematic characters to get back original gene name for further processing
				query_gene = query_gene.strip("'")
				ind_orth = test_group.index(ortholog_hit)
				if ind_orth > top_rank_hits:
					score_reduction_factor = float(score_reduction_factor / 2)
				return query_gene, ortholog_hit, candidate_orthologs, score_reduction_factor
			else:
				query_gene = query_gene.strip("'")
				return query_gene, None, candidate_orthologs, score_reduction_factor
		else:
			sys.stdout.write(f"Tree file {tree_path} not found.\n")
			query_gene = query_gene.strip("'")
			sys.stdout.write(f"{query_gene}\tNone")
			return query_gene, None, candidate_orthologs, score_reduction_factor
	finally:
		# Clean up of the tmp fasta, alignment and tree files
		for path in [fasta_path, aligned_path, cleaned_aligned_path, tree_path]:
			try:
				os.remove(path)
			except Exception:
				pass

#v13 function for collecting best hits from forward blast/ diamond/ mmseqs2 result files
def besthits(result, best_hits_file,bitscore_dic,pep_query,pep_ref,number_of_hits,tmp_dir,mafft,occupancy,fasttree,cpu,logger,fwd_similarity_cutoff,score_ratio_cutoff, pos_dic_query, pos_dic_ref,orthologs_file, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree):
	best_hits = {}
	all_fwd_hits = {}
	top_local_hits = {}
	all_ref_genes_bitscores = {}
	max_bitscore = 0.0
	with open(result, "r") as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] in bitscore_dic.keys():
				max_bitscore = float(bitscore_dic[parts[0]])#bit score of the query gene's hit against itself in self blast, making this the maximum possible bit score or limit of biological similarity
			elif parts[0] not in bitscore_dic.keys():
				line = f.readline()
				continue
			if float(parts[2]) > fwd_similarity_cutoff:#fwd hit should have some amount of good similarity percentage
				if float((float(parts[-1])/max_bitscore)) > score_ratio_cutoff:#fwd hit's bit score should be significant when compared to the query's self hit bit score to be biologically significant
					try:
						best_hits[parts[0]].append({'score': float((float(parts[-1])/max_bitscore)), 'evalue':float(parts[-2]), 'line': line})
					except:
						best_hits.update({parts[0]: [{'score': float((float(parts[-1])/max_bitscore)), 'evalue':float(parts[-2]), 'line': line}]})
			#making a dictionary of ref genes and the bit scores of the forward hits
			all_ref_genes_bitscores[str(parts[0])+'\t'+str(parts[1])] = float(parts[-1])
			#making a dictionary of all forward hits for each query gene for the sake of comprehensive phylogeny step w/o filters
			try:
				all_fwd_hits[parts[0]].append({'score': float((float(parts[-1]) / max_bitscore)), 'evalue': float(parts[-2]), 'line': line})
			except:
				all_fwd_hits.update({parts[0]: [{'score': float((float(parts[-1]) / max_bitscore)), 'evalue': float(parts[-2]), 'line': line}]})
			line = f.readline()
	final_sorted_best_hits = {}
	for key in sorted(best_hits.keys()):
		data = sorted(best_hits[key], key=lambda x: (-x['score'], x['evalue']))#sorting the hits based on descending order of bit scores; in case of bit score tie, go for second level sorting of hits based on ascending order of e values
		final_sorted_best_hits.update({key: data})
	#making a sorted dictionary of all the forward hits without filters for the phylogeny step
	all_sorted_hits = {}
	for key in sorted(all_fwd_hits.keys()):
		data = sorted(all_fwd_hits[key], key=lambda x: (-x['score'], x['evalue']))#sorting the hits based on descending order of bit scores; in case of bit score tie, go for second level sorting of hits based on ascending order of e values
		all_sorted_hits.update({key: data})
	#making a dictionary of top BLAST hits for each query gene
	for query_gene, hits in sorted(final_sorted_best_hits.items()):
		parts = hits[0]['line'].strip().split('\t')
		top_local_hits[query_gene]=parts[1]
	#making a dictionary of the query pep sequences and the headers
	query_sequences = {}
	with open(pep_query) as f:
		gene_id = None
		seq_lines = []
		for line in f:
			line = line.strip()
			if line.startswith(">"):
				if gene_id:
					query_sequences[gene_id] = "".join(seq_lines)
				gene_id = line[1:].split()[0]
				seq_lines = []
			else:
				seq_lines.append(line)
		if gene_id:
			query_sequences[gene_id] = "".join(seq_lines)
	# making a dictionary of the reference pep sequences and the headers
	ref_sequences = {}
	with open(pep_ref) as f:
		gene_id = None
		seq_lines = []
		for line in f:
			line = line.strip()
			if line.startswith(">"):
				if gene_id:
					ref_sequences[gene_id] = "".join(seq_lines)
				gene_id = line[1:].split()[0]
				seq_lines = []
			else:
				seq_lines.append(line)
		if gene_id:
			ref_sequences[gene_id] = "".join(seq_lines)
	#combining the two dictionaries
	seqs_dic = query_sequences | ref_sequences
	best_hits_dic = {}
	logger.info(f"Allowing up to {cpu} inner workers per organism's ortholog detection steps.")
	#collecting top n hits along with the query iteratively into a list and parallelizing further steps with a nested parallelization approach using a worker function
	with ProcessPoolExecutor(max_workers=cpu)as executor:
		futures=[]
		for query_gene, hits in sorted(final_sorted_best_hits.items()):
			test_group = []
			test_group.append(query_gene)
			if len(final_sorted_best_hits[query_gene]) > number_of_hits:
				top_hits = hits[:number_of_hits]
			elif len(final_sorted_best_hits[query_gene]) <= number_of_hits:
				top_hits = hits[:len(final_sorted_best_hits[query_gene])]
			if len(top_hits) == 1:
				parts = hits[0]['line'].strip().split('\t')
				similarity_diff = ((float(parts[2]) - fwd_similarity_cutoff) / float(parts[2])) * 100
				score_ratio_diff = ((float(parts[-1]) - score_ratio_cutoff) / float(parts[-1])) * 100
				# High confidence only if BOTH conditions are met
				if similarity_diff > 30 and score_ratio_diff > 30:
					best_hits_dic[query_gene] = [parts[1], 1, None, 'High confidence ortholog']
				else:
					# Run synteny check for all other cases
					synteny_found = False
					sublist = [parts[0]]
					anchor_list = [parts[1]]
					synteny = synteny_checker(anchor_list, sublist, pos_dic_query, pos_dic_ref,synteny_cutoff, flank_number, best_side, top_local_hits)
					if synteny == 'yes':
						best_hits_dic[query_gene] = [parts[1], 0.5, None, 'Moderate confidence ortholog']
					else:
						best_hits_dic[query_gene] = [parts[1], 0.1, None, 'Low confidence ortholog - potential false positive']
			elif len(top_hits) != 1:
				bscore1 = top_hits[0]['score']
				bscore2 = top_hits[1]['score']
				escore1 = top_hits[0]['evalue']
				escore2 = top_hits[1]['evalue']
				Br=float(bscore1/bscore2)
				Er=0
				significant_escore = False
				if escore1!=0:
					Er=(escore2//escore1)
				elif escore1==0 and escore2!=0:
					significant_escore = True
				if Br>1.5 or Er>10**20 or significant_escore == True:#code block to decide if top BLAST hit is a clear winner or if phylogenetic analysis is needed
					parts = top_hits[0]['line'].strip().split('\t')
					best_hits_dic[query_gene] = [parts[1],1, None,'High confidence ortholog']
				else:
					parts = top_hits[0]['line'].strip().split('\t')
					sublist = [parts[0]]
					anchor_list_all = []
					synteny_found = False  # Flag to track if synteny was found
					if len(final_sorted_best_hits[query_gene]) > 1: #1 is here to check synteny of only the top BLAST hit
						potential_orthologs_info = hits[:1]
					elif len(final_sorted_best_hits[query_gene]) <= 1:
						potential_orthologs_info = hits[:len(final_sorted_best_hits[query_gene])]
					# Check synteny for each potential ortholog
					for each in potential_orthologs_info:
						potential_ortholog = each['line'].strip().split('\t')[1]
						anchor_list_all.append(potential_ortholog)
						for anchor in anchor_list_all:
							anchor_list = [anchor]
							synteny = synteny_checker(anchor_list, sublist, pos_dic_query, pos_dic_ref,synteny_cutoff, flank_number, best_side, top_local_hits)
							if synteny == 'yes':
								best_hits_dic[query_gene] = [anchor, 0.5, None, 'Moderate confidence ortholog']
								synteny_found = True
								break
						if synteny_found:
							break
					# Only proceed with tree analysis if top hit was NOT syntenic
					if not synteny_found:
						if do_tree == 'yes':
							# Build test_group only once, outside the synteny loop
							if len(final_sorted_best_hits[query_gene]) > number_of_hits:
								top_hits = final_sorted_best_hits[query_gene][:number_of_hits]
							elif len(final_sorted_best_hits[query_gene]) <= number_of_hits:
								top_hits = final_sorted_best_hits[query_gene][:len(final_sorted_best_hits[query_gene])]
							# Add hits to test_group only once
							for hit in top_hits:
								parts = hit['line'].strip().split('\t')
								hit_gene = parts[1]
								test_group.append(hit_gene)
							outgroup = test_group[-1]
							top_fwd_hit = test_group[1]
							top_gene_pair = (str(query_gene) + '\t' + str(top_fwd_hit))
							out_gene_pair = (str(query_gene) + '\t' + str(outgroup))
							top_bit_score = all_ref_genes_bitscores[top_gene_pair]
							outgroup_bit_score = all_ref_genes_bitscores[out_gene_pair]
							futures.append(executor.submit(do_mafft_clean_tree, query_gene, test_group,query_sequences, ref_sequences, tmp_dir, mafft,occupancy, fasttree, outgroup, possible_orthologs,logger,top_bit_score,outgroup_bit_score))
						elif do_tree == 'no':
							candidate_orthologs = []
							if len(final_sorted_best_hits[query_gene]) > possible_orthologs:
								potential_orthologs_info = hits[:possible_orthologs]
							elif len(final_sorted_best_hits[query_gene]) <= possible_orthologs:
								potential_orthologs_info = hits[:len(final_sorted_best_hits[query_gene])]
							for each in potential_orthologs_info:
								potential_ortholog = each['line'].strip().split('\t')[1]
								candidate_orthologs.append(potential_ortholog)
							best_hits_dic[query_gene] = [parts[1], 0.1, candidate_orthologs, 'Low confidence ortholog']
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		for future in iterator:
			try:
				query, ortholog, candidate_orthologs, score_reduction_factor = future.result()
				if ortholog:
					conf_score = float(0.25 * score_reduction_factor)
					best_hits_dic[query] = [ortholog, conf_score, candidate_orthologs,'Low confidence ortholog']
				else:
					best_hits_dic[query] = [candidate_orthologs[0], 0.05, candidate_orthologs,'Phylogeny failed; Fallback to top BLAST hit']
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise

	#writing out the orthologs file
	with open(orthologs_file, "w") as out:
		out.write(f"Sample_gene\tOrtholog\tOrtholog_score\tCandidate_orthologs\tOrtholog_reliability\n")
		for query, values in best_hits_dic.items():
			ortholog   = values[0] if values[0] is not None else "--"
			score      = values[1] if values[1] is not None else "--"
			candidates = ", ".join(str(x) for x in values[2]) if values[2] is not None else "--"
			confidence = values[3] if values[3] is not None else "--"
			out.write(f"{query}\t{ortholog}\t{score}\t{candidates}\t{confidence}\n")

	#writing out the fwd best hits file
	with open(best_hits_file, "w") as out:
		for query_gene, hits in final_sorted_best_hits.items():
			for hit in hits:
				line = hit['line']
				parts = line.strip().split("\t")
				if top_local_hits[query_gene] == parts[1]:
					out.write(hit['line'])
	return best_hits_file

#functions to check for presence of working BUSCO installation in the user's system
def check_native_busco(user_path=None):
	# If user provides an explicit path, use that
	if user_path:
		if os.path.isfile(user_path) and os.access(user_path, os.X_OK):
			try:
				subprocess.run([user_path, "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
				return user_path  # Valid and executable
			except Exception:
				return None
		else:
			return None  # Invalid path

	# Else fallback to checking PATH
	busco_path = shutil.which("busco")
	if busco_path:
		try:
			subprocess.run(["busco", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			return "busco"
		except Exception:
			return None

	return None  # BUSCO not found


def check_docker_busco():
	# CHeck docker installation
	if shutil.which("docker") is None:
		return None

	# Check for presence of 'busco' in the name of the local image?
	try:
		images = subprocess.check_output(["docker", "images", "--format", "{{.Repository}}:{{.Tag}}"]).decode().splitlines()
		busco_images = [img for img in images if "busco" in img.lower()]
		if not busco_images:
			return None
	except Exception:
		return None

	# Try running 'busco --version' inside one of them
	for image in busco_images:
		try:
			subprocess.run(
				["docker", "run", "--rm", image, "busco", "--version"],
				check=True,
				stdout=subprocess.DEVNULL,
				stderr=subprocess.DEVNULL
			)
			return image  # BUSCO is usable via this Docker image
		except Exception:
			continue

	return None  # BUSCO not available via Docker

def detect_busco(busco_user_path):

	if busco_user_path == 'busco_docker':
		docker_busco = check_docker_busco()
		if docker_busco:
			return docker_busco

	else:
		#try native busco
		native = check_native_busco(busco_user_path)
		if native:
			return native

	return None

#function to collect self max bit scores from self blast for normalization in the forward blast step
def get_self_max_bit_scores(result, max_bit_score_dic_dir):
	file_name = get_basename(str(result))
	org_name = str(file_name[0]).strip().split("_self")
	best_hits = {}
	with open(result, "r") as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				best_hits[parts[0]].append({'score': float(parts[-1]), 'evalue': float(parts[-2]), 'line': line})
			except:
				best_hits.update({parts[0]: [{'score': float(parts[-1]), 'evalue': float(parts[-2]), 'line': line}]})
			line = f.readline()
	final_sorted_hits = {}
	singletons_list_initial = []
	for key in sorted(best_hits.keys()):
		data = sorted(best_hits[key], key=lambda x: (-x['score'], x['evalue'])) #sorting the hits based on descending order of bit scores; in case of bit scores tie go for second step sorting based on ascending order of e-value
		final_sorted_hits.update({key: data})
	#code block to get the maximum possible bit score of each gene i.e. the bitscore of the genes' self hit against itself
	self_max_bit_scores = {}
	gene_bit_scores = {}
	for key in sorted(final_sorted_hits.keys()):
		if (final_sorted_hits[key][0]['line'].split('\t')[0]) == (final_sorted_hits[key][0]['line'].split('\t')[1]):
			gene_bit_scores[final_sorted_hits[key][0]['line'].split('\t')[0]]=(final_sorted_hits[key][0]['line']).strip().split('\t')[-1]
		elif (final_sorted_hits[key][1]['line'].split('\t')[0]) == (final_sorted_hits[key][1]['line'].split('\t')[1]):
			gene_bit_scores[final_sorted_hits[key][1]['line'].split('\t')[0]] = (final_sorted_hits[key][1]['line']).strip().split('\t')[-1]
	self_max_bit_scores[org_name[0]]=gene_bit_scores
	pickle_dic_file = os.path.join(max_bit_score_dic_dir,str(org_name[0])+".pkl")
	with open(pickle_dic_file,'wb') as f:
		pickle.dump(self_max_bit_scores,f)

#Collecting second best hits from self blast/ diamond/ mmseqs2 result files and segregating singletons and duplicates
def sec_besthits_get_singletons(result, best_hits_file, singletons, busco_genes_list,self_similarity_cutoff,normalized_bit_score,tmp_dir,dpi_no,norm_bit_score_plots_dir):
	#removing placeholders from gene names used for busco organism_analysis
	original_busco_genes_list=[]
	reverse_replacements = {v: k for k, v in replacements.items()}
	for gene in busco_genes_list:
		original_gene = gene
		for replacement, char in reverse_replacements.items():
			original_gene=original_gene.replace(replacement,char)
		original_busco_genes_list.append(original_gene)
	messages = []
	file_name = get_basename(str(result))
	org_name = str(file_name[0]).strip().split("_self")[0]
	best_hits = {}
	with open(result, "r") as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				best_hits[parts[0]].append({'score': float(parts[-1]), 'evalue':float(parts[-2]), 'line': line})
			except:
				best_hits.update({parts[0]: [{'score': float(parts[-1]), 'evalue':float(parts[-2]), 'line': line}]})
			line = f.readline()
	final_sorted_hits = {}
	singletons_list_initial=[]
	for key in sorted(best_hits.keys()):
		data = sorted(best_hits[key], key=lambda x: (-x['score'], x['evalue'])) #sorting the hits based on descending order of bit scores; in case of bit scores tie go for second step sorting based on ascending order of e-value
		final_sorted_hits.update({key: data})
	#code to get normalized bit scores for automatic threshold determination
	all_norm_bit_scores = []
	busco_norm_bit_scores=[]
	for key in sorted(final_sorted_hits.keys()):
		if (len(final_sorted_hits[key])) > 1:
			if (final_sorted_hits[key][1]['line'].split('\t')[0]) == (final_sorted_hits[key][1]['line'].split('\t')[1]):
				individual_norm_bit_scores = (float(final_sorted_hits[key][0]['score']) / float(final_sorted_hits[key][1]['score']))
			elif (final_sorted_hits[key][1]['line'].split('\t')[0]) != (final_sorted_hits[key][1]['line'].split('\t')[1]):
				individual_norm_bit_scores = (float(final_sorted_hits[key][1]['score']) / float(final_sorted_hits[key][0]['score']))
			all_norm_bit_scores.append(individual_norm_bit_scores)
			# Check if this gene is a BUSCO single copy Complete gene
			if key in original_busco_genes_list:
				busco_norm_bit_scores.append(individual_norm_bit_scores)
	#code to plot normalized bit score distribution
	# Create figure with larger size for better visibility
	plt.figure(figsize=(12, 8))
	# Calculate optimal number of bins using Freedman-Diaconis rule
	data_array = np.array(all_norm_bit_scores)
	q75, q25 = np.percentile(data_array, [75, 25])
	iqr = q75 - q25
	bin_width = 2 * iqr / (len(data_array) ** (1 / 3))
	n_bins = int((data_array.max() - data_array.min()) / bin_width)
	n_bins = min(max(n_bins, 30), 100)  # Ensure between 30-100 bins
	# Create histogram
	counts, bins, patches = plt.hist(all_norm_bit_scores, bins=n_bins, density=True,alpha=0.6, color='skyblue', edgecolor='navy',linewidth=0.5, label='Histogram')
	#Bounded KDE using reflection (most accurate smooth curve)
	# Convert list to numpy array for arithmetic operations
	data_array = np.array(all_norm_bit_scores)
	data_min, data_max = data_array.min(), data_array.max()
	# Create reflected data
	reflected_data = np.concatenate([2 * data_min - data_array,data_array,2 * data_max - data_array ])  # Reflect above maximum
	# Create KDE and evaluate only in the data range
	kde = stats.gaussian_kde(reflected_data)
	x_range = np.linspace(data_min, data_max, 1000)
	density = kde(x_range)
	# Plot the bounded curve
	plt.plot(x_range, density, color='red', linewidth=2.5, label='Bounded KDE', alpha=0.8)
	# Formatting
	plt.xlabel('Normalized bit scores', fontsize=12)
	plt.ylabel('Density', fontsize=12)
	plt.title(f'Distribution of Normalized Bit Scores\n(n={len(all_norm_bit_scores):,}, bins={n_bins})',fontsize=14)
	plt.xlim(0, 1)  # Set limits to meaningful range
	plt.grid(True, alpha=0.3)
	plt.legend()
	# Add summary statistics as text
	stats_text = f"""Statistics:
	Min: {data_array.min():.4f}
	Max: {data_array.max():.4f}
	Mean: {data_array.mean():.4f}
	Median: {np.median(data_array):.4f}"""
	stats_text = stats_text.replace("\t", " ")
	plt.text(0.01, 0.9, stats_text, transform=plt.gca().transAxes,verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
	# Save the plot
	plot_destination = os.path.join(norm_bit_score_plots_dir,f"{org_name}.png")
	plt.savefig(plot_destination,dpi=dpi_no, bbox_inches='tight')
	plt.close()
	if normalized_bit_score == 'auto':
		if len (busco_genes_list) > 50:
			threshold = np.percentile(busco_norm_bit_scores, 95)
			method_used = 'busco'
			self_similarity_cutoff = 0
		elif len(busco_genes_list) < 50:
			threshold = 0
			method_used = 'default fall back'
	else:
		threshold = normalized_bit_score
		method_used = 'manual'
	messages.append(f"For segregating singletons and duplicates in '{org_name}', {threshold} was used as threshold value. Method used to determine the threshold was {method_used}.")
	with open(best_hits_file, "w") as out:
		for key in sorted(final_sorted_hits.keys()):
			if (len(final_sorted_hits[key]) > 1 and ((float(final_sorted_hits[key][1]['score'])/float(final_sorted_hits[key][0]['score'])) >= float(threshold))):  # Ensure there are at least 2 hits
				if (final_sorted_hits[key][1]['line'].split('\t')[0]) == (final_sorted_hits[key][1]['line'].split('\t')[1]):
					if float(final_sorted_hits[key][0]['line'].split('\t')[2])>self_similarity_cutoff:#second level of cutoff check to ensure enough similarity among gene duplicates
						out.write(final_sorted_hits[key][0]['line'])
					else:
						singletons_list_initial.append(final_sorted_hits[key][0]['line'].split('\t')[0])
				else:
					if float(final_sorted_hits[key][1]['line'].split('\t')[2]) > self_similarity_cutoff:
						out.write(final_sorted_hits[key][1]['line'])
					else:
						singletons_list_initial.append(final_sorted_hits[key][1]['line'].split('\t')[0])
			elif ((len(final_sorted_hits[key]) > 1 and ((float(final_sorted_hits[key][1]['score'])/float(final_sorted_hits[key][0]['score'])) < float(threshold)) or (len(final_sorted_hits[key]) == 1))):# if the gene's second best self hit is lower than the normalized bit score cut-off or if the gene is its own and only best hit, write out the gene line to singleton file
				singletons_list_initial.append(final_sorted_hits[key][0]['line'].split('\t')[0])
	duplicates = []
	duplicate_pair = []
	with open(best_hits_file, 'r') as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			duplicate_pair.append(parts[0])
			duplicate_pair.append(parts[1])
			duplicates.append(duplicate_pair)
			duplicate_pair = []
			line = f.readline()
	duplicates_copy = duplicates.copy()
	stratify = []
	for subgroup in duplicates_copy:
		for each in subgroup:
			stratify.append(str(each))
	stratify_set = set(stratify)
	singletons_set = set(singletons_list_initial)
	common_elements = singletons_set.intersection(stratify_set)
	new_singleton_set = singletons_set.difference(common_elements)
	singletons_list=list(new_singleton_set)
	if singletons_list:
		singletons_count = len(singletons_list)
		with open(singletons, 'w') as outs:
			for each in singletons_list:
				outs.write(str(each)+'\n')  # writes only the gene IDs in the first column into the output file
	return messages

#function to check if pickle files in a folder have empty dictionaires
def all_pickles_empty(folder_path):
	pickle_files = glob.glob(os.path.join(folder_path, '*.pkl'))

	if not pickle_files:  # No pickle files found
		return False  # or True, depending on your logic

	for f in pickle_files:
		with open(f, 'rb') as file:
			data = pickle.load(file)
			if data != {}:  # Found non-empty dictionary
				return False
	return True  # All dictionaries are empty

#master function to parallelize singleton duplicate segregation
def parallelize_singleton_duplicate_segregation(self_hits_dir,self_best_hits_dir,singleton_dir,tmp_dir,self_similarity_cutoff,normalized_bit_score,cores,logger,dpi_no,norm_bit_score_plots_dir,busco_dic_file_list,busco_dic_dir):
	# Dynamically detect number of CPU cores
	total_cores = cores
	selfblast_folder_path = Path(self_hits_dir)
	selfblast_file_list = list(selfblast_folder_path.iterdir())
	num_tasks = len(selfblast_file_list)
	max_parallel_organisms = distribute_cores_single_level(total_cores, num_tasks)
	logger.info(f" Parallelizing {max_parallel_organisms} organism-level processes")

	with ProcessPoolExecutor(max_workers=max_parallel_organisms) as executor:
		futures = []
		for each in selfblast_file_list:
			orgs_name = get_basename(str(each))[0]
			organism = orgs_name.split('_self')[0]
			blast_sec_best_hits_file = os.path.join(self_best_hits_dir, organism + "_self_sec_best_hits.tsv")
			singleton_file = os.path.join(singleton_dir, organism + "_singletons.tsv")
			if all_pickles_empty(busco_dic_dir):
				busco_genes_list = []
			else:
				org_found =False
				for busco_dic_file in busco_dic_file_list:
					busco_orgname = get_basename(str(busco_dic_file))[0]
					if busco_orgname == organism:
						with open(busco_dic_file, 'rb') as f:
							busco_dic = pickle.load(f)
							if 'single_copy_genes' in busco_dic:
								busco_genes_list = busco_dic['single_copy_genes']
							else:
								busco_genes_list = []
							org_found = True
				if org_found == False:
					logger.info(f"Single copy BUSCO gene lists of organism {organism} are not available.")
					busco_genes_list = []
			futures.append(executor.submit(sec_besthits_get_singletons, str(each),blast_sec_best_hits_file,singleton_file,busco_genes_list,self_similarity_cutoff,normalized_bit_score,tmp_dir,dpi_no,norm_bit_score_plots_dir))

		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures),desc="Organisms processed for singleton duplicate genes segregation step.")
		# Add tqdm to monitor futures
		for future in iterator:
			try:
				messages = future.result(timeout = 60)
				if messages:
					for msg in messages:
						logger.info(msg)
			except TimeoutError:
				logger.error(f"[TIMEOUT] Processing took too long.")
	logger.info(f"Completed singleton and duplicate genes segregation steps.")

#Executing blast/ diamond/ mmseqs2 and obtaining the best hits file by calling the besthits function from within this function
def align_func(no_trans_pep, database, result, best_hits_file, cpu, tmp_dir, tool, blastp, diamond, mmseqs2, eval, self, singletons, normalized_bit_score,messages,bitscore_dic,pep_ref,number_of_hits,mafft,occupancy,fasttree,logger,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,organism_name,org_db, pos_dic_query, pos_dic_ref,orthologs_file, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir):
	if tool == 'blast':
		if not os.path.isfile(result):
			cmd = blastp + ' -query ' + no_trans_pep + ' -db ' + database + ' -evalue ' + str(eval) + ' -out ' + result + " -outfmt 6 " + " -num_threads " + str(cpu)
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()
		if not self and singletons is None:
			if not os.path.isfile(orthologs_file):
				aligner_final_result = besthits(result, best_hits_file,bitscore_dic,no_trans_pep,pep_ref,number_of_hits,tmp_dir,mafft,occupancy,fasttree,cpu,logger,fwd_similarity_cutoff,score_ratio_cutoff, pos_dic_query, pos_dic_ref,orthologs_file, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree)
			return aligner_final_result
		elif self and singletons is not None:
			get_self_max_bit_scores(result, max_bit_score_dic_dir)

	elif tool == 'diamond':
		if not os.path.isfile(result):
			cmd = diamond + ' blastp --db ' + database + ' --evalue ' + str(eval) + ' --query ' + no_trans_pep + ' --out ' + result + ' --outfmt 6 --quiet --ultra-sensitive --threads ' + str(cpu)
			p = subprocess.Popen(args = cmd, shell = True)
			p.communicate()
		if not self and singletons is None:
			if not os.path.isfile(orthologs_file):
				aligner_final_result = besthits(result, best_hits_file,bitscore_dic,no_trans_pep,pep_ref,number_of_hits,tmp_dir,mafft,occupancy,fasttree,cpu,logger,fwd_similarity_cutoff,score_ratio_cutoff, pos_dic_query, pos_dic_ref,orthologs_file, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree)
			return aligner_final_result
		elif self and singletons is not None:
			get_self_max_bit_scores(result, max_bit_score_dic_dir)

	elif tool == 'mmseqs2':
		result_db = os.path.join(tmp_dir,organism_name+'_mmseqs_db')
		mmseqs_tmp = os.path.join(tmp_dir, str('MMseqs2_tmp'))
		os.makedirs(mmseqs_tmp, exist_ok=True)
		if self and singletons is not None:
			if not os.path.isfile(result):
				cmd = mmseqs2 + ' search ' + database + ' ' + database + ' ' + result_db + ' ' + mmseqs_tmp + ' --threads ' + str(cpu)
				p = subprocess.Popen(args = cmd, shell = True)
				p.communicate()
				cmd = mmseqs2 + ' convertalis ' + database + ' ' + database + ' ' + result_db + ' ' + result + " --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'"
				p = subprocess.Popen(args=cmd, shell=True)
				p.communicate()
			get_self_max_bit_scores(result, max_bit_score_dic_dir)
		elif not self and singletons is None:
			if not os.path.isfile(result):
				cmd = mmseqs2 + ' search ' + org_db + ' ' + database + ' ' + result_db + ' ' + mmseqs_tmp + ' --threads ' + str(cpu)
				p = subprocess.Popen(args=cmd, shell=True)
				p.communicate()
				cmd = mmseqs2 + ' convertalis ' + org_db + ' ' + database + ' ' + result_db + ' ' + result + " --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'"
				p = subprocess.Popen(args=cmd, shell=True)
				p.communicate()
			if not os.path.isfile(orthologs_file):
				aligner_final_result = besthits(result, best_hits_file,bitscore_dic,no_trans_pep,pep_ref,number_of_hits,tmp_dir,mafft,occupancy,fasttree,cpu,logger,fwd_similarity_cutoff,score_ratio_cutoff, pos_dic_query, pos_dic_ref,orthologs_file, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree)
			return aligner_final_result

#function to replace transcript names with parent gene names for prokaryote organisms since they lack alternate transcripts
def convert_child_to_parent_attribute (gff3_file,child_attribute,parent_attribute,child_parent_linker):
	child_parent_linker_dic = {}
	if gff3_file[-2:].lower() != 'gz':
		with open(gff3_file, "r") as f:
			gff_lines = f.readlines()
			# checking if cds feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if
						   not line.startswith('#') and len(line.split('\t')) >= 3)
	else:
		with gzip.open(gff3_file, "rt") as f:
			gff_lines = f.readlines()
			# checking if cds feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if
						   not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_mrna:
		coding_feature = 'MRNA'
		has_parent = any(
			line.split('\t')[2].upper() == 'MRNA' and str(child_parent_linker) + '=' in line.split('\t')[8]
			for line in gff_lines
			if not line.startswith('#') and len(line.split('\t')) >= 9
		)
	else:
		has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
							 if not line.startswith('#') and len(line.split('\t')) >= 3)
		if has_transcript:
			coding_feature = 'TRANSCRIPT'
			has_parent = any(
				line.split('\t')[2].upper() == 'TRANSCRIPT' and str(child_parent_linker) + '=' in line.split('\t')[8]
				for line in gff_lines
				if not line.startswith('#') and len(line.split('\t')) >= 9
			)
		else:
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)
			if has_cds:
				coding_feature = 'CDS'
				has_parent = any(
					line.split('\t')[2].upper() == 'CDS' and str(child_parent_linker) + '=' in line.split('\t')[8]
					for line in gff_lines
					if not line.startswith('#') and len(line.split('\t')) >= 9
				)
			else:
				coding_feature = 'EXON'
				has_parent = any(
					line.split('\t')[2].upper() == 'EXON' and str(child_parent_linker) + '=' in line.split('\t')[8]
					for line in gff_lines
					if not line.startswith('#') and len(line.split('\t')) >= 9
				)
	if gff3_file[-2:].lower() != 'gz':
		with open(gff3_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					parts = line.strip().split('\t')
					if len(parts) > 2:
						if has_gene and has_parent:
							if parts[2].upper() == coding_feature:
								attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
								id = attributes.get(child_attribute, "")
								child_parent_linker_id = attributes.get(child_parent_linker, "")
								child_parent_linker_dic[id]=child_parent_linker_id
				line = f.readline()
	else:
		with gzip.open(gff3_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					parts = line.strip().split('\t')
					if len(parts) > 2:
						if has_gene and has_parent:
							if parts[2].upper() == coding_feature:
								attributes = dict(attr.split("=", 1) for attr in parts[8].split(";") if "=" in attr)
								id = attributes.get(child_attribute, "")
								child_parent_linker_id = attributes.get(child_parent_linker, "")
								child_parent_linker_dic[id]=child_parent_linker_id
				line = f.readline()
	return child_parent_linker_dic

#Converting the position information in the GFF3 file to dictionary of lists with keys being the contig names and the lists having the genes ordered according to their positions
def GFF3_position_to_dic(no_trans_pep, gff3_file,logger,process_pseudos,gene_name_errors,child_parent_linker,parent_attribute):
	gene_name_errors=[]
	file_name = get_basename(str(no_trans_pep))
	final_file_name = str(file_name[0]).strip().split("_no_alt_trans")
	no_gene_no_parent = False
	has_gene = False
	has_parent = False
	if gff3_file[-2:].lower() != 'gz':
		with open(gff3_file, "r") as f:
			gff_lines = f.readlines()
			# checking if cds feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	else:
		with gzip.open(gff3_file, "rt") as f:
			gff_lines = f.readlines()
			# checking if cds feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_mrna:
		coding_feature = 'MRNA'
		has_parent = any(
			line.split('\t')[2].upper() == 'MRNA' and str(child_parent_linker)+'=' in line.split('\t')[8]
			for line in gff_lines
			if not line.startswith('#') and len(line.split('\t')) >= 9
		)
	else:
		has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
		if has_transcript:
			coding_feature = 'TRANSCRIPT'
			has_parent = any(
				line.split('\t')[2].upper() == 'TRANSCRIPT' and str(child_parent_linker)+'=' in line.split('\t')[8]
				for line in gff_lines
				if not line.startswith('#') and len(line.split('\t')) >= 9
			)
		else:
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			if has_cds:
				coding_feature = 'CDS'
				has_parent = any(
					line.split('\t')[2].upper() == 'CDS' and str(child_parent_linker)+'=' in line.split('\t')[8]
					for line in gff_lines
					if not line.startswith('#') and len(line.split('\t')) >= 9
				)
			else:
				coding_feature = 'EXON'
				has_parent = any(
					line.split('\t')[2].upper() == 'EXON' and str(child_parent_linker)+'=' in line.split('\t')[8]
					for line in gff_lines
					if not line.startswith('#') and len(line.split('\t')) >= 9
				)
	#making a list of coding genes in organism
	coding_genes=[]
	with open(no_trans_pep,'r') as f:
		line=f.readline()
		while line:
			if '>' in line:
				coding_genes.append(str(line).replace('>','').replace('\n',''))
			line=f.readline()
	if process_pseudos == 'no':
		feature = ['GENE']
	else:
		feature = ['GENE','PSEUDOGENE','MRNA','TRANSCRIPT']
	#making dictionary of lists for the organism
	pos_nos = {}
	gene_pos_per_chr = {}
	if gff3_file[-2:].lower() != 'gz':  # uncompressed gff file
		with open( gff3_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					parts = line.strip().split('\t')
					if len( parts ) > 2:
						if has_gene and has_parent:
							if parts[2].upper() in feature:
								partsnew = parts[-1].strip().split(';')
								for each in partsnew:
									if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
										if 'PSEUDO=TRUE' in parts[8].upper() or 'PSEUDOGENE' in parts[8].upper():
											pattern_ID = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
											if re.match(pattern_ID, each):
												partsnew1 = str(each).replace(str(child_parent_linker) + '=', "")
									pattern_ID = r'^' + re.escape(parent_attribute + '=') + r'.*$'
									if re.match(pattern_ID, each):
										partsnew1 = str(each).replace(str(parent_attribute) + '=', "")
								try:  # add the tuple of gene names, their corresponding start positions to the list in the dictionary of lists
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]].add((partsnew1, start_position))
									pos_nos[partsnew1] = start_position
								except KeyError:
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]] = {(partsnew1, start_position)}
									pos_nos[partsnew1] = start_position
						else:
							if parts[2].upper() == coding_feature:
								partsnew = parts[-1].strip().split(';')
								for each in partsnew:
									pattern_ID = r'^' + re.escape(parent_attribute + '=') + r'.*$'
									if re.match(pattern_ID, each):
										partsnew1 = str(each).replace(str(parent_attribute) + '=', "")
								try:  # add the tuple of gene names, their corresponding start positions to the list in the dictionary of lists
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]].add((partsnew1, start_position))
									pos_nos[partsnew1] = start_position
								except KeyError:
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]] = {(partsnew1, start_position)}
									pos_nos[partsnew1] = start_position
				line = f.readline()
		for chr, genes in gene_pos_per_chr.items():#sort the genes in each contig based on their start positions and give the list with only gene ids and not positions
			gene_pos_per_chr[chr]=[gene for gene, _ in sorted(genes, key=lambda x: x[1])]
		contig_names = gene_pos_per_chr.keys()
	else:
		with gzip.open(gff3_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					parts = line.strip().split('\t')
					if len( parts ) > 2:
						if has_gene and has_parent:
							if parts[2].upper() in feature:
								partsnew = parts[-1].strip().split(';')
								for each in partsnew:
									if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
										if 'PSEUDO=TRUE' in parts[8].upper() or 'PSEUDOGENE' in parts[8].upper():
											pattern_ID = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
											if re.match(pattern_ID, each):
												partsnew1 = str(each).replace(str(child_parent_linker) + '=', "")
									pattern_ID = r'^' + re.escape(parent_attribute + '=') + r'.*$'
									if re.match(pattern_ID, each):
										partsnew1 = str(each).replace(str(parent_attribute) + '=', "")
								try:  # add the tuple of gene names, their corresponding start positions to the list in the dictionary of lists
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]].add((partsnew1, start_position))
									pos_nos[partsnew1] = start_position
								except KeyError:
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]] = {(partsnew1, start_position)}
									pos_nos[partsnew1] = start_position
						else:
							if parts[2].upper() == coding_feature:
								partsnew = parts[-1].strip().split(';')
								for each in partsnew:
									pattern_ID = r'^' + re.escape(parent_attribute + '=') + r'.*$'
									if re.match(pattern_ID, each):
										partsnew1 = str(each).replace(str(parent_attribute) + '=', "")
								try:  # add the tuple of gene names, their corresponding start positions to the list in the dictionary of lists
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]].add((partsnew1, start_position))
									pos_nos[partsnew1] = start_position
								except KeyError:
									start_position = int(parts[3])
									gene_pos_per_chr[parts[0]] = {(partsnew1, start_position)}
									pos_nos[partsnew1] = start_position
				line = f.readline()
		for chr, genes in gene_pos_per_chr.items():#sort the genes in each contig based on their start positions and give the list with only gene ids and not positions
			gene_pos_per_chr[chr]=[gene for gene, _ in sorted(genes, key=lambda x: x[1])]
		contig_names = gene_pos_per_chr.keys()
	gene_pos_per_chr_coding_non_coding=copy.deepcopy(gene_pos_per_chr)
	gene_pos_per_chr_coding = copy.deepcopy(gene_pos_per_chr)
	#removing non-coding genes from the dictionary of lists
	# Flatten all gene names from the contig_dict values
	all_genes_in_contigs = set(gene for genes in gene_pos_per_chr.values() for gene in genes)
	# Check if all genes in gene_list are in the contig dictionary values
	all_present = all(gene in all_genes_in_contigs for gene in coding_genes)

	if all_present==True:
		for each in contig_names:
			for every in gene_pos_per_chr[each]:
				if every not in coding_genes:
					gene_pos_per_chr_coding[each].remove(every)
	elif all_present==False:
		logger.error(f"The gene and protein names are different for {final_file_name[0]}. Please check your GFF3 and PEP files.")
		gene_name_errors.append(f"The gene and protein names are different for {final_file_name[0]}.")
	return gene_pos_per_chr_coding, pos_nos, gene_pos_per_chr_coding_non_coding, gene_name_errors

#Function to clean alternate transcripts in TPM counts file
def tpmclean(gff3_input_file,output_file,child_parent_linker,child_attribute):
	no_gene_no_parent = False
	has_gene = False
	if gff3_input_file[-2:].lower() != 'gz':
		with open(gff3_input_file, "r") as f:
			gff_lines = f.readlines()
			# checking if cds feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
	else:
		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			# checking if cds feature is present in the GFF file
			has_gene = any(line.split('\t')[2].upper() == 'GENE' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
						   if not line.startswith('#') and len(line.split('\t')) >= 3)
	if has_mrna:
		coding_feature = 'MRNA'
	else:
		has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
							 if not line.startswith('#') and len(line.split('\t')) >= 3)
		if has_transcript:
			coding_feature = 'TRANSCRIPT'
		else:
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)
			if has_cds:
				coding_feature = 'CDS'
			else:
				coding_feature = 'EXON'

	nogene_noparent_counter = 0
	if gff3_input_file[-2:].lower() != 'gz':#uncompressed gff file
		transcripts_per_gene = {}
		with open(gff3_input_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != "#":
					no_gene_no_parent = False
					parts = line.strip().split('\t')
					if len(parts) > 2:
						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker)+'=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								if child_attribute in parts[-1]:
									for each in partsnew:
										pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
										if re.match(pattern_par, each):
											partsnew1 = str(each).replace(str(child_parent_linker) + '=', "")
									for every in partsnew:
										pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
										if re.match(pattern_ID, every):
											partsnew0 = str(every).replace(str(child_attribute) + '=', "")
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
					no_gene_no_parent = False
					parts = line.strip().split('\t')
					if len(parts) > 2:
						if parts[2].upper() == coding_feature:
							partsnew = parts[-1].strip().split(';')
							# Check if any attribute starts with 'Parent='
							has_parent = any(attr.startswith(str(child_parent_linker)+'=') for attr in partsnew)
							if has_gene and has_parent:
								nogene_noparent_counter += 1
								if child_attribute in parts[-1]:
									for each in partsnew:
										pattern_par = r'^' + re.escape(child_parent_linker + '=') + r'.*$'
										if re.match(pattern_par, each):
											partsnew1 = str(each).replace(str(child_parent_linker) + '=', "")
									for every in partsnew:
										pattern_ID = r'^' + re.escape(child_attribute + '=') + r'.*$'
										if re.match(pattern_ID, every):
											partsnew0 = str(every).replace(str(child_attribute) + '=', "")
							try:
								transcripts_per_gene[partsnew1].append(partsnew0)
							except KeyError:
								transcripts_per_gene.update({partsnew1: [partsnew0]})
				line = f.readline()
	gene_names = list(transcripts_per_gene.keys())
	return transcripts_per_gene, gene_names

#function for tpm file cleaning step
def tpm_cleaning_func(each,every,validated_tpm_files, exp_name, pep_file,output_tpm,child_parent_linker,child_attribute):
	messages = []
	# checking for enough datasets in the counts file for gene expression analysis
	if str(every)[-2:].lower() != 'gz':  # dealing with uncompressed exp file
		with open(str(every)) as file:
			num_columns = len(file.readline().strip().split("\t"))
	else: # compressed exp file
		with gzip.open(str(every),'rt') as file:
			num_columns = len(file.readline().strip().split("\t"))
	if num_columns <= 10:
		messages.append(str(every) + ' is not a valid counts file or does not have enough samples for reliable gene expression analysis. Skipping it...'+'\n')
		return messages
	else:
		transcripts_per_gene, gene_names = tpmclean(str(each),pep_file,child_parent_linker,child_attribute)
	# Step 1: Read the TPM file and replace transcript names with gene names
	# Create a mapping of transcript to gene for quick lookup
	transcript_to_gene = {transcript: gene for gene, transcripts in transcripts_per_gene.items() for transcript in transcripts}
	# Initialize dictionary to store gene expression values
	gene_expression = {}
	gene_name = None
	transcript_name_mismatch = 0
	if str(every)[-2:].lower() != 'gz':  # dealing with uncompressed exp file
		# Read the count table file and process the data
		with open(str(every), "r") as f:
			reader = csv.reader(f, delimiter="\t")  # Assuming tab-separated file
			header = next(reader)  # Read the first line (header)
			sample_ids = header[1:]  # Extract sample IDs (excluding 'gene' column)
			for row in reader:
				transcript = row[0]  # First column: transcript name
				expression_values = row[1:]  # Other columns: expression values
				# Find the corresponding gene name
				if transcript in gene_names:
					gene_name = transcript #when gene names are present in counts file instead of transcript names retain gene names
				else:
					if transcript in transcript_to_gene:
						gene_name = transcript_to_gene.get(transcript, transcript)  # Default to transcript name if not found
					else:
						gene_name = None
						transcript_name_mismatch += 1
				if gene_name != None:
					# Store expression values against the gene
					if gene_name in gene_expression:
						# Sum expression values if multiple transcripts map to the same gene
						gene_expression[gene_name] = [
							float(gene_expression[gene_name][i]) + float(expression_values[i])
							for i in range(len(expression_values))]
					else:
						# Store the first transcript’s values
						gene_expression[gene_name] = list(map(float, expression_values))
	else:  # dealing with compressed TPM file
		# Read the count table file and process the data
		with gzip.open(str(every), "rt") as f:
			reader = csv.reader(f, delimiter="\t")  # Assuming tab-separated file
			header = next(reader)  # Read the first line (header)
			sample_ids = header[1:]  # Extract sample IDs (excluding 'gene' column)
			for row in reader:
				transcript = row[0]  # First column: transcript name
				expression_values = row[1:]  # Other columns: expression values
				# Find the corresponding gene name
				if transcript in gene_names:
					gene_name = transcript  # when gene names are present in counts file instead of transcript names retain gene names
				else:
					if transcript in transcript_to_gene:
						gene_name = transcript_to_gene.get(transcript, transcript) # Default to transcript name if not found
					else:
						gene_name = None
						transcript_name_mismatch += 1
				if gene_name != None:
					# Store expression values against the gene
					if gene_name in gene_expression:
						# Sum expression values if multiple transcripts map to the same gene
						gene_expression[gene_name] = [
							float(gene_expression[gene_name][i]) + float(expression_values[i])
							for i in range(len(expression_values))]
					else:
						# Store the first transcript’s values
						gene_expression[gene_name] = list(map(float, expression_values))
	if transcript_name_mismatch == 0:
		# Step 2: Write the cleaned data to a new file
		# Write the processed data to the gzipped output file
		with gzip.open(output_tpm, "wt", newline="") as f_out:
			writer = csv.writer(f_out, delimiter="\t")
			# Write the header (gene + sample IDs)
			writer.writerow(["gene"] + sample_ids)
			# Write each gene and its expression values
			for gene, values in gene_expression.items():
				writer.writerow([gene] + values)
		messages.append('expression file written for ' + str(every))
	else:
		messages.append('expression file not processed for ' + str(every) + ' due to transcript and gene names mismatch. Expression analysis skipped for this organism.')
	return messages

#set of functions to parallelize the GeMoMa steps
def GeMoMa_executor (outdir_gemoma,input,inputs,each,GeMoMa,gff3_folder,orgs_name,threads_per_organism):
	destination_file_path = os.path.join(gff3_folder, str(orgs_name) + '.gff')
	if not os.path.isfile(destination_file_path):
		# --- filtering out organisms without annotation and running GeMoMa to obtain their annotations based on the reference annotation--- #
		cmd = 'java -jar ' + GeMoMa + ' CLI GeMoMaPipeline threads=' + str(threads_per_organism) + ' outdir=' + outdir_gemoma + ' GeMoMa.Score=ReAlign AnnotationFinalizer.r=COMPOSED AnnotationFinalizer.p=g AnnotationFinalizer.n=false t=' + str(each) + ' a=' + str(inputs) + ' g=' + str(input)
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
		source_file_path = os.path.join(outdir_gemoma, 'final_annotation.gff')
		shutil.move(source_file_path,destination_file_path)  # renaming the final_annotation.gff file to the specific orgname.gff and moving it into the GFF files containing folder from the TMP folder of GeMoMa outputs for further processing

#master function to parallelize GeMoMa annotation runs
def parallelize_GeMoMa(fasta_input_file,queries_to_annotate,tmp_dir,gff3_input_files,GeMoMa,gff3_folder,cores,logger):
	# Dynamically detect number of CPU cores
	total_cores = cores
	num_files = len(queries_to_annotate)

	max_parallel_organisms, threads_per_organism = distribute_cores(total_cores, num_files)

	logger.info(f" Parallelizing {max_parallel_organisms} organism-level processes")
	logger.info(f"{threads_per_organism} threads inside each organism for GeMoMa annotation")

	with ProcessPoolExecutor(max_workers=max_parallel_organisms) as executor:
		futures = []
		for each in fasta_input_file:
			orgs_name = get_basename(str(each))
			if orgs_name[0] in queries_to_annotate:
				outdir_gemoma = os.path.join(tmp_dir, str(orgs_name[0]) + '_GeMoMa')  # creating output folder for GeMoMa output folder of the particular query organism
				os.makedirs(outdir_gemoma, exist_ok=True)
				for input in fasta_input_file:
					ref = get_basename(str(input))
					if ref[0] == queries_to_annotate[orgs_name[0]]:
						for inputs in gff3_input_files:
							gff = get_basename(str(inputs))
							if gff[0] == queries_to_annotate[orgs_name[0]]:
								futures.append(executor.submit(GeMoMa_executor,outdir_gemoma,str(input),str(inputs),str(each),GeMoMa,gff3_folder,orgs_name[0],threads_per_organism))
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures),desc="Organisms processed for the GeMoMa annotation step")
		for future in iterator:
			try:
				result = future.result()
				if result:
					logger.info(result)
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise
	logger.info(f"Completed annotation of {len(fasta_input_file)} organism files.")

#functions to parallelize the steps involved in self blast and no alternate transcripts pep files
def process_organism_gff_fasta(fasta_file, gff3_file,tmp_dir,no_alt_trans_dir,no_alt_trans_cds_dir,ref_name,blast_dir,database_dir, self_hits_dir,self_best_hits_dir,singleton_dir,strict_start,strict_end,blast_pep,blast_db,tool, makeblastdb, blastp, diamond, mmseqs2,eval,normalized_bit_score, evo_analysis, file_type,logger,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_type,cds_dir,pep_dir,process_pseudos, validated_tpm_files, exp_input_file,validated_fasta_files,clean_cds_dir,gff_dic_dir,gff_config_parameters,attribute_dic_dir):
	messages=[]
	message2 = []
	message = []
	gene_name_errors = []
	gff_pos_details = {}
	self_max_bit_scores = {}
	peptide_to_cds_map={}
	if file_type=='fasta only':
		seqs = load_sequences(str(fasta_file))
	fasta_name = get_basename(str(fasta_file))
	gff_name = get_basename(str(gff3_file))
	if 'all' in gff_config_parameters:
		child_attribute = gff_config_parameters['all']['child_attribute']
		child_parent_linker = gff_config_parameters['all']['child_parent_linker']
		parent_attribute = gff_config_parameters['all']['parent_attribute']
	elif 'default' in gff_config_parameters:
		child_attribute = gff_config_parameters['default']['child_attribute']
		child_parent_linker = gff_config_parameters['default']['child_parent_linker']
		parent_attribute = gff_config_parameters['default']['parent_attribute']
	else:
		child_attribute = gff_config_parameters[str(fasta_name[0])]['child_attribute']
		child_parent_linker = gff_config_parameters[str(fasta_name[0])]['child_parent_linker']
		parent_attribute = gff_config_parameters[str(fasta_name[0])]['parent_attribute']
	if fasta_name != gff_name:
		return None
	if file_type == 'fasta only':
		if org_type == 'eukaryote':
			cds_out = os.path.join(cds_dir, fasta_name[0] + ".cds.fasta")
			pep_out = os.path.join(pep_dir, fasta_name[0] + ".pep.fasta")
		elif org_type == 'prokaryote':
			cds_out = os.path.join(cds_dir, fasta_name[0] + "_no_alt_trans" + ".cds.fasta")
			pep_out = os.path.join(pep_dir, fasta_name[0] + "_no_alt_trans" + ".pep.fasta")
	if file_type == 'cds fasta':
		if org_type == 'eukaryote':
			pep_out = os.path.join(pep_dir, fasta_name[0] + ".pep.fasta")
		elif org_type == 'prokaryote':
			pep_out = os.path.join(pep_dir, fasta_name[0] + "_no_alt_trans" + ".pep.fasta")
	if org_type == 'eukaryote':
		no_alt_trans_pep_out = os.path.join(no_alt_trans_dir, fasta_name[0] + "_no_alt_trans" + ".pep.fasta")
		no_alt_trans_cds_out = os.path.join(no_alt_trans_cds_dir, fasta_name[0] + "_no_alt_trans" + ".cds.fasta")
	elif org_type == 'prokaryote':
		no_alt_trans_pep_out = pep_out
		no_alt_trans_cds_out = os.path.join(clean_cds_dir, fasta_name[0] +".cds.fasta")#containing cds sequences that gave rise to the pep
	if file_type == 'fasta only':
		if not os.path.isfile(cds_out):
			transcript_information, message = load_transcript_information_from_gff3(str(gff3_file),process_pseudos,child_attribute,child_parent_linker)
			construct_CDS_file(transcript_information, cds_out, seqs, child_parent_linker)
		if not os.path.isfile(pep_out):
			message2, peptide_to_cds_map = transeq(cds_out, pep_out, strict_start, strict_end,peptide_to_cds_map)
		if org_type == 'eukaryote':
			if not os.path.isfile(no_alt_trans_pep_out):
				pepclean(str(gff3_file), pep_out, no_alt_trans_pep_out,no_alt_trans_cds_out,'fasta',peptide_to_cds_map,child_attribute,child_parent_linker)
		if org_type == 'prokaryote':
			attribute_dic_file = os.path.join(attribute_dic_dir, str(fasta_name[0]) + '.pkl')
			if not os.path.isfile(attribute_dic_file):
				replacement_dic = convert_child_to_parent_attribute(gff3_file, child_attribute, parent_attribute,child_parent_linker)
				with open(attribute_dic_file, "wb") as f:
					pickle.dump(replacement_dic, f)
			with open(attribute_dic_file,'rb') as f:
				replacement_dic = pickle.load(f)
			with open(no_alt_trans_cds_out,'w')as out:
				for header in peptide_to_cds_map.keys():
					out.write(">"+str(replacement_dic[header])+"\n")
					out.write(str(peptide_to_cds_map[header])+"\n")
			with fileinput.FileInput(pep_out,inplace=True) as file:
				for line in file:
					if line.startswith(">"):
						old_header = line[1:].strip()
						new_header = replacement_dic.get(old_header, old_header)  # keep same if not in map
						print(f">{new_header}")
					else:
						print(line, end="")

	if file_type == 'cds fasta':
		if not os.path.isfile(pep_out):
			message2, peptide_to_cds_map = transeq(str(fasta_file), pep_out, strict_start, strict_end,peptide_to_cds_map)
		if org_type == 'eukaryote':
			if not os.path.isfile(no_alt_trans_pep_out):
				pepclean(str(gff3_file), pep_out, no_alt_trans_pep_out, no_alt_trans_cds_out,'cds', peptide_to_cds_map,child_attribute,child_parent_linker)
		if org_type == 'prokaryote':
			attribute_dic_file = os.path.join(attribute_dic_dir, str(fasta_name[0]) + '.pkl')
			if not os.path.isfile(attribute_dic_file):
				replacement_dic = convert_child_to_parent_attribute(gff3_file, child_attribute, parent_attribute,child_parent_linker)
			with open(attribute_dic_file, "wb") as f:
				pickle.dump(replacement_dic, f)
			with open(attribute_dic_file, 'rb') as f:
				replacement_dic = pickle.load(f)
			with open(no_alt_trans_cds_out, 'w') as out:
				for header in peptide_to_cds_map.keys():
					out.write(">" + str(replacement_dic[header]) + "\n")
					out.write(str(peptide_to_cds_map[header]) + "\n")
			with fileinput.FileInput(pep_out,inplace=True) as file:
				for line in file:
					if line.startswith(">"):
						old_header = line[1:].strip()
						new_header = replacement_dic.get(old_header, old_header)  # keep same if not in map
						print(f">{new_header}")
					else:
						print(line, end="")

	if file_type == 'pep fasta':
		if org_type == 'eukaryote':
			if not os.path.isfile(no_alt_trans_pep_out):
				try:
					pepclean(str(gff3_file), str(fasta_file), no_alt_trans_pep_out, no_alt_trans_cds_out, 'pep',peptide_to_cds_map,child_attribute,child_parent_linker)
				except Exception as e:
					logger.error(f"Worker crashed with error: {e}")
					logger.exception("The error is as follows:")
					raise
		elif org_type == 'prokaryote':
			attribute_dic_file = os.path.join(attribute_dic_dir, str(fasta_name[0]) + '.pkl')
			if not os.path.isfile(attribute_dic_file):
				replacement_dic = convert_child_to_parent_attribute(gff3_file, child_attribute, parent_attribute,child_parent_linker)
			with open(attribute_dic_file, "wb") as f:
				pickle.dump(replacement_dic, f)
			with open(attribute_dic_file, 'rb') as f:
				replacement_dic = pickle.load(f)
			if fasta_file[-2].lower() != 'gz':
				pep_renamed_file = os.path.join(pep_dir, fasta_name[0] + "_no_alt_trans" + ".pep.fasta")
				shutil.move(fasta_file,pep_renamed_file)
				no_alt_trans_pep_out = pep_renamed_file
				with fileinput.FileInput(no_alt_trans_pep_out,inplace=True) as file:
					for line in file:
						if line.startswith(">"):
							old_header = line[1:].strip()
							new_header = replacement_dic.get(old_header, old_header)  # keep same if not in map
							print(f">{new_header}")
						else:
							print(line, end="")
			else:
				no_alt_trans_pep_out = os.path.join(pep_dir, fasta_name[0] + "_no_alt_trans" + ".pep.fasta")
				with gzip.open(fasta_file,'rt') as fin, open (no_alt_trans_pep_out,'w') as fout:
					for line in fin:
						if line.startswith(">"):
							header = line[1:].strip()
							new_header = replacement_dic.get(header, header)
							fout.write(f">{new_header}\n")
						else:
							fout.write(line)
				# remove the compressed fasta file
				os.remove(fasta_file)

	#doing the gff3_position_dic conversion here as first to just check if the coding gene names and the ID attribute names match at the start of the tool instead of putting it off to the mid step
	peptides_folder_path = Path(pep_dir)
	peptides_file_list = list(peptides_folder_path.iterdir())
	pickle_dic_file = os.path.join(gff_dic_dir,str(fasta_name[0])+'.pkl')
	if not os.path.isfile(pickle_dic_file):
		#validating and cleaning tpm files in the first few steps to avoid later errors
		gene_pos_per_chr_coding, pos_nos, gene_pos_per_chr_coding_non_coding, gene_name_errors = GFF3_position_to_dic(no_alt_trans_pep_out, gff3_file, logger, process_pseudos,gene_name_errors,child_parent_linker,parent_attribute)
		gff_pos_details[str(fasta_name[0])] = (gene_pos_per_chr_coding, pos_nos, gene_pos_per_chr_coding_non_coding)
		with open(pickle_dic_file, "wb") as f:
			pickle.dump(gff_pos_details, f)
	if len(exp_input_file)!=0 and org_type == 'eukaryote':
		output_tpm = os.path.join(validated_tpm_files, str(fasta_name[0])+'.txt.gz')
		if not os.path.isfile(output_tpm):
			for every in exp_input_file:
				exp_name = get_basename(every)[0]
				for pepfile in peptides_file_list:
					if org_type == 'eukaryote':
						pepname = get_basename(str(pepfile))[0]
					elif org_type == 'prokaryote':
						pepname = (get_basename(str(pepfile))[0]).strip().split('no_alt_trans')[0]
					if exp_name == fasta_name[0] and exp_name == pepname:
						messages_tpm = tpm_cleaning_func(gff3_file, every, validated_tpm_files, exp_name,str(pepfile),output_tpm,child_parent_linker,child_attribute)
						for msg in messages_tpm:
							logger.info(msg)
	combined_messages = message + message2 + messages
	if fasta_name[0]==ref_name:
		ref_database = os.path.join(database_dir, ref_name + "_db")
		if not os.path.isfile(ref_database):
			makedb(no_alt_trans_pep_out, ref_database, tmp_dir, tool, makeblastdb, blastp, diamond, mmseqs2, logger)
	return combined_messages, gene_name_errors

#master function to parallelize PEP file cleaning by removing alternate transcripts
def parallel_process_fasta_gff(fasta_input_file, gff3_input_file, tmp_dir,no_alt_trans_dir,no_alt_trans_cds_dir,ref_name,blast_dir,database_dir,self_hits_dir,self_best_hits_dir,singleton_dir,strict_start,strict_end,blast_pep,blast_db,tool, makeblastdb, blastp, diamond, mmseqs2,eval,normalized_bit_score, evo_analysis, file_type,cores,logger,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_type,cds_dir,pep_dir,process_pseudos, validated_tpm_files, exp_input_file,validated_fasta_files,clean_cds_dir,gff_dic_dir,gff_config_parameters,attribute_dic_dir):
	total_cores = cores
	num_files = len(fasta_input_file)
	results = []
	correct_exp_input_file = []
	final_correct_exp_input_file=[]
	pep_gene_name_errors = []
	if len(exp_input_file)!=0:
		num_exp_files = (len(exp_input_file))
		# checking for file extensions
		allowed_patterns = [
			'*.tsv', '*.txt', '*.txt.gz', '*.tpms.txt', '*.tpms.txt.gz'
		]
		tpm_ext_matches = []
		for every in exp_input_file:
			tpm_path = Path(every)
			tpm_ext_match = [pattern for pattern in allowed_patterns if fnmatch.fnmatch(tpm_path.name, pattern)]
			if not tpm_ext_match:  # empty list = no matches
				logger.error(f"Please check the extensions of your input expression file {every}. It does not match the allowed extensions.")
				tpm_ext_matches.append(tpm_ext_match)
			else:  # non-empty list = at least one match
				correct_exp_input_file.append(every)

		if len(correct_exp_input_file)==0:
			logger.error('Expression analysis not possible as input files do not have allowed extesnions.')
		else:
			i = 0
			for every in correct_exp_input_file:
				exp_name = get_basename(str(every))
				for each in gff3_input_file:
					gff_name = get_basename(str(each))
					if exp_name[0] == gff_name[0]:
						final_correct_exp_input_file.append(every)
			if len(final_correct_exp_input_file)==0:
				logger.error(" File name mismatch in one or more of your expression input files. Expression analysis not possible")

	max_parallel_organisms = distribute_cores_single_level(total_cores, num_files)

	logger.info(f"Parallelizing {max_parallel_organisms} organism-level processes")
	futures = []
	with ProcessPoolExecutor(max_workers=max_parallel_organisms) as executor:
		for fasta_file in fasta_input_file:
			name = get_basename(str(fasta_file))
			for gff_file in gff3_input_file:
				orgname = get_basename(str(gff_file))
				if name == orgname:
					futures.append(executor.submit(process_organism_gff_fasta, fasta_file, gff_file, tmp_dir,no_alt_trans_dir,no_alt_trans_cds_dir,ref_name,blast_dir,database_dir,self_hits_dir,self_best_hits_dir,singleton_dir,strict_start,strict_end,blast_pep,blast_db,tool, makeblastdb, blastp, diamond, mmseqs2,eval,normalized_bit_score, evo_analysis, file_type,logger,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_type,cds_dir,pep_dir,process_pseudos, validated_tpm_files, final_correct_exp_input_file,validated_fasta_files,clean_cds_dir,gff_dic_dir,gff_config_parameters,attribute_dic_dir))

		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures),desc="Organisms processed for PEP processing step")
		# Add tqdm to monitor futures
		for future in iterator:
			try:
				combined_messages, gene_name_errors = future.result(timeout = 60)
				if combined_messages:
					for msg in combined_messages:
						logger.info(msg)
				if gene_name_errors:
					pep_gene_name_errors.append(gene_name_errors)
			except TimeoutError:
				logger.error(f"[TIMEOUT] Processing took too long.")
			#except Exception as e:
				#logger.error(f"Worker crashed with error: {e}")
				#logger.exception("The error is as follows:")
				#raise
	logger.info(f"Obtained PEP files without alternate transcripts.")
	return pep_gene_name_errors

#function to carry out self local alignment runs
def process_organism_self_blast(no_alt_trans_pep_out, database, blast_output_file,blast_sec_best_hits_file, blast_threads, tmp_dir, tool, blastp, diamond,mmseqs2, eval, self, singletons, normalized_bit_score, messages, dummy_dic,pep_ref, number_of_hits, mafft, occupancy, fasttree, logger,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff,pep_name, pos_dic_query, pos_dic_ref, max_bit_score_dic_dir):
	align_func(no_alt_trans_pep_out, database, blast_output_file,blast_sec_best_hits_file, blast_threads, tmp_dir, tool, blastp, diamond,mmseqs2, eval, self, singletons, normalized_bit_score, messages, dummy_dic,pep_ref, number_of_hits, mafft, occupancy, fasttree, logger,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff,pep_name, None, pos_dic_query, pos_dic_ref, None, None, None, None,None, None, max_bit_score_dic_dir)

#function to parallelize self local alignment runs
def parallel_process_self_blast(ref_name, cores,peptides_file_list,database_dir, tmp_dir, tool, makeblastdb, blastp, diamond, mmseqs2, logger,self_hits_dir,self_best_hits_dir,singleton_dir,eval,normalized_bit_score, max_bit_score_dic_dir,mafft,occupancy,fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff):
	results = []
	total_cores = cores
	if ref_name != 'NA':
		num_files = (len(peptides_file_list)-1)
	else:
		num_files = (len(peptides_file_list))
	max_parallel_organisms, blast_threads_per_organism = distribute_cores(total_cores, num_files)
	logger.info(f"Parallelizing {max_parallel_organisms} organism-level processes")
	logger.info(f"{blast_threads_per_organism} threads inside each organism for self alignment")
	futures = []
	with ProcessPoolExecutor(max_workers=max_parallel_organisms) as executor:
		for file in peptides_file_list:
			pep_name = get_basename(str(file))[0].strip().split("_no_alt_trans")[0]
			if pep_name != ref_name:
				database = os.path.join(database_dir, pep_name + "_db")
				if not os.path.isfile(database):
					makedb(str(file), database, tmp_dir, tool, makeblastdb, blastp, diamond, mmseqs2, logger)
					messages = []
					self = True
					dummy_dic = {}
					pep_ref = None
					number_of_hits = 0
					pos_dic_query = {}
					pos_dic_ref = {}
					blast_output_file = os.path.join(self_hits_dir, pep_name + "_self.tsv")
					blast_sec_best_hits_file = os.path.join(self_best_hits_dir,pep_name + "_self_sec_best_hits.tsv")
					singletons = os.path.join(singleton_dir, pep_name + "_singletons.tsv")
					files_to_check = [blast_output_file,blast_sec_best_hits_file,singletons]
					all_files_exist = all(os.path.isfile(f) for f in files_to_check)
					if all_files_exist:
						pass
					else:
						futures.append(executor.submit(process_organism_self_blast,str(file), database, blast_output_file,blast_sec_best_hits_file, blast_threads_per_organism, tmp_dir, tool, blastp, diamond,mmseqs2, eval, self, singletons, normalized_bit_score, messages, dummy_dic,pep_ref, number_of_hits, mafft, occupancy, fasttree, logger,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff,pep_name, pos_dic_query, pos_dic_ref, max_bit_score_dic_dir))
		# Wrap the iterator conditionally
		iterator=as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures),desc="Organisms processed for the self BLAST step")
		#Add tqdm to monitor futures
		for future in iterator:
			try:
				result=future.result()
				results.append(result)
			except Exception as e:
				logger.error(f"Error in future: {e}")
	logger.info(f"Completed self alignment of query organism(s).")

#function to predownload BUSCO datasets
def pre_download_databases(busco_path, org_type, busco_dir, busco_db_dir, logger):
	"""Pre-download BUSCO databases to avoid concurrent download issues"""
	logger.info("Pre-downloading BUSCO databases...")

	# Ensure the database directory exists
	os.makedirs(busco_db_dir, exist_ok=True)

	dummy_fasta = os.path.join(busco_dir, "dummy.faa")
	with open(dummy_fasta, 'w') as f:
		f.write(">dummy\nMVKIILFVGLLFSSVTYGC\n")

	if org_type == 'eukaryote':
		cmd = f"{busco_path} -i {dummy_fasta} -m proteins --auto-lineage-euk -q -o dummy_download --out_path {busco_dir} --download_path {busco_db_dir}"
	elif org_type == 'prokaryote':
		cmd = f"{busco_path} -i {dummy_fasta} -m proteins --auto-lineage-prok -q -o dummy_download --out_path {busco_dir} --download_path {busco_db_dir}"

	logger.info(f"Downloading BUSCO databases to: {busco_db_dir}")
	# Redirect stderr to devnull to suppress busco error messages for the above dummy dataset
	with open(os.devnull, 'w') as devnull:
		p = subprocess.Popen(args=cmd, shell=True, stderr=devnull, stdout=devnull)
		p.communicate()

	# Cleanup
	os.remove(dummy_fasta)
	dummy_dir = os.path.join(busco_dir, "dummy_download")
	if os.path.exists(dummy_dir):
		shutil.rmtree(dummy_dir)

	logger.info("Database pre-download completed.")

#function to clean gene names for BUSCO analysis
def clean_headers_for_busco(gene_id):
	for char, replacement in replacements.items():
		safe_gene_id = gene_id.replace(char, replacement)
	return safe_gene_id

#function to run BUSCO
def run_busco(orgname, fasta_file, busco_path, busco_dir, busco_dir_final, busco_threads_per_organism, org_type,host_cache_dir, busco_db_dir, busco_pickle_file,logger,buscolineage):
	#code to extract lineages to use from busco lineage config file
	buscolineage_dic={}
	if buscolineage!='auto':
		with open (buscolineage,'r')as f:
			line=f.readline()
			while line:
				parts=line.strip().split()
				buscolineage_dic[parts[0]]=parts[1].strip()
				line=f.readline()
	#code to create a temp protein fasta file to clean headers and use this clean protein file for busco analysis - this temp file will be deleted after the busco run completes
	tmp_fasta = os.path.join(busco_dir, f"{orgname}_temp.pep.fasta")
	try:
		# Process and write headers + sequences
		with open(fasta_file, "r") as infile, open(tmp_fasta, "w") as outfile:
			for line in infile:
				if line.startswith(">"):
					new_header = clean_headers_for_busco(line.strip()[1:])  # remove ">"
					outfile.write(f">{new_header}\n")
				else:
					outfile.write(line)
		if org_type == 'eukaryote':
			if buscolineage=='auto':
				cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --auto-lineage-euk -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
			else:
				if orgname in buscolineage_dic.keys():
					lineage=buscolineage_dic[orgname]
					cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --l '+ lineage + ' -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
				else:
					logger.warning(f"BUSCO lineage not specified for {orgname}. Reverting back to auto lineage detection.")
					cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --auto-lineage-euk -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
		elif org_type == 'prokaryote':
			if buscolineage=='auto':
				cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --auto-lineage-prok -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
			else:
				if orgname in buscolineage_dic.keys():
					lineage=buscolineage_dic[orgname]
					cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --l ' + lineage + ' -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir
				else:
					logger.warning(f"BUSCO lineage not specified for {orgname}. Reverting back to auto lineage detection.")
					cmd = busco_path + ' -i ' + tmp_fasta + ' -m proteins --auto-lineage-prok -q -o ' + orgname + ' -c ' + busco_threads_per_organism + ' --out_path ' + busco_dir + ' --download_path ' + busco_db_dir

		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
	finally:
		#remove the temporary pep fasta file
		if os.path.exists(tmp_fasta):
			os.remove(tmp_fasta)
	species_dir = Path(busco_dir) / orgname
	# Get all run_* folders (ignore auto_lineage and subdirs of auto_lineage)
	run_dirs = [
		p for p in species_dir.iterdir()
		if p.is_dir() and p.name.startswith("run_") and "auto_lineage" not in p.name
	]

	# Decide which run folder to use
	if not run_dirs:
		logger.error(f"BUSCO error for {orgname}")
		feff_result = str(orgname) + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA'
		busco_dic = {}
		with open(busco_pickle_file, 'wb') as f:
			pickle.dump(busco_dic, f)
		selected_dir = ''

	elif len(run_dirs) == 1:
		selected_dir = run_dirs[0]

	else:
		# Prefer a folder that does NOT include "eukaryota"
		if org_type == 'eukaryote':
			filtered = [d for d in run_dirs if "eukaryota" not in d.name.lower()]
			if not filtered:
				logger.warning(f"Multiple run_* folders found for {orgname}, but all are eukaryota-based.")
				selected_dir = run_dirs[0] # use the first eukaryota run folder
			else:
				selected_dir = filtered[0]  # Use the first non-eukaryota run folder
		elif org_type == 'prokaryote':
			filtered = [d for d in run_dirs if "prokaryota" not in d.name.lower()]
			if not filtered:
				logger.warning(f"Multiple run_* folders found for {orgname}, but all are prokaryota-based..")
				selected_dir = run_dirs[0] # use the first prokaryota run folder
			else:
				selected_dir = filtered[0]  # Use the first non-prokaryota run folder

	#Find full_table.tsv inside selected folder (ignore subfolders of auto_lineage if any)
	if selected_dir:
		full_tables = list(selected_dir.rglob("full_table.tsv"))
	else:
		full_tables = []

	if not full_tables:
		logger.error(f"No 'full_table.tsv' found for {orgname}")
		if not os.path.isfile(busco_pickle_file):
			busco_dic = {}
			with open(busco_pickle_file, 'wb') as f:
				pickle.dump(busco_dic, f)

	else:
		# Get the path to the full_table.tsv
		out_path = full_tables[0]
		full_table = str(out_path)
		complete_genes = []
		all_entries = []
		busco_list = []
		with open(full_table, 'r') as f:
			for line in f:
				if not line.startswith('#'):
					columns = line.strip().split('\t')
					if columns and  len(columns) >= 3:  # Ensure the line is not empty
						busco_list.append(columns[0])  # First column contains the BUSCO ID
						busco_id = columns[0]
						status = columns[1]
						gene_id = columns[2]

						all_entries.append({
							'busco_id': busco_id,
							'status': status,
							'gene_id': gene_id
						})
		# Count frequencies
		busco_id_counts = Counter([entry['busco_id'] for entry in all_entries])
		status_counts = Counter([entry['status'] for entry in all_entries])

		# Extract single-copy genes (Strategy 2)
		complete_from_single = []
		complete_from_multi = []

		for entry in all_entries:
			if entry['status'] == 'Complete' and entry['gene_id']:
				if busco_id_counts[entry['busco_id']] == 1:
					complete_from_single.append(entry['gene_id'])
				else:
					complete_from_multi.append(entry['gene_id'])

		single_copy_busco_genes = complete_from_single

		gene_counts = Counter(busco_list)  # Count occurrences of each gene
		frequency_counts = Counter(gene_counts.values())  # Count occurrences of each frequency

		x_values = sorted(frequency_counts.keys())  # Unique occurrence frequencies
		y_values = [frequency_counts[x] for x in x_values]  # Count of genes occurring those many times

		#code to extract C, D percentages from busco log file
		busco_log_path = next(Path(species_dir).rglob("busco.log"))
		busco_log = str(busco_log_path)
		with open (busco_log,'r') as f:
			contents = f.read()
			# Find all matches of the C/S/D/F/M pattern to take the C, D values of the specific BUSCO lineage instead of the general one that is presented first
			matches = re.findall(r"C:(\d+(?:\.\d+)?)%\[S:\d+(?:\.\d+)?%,D:(\d+(?:\.\d+)?)%", contents)
		if not matches:
			logger.error(f"No C/D values found in BUSCO log. BUSCO run might be incomplete for {orgname}")
			if not os.path.isfile(busco_pickle_file):
				busco_dic = {}
				with open(busco_pickle_file, 'wb') as f:
					pickle.dump(busco_dic, f)
		else:
			# If multiple matches, use the last one (usually the lineage-specific one)
			C_value = float(matches[-1][0])
			D_value = float(matches[-1][1])
			plt.bar(x_values, y_values, color='blue', alpha=0.7)
			plt.xlabel("Occurrence Frequency")
			plt.ylabel("Number of BUSCO Genes")
			plt.title("Gene Occurrence Frequency Distribution")
			plt.xticks(x_values)  # Ensure all x-axis values are labeled
			plt.grid(axis='y', linestyle='--', alpha=0.7)
			plot_file = os.path.join(busco_dir_final, orgname+'_gene_frequency_distribution.png')
			plt.savefig(plot_file)
			plt.close()
			total_genes = sum(frequency_counts.values())  # Total number of BUSCO genes
			Feff_numerator = sum(freq * count for freq, count in frequency_counts.items())  # Σ (C_i * N_i)
			Feff = Feff_numerator / total_genes if total_genes > 0 else 1  # Avoid division by zero
			feff_result = str(orgname)+'\t'+str(Feff)+'\t'+str(C_value)+'\t'+str(D_value)
			if not os.path.isfile(busco_pickle_file):
				busco_dic={}
				busco_dic['feff'] = feff_result
				busco_dic['single_copy_genes'] = single_copy_busco_genes
				with open (busco_pickle_file,'wb')as f:
					pickle.dump(busco_dic,f)

#functions to parallelize the steps involved in forward blast files
#set of functions to parallelize the fwd blast steps
def process_organisms_fwd_blast(blast_threads_per_organism, every,blast_pep_org_name,fwd_blast_db,ref_name, tmp_dir, blast_dir, fwd_best_hits_dir, tool, blastp, diamond, mmseqs2, eval, self, singletons, normalized_bit_score,bitscore_dic,pep_ref,number_of_hits,mafft,occupancy,fasttree,logger,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_db, orthologs_dir, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir, gff_dic_dir,blast_output_file_fwd, blast_best_hits_file, orthologs_file):
	# --- Forward blast of queries vs reference ---#
	gff_posdic_folder_path = Path(gff_dic_dir)
	gff_posdic_file_list = list(gff_posdic_folder_path.iterdir())
	messages = []
	for gffdic_file in gff_posdic_file_list:
		gffdic_org_name = get_basename(str(gffdic_file))[0]
		with open(gffdic_file, 'rb') as f:
			gff_pos_details = pickle.load(f)
		if gffdic_org_name == blast_pep_org_name and gffdic_org_name != ref_name:
			pos_dic_query, pos_nos_query, coding_non_coding_query = gff_pos_details[gffdic_org_name]
		elif gffdic_org_name == ref_name:
			pos_dic_ref, pos_nos_ref, coding_non_coding_ref = gff_pos_details[gffdic_org_name]
	self=False
	singletons=None
	align_func(str(every), fwd_blast_db, blast_output_file_fwd, blast_best_hits_file, blast_threads_per_organism, tmp_dir, tool, blastp, diamond, mmseqs2, eval, self, singletons, normalized_bit_score,messages,bitscore_dic,pep_ref,number_of_hits,mafft,occupancy,fasttree,logger,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,blast_pep_org_name,org_db, pos_dic_query, pos_dic_ref,orthologs_file, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir)

#master function to parallelize forward local alignment runs
def parallel_process_organism_fwd_blast(fwd_blast_db,peptides_file_list,ref_name, tmp_dir, blast_dir, fwd_best_hits_dir, tool, blastp, diamond, mmseqs2, eval, self_dummy, singletons_dummy, normalized_bit_score,cores,logger,pep_ref,number_of_hits,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,db_names,database_folder_path,orthologs_dir, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir, gff_dic_dir):
	# Dynamically detect number of CPU cores
	total_cores = cores
	num_files = len(peptides_file_list)-1

	max_parallel_organisms, blast_threads_per_organism = distribute_cores(total_cores, num_files)

	logger.info(f"Parallelizing {max_parallel_organisms} organism-level processes")
	logger.info(f"{blast_threads_per_organism} BLAST threads inside each organism for forward blast")

	bitscores_folder_path = Path(max_bit_score_dic_dir)
	bitscores_file_list = list(bitscores_folder_path.iterdir())

	futures = []
	with ProcessPoolExecutor(max_workers=max_parallel_organisms) as executor:
		for every in peptides_file_list:
			blast_pep_org_name = (os.path.basename(str(every))).split('_no_alt_trans')[0]
			if blast_pep_org_name != ref_name:
				for bitscore_file in bitscores_file_list:
					bitscore_org_name = get_basename(str(bitscore_file))[0]
					if bitscore_org_name == blast_pep_org_name:
						with open(bitscore_file, "rb") as f:
							bitscore_full_dic = pickle.load(f)
						bitscore_dic = bitscore_full_dic[bitscore_org_name]
				for each in db_names:
					org_db_name = each.replace('_db', '')
					if blast_pep_org_name == org_db_name:
						org_db = os.path.join(database_folder_path, each)
				blast_output_file_fwd = os.path.join(blast_dir, blast_pep_org_name + "_fwd.tsv")
				blast_best_hits_file = os.path.join(fwd_best_hits_dir, blast_pep_org_name + "_best_hits.tsv")
				orthologs_file = os.path.join(orthologs_dir, blast_pep_org_name + "_orthologs.tsv")
				files_to_check = [blast_output_file_fwd, blast_best_hits_file, orthologs_file]
				all_files_exist = all(os.path.isfile(f) for f in files_to_check)
				if all_files_exist:
					pass
				else:
					futures.append(executor.submit(process_organisms_fwd_blast,blast_threads_per_organism, str(every), blast_pep_org_name,fwd_blast_db,ref_name, tmp_dir, blast_dir, fwd_best_hits_dir, tool, blastp, diamond, mmseqs2, eval, self_dummy, singletons_dummy, normalized_bit_score,bitscore_dic,pep_ref,number_of_hits,mafft,occupancy,fasttree,logger,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_db, orthologs_dir, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir, gff_dic_dir,blast_output_file_fwd, blast_best_hits_file, orthologs_file))
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures),desc="Organisms processed for the forward BLAST step")
		# Add tqdm to monitor futures
		for future in iterator:
			try:
				result = future.result()
				if result:
					logger.info(result)
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise
	logger.info("Forward local alignment completed.")

#function for appending pairs of genes in self blast file to a list for cluster confidence scoring
def self_to_set(file):
	self_pairs = set()
	with open(file, 'r')as f:
		line=f.readline()
		while line:
			parts = line.strip().split('\t')
			pair = (parts[0],parts[1])
			self_pairs.add(pair)
			line=f.readline()
	return self_pairs

#function to convert a self local alignment file to dictioanry for proximal-tandem reclassification step
def self_to_dic_for_reclassification(file):
	gene_pairs = {}
	max_scores = defaultdict(float)
	org = get_basename(file)
	all_genes = set()
	#initilaizing an empty list to append the unclassified genes i.e. genes that do not show self hit in blast, but show self hit against other genes in self blast - such genes are put in the unclassified category
	unclassified = []
	# First pass: Determine the highest bit score for each gene (self-hit)
	with open(file, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			gene1, gene2, bit_score = row[0], row[1], float(row[-1])  # Last column is bit score
			all_genes.update([gene1, gene2])
			if gene1 == gene2:  # Self-hit
				max_scores[gene1] = max(max_scores[gene1], bit_score)

	# to catch genes that do not show self hit in blast, but show self hit against other genes in self blast
	unclassified = [gene for gene in all_genes if gene not in max_scores]
	# Second pass: Normalize bit scores and store in dictionary
	with open(file, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		for row in reader:
			gene1, gene2, bit_score = row[0], row[1], float(row[-1])  # Last column is bit score
			if max_scores[gene1]!=0:
				norm_score = bit_score / max_scores[gene1]
			else:
				norm_score = 0
			# Apply the 50% threshold condition
			if norm_score >= 0.5:
				gene_pairs[f"{gene1}\t{gene2}"] = norm_score
			else:
				pass
	return gene_pairs,unclassified

#function for converting the genes in the first two columns of the self best hits into lists within a list called duplicates
def sec_besthits_to_duplicates_list(best_hits_file):
	duplicates=[]
	duplicate_pair=[]
	with open(best_hits_file,'r') as f:
		line=f.readline()
		while line:
			parts=line.strip().split('\t')
			duplicate_pair.append(parts[0])
			duplicate_pair.append(parts[1])
			duplicates.append(duplicate_pair)
			duplicate_pair = []
			line = f.readline()
	#code to remove reciprocal hits in the gene pairs in the masters duplicates list to reduce list comprehension time in the further steps
	# Use set of frozensets to remove duplicates regardless of order
	unique_pairs = {frozenset(pair) for pair in duplicates}
	# Convert frozensets back to sorted lists (optional)
	cleaned_duplicates = [list(pair) for pair in unique_pairs]

	duplicates_copy = copy.deepcopy(cleaned_duplicates)
	stratify = []
	for subgroup in duplicates_copy:
		for each in subgroup:
			stratify.append(str(each))
	stratify_set = set(stratify)
	return cleaned_duplicates

#function to add those genes that didn't show up in the self blast hits file into the temporary singletons list
def get_singletons_from_self_blast_pep (pep_file,self_file):
	set_self = set()
	with open(self_file, 'r') as file1:
		for line in file1:
			if line.strip():  # Skip empty lines
				gene_name = line.strip().split('\t')[0]
				set_self.add(gene_name)
	set_pep = set()
	with open(pep_file, 'r') as file2:
		for line in file2:
			if line.startswith('>'):
				header = line[1:].strip().split()[0]  # Extract header up to the first whitespace
				set_pep.add(header)
	remaining_singletons_set = set_pep - set_self
	remaining_singletons = list(remaining_singletons_set)
	return remaining_singletons

#grouping potential tandem gene arrays using the duplicates list and the dictionary of lists with the gene position information using the depth first search (dfs) graph-based approach
def tandems_grouping(duplicates, gene_pos_per_chr):
	# Create an adjacency graph from the duplicates list
	graph = {}
	for pair in duplicates:
		gene1, gene2 = pair
		if gene1 not in graph:
			graph[gene1] = []
		if gene2 not in graph:
			graph[gene2] = []
		graph[gene1].append(gene2)
		graph[gene2].append(gene1)

	def dfs(gene, visited, current_group):
		visited.add(gene)
		current_group.append(gene)
		for neighbor in graph.get(gene, []):
			if neighbor not in visited:
				dfs(neighbor, visited, current_group)

	tandem_arrays = []
	tandems_potential = []

	for key, genes_list in gene_pos_per_chr.items():
		visited = set()

		for gene in genes_list:
			if gene not in visited:
				current_group = []
				dfs(gene, visited, current_group)

				# Filter to keep only genes present in the current genes_list
				current_group = [g for g in current_group if g in genes_list]

				# Sort the group by their positions in the gene list
				current_group.sort(key=lambda g: genes_list.index(g))

				# Check if all genes are consecutive
				# Find all consecutive subgroups inside current_group
				indices = [genes_list.index(g) for g in current_group]
				subgroup = [current_group[0]]

				for i in range(1, len(current_group)):
					if indices[i] - indices[i - 1] == 1:
						subgroup.append(current_group[i])
					else:
						if len(subgroup) > 1:
							tandem_arrays.append(subgroup)
						subgroup = [current_group[i]]

				# Add the last subgroup
				if len(subgroup) > 1:
					tandem_arrays.append(subgroup)

	#code block to remove tandem gene duplicates from the duplicates list
	tandem_arrays_copy = copy.deepcopy(tandem_arrays)
	duplicates_copy = copy.deepcopy(duplicates)
	flattened_tandem_array = [item for sublist in tandem_arrays_copy for item in sublist]
	tandem_set = set(flattened_tandem_array)

	# Separate into kept and rejected list of lists
	no_tandems_duplicates = []
	rejected_duplicates_from_tand = []

	# Reject pairs if both the genes of a pair are in the tandem_set; Keep the other pairs in the master duplicates list
	for pair in duplicates:
		if pair[0] in tandem_set and pair[1] in tandem_set:
			rejected_duplicates_from_tand.append(pair)
		else:
			no_tandems_duplicates.append(pair)

	no_tand_dups=no_tandems_duplicates.copy()
	# Create sets for unique gene counts
	rej_set = set(
		element.strip()
		for pair in rejected_duplicates_from_tand
		for item in pair
		for element in item.split(','))

	notand_set = set(
		element.strip()
		for pair in no_tand_dups
		for item in pair
		for element in item.split(','))

	return tandem_arrays, no_tandems_duplicates, tandem_set

#grouping potential proximal gene arrays using the duplicates list and the dictionary of lists with the gene position information using the depth first search (dfs) graph-based approach
def proximal_grouping(no_tandems_duplicates, gene_pos_per_chr, proximity):
	# Create adjacency graph from the duplicates list
	graph = {}
	for gene1, gene2 in no_tandems_duplicates:
		graph.setdefault(gene1, []).append(gene2)
		graph.setdefault(gene2, []).append(gene1)

	def dfs(gene, visited, current_group):
		visited.add(gene)
		current_group.append(gene)
		for neighbor in graph.get(gene, []):
			if neighbor not in visited:
				dfs(neighbor, visited, current_group)

	proximal_arrays = []

	for key, genes_list in gene_pos_per_chr.items():
		visited = set()

		for gene in genes_list:
			if gene not in visited:
				current_group = []
				dfs(gene, visited, current_group)

				# Filter only genes in the current chromosome list
				current_group = [g for g in current_group if g in genes_list]

				# Sort the group by their positions in the gene list
				current_group.sort(key=lambda g: genes_list.index(g))

				indices = [genes_list.index(g) for g in current_group]
				gene_pos = sorted(zip(current_group, indices), key=lambda x: x[1])

				temp_group = [gene_pos[0][0]]
				prev_index = gene_pos[0][1]

				for i in range(1, len(gene_pos)):
					gene, index = gene_pos[i]
					gap = index - prev_index

					if 1 < gap <= (proximity+1):
						temp_group.append(gene)
					else:
						if len(temp_group) >= 2:
							proximal_arrays.append(temp_group)
						temp_group = [gene]

					prev_index = index

				# Final check
				if len(temp_group) >= 2:
					proximal_arrays.append(temp_group)

	#code block to remove proximal gene duplicates from the duplicates list
	proximal_arrays_copy = copy.deepcopy(proximal_arrays)
	no_tandems_duplicates_copy = copy.deepcopy(no_tandems_duplicates)
	flattened_proximal_arrays = [item for sublist in proximal_arrays for item in sublist]
	proximal_set = set(flattened_proximal_arrays)

	# Separate into kept and rejected list of lists
	no_tandems_proximals_duplicates = []
	rejected_duplicates_from_prox = []

	# Reject pairs if both the genes of a pair are in the proximal_set; Keep the other pairs in the master no tandems duplicates list
	for pair in no_tandems_duplicates:
		if pair[0] in proximal_set and pair[1] in proximal_set:
			rejected_duplicates_from_prox.append(pair)
		else:
			no_tandems_proximals_duplicates.append(pair)
	# Create sets for unique gene counts
	rej_set = set(
		element.strip()
		for pair in rejected_duplicates_from_prox
		for item in pair
		for element in item.split(','))

	noprox_set = set(
		element.strip()
		for pair in no_tandems_proximals_duplicates
		for item in pair
		for element in item.split(','))

	return proximal_arrays,no_tandems_proximals_duplicates, proximal_set

#function for reclassification of proximals and tandems
def proximal_tandem_reclassification(tandems, proximals, no_tandems_proximals_mixed_duplicates_list, pos_dic, self_dic, duplicates_master_list):
	tandems_final = []
	new_tandems = []
	tandems_copy = copy.deepcopy(tandems)  # Ensure no shared references
	proximals_copy = copy.deepcopy(proximals)  # Ensure nested lists are not shared
	tandems_copy_prox_check = copy.deepcopy(tandems)
	tandems_copy_prox_check_set=set()
	for every in tandems_copy_prox_check:
		for each in every:
			tandems_copy_prox_check_set.add(str(each))

	#collecting all the genes enclosed within a proximal array group for potential reclassification
	for subgroup in proximals:
		# Create a set of all genes in the dictionary for quick look-up
		all_genes = [gene for genes in pos_dic.values() for gene in genes]
		# Create a dictionary to map each gene to its index in the contig order
		gene_order = {gene: i for i, gene in enumerate(all_genes)}
		# Sort genes_list based on their position in gene_order
		subgroup.sort(key=lambda gene: gene_order.get(gene, float('inf')))
		subgroup_copy = copy.deepcopy(subgroup)  # Use deepcopy to isolate changes
		start = subgroup[0]
		end = subgroup[-1]
		for genes in pos_dic.values():
			if start in genes and end in genes:
				start_index = genes.index(start)
				end_index = genes.index(end)
				if start_index > end_index:
					start_index, end_index = end_index, start_index  # Ensure correct order
				focus_array = genes[start_index:end_index + 1]
				focus_array_copy = focus_array.copy()
				for each in focus_array:
					if each in subgroup:
						focus_array_copy.remove(each)
				for every in focus_array_copy:
					for each in subgroup:
						test_pair1 = str(every) + '\t' + str(each)
						test_pair2 = str(each) + '\t' + str(every)
						if test_pair1 in self_dic and test_pair2 in self_dic:
							test_pair = test_pair1  # Prioritize test_pair1
						elif test_pair1 in self_dic and test_pair2 not in self_dic:
							test_pair = test_pair1
						elif test_pair2 in self_dic and test_pair1 not in self_dic:
							test_pair = test_pair2
						else:
							test_pair = None
						if test_pair and test_pair in self_dic and self_dic[test_pair] != 0:#short circuiting based complex check to sequentially evaluate conditions
							subgroup_copy.append(str(every))
							break
				# Create a set of all genes in the dictionary for quick look-up
				all_genes = [gene for genes in pos_dic.values() for gene in genes]
				# Create a dictionary to map each gene to its index in the contig order
				gene_order = {gene: i for i, gene in enumerate(all_genes)}
				# Sort genes_list based on their position in gene_order
				subgroup_copy.sort(key=lambda gene: gene_order.get(gene, float('inf')))

				if subgroup_copy == subgroup:
					continue  # No change in the group

				gene_positions = {}
				for contig_genes in pos_dic.values():
					for i, gene in enumerate(contig_genes):
						gene_positions[gene] = i

				adjacent_groups = []
				current_group = []

				for i in range(len(subgroup_copy) - 1):
					gene1 = subgroup_copy[i]
					gene2 = subgroup_copy[i + 1]
					if abs(gene_positions[gene1] - gene_positions[gene2]) == 1:
						if not current_group:
							current_group.append(gene1)
						current_group.append(gene2)
					else:
						if current_group:
							adjacent_groups.append(current_group)
							current_group = []

				if current_group:
					adjacent_groups.append(current_group)
				# Create a set of all genes in the dictionary for quick look-up
				all_genes = [gene for genes in pos_dic.values() for gene in genes]
				# Create a dictionary to map each gene to its index in the contig order
				gene_order = {gene: i for i, gene in enumerate(all_genes)}
				adjacent = [adj for group in adjacent_groups for adj in group]
				adjacent.sort(key=lambda gene: gene_order.get(gene, float('inf')))
				if adjacent_groups:
					new_tandems.extend(adjacent_groups)
					subgroup_copy_no_adjacent = [g for g in subgroup_copy if g not in adjacent]
					if len(subgroup_copy_no_adjacent)>1:
						# Create a set of all genes in the dictionary for quick look-up
						all_genes = [gene for genes in pos_dic.values() for gene in genes]
						# Create a dictionary to map each gene to its index in the contig order
						gene_order = {gene: i for i, gene in enumerate(all_genes)}
						# Sort genes_list based on their position in gene_order
						subgroup_copy_no_adjacent.sort(key=lambda gene: gene_order.get(gene, float('inf')))
					if subgroup_copy == adjacent:
						if subgroup in proximals_copy:
							proximals_copy.remove(subgroup)
						continue#entire array has become a tandem array
					else:
						for sublist in adjacent_groups:#code for maintaining the relation link between tandem and proximal genes
							subgroup_copy_no_adjacent.append(sublist[-1])#Retain the last member of each converted tandem array in the proximal sublist to maintain the relation between the proximal and tandem sublists
						# Create a set of all genes in the dictionary for quick look-up
						all_genes = [gene for genes in pos_dic.values() for gene in genes]
						# Create a dictionary to map each gene to its index in the contig order
						gene_order = {gene: i for i, gene in enumerate(all_genes)}
						# Sort genes_list based on their position in gene_order
						subgroup_copy_no_adjacent.sort(key=lambda gene: gene_order.get(gene, float('inf')))

						if subgroup in proximals_copy:
							proximals_copy.remove(subgroup)
						proximals_copy.append(subgroup_copy_no_adjacent)

				else:
					# No tandems, but subgroup_copy has new non-adjacent genes
					if subgroup_copy != subgroup:
						# Create a set of all genes in the dictionary for quick look-up
						all_genes = [gene for genes in pos_dic.values() for gene in genes]
						# Create a dictionary to map each gene to its index in the contig order
						gene_order = {gene: i for i, gene in enumerate(all_genes)}
						# Sort genes_list based on their position in gene_order
						subgroup_copy.sort(key=lambda gene: gene_order.get(gene, float('inf')))
						if subgroup in proximals_copy:
							proximals_copy.remove(subgroup)
						proximals_copy.append(subgroup_copy)
	new_tandems_copy = copy.deepcopy(new_tandems)

	#code block to identify tandem arrays with common genes between the newly identified tandems and the old tandems list of lists using union find method
	def find(parent, x):
		if parent[x] != x:
			parent[x] = find(parent, parent[x])
		return parent[x]

	def union(parent, rank, x, y):
		root_x = find(parent, x)
		root_y = find(parent, y)
		if root_x != root_y:
			if rank[root_x] > rank[root_y]:
				parent[root_y] = root_x
			elif rank[root_x] < rank[root_y]:
				parent[root_x] = root_y
			else:
				parent[root_y] = root_x
				rank[root_x] += 1

	def merge_tandem_lists(new_tandems, tandems, pos_dic):
		sets = new_tandems + tandems
		parent = {}
		rank = {}

		# Build gene_positions dictionary from pos_dic
		gene_positions = {}
		for contig, contig_genes in pos_dic.items():
			for i, gene in enumerate(contig_genes):
				gene_positions[gene] = i  # Just position, or use (contig, i) if you need

		# Initialize union-find structure
		for group in sets:
			for gene in group:
				if gene not in parent:
					parent[gene] = gene
					rank[gene] = 0

		# Perform union operation for each group
		for group in sets:
			for i in range(len(group) - 1):
				union(parent, rank, group[i], group[i + 1])

		# Group genes by their representative parent
		merged_tandems = defaultdict(list)
		for gene in parent:
			root = find(parent, gene)
			merged_tandems[root].append(gene)

		# Sort individual groups based on their position in the contig
		sorted_groups = [
			sorted(group, key=lambda x: gene_positions.get(x, float('inf')))
			for group in merged_tandems.values()
		]

		# Sort the list of groups based on the position of the first gene in each group
		return sorted(sorted_groups, key=lambda x: gene_positions.get(x[0], float('inf')))

	merged_new_tandems_with_common_nodes = merge_tandem_lists(new_tandems, tandems, pos_dic)
	prefinal_tandems = merged_new_tandems_with_common_nodes
	#code block to identify tandem sublists with succeeding or preceding elements in other sublists
	def merge_tandem_groups(pos_dic, tandem_lists):
		# Step 1: Build the gene_positions dictionary using pos_dic
		gene_positions = {}
		for contig, contig_genes in pos_dic.items():
			for i, gene in enumerate(contig_genes):
				gene_positions[gene] = (contig, i)  # Store (contig, position) tuple

		# Step 2: Sort all sublists by the first element's position
		tandem_lists.sort(key=lambda x: gene_positions[x[0]][1])  # Sort by first gene position

		merged = []  # This will store the final merged sublists
		current_group = tandem_lists[0]  # Start with the first sublist

		# Step 3: Merging adjacent sublists
		for sublist in tandem_lists[1:]:
			# Get the last gene of the current group and the first gene of the next sublist
			last_gene_in_group = current_group[-1]
			first_gene_in_sublist = sublist[0]

			# Check if the last gene in current group is adjacent to the first gene in the next sublist
			last_gene_contig, last_gene_pos = gene_positions[last_gene_in_group]
			first_gene_contig, first_gene_pos = gene_positions[first_gene_in_sublist]

			# Merge them if they are adjacent on the same contig
			if last_gene_contig == first_gene_contig and last_gene_pos + 1 == first_gene_pos:
				current_group.extend(sublist)  # Merge sublists
			else:
				merged.append(current_group)  # Save the previous group
				current_group = sublist  # Start a new group

		# Append the last processed group to merged
		merged.append(current_group)
		return merged
	tandems_final = merge_tandem_groups(pos_dic, prefinal_tandems)#reclassified final tandems list of lists
	proximals_final = proximals_copy#reclassified final proximals list of lists
	tandem_set=set()
	proximal_set=set()

	for sub in tandems_final:
		for every in sub:
			tandem_set.add(str(every))
	for sub in proximals_final:
		for every in sub:
			proximal_set.add(str(every))
	final_set = tandem_set | proximal_set#final set of genes in the tandem and proximal arrays
	#code to produce a new list of lists of gene duplicates without the tandem and proximal genes
	# Separate into kept and rejected list of lists
	no_tandems_proximals_duplicates_final = []
	rejected_duplicates = []

	# Reject pairs if both the genes of a pair are in the final_set; Keep the other pairs in the master no tandems proximals final duplicates list
	for pair in duplicates_master_list:
		if pair[0] in final_set and pair[1] in final_set:
			rejected_duplicates.append(pair)
		else:
			no_tandems_proximals_duplicates_final.append(pair)
	no_tandems_proximals_duplicates_final_set = set()
	for sublist in no_tandems_proximals_duplicates_final:
		for each in sublist:
			no_tandems_proximals_duplicates_final_set.add(each)

	return tandems_final, proximals_final, no_tandems_proximals_duplicates_final, final_set, no_tandems_proximals_duplicates_final_set

#grouping dispersed duplicates
def dispersed_dups_grouping(no_tandems_no_proximals_duplicates):
	dispersed_dup_copy = no_tandems_no_proximals_duplicates.copy()
	# Create adjacency graph from the duplicates list
	# Build the adjacency graph, ensuring all genes are included
	graph = {}
	for gene1, gene2 in no_tandems_no_proximals_duplicates:
		graph.setdefault(gene1, []).append(gene2)
		graph.setdefault(gene2, []).append(gene1)

	# Perform DFS to identify connected components
	def dfs(gene, visited, current_group):
		visited.add(gene)
		current_group.append(gene)
		for neighbor in graph.get(gene, []):
			if neighbor not in visited:
				dfs(neighbor, visited, current_group)

	# Collapse connected pairs and include unconnected pairs/singletons
	visited = set()
	dispersed_arrays = []

	for gene in graph:
		if gene not in visited:
			current_group = []
			dfs(gene, visited, current_group)
			dispersed_arrays.append(sorted(set(current_group)))  # Remove duplicates and sort
	return dispersed_arrays

#Functions to merge sublists with genes having connections across the different gene duplication groups
def merge_sublists(*lists):
	seen = set()
	merged = []
	for lst in lists:
		for item in lst:
			if item not in seen:
				seen.add(item)
				merged.append(item)
	return merged

#pulling together gene duplicates with connections across different duplicates classifications
def link_tandems_proximals_dispersed (tandems, proximals, dispersed,classify):
	def find_merge_groups(all_sublists):
		merged = []
		used = set()

		for i in range(len(all_sublists)):
			if i in used:
				continue
			group = list(all_sublists[i][0])
			sources = [all_sublists[i][1]]
			indices = {i}
			used.add(i)
			merged_this_round = True

			while merged_this_round:
				merged_this_round = False
				for j in range(len(all_sublists)):
					if j not in used and any(e in group for e in all_sublists[j][0]):
						for e in all_sublists[j][0]:
							if e not in group:
								group.append(e)
						if all_sublists[j][1] not in sources:
							sources.append(all_sublists[j][1])
						indices.add(j)
						used.add(j)
						merged_this_round = True

			if len(indices) > 1:
				merged.append((group, sources, indices))
		return merged

	all_sublists = []
	for idx, sublist in enumerate(tandems):
		all_sublists.append((sublist, 'Tandems', idx))
	for idx, sublist in enumerate(proximals):
		all_sublists.append((sublist, 'Proximals', idx))
	for idx, sublist in enumerate(dispersed):
		all_sublists.append((sublist, 'Dispersed duplicates', idx))

	merged_groups = merged_groups = find_merge_groups(all_sublists)

	mixed = []
	merged_indices = set()
	for elements, sources, indices in merged_groups:
		source_str = "-".join(sources)
		mixed.append(elements + [source_str])
		merged_indices.update(indices)

	if classify == 'strict':
		tandems_new, proximals_new, dispersed_new = [], [], []
		for idx, (elems, src, old_idx) in enumerate(all_sublists):
			if idx not in merged_indices:
				if src == 'Tandems':
					tandems_new.append(elems)
				elif src == 'Proximals':
					proximals_new.append(elems)
				elif src == 'Dispersed duplicates':
					dispersed_new.append(elems)
		return tandems_new, proximals_new, dispersed_new, mixed
	else:
		return tandems, proximals, dispersed, mixed

#function to find the number of genes between genes in a list based on a dictionary
def build_gene_to_position(pos_dic):
	gene_to_position = {}
	for contig, genes in pos_dic.items():
		for idx, gene in enumerate(genes):
			gene_to_position[gene] = (contig, idx)
	return gene_to_position

def find_intervening_genes(gene_to_position,gene_list):
	# For each consecutive pair, calculate the distance
	gene_nos = []
	for i in range(len(gene_list) - 1):
		gene1 = gene_list[i]
		gene2 = gene_list[i + 1]

		if gene1 in gene_to_position and gene2 in gene_to_position:
			contig1, idx1 = gene_to_position[gene1]
			contig2, idx2 = gene_to_position[gene2]

			if contig1 == contig2:
				genes = abs(idx1 - idx2) - 1
			else:
				genes = "different_contig"  # or some special value if they are on different contigs

		gene_nos.append(genes)
	return gene_nos

#function to filter proximals output before writing them to the output file in case of mode being overlap
def output_subset_filter(duplicates_final):
	# Sort by length descending so we check bigger sets first
	sorted_lists = sorted(duplicates_final, key=lambda x: -len(x))

	seen = []  # Stores sets that we want to keep
	filtered = []

	for sublist in sorted_lists:
		sub_set = set(sublist)
		is_subset = False

		# Only check against larger or equal sets already accepted
		for seen_set in seen:
			if sub_set < seen_set:  # proper subset
				is_subset = True
				break

		if not is_subset:
			seen.append(sub_set)
			filtered.append(sublist)

	return filtered

#function to remove repeating tandems and proximals from dispersed list of lists
def short_distance_filter (duplicates_list, pos_dic,proximity,dup_type):
	cleaned_lists = []
	gene_to_position = build_gene_to_position(pos_dic)

	for sublist in duplicates_list:

		genes = sublist[:]#making a copy of sublist as genes

		while True:
			intervening = find_intervening_genes(gene_to_position, genes)
			to_remove = None

			for i, count in enumerate(intervening):
				if isinstance(count, str):
					continue  # Skip 'different' or any non-integer entry

				if dup_type == 'dispersed':
					if count < (proximity+1):
						to_remove = i + 1  # Remove the second gene in the pair
						break
				if dup_type == 'proximal':
					if count < 1:
						to_remove = i + 1  # Remove the second gene in the pair
						break

			if to_remove is not None:
				del genes[to_remove]
			else:
				break  # No more deletions needed

		cleaned_lists.append(genes)

	return cleaned_lists


#Converts the gene pairs in the forward blast best hits file to a dictionary such that the query genes in the first columns are the keys and the reference genes in the second column are the values
def fwd_best_hits_to_dic(best_hits_file):
	fwd_best_hits_dic={}
	with open(best_hits_file,'r')as f:
		line=f.readline()
		while line:
			parts=line.strip().split('\t')
			fwd_best_hits_dic.update({parts[0]:parts[1]})
			line=f.readline()
	return fwd_best_hits_dic

#Converts the gene pairs in the orthologs file to a dictionary such that the query genes in the first columns are the keys and the orthologous reference genes in the second column are the values
def orthologs_hits_to_dic(orthologs_file):
	orth_hits_dic={}
	with open(orthologs_file,'r')as f:
		f.readline()
		line=f.readline()
		while line:
			parts=line.strip().split('\t')
			if parts[1]!='--':
				orth_hits_dic.update({parts[0]:parts[1]})
			line=f.readline()
	return orth_hits_dic

#function to get the downstream flanking genes of a specified gene
def left_flank(target_gene, gene_pos_per_chr,flank_number):
	l=1
	org_left=[]
	for contig in gene_pos_per_chr:
		if target_gene in gene_pos_per_chr[contig]:
			gene = target_gene
			ind = gene_pos_per_chr[contig].index(gene)
			length = len(gene_pos_per_chr[contig])
			pos = ind + 1
			p = pos - 1
			q = length - pos
			if ind == 0:  # if the target gene is at the start of a contig
				looper = 0
				org_left = []
			elif ind == (length - 1) and length != 1:  # if the target gene is at the end of a contig
				if length > int(flank_number):
					looper = int(flank_number)
					while l <= looper:
						org_left.append(str(gene_pos_per_chr[contig][ind - l]))
						l += 1
				else:
					looper = length
					while l < looper:
						org_left.append(str(gene_pos_per_chr[contig][ind - l]))
						l += 1
			elif ind != 0 and ind != (length - 1):  # if the target gene occupies a non-end position in a contig
				if p > int(flank_number):
					looper = int(flank_number)
				elif p <= int(flank_number):
					looper = p
				while l <= looper:
					org_left.append(str(gene_pos_per_chr[contig][ind - l]))
					l += 1
	return org_left

#function to get the upstream flanking genes of a specified gene
def right_flank(target_gene, gene_pos_per_chr,flank_number):
	r=1
	org_right=[]
	for contig in gene_pos_per_chr:
		if target_gene in gene_pos_per_chr[contig]:
			gene = target_gene
			ind = gene_pos_per_chr[contig].index(gene)
			length = len(gene_pos_per_chr[contig])
			pos = ind + 1
			p = pos - 1
			q = length - pos
			if ind == 0 and length != 1:  # if the target gene is at the start of a contig
				if length > int(flank_number):
					looper = int(flank_number)
					while r <= looper:
						org_right.append(str(gene_pos_per_chr[contig][ind + r]))
						r += 1
				else:
					looper = length
					while r < looper:
						org_right.append(str(gene_pos_per_chr[contig][ind + r]))
						r += 1
			elif ind == (length - 1):  # if the target gene is at the end of a contig
				looper = 0
				org_right = []
			elif ind != 0 and ind != (length - 1):  # if the target gene occupies a non-end position in a contig
				if q > int(flank_number):
					looper = int(flank_number)
				elif q <= int(flank_number):
					looper = q
				while r <= looper:
					org_right.append(str(gene_pos_per_chr[contig][ind + r]))
					r += 1
	return org_right

#Checking for synteny in the specified duplicate and reference gene regions
def synteny_checker(anchor_list,sublist,gene_pos_per_chr_q, gene_pos_per_chr_r, synteny_cutoff, flank_number, best_side, fwd_best_hits_dic):
	orgq_left = []
	orgq_right = []
	orgr_left = []
	orgr_right = []
	if len(sublist)==1:
		only = sublist[0]
		# collecting the left and right flanking genes of the given gene
		orgq_left = left_flank(only, gene_pos_per_chr_q, flank_number)
		orgq_right = right_flank(only, gene_pos_per_chr_q, flank_number)
	elif len(sublist)>1:
		start = sublist[0]
		end = sublist[-1]
		#collecting the left and right flanking genes of the given duplicates list
		orgq_left = left_flank(start,gene_pos_per_chr_q,flank_number)
		orgq_right=right_flank(end,gene_pos_per_chr_q,flank_number)
	anchor_nos=len(anchor_list)
	if anchor_nos==1:
		anchor = anchor_list[0]
		#collecting the left and right flanking genes of a given single anchor gene in reference
		orgr_left = left_flank(anchor,gene_pos_per_chr_r,flank_number)
		orgr_right = right_flank(anchor,gene_pos_per_chr_r,flank_number)
	elif anchor_nos>1:
		anchor_start=anchor_list[0]
		anchor_end=anchor_list[-1]
		#collecting the left and right flanking genes of multiple anchor genes in the reference
		orgr_left = left_flank(anchor_start, gene_pos_per_chr_r,flank_number)
		orgr_right = right_flank(anchor_end, gene_pos_per_chr_r,flank_number)
	#counting the number of flanking gene pairs that show up as best hits in the forward blast of the query vs reference
	orgq_flank=orgq_left+orgq_right
	orgr_flank=orgr_left+orgr_right
	best_left = 0
	best_right = 0
	cross_flank_left = 0
	cross_flank_right = 0
	for flank in orgq_flank:
		if flank in fwd_best_hits_dic:
			if fwd_best_hits_dic[flank] in orgr_flank:
				#code segment to account hits for synteny in case of inversions between gene duplicate groups and the corresponding reference gene
				if (flank in orgq_left and fwd_best_hits_dic[flank] in orgr_left): #to ensure that the flanking genes forward best hits correspond to the left flanking genes correspondingly in the reference organism
					best_left += 1
				elif (flank in orgq_right and fwd_best_hits_dic[flank] in orgr_right): #to ensure that the flanking genes forward best hits correspond to the right flanking genes correspondingly in the reference organism
					best_right += 1
				elif (flank in orgq_left and fwd_best_hits_dic[flank] in orgr_right):
					cross_flank_left += 1
				elif (flank in orgq_right and fwd_best_hits_dic[flank] in orgr_left):
					cross_flank_right += 1
	total_best_left = best_left + cross_flank_left
	total_best_right = best_right + cross_flank_right
	total_best = total_best_left + total_best_right

	if len(orgq_flank) > 0:
		ratio = float(total_best)/float(len(orgq_flank))
	else:
		ratio = 0
	if (ratio >= synteny_cutoff) and (total_best_left >= int(best_side)) and (total_best_right >= int(best_side)):
		synteny ='yes'
	else:
		synteny ='no'
	return synteny


#Functions to find the ka,ks values of gene pairs

# Constants
CODON = 64
XSIZE = 16
DNASIZE = 4

# Genetic code table (standard genetic code)
codon_table = {
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
	'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def convertChar(ch):
	"""Convert nucleotide character to number: T/U->0, C->1, A->2, G->3"""
	if ch in ['T', 't', 'U', 'u']:
		return 0
	elif ch in ['C', 'c']:
		return 1
	elif ch in ['A', 'a']:
		return 2
	elif ch in ['G', 'g']:
		return 3
	else:
		return -1


def getID(codon):
	"""Get codon ID (0-63)"""
	if len(codon) != 3:
		return -1
	id = 0
	for i in range(3):
		n = convertChar(codon[i])
		if n < 0:
			return -1
		id = id * 4 + n
	return id


def getAminoAcid(codon_id):
	"""Get amino acid from codon ID"""
	if codon_id < 0 or codon_id >= 64:
		return '*'
	# Convert ID back to codon string
	nucleotides = ['T', 'C', 'A', 'G']
	codon = ''
	temp_id = codon_id
	for i in range(3):
		codon = nucleotides[temp_id % 4] + codon
		temp_id //= 4
	return codon_table.get(codon, '*')


def initArray(arr, size):
	"""Initialize array with zeros"""
	return [0.0] * size


def sumArray(arr, size):
	"""Sum array elements"""
	return sum(arr[:size])


def square(x):
	"""Square a number"""
	return x * x


def SIGN(a, b):
	"""Return abs(a) with sign of b"""
	return abs(a) if b >= 0.0 else -abs(a)

class MYN:
	# Class-level constants
	CODON = 64
	XSIZE = 16
	DNASIZE = 4

	def __init__(self):
		self.name = "MYN"
		self.GAMMA = 0
		self.kappa = 0
		self.kappatc = 0
		self.kappaag = 0
		self.KAPPA = [0, 0]
		self.omega = 0
		self.t = 0
		self.S = 0
		self.N = 0
		self.Sd = 0
		self.Nd = 0
		self.snp = 0
		self.pi = initArray([], 64)
		self.pi_sqrt = initArray([], 64)
		self.npi0 = 0
		self.f12pos = initArray([], 12)
		self.iteration = 1  # From YN00 base class

	def GetKappa(self, seq1, seq2):
		"""Get the two kappas between purines and between pyrimidines"""
		F = [[0.0] * self.XSIZE for _ in range(2)]
		S = [0.0, 0.0]
		wk = [0.0, 0.0]
		pi4 = [0.0] * 4
		kappatc_TN = [0.0, 0.0]
		kappaag_TN = [0.0, 0.0]
		kappa_TN = [0.0, 0.0]
		kdefault = 2

		by = [16, 4, 1]

		# Get Pi[] of A,C,G,T
		for h in range(0, len(seq1), 3):
			# c[]: amino acid(0--63)
			c = [getID(seq1[h:h + 3]), getID(seq2[h:h + 3])]
			# aa[]: amino acid
			aa = [getAminoAcid(c[0]), getAminoAcid(c[1])]
			# b[][]: 0--3
			b = [[0, 0, 0], [0, 0, 0]]
			for j in range(3):
				b[0][j] = convertChar(seq1[h + j])
				b[1][j] = convertChar(seq2[h + j])

			# Find non-degenerate sites
			for pos in range(3):
				nondeg = 0
				for k in range(2):
					found = False
					for i in range(4):
						if i != b[k][pos]:
							if getAminoAcid(c[k] + (i - b[k][pos]) * by[pos]) == aa[k]:
								found = True
								break
					if not found:
						nondeg += 1

				# F[0][]: 0-fold
				if nondeg == 2:
					F[0][b[0][pos] * 4 + b[1][pos]] += 0.5
					F[0][b[1][pos] * 4 + b[0][pos]] += 0.5

			# Find 4-fold degenerate sites at 3rd position
			fourdeg = 0
			for k in range(2):
				i = c[k] - b[k][2]
				all_same = True
				for j in range(4):
					if j != b[k][2] and getAminoAcid(i + j) != aa[k]:
						all_same = False
						break
				if aa[0] == aa[1] and all_same:
					fourdeg += 1

			# F[1][]: 4-fold
			if fourdeg == 2:
				F[1][b[0][2] * 4 + b[1][2]] += 0.5
				F[1][b[1][2] * 4 + b[0][2]] += 0.5

		for k in range(2):  # two kinds of sites
			S[k] = sumArray(F[k], 16)
			if S[k] <= 0:
				wk[k] = 0
				continue

			for j in range(16):
				F[k][j] /= S[k]

			# Transitions between purines
			T1 = 2 * F[k][2 * self.DNASIZE + 3]
			# Transitions between pyrimidines
			T2 = 2 * F[k][0 * self.DNASIZE + 1]
			# Transversions
			V = 1 - T1 - T2 - (F[k][0 * self.DNASIZE + 0] + F[k][1 * self.DNASIZE + 1] +
							   F[k][2 * self.DNASIZE + 2] + F[k][3 * self.DNASIZE + 3])

			# pi[]: the sum probability of T, C, A, G, respectively
			for j in range(4):
				pi4[j] = sumArray(F[k][j * 4:], 4)

			kappatc_TN[k], kappaag_TN[k] = self.CorrectKappaTN93(S[k], T1, T2, V, pi4)
			wk[k] = S[k] if (kappatc_TN[k] > 0 and kappaag_TN[k] > 0) else 0

			# R = (πT*πC*κ1 + πA*πG*κ2)/(πY*πR), kappa = 2R in PAML's DOC
			kappa_TN[k] = 2 * (kappatc_TN[k] * pi4[0] * pi4[1] + kappaag_TN[k] * pi4[2] * pi4[3]) / (
						(pi4[0] + pi4[1]) * (pi4[2] + pi4[3]))

		if wk[0] + wk[1] == 0:
			self.kappatc = self.kappaag = self.kappa = kdefault
		else:
			self.kappatc = (kappatc_TN[0] * wk[0] + kappatc_TN[1] * wk[1]) / (wk[0] + wk[1])
			self.kappaag = (kappaag_TN[0] * wk[0] + kappaag_TN[1] * wk[1]) / (wk[0] + wk[1])
			self.kappa = (kappa_TN[0] * wk[0] + kappa_TN[1] * wk[1]) / (wk[0] + wk[1])

		self.KAPPA[0] = self.kappatc
		self.KAPPA[1] = self.kappaag

		return 0

	def GetPMatCodon(self, PMatrix, kappa, omega):
		"""Calculate transition probability matrix(64*64)"""
		by = [16, 4, 1]
		U = [0.0] * (self.CODON * self.CODON)
		V = [0.0] * (self.CODON * self.CODON)
		Root = [0.0] * (self.CODON * self.CODON)

		# Initialize PMatrix
		for i in range(self.CODON * self.CODON):
			PMatrix[i] = 0.0

		for i in range(self.CODON):
			for j in range(i):
				# codon 'from'
				from_codon = [i // 16, (i // 4) % 4, i % 4]
				# codon 'to'
				to_codon = [j // 16, (j // 4) % 4, j % 4]
				# amino acid of 'from' and 'to'
				c = [getAminoAcid(i), getAminoAcid(j)]
				# stop codon
				if c[0] == '*' or c[1] == '*':
					continue

				# whether two codons only have one difference
				ndiff = 0
				pos = 0
				for k in range(3):
					if from_codon[k] != to_codon[k]:
						ndiff += 1
						pos = k

				if ndiff == 1:
					# only have one difference
					PMatrix[i * self.CODON + j] = 1
					# transition
					if (from_codon[pos] + to_codon[pos] - 1) * (from_codon[pos] + to_codon[pos] - 5) == 0:
						if from_codon[pos] + to_codon[pos] == 1:
							PMatrix[i * self.CODON + j] *= self.kappatc  # T<->C
						else:
							PMatrix[i * self.CODON + j] *= self.kappaag  # A<->G

					# nonsynonymous
					if c[0] != c[1]:
						PMatrix[i * self.CODON + j] *= omega

					# diagonal element is equal
					PMatrix[j * self.CODON + i] = PMatrix[i * self.CODON + j]

		# PMatrix[](*Q): transition probability matrix
		for i in range(self.CODON):
			for j in range(self.CODON):
				PMatrix[i * self.CODON + j] *= self.pi[j]

		# scale the sum of PMat[][j](j=0-63) to zero
		mr = 0
		for i in range(self.CODON):
			PMatrix[i * self.CODON + i] = -sumArray(PMatrix[i * self.CODON:], self.CODON)
			# The sum of transition probability of main diagonal elements
			mr -= self.pi[i] * PMatrix[i * self.CODON + i]

		# Safety check for potential error:
		if abs(mr) < 1e-10:
			return -1

		# calculate exp(PMatrix*t)
		status = self.eigenQREV(PMatrix, self.pi, self.pi_sqrt, self.CODON, self.npi0, Root, U, V)
		if status != 0:
			return -1
		for i in range(self.CODON):
			Root[i] /= mr
		self.PMatUVRoot(PMatrix, self.t, 64, U, V, Root)

		return 0

	def CorrectKappaTN93(self, n, P1, P2, Q, pi4):
		"""Correct kappas"""
		Qsmall = min(1e-10, 0.1 / n)
		default_kappa = 2
		maxkappa = 99

		kappatc_TN93 = -1
		kappaag_TN93 = -1

		failTN93 = 0

		Y = pi4[0] + pi4[1]
		R = pi4[2] + pi4[3]

		tc = pi4[0] * pi4[1]
		ag = pi4[2] * pi4[3]

		if (P1 + P2 + Q) > 1.0 or (P1 + P2) < -1e-10 or Q < -1e-10 or abs(Y + R - 1) > 1e-8:
			return default_kappa, default_kappa  # Return default values instead of -1

		if Q < Qsmall:
			failTN93 = 1
		elif Y <= 0 or R <= 0 or (tc <= 0 and ag <= 0):
			failTN93 = 1
		else:  # TN93 for multiple substitutions
			A = tc / Y + ag / R
			B = tc + ag
			C = Y * R

			# Add safety checks for division by zero
			if ag > 0 and tc > 0 and R > 0 and Y > 0:
				a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R)
				a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y)
				b = 1 - Q / (2 * C)

				if a1 < 0 or a2 < 0 or b < 0:
					failTN93 = 1
				else:
					a1 = math.log(a1)
					a2 = math.log(a2)
					b = math.log(b)
					# Kappa
					if R * b != 0:  # Additional safety check
						kappaag_TN93 = (Y * b - a1) / (-R * b)
					else:
						failTN93 = 1

					if Y * b != 0:  # Additional safety check
						kappatc_TN93 = (R * b - a2) / (-Y * b)
					else:
						failTN93 = 1
			else:
				failTN93 = 1

		if failTN93:  # Fail to correct kappa
			kappatc_TN93 = default_kappa
			kappaag_TN93 = default_kappa

		return kappatc_TN93, kappaag_TN93

	def CorrectKaksTN93(self, n, P1, P2, Q, pi4):
		"""Correct Ka and Ks"""
		Qsmall = 1e-10
		failTN93 = 0
		kaks = 0
		SEkaks = 0

		Y = pi4[0] + pi4[1]
		R = pi4[2] + pi4[3]

		tc = pi4[0] * pi4[1]
		ag = pi4[2] * pi4[3]

		if P1 + P2 + Q > 1 or abs(Y + R - 1) > Qsmall or Y <= 0 or R <= 0 or (tc <= 0 and ag <= 0):
			failTN93 = 1
		else:  # TN93 for multiple substitutions
			A = tc / Y + ag / R
			B = tc + ag
			C = Y * R

			# Add safety checks for division by zero
			if ag > 0 and tc > 0 and R > 0 and Y > 0:
				a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R)
				a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y)
				b = 1 - Q / (2 * C)

				if a1 < 0 or a2 < 0 or b < 0:
					failTN93 = 1
				else:
					a1 = math.log(a1)
					a2 = math.log(a2)
					b = math.log(b)
					# Ka or Ks
					kaks = (-2 * ag * a1 / R) + (-2 * tc * a2 / Y) + (-2 * (C - ag * Y / R - tc * R / Y) * b)

					# Calculate standard errors with safety checks
					denom1 = 2 * ag * R - R * R * P1 - ag * Q
					denom2 = 2 * tc * Y - Y * Y * P2 - tc * Q
					denom3 = 2 * R * R * Y * Y - R * Y * Q

					if denom1 != 0 and denom2 != 0 and denom3 != 0:
						cc1 = 2 * ag * R / denom1
						cc2 = 2 * tc * Y / denom2
						cc3 = 2 * ag * ag / (R * denom1)
						cc3 += 2 * tc * tc / (Y * denom2)
						cc3 += (R * R * (Y * Y - 2 * tc) + Y * Y * (R * R - 2 * ag)) / denom3
						SEkaks = (square(cc1) * P1 + square(cc2) * P2 + square(cc3) * Q - square(
							cc1 * P1 + cc2 * P2 + cc3 * Q)) / n
					else:
						failTN93 = 1
			else:
				failTN93 = 1

		if failTN93 == 1:  # Use YN00's correction for Ka, Ks
			kaks, t, SEkaks = self.DistanceF84(n, P1 + P2, Q, pi4, Qsmall)

		return kaks, SEkaks

	def CountDiffs(self, seq1, seq2, PMatrix):
		"""Count differences, considering different transitional pathways between purines and between pyrimidines"""
		by = [16, 4, 1]
		Sdts1 = Sdts2 = Sdtv = Ndts1 = Ndts2 = Ndtv = 0
		self.snp = 0

		for h in range(0, len(seq1), 3):
			c = [getID(seq1[h:h + 3]), getID(seq2[h:h + 3])]
			# Difference?
			if c[0] == c[1]:
				continue

			b = [[0, 0, 0], [0, 0, 0]]
			aa = ['', '']
			for i in range(2):
				b[i][0] = c[i] // 16
				b[i][1] = (c[i] % 16) // 4
				b[i][2] = c[i] % 4
				aa[i] = getAminoAcid(c[i])

			# ndiff: differences of two codons
			ndiff = 0
			sts1 = sts2 = stv = nts1 = nts2 = ntv = 0
			# dmark[]: position of different codon
			dmark = [-1, -1, -1]
			for k in range(3):
				if b[0][k] != b[1][k]:
					dmark[ndiff] = k
					ndiff += 1

			self.snp += ndiff

			npath = 1
			if ndiff > 1:
				npath = 2 if ndiff == 2 else 6

			if ndiff == 1:
				transi = b[0][dmark[0]] + b[1][dmark[0]]
				if aa[0] == aa[1]:
					if transi == 5:
						sts1 += 1
					elif transi == 1:
						sts2 += 1
					else:
						stv += 1
				else:
					if transi == 5:
						nts1 += 1
					elif transi == 1:
						nts2 += 1
					else:
						ntv += 1
			else:  # ndiff=2 or 3
				nstop = 0
				sts1path = [0] * 6
				sts2path = [0] * 6
				stvpath = [0] * 6
				nts1path = [0] * 6
				nts2path = [0] * 6
				ntvpath = [0] * 6
				ppath = [0.0] * 6

				for k in range(npath):
					# set the step[]
					step = [-1, -1, -1]
					if ndiff == 2:
						step[0] = dmark[k]
						step[1] = dmark[1 - k]
					else:
						step[0] = k // 2
						step[1] = k % 2
						if step[0] <= step[1]:
							step[1] += 1
						step[2] = 3 - step[0] - step[1]

					bt1 = b[0][:]
					bt2 = b[0][:]

					sts1path[k] = sts2path[k] = stvpath[k] = nts1path[k] = nts2path[k] = ntvpath[k] = 0

					# ppath[]: probability of each path
					ppath[k] = 1
					for i1 in range(ndiff):
						bt2[step[i1]] = b[1][step[i1]]

						# ct[]: mutated codon's ID(0--63)
						ct = [0, 0]
						for i2 in range(3):
							ct[0] += bt1[i2] * by[i2]
							ct[1] += bt2[i2] * by[i2]

						# ppath[k]: probability of path k
						ppath[k] *= PMatrix[ct[0] * self.CODON + ct[1]]
						aa_temp = [getAminoAcid(ct[0]), getAminoAcid(ct[1])]

						if aa_temp[1] == '*':
							nstop += 1
							ppath[k] = 0
							break

						transi = b[0][step[i1]] + b[1][step[i1]]

						# ts & tr when syn & nonsyn in path k
						if aa_temp[0] == aa_temp[1]:
							if transi == 5:
								sts1path[k] += 1
							elif transi == 1:
								sts2path[k] += 1
							else:
								stvpath[k] += 1
						else:
							if transi == 5:
								nts1path[k] += 1
							elif transi == 1:
								nts2path[k] += 1
							else:
								ntvpath[k] += 1

						bt1 = bt2[:]

				if npath == nstop:  # all paths through stop codons
					if ndiff == 2:
						nts1 = 0.25
						nts2 = 0.25
						ntv = 1.5
					else:
						nts1 = 0.25
						nts2 = 0.25
						ntv = 2.5
				else:
					# sum probability of all paths
					sump = sumArray(ppath, npath)
					if sump > 1e-20:
						for k in range(npath):
							p = ppath[k] / sump
							sts1 += sts1path[k] * p
							sts2 += sts2path[k] * p
							stv += stvpath[k] * p
							nts1 += nts1path[k] * p
							nts2 += nts2path[k] * p
							ntv += ntvpath[k] * p

			Sdts1 += sts1
			Sdts2 += sts2
			Sdtv += stv
			Ndts1 += nts1
			Ndts2 += nts2
			Ndtv += ntv

		return Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv

	# Main function for calculating Ka and Ks
	def DistanceYN00(self, seq1, seq2):
		nround = 100
		status = 1
		fbS = [0.0] * 4
		fbN = [0.0] * 4
		fbSt = [0.0] * 4
		fbNt = [0.0] * 4
		accu = 5e-8
		minomega = 1e-5
		maxomega = 99
		PMatrix = [0.0] * (self.CODON * self.CODON)

		# initial values for t and omega(Ka/Ks)
		self.t = 0.09
		self.omega = 0.5
		self.S = self.N = 0.0

		# Get frequency
		self.getFrequency(seq1, seq2)

		# Count sites of sequence 1
		St, Nt, fbSt, fbNt = self.CountSites(seq1)
		self.S += St / 2
		self.N += Nt / 2
		for j in range(4):
			fbS[j] += fbSt[j] / 2
			fbN[j] += fbNt[j] / 2

		# Count sites of sequence 2
		St, Nt, fbSt, fbNt = self.CountSites(seq2)
		self.S += St / 2
		self.N += Nt / 2
		for j in range(4):
			fbS[j] += fbSt[j] / 2
			fbN[j] += fbNt[j] / 2

		# Iterative loop
		w0 = 0
		dS0 = 0
		dN0 = 0

		for ir in range(nround):
			# Get transition probability matrix from one codon to another
			matrix_result = self.GetPMatCodon(PMatrix, self.kappa, self.omega)

			# Count differences
			Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv = self.CountDiffs(seq1, seq2, PMatrix)

			# Synonymous(Sd) and nonsynonymous(Nd) differences
			self.Sd = Sdts1 + Sdts2 + Sdtv
			self.Nd = Ndts1 + Ndts2 + Ndtv

			# Seldom happen
			if self.Sd > self.S:
				Sdts1 *= (self.S / self.Sd)
				Sdts2 *= (self.S / self.Sd)
				Sdtv *= (self.S / self.Sd)
			if self.Nd > self.N:
				Ndts1 *= (self.N / self.Nd)
				Ndts2 *= (self.N / self.Nd)
				Ndtv *= (self.N / self.Nd)

			# Ks
			dS, SEdS = self.CorrectKaksTN93(self.S, Sdts1 / self.S, Sdts2 / self.S, Sdtv / self.S, fbS)
			# Ka
			dN, SEdN = self.CorrectKaksTN93(self.N, Ndts1 / self.N, Ndts2 / self.N, Ndtv / self.N, fbN)

			status = -1

			if dS < 1e-9:
				status = -1
				self.omega = maxomega
			else:
				self.omega = max(minomega, dN / dS)

			self.t = dS * 3 * self.S / (self.S + self.N) + dN * 3 * self.N / (self.S + self.N)

			if abs(dS - dS0) < accu and abs(dN - dN0) < accu and abs(self.omega - w0) < accu:
				break

			dS0 = dS
			dN0 = dN
			w0 = self.omega

		#Check for convergence failure
		if ir == nround - 1:  # C code checks ir==nround
			status = -2

		return dS, dN, SEdS, SEdN, status
# function to count the synonymous and nonsynonymous sites of a sequence
	def CountSites(self, seq):
		by = [16, 4, 1]
		Stot = Ntot = 0
		fbS = [0.0] * 4
		fbN = [0.0] * 4

		for h in range(0, len(seq), 3):
			# Get codon id and amino acid
			c = [getID(seq[h:h + 3]), 0]
			aa = [getAminoAcid(c[0]), '']
			b = [0, 0, 0]
			for i in range(3):
				b[i] = convertChar(seq[h + i])

			S = N = 0
			for j in range(3):
				for k in range(4):  # b[j] changes to k
					if k == b[j]:
						continue
					# c[0] change at position j
					c[1] = c[0] + (k - b[j]) * by[j]
					aa[1] = getAminoAcid(c[1])

					if aa[1] == '*':
						continue

					r = self.pi[c[1]]
					if k + b[j] == 1 or k + b[j] == 5:  # transition
						if k + b[j] == 1:
							r *= self.kappatc
						else:
							r *= self.kappaag  # (k+b[j]==5)

					if aa[0] == aa[1]:  # synonymous
						S += r
						fbS[b[j]] += r  # syn probability of A,C,G,T
					else:  # nonsynonymous
						N += r
						fbN[b[j]] += r  # nonsyn probability of A,C,G,T

			Stot += S
			Ntot += N

		# Scale Stot+Ntot to seq.length()
		r = len(seq) / (Stot + Ntot)
		Stot *= r
		Ntot *= r

		# get probability of syn of four nul.
		r = sumArray(fbS, 4)
		for k in range(4):
			fbS[k] /= r

		# get probability of nonsyn of four nul.
		r = sumArray(fbN, 4)
		for k in range(4):
			fbN[k] /= r

		return Stot, Ntot, fbS, fbN

	# Additional methods from YN00 that MYN uses
	def getFrequency(self, seq1, seq2):
		"""Get A,C,G,T frequency at three positions"""
		fstop = 0.0

		# Initialize f12pos
		self.f12pos = [0.0] * 12

		# Get A,C,G,T frequency at three positions
		for i in range(len(seq1)):
			self.f12pos[(i % 3) * 4 + convertChar(seq1[i])] += 1
			self.f12pos[(i % 3) * 4 + convertChar(seq2[i])] += 1

		for i in range(12):
			self.f12pos[i] /= (len(seq1) + len(seq2)) / 3

		# Get 64 amino acid probability
		for i in range(self.CODON):
			self.pi[i] = self.f12pos[i // 16] * self.f12pos[4 + (i % 16) // 4] * self.f12pos[8 + i % 4]
			if getAminoAcid(i) == '*':
				fstop += self.pi[i]
				self.pi[i] = 0

		# Scale the sum of pi[] to 1
		for i in range(self.CODON):
			self.pi[i] /= (1.0 - fstop)

		if abs(1 - sumArray(self.pi, self.CODON)) > 1e-6:
			sys.stdout.write("Warning: error in get codon frequency.")
			sys.stdout.flush()

		# Create compressed pi_sqrt array
		epsilon = 1e-10
		compressed_index = 0
		for i in range(self.CODON):
			if self.pi[i] > epsilon:
				self.pi_sqrt[compressed_index] = math.sqrt(self.pi[i])
				compressed_index += 1
			elif self.pi[i] < 0:
				self.pi[i] = 0  # Ensure no negative values

		# npi0 is the number of codons with zero frequency
		self.npi0 = self.CODON - compressed_index

	#Correct for multiple substitutions
	def DistanceF84(self, n, P, Q, pi4, Qsmall):
		failF84 = 0
		failK80 = 0
		failJC69 = 0
		maxkappa = 2
		maxt = 99

		k_HKY = -1
		t = 0
		SEt = 0

		Y = pi4[0] + pi4[1]
		R = pi4[2] + pi4[3]

		tc = pi4[0] * pi4[1]
		ag = pi4[2] * pi4[3]

		if self.GAMMA == 4 or self.GAMMA == -1:
			self.name = "GYN"

		if P + Q > 1:
			t = maxt
			k_HKY = 1
			return k_HKY, t, SEt

		if P < -1e-10 or Q < -1e-10 or abs(Y + R - 1) > 1e-8:
			return k_HKY, t, SEt

		# HKY85
		if Q < Qsmall:
			failF84 = failK80 = 1
		elif Y <= 0 or R <= 0 or (tc <= 0 and ag <= 0):
			failF84 = 1
		else:
			A = tc / Y + ag / R
			B = tc + ag
			C = Y * R
			a = (2 * B + 2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * C)) - P) / (2 * A)
			b = 1 - Q / (2 * C)
			if a <= 0 or b <= 0:
				failF84 = 1

		if not failF84:
			if self.GAMMA == 4 or self.GAMMA == 20:
				a = 0.5 * self.GAMMA * (pow(a, -1.0 / self.GAMMA) - 1)
				b = 0.5 * self.GAMMA * (pow(b, -1.0 / self.GAMMA) - 1)
			else:
				a = -0.5 * math.log(a)
				b = -0.5 * math.log(b)

			if b <= 0:
				failF84 = 1
			else:
				k_F84 = a / b - 1
				t = 4 * b * (tc * (1 + k_F84 / Y) + ag * (1 + k_F84 / R) + C)
				k_HKY = (B + (tc / Y + ag / R) * k_F84) / B  # k_F84=>k_HKY85

				# Standard errors
				a = A * C / (A * C - C * P / 2 - (A - B) * Q / 2)
				b = A * (A - B) / (A * C - C * P / 2 - (A - B) * Q / 2) - (A - B - C) / (C - Q / 2)
				SEt = math.sqrt((a * a * P + b * b * Q - square(a * P + b * Q)) / n)

		# K80
		if failF84 and not failK80:
			a = 1 - 2 * P - Q
			b = 1 - 2 * Q
			if a <= 0 or b <= 0:
				failK80 = 1
			else:
				a = -math.log(a)
				b = -math.log(b)
				if b <= 0:
					failK80 = 1
				else:
					k_HKY = (0.5 * a - 0.25 * b) / (0.25 * b)
					t = 0.5 * a + 0.25 * b
				if SEt:
					a = 1 / (1 - 2 * P - Q)
					b = (a + 1 / (1 - 2 * Q)) / 2
					SEt = math.sqrt((a * a * P + b * b * Q - square(a * P + b * Q)) / n)

		if failK80:  # try JC69
			P += Q
			if P >= 0.75:
				failJC69 = 1
				P = 0.75 * (n - 1) / n
			t = -0.75 * math.log(1 - P * 4 / 3)

			if t > 99:
				t = maxt
			if SEt:
				SEt = math.sqrt(9 * P * (1 - P) / n) / (3 - 4 * P)

		if k_HKY > 99:
			k_HKY = maxkappa

		return t, t, SEt  # Return t for kaks as well

	# function for matrix operations from YN00
	def eigenQREV(self, Q, pi, pi_sqrt, n, npi0, Root, U, V):
		"""Find eigen solution of rate matrix Q for time-reversible Markov process"""
		nnew = n - npi0

		if npi0 == 0:  # All codons have non-zero frequency
			# Set U[64*64]
			for i in range(n):
				U[i * n + i] = Q[i * n + i]
				for j in range(i):
					U[i * n + j] = U[j * n + i] = Q[i * n + j] * pi_sqrt[i] / pi_sqrt[j]

			# Set U[64*64]
			status = self.eigenRealSym(U, n, Root, V)

			for i in range(n):
				for j in range(n):
					V[i * n + j] = U[j * n + i] * pi_sqrt[j]

			for i in range(n):
				for j in range(n):
					U[i * n + j] /= pi_sqrt[i]
		else:
			# Compressed case - pi_sqrt contains only non-zero values
			inew = 0
			for i in range(n):
				if pi[i] > 0:
					jnew = 0
					for j in range(i):
						if pi[j] > 0:
							U[inew * nnew + jnew] = U[jnew * nnew + inew] = Q[i * n + j] * pi_sqrt[inew] / pi_sqrt[jnew]
							jnew += 1
					U[inew * nnew + inew] = Q[i * n + i]
					inew += 1

			status = self.eigenRealSym(U, nnew, Root, V)

			# construct Root
			inew = nnew - 1
			for i in range(n - 1, -1, -1):
				Root[i] = Root[inew] if pi[i] > 0 else 0
				if pi[i] > 0:
					inew -= 1

			# construct V
			inew = nnew - 1
			for i in range(n - 1, -1, -1):
				if pi[i] > 0:
					jnew = nnew - 1
					for j in range(n - 1, -1, -1):
						if pi[j] > 0:
							V[i * n + j] = U[jnew * nnew + inew] * pi_sqrt[jnew]
							jnew -= 1
						else:
							V[i * n + j] = 1 if i == j else 0
					inew -= 1
				else:
					for j in range(n):
						V[i * n + j] = 1 if i == j else 0

			# construct U
			inew = nnew - 1
			for i in range(n - 1, -1, -1):
				if pi[i] > 0:
					jnew = nnew - 1
					for j in range(n - 1, -1, -1):
						if pi[j] > 0:
							U[i * n + j] = U[inew * nnew + jnew] / pi_sqrt[inew]
							jnew -= 1
						else:
							U[i * n + j] = 1 if i == j else 0
					inew -= 1
				else:
					for j in range(n):
						U[i * n + j] = 1 if i == j else 0

		return status

	def PMatUVRoot(self, P, t, n, U, V, Root):
		"""P(t) = U * exp{Root*t} * V"""
		smallp = 0

		# Initialize P with zeros
		for i in range(n * n):
			P[i] = 0

		for k in range(n):
			expt = math.exp(t * Root[k])
			for i in range(n):
				uexpt = U[i * n + k] * expt
				for j in range(n):
					P[i * n + j] += uexpt * V[k * n + j]

		for i in range(n * n):
			if P[i] < smallp:
				P[i] = 0

		return 0

	#function to find eigen solution of real symmetrical matrix A[n*n]
	def eigenRealSym(self, A, n, Root, work):
		status = 0
		self.HouseholderRealSym(A, n, Root, work)
		status = self.EigenTridagQLImplicit(Root, work, n, A)

		# Check for negative eigenvalues which shouldn't happen for transition matrices
		for i in range(n):
			if Root[i] < -1e-10:  # Significant negative eigenvalue
				return -1
			
		self.EigenSort(Root, A, n)
		return status

	#function for householder transformation to reduce matrix to tridiagonal form
	def HouseholderRealSym(self, a, n, d, e):
		for i in range(n - 1, 0, -1):
			m = i - 1
			h = scale = 0
			if m > 0:
				for k in range(m + 1):
					scale += abs(a[i * n + k])
				if scale == 0:
					e[i] = a[i * n + m]
				else:
					for k in range(m + 1):
						a[i * n + k] /= scale
						h += a[i * n + k] * a[i * n + k]
					f = a[i * n + m]
					g = -math.sqrt(h) if f >= 0 else math.sqrt(h)
					e[i] = scale * g
					h -= f * g
					a[i * n + m] = f - g
					f = 0
					for j in range(m + 1):
						a[j * n + i] = a[i * n + j] / h
						g = 0
						for k in range(j + 1):
							g += a[j * n + k] * a[i * n + k]
						for k in range(j + 1, m + 1):
							g += a[k * n + j] * a[i * n + k]
						e[j] = g / h
						f += e[j] * a[i * n + j]
					hh = f / (h * 2)
					for j in range(m + 1):
						f = a[i * n + j]
						e[j] = g = e[j] - hh * f
						for k in range(j + 1):
							a[j * n + k] -= (f * e[k] + g * a[i * n + k])
			else:
				e[i] = a[i * n + m]
			d[i] = h

		d[0] = e[0] = 0

		# Get eigenvectors
		for i in range(n):
			m = i - 1
			if d[i]:
				for j in range(m + 1):
					g = 0
					for k in range(m + 1):
						g += a[i * n + k] * a[k * n + j]
					for k in range(m + 1):
						a[k * n + j] -= g * a[k * n + i]
			d[i] = a[i * n + i]
			a[i * n + i] = 1
			for j in range(m + 1):
				a[j * n + i] = a[i * n + j] = 0

	#function to find eigen solution of tridiagonal matrix
	def EigenTridagQLImplicit(self, d, e, n, z):
		niter = 30
		status = 0

		for i in range(1, n):
			e[i - 1] = e[i]
		e[n - 1] = 0

		for j in range(n):
			iter = 0
			while True:
				m = j
				for m in range(j, n - 1):
					dd = abs(d[m]) + abs(d[m + 1])
					if abs(e[m]) + dd == dd:
						break

				if m != j:
					if iter == niter:
						status = -1
						break
					iter += 1

					g = (d[j + 1] - d[j]) / (2 * e[j])

					# r=pythag(g,1)
					aa = abs(g)
					if aa > 1:
						r = aa * math.sqrt(1 + 1 / (g * g))
					else:
						r = math.sqrt(1 + g * g)

					g = d[m] - d[j] + e[j] / (g + SIGN(r, g))
					s = c = 1
					p = 0

					for i in range(m - 1, j - 1, -1):
						f = s * e[i]
						b = c * e[i]

						# r=pythag(f,g)
						aa = abs(f)
						bb = abs(g)
						if aa > bb:
							bb /= aa
							r = aa * math.sqrt(1 + bb * bb)
						elif bb == 0:
							r = 0
						else:
							aa /= bb
							r = bb * math.sqrt(1 + aa * aa)

						e[i + 1] = r
						if r == 0:
							d[i + 1] -= p
							e[m] = 0
							break

						s = f / r
						c = g / r
						g = d[i + 1] - p
						r = (d[i] - g) * s + 2 * c * b
						p = s * r
						d[i + 1] = g + p
						g = c * r - b

						for k in range(n):
							f = z[k * n + i + 1]
							z[k * n + i + 1] = s * z[k * n + i] + c * f
							z[k * n + i] = c * z[k * n + i] - s * f

					if r == 0 and i >= j:
						continue
					d[j] -= p
					e[j] = g
					e[m] = 0
				else:
					break

		return status


	#function to sort eigenvalues and eigenvectors
	def EigenSort(self, d, U, n):
		for i in range(n - 1):
			k = i
			p = d[i]
			for j in range(i + 1, n):
				if d[j] >= p:
					p = d[k := j]
			if k != i:
				d[k] = d[i]
				d[i] = p
				for j in range(n):
					p = U[j * n + i]
					U[j * n + i] = U[j * n + k]
					U[j * n + k] = p

#function to calculate Ka/Ks using the MYN method
def calculate_myn_kaks(seq1, seq2,seq1_header,seq2_header):
	# Convert to uppercase
	seq1 = seq1.upper()
	seq2 = seq2.upper()

	# Validate input
	if len(seq1) != len(seq2):
		raise ValueError("Sequences must be of same length")

	if len(seq1) % 3 != 0:
		raise ValueError("Sequence length must be divisible by 3")

	# Check for valid nucleotides
	valid_nucleotides = set('ACGT')
	if not set(seq1).issubset(valid_nucleotides) or not set(seq2).issubset(valid_nucleotides):
		raise ValueError("Sequences should only contain A, C, G, T")

	# Create MYN instance
	myn = MYN()

	# Get kappa values
	myn.GetKappa(seq1, seq2)

	# Calculate Ka and Ks
	dS, dN, SEdS, SEdN, status = myn.DistanceYN00(seq1, seq2)

	# Calculate Ka/Ks ratio
	ka_ks_ratio = dN / dS if dS > 1e-9 else float('inf')
	S = myn.S
	N = myn.N
	Sd = myn.Sd
	Nd = myn.Nd
	return dN, dS, Nd, Sd, N, S

# function to align two protein sequences using MAFFT
def align_proteins(seq1, seq2, mafft,seq1_header,seq2_header):
	with tempfile.TemporaryDirectory() as tmpdir:
		in_fasta = os.path.join(tmpdir, "input.fasta")
		out_fasta = os.path.join(tmpdir, "aligned.fasta")

		# Write input protein sequences
		with open(in_fasta, "w") as f:
			f.write(">seq1\n" + seq1 + "\n")
			f.write(">seq2\n" + seq2 + "\n")

		# Run MAFFT
		try:
			subprocess.run([mafft, "--quiet", "--globalpair", "--maxiterate", "1000","--amino", in_fasta], stdout=open(out_fasta, "w"), stderr=subprocess.PIPE,check = False)
		except subprocess.CalledProcessError as e:
			raise RuntimeError(f"MAFFT alignment failed: {e}")

		# Read the aligned sequences
		aligned_seqs = []
		with open(out_fasta) as f:
			seq = ''
			for line in f:
				if line.startswith('>'):
					if seq:
						aligned_seqs.append(seq)
						seq = ''
				else:
					seq += line.strip()
			if seq:
				aligned_seqs.append(seq)

	return aligned_seqs[0], aligned_seqs[1]

# function to project gaps from protein alignment into CDS
def project_gaps_nogap_cds_aligner(protein_aligned_list, cds_seq_list):
	if len(protein_aligned_list) != len(cds_seq_list):
		raise ValueError("Number of protein and CDS sequences must match")

	if not protein_aligned_list:
		return []

	# All aligned proteins should have the same length
	alignment_length = len(protein_aligned_list[0])
	for i, prot in enumerate(protein_aligned_list):
		if len(prot) != alignment_length:
			raise ValueError(f"All aligned proteins must have same length. Sequence {i} has length {len(prot)}, expected {alignment_length}")

	# Identify positions without gaps in ANY sequence
	gap_free_positions = []
	for pos in range(alignment_length):
		has_gap = any(protein_seq[pos] == '-' for protein_seq in protein_aligned_list)
		if not has_gap:
			gap_free_positions.append(pos)

	# Project gap-free positions to CDS
	codon_aligned_list = []
	for i, (protein_aligned, cds_seq) in enumerate(zip(protein_aligned_list, cds_seq_list)):
		codon_aligned = ''
		cds_pos = 0

		for pos in range(alignment_length):
			if pos in gap_free_positions:
				# This position has no gaps in any sequence, so include the codon
				if cds_pos + 3 <= len(cds_seq):
					codon_aligned += cds_seq[cds_pos:cds_pos + 3]
				else:
					# Handle case where CDS is shorter than expected
					codon_aligned += 'NNN'  # or handle as appropriate

			# Always advance CDS position if this protein doesn't have a gap here
			if protein_aligned[pos] != '-':
				cds_pos += 3

		codon_aligned_list.append(codon_aligned)

	return codon_aligned_list

#Get amino acid from codon, return '*' for stop codons
def get_amino_acid(codon):
	return codon_table.get(codon.upper(), '*')

#Convert nucleotide character to integer (0=A, 1=C, 2=G, 3=T)
def convert_char_to_int(char):
	char_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
	return char_map.get(char.upper(), -1)

#Convert integer to nucleotide character
def convert_int_to_char(i):
	int_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
	return int_map.get(i, 'N')

#Count synonymous(S) and nonsynonymous(N) sites for a single codon
def get_codon_site(codon):
	if get_amino_acid(codon) == '*':
		return 0.0, 0.0

	syn = 0.0
	stop = 0

	# Only check positions 0 and 2 (first and third position in codon)
	for i in [0, 2]:
		for j in range(4):
			temp = list(codon.upper())
			if j != convert_char_to_int(temp[i]):
				temp[i] = convert_int_to_char(j)
				temp_codon = ''.join(temp)

				if get_amino_acid(temp_codon) == '*':
					stop += 1
				else:
					if get_amino_acid(temp_codon) == get_amino_acid(codon):
						syn += 1

	S = syn / 3.0
	N = 3 - stop / 3.0 - syn / 3.0

	return S, N

#Count synonymous(Sd) and nonsynonymous(Nd) differences between two codons
def get_codon_difference(codon1, codon2):
	if get_amino_acid(codon1) == '*' or get_amino_acid(codon2) == '*':
		return 0.0, 0.0, 0

	codon1 = codon1.upper()
	codon2 = codon2.upper()

	# Find positions that differ
	diff = []
	for i in range(3):
		if codon1[i] != codon2[i]:
			diff.append(i)

	num = len(diff)

	# Two codons are same
	if num == 0:
		return 0.0, 0.0, 0

	snp = num
	path = 1
	for i in range(1, num + 1):
		path *= i

	sd_temp = 0.0
	nd_temp = 0.0
	stop = 0

	# Only one difference between two codons
	if num == 1:
		if get_amino_acid(codon1) == get_amino_acid(codon2):
			sd_temp = 1.0
		else:
			nd_temp = 1.0

	# Two differences between two codons
	elif num == 2:
		for i in range(num):
			for j in range(num):
				if i != j:
					temp1 = list(codon1)
					temp1[diff[i]] = codon2[diff[i]]
					temp1_str = ''.join(temp1)

					if get_amino_acid(temp1_str) != '*':
						# codon1 <-> temp1
						if get_amino_acid(temp1_str) == get_amino_acid(codon1):
							sd_temp += 1
						else:
							nd_temp += 1

						# temp1 <-> codon2
						if get_amino_acid(temp1_str) == get_amino_acid(codon2):
							sd_temp += 1
						else:
							nd_temp += 1
					else:
						stop += 1

	# Three differences between two codons
	elif num == 3:
		for i in range(3):
			for j in range(3):
				for k in range(3):
					if (i != j) and (i != k) and (j != k):
						temp1 = list(codon1)
						temp1[diff[i]] = codon2[diff[i]]
						temp1_str = ''.join(temp1)

						temp2 = list(temp1_str)
						temp2[diff[j]] = codon2[diff[j]]
						temp2_str = ''.join(temp2)

						if get_amino_acid(temp1_str) != '*' and get_amino_acid(temp2_str) != '*':
							# codon1 <-> temp1
							if get_amino_acid(temp1_str) == get_amino_acid(codon1):
								sd_temp += 1
							else:
								nd_temp += 1

							# temp1 <-> temp2
							if get_amino_acid(temp2_str) == get_amino_acid(temp1_str):
								sd_temp += 1
							else:
								nd_temp += 1

							# temp2 <-> codon2
							if get_amino_acid(codon2) == get_amino_acid(temp2_str):
								sd_temp += 1
							else:
								nd_temp += 1
						else:
							stop += 1

	if path == stop:
		# All pathways are through stop codons
		if num == 2:
			Sd = 0.5
			Nd = 1.5
		else:
			Sd = 1.0
			Nd = 2.0
	else:
		Sd = sd_temp / (path - stop)
		Nd = nd_temp / (path - stop)

	return Sd, Nd, snp

#Process sequences to count sites and differences
def preprocess_sequences(seq1, seq2):
	S = 0.0
	N = 0.0
	Sd = 0.0
	Nd = 0.0
	snp_total = 0

	# Count sites and differences
	for i in range(0, len(seq1), 3):
		codon1 = seq1[i:i + 3]
		codon2 = seq2[i:i + 3]

		if len(codon1) == 3 and len(codon2) == 3:
			# Skip codons with gaps or ambiguous bases
			if '-' in codon1 or '-' in codon2 or 'N' in codon1 or 'N' in codon2:
				continue

			s1, n1 = get_codon_site(codon1)
			s2, n2 = get_codon_site(codon2)
			S += s1 + s2
			N += n1 + n2

			sd, nd, snp = get_codon_difference(codon1, codon2)
			Sd += sd
			Nd += nd
			snp_total += snp

	# Average sites from both sequences
	S /= 2.0
	N /= 2.0

	# Scale the sum of S+N to the length of sequence
	total_length = len(seq1)  # Should be same as len(seq2)
	if S + N > 0:
		y = total_length / (S + N)
		S *= y
		N *= y

	return S, N, Sd, Nd, snp_total


#Jukes & Cantor's one-parameter formula for correction
def kaks_formula(p):
	d = 1 - (4 * p) / 3
	if d <= 0.0:
		return float('inf')  # NA in C++ code
	else:
		d = math.log(d)
		if d > 0.0:
			return float('inf')  # NA in C++ code
		else:
			return (-3.0) * d / 4.0

#Calculate Ka and Ks using Modified YN00 (MYN) method
def calculate_ka_ks_MYN(codon_seq1, codon_seq2, seq1_header, seq2_header):
	if len(codon_seq1) != len(codon_seq2):
		raise ValueError("Sequences must be of equal length")

	if len(codon_seq1) % 3 != 0:
		raise ValueError("Sequence length must be divisible by 3")

	try:
		# Create MYN instance
		myn = MYN()

		# Get kappa values - wrap in try-catch for math domain errors
		try:
			myn.GetKappa(codon_seq1, codon_seq2)
		except (ValueError, ZeroDivisionError) as e:
			if "math domain error" in str(e) or isinstance(e, ZeroDivisionError):
				# Use default kappa values and continue
				myn.kappatc = myn.kappaag = myn.kappa = 2.0
				myn.KAPPA = [2.0, 2.0]
			else:
				raise e

		# Calculate Ka and Ks - wrap in try-catch
		try:
			dS, dN, SEdS, SEdN, status = myn.DistanceYN00(codon_seq1, codon_seq2)
		except (ValueError, ZeroDivisionError) as e:
			if "math domain error" in str(e) or isinstance(e, ZeroDivisionError):
				# Return NaN values for problematic pairs
				return float('nan'), float('nan'), 0, 0, 0, 0
			else:
				raise e

		S = myn.S
		N = myn.N
		Sd = myn.Sd
		Nd = myn.Nd
		return dN, dS, Nd, Sd, N, S

	except Exception as e:
		print(f"Error in calculate_ka_ks_MYN: {e}")
		import traceback
		traceback.print_exc()
		return float('nan'), float('nan'), 0, 0, 0, 0

#function to calculate Ka and Ks using Nei-Gojobori method
def calculate_ka_ks_nei_gojobori(codon_seq1, codon_seq2):
	if len(codon_seq1) != len(codon_seq2):
		raise ValueError("Sequences must be of equal length")

	if len(codon_seq1) % 3 != 0:
		raise ValueError("Sequence length must be divisible by 3")

	S, N, Sd, Nd, snp_total = preprocess_sequences(codon_seq1, codon_seq2)

	if S <= 0 or N <= 0:
		return float('nan'), float('nan'), Nd, Sd, N, S

	pS = Sd / S
	pN = Nd / N

	Ks = kaks_formula(pS)
	Ka = kaks_formula(pN)

	return Ka, Ks, Nd, Sd, N, S

# master function controlling the different Ka/Ks calculation algorithm methods' execution
def kaks_finder(pep1withstop,pep2withstop,cds1withstop,cds2withstop, mafft,seq1_header,seq2_header,method,logger):
	# Check length of PEP and corresponding CDS sequences
	# Check if protein 1 actually has a stop codon
	if pep1withstop.endswith('*'):
		# Protein has stop codon - remove it and corresponding CDS nucleotides
		pep1 = pep1withstop.rstrip('*')
		if len(cds1withstop) == len(pep1withstop) * 3:
			# Perfect match with stop codon
			cds1 = cds1withstop[:-3]
		else:
			# Length mismatch - use calculated length
			expected_cds1_length = len(pep1) * 3
			cds1 = cds1withstop[:expected_cds1_length]
	else:
		# Protein has no stop codon - don't remove anything
		pep1 = pep1withstop
		expected_cds1_length = len(pep1) * 3
		cds1 = cds1withstop[:expected_cds1_length]

	# Check if protein 2 actually has a stop codon
	if pep2withstop.endswith('*'):
		# Protein has stop codon - remove it and corresponding CDS nucleotides
		pep2 = pep2withstop.rstrip('*')
		if len(cds2withstop) == len(pep2withstop) * 3:
			# Perfect match with stop codon
			cds2 = cds2withstop[:-3]
		else:
			# Length mismatch - use calculated length
			expected_cds2_length = len(pep2) * 3
			cds2 = cds2withstop[:expected_cds2_length]
	else:
		# Protein has no stop codon - don't remove anything
		pep2 = pep2withstop
		expected_cds2_length = len(pep2) * 3
		cds2 = cds2withstop[:expected_cds2_length]

	# Step 1: Align protein sequences
	aligned_protein1, aligned_protein2 = align_proteins(pep1,pep2,mafft,seq1_header,seq2_header)

	# Step 2: Project gaps onto CDS and give CDS sequences without gaps like pal2nal nogap flag
	aligned_cds_list = project_gaps_nogap_cds_aligner([aligned_protein1, aligned_protein2], [cds1, cds2])
	aligned_cds1, aligned_cds2 = aligned_cds_list[0], aligned_cds_list[1]

	# Step 3: Calculate Ka and Ks
	if method == 'NG':
		ka, ks, nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites = calculate_ka_ks_nei_gojobori(aligned_cds1, aligned_cds2)
	else:
		ka, ks, nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites = calculate_ka_ks_MYN(aligned_cds1, aligned_cds2,seq1_header,seq2_header)

	# Step 4: Compute the ka/ks ratio
	if ks > 0 and ks != float('inf'):
		kaks_ratio = ka / ks
	else:
		kaks_ratio = float('inf') if ka > 0 else float('nan')
	return kaks_ratio, nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites

#second level function for kaks analysis processing the CDS, PEP sequences correctly
def kaks_analyzer(pair,pep_sequences, cds_sequences,mafft,logger,method):
	try:
		pep1withstop = pep_sequences[pair[0]]
		pep2withstop = pep_sequences[pair[1]]
		cds1withstop = cds_sequences[pair[0]]
		cds2withstop = cds_sequences[pair[1]]
		kaks, nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites = kaks_finder(pep1withstop, pep2withstop, cds1withstop, cds2withstop, mafft,str(pair[0]),str(pair[1]),method,logger)
		return (tuple(pair),kaks,nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites)
	except Exception as  e:
		logger.error(f"Worker crashed with error: {e}")
		logger.exception("The error is as follows:")
		return(tuple(pair), float('nan'),0,0,0,0)

#function to find weighted median
def find_weighted_median (kaks_list,weights_list):
	kaks_array = np.array(kaks_list)
	weights_array = np.array(weights_list)
	sorted_indices = np.argsort(kaks_array)
	values_sorted = kaks_array[sorted_indices]
	weights_sorted = weights_array[sorted_indices]
	cumulative_weight = np.cumsum(weights_sorted)
	total_weight = cumulative_weight[-1]
	cutoff = 0.5* total_weight
	return values_sorted[np.searchsorted(cumulative_weight, cutoff)]

#function to load duplicate groups
def load_duplicate_groups(file_paths, pepname, classify,ref_name):
	"""Load duplicate groups from files"""
	duplicate_groups = []
	for file_path in file_paths:
		if not os.path.exists(file_path):
			continue
		basename = os.path.basename(file_path)
		# Extract organism name based on file type
		if '_tandems' in basename:
			organism_name = basename.split('_tandems')[0]
		elif '_proximals' in basename:
			organism_name = basename.split('_proximals')[0]
		elif '_dispersed_duplicates' in basename:
			organism_name = basename.split('_dispersed_duplicates')[0]
		elif '_mixed_duplicates' in basename:
			organism_name = basename.split('_mixed_duplicates')[0]
		else:
			continue
		if organism_name == pepname:
			with open(file_path, 'r') as f:
				next(f)  # Skip the header line
				for line in f:
					parts = line.strip().split('\t')
					if len(parts) >= 2:
						if ref_name != 'NA':
							gene_list = [gene.strip() for gene in parts[1].strip().split(',')]
						else:
							gene_list = [gene.strip() for gene in parts[0].strip().split(',')]
						duplicate_groups.append(gene_list)
	return duplicate_groups

#function to make dictionaries from CDS, PEP files
def load_sequences_to_dic(filepath):
	"""Load sequences from FASTA file"""
	sequences = {}
	if filepath[-2:].lower() != 'gz':
		with open(filepath, 'r') as f:
			gene_id = None
			seq_lines = []
			for line in f:
				line = line.strip()
				if line.startswith(">"):
					if gene_id:
						sequences[gene_id] = "".join(seq_lines)
					gene_id = line[1:].strip().split()[0]
					seq_lines = []
				else:
					seq_lines.append(line)
			if gene_id:
				sequences[gene_id] = "".join(seq_lines)
	else:
		with gzip.open(filepath, 'rt') as f:
			gene_id = None
			seq_lines = []
			for line in f:
				line = line.strip()
				if line.startswith(">"):
					if gene_id:
						sequences[gene_id] = "".join(seq_lines)
					gene_id = line[1:].strip().split()[0]
					seq_lines = []
				else:
					seq_lines.append(line)
			if gene_id:
				sequences[gene_id] = "".join(seq_lines)
	return sequences

#first level worker function for kaks analysis that finally computes individual gene-level Ka/Ks values and their statistical significance
def collect_groups_for_kaks(cds,pep,pepname,tandems_dir,proximals_dir,dispersed_dir,mixed_dir,classify,logger,mafft,threads_per_organism,kaks_dir,ref_name,method,kaks_output_file,kaks_plot_file,dpi_no):
	all_files = []
	if classify == 'strict':
		all_files.extend(tandems_dir + proximals_dir + dispersed_dir + mixed_dir)
	elif classify == 'overlap':
		all_files.extend(tandems_dir + proximals_dir + dispersed_dir)
	# Collect duplicate groups
	duplicate_groups = load_duplicate_groups(all_files, pepname, classify,ref_name)
	# Load sequences onto dictionaries
	pep_sequences =load_sequences_to_dic(pep)
	cds_sequences =load_sequences_to_dic(cds)
	#obtaining a set of all the genes in the gene duplicates groups
	duplicate_genes_set = set()
	for gene_list in duplicate_groups:
		for gene in gene_list:
			duplicate_genes_set.add(gene)
	group_occurrence_normalizer = {}
	if classify == 'strict':
		pass
	elif classify == 'overlap':
		#making a dictionary where key is the gene and the value is inverse of number of duplicate groups in which it occurs
		for gene in duplicate_genes_set:
			count = 0
			for gene_list in duplicate_groups:
				if gene in gene_list:
					count+=1
			group_occurrence_normalizer[gene] = float(1/count)
	group_length_normalizer = {}
	kaks_pairs = []
	for gene_list in duplicate_groups:
		pairs=[]
		for i in range(len(gene_list)):
			for j in range(i+1,len(gene_list)):
				pair = (gene_list[i],gene_list[j])
				pairs.append(pair)
		group_pair_length = len(pairs)
		for each in pairs:
			group_length_normalizer[each]=float(1/group_pair_length)
			kaks_pairs.append(list(each))
	logger.info(f"Processing {len(kaks_pairs)} gene pairs for {pepname}")
	kaks_results = {}
	with concurrent.futures.ProcessPoolExecutor (max_workers = threads_per_organism) as executor:
		futures = {}
		for pair in kaks_pairs:
			future = executor.submit(kaks_analyzer, pair, pep_sequences, cds_sequences, mafft, logger,method)
			futures[future] = pair  # map future to the pair
		#collect results of inner parallelization with kaks_analyzer
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		for future in iterator:
			pair = futures[future]
			try:
				pair_tuple, kaks_value, nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites = future.result()
				kaks_results[pair_tuple] = (kaks_value, nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites)
			except Exception as e:
				logger.error(f"Error processing pair {pair}: {e}")
				kaks_results[tuple(pair)] = (float('nan'), 0,0,0,0)
	# Consolidate results gene wise for each organism
	gene_counts = defaultdict(lambda: {"N": 0, "S": 0, "N_sites": 0, "S_sites": 0})
	kaks_consolidator = {}
	fisher_counter = 0
	target_genes_set = duplicate_genes_set if classify == 'strict' else group_occurrence_normalizer.keys()
	for each in target_genes_set:
		specific_gene_containing_pairs = []
		kaks_list = []
		weights_list = []

		for pairs in kaks_pairs:
			if each in pairs:
				specific_gene_containing_pairs.append(pairs)

		for pair in specific_gene_containing_pairs:
			kaks_value, nonsyn_differences, syn_differences, total_nonsyn_sites, total_syn_sites = kaks_results[tuple(pair)]

			# Calculate weight based on classification method
			if classify == 'strict':
				weight_value = float(group_length_normalizer[tuple(pair)])
			elif classify == 'overlap':
				weight_value = float(group_length_normalizer[tuple(pair)]) * float(group_occurrence_normalizer[each])

			if not (math.isnan(kaks_value) or math.isinf(kaks_value)):
				kaks_list.append(kaks_value)
				weights_list.append(weight_value)
				# Accumulate substitution counts for this gene
				gene_counts[each]["N"] += nonsyn_differences
				gene_counts[each]["S"] += syn_differences
				gene_counts[each]["N_sites"] += total_nonsyn_sites
				gene_counts[each]["S_sites"] += total_syn_sites

		if len(kaks_list) > 2:
			representative_kaks = find_weighted_median(kaks_list, weights_list)
		elif len(kaks_list) == 2:
			representative_kaks = np.mean(kaks_list)
		elif len(kaks_list) == 1:
			representative_kaks = kaks_list[0]
		elif len(kaks_list) == 0:
			representative_kaks = float('nan')
		if not (math.isnan(representative_kaks) or math.isinf(representative_kaks)):
			# Build contingency table
			gc = gene_counts[each]
			observed_nonsyn = max(0, int(gc["N"]))
			expected_nonsyn = max(0, int(gc["N_sites"] - gc["N"]))
			observed_syn = max(0, int(gc["S"]))
			expected_syn = max(0, int(gc["S_sites"] - gc["S"]))
			if (observed_nonsyn + expected_nonsyn > 0) and (observed_syn + expected_syn > 0):
				contingency_table = [[observed_nonsyn, expected_nonsyn],
									 [observed_syn, expected_syn]]
				try:
					_, pval = fisher_exact(contingency_table)
					fisher_counter += 1
				except Exception as e:
					logger.warning(f"Fisher's exact test failed for gene {each}: {e}")
					pval = float('nan')
			else:
				pval = float('nan')
		else:
			pval = float('nan')
		kaks_consolidator[each] = [representative_kaks, pval]

	#Bonferroni correction and writing the kaks values to output summary file
	if fisher_counter > 0:
		for gene in kaks_consolidator.keys():
			if not math.isnan(kaks_consolidator[gene][1]):
				pval = float(kaks_consolidator[gene][1])
				adjusted_pval = min(pval * fisher_counter, 1.0)  # Bonferroni correction
				kaks_consolidator[gene][1] = adjusted_pval
	with open(kaks_output_file, 'w') as out:
		out.write(f"Gene\tka/ks\tFisher's exact test p-value\tSelection pressure\n")
		for gene in kaks_consolidator.keys():
			kaks_value, pval = kaks_consolidator[gene]
			if math.isnan(kaks_value):
				out.write(f"{gene}\tNA\tNA\tNot determined\n")
				continue
			if math.isnan(pval):
				selection_type = "Inconclusive - statistical support not determined"
			elif pval <= 0.05:
				if kaks_value > 1:
					selection_type = "Positive selection"
				elif kaks_value < 1:
					selection_type = "Negative selection"
				else:
					selection_type = "Neutral selection"
			else:
				selection_type = "Inconclusive - no statistical support"
			out.write(f"{gene}\t{kaks_value:.6f}\t{pval if not math.isnan(pval) else 'NA'}\t{selection_type}\n")
	#code for ka/ks plot
	df =pd.read_csv(kaks_output_file,sep='\t')
	#remove NA values
	df = df[df['ka/ks'] != 'NA'].copy()
	df['ka/ks'] = pd.to_numeric(df['ka/ks'])
	df['pval'] = pd.to_numeric(df["Fisher's exact test p-value"], errors='coerce')
	#only take ka/ks values less than or equal to 0.05
	df_sig = df[df['pval'] <= 0.05]
	# Filter extreme Ka/Ks values
	df_filtered = df_sig[(df_sig['ka/ks'] > 0) & (df_sig['ka/ks'] < 5)]
	# Plot histogram
	fig, ax = plt.subplots(figsize=(10, 6))
	ax.hist(df_filtered['ka/ks'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
	ax.set_xlabel('Ka/Ks (ω)', fontsize=12)
	ax.set_ylabel('Number of genes', fontsize=12)
	ax.set_title(f'Ka/Ks Distribution (n={len(df_filtered)} genes, p ≤ 0.05)', fontsize=14)
	ax.legend()
	plt.tight_layout()
	plt.savefig(kaks_plot_file,dpi=dpi_no)
	logger.info(f"Ka/Ks analysis completed for {pepname}.")

#master function to control parallelization of ka/ks calculations
def parallelize_kaks_calculations(cds_dir,pep_dir,tandems_dir,proximals_dir,dispersed_dir,mixed_dir,mafft,cores,classify,logger,kaks_dir,ref_name,method,org_type,dpi_no):
	if ref_name != 'NA':
		num_files = len(cds_dir)-1
	else:
		num_files = len(cds_dir)
	total_cores = cores
	max_parallel_organisms, threads_per_organism = distribute_cores(total_cores, num_files)
	logger.info(f"Detected {total_cores} CPU cores.")
	logger.info(f"Parallelizing ka/ks calculation processes at outer level using {max_parallel_organisms} processes...")
	logger.info(f"Each organism will use {threads_per_organism} threads for inner parallelization.")
	with concurrent.futures.ProcessPoolExecutor (max_workers = max_parallel_organisms) as executor:
		futures=[]
		for each in cds_dir:
			if org_type == 'eukaryote':
				cdsname = (os.path.basename(str(each))).split('_no_alt_trans')[0]
			elif org_type == 'prokaryote':
				cdsname = (os.path.basename(str(each)))
			for every in pep_dir:
				pepname = (os.path.basename(str(every))).split('_no_alt_trans')[0]
				if cdsname == pepname and cdsname!=ref_name and pepname!=ref_name:
						kaks_output_file = os.path.join(kaks_dir, str(pepname) + "_kaks_summary.tsv")
						kaks_plot_file=os.path.join(kaks_dir, str(pepname) + "_kaks_plot.png")
						if not os.path.isfile(kaks_output_file):
							futures.append(executor.submit(collect_groups_for_kaks, str(each), str(every), pepname,tandems_dir,proximals_dir,dispersed_dir,mixed_dir,classify,logger,mafft,threads_per_organism,kaks_dir,ref_name,method,kaks_output_file,kaks_plot_file,dpi_no))
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures),desc="Organisms processed for ka/ks calculation")
		for future in iterator:
			try:
				messages = future.result()
				if messages:
					for msg in messages:
						logger.info(msg)
			# Wait and capture exceptions
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise
	logger.info("All Ka/Ks calculations completed successfully.")

#function to assign orthologs to the genes in the duplicates groups and singleton lists, determine genetic distance between duplicates, and output the nature of the tandem, and proximal duplicate groups
#this also helps in finding if the given duplicate array or pair has one (or multiple) reference anchor gene(s) in the reference organism or not
def ref_anchor_finder_output_writer(fwd_best_hits_dic, orthologs_dic, gene_duplicates, duplicates_output_file, gene_pos_per_chr_q, gene_pos_per_chr_r, gene_pos_nos_q, coding_non_coding_query, synteny_cutoff, flank_number, best_side, self_pairs_set):
	gene_to_contig = {
		gene: contig
		for contig, genes in gene_pos_per_chr_q.items()
		for gene in genes
	}
	gene_pos_q_lookup = build_gene_to_position(gene_pos_per_chr_q)
	coding_nc_lookup = build_gene_to_position(coding_non_coding_query)
	with open(duplicates_output_file, 'w') as out:
		if not gene_pos_per_chr_r:
			out.write("Query" + "\t" + "Pairwise gene distance (bp)" + "\t" + "Actual intervening gene number" + "\t" + "Apparent intervening gene number" + "\t" + "Group confidence" + "\n")
		else:
			out.write("Reference"+ "\t" +"Query"+ "\t" + "Pairwise gene distance (bp)" + "\t"+"Actual intervening gene number"+"\t"+"Apparent intervening gene number"+"\t"+"Synteny information"+"\t"+"Nature of duplicate group"+"\t"+"Group confidence"+"\n")
		low_confidence = 0
		moderate_confidence = 0
		high_confidence = 0
		for sublist in gene_duplicates:
			#code to perform confidence scoring of gene duplicate groups/ clusters
			pairs = []
			hit_pairs = 0
			for i in range(len(sublist)):
				for j in range(i + 1, len(sublist)):
					pair = (sublist[i], sublist[j])
					pairs.append(pair)
					inverse_pair = (sublist[j], sublist[i])
					if pair in self_pairs_set or inverse_pair in self_pairs_set:
						hit_pairs += 1
			total_pairs = len(pairs)
			confidence = float(hit_pairs/total_pairs)
			if confidence <= 0.3:
				group_confidence = 'Low confidence group'
				low_confidence += 1
			elif 0.3 < confidence <= 0.5:
				group_confidence = 'Moderate confidence group'
				moderate_confidence += 1
			elif confidence > 0.5:
				group_confidence = 'High confidence group'
				high_confidence += 1

			pairwise_distances = []
			for i in range(len(sublist) - 1):
				g1, g2 = sublist[i], sublist[i + 1]
				dist = abs(gene_pos_nos_q[g1] - gene_pos_nos_q[g2])
				pairwise_distances.append(dist)
			actual_gene_nos = find_intervening_genes(coding_nc_lookup, sublist)
			apparent_gene_nos = find_intervening_genes(gene_pos_q_lookup, sublist)

			if not gene_pos_per_chr_r:
				out.write((", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + str(group_confidence) + '\n')
			else:
				anchors_original = []
				sublist_len = len(sublist)
				for each in sublist:
					anchors_original.append(orthologs_dic.get(str(each)))
				anchors = list(dict.fromkeys(anchors_original))
				if len(anchors) == 1 and None in anchors:  # case1: the duplicates are not present in the forward best hits and hence does not have a best hit in the reference
					out.write('-' + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Missing reference gene' + '\t' + 'De novo duplication' + '\t' + str(group_confidence) + '\n')
				elif len(anchors) == 1 and None not in anchors:  # case2: the duplicates all have one and the same best hit in the reference
					synteny = synteny_checker(anchors, sublist, gene_pos_per_chr_q, gene_pos_per_chr_r, synteny_cutoff,flank_number, best_side, fwd_best_hits_dic)
					if synteny == 'yes':
						out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Gene group expansion' + '\t' + str(group_confidence) + '\n')
					elif synteny == 'no':
						out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Non-syntenic' + '\t' + 'Gene group expansion' + '\t' + str(group_confidence) + '\n')
				elif len(anchors) > 1 and None in anchors:  # case3: the duplicates have multiple best hits in the reference while some of them do not have best hits
					None_count = anchors.count(None)
					anchors_copy = copy.deepcopy(anchors)
					anchors_copy.remove(None)
					anchors_len = len(anchors_copy)
					synteny = synteny_checker(anchors_copy, sublist, gene_pos_per_chr_q, gene_pos_per_chr_r, synteny_cutoff, flank_number, best_side, fwd_best_hits_dic)
					if synteny == 'yes':
						if anchors_len > 1:
							begin = anchors_copy[0]
							end = anchors_copy[-1]
							contig_names = list(gene_pos_per_chr_r.keys())
							for contig in contig_names:
								for gene in gene_pos_per_chr_r[contig]:
									if begin == str(gene):
										ind_begin = gene_pos_per_chr_r[contig].index(gene)
									elif end == str(gene):
										ind_end = gene_pos_per_chr_r[contig].index(gene)
							if ind_begin > ind_end:
								index = anchors.index(None)  # Find the index of None
								anchors[index] = '-'  # Replace None with '-'
								if anchors_len == (sublist_len - None_count):
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic; Gene inversion' + '\t' + 'Conserved gene group and de novo duplication' + '\t' + str(group_confidence) + '\n')
								elif anchors_len != (sublist_len - None_count):
									if anchors_len < (sublist_len - None_count):
										out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic; Gene inversion' + '\t' + 'Gene group expansion and de novo duplication' + '\t' + str(group_confidence) + '\n')
									else:
										out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic; Gene inversion' + '\t' + 'Gene group contraction and de novo duplication' + '\t' + str(group_confidence) + '\n')
							elif ind_begin < ind_end:
								index = anchors.index(None)  # Find the index of None
								anchors[index] = '-'  # Replace None with '-'
								if anchors_len == (sublist_len - None_count):
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Conserved gene group and de novo duplication' + '\t' + str(group_confidence) + '\n')
								elif anchors_len != (sublist_len - None_count):
									if anchors_len < (sublist_len - None_count):
										out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Gene group expansion and de novo duplication' + '\t' + str(group_confidence) + '\n')
									else:
										out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Gene group contraction and de novo duplication' + '\t' + str(group_confidence) + '\n')
						elif anchors_len == 1:
							index = anchors.index(None)  # Find the index of None
							anchors[index] = '-'  # Replace None with '-'
							if anchors_len == (sublist_len - None_count):
								out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Conserved gene group and de novo duplication' + '\t' + str(group_confidence) + '\n')
							elif anchors_len != (sublist_len - None_count):
								if anchors_len < (sublist_len - None_count):
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Gene group expansion and de novo duplication' + '\t' + str(group_confidence) + '\n')
								else:
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Gene group contraction and de novo duplication' + '\t' + str(group_confidence) + '\n')
					elif synteny == 'no':
						index = anchors.index(None)  # Find the index of None
						anchors[index] = '-'  # Replace None with '-'
						if anchors_len == (sublist_len - None_count):
							out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Non-syntenic' + '\t' + 'Conserved gene group and de novo duplication' + '\t' + str(group_confidence) + '\n')
						elif anchors_len != (sublist_len - None_count):
							if anchors_len < (sublist_len - None_count):
								out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Non-syntenic' + '\t' + 'Gene group expansion and de novo duplication' + '\t' + str(group_confidence) + '\n')
							else:
								out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Non-syntenic' + '\t' + 'Gene group contraction and de novo duplication' + '\t' + str(group_confidence) + '\n')
				elif len(anchors) > 1 and None not in anchors:  # case4: the duplicates have multiple best hits in the reference
					synteny = synteny_checker(anchors, sublist, gene_pos_per_chr_q, gene_pos_per_chr_r, synteny_cutoff,flank_number, best_side, fwd_best_hits_dic)
					if synteny == 'yes':
						begin = anchors[0]
						end = anchors[-1]
						contig_names = list(gene_pos_per_chr_r.keys())
						for contig in contig_names:
							for gene in gene_pos_per_chr_r[contig]:
								if begin == str(gene):
									ind_begin = gene_pos_per_chr_r[contig].index(gene)
								elif end == str(gene):
									ind_end = gene_pos_per_chr_r[contig].index(gene)
						if ind_begin > ind_end:
							if len(anchors) == sublist_len:
								out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic; Gene inversion' + '\t' + 'Conserved gene group' + '\t' + str(group_confidence) + '\n')
							elif len(anchors) != sublist_len:
								if len(anchors) < sublist_len:
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic; Gene inversion' + '\t' + 'Gene group expansion and conserved gene group' + '\t' + str(group_confidence) + '\n')
								else:
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic; Gene inversion' + '\t' + 'Gene group contraction and conserved gene group' + '\t' + str(group_confidence) + '\n')
						elif ind_begin < ind_end:
							if len(anchors) == sublist_len:
								out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Conserved gene group' + '\t' + str(group_confidence) + '\n')
							elif len(anchors) != sublist_len:
								if len(anchors) < sublist_len:
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Gene group expansion and conserved gene group' + '\t' + str(group_confidence) + '\n')
								else:
									out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Syntenic' + '\t' + 'Gene group contraction and conserved gene group' + '\t' + str(group_confidence) + '\n')
					elif synteny == 'no':
						if len(anchors) == sublist_len:
							out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Non-syntenic' + '\t' + 'Conserved gene group' + '\t' + str(group_confidence) + '\n')
						elif len(anchors) != sublist_len:
							if len(anchors) < sublist_len:
								out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Non-syntenic' + '\t' + 'Gene group expansion and conserved gene group' + '\t' + str(group_confidence) + '\n')
							else:
								out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + 'Non-syntenic' + '\t' + 'Gene group contraction and conserved gene group' + '\t' + str(group_confidence) + '\n')
	return low_confidence, moderate_confidence, high_confidence

#function to assign orthologs to the genes in the duplicates groups and singleton lists, determine genetic distance between duplicates, and output the nature of the dispersed duplicate groups
#Finding if the given duplicate array or pair has one (or multiple) reference anchor gene(s) in the reference organism or not and write the output file
def ref_anchor_finder_output_writer_dispersed_mixed(fwd_best_hits_dic, orthologs_dic, gene_duplicates, duplicates_output_file, gene_pos_per_chr_q, gene_pos_per_chr_r, gene_pos_nos_q, coding_non_coding_query, dups_type, self_pairs_set):
	# Build a reverse lookup dict: gene → contig
	gene_to_contig = {
		gene: contig
		for contig, genes in gene_pos_per_chr_q.items()
		for gene in genes
	}
	gene_pos_q_lookup = build_gene_to_position(gene_pos_per_chr_q)
	coding_nc_lookup = build_gene_to_position(coding_non_coding_query)
	with open(duplicates_output_file, 'w') as out:
		if dups_type == 'mixed':
			if not gene_pos_per_chr_r:
				out.write("Query" + "\t" + "Pairwise gene distance (bp)" + "\t" + "Actual intervening gene number" + "\t" + "Apparent intervening gene number" + "\t" + "Mixed group duplicate types"+ "\t" + "Group confidence" + "\n")
			else:
				out.write("Reference" + "\t" + "Query" + "\t" + "Pairwise gene distance (bp)" + "\t" + "Actual intervening gene number" + "\t" + "Apparent intervening gene number" + "\t" + "Mixed group duplicate types" + "\t" + "Group confidence" + "\n")
		else:
			if not gene_pos_per_chr_r:
				out.write("Query" + "\t" + "Pairwise gene distance (bp)" + "\t" + "Actual intervening gene number" + "\t" + "Apparent intervening gene number" + "\t" + "Group confidence" + "\n")
			else:
				out.write("Reference" + "\t" + "Query" + "\t" + "Pairwise gene distance (bp)" + "\t" + "Actual intervening gene number" + "\t" + "Apparent intervening gene number" + "\t" + "Group confidence" + "\n")
		low_confidence = 0
		moderate_confidence = 0
		high_confidence = 0
		for sublist in gene_duplicates:
			# code to perform confidence scoring of gene duplicate groups/ clusters
			pairs = []
			hit_pairs = 0
			for i in range(len(sublist)):
				for j in range(i + 1, len(sublist)):
					pair = (sublist[i], sublist[j])
					pairs.append(pair)
					inverse_pair = (sublist[j], sublist[i])
					if pair in self_pairs_set or inverse_pair in self_pairs_set:
						hit_pairs += 1
			total_pairs = len(pairs)
			confidence = float(hit_pairs / total_pairs)
			if confidence <= 0.3:
				group_confidence = 'Low confidence group'
				low_confidence += 1
			elif 0.3 < confidence <= 0.5:
				group_confidence = 'Moderate confidence group'
				moderate_confidence += 1
			elif confidence > 0.5:
				group_confidence = 'High confidence group'
				high_confidence += 1

			#code to check if two genes in the sublist are present on the same contig and then give their pairwise distance
			pairwise_distances = []
			for i in range(len(sublist) - 1):
				g1, g2 = sublist[i], sublist[i + 1]
				if gene_to_contig.get(g1) == gene_to_contig.get(g2):
					dist = abs(gene_pos_nos_q[g1] - gene_pos_nos_q[g2])
				else:
					dist = "different_contig"
				pairwise_distances.append(dist)
			actual_gene_nos = find_intervening_genes(coding_nc_lookup, sublist)
			apparent_gene_nos = find_intervening_genes(gene_pos_q_lookup, sublist)

			if not gene_pos_per_chr_r:
				if dups_type == 'mixed':
					comment = str(sublist[-1])
					sublist = sublist[:-1]
					out.write( (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + str(comment).strip('[').strip(']') + '\t' + str(group_confidence) + '\n')
				else:
					out.write((", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + str(group_confidence) + '\n')
			else:
				anchors_original = []
				sublist_len = len(sublist)
				if dups_type == 'mixed':
					comment = str(sublist[-1])
					sublist = sublist[:-1]
					for each in sublist:
						anchors_original.append(orthologs_dic.get(str(each)))
				else:
					for each in sublist:
						anchors_original.append(orthologs_dic.get(str(each)))
				if dups_type=='mixed':
					out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + str(comment).strip('[').strip(']') + '\t' + str(group_confidence) + '\n')
				else:
					out.write((", ".join(str(x) for x in anchors_original)) + '\t' + (", ".join(str(x) for x in sublist)) + '\t' + str(pairwise_distances).strip('[').strip(']') + '\t' + str(actual_gene_nos).strip('[').strip(']') + '\t' + str(apparent_gene_nos).strip('[').strip(']') + '\t' + str(group_confidence) + '\n')

	return low_confidence, moderate_confidence, high_confidence

#functions to parallelize the gene duplicates identification, classification, ka/ks analysis, singletons and output writing steps
def identify_classify_duplicates_singletons(self_pairs_set,component,pos_dic_query,proximity,self_dic,orgname,each,no_alt_trans_dir,evo_analysis,files,pos_dic_ref, pos_nos_query, coding_non_coding_query,synteny_cutoff, flank_number, best_side, mafft,kaks_dir,singleton_dir_final,base_name,remaining_singletons, classify,tandems_dir, proximals_dir, dispersed_dir,mixed_dir,orth_file,tandems_output,proximals_output,dispersed_dups_output,mixed_dups_output,singletons_final, frequency_plots_dir, dpi_no):
	# code lines to get the final gene duplicates outputs
	t_master_start = time.perf_counter()
	duplicates_master_list = sec_besthits_to_duplicates_list(str(component))
	t_master_end = time.perf_counter()
	duplicates_master_list_copy = copy.deepcopy(duplicates_master_list)
	t_tand_start = time.perf_counter()
	tandems, no_tandems_master_duplicates_list, tandem_set = tandems_grouping(duplicates_master_list, pos_dic_query)
	t_prox_start = time.perf_counter()
	proximals, no_tandems_proximals_mixed_duplicates_list, proximal_set = proximal_grouping(no_tandems_master_duplicates_list, pos_dic_query, proximity)
	t_reclass_start = time.perf_counter()
	prefinal_tandems, prefinal_proximals, final_no_tandems_proximals_dispersed_duplicates_list, tandem_proximal_set, no_tandems_proximals_mixed_duplicates_set = proximal_tandem_reclassification(tandems, proximals, no_tandems_proximals_mixed_duplicates_list, pos_dic_query, self_dic,duplicates_master_list_copy)
	t_dispersed_start = time.perf_counter()
	dispersed_duplicates = dispersed_dups_grouping(final_no_tandems_proximals_dispersed_duplicates_list)
	t_connection_start = time.perf_counter()
	# Code to do final regrouping of different duplicate types to maintain the relationships between different duplicates as well as to attribute correct classification to the genes
	tandems_final, proximals_final, dispersed_final, mixed_final = link_tandems_proximals_dispersed(prefinal_tandems,prefinal_proximals,dispersed_duplicates,classify)
	t_connection_end = time.perf_counter()
	# sorting the genes in the duplicates outputs according to their order in the contig
	t_tandsort_start = time.perf_counter()
	# Build gene order mapping once
	all_genes = [gene for genes in pos_dic_query.values() for gene in genes]
	gene_order = {gene: i for i, gene in enumerate(all_genes)}
	for subgroup in tandems_final:
		subgroup.sort(key=lambda gene: gene_order.get(gene, float('inf')))
	t_proxsort_start = time.perf_counter()

	for subgroup in proximals_final:
		subgroup.sort(key=lambda gene: gene_order.get(gene, float('inf')))
	t_prox_filter_start = time.perf_counter()
	if classify == 'strict':
		pass
	else:
		proximals_filtered = output_subset_filter(proximals_final)
		proximals_cleaned = short_distance_filter(proximals_filtered, pos_dic_query, proximity, 'proximal')
	t_prox_filter_end = time.perf_counter()
	for subgroup in dispersed_final:
		subgroup.sort(key=lambda gene: gene_order.get(gene, float('inf')))
	t_disp_sort_end = time.perf_counter()
	for subgroup in mixed_final:
		subgroup.sort(key=lambda gene: gene_order.get(gene, float('inf')))
	t_mix_sort_end = time.perf_counter()
	if classify == 'strict':
		pass
	else:
		dispersed_filtered = output_subset_filter(dispersed_final)
		dispersed_cleaned = short_distance_filter(dispersed_filtered, pos_dic_query, proximity, 'dispersed')
	t_disp_filter_end = time.perf_counter()
	final_tandem_set = set()
	final_proximal_set = set()
	final_dispersed_set = set()
	final_mixed_set = set()
	messages = []
	i = 0
	for subgroup in tandems_final:
		for every in subgroup:
			final_tandem_set.add(every)
			i += 1
	t_tandcount_end = time.perf_counter()
	j=0
	if classify == 'strict':
		for subgroup in proximals_final:
			for every in subgroup:
				final_proximal_set.add(every)
				j += 1
	else:
		for subgroup in proximals_cleaned:
			for every in subgroup:
				final_proximal_set.add(every)
		j = len(final_proximal_set)
	t_proxcount_end = time.perf_counter()

	k = 0
	if classify == 'strict':
		for subgroup in dispersed_final:
			for every in subgroup:
				final_dispersed_set.add(every)
				k += 1
	else:
		for subgroup in dispersed_cleaned:
			for every in subgroup:
				final_dispersed_set.add(every)
				k += 1
	t_dispcount_end = time.perf_counter()
	l = 0
	for subgroup in mixed_final:
		subgroup = subgroup[:-1]
		for every in subgroup:
			final_mixed_set.add(every)
			l += 1
	t_mixcount_end = time.perf_counter()

	total_duplicates_set = final_tandem_set | final_proximal_set | final_dispersed_set | final_mixed_set

	if files != 'NA':
		fwd_dic = fwd_best_hits_to_dic(str(files))
		orthologs_dic = orthologs_hits_to_dic(orth_file)
	else:
		fwd_dic = {}
		orthologs_dic = {}
	low_confidence_tand, moderate_confidence_tand, high_confidence_tand = ref_anchor_finder_output_writer(fwd_dic, orthologs_dic, tandems_final, tandems_output, pos_dic_query,pos_dic_ref, pos_nos_query, coding_non_coding_query,synteny_cutoff, flank_number, best_side,self_pairs_set)
	t_tandwrite_end = time.perf_counter()
	if classify == 'strict':
		low_confidence_prox, moderate_confidence_prox, high_confidence_prox = ref_anchor_finder_output_writer(fwd_dic, orthologs_dic, proximals_final, proximals_output, pos_dic_query,pos_dic_ref, pos_nos_query, coding_non_coding_query,synteny_cutoff, flank_number, best_side,self_pairs_set)
	else:
		low_confidence_prox, moderate_confidence_prox, high_confidence_prox = ref_anchor_finder_output_writer(fwd_dic, orthologs_dic, proximals_cleaned, proximals_output, pos_dic_query,pos_dic_ref, pos_nos_query, coding_non_coding_query,synteny_cutoff, flank_number, best_side,self_pairs_set)
	t_proxwrite_end = time.perf_counter()
	if classify == 'strict':
		low_confidence_disp, moderate_confidence_disp, high_confidence_disp = ref_anchor_finder_output_writer_dispersed_mixed(fwd_dic, orthologs_dic, dispersed_final,dispersed_dups_output, pos_dic_query,pos_dic_ref, pos_nos_query,coding_non_coding_query, 'dispersed',self_pairs_set)
	else:
		low_confidence_disp, moderate_confidence_disp, high_confidence_disp = ref_anchor_finder_output_writer_dispersed_mixed(fwd_dic, orthologs_dic, dispersed_cleaned,dispersed_dups_output, pos_dic_query,pos_dic_ref, pos_nos_query,coding_non_coding_query, 'dispersed',self_pairs_set)
	t_dispwrite_end = time.perf_counter()
	low_confidence_mix, moderate_confidence_mix, high_confidence_mix = ref_anchor_finder_output_writer_dispersed_mixed(fwd_dic, orthologs_dic, mixed_final, mixed_dups_output,pos_dic_query, pos_dic_ref, pos_nos_query,coding_non_coding_query, 'mixed',self_pairs_set)

	if classify == 'strict':
		tot_low_confidence_groups = low_confidence_tand + low_confidence_prox + low_confidence_disp + low_confidence_mix
		tot_moderate_confidence_groups = moderate_confidence_tand + moderate_confidence_prox + moderate_confidence_disp + moderate_confidence_mix
		tot_high_confidence_groups = high_confidence_tand + high_confidence_prox + high_confidence_disp + high_confidence_mix
	else:
		tot_low_confidence_groups = low_confidence_tand + low_confidence_prox + low_confidence_disp
		tot_moderate_confidence_groups = moderate_confidence_tand + moderate_confidence_prox + moderate_confidence_disp
		tot_high_confidence_groups = high_confidence_tand + high_confidence_prox + high_confidence_disp 

	t_mixwrite_end = time.perf_counter()
	singletons_set_final = set()
	t_singletons_start = time.perf_counter()
	with open(str(each), 'r') as f:
		temp_singletons = f.read().splitlines()
	singletons_master_checklist = temp_singletons + remaining_singletons
	t_singleton_start = time.perf_counter()
	with open(singletons_final, 'w') as out:
		if not pos_dic_ref:
			out.write("Singletons" + '\n')
			for element in singletons_master_checklist:
				if element not in total_duplicates_set:
					sublist = [element]
					out.write((", ".join(sublist)) + '\n')
		else:
			out.write("Reference" + '\t' + "Query" + '\t' + 'Synteny information' + '\t' + 'Nature of singleton' + '\n')
			for element in singletons_master_checklist:
				if element not in total_duplicates_set:
					singletons_set_final.add(element)
				if element not in total_duplicates_set:
					if element in fwd_dic.keys():
						sublist = [element]
						anchor = orthologs_dic[element]
						anchor_list = [anchor]
						synteny = synteny_checker(anchor_list, sublist, pos_dic_query, pos_dic_ref, synteny_cutoff,flank_number, best_side, fwd_dic)
						if synteny == 'yes':
							out.write((", ".join(anchor_list)) + '\t' + (", ".join(sublist)) + '\t' + 'Syntenic' + '\t' + 'Conserved ortholog' + '\n')
						elif synteny == 'no':
							out.write((", ".join(anchor_list)) + '\t' + (", ".join(sublist)) + '\t' + 'Non-syntenic' + '\t' + 'Migrated ortholog' + '\n')
					else:
						sublist = [element]
						out.write(str('-') + '\t' + (", ".join(sublist)) + '\t' + 'NA' + '\t' + 'Ortholog missing in reference' + '\n')
	with open(singletons_final, 'r') as f:
		singles = f.readlines()[1:]
	t_singleton_end = time.perf_counter()

	if classify == 'strict':
		messages.append(str(orgname) + "\t" + str(i) + "\t" + str(j)+"\t" + str(k)+"\t" + str(l)+"\t"+str(len(singles)) + "\t" + str(tot_low_confidence_groups) + "\t" +str(tot_moderate_confidence_groups) + "\t" + str(tot_high_confidence_groups))
		#code to plot absolute and relative frequency plots
		category_num_data={"Tandem": i,
		"Proximal": j,
		"Dispersed": k,
		"Mixed": l,
		"Singletons": singles}
		# Convert to pandas Series
		s = pd.Series(category_num_data)

		# Compute relative frequencies (%)
		relative = s / s.sum() * 100
		# Absolute frequency plot
		plt.figure(figsize=(8, 5))
		s.plot(kind="bar")
		plt.ylabel("Absolute gene frequency")
		plt.xlabel("Classification category")
		plt.xticks(rotation=0)
		plt.tight_layout()
		abs_frequency_plot = os.path.join(frequency_plots_dir,orgname+'_absolute_frequency_plot.png')
		plt.savefig(abs_frequency_plot,dpi=dpi_no)

		# Relative frequency plot
		plt.figure(figsize=(8, 5))
		relative.plot(kind="bar")
		plt.ylabel("Relative gene frequency (%)")
		plt.xlabel("Classification category")
		plt.xticks(rotation=0)
		plt.tight_layout()
		rel_frequency_plot = os.path.join(frequency_plots_dir, orgname + '_relative_frequency_plot.png')
		plt.savefig(rel_frequency_plot, dpi=dpi_no)
	else:
		messages.append(str(orgname) + "\t" + str(i) + "\t" + str(j)+"\t" + str(k)+"\t"+str(len(singles))+ "\t" + str(tot_low_confidence_groups) + "\t" +str(tot_moderate_confidence_groups) + "\t" + str(tot_high_confidence_groups))
		# code to plot absolute and relative frequency plots
		category_num_data = {"Tandem": i,
							 "Proximal": j,
							 "Dispersed": k,
							 "Singletons": singles}
		# Convert to pandas Series
		s = pd.Series(category_num_data)

		# Compute relative frequencies (%)
		relative = s / s.sum() * 100
		# Absolute frequency plot
		plt.figure(figsize=(8, 5))
		s.plot(kind="bar")
		plt.ylabel("Absolute gene frequency")
		plt.xlabel("Classification category")
		plt.xticks(rotation=0)
		plt.tight_layout()
		abs_frequency_plot = os.path.join(frequency_plots_dir, orgname + '_absolute_frequency_plot.png')
		plt.savefig(abs_frequency_plot, dpi=dpi_no)

		# Relative frequency plot
		plt.figure(figsize=(8, 5))
		relative.plot(kind="bar")
		plt.ylabel("Relative gene frequency (%)")
		plt.xlabel("Classification category")
		plt.xticks(rotation=0)
		plt.tight_layout()
		rel_frequency_plot = os.path.join(frequency_plots_dir, orgname + '_relative_frequency_plot.png')
		plt.savefig(rel_frequency_plot, dpi=dpi_no)
	return messages

#master function to parallelize gene duplicates classification and ortholog assignments in the presence of a reference organism
def parallelize_duplicates_processing(gff3_input_file,peptide_file_list,self_list,self_file_list,fwd_file_list,orthologs_file_list, tmp_singleton_files,proximity,no_alt_trans_dir,evo_analysis,pos_dic_ref, synteny_cutoff, flank_number, best_side, mafft,kaks_dir,singleton_dir_final,classify,cores,ref_name,tandems_dir, proximals_dir, dispersed_dir,mixed_dir,unclassified_dir_final,output_dir,logger,process_pseudos, gff_posdic_file_list, frequency_plots_dir, dpi_no):
	# Dynamically detect number of CPU cores
	total_cores = cores
	if ref_name != 'NA':
		num_fasta_files = len(peptide_file_list)-1
	else:
		num_fasta_files = len(peptide_file_list)

	max_workers =distribute_cores_single_level(total_cores, num_fasta_files)
	logger.info(f"Detected {total_cores} CPU cores.")
	logger.info(f"Parallelizing gene duplicates identification, classification and singleton file generation using {max_workers} processes...")
	dups_results = []
	summary_file = os.path.join(output_dir,'Summary.tsv')

	with ProcessPoolExecutor(max_workers=max_workers) as executor:
		futures = []
		# code to get the positional dictionary of the query organism and final grouping and writing of gene duplicates
		for every in gff3_input_file:
			orgname = get_basename(str(every))
			unclassified_file = os.path.join(unclassified_dir_final,orgname[0]+'_unclassified_genes.tsv')
			for each in peptide_file_list:
				pep_name = (os.path.basename(str(each))).split('_no_alt_trans')[0]
				if orgname[0] == pep_name and orgname[0] != ref_name:
					for gffdic_file in gff_posdic_file_list:
						gffdic_org_name = get_basename(str(gffdic_file))[0]
						if gffdic_org_name == pep_name:
							with open(gffdic_file, 'rb') as f:
								gff_pos_details = pickle.load(f)
								pos_dic_query, pos_nos_query, coding_non_coding_query = gff_pos_details[gffdic_org_name]
					for file in self_list:
						file_name = (os.path.basename(str(file))).split('_self')[0]
						if file_name == orgname[0]:
							self_dic,unclassified_list = self_to_dic_for_reclassification(str(file))
							self_pairs_set = self_to_set(str(file))
							if unclassified_list:
								with open(unclassified_file,'w')as out:
									for unclassified in unclassified_list:
										out.write(str(unclassified)+'\n')
							remaining_singletons = get_singletons_from_self_blast_pep(str(each), str(file))
					for component in self_file_list:
						self_name = (os.path.basename(str(component))).split('_self_sec_best_hits')[0]
						if self_name == orgname[0]:
							if len(fwd_file_list)!=0:
								for files in fwd_file_list:
									fwd_name = (os.path.basename(str(files))).split('_best_hits')[0]
									if fwd_name == orgname[0]:
										fwd_file = files
								for orth in orthologs_file_list:
									orth_name = (os.path.basename(str(orth))).split('_orthologs')[0]
									if orth_name == orgname[0]:
										orth_file = orth
							else:
								fwd_file = 'NA'
								orth_file = '' #dummy ortholog file variable
							for each in tmp_singleton_files:
								base_name = (os.path.basename(str(each))).split('_singletons')[0]
								if base_name == orgname[0]:
									tandems_output = os.path.join(tandems_dir, str(orgname[0]) + "_tandems.tsv")
									proximals_output = os.path.join(proximals_dir, str(orgname[0]) + "_proximals.tsv")
									dispersed_dups_output = os.path.join(dispersed_dir, str(orgname[0]) + "_dispersed_duplicates.tsv")
									if classify == 'strict':
										mixed_dups_output = os.path.join(mixed_dir, str(orgname[0]) + "_mixed_duplicates.tsv")
									else:
										mixed_dups_output = os.path.join(mixed_dir, str(orgname[0]) + "_duplicates_relationships.tsv")
									singletons_final = os.path.join(singleton_dir_final, base_name + "_singletons.tsv")
									files_to_be_checked = [tandems_output,proximals_output,dispersed_dups_output,mixed_dups_output,singletons_final]
									all_files_exist = all(os.path.isfile(f) for f in files_to_be_checked)
									if all_files_exist:
										pass
									else:
										futures.append(executor.submit(identify_classify_duplicates_singletons, self_pairs_set,str(component), pos_dic_query, proximity, self_dic,str(orgname[0]), each, no_alt_trans_dir, evo_analysis,fwd_file, pos_dic_ref, pos_nos_query,coding_non_coding_query, synteny_cutoff, flank_number,best_side, mafft, kaks_dir, singleton_dir_final, base_name,remaining_singletons, classify, tandems_dir, proximals_dir,dispersed_dir, mixed_dir,orth_file,tandems_output,proximals_output,dispersed_dups_output,mixed_dups_output,singletons_final, frequency_plots_dir, dpi_no))
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures), desc="Organisms processed in the gene duplicates identification, classification and singleton file generation step")
		# Add tqdm to monitor futures
		for future in iterator:
			try:
				dups_stats = future.result()
				dups_results.append(dups_stats)
			# Wait and capture exceptions
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise
	if dups_results:
		with open(summary_file,'w')as out:
			if classify == 'overlap':
				out.write("Organism" + "\t" + "Tandem duplicates" + "\t" + "Proximal duplicates" + "\t" + "Dispersed duplicates" + "\t" + "Singletons" + "\t" + "Low confidence groups" + "\t" + "Moderate confidence groups" + "\t" + "High confidence groups" + "\n")
				for each in dups_results:
					out.write('\t'.join(str(x) for x in each)+'\n')
			elif classify == 'strict':
				out.write("Organism" + "\t" + "Tandem duplicates" + "\t" + "Proximal duplicates" + "\t" + "Dispersed duplicates"+ "\t" + "Mixed duplicates" + "\t" + "Singletons" + "\t" + "Low confidence groups" + "\t" + "Moderate confidence groups" + "\t" + "High confidence groups" + "\n")
				for each in dups_results:
					out.write('\t'.join(str(x) for x in each)+'\n')

#function to merge the gene duplicates outputs across query species
def merge_files(folder_path,ref_name):
	file_list = sorted(glob.glob(os.path.join(folder_path, '*.tsv')))
	merged_df = None
	for file in file_list:
		if 'tandems' in file:
			orgname = (os.path.basename(str(file))).split('_tandems')[0]
		elif 'proximals' in file:
			orgname = (os.path.basename(str(file))).split('_proximals')[0]
		elif 'dispersed_duplicates' in file:
			orgname = (os.path.basename(str(file))).split('_dispersed_duplicates')[0]
		elif 'mixed_duplicates' in file:
			orgname = (os.path.basename(str(file))).split('_mixed_duplicates')[0]
		elif 'duplicates_relationships' in file:
			orgname = (os.path.basename(str(file))).split('_duplicates_relationships')[0]
		elif 'singletons' in file:
			orgname = (os.path.basename(str(file))).split('_singletons')[0]
		df = pd.read_csv(file, sep='\t', usecols=[0, 1], header=None, names=[ref_name, orgname])
		if merged_df is None:
			merged_df = df
		else:
			merged_df = pd.merge(merged_df, df, on=ref_name, how='outer')
	merged_df.fillna('--', inplace=True)
	return merged_df

#function to make lists of specific lines containing the user specified desired genes
def lists_for_specified_results(duplicates_dir,desired_genes_set,ref_name):
	results = defaultdict(dict)  # {ref_gene_group: {organism: target_genes_only}}
	all_organisms = set()  # To track organism (file) names
	merged_rows = []#master list to append all subsequent rows to be printed in the final specified_genes_results file
	file_list = sorted(glob.glob(os.path.join(duplicates_dir, '*.tsv')))
	for file in file_list:
		base_filename = (os.path.basename(str(file)))
		if 'tandems' in base_filename:
			orgname = base_filename.split('_tandems')[0]
		elif 'proximals' in base_filename:
			orgname = base_filename.split('_proximals')[0]
		elif 'dispersed_duplicates' in base_filename:
			orgname = base_filename.split('_dispersed_duplicates')[0]
		elif 'mixed_duplicates' in base_filename:
			orgname = base_filename.split('_mixed_duplicates')[0]
		elif 'duplicates_relationships' in base_filename:
			orgname = base_filename.split('_duplicates_relationships')[0]
		elif 'singletons' in base_filename:
			orgname = base_filename.split('_singletons')[0]
		all_organisms.add(orgname)
		with open(file, newline='') as f:
			reader = csv.reader(f, delimiter='\t')
			next(reader)  # Skip the header row
			for row in reader:
				ref_genes = [g.strip().strip("'") for g in row[0].split(',')]
				query_genes = row[1].strip()
				# Check for match between reference and desired genes
				if any(g in desired_genes_set for g in ref_genes):
					ref_key = ','.join([f"{g}" for g in ref_genes])  #keep single quotes
					if orgname in results[ref_key]:
						if results[ref_key][orgname] != '--':
							results[ref_key][orgname] += f", {query_genes}"
						else:
							results[ref_key][orgname] = query_genes
					else:
						results[ref_key][orgname] = query_genes

	#make file headers
	organism_headers = sorted(all_organisms)
	merged_rows.append([ref_name] + organism_headers + ['classification'])
	#storing the further data rows to the list
	for ref_key in sorted(results.keys()):
		row = [ref_key]
		for org in organism_headers:
			value = results[ref_key].get(org, '--')
			row.append(value)
		if 'Tandem' in duplicates_dir:
			row.append('Tandem duplication class')
			merged_rows.append(row)
		elif 'Proximal' in duplicates_dir:
			row.append('Proximal duplication class')
			merged_rows.append(row)
		elif 'Dispersed' in duplicates_dir:
			row.append('Dispersed duplication class')
			merged_rows.append(row)
		elif 'Mixed' in duplicates_dir:
			row.append('Mixed duplication class')
			merged_rows.append(row)
		elif 'Singletons' in duplicates_dir:
			row.append('Singleton')
			merged_rows.append(row)
	return merged_rows

#function to remove perceived pseudogenes from the TPM file and randomly pair the genes
def find_pseudogenes_random_pairing(exp_table,output,rand_genes,avg_exp_cutoff):
	if exp_table[-2:].lower() != 'gz':  # dealing with uncompressed expression table file
		with open(exp_table,'r') as f:
			gene_expression_full={}
			pseudos=[]
			gene_names = []
			line = f.readline()
			while line:
				if 'gene' not in line:
					parts = line.strip().split('\t')
					# Get the expression data (all columns except the first one)
					expression_values = list(map(float, parts[1:]))  # Convert to float
					#adding new filters to remove pseudogenes
					expression_values_copy=expression_values.copy()
					# Convert to a NumPy array for easier handling of NaNs
					expression_array = np.array(expression_values_copy, dtype=float)
					# Remove NaN values
					valid_values = expression_array[~np.isnan(expression_array)]  # Keeps only non-NaN values
					if len(valid_values) > 0:  # Ensure there is at least one valid value
						average_expression = np.mean(valid_values)
					# Check if the expression data is constant (all values are the same)
					if (average_expression > avg_exp_cutoff):  # If there are more than 1 unique value, it's not constant, also the gene should have significant non-zero expression across samples
						gene_names.append(str(parts[0]))
						gene_expression_full[parts[0]] = parts[1:]
					elif (average_expression <= avg_exp_cutoff):#collecting the perceived pseudogenes
						pseudos.append(str(parts[0]))
				line = f.readline()
	else:
		with gzip.open(exp_table,'rt') as f:
			gene_expression_full={}
			pseudos=[]
			gene_names = []
			line = f.readline()
			while line:
				if 'gene' not in line:
					parts = line.strip().split('\t')
					# Get the expression data (all columns except the first one)
					expression_values = list(map(float, parts[1:]))  # Convert to float
					# adding new filters to remove pseudogenes
					expression_values_copy = expression_values.copy()
					# Convert to a NumPy array for easier handling of NaNs
					expression_array = np.array(expression_values_copy, dtype=float)
					# Remove NaN values
					valid_values = expression_array[~np.isnan(expression_array)]  # Keeps only non-NaN values
					if len(valid_values) > 0:  # Ensure there is at least one valid value
						average_expression = np.mean(valid_values)
					# Check if the expression data is constant (all values are the same)
					if (average_expression > avg_exp_cutoff):  # the gene should have significant non-zero expression across samples
						gene_names.append(str(parts[0]))
						gene_expression_full[parts[0]] = parts[1:]
					elif (average_expression <= avg_exp_cutoff):  # collecting the perceived pseudogenes
						pseudos.append(str(parts[0]))
				line = f.readline()

	#random pairing of genes
	if len(gene_names) <= rand_genes:
		randoms = len(gene_names)
	elif len(gene_names) > rand_genes:
		randoms = rand_genes
	# Ensuring randoms is even for proper pairing and avoiding index issues
	if randoms % 2 != 0:
		randoms -= 1

	if randoms < 2:
		print("Warning: Not enough genes for random pairing. Expression divergence threshold determination not possible.")
		return [], {}, [], 0
	random_rows = random.sample(gene_names, randoms)
	random.shuffle(random_rows)
	# Ensure no self-pairing by checking the list
	while any(random_rows[i] == random_rows[i + 1] for i in range(0, len(random_rows), 2)):
		random.shuffle(random_rows)  # Re-shuffle if self-pairing is detected
	pairs = [(random_rows[i], random_rows[i + 1]) for i in range(0, len(random_rows), 2)]
	#writing perceived pseudogenes to the output file
	pseudos_output_file_path = os.path.join(output, "Perceived_pseudogenes.txt")
	with open(pseudos_output_file_path, 'w') as out:
		for each in pseudos:
			out.write(str(each)+'\n')
	#debug
	return pseudos, gene_expression_full, pairs, randoms

#finding the r value threshold for expression divergence
def find_r(pairs,gene_expression_full,plots, dpi_no):
	r_list=[]
	for pair in pairs:
		gene_expression = {}
		for each in pair:
			if str(each) in gene_expression_full:
				gene_expression[str(each)] = gene_expression_full[str(each)]
		df = pd.DataFrame(gene_expression)  # Create a DataFrame
		# Convert to numeric (if necessary) and handle non-numeric values
		df = df.apply(pd.to_numeric, errors='coerce')
		# calculate the column-wise NaN percentage for each gene in the given gene pair
		NaN_list = []
		for every in pair:
			nan = df[str(every)].isna().mean() * 100
			NaN_list.append(nan)
		if NaN_list[0] > 10 or NaN_list[1] > 10:
			df = df.dropna()  # Drop rows with NaNs in either gene
		else:
			# If NaN percentage < 10% for both genes, fill NaNs with the median of each column
			df[str(pair[0])] = df[str(pair[0])].fillna(df[str(pair[0])].median())
			df[str(pair[1])] = df[str(pair[1])].fillna(df[str(pair[1])].median())
		# Log-transform the data (adding a small constant to avoid log(0))
		df = np.log1p(df)  # log1p is equivalent to log(1 + x)
		gene1data = df.iloc[:, 0]  # gets the first column of the dataframe
		gene2data = df.iloc[:, 1]  # gets the second column of the dataframe
		r, pval = spearmanr(gene1data, gene2data)
		r_list.append(r)
	rthresh = np.percentile(r_list, 95)
	# Create the density distribution plot of the Spearman r values
	plt.figure(figsize=(8, 6))
	sns.kdeplot(r_list, fill=False, color="purple", label="Density Distribution")

	# Add a vertical dotted line at the threshold
	plt.axvline(x=rthresh, color='black', linestyle='--', label=f"Expression Divergence Threshold")

	# Add a label for the threshold line
	plt.text(rthresh + 0.5, 0.1, f"r: {rthresh}", color='black', fontsize=12)

	# Customize the plot
	plt.title("Density Distribution of Spearman correlation coefficients obtained from random gene pairing")
	plt.xlabel("Spearman r")
	plt.ylabel("Density")
	plt.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), borderaxespad=0.0)
	# Save the density distribution plot
	plt.savefig(os.path.join(plots, 'Correlation_cofficients_KDE_plot.png'), bbox_inches='tight', dpi=dpi_no)
	# Close the main plot figure to free up memory
	plt.close()
	return rthresh

#function to confert the count table expression file to a dataframe
def tpm_to_df (gene_duplicates_list, exp_table,plots):
	if exp_table[-2:].lower() != 'gz':  # dealing with uncompressed expression table file
		with open(exp_table,'r') as f:
			gene_expression = {}
			line = f.readline()
			while line:
				for each in gene_duplicates_list:
					if line.startswith(str(each)):
						parts = line.strip().split('\t')
						gene_expression[parts[0]] = parts[1:]
				line = f.readline()
		df = pd.DataFrame(gene_expression)# Create a DataFrame
	else:#dealing with compressed expression table file
		with gzip.open(exp_table, 'rt') as f:
			gene_expression = {}
			line = f.readline()
			while line:
				for each in gene_duplicates_list:
					if line.startswith(str(each)):
						parts = line.strip().split('\t')
						gene_expression[parts[0]] = parts[1:]
				line = f.readline()
		df = pd.DataFrame(gene_expression)# Create a DataFrame
	return df

# Preprocessing function to handle NaNs and log-transform data for pairwise comparisons
def preprocess_data_pairwise(df, x_gene, y_gene, pseudos):
	# Create a temporary DataFrame with the two genes being compared
	temp_df = df[[x_gene, y_gene]].copy()
	# Convert to numeric (if necessary) and handle non-numeric values
	temp_df = temp_df.apply(pd.to_numeric, errors='coerce')
	# Calculate the percentage of NaNs
	nan_x=temp_df[x_gene].isna().mean() * 100
	nan_y=temp_df[y_gene].isna().mean() * 100
	if nan_x > 10 or nan_y > 10:
		temp_df = temp_df.dropna()  # Drop rows with NaNs in either gene
	else:
		# If NaN percentage < 10% for both genes, fill NaNs with the median of each column
		temp_df[x_gene] = temp_df[x_gene].fillna(temp_df[x_gene].median())
		temp_df[y_gene] = temp_df[y_gene].fillna(temp_df[y_gene].median())
	# Log-transform the data (adding a small constant to avoid log(0))
	temp_df = np.log1p(temp_df)  # log1p is equivalent to log(1 + x)
	return temp_df

#function to plot gene expression as 2-D histogram
def plot_gene_expression(df,gene1, gene2, ax, correlation_info, pseudo_genes, temp_output_list):
	# Preprocess the data pairwise for the current comparison
	temp_df = preprocess_data_pairwise(df, gene1, gene2, pseudo_genes)
	x=temp_df[gene2]
	y=temp_df[gene1]

	# Plot the 2D histogram on the provided axis
	hist = ax.hist2d(x, y, norm=mcolors.LogNorm(), cmap='gnuplot2_r', bins=25)
	# Add a sample size text
	sample_size_text = f'n = {len(x)}'
	ax.text(0.5, 0.90, sample_size_text, transform=ax.transAxes, ha='center', va='center', fontsize=8)
	# Let Matplotlib automatically set the ticks
	ax.tick_params(axis='both', which='major', labelsize=8)  # Adjust tick label size
	# Add axis labels to subplots
	ax.set_xlabel(gene2, fontsize=8)
	ax.set_ylabel(gene1, fontsize=8)

	if gene1 not in pseudo_genes and gene2 not in pseudo_genes:
		# Perform Spearman correlation
		corr, p_value = spearmanr(x, y)
		# Add correlation and p-value to the list for later use
		correlation_info.append((gene1, gene2, corr, p_value))
	return hist  # Return the histogram object to add a colorbar later

#function to create matrix gene expression plots
def create_matrix_plot_do_spearman(df, gene_duplicates_list, pseudo_genes, plots, temp_output_list,dpi_no):
	# Number of genes
	n = len(gene_duplicates_list)

	# Create a list to store correlation information
	correlation_info = []

	# Dynamically adjust figure size
	fig_size = max(2 * n, 10)  # Ensure a minimum size for readability
	fig, axes = plt.subplots(n, n, figsize=(fig_size, fig_size), constrained_layout=False)

	# Loop through the gene list and plot only upper triangle (excluding diagonal)
	for i in range(n):
		for j in range(n):
			if i < j:
				gene1 = gene_duplicates_list[i]
				gene2 = gene_duplicates_list[j]

				# Plot pairwise expression in the corresponding axis (upper triangle of matrix)
				ax = axes[i, j]
				hist = plot_gene_expression(df, gene1, gene2, ax, correlation_info, pseudo_genes, temp_output_list)
				# Turn off axes for unused subplots
				axes[j, i].axis('off')

			elif i == j:
				# Turn off diagonal plots
				axes[i, j].axis('off')

	# Add gene names to the top row and left column
	for i, gene in enumerate(gene_duplicates_list):
		# Add gene names to the left side of the matrix (y-axis labels for rows)
		axes[i, 0].annotate(gene, xy=(0, 0.5), xycoords='axes fraction', fontsize=15, ha='right', va='center',rotation=0, annotation_clip=False)
		# Add gene names to the top of the matrix (x-axis labels for columns)
		axes[0, i].set_title(gene_duplicates_list[i], fontsize=15, pad=20, loc='center')
	# Add a single colorbar in the blank lower triangle space
	cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.02])  # Position at the bottom of the figure
	cbar = fig.colorbar(hist[3], cax=cbar_ax, orientation='horizontal')
	cbar.set_label(r'$\log_{10}$ density of points', fontsize=12)

	# Add a note near the colorbar
	fig.text(0.5, 0.15, "NOTE_1: 'n' represents the number of samples used to create the gene expression plot",fontsize=20, ha='center', va='center')
	fig.text(0.5, 0.1, "NOTE_2: The x and y axes show the gene expression in log(1+TPM)",fontsize=20, ha='center', va='center')

	# Adjust layout using subplots_adjust instead of tight_layout
	fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15, wspace=0.3, hspace=0.3)

	# Save the matrix figure plot
	plt.savefig(os.path.join(plots, 'Pairwise_gene_expression_matrix.png'), bbox_inches='tight', dpi=dpi_no)
	# Close the main plot figure to free up memory
	plt.close()
	return correlation_info  # Return the list of correlation info

# function to perform Spearman correlation tests alone
def calculate_spearman_correlation(gene_duplicates_list, df, pseudo_genes, temp_output_list):

	# Number of genes
	n = len(gene_duplicates_list)

	# Create a list to store correlation information
	correlation_info = []

	# Loop through the gene pairs (upper triangle, excluding diagonal)
	for i in range(n):
		for j in range(i + 1, n):
			gene1 = gene_duplicates_list[i]
			gene2 = gene_duplicates_list[j]
			if gene1 not in pseudo_genes and gene2 not in pseudo_genes:
				# Preprocess the data pairwise for the current comparison
				temp_df = preprocess_data_pairwise(df, gene1, gene2, pseudo_genes)
				x = temp_df[gene1]
				y = temp_df[gene2]
				# Perform Spearman correlation
				corr, p_value = spearmanr(x, y)
				# Store results in the list
				correlation_info.append((gene1, gene2, corr, p_value))
	return correlation_info


#set of fucntions to parallelize rthreshold and pseudogenes finding steps
#functions to parallelize the tpm file cleaning step
def rthreshold_pseudogenes_func(rand_genes, org_folder, org_name, exp_table,dpi_no,avg_exp_cutoff,rthresh_pseudo_file):
	pseudo_genes_dic={}
	r_threshold_dic ={}
	# function call to calculate the r threshold value for classifying into expression divergence
	pseudo_genes, gene_expression_dic, random_gene_pairs, random_gene_nos = find_pseudogenes_random_pairing(exp_table,org_folder,int(rand_genes),avg_exp_cutoff)
	r_threshold = find_r(random_gene_pairs, gene_expression_dic, org_folder, dpi_no)
	pseudo_genes_dic[org_name]=pseudo_genes
	r_threshold_dic[org_name]=r_threshold
	with open(rthresh_pseudo_file,'wb')as f:
		pickle.dump((pseudo_genes_dic,r_threshold_dic),f)

#master function to parallelize finding the r value for divergent expression determination
def parallelize_rthresh_pseudogenes (analyses,validated_tpm_input_file,rand_genes,dpi_no,cores,logger,avg_exp_cutoff,rthreshold_pseudogenes_dic_dir):
	# Dynamically detect number of CPU cores
	total_cores = cores
	num_exp_files = (len(validated_tpm_input_file))
	max_workers = distribute_cores_single_level(total_cores, num_exp_files)

	logger.info(f"Detected {total_cores} CPU cores.")
	logger.info(f"Parallelizing rthreshold and pseudogenes identification steps using {max_workers} processes...")
	with concurrent.futures.ProcessPoolExecutor (max_workers = max_workers) as executor:
		futures=[]
		for every in validated_tpm_input_file:
			exp_name = get_basename(str(every))
			org_folder = os.path.join(analyses, str(exp_name[0]))  # Creating the path for the folder for storing the data of particular organism
			os.makedirs(org_folder, exist_ok=True)
			org_name=exp_name[0]
			exp_table = every
			rthresh_pseudo_file = os.path.join(rthreshold_pseudogenes_dic_dir,str(org_name)+'.pkl')
			if not os.path.isfile(rthresh_pseudo_file):
				futures.append(executor.submit(rthreshold_pseudogenes_func, rand_genes, org_folder,org_name,exp_table,dpi_no,avg_exp_cutoff,rthresh_pseudo_file))
		results = []
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures), desc="Organisms processed for the rthreshold and pseudogenes finding step")

#functions to parallelize the specific gene duplicate analysis steps
def specific_sublist_analysis (dups,pseudo_genes,exp_table,plots_org,stats_org,dpi_no,r_threshold,specific_pseudos):
	# function call to convert the gene expression data of the specified gene duplicates into a dataframe
	duplicates_list_copy = copy.deepcopy(dups)
	for gene in dups:
		if gene in pseudo_genes:
			specific_pseudos.append(gene)
			duplicates_list_copy.remove(gene)
	gene_exp_df = tpm_to_df(duplicates_list_copy, exp_table, plots_org)
	temp_output_list = []
	if 1 < len(duplicates_list_copy) <= 10:
		# function call to make the Spearman correlation tests and the multi plot matrix figure
		correlation_possibility = 'yes'
		correlation_information = create_matrix_plot_do_spearman(gene_exp_df, duplicates_list_copy,pseudo_genes, plots_org, temp_output_list,dpi_no)
	elif len(duplicates_list_copy) > 10:
		# function call to make the Spearman correlation tests
		correlation_possibility = 'yes'
		correlation_information = calculate_spearman_correlation(duplicates_list_copy, gene_exp_df,pseudo_genes, temp_output_list)
	elif len(duplicates_list_copy) == 1:
		correlation_possibility = 'no'
		correlation_information = []
	elif len(duplicates_list_copy) == 0:
		correlation_possibility = 'no as all pseudogenes'
		correlation_information = []
	# Bonferroni correction (adjust p-values)
	num_tests = len(correlation_information)
	corrected_p_values = [min(p * num_tests, 1.0) for _, _, _, p in correlation_information]
	# Write the Bonferroni adjusted p-values and correlation coefficients as a stats table
	output_file_path = os.path.join(stats_org, str(dups[0]) + "_group_stats.tsv")
	with open(output_file_path, 'w') as out:
		if correlation_possibility == 'yes':
			for (gene1, gene2, corr, p_val), corrected_p in zip(correlation_information, corrected_p_values):
				if corr >= r_threshold:
					if 0.0 <= abs(corr) < 0.1:
						out.write(f"{gene1}\t{gene2}\t{float(p_val):.3e}\t{float(corr):.3f}\tNegligible correlation\n")
					elif 0.1 <= abs(corr) < 0.4:
						out.write(f"{gene1}\t{gene2}\t{float(p_val):.3e}\t{float(corr):.3f}\tWeak correlation\n")
					elif 0.4 <= abs(corr) < 0.7:
						out.write(f"{gene1}\t{gene2}\t{float(p_val):.3e}\t{float(corr):.3f}\tModerate correlation\n")
					elif 0.7 <= abs(corr) < 0.9:
						out.write(f"{gene1}\t{gene2}\t{float(p_val):.3e}\t{float(corr):.3f}\tStrong correlation\n")
					elif 0.9 <= abs(corr) <= 1.0:
						out.write(f"{gene1}\t{gene2}\t{float(p_val):.3e}\t{float(corr):.3f}\tVery strong correlation\n")
				elif corr < r_threshold:
					out.write(f"{gene1}\t{gene2}\t{float(p_val):.3e}\t{float(corr):.3f}\tPotential expression divergence\n")
				if len(specific_pseudos) > 1 and len(specific_pseudos) != len(dups) and len(specific_pseudos) != (len(dups) - 1):
					out.write(str(specific_pseudos).strip("[]'").replace("'", "") + ' are pseudogenes in the duplicates array.'+'\n')
				elif len(specific_pseudos) == 1 and len(specific_pseudos) != len(dups) and len(specific_pseudos) != (len(dups) - 1):
					out.write(str(specific_pseudos).strip("[]'").replace("'", "") + ' is a pseudogene in the duplicates array.'+'\n')
		elif correlation_possibility == 'no' and len(specific_pseudos) == 1:  # one pseudogene in a gene duplicate pair array
			out.write(str(specific_pseudos).strip("[]'").replace("'", "") + ' is a pseudogene in the duplicate pair. ' +str(duplicates_list_copy[0]) + ' is a non-pseudogene. Potential expression divergence. No correlation analysis possible'+'\n')
		elif correlation_possibility == 'no' and len(specific_pseudos) > 1:  # multiple pseudogenes
			out.write(str(specific_pseudos).strip("[]'").replace("'", "") + ' are pseudogenes.'+str(duplicates_list_copy[0]) + 'is the only non-pseudogene. Potential expression divergence. No correlation analysis possible.'+'\n')
		elif correlation_possibility == 'no as all pseudogenes' and len(specific_pseudos) == len(dups):  # all genes of the duplicates array are pseudogenes
			out.write(str(specific_pseudos).strip("[]'").replace("'", "") + ' are all pseudogenes. No correlation analysis possible.'+'\n')
	with open(output_file_path, "r") as f:
		lines = f.readlines()
	unique_lines = list(dict.fromkeys(lines))  # removes duplicates and preserves order
	with open(output_file_path, "w") as f:
		f.writelines(unique_lines)

#function to perform gene expression analysis of gene duplicates in an organism-wise manner in the case of specific gene uplicates expression analysis
def specific_organism_analysis (df,dpi_no,rand_genes,org_col,every,org_folder,classification_col,chunk_size,pseudo_genes,r_threshold,logger,max_inner_workers):
	plots_folder = os.path.join(org_folder,'Plots')  # Creating the path for the folder for storing the data of particular organism
	os.makedirs(plots_folder, exist_ok=True)  # Creates an organism specific folder in the specified path location
	stats_folder = os.path.join(org_folder,'Stats')  # Creating the path for the folder for storing the data of particular organism
	os.makedirs(stats_folder, exist_ok=True)  # Creates an organism specific folder in the specified path location
	exp_table = str(every)
	duplicates_lists_per_column = []
	for index, row in df.iterrows():
		gene_name = row[org_col]  # Get gene name
		classification = row[classification_col]  # Get classification
		# Skip empty values or "--" entries
		if gene_name != "--":
			if classification != 'Singleton' and 'No hits in queries' not in classification:
				duplicates_list = [gene.strip().strip("'") for gene in str(gene_name).strip().split(',')]
				duplicates_lists_per_column.append(duplicates_list)
	total_duplicates_lists = len(duplicates_lists_per_column)
	if total_duplicates_lists<=500:
		with concurrent.futures.ProcessPoolExecutor(max_workers=max_inner_workers) as duplicates_list_executor:
			futures = []
			for dups in duplicates_lists_per_column:
				specific_pseudos = []  # list definition to collect pseudogenes in a duplicates group
				plots_org = os.path.join(plots_folder, str(dups[0]) + "_group_plots")  # Creating the path for the folder for storing the stats of particular organism
				counter = 1
				while os.path.exists(plots_org):
					plots_org = os.path.join(plots_folder, f"{dups[0]}_group{counter}_plots")
					counter+=1 #to avoid folder name collisions and file overwriting
				os.makedirs(plots_org,exist_ok=True)  # Creates an organism specific stats folder in the specified path location
				stats_org = os.path.join(stats_folder, str(dups[0]) + "_group_stats")  # Creating the path for the folder for storing the stats of particular organism
				stats_counter=1 #to avoid folder name collisions and file overwriting
				while os.path.exists(stats_org):
					stats_org = os.path.join(stats_folder,f"{dups[0]}_group{stats_counter}_stats")
					stats_counter += 1
				os.makedirs(stats_org,exist_ok=True)  # Creates an organism specific stats folder in the specified path location
				futures.append(duplicates_list_executor.submit(specific_sublist_analysis,dups,pseudo_genes,exp_table,plots_org,stats_org,dpi_no,r_threshold,specific_pseudos))
			# Add tqdm for inner progress
			# Prepare the iterator
			iterator = concurrent.futures.as_completed(futures)
			if tqdm_available:
				iterator = tqdm(iterator, total=len(futures), desc=f"Processing duplicates groups ({org_col})")
			for future in iterator:
				try:
					messages = future.result()
					if messages:
						for msg in messages:
							logger.info(msg)
				# Wait and capture exceptions
				except Exception as e:
					logger.error(f"Worker crashed with error: {e}")
					logger.exception("The error is as follows:")
					raise

	else:
		num_chunks = ceil(total_duplicates_lists / chunk_size)
		for i in range(num_chunks):
			start = i * chunk_size
			end = min((i + 1) * chunk_size, total_duplicates_lists)
			batch = duplicates_lists_per_column[start:end]
			with concurrent.futures.ProcessPoolExecutor(max_workers=max_inner_workers) as batch_executor:
				batch_futures = []
				for dups in batch:
					specific_pseudos=[]
					plots_org = os.path.join(plots_folder, str(dups[0]) + "_group_plots")  # Creating the path for the folder for storing the stats of particular organism
					counter = 1
					while os.path.exists(plots_org):
						plots_org = os.path.join(plots_folder, f"{dups[0]}_group{counter}_plots")
						counter += 1  # to avoid folder name collisions and file overwriting
					os.makedirs(plots_org, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
					stats_org = os.path.join(stats_folder, str(dups[0]) + "_group_stats")  # Creating the path for the folder for storing the stats of particular organism
					stats_counter = 1  # to avoid folder name collisions and file overwriting
					while os.path.exists(stats_org):
						stats_org = os.path.join(stats_folder, f"{dups[0]}_group{stats_counter}_stats")
						stats_counter+=1
					os.makedirs(stats_org, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
					batch_futures.append(batch_executor.submit(specific_sublist_analysis,dups,pseudo_genes,exp_table,plots_org,stats_org,dpi_no,r_threshold,specific_pseudos))
				# Add tqdm for inner progress
				# Prepare the iterator
				iterator = concurrent.futures.as_completed(batch_futures)
				if tqdm_available:
					iterator = tqdm(iterator, total=len(batch_futures),desc=f" Batch processing in progress - Batch {i + 1}/{num_chunks} ({org_col})")
				for future in iterator:
					try:
						messages = future.result()
						if messages:
							for msg in messages:
								logger.info(msg)
					# Wait and capture exceptions
					except Exception as e:
						logger.error(f"Worker crashed with error: {e}")
						logger.exception("The error is as follows:")
						raise

#master function to parallelize gene expression plotting, and statistical analysis of specific gene duplicates organism-wise in the specific gene duplicates analysis mode
def specific_two_step_parallelization_master(validated_tpm_input_file,master_folder, desired_genes_output_file,dpi_no,rand_genes,rthreshold_pseudogenes_file_list,logger,cores):
	df = pd.read_csv(desired_genes_output_file, sep='\t')
	# Extract column names (headers)
	columns = df.columns.tolist()
	# Skip the first and last columns
	org_columns = columns[1:-1]  # Organism name columns
	org_columns_copy = copy.deepcopy(org_columns)
	classification_col = columns[-1]  # Last column (classification)
	num_organisms = len(org_columns)
	total_cores = cores
	logger.info(f"Detected {total_cores} CPU cores.")

	outer_cores, inner_cores = distribute_cores(total_cores, num_organisms)

	logger.info(f"Allocating {outer_cores} cores for organism processing.")
	logger.info(f"Allowing up to {inner_cores} inner workers per organism.")
	chunk_size = 100
	#code line to avoid processing the column of organism whose expresio file is inappropriate or not given
	exp_set = set()
	for every in validated_tpm_input_file:
		exp_name = get_basename(str(every))
		exp_set.add(exp_name[0])
	for org_col in org_columns:
		if org_col not in exp_set:
			org_columns_copy.remove(org_col)
	# Parallelize across organisms
	with concurrent.futures.ProcessPoolExecutor(max_workers=outer_cores) as org_executor:
		org_futures = []
		for org_col in org_columns_copy:
			for every in validated_tpm_input_file:
				exp_name = get_basename(str(every))
				if exp_name[0] == org_col:
					target_name = exp_name[0]
					exp_file = str(every)
			org_folder = os.path.join(master_folder, str(org_col))  # Creating the path for the folder for storing the data of particular organism
			os.makedirs(org_folder,exist_ok=True)  # Creates an organism specific folder in the specified path location
			for rthresh_pseudo_file in rthreshold_pseudogenes_file_list:
				rthresh_orgname = get_basename(str(rthresh_pseudo_file))[0]
				if rthresh_orgname == target_name:
					with open(rthresh_pseudo_file,'rb')as f:
						pseudo_dic, r_dic = pickle.load(f)
					pseudo_genes = pseudo_dic[target_name]
					r_threshold = r_dic[target_name]
			org_futures.append(org_executor.submit(specific_organism_analysis,df,dpi_no,rand_genes,org_col,exp_file,org_folder,classification_col,chunk_size,pseudo_genes,r_threshold,logger,max_inner_workers=inner_cores))
		# Prepare the iterator
		iterator = concurrent.futures.as_completed(org_futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(org_futures), desc="Organisms processed for specific duplicates analysis")
		for future in iterator:
			try:
				messages = future.result()
				if messages:
					for msg in messages:
						logger.info(msg)
			# Wait and capture exceptions
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise
	logger.info('Completed expression analysis of the specific duplicates.')

#set of functions to parallelize the general gene duplications analysis steps for each organism in the regular expression analysis mode
def ref_independent_organism_analysis(exp_table,org_folder,rand_genes,dpi_no,each,chunk_size,folder_type,pseudo_genes, r_threshold, org_name, specific_dup_list, ref_free_folder, orgname, singleton, logger, max_inner_workers):
	singleton_list=[]
	with open(singleton,'r')as f:
		f.readline()
		line=f.readline()
		while line:
			gene=line.strip()
			singleton_list.append(gene)
			line=f.readline()
	specific_dup_list_copy = copy.deepcopy(specific_dup_list)
	for genedup in specific_dup_list:
		if genedup in singleton_list:
			specific_dup_list_copy.remove(genedup)
			logger.info(f'{genedup} is a singleton. Expression analysis and statistical analysis not possible.')
	orgfolder = os.path.join(ref_free_folder,str(orgname))
	os.makedirs(orgfolder,exist_ok=True)
	plots = os.path.join(orgfolder,"Plots")  # Creating the path for the folder for storing the plots of particular organism
	os.makedirs(plots, exist_ok=True)  # Creates an organism specific plots folder in the specified path location

	stats = os.path.join(orgfolder,"Stats")  # Creating the path for the folder for storing the stats of particular organism
	os.makedirs(stats, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
	tot_dups = []
	with open(str(each), 'r') as f:
		f.readline()
		line=f.readline()
		while line:
			parts = line.strip().split('\t')
			duplicates_list = [gene.strip().strip("'") for gene in str(parts[0]).strip().split(',')]
			tot_dups.append(duplicates_list)
			line=f.readline()
	tot_dups_no = len(tot_dups)
	if tot_dups_no<=1000 and len(specific_dup_list_copy)!=0:
		with concurrent.futures.ProcessPoolExecutor(max_workers=max_inner_workers) as duplicates_list_executor:
			futures = []
			for gene_dup in specific_dup_list_copy:
				for dups in tot_dups:
					if gene_dup in dups:
						specific_pseudos = []  # list definition to collect pseudogenes in a duplicates group
						plots_org = os.path.join(plots, str(gene_dup) + "_group_plots")  # Creating the path for the folder for storing the stats of particular organism
						counter = 1
						while os.path.exists(plots_org):
							plots_org = os.path.join(plots, f"{gene_dup}_group{counter}_plots")
							counter+=1 #to avoid folder name collisions and file overwriting
						os.makedirs(plots_org,exist_ok=True)  # Creates an organism specific stats folder in the specified path location
						stats_org = os.path.join(stats, str(gene_dup) + "_group_stats")  # Creating the path for the folder for storing the stats of particular organism
						stats_counter=1 #to avoid folder name collisions and file overwriting
						while os.path.exists(stats_org):
							stats_org = os.path.join(stats,f"{gene_dup}_group{stats_counter}_stats")
							stats_counter+=1
						os.makedirs(stats_org,exist_ok=True)  # Creates an organism specific stats folder in the specified path location
						futures.append(duplicates_list_executor.submit(specific_sublist_analysis,dups,pseudo_genes,exp_table,plots_org,stats_org,dpi_no,r_threshold,specific_pseudos))
			# Add tqdm for inner progress
			# Prepare the iterator
			iterator = concurrent.futures.as_completed(futures)
			if tqdm_available:
				iterator = tqdm(iterator, total=len(futures),desc=f"Processing duplicates groups ({org_name})")
			# Wait and capture exceptions
			for future in iterator:
				try:
					result = future.result()
				except Exception as e:
					logger.error(f"Worker crashed with error: {e}")
					logger.exception("The error is as follows:")
					raise
	elif tot_dups_no>1000 and len(specific_dup_list_copy)!=0:
		num_chunks = ceil(tot_dups_no / chunk_size)
		for i in range(num_chunks):
			start = i * chunk_size
			end = min((i + 1) * chunk_size, tot_dups_no)
			batch = tot_dups[start:end]
			with concurrent.futures.ProcessPoolExecutor(max_workers=max_inner_workers) as batch_executor:
				batch_futures = []
				for gene_dup in specific_dup_list_copy:
					for dups in batch:
						if gene_dup in dups:
							specific_pseudos = []
							plots_org = os.path.join(plots, str(dups[0]) + "_group_plots")  # Creating the path for the folder for storing the stats of particular organism
							counter = 1
							while os.path.exists(plots_org):
								plots_org = os.path.join(plots, f"{dups[0]}_group{counter}_plots")
								counter += 1
							os.makedirs(plots_org, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
							stats_org = os.path.join(stats, str(dups[0]) + "_group_stats")  # Creating the path for the folder for storing the stats of particular organism
							stats_counter = 1
							while os.path.exists(stats_org):
								stats_org = os.path.join(stats, f"{dups[0]}_group{stats_counter}_stats")
								stats_counter+=1
							os.makedirs(stats_org, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
							batch_futures.append(batch_executor.submit(specific_sublist_analysis,dups,pseudo_genes,exp_table,plots_org,stats_org,dpi_no,r_threshold,specific_pseudos))
				# Add tqdm for inner progress
				# Prepare the iterator
				iterator = concurrent.futures.as_completed(batch_futures)
				if tqdm_available:
					iterator = tqdm(iterator, total=len(batch_futures), desc=f" Batch processing in progress - Batch {i+1}/{num_chunks} ({org_name})")
				# Wait and capture exceptions
				for future in iterator:
					try:
						result = future.result()
					except Exception as e:
						logger.error(f"Worker crashed with error: {e}")
						logger.exception("The error is as follows:")
						raise

#set of functions to parallelize the general gene duplications analysis steps for each organism in the regular expression analysis mode
def organism_analysis(exp_table,org_folder,rand_genes,dpi_no,each,chunk_size,folder_type,pseudo_genes, r_threshold, org_name, logger, ref_name,max_inner_workers):

	org_duplicates_folder = os.path.join(org_folder,folder_type+"_duplicates")  # Creating the path for the folder for storing the data of tandems of particular organism
	os.makedirs(org_duplicates_folder,exist_ok=True)  # Creates an organism specific tandem duplicates analysis folder in the specified path location

	plots = os.path.join(org_duplicates_folder,"Plots")  # Creating the path for the folder for storing the plots of particular organism
	os.makedirs(plots, exist_ok=True)  # Creates an organism specific plots folder in the specified path location

	stats = os.path.join(org_duplicates_folder,"Stats")  # Creating the path for the folder for storing the stats of particular organism
	os.makedirs(stats, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
	tot_dups = []
	if ref_name != 'NA':
		i=1
	elif ref_name =='NA':
		i=0
	with open(str(each), 'r') as f:
		f.readline()
		line=f.readline()
		while line:
			parts = line.strip().split('\t')
			duplicates_list = [gene.strip().strip("'") for gene in str(parts[i]).strip().split(',')]
			tot_dups.append(duplicates_list)
			line=f.readline()
	tot_dups_no = len(tot_dups)
	if tot_dups_no<=1000:
		with concurrent.futures.ProcessPoolExecutor(max_workers=max_inner_workers) as duplicates_list_executor:
			futures = []
			for dups in tot_dups:
				specific_pseudos = []  # list definition to collect pseudogenes in a duplicates group
				plots_org = os.path.join(plots, str(dups[0]) + "_group_plots")  # Creating the path for the folder for storing the stats of particular organism
				counter = 1
				while os.path.exists(plots_org):
					plots_org = os.path.join(plots, f"{dups[0]}_group{counter}_plots")
					counter+=1 #to avoid folder name collisions and file overwriting
				os.makedirs(plots_org,exist_ok=True)  # Creates an organism specific stats folder in the specified path location
				stats_org = os.path.join(stats, str(dups[0]) + "_group_stats")  # Creating the path for the folder for storing the stats of particular organism
				stats_counter=1 #to avoid folder name collisions and file overwriting
				while os.path.exists(stats_org):
					stats_org = os.path.join(stats,f"{dups[0]}_group{stats_counter}_stats")
					stats_counter+=1
				os.makedirs(stats_org,exist_ok=True)  # Creates an organism specific stats folder in the specified path location
				futures.append(duplicates_list_executor.submit(specific_sublist_analysis,dups,pseudo_genes,exp_table,plots_org,stats_org,dpi_no,r_threshold,specific_pseudos))
			# Add tqdm for inner progress
			# Prepare the iterator
			iterator = concurrent.futures.as_completed(futures)
			if tqdm_available:
				iterator = tqdm(iterator, total=len(futures),desc=f"Processing duplicates groups ({org_name})")
			# Wait and capture exceptions
			for future in iterator:
				try:
					result = future.result()
				except Exception as e:
					logger.error(f"Worker crashed with error: {e}")
					logger.exception("The error is as follows:")
					raise
	else:
		num_chunks = ceil(tot_dups_no / chunk_size)
		for i in range(num_chunks):
			start = i * chunk_size
			end = min((i + 1) * chunk_size, tot_dups_no)
			batch = tot_dups[start:end]
			with concurrent.futures.ProcessPoolExecutor(max_workers=max_inner_workers) as batch_executor:
				batch_futures = []
				for dups in batch:
					specific_pseudos = []
					plots_org = os.path.join(plots, str(dups[0]) + "_group_plots")  # Creating the path for the folder for storing the stats of particular organism
					counter = 1
					while os.path.exists(plots_org):
						plots_org = os.path.join(plots, f"{dups[0]}_group{counter}_plots")
						counter += 1
					os.makedirs(plots_org, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
					stats_org = os.path.join(stats, str(dups[0]) + "_group_stats")  # Creating the path for the folder for storing the stats of particular organism
					stats_counter = 1
					while os.path.exists(stats_org):
						stats_org = os.path.join(stats, f"{dups[0]}_group{stats_counter}_stats")
						stats_counter+=1
					os.makedirs(stats_org, exist_ok=True)  # Creates an organism specific stats folder in the specified path location
					batch_futures.append(batch_executor.submit(specific_sublist_analysis,dups,pseudo_genes,exp_table,plots_org,stats_org,dpi_no,r_threshold,specific_pseudos))
				# Add tqdm for inner progress
				# Prepare the iterator
				iterator = concurrent.futures.as_completed(batch_futures)
				if tqdm_available:
					iterator = tqdm(iterator, total=len(batch_futures), desc=f" Batch processing in progress - Batch {i+1}/{num_chunks} ({org_name})")
				# Wait and capture exceptions
				for future in iterator:
					try:
						result = future.result()
					except Exception as e:
						logger.error(f"Worker crashed with error: {e}")
						logger.exception("The error is as follows:")
						raise

#master function to parallelize gene expression plotting and statistical analysis
def two_step_parallelization_master(duplicates_analysis,validated_tpm_input_file,analyses,rand_genes,dpi_no,folder_type,rthreshold_pseudogenes_file_list,cores,logger,specific_dup_list,ref_free_folder,singleton_files,ref_name):
	total_cores = cores
	logger.info(f"Detected {total_cores} CPU cores.")
	num_organisms = len(validated_tpm_input_file)
	logger.info(f"Detected {num_organisms} to process.")
	outer_cores, inner_cores = distribute_cores(total_cores, num_organisms)

	logger.info(f"Allocating {outer_cores} cores for organism processing.")
	logger.info(f"Allowing up to {inner_cores} inner workers per organism.")
	chunk_size=100

	with concurrent.futures.ProcessPoolExecutor(max_workers=outer_cores) as org_executor:
		org_futures=[]
		for each in duplicates_analysis:
			orgname = get_basename(str(each))[0]
			if 'tandems' in orgname:
				final_orgname = str(orgname).strip().split('_tandems')
			elif 'proximals' in orgname:
				final_orgname = str(orgname).strip().split('_proximals')
			elif 'dispersed_duplicates' in orgname:
				final_orgname = str(orgname).strip().split('_dispersed_duplicates')
			elif 'mixed_duplicates' in orgname:
				final_orgname = str(orgname).strip().split('_mixed_duplicates')
			elif 'duplicates_relationships' in orgname:
				final_orgname = str(orgname).strip().split('_duplicates_relationships')
			for every in validated_tpm_input_file:
				exp_name = get_basename(str(every))
				if final_orgname[0] == exp_name[0]:
					exp_table = str(every)
			org_folder = os.path.join(analyses, str(final_orgname[0]))  # Creating the path for the folder for storing the data of particular organism
			os.makedirs(org_folder, exist_ok=True)  # Creates an organism specific folder in the specified path location
			for rthresh_pseudo_file in rthreshold_pseudogenes_file_list:
				rthresh_orgname = get_basename(str(rthresh_pseudo_file))[0]
				if rthresh_orgname == final_orgname[0]:
					with open(rthresh_pseudo_file, 'rb') as f:
						pseudo_dic, r_dic = pickle.load(f)
					pseudo_genes = pseudo_dic[final_orgname[0]]
					r_threshold = r_dic[final_orgname[0]]
			org_name = final_orgname[0]

			if specific_dup_list:
				for every_single in singleton_files:
					orgname_singleton = get_basename(str(every_single))[0]
					final_org_singleton_name = str(orgname_singleton).strip().split('_singletons')[0]
					if final_org_singleton_name==org_name:
						singleton=every_single
				org_futures.append(org_executor.submit(ref_independent_organism_analysis, exp_table,org_folder,rand_genes,dpi_no,each,chunk_size,folder_type,pseudo_genes, r_threshold, org_name, specific_dup_list, ref_free_folder,str(final_orgname[0]),singleton,logger,max_inner_workers=inner_cores))
			else:
				org_futures.append(org_executor.submit(organism_analysis, exp_table, org_folder, rand_genes, dpi_no, each,chunk_size, folder_type, pseudo_genes, r_threshold, org_name,logger,ref_name,max_inner_workers=inner_cores))
		# Prepare the iterator
		iterator = concurrent.futures.as_completed(org_futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(org_futures),desc="Organisms processed for duplicates analysis")
		# Wait and capture exceptions
		for future in iterator:
			try:
				result = future.result()
			except Exception as e:
				logger.error(f"Worker crashed with error: {e}")
				logger.exception("The error is as follows:")
				raise

#main master function that runs all the parts of this script
def main( arguments ):
	"""! @brief runs all parts of this script """
	if '--prokaryote' in arguments:
		org_type = 'prokaryote'
	else:
		org_type = 'eukaryote'
	if '--ref' in arguments:#reference organism name for duplication anlysis
		ref_name = arguments[arguments.index('--ref') + 1]#name of the reference organism
	else:
		ref_name = 'NA'

	output_dir = arguments[ arguments.index('--out')+1 ] #full path to output folder
	if output_dir[-1] != "/":
		output_dir += "/"
	if not os.path.exists(output_dir):
		os.makedirs( output_dir )

	frequency_plots_dir=os.path.join(output_dir, "Duplication_frequency_plots")#creating the path for the folder for storing frequency plots of classification category orgnaism-wise
	os.makedirs(frequency_plots_dir, exist_ok=True)#Creates Duplication_freuqncy_plots folder in the specified location

	norm_bit_score_plots_dir = os.path.join(output_dir, "Duplication_landscape_plots")#Creating the path for the folder for storing normalized bit score plots for each query species
	os.makedirs(norm_bit_score_plots_dir, exist_ok=True)#Creates a Duplication_landscape_plots folder in the specified path location

	singleton_dir_final = os.path.join(output_dir, "Singletons")#Creating the path for the folder for storing singleton genes for each query species
	os.makedirs(singleton_dir_final, exist_ok=True)#Creates a Singletons folder in the specified path location

	unclassified_dir_final = os.path.join(output_dir, "Unclassified")#Creating the path for the folder for storing unclassified genes for each query species
	os.makedirs(unclassified_dir_final, exist_ok=True)#Creates a Unclassified folder in the specified path location

	tandems_dir = os.path.join(output_dir,"Tandem_duplicates")  # Creating the path for the folder for storing tandem gene duplicates for each query species
	os.makedirs(tandems_dir, exist_ok=True)  # Creates a TANDEM_DUPLICATES folder in the specified path location

	proximals_dir = os.path.join(output_dir,"Proximal_duplicates")  # Creating the path for the folder for storing proximal gene duplicates for each query species
	os.makedirs(proximals_dir, exist_ok=True)  # Creates a PROXIMAL_DUPLICATES folder in the specified path location

	dispersed_dir = os.path.join(output_dir, "Dispersed_duplicates")  # Creating the path for the folder for storing dispersed gene duplicates for each query species
	os.makedirs(dispersed_dir, exist_ok=True)  # Creates a Dispersed_DUPLICATES folder in the specified path location

	tmp_dir = os.path.join(output_dir, "Tmp")#Creating the path for the folder for storing all auxiliary or intermediate files of the Dupylicate analysis
	os.makedirs(tmp_dir, exist_ok=True)#Creates a TMP folder in the specified path location

	validated_tpm_files = os.path.join(tmp_dir, "Cleaned_expression_files")
	os.makedirs(validated_tpm_files, exist_ok=True)

	busco_dir = os.path.join(tmp_dir,"BUSCO_results")  # Creating the path for the folder for storing BUSCO outputs for each query species
	os.makedirs(busco_dir, exist_ok=True)  # Creates a busco output folder in the specified path location

	host_cache_dir = os.path.join(os.path.dirname(busco_dir), "busco_cache")
	os.makedirs(host_cache_dir, exist_ok=True)

	orthologs_dir = os.path.join(output_dir,"Orthologs")  # Creating the path for the folder for storing orthologs result for each query species
	os.makedirs(orthologs_dir, exist_ok=True)  # Creates an orthologs folder in the specified path location

	busco_dir_final = os.path.join(output_dir,"BUSCO_QC")  # Creating the path for the folder for storing BUSCO qc result for each query species
	os.makedirs(busco_dir_final, exist_ok=True)  # Creates a busco qc folder in the specified path location

	kaks_dir = os.path.join(output_dir, "Ka_Ks_analysis")  # Creating the path for the folder for storing the kaks summary files of the Dupylicate analysis
	os.makedirs(kaks_dir, exist_ok=True)  # Creates a kaks folder in the specified path location

	singleton_dir = os.path.join(tmp_dir, "Singletons_Tmp")  # Creating the path for the folder for storing singleton genes for each query species in a TMP file
	os.makedirs(singleton_dir, exist_ok=True)  # Creates a SINGLETONS_TMP folder in the specified path location

	blast_dir = os.path.join(tmp_dir, "Aln_outputs") #Creating the path for the folder for storing all BLAST output files to make further processing simpler
	os.makedirs(blast_dir, exist_ok=True)  # Creates a TMP folder in the specified path location

	fwd_best_hits_dir = os.path.join(blast_dir, "Fwd_best_hits")#Creating the path for the folder for storing forward best hits from BLAST/Diamond/mmseqs2
	os.makedirs(fwd_best_hits_dir, exist_ok=True)#Creates the FORWARD_BEST_HITS folder in the specified location

	self_best_hits_dir = os.path.join(blast_dir, "Self_best_hits")#Creating the path for the folder for storing self best hits from BLAST/Diamond/mmseqs2
	os.makedirs(self_best_hits_dir, exist_ok=True)#Creates the SELF_BEST_HITS folder in the specified location

	self_hits_dir = os.path.join(blast_dir,"Self_hits")  # Creating the path for the folder for storing self best hits from BLAST/Diamond/mmseqs2
	os.makedirs(self_hits_dir, exist_ok=True)  # Creates the SELF_BEST_HITS folder in the specified location

	database_dir = os.path.join(blast_dir,"Databases")  # Creating the path for the folder for storing database files from BLAST/Diamond/mmseqs2
	os.makedirs(database_dir, exist_ok=True)  # Creates the database folder in the specified location

	cds_dir = os.path.join(tmp_dir,"CDS_FASTAs")  # Creating the path for storing CDS FASTA files
	os.makedirs(cds_dir,exist_ok=True)  # Creates the cds_dir folder in the specified path location

	pep_dir = os.path.join(tmp_dir,"PEP_FASTAs")  # Creating the path for storing PEP FASTA files
	os.makedirs(pep_dir,exist_ok=True)  # Creates the pep_dir folder in the specified path location

	no_alt_trans_dir = os.path.join(tmp_dir, "PEP_FASTAs_no_alt_trans") #Creating the path for storing no alternate transcripts peptide FASTA files
	os.makedirs(no_alt_trans_dir, exist_ok=True) #Creates the PEPTIDE_FASTAs_SANS_ALTERNATE_TRANSCRIPTS folder in the specified path location

	no_alt_trans_cds_dir = os.path.join(tmp_dir,"CDS_FASTAs_no_alt_trans")  # Creating the path for storing no alternate transcripts peptide FASTA files
	os.makedirs(no_alt_trans_cds_dir,exist_ok=True)  # Creates the PEPTIDE_FASTAs_SANS_ALTERNATE_TRANSCRIPTS folder in the specified path location

	clean_cds_dir = os.path.join(tmp_dir,"Clean_CDS")  # Creating the path for the folder for storing cleaned cds seqeunces that were used to give rise to the pep sequences in the translate function
	os.makedirs(clean_cds_dir, exist_ok=True)

	gff_dic_dir = os.path.join(tmp_dir,"GFF_position_dictionaries")
	os.makedirs(gff_dic_dir, exist_ok=True)

	attribute_dic_dir = os.path.join(tmp_dir, "Attribute_dictionaries")
	os.makedirs(attribute_dic_dir, exist_ok=True)

	max_bit_score_dic_dir = os.path.join(tmp_dir,'Max_bit_scores_dictionaries')
	os.makedirs(max_bit_score_dic_dir, exist_ok=True)

	busco_dic_dir = os.path.join(tmp_dir,"BUSCO_dictionaries")
	os.makedirs(busco_dic_dir, exist_ok=True)

	rthreshold_pseudogenes_dic_dir = os.path.join(tmp_dir,'rthreshold_pseudogenes_dictionaries')
	os.makedirs(rthreshold_pseudogenes_dic_dir, exist_ok=True)

	logfile = os.path.join(output_dir,'Dupylicate.log')

	#Create a logger
	logger = logging.getLogger("dupylicate_logger")
	logger.setLevel(logging.DEBUG)  # Capture all levels of logs

	#Create a file handler to write logs to a file
	file_handler = logging.FileHandler(logfile)
	file_handler.setLevel(logging.DEBUG)  # Write all logs to file

	#Create a stream handler to print logs to console
	console_handler = logging.StreamHandler()
	console_handler.setLevel(logging.DEBUG)  # Show only INFO and above in console

	#Create a common log format
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
	file_handler.setFormatter(formatter)
	console_handler.setFormatter(formatter)

	# Step 5: Add both handlers to the logger
	logger.addHandler(file_handler)
	logger.addHandler(console_handler)

	if '--pseudos_gff' in arguments:#option to choose inclusion or exclusion of pseuodgenes from the gff file
		process_pseudos = arguments[arguments.index('--pseudos_gff')+1]#<yes or no>
	else:
		process_pseudos = 'no'
	gff_config_parameters = {}
	if '--gff_config' in arguments:
		gff_config_file = arguments[arguments.index('--gff_config')+1]
		with open(gff_config_file,'r')as f:
			for line in f:
				parts=line.strip().split()
				if len(parts)==4:
					gff_config_parameters[parts[0]] = {
							'child_attribute': parts[1],
							'child_parent_linker': parts[2],
							'parent_attribute': parts[3]
						}
	else:
		gff_config_parameters['default'] = {
			'child_attribute': 'ID',
			'child_parent_linker': 'Parent',
			'parent_attribute': 'ID'
		}
	# Option for user to give full path to busco
	if '--busco' in arguments:
		busco_path = arguments[arguments.index('--busco')+1]# busco_docker or default busco
	else:
		busco_path = 'busco'

	# Option for user to specify BUSCO lineage per input organism
	if '--busco_lineage' in arguments:
		buscolineage=arguments[arguments.index('--busco_lineage')+1]
	else:
		buscolineage='auto'

	if '--busco_version' in arguments:
		busco_version = arguments[arguments.index('--busco_version')+1]
	else:
		busco_version = "v6.0.0"#the most recent version of BUSCO at the time of writing this script

	if '--container_version' in arguments:
		container_version = arguments[arguments.index('--container_version')+1]
	else:
		container_version = 'cv1'

	if '--docker_host_path' in arguments:
		host_path = arguments[arguments.index('--docker_host_path')+1]
	else:
		host_path = "NA"

	if '--docker_container_path' in arguments:
		container_path = arguments[arguments.index('--docker_container_path')+1]
	else:
		container_path = "NA"

	#option for user to choose between BLAST, Diamond and mmseqs2
	if '--seq_aligner' in arguments:
		tool = arguments[arguments.index('--seq_aligner')+1]#<blast | diamond | mmseqs2 >
	else:
		tool = 'diamond'

	#Option for user to give full path to blast
	if '--blast' in arguments:
		aligner = arguments[arguments.index('--blast')+1]
		aligner = os.path.join(aligner, '')#to ensure that / is automatically added after the path is mentioned by the user if it is not already there
		aligner = aligner.replace(os.sep, '/')#to ensure adding only forward slash for ubuntu and windows systems
		makeblastdb = os.path.join(aligner + 'makeblastdb')
		blastp = os.path.join(aligner + 'blastp')
	else:
		makeblastdb = 'makeblastdb'
		blastp = 'blastp'

	#Option for user to give full path to diamond
	if '--diamond' in arguments:
		aligner = arguments[arguments.index('--diamond')+1]
		aligner = os.path.join(aligner, '')#to ensure that / is automatically added after the path is mentioned by the user if it is not already there
		aligner = aligner.replace(os.sep, '/')#to ensure adding only forward slash for ubuntu and windows systems
		diamond = os.path.join(aligner + 'diamond')
	else:
		diamond = 'diamond'

	#Option for user to give full path to mmseqs2
	if '--mmseqs' in arguments:
		aligner = arguments[arguments.index('--mmseqs')+1]
		aligner = os.path.join(aligner, '')#to ensure that / is automatically added after the path is mentioned by the user if it is not already there
		aligner = aligner.replace(os.sep, '/')#to ensure adding only forward slash for ubuntu and windows systems
		mmseqs2 = os.path.join(aligner + 'mmseqs')
	else:
		mmseqs2 = 'mmseqs'

	if '--tree' in arguments:
		do_tree = arguments[arguments.index('--tree')+1]#yes or no for including phylogeny during ortholog finding when reference is given
	else:
		do_tree='yes'

	#Option for user to give full path to MAFFT
	if '--mafft' in arguments:
		mafft = arguments[arguments.index('--mafft')+1]#full path to mafft including the mafft.bat file
	else:
		mafft = 'mafft' #defaults to v7.526 of MAFFT - the most recent version of MAFFT while developing this script

	#Option for user to specify the occupancy for alignment cleaning of the MAFFT file
	if '--occupancy' in arguments:
		occupancy = arguments[arguments.index('--occupancy')+1]
	else:
		occupancy = 0.1

	#Option for user to give full path to FastTree
	if '--fasttree' in arguments:
		fasttree = arguments[arguments.index('--fasttree')+1]
	else:
		fasttree = 'FastTree'

	#Option for user to give full path to GeMoMa
	if '--gemoma' in arguments:
		GeMoMa = arguments[arguments.index('--gemoma')+1]#full path including the name of the GeMoMa jar file

	if '--cores' in arguments:#total number of cores specified to run the Dupylicate analysis
		cores = int(arguments[arguments.index('--cores')+1])
	else:
		cores = 4

	if '--evalue' in arguments:
		eval = arguments[arguments.index('--evalue')+1]#user-specified string input for e value
	else:
		eval = '1e-5'

	if '--hits' in arguments:
		number_of_hits = int(arguments[arguments.index('--hits')+1])#user specified integer value for taking the top n best forward hits
	else:
		number_of_hits = 10

	if '--ortho_candidates' in arguments:
		possible_orthologs = int(arguments[arguments.index('--ortho_candidates')+1])#user specified integer value for listing the potential ortholog candidates for low confidence ortholog cases
	else:
		possible_orthologs = 3

	if '--score' in arguments:
		normalized_bit_score = arguments[arguments.index('--score')+1] #auto for automatic threshold finding and <float number> for manual threshold finding
	else:
		normalized_bit_score = 'auto'

	if '--qc' in arguments:
		quality_check = arguments[arguments.index('--qc')+1]#yes or no
	else:
		quality_check = 'no'

	if '--scoreratio' in arguments:#ratio of forward blast bit score and self blast bit score of query to assess if the forward hit is valid
		score_ratio_cutoff = float( arguments[ arguments.index('--scoreratio')+1 ] )
	else:
		score_ratio_cutoff = 0.3

	if '--fwd_simcut' in arguments:
		fwd_similarity_cutoff = float( arguments[ arguments.index('--fwd_simcut')+1 ] )#percentage value for similarity cutoff in forward alignment hits
	else:
		fwd_similarity_cutoff = 40.0#value in percent

	if '--self_simcut' in arguments:
		self_similarity_cutoff = float( arguments[ arguments.index('--self_simcut')+1 ] )#percentage value for similarity cutoff in self alignment hits
	else:
		self_similarity_cutoff = 50.0#value in percent

	if '--proximity' in arguments:
		proximity=int(arguments[arguments.index('--proximity')+1])#user specified integer value for the number of intervening genes to detect proximal duplications
	else:
		proximity = 10

	if '--synteny_score' in arguments:
		synteny_cutoff = float(arguments[arguments.index('--synteny_score')+1])#user specified float value for synteny analysis
	else:
		synteny_cutoff = 0.5

	if '--flank' in arguments:
		flank_number=int(arguments[arguments.index('--flank')+1])#user specified integer value for flanking gene count in synteny analysis
	else:
		flank_number=5

	if '--side' in arguments:
		best_side=int(arguments[arguments.index('--side')+1])#user specified integer value for synteny support from either side of a flanking region
	else:
		best_side = 1

	if '--ka_ks' in arguments:
		evo_analysis = arguments[arguments.index('--ka_ks')+1]#yes or no to calculate ka, ks values
	else:
		evo_analysis = 'no'

	if '--ka_ks_method' in arguments:
		method = arguments[arguments.index('--ka_ks_method')+1] # <MYN or NG>
	else:
		method = 'NG'

	if '--duplicates_analysis' in arguments:# to choose or not choose further statistical analysis of identified gene duplicates
		analysis = arguments[arguments.index('--duplicates_analysis')+1]
	else:
		analysis = 'no'

	if '--specific_duplicates_analysis' in arguments:#to choose or not choose further statistical analysis of specified ref genes' gene duplicates
		specific_analysis = arguments[arguments.index('--specific_duplicates_analysis')+1]
	else:
		specific_analysis = 'no'

	if '--ref_free_specific_duplicates_analysis' in arguments:#user given list of specific sample genes with one transcript name per line, for analysing gene duplication outputs
		ref_free_specific_analysis = arguments[arguments.index('--ref_free_specific_duplicates_analysis')+1]
		with open(ref_free_specific_analysis, 'r') as f:
			ref_free_desired_genes = [line.strip() for line in f]
	else:
		ref_free_desired_genes = []

	if '--dpi' in arguments:#resolution value desired by the user for plots
		dpi_level = arguments[arguments.index('--dpi')+1]#low/moderate/high/very high
		if dpi_level == 'low':
			dpi_no = 75
		elif dpi_level == 'moderate':
			dpi_no = 150
		elif dpi_level == 'high':
			dpi_no = 300
		elif dpi_level == 'very high':
			dpi_no = 600
	else:
		dpi_level = 'moderate'
		dpi_no = 150

	if '--analyse_disperse_duplicates' in arguments:#to choose or not choose statistical analysis of mixed gene duplicates
		dispersed_analysis_choice = arguments[arguments.index('--analyse_disperse_duplicates')+1]
	else:
		dispersed_analysis_choice = 'no'

	if '--exp' in arguments:
		exp_folder = arguments[arguments.index('--exp') + 1] # getting the full path to the counts TPM folder/file(s) from the user
		if exp_folder[-1] == "/":
			exp_input_file = sorted(glob.glob(exp_folder + "*.txt") + glob.glob(exp_folder + "*.txt.gz")+ glob.glob(exp_folder + "*.tpms.txt.gz")+ glob.glob(exp_folder + "*.tpms.txt"))
			exp_input_file = sorted(set(exp_input_file))
		else:
			exp_input_file = [exp_folder]
	else:
		exp_input_file = []

	if '--avg_exp_cut' in arguments:
		avg_exp_cutoff = arguments[arguments.index('--avg_exp_cut')+1]#float value from the user to cutoff the average expression value across samples for classifying a gene as pseudogene
	else:
		avg_exp_cutoff = 1

	if '--genes' in arguments:
		rand_genes=int(arguments[arguments.index('--genes')+1])#user specified integer value of genes to be taken for random pairing to determine the r value threshold
	else:
		rand_genes = 10000

	if '--specific_genes' in arguments:
		txt_file = arguments[arguments.index('--specific_genes') + 1]#user given list of specific genes with one transcript name per line, for analysing gene duplication outputs
		with open(txt_file, 'r') as f:
			desired_genes = [line.strip() for line in f]
	else:
		desired_genes = []

	if '--mode' in arguments:
		classify = arguments[arguments.index('--mode')+1]#overlap or strict
	else:
		classify = 'overlap'

	if '--clean_up' in arguments:
		clean_remnants = arguments[arguments.index('--clean_up')+1]# yes/ no
	else:
		clean_remnants = 'yes'

	if classify == 'strict':
		mixed_dir = os.path.join(output_dir,"Mixed_duplicates")  # Creating the path for the folder for storing mixed gene duplicates for each query species
		os.makedirs(mixed_dir, exist_ok=True)  # Creates a MIXED_DUPLICATES folder in the specified path location
	else:
		mixed_dir = os.path.join(output_dir, "Duplicates_relationships")  # Creating the path for the folder for storing mixed gene duplicates for each query species
		os.makedirs(mixed_dir, exist_ok=True)  # Creates a MIXED_DUPLICATES folder in the specified path location

	analyses = os.path.join(output_dir,"Duplicates_analysis")  # Creating the path for the folder for storing the plots and stats
	os.makedirs(analyses, exist_ok=True)  # Creates a PLOTS folder in the specified path location

	t_start = time.perf_counter()
	logger.info('Welcome to Dupylicate!!! Your analysis has started!!!')

	# --- if GFF3 and FASTA files of reference and query organisms are provided--- #

	if '--gff' in arguments and '--fasta' in arguments and '--to_annotate' not in arguments:
		gff3_folder = arguments[arguments.index('--gff')+1] #full path to GFF3 file or folder containing GFF3 annotation files
		if gff3_folder[-1] == "/":
			gff3_input_file = sorted(glob.glob(gff3_folder + "*.gff") + glob.glob(gff3_folder + "*.gff3") + glob.glob(gff3_folder + "*.gff.gz") + glob.glob(gff3_folder + "*.gff3.gz"))
		else:
			gff3_input_file = [gff3_folder]
		fasta_folder = arguments[arguments.index('--fasta')+1] #full path to folder containing assembly FASTA files
		if fasta_folder[-1] == "/":
			fasta_input_file = sorted(glob.glob(fasta_folder + "*.fa") + glob.glob(fasta_folder + "*.fna") + glob.glob(fasta_folder + "*.fna.gz") + glob.glob(fasta_folder + "*.fa.gz") + glob.glob(fasta_folder + "*.fasta.gz") + glob.glob(fasta_folder + "*.fasta") + glob.glob(fasta_folder + "*.genome.fa") + glob.glob(fasta_folder + "*.genome.fna") + glob.glob(fasta_folder + "*.genome.fna.gz") + glob.glob(fasta_folder + "*.genome.fa.gz") + glob.glob(fasta_folder + "*.genome.fasta.gz") + glob.glob(fasta_folder + "*.genome.fasta"))
		else:
			fasta_input_file = [fasta_folder]

		#Code to get md5sums of input files and write run parameters to a json file
		md5sum_dic ={}
		for each in fasta_input_file:
			file_md5sum = calculate_md5sum(str(each))
			md5sum_dic[str(each)] = file_md5sum
		for every in gff3_input_file:
			file_md5sum = calculate_md5sum(str(every))
			md5sum_dic[str(every)] = file_md5sum
		md5sum_file_path = os.path.join(output_dir,'md5sum.tsv')
		with open (md5sum_file_path,'w')as out:
			for file_path, checksum in md5sum_dic.items():
				out.write(f"{file_path}\t{checksum}\n")
		parameters = parse_arguments(arguments)
		with open(os.path.join(parameters["out"], "run_parameters.json"), "w") as f:json.dump(parameters, f, indent=4)

		strict_start = False  # option not implemented yet
		strict_end = False  # option not implemented yet
		blast_pep = []
		blast_db = []

		validated_fasta_files = os.path.join(tmp_dir,"Validated_FASTA")
		os.makedirs(validated_fasta_files, exist_ok=True)

		error_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(error_files, exist_ok=True)

		warning_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(warning_files, exist_ok=True)
		gemoma=False
		#Calling functions to check and validate the GFF and FASTA files
		total_errors_final = parallel_clean_and_validate(fasta_input_file, gff3_input_file, validated_fasta_files, error_files, warning_files, 'fasta',cores,logger, gff_config_parameters,gemoma)
		for everything in total_errors_final:
			if len(everything) != 0:
				logger.error('There are errors in your input file(s). Dupylicate analysis exiting...')
				sys.exit()
		validated_fasta_input_file = glob.glob(os.path.join(validated_fasta_files, "*.fa")) + glob.glob(os.path.join(validated_fasta_files, "*.fa.gz")) + glob.glob(os.path.join(validated_fasta_files, "*.fasta.gz")) + glob.glob(os.path.join(validated_fasta_files, "*.fasta"))

		# --- Calling the functions relevant for getting CDS, PEP, No alternate transcripts PEP --- #
		pep_gene_name_errors = parallel_process_fasta_gff(validated_fasta_input_file, gff3_input_file, tmp_dir, no_alt_trans_dir,no_alt_trans_cds_dir, ref_name, blast_dir, database_dir, self_hits_dir,self_best_hits_dir, singleton_dir, strict_start, strict_end, blast_pep,blast_db, tool, makeblastdb, blastp, diamond, mmseqs2, eval,normalized_bit_score, evo_analysis, 'fasta only',cores,logger,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_type,cds_dir,pep_dir,process_pseudos, validated_tpm_files, exp_input_file,validated_fasta_files,clean_cds_dir,gff_dic_dir,gff_config_parameters,attribute_dic_dir)
		if len(pep_gene_name_errors)!=0:
			logger.error('Please check the protein coding gene names in the GFF file and PEP file. Exiting ...')
			sys.exit()
		if org_type == 'eukaryote':
			peptides_folder_path = Path(no_alt_trans_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		elif org_type == 'prokaryote':
			peptides_folder_path = Path(pep_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		parallel_process_self_blast(ref_name, cores, peptides_file_list, database_dir, tmp_dir, tool, makeblastdb, blastp,diamond, mmseqs2, logger, self_hits_dir, self_best_hits_dir, singleton_dir, eval,normalized_bit_score, max_bit_score_dic_dir, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff)
		# code to do BUSCO based QC
		if normalized_bit_score == 'auto' or quality_check == 'yes':
			busco_cmd = detect_busco(busco_path)
			if busco_cmd:
				busco_qc_result = os.path.join(busco_dir_final, 'BUSCO_QC.tsv')
				busco_db_dir = os.path.join(busco_dir, 'busco_databases')
				total_cores = cores
				num_files = len(peptides_file_list)
				busco_threads_per_organism = total_cores
				logger.info(f"Serializing {num_files} organism-level processes for BUSCO analysis")
				logger.info(f"{busco_threads_per_organism} threads inside each organism for BUSCO analysis")
				busco_path_generated = "NA"
				container_workdir = busco_dir
				container_cache_dir = host_cache_dir

				if busco_path == 'busco_docker':
					container_db_dir = busco_db_dir
					docker_image = f"ezlabgva/busco:{busco_version}_{container_version}"
					# writing bash script for executing busco installed with docker
					bash_script = f"""#!/bin/bash
							# BUSCO Docker wrapper script with proper cleanup

							# Container name based on process ID and timestamp for uniqueness
							CONTAINER_NAME="busco_$$_$(date +%s)_$(shuf -i 1000-9999 -n 1)"

							# Cleanup function
							cleanup() {{
								echo "Cleaning up Docker container: $CONTAINER_NAME" >&2
								docker stop "$CONTAINER_NAME" 2>/dev/null || true
								docker rm "$CONTAINER_NAME" 2>/dev/null || true
								exit $1
							}}

							# Set up signal handlers for proper cleanup
							trap 'cleanup 130' INT    # Ctrl+C (SIGINT)
							trap 'cleanup 143' TERM   # Termination (SIGTERM)
							trap 'cleanup 1' EXIT     # Any exit
							trap 'cleanup 1' ERR      # Any error

							# Check if Docker is available
							if ! command -v docker &>/dev/null; then
								echo "Error: Docker is not installed or not in PATH." >&2
								exit 1
							fi

							# Run BUSCO in Docker with automatic cleanup
							docker run --rm \\
								--name "$CONTAINER_NAME" \\
								-u $(id -u) \\
								-v "{host_path}:{container_path}" \\
								-v "{host_cache_dir}:{container_cache_dir}" \\
								-v "{busco_db_dir}:{container_db_dir}" \\
								-e XDG_CONFIG_HOME="{container_cache_dir}" \\
								-w "{container_workdir}" \\
								{docker_image} \\
								busco "$@"

							# Capture exit code and exit cleanly
							EXIT_CODE=$?
							exit $EXIT_CODE
							"""
					output_filename = os.path.join(busco_dir_final, "run_busco_docker.sh")
					with open(output_filename, "w") as f:
						f.write(bash_script)
					os.chmod(output_filename, 0o755)
					busco_path_generated = output_filename
					# Pre-download for Docker
					if not os.path.isdir(busco_db_dir):
						pre_download_databases(busco_path_generated, org_type, busco_dir, busco_db_dir, logger)
				else:
					if not os.path.isdir(busco_db_dir):
						# Pre-download for normal BUSCO
						pre_download_databases(busco_path, org_type, busco_dir, busco_db_dir, logger)
				logger.info("starting actual BUSCO runs")
				ploidy_results = []
				busco_single_copy_lists = {}
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir,str(orgname)+'.pkl')
					if not os.path.isfile(busco_pickle_file):
						if orgname != ref_name:
							if busco_path != 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path, busco_dir, busco_dir_final,str(busco_threads_per_organism), org_type, host_cache_dir,busco_db_dir, busco_pickle_file,logger,buscolineage)
							elif busco_path == 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path_generated, busco_dir,busco_dir_final, str(busco_threads_per_organism), org_type,container_cache_dir, container_db_dir, busco_pickle_file,logger,buscolineage)
				busco_dic_file_path = Path(busco_dic_dir)
				busco_dic_file_list = list(busco_dic_file_path.iterdir())
				with open(busco_qc_result, 'w') as out:
					out.write('Organism' + '\t' + 'Pseudo ploidy number' + '\t' + 'BUSCO Completeness (%)' + '\t' + 'BUSCO Duplication (%)' + '\n')
					for busco_dic_file in busco_dic_file_list:
						with open(busco_dic_file,'rb') as f:
							busco_dic = pickle.load(f)
							if 'feff' in busco_dic:
								out.write(str(busco_dic['feff']) + '\n')
				logger.info(f"BUSCO QC check completed.")
			else:
				logger.info("BUSCO not found. No QC check possible. Using default method to find the threshold to segregate singletons and duplicated genes in the query organism(s)")
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
					busco_single_copy_lists = {}
					with open (busco_pickle_file,'wb') as f:
						pickle.dump(busco_single_copy_lists,f)

		else:
			for fasta_file in peptides_file_list:
				name = get_basename(fasta_file)
				orgname = str(name[0]).split('_no_alt_trans')[0]
				busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
				busco_single_copy_lists = {}
				with open(busco_pickle_file, 'wb') as f:
					pickle.dump(busco_single_copy_lists, f)
		busco_dic_file_path = Path(busco_dic_dir)
		busco_dic_file_list = list(busco_dic_file_path.iterdir())
		parallelize_singleton_duplicate_segregation(self_hits_dir,self_best_hits_dir,singleton_dir,tmp_dir,self_similarity_cutoff,normalized_bit_score,cores,logger,dpi_no,norm_bit_score_plots_dir,busco_dic_file_list,busco_dic_dir)

		if ref_name != 'NA':
			database_folder_path = Path(database_dir)
			database_file_list = os.listdir(database_folder_path)
			db_names = set()
			for files in database_file_list:
				if "_db" in files:
					orgname = files.split("_db")[0]
					db_names.add(orgname)

			# --- Forward blast of queries vs reference ---#
			for each in db_names:
				if ref_name == str(each):
					fwd_blast_db = os.path.join(database_folder_path, str(each) + '_db')

			self_dummy = False
			singletons_dummy = None
			for every in peptides_file_list:
				blast_pep_org_name = (os.path.basename(str(every))).split('_no_alt_trans')[0]
				if blast_pep_org_name == ref_name:
					pep_ref = every
			parallel_process_organism_fwd_blast(fwd_blast_db, peptides_file_list, ref_name, tmp_dir, blast_dir,fwd_best_hits_dir, tool, blastp, diamond, mmseqs2, eval, self_dummy,singletons_dummy, normalized_bit_score, cores, logger,pep_ref, number_of_hits, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff,db_names,database_folder_path, orthologs_dir, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir, gff_dic_dir)

	# --- if GFF3 and FASTA files of reference organism are provided with one or more of the query organisms missing annotation --- #
	if '--gff' in arguments and '--fasta' in arguments and '--to_annotate' in arguments:
		gff3_folder = arguments[arguments.index('--gff')+1] #full path to GFF3 file or folder containing GFF3 annotation files
		if gff3_folder[-1] == "/":
			gff3_input_files = sorted(glob.glob(gff3_folder + "*.gff") + glob.glob(gff3_folder + "*.gff3") + glob.glob(gff3_folder + "*.gff.gz") + glob.glob(gff3_folder + "*.gff3.gz"))
		else:
			gff3_input_files = [gff3_folder]
		fasta_folder = arguments[arguments.index('--fasta')+1] #full path to folder containing assembly FASTA files
		if fasta_folder[-1] == "/":
			fasta_input_file = fasta_input_file = sorted(glob.glob(fasta_folder + "*.fa") + glob.glob(fasta_folder + "*.fna") + glob.glob(fasta_folder + "*.fna.gz") + glob.glob(fasta_folder + "*.fa.gz") + glob.glob(fasta_folder + "*.fasta.gz") + glob.glob(fasta_folder + "*.fasta") + glob.glob(fasta_folder + "*.genome.fa") + glob.glob(fasta_folder + "*.genome.fna") + glob.glob(fasta_folder + "*.genome.fna.gz") + glob.glob(fasta_folder + "*.genome.fa.gz") + glob.glob(fasta_folder + "*.genome.fasta.gz") + glob.glob(fasta_folder + "*.genome.fasta"))
		else:
			fasta_input_file = [fasta_folder]
		annot_list = arguments[arguments.index('--to_annotate') + 1]  # full path to file containing names of queries to be annotated
		queries = []
		if annot_list[-3:].lower() == "txt":
			queries_to_annotate = {}
			with open(annot_list, 'r') as file:
				line = file.readline()
				while line:
					parts = line.strip().split(',')
					queries_to_annotate[str(parts[0])] = str(parts[1])
					line = file.readline()

		strict_start = False  # option not implemented yet
		strict_end = False  # option not implemented yet
		blast_pep = []
		blast_db = []

		# --- filtering out organisms without annotation and running GeMoMa to obtain their annotations based on the reference annotation--- #
		parallelize_GeMoMa(fasta_input_file, queries_to_annotate, tmp_dir, gff3_input_files,GeMoMa,gff3_folder, cores,logger)

		gff3_input_file = sorted(glob.glob(gff3_folder + "*.gff") + glob.glob(gff3_folder + "*.gff.gz"))

		# Code to get md5sums of input files and write run parameters to a json file

		md5sum_dic = {}
		for each in fasta_input_file:
			file_md5sum = calculate_md5sum(str(each))
			md5sum_dic[str(each)] = file_md5sum
		for every in gff3_input_file:
			file_md5sum = calculate_md5sum(str(every))
			md5sum_dic[str(every)] = file_md5sum
		md5sum_file_path = os.path.join(output_dir, 'md5sum.tsv')
		with open(md5sum_file_path, 'w') as out:
			for file_path, checksum in md5sum_dic.items():
				out.write(f"{file_path}\t{checksum}\n")
		parameters = parse_arguments(arguments)
		with open(os.path.join(parameters["out"], "run_parameters.json"), "w") as f:
			json.dump(parameters, f, indent=4)

		#Calling functions to check and validate the GFF and FASTA files
		validated_fasta_files = os.path.join(tmp_dir, "Validated_FASTA")
		os.makedirs(validated_fasta_files, exist_ok=True)

		error_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(error_files, exist_ok=True)

		warning_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(warning_files, exist_ok=True)
		gemoma=True
		total_errors_final = parallel_clean_and_validate(fasta_input_file, gff3_input_file, validated_fasta_files, error_files, warning_files, 'fasta',cores,logger, gff_config_parameters,gemoma)
		for everything in total_errors_final:
			if len(everything) != 0:
				logger.error('There are errors in your input file(s). Dupylicate analysis exiting.')
				sys.exit()
		validated_fasta_input_file = sorted(glob.glob(os.path.join(validated_fasta_files , "*.fa")) + glob.glob(os.path.join(validated_fasta_files , "*.fa.gz")) + glob.glob(os.path.join(validated_fasta_files , "*.fasta.gz")) + glob.glob(os.path.join(validated_fasta_files , "*.fasta")))

		# --- Calling the functions relevant for getting CDS, PEP, No alternate transcripts PEP --- #
		pep_gene_name_errors = parallel_process_fasta_gff(validated_fasta_input_file, gff3_input_file, tmp_dir, no_alt_trans_dir,no_alt_trans_cds_dir, ref_name, blast_dir, database_dir, self_hits_dir,self_best_hits_dir, singleton_dir, strict_start, strict_end, blast_pep,blast_db, tool, makeblastdb, blastp, diamond, mmseqs2, eval,normalized_bit_score, evo_analysis, 'fasta only',cores,logger,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_type,cds_dir,pep_dir,process_pseudos, validated_tpm_files, exp_input_file,validated_fasta_files,clean_cds_dir,gff_dic_dir,gff_config_parameters,attribute_dic_dir)
		if len(pep_gene_name_errors) != 0:
			logger.error('Please check the protein coding gene names in the GFF file and PEP file. Exiting ...')
			sys.exit()
		if org_type == 'eukaryote':
			peptides_folder_path = Path(no_alt_trans_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		elif org_type == 'prokaryote':
			peptides_folder_path = Path(pep_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		parallel_process_self_blast(ref_name, cores, peptides_file_list, database_dir, tmp_dir, tool, makeblastdb, blastp,diamond, mmseqs2, logger, self_hits_dir, self_best_hits_dir, singleton_dir, eval,normalized_bit_score, max_bit_score_dic_dir, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff)
		# code to do BUSCO based QC
		if normalized_bit_score == 'auto' or quality_check == 'yes':
			busco_cmd = detect_busco(busco_path)
			if busco_cmd:
				logger.info(f"Using BUSCO for QC check of the query organism(s).")
				busco_qc_result = os.path.join(busco_dir_final, 'BUSCO_QC.tsv')
				busco_db_dir = os.path.join(busco_dir, 'busco_databases')
				total_cores = cores
				num_files = len(peptides_file_list)
				busco_threads_per_organism = total_cores
				logger.info(f"Serializing {num_files} organism-level processes for BUSCO analysis")
				logger.info(f"{busco_threads_per_organism} threads inside each organism for BUSCO analysis")
				busco_path_generated = "NA"
				container_workdir = busco_dir
				container_cache_dir = host_cache_dir

				if busco_path == 'busco_docker':
					container_db_dir = busco_db_dir
					docker_image = f"ezlabgva/busco:{busco_version}_{container_version}"
					# writing bash script for executing busco installed with docker
					bash_script = f"""#!/bin/bash
							# BUSCO Docker wrapper script with proper cleanup

							# Container name based on process ID and timestamp for uniqueness
							CONTAINER_NAME="busco_$$_$(date +%s)_$(shuf -i 1000-9999 -n 1)"

							# Cleanup function
							cleanup() {{
								echo "Cleaning up Docker container: $CONTAINER_NAME" >&2
								docker stop "$CONTAINER_NAME" 2>/dev/null || true
								docker rm "$CONTAINER_NAME" 2>/dev/null || true
								exit $1
							}}

							# Set up signal handlers for proper cleanup
							trap 'cleanup 130' INT    # Ctrl+C (SIGINT)
							trap 'cleanup 143' TERM   # Termination (SIGTERM)
							trap 'cleanup 1' EXIT     # Any exit
							trap 'cleanup 1' ERR      # Any error

							# Check if Docker is available
							if ! command -v docker &>/dev/null; then
								echo "Error: Docker is not installed or not in PATH." >&2
								exit 1
							fi

							# Run BUSCO in Docker with automatic cleanup
							docker run --rm \\
								--name "$CONTAINER_NAME" \\
								-u $(id -u) \\
								-v "{host_path}:{container_path}" \\
								-v "{host_cache_dir}:{container_cache_dir}" \\
								-v "{busco_db_dir}:{container_db_dir}" \\
								-e XDG_CONFIG_HOME="{container_cache_dir}" \\
								-w "{container_workdir}" \\
								{docker_image} \\
								busco "$@"

							# Capture exit code and exit cleanly
							EXIT_CODE=$?
							exit $EXIT_CODE
							"""
					output_filename = os.path.join(busco_dir_final, "run_busco_docker.sh")
					with open(output_filename, "w") as f:
						f.write(bash_script)
					os.chmod(output_filename, 0o755)
					busco_path_generated = output_filename
					# Pre-download for Docker
					if not os.path.isdir(busco_db_dir):
						pre_download_databases(busco_path_generated, org_type, busco_dir, busco_db_dir, logger)
				else:
					if not os.path.isdir(busco_db_dir):
						# Pre-download for normal BUSCO
						pre_download_databases(busco_path, org_type, busco_dir, busco_db_dir, logger)
				logger.info("starting actual BUSCO runs")
				ploidy_results = []
				busco_single_copy_lists = {}
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir,str(orgname)+'.pkl')
					if not os.path.isfile(busco_pickle_file):
						if orgname != ref_name:
							if busco_path != 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path, busco_dir, busco_dir_final,str(busco_threads_per_organism), org_type, host_cache_dir,busco_db_dir, busco_pickle_file,logger,buscolineage)
							elif busco_path == 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path_generated, busco_dir,busco_dir_final, str(busco_threads_per_organism), org_type,container_cache_dir, container_db_dir, busco_pickle_file,logger,buscolineage)
				busco_dic_file_path = Path(busco_dic_dir)
				busco_dic_file_list = list(busco_dic_file_path.iterdir())
				with open(busco_qc_result, 'w') as out:
					out.write('Organism' + '\t' + 'Pseudo ploidy number' + '\t' + 'BUSCO Completeness (%)' + '\t' + 'BUSCO Duplication (%)' + '\n')
					for busco_dic_file in busco_dic_file_list:
						with open(busco_dic_file,'rb') as f:
							busco_dic = pickle.load(f)
							if 'feff' in busco_dic:
								out.write(str(busco_dic['feff']) + '\n')
				logger.info(f"BUSCO QC check completed.")
			else:
				logger.info("BUSCO not found. No QC check possible. Using default method to find the threshold to segregate singletons and duplicated genes in the query organism(s)")
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
					busco_single_copy_lists = {}
					with open (busco_pickle_file,'wb') as f:
						pickle.dump(busco_single_copy_lists,f)

		else:
			for fasta_file in peptides_file_list:
				name = get_basename(fasta_file)
				orgname = str(name[0]).split('_no_alt_trans')[0]
				busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
				busco_single_copy_lists = {}
				with open(busco_pickle_file, 'wb') as f:
					pickle.dump(busco_single_copy_lists, f)
		busco_dic_file_path = Path(busco_dic_dir)
		busco_dic_file_list = list(busco_dic_file_path.iterdir())
		parallelize_singleton_duplicate_segregation(self_hits_dir,self_best_hits_dir,singleton_dir,tmp_dir,self_similarity_cutoff,normalized_bit_score,cores,logger,dpi_no,norm_bit_score_plots_dir,busco_dic_file_list,busco_dic_dir)

		if ref_name != 'NA':
			# --- Forward blast of queries vs reference ---#

			database_folder_path = Path(database_dir)
			database_file_list = os.listdir(database_folder_path)
			db_names = set()
			for files in database_file_list:
				if "_db" in files:
					orgname = files.split("_db")[0]
					db_names.add(orgname)

			# --- Forward blast of queries vs reference ---#
			for each in db_names:
				if ref_name == str(each):
					fwd_blast_db = os.path.join(database_folder_path, str(each) + '_db')

			self_dummy = False
			singletons_dummy = None
			for every in peptides_file_list:
				blast_pep_org_name = (os.path.basename(str(every))).split('_no_alt_trans')[0]
				if blast_pep_org_name == ref_name:
					pep_ref = every
			parallel_process_organism_fwd_blast(fwd_blast_db, peptides_file_list, ref_name, tmp_dir, blast_dir,fwd_best_hits_dir, tool, blastp, diamond, mmseqs2, eval, self_dummy,singletons_dummy, normalized_bit_score, cores, logger,pep_ref, number_of_hits, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff,db_names,database_folder_path, orthologs_dir, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir, gff_dic_dir)

	# --- if GFF3 and CDS FASTA files of reference and query organisms are provided --- #
	if '--gff' in arguments and '--cds' in arguments:
		gff3_folder = arguments[arguments.index('--gff')+1] #full path to GFF3 file or folder containing GFF3 annotation files
		if gff3_folder[-1] == "/":
			gff3_input_file = sorted(glob.glob(gff3_folder + "*.gff") + glob.glob(gff3_folder + "*.gff3") + glob.glob(gff3_folder + "*.gff.gz") + glob.glob(gff3_folder + "*.gff3.gz"))
		else:
			gff3_input_file = [gff3_folder]

		cds_folder = arguments[arguments.index('--cds')+1] #full path to CDS FASTA file or folder containing CDS FASTA files
		if cds_folder[-1] == "/":
			cds_input_file = sorted(glob.glob(cds_folder + "*.cds.fa") + glob.glob(cds_folder + "*.cds.fa.gz") + glob.glob(cds_folder + "*.cds.fasta") + glob.glob(cds_folder + "*.cds.fasta.gz") + glob.glob(cds_folder + "*.cds.fna") + glob.glob(cds_folder + "*.cds.fna.gz"))
		else:
			cds_input_file = [cds_folder]
		# Code to get md5sums of input files and write run parameters to a json file
		md5sum_dic = {}
		for each in cds_input_file:
			file_md5sum = calculate_md5sum(str(each))
			md5sum_dic[str(each)] = file_md5sum
		for every in gff3_input_file:
			file_md5sum = calculate_md5sum(str(every))
			md5sum_dic[str(every)] = file_md5sum
		md5sum_file_path = os.path.join(output_dir, 'md5sum.tsv')
		with open(md5sum_file_path, 'w') as out:
			for file_path, checksum in md5sum_dic.items():
				out.write(f"{file_path}\t{checksum}\n")
		parameters = parse_arguments(arguments)
		with open(os.path.join(parameters["out"], "run_parameters.json"), "w") as f:
			json.dump(parameters, f, indent=4)

		strict_start = False  # option not implemented yet
		strict_end = False  # option not implemented yet
		blast_pep = []
		blast_db = []

		validated_fasta_files = os.path.join(tmp_dir, "Validated_CDS")
		os.makedirs(validated_fasta_files, exist_ok=True)

		error_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(error_files, exist_ok=True)

		warning_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(warning_files, exist_ok=True)
		gemoma=False
		# Calling functions to check and validate the GFF and FASTA files
		total_errors_final = parallel_clean_and_validate(cds_input_file, gff3_input_file, validated_fasta_files, error_files, warning_files, 'cds',cores,logger, gff_config_parameters,gemoma)
		for everything in total_errors_final:
			if len(everything) != 0:
				logger.error('There are errors in your input file(s). Dupylicate analysis exiting.')
				sys.exit()
		validated_cds_input_file = sorted(glob.glob(os.path.join(validated_fasta_files , "*.fa")) + glob.glob(os.path.join(validated_fasta_files , "*.fa.gz")) + glob.glob(os.path.join(validated_fasta_files , "*.fasta.gz")) + glob.glob(os.path.join(validated_fasta_files , "*.fasta")))
		# --- Calling the functions relevant for getting PEP, No alternate transcripts PEP --- #
		pep_gene_name_errors = parallel_process_fasta_gff(validated_cds_input_file, gff3_input_file, tmp_dir, no_alt_trans_dir,no_alt_trans_cds_dir, ref_name, blast_dir, database_dir, self_hits_dir,self_best_hits_dir, singleton_dir, strict_start, strict_end, blast_pep,blast_db, tool, makeblastdb, blastp, diamond, mmseqs2, eval,normalized_bit_score, evo_analysis, 'cds fasta',cores,logger,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_type,cds_dir,pep_dir,process_pseudos, validated_tpm_files, exp_input_file,validated_fasta_files,clean_cds_dir,gff_dic_dir,gff_config_parameters,attribute_dic_dir)
		if len(pep_gene_name_errors) != 0:
			logger.error('Please check the protein coding gene names in the GFF file and PEP file. Exiting ...')
			sys.exit()
		if org_type == 'eukaryote':
			peptides_folder_path = Path(no_alt_trans_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		elif org_type == 'prokaryote':
			peptides_folder_path = Path(pep_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		parallel_process_self_blast(ref_name, cores, peptides_file_list, database_dir, tmp_dir, tool, makeblastdb, blastp,diamond, mmseqs2, logger, self_hits_dir, self_best_hits_dir, singleton_dir, eval,normalized_bit_score, max_bit_score_dic_dir, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff)
		# code to do BUSCO based QC
		if normalized_bit_score == 'auto' or quality_check == 'yes':
			busco_cmd = detect_busco(busco_path)
			if busco_cmd:
				logger.info(f"Using BUSCO for QC check of the query organism(s).")
				busco_qc_result = os.path.join(busco_dir_final, 'BUSCO_QC.tsv')
				busco_db_dir = os.path.join(busco_dir, 'busco_databases')
				total_cores = cores
				num_files = len(peptides_file_list)
				busco_threads_per_organism = total_cores
				logger.info(f"Serializing {num_files} organism-level processes for BUSCO analysis")
				logger.info(f"{busco_threads_per_organism} threads inside each organism for BUSCO analysis")
				busco_path_generated = "NA"
				container_workdir = busco_dir
				container_cache_dir = host_cache_dir

				if busco_path == 'busco_docker':
					container_db_dir = busco_db_dir
					docker_image = f"ezlabgva/busco:{busco_version}_{container_version}"
					# writing bash script for executing busco installed with docker
					bash_script = f"""#!/bin/bash
							# BUSCO Docker wrapper script with proper cleanup

							# Container name based on process ID and timestamp for uniqueness
							CONTAINER_NAME="busco_$$_$(date +%s)_$(shuf -i 1000-9999 -n 1)"

							# Cleanup function
							cleanup() {{
							echo "Cleaning up Docker container: $CONTAINER_NAME" >&2
							docker stop "$CONTAINER_NAME" 2>/dev/null || true
							docker rm "$CONTAINER_NAME" 2>/dev/null || true
							exit $1
							}}

							# Set up signal handlers for proper cleanup
							trap 'cleanup 130' INT    # Ctrl+C (SIGINT)
							trap 'cleanup 143' TERM   # Termination (SIGTERM)
							trap 'cleanup 1' EXIT     # Any exit
							trap 'cleanup 1' ERR      # Any error

							# Check if Docker is available
							if ! command -v docker &>/dev/null; then
							echo "Error: Docker is not installed or not in PATH." >&2
							exit 1
							fi

							# Run BUSCO in Docker with automatic cleanup
							docker run --rm \\
							--name "$CONTAINER_NAME" \\
							-u $(id -u) \\
							-v "{host_path}:{container_path}" \\
							-v "{host_cache_dir}:{container_cache_dir}" \\
							-v "{busco_db_dir}:{container_db_dir}" \\
							-e XDG_CONFIG_HOME="{container_cache_dir}" \\
							-w "{container_workdir}" \\
							{docker_image} \\
							busco "$@"

							# Capture exit code and exit cleanly
							EXIT_CODE=$?
							exit $EXIT_CODE
							"""
					output_filename = os.path.join(busco_dir_final, "run_busco_docker.sh")
					with open(output_filename, "w") as f:
						f.write(bash_script)
					os.chmod(output_filename, 0o755)
					busco_path_generated = output_filename
					# Pre-download for Docker
					if not os.path.isdir(busco_db_dir):
						pre_download_databases(busco_path_generated, org_type, busco_dir, busco_db_dir, logger)
				else:
					if not os.path.isdir(busco_db_dir):
						# Pre-download for normal BUSCO
						pre_download_databases(busco_path, org_type, busco_dir, busco_db_dir, logger)
				logger.info("starting actual BUSCO runs")
				ploidy_results = []
				busco_single_copy_lists = {}
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir,str(orgname)+'.pkl')
					if not os.path.isfile(busco_pickle_file):
						if orgname != ref_name:
							if busco_path != 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path, busco_dir, busco_dir_final,str(busco_threads_per_organism), org_type, host_cache_dir,busco_db_dir, busco_pickle_file,logger,buscolineage)
							elif busco_path == 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path_generated, busco_dir,busco_dir_final, str(busco_threads_per_organism), org_type,container_cache_dir, container_db_dir, busco_pickle_file,logger,buscolineage)
				busco_dic_file_path = Path(busco_dic_dir)
				busco_dic_file_list = list(busco_dic_file_path.iterdir())
				with open(busco_qc_result, 'w') as out:
					out.write('Organism' + '\t' + 'Pseudo ploidy number' + '\t' + 'BUSCO Completeness (%)' + '\t' + 'BUSCO Duplication (%)' + '\n')
					for busco_dic_file in busco_dic_file_list:
						with open(busco_dic_file,'rb') as f:
							busco_dic = pickle.load(f)
							if 'feff' in busco_dic:
								out.write(str(busco_dic['feff']) + '\n')
				logger.info(f"BUSCO QC check completed.")
			else:
				logger.info("BUSCO not found. No QC check possible. Using default method to find the threshold to segregate singletons and duplicated genes in the query organism(s)")
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
					busco_single_copy_lists = {}
					with open (busco_pickle_file,'wb') as f:
						pickle.dump(busco_single_copy_lists,f)

		else:
			for fasta_file in peptides_file_list:
				name = get_basename(fasta_file)
				orgname = str(name[0]).split('_no_alt_trans')[0]
				busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
				busco_single_copy_lists = {}
				with open(busco_pickle_file, 'wb') as f:
					pickle.dump(busco_single_copy_lists, f)
		busco_dic_file_path = Path(busco_dic_dir)
		busco_dic_file_list = list(busco_dic_file_path.iterdir())
		parallelize_singleton_duplicate_segregation(self_hits_dir,self_best_hits_dir,singleton_dir,tmp_dir,self_similarity_cutoff,normalized_bit_score,cores,logger,dpi_no,norm_bit_score_plots_dir,busco_dic_file_list,busco_dic_dir)

		if ref_name != 'NA':
			# --- Forward blast of queries vs reference ---#

			database_folder_path = Path(database_dir)
			database_file_list = os.listdir(database_folder_path)
			db_names = set()
			for files in database_file_list:
				if "_db" in files:
					orgname = files.split("_db")[0]
					db_names.add(orgname)
			# --- Forward blast of queries vs reference ---#
			for each in db_names:
				if ref_name == str(each):
					fwd_blast_db = os.path.join(database_folder_path, str(each) + '_db')
			self_dummy = False
			singletons_dummy = None
			for every in peptides_file_list:
				blast_pep_org_name = (os.path.basename(str(every))).split('_no_alt_trans')[0]
				if blast_pep_org_name == ref_name:
					pep_ref = every
			parallel_process_organism_fwd_blast(fwd_blast_db, peptides_file_list, ref_name, tmp_dir, blast_dir,fwd_best_hits_dir, tool, blastp, diamond, mmseqs2, eval, self_dummy,singletons_dummy, normalized_bit_score, cores, logger,pep_ref, number_of_hits, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff,db_names,database_folder_path, orthologs_dir, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir, gff_dic_dir)

	# --- if GFF3 and PEP FASTA files of reference and query organisms are provided --- #
	if '--gff' in arguments and '--pep' in arguments:
		gff3_folder = arguments[arguments.index('--gff')+1] #full path to GFF3 file or folder containing GFF3 annotation files
		if gff3_folder[-1] == "/":
			gff3_input_file = sorted(glob.glob(gff3_folder + "*.gff") + glob.glob(gff3_folder + "*.gff3") + glob.glob(gff3_folder + "*.gff.gz") + glob.glob(gff3_folder + "*.gff3.gz"))
		else:
			gff3_input_file = [gff3_folder]
		pep_folder = arguments[arguments.index('--pep')+1] #full path to folder containing PEP FASTA files
		if pep_folder[-1] == "/":
			pep_input_file = sorted(glob.glob(pep_folder + "*.pep.fa") + glob.glob(pep_folder + "*.pep.fa.gz") + glob.glob(pep_folder + "*.pep.fasta") + glob.glob(pep_folder + "*.pep.fasta.gz") + glob.glob(pep_folder + "*.pep.fna") + glob.glob(pep_folder + "*.pep.fna.gz"))
		else:
			pep_input_file = [pep_folder]

		# Code to get md5sums of input files and write run parameters to a json file
		md5sum_dic = {}
		for each in pep_input_file:
			file_md5sum = calculate_md5sum(str(each))
			md5sum_dic[str(each)] = file_md5sum
		for every in gff3_input_file:
			file_md5sum = calculate_md5sum(str(every))
			md5sum_dic[str(every)] = file_md5sum
		md5sum_file_path = os.path.join(output_dir, 'md5sum.tsv')
		with open(md5sum_file_path, 'w') as out:
			for file_path, checksum in md5sum_dic.items():
				out.write(f"{file_path}\t{checksum}\n")
		parameters = parse_arguments(arguments)
		with open(os.path.join(parameters["out"], "run_parameters.json"), "w") as f:
			json.dump(parameters, f, indent=4)
		
		strict_start = False  # option not implemented yet
		strict_end = False  # option not implemented yet
		blast_pep = []
		blast_db = []

		validated_fasta_files = os.path.join(tmp_dir, "Validated_PEP")
		os.makedirs(validated_fasta_files, exist_ok=True)

		error_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(error_files, exist_ok=True)

		warning_files = os.path.join(tmp_dir, "Errors_warnings")
		os.makedirs(warning_files, exist_ok=True)
		gemoma=False
		# Calling functions to check and validate the GFF and FASTA files
		total_errors_final = parallel_clean_and_validate(pep_input_file, gff3_input_file, validated_fasta_files, error_files, warning_files, 'pep',cores,logger, gff_config_parameters,gemoma)
		for everything in total_errors_final:
			if len(everything) != 0:
				logger.error('There are errors in your input file(s). Dupylicate analysis exiting.')
				sys.exit()
		if org_type != 'prokaryote':
			# copying the files in validated PEP files folder to pep_dir folder for the later tpmclean steps when user has specific genes to analyse
			for filename in os.listdir(validated_fasta_files):
				src_path = os.path.join(validated_fasta_files, filename)
				dst_path = os.path.join(pep_dir, filename)

				# Copy only files (not subfolders)
				if os.path.isfile(src_path):
					shutil.copy2(src_path, dst_path)  # copy2 preserves metadata
		validated_pep_input_file = sorted(glob.glob(os.path.join(validated_fasta_files, "*.fa")) + glob.glob(os.path.join(validated_fasta_files , "*.fa.gz")) + glob.glob(os.path.join(validated_fasta_files , "*.fasta.gz")) + glob.glob(os.path.join(validated_fasta_files , "*.fasta")))
		# --- Calling the functions relevant for getting PEP, No alternate transcripts PEP --- #
		pep_gene_name_errors = parallel_process_fasta_gff(validated_pep_input_file, gff3_input_file, tmp_dir, no_alt_trans_dir,no_alt_trans_cds_dir, ref_name, blast_dir, database_dir, self_hits_dir,self_best_hits_dir, singleton_dir, strict_start, strict_end, blast_pep,blast_db, tool, makeblastdb, blastp, diamond, mmseqs2, eval,normalized_bit_score, evo_analysis, 'pep fasta',cores,logger,mafft,occupancy,fasttree,fwd_similarity_cutoff,score_ratio_cutoff,self_similarity_cutoff,org_type,cds_dir,pep_dir,process_pseudos, validated_tpm_files, exp_input_file,validated_fasta_files,clean_cds_dir,gff_dic_dir,gff_config_parameters,attribute_dic_dir)
		if len(pep_gene_name_errors) != 0:
			logger.error('Please check the protein coding gene names in the GFF file and PEP file. Exiting ...')
			sys.exit()
		if org_type == 'eukaryote':
			peptides_folder_path = Path(no_alt_trans_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		elif org_type == 'prokaryote':
			peptides_folder_path = Path(pep_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
		parallel_process_self_blast(ref_name, cores, peptides_file_list, database_dir, tmp_dir, tool, makeblastdb, blastp,diamond, mmseqs2, logger, self_hits_dir, self_best_hits_dir, singleton_dir, eval,normalized_bit_score, max_bit_score_dic_dir, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff)
		# code to do BUSCO based QC
		if normalized_bit_score == 'auto' or quality_check == 'yes':
			busco_cmd = detect_busco(busco_path)
			if busco_cmd:
				logger.info(f"Using BUSCO for QC check of the query organism(s).")
				busco_qc_result = os.path.join(busco_dir_final, 'BUSCO_QC.tsv')
				busco_db_dir = os.path.join(busco_dir, 'busco_databases')
				total_cores = cores
				num_files = len(peptides_file_list)
				busco_threads_per_organism = total_cores
				logger.info(f"Serializing {num_files} organism-level processes for BUSCO analysis")
				logger.info(f"{busco_threads_per_organism} threads inside each organism for BUSCO analysis")
				busco_path_generated = "NA"
				container_workdir = busco_dir
				container_cache_dir = host_cache_dir

				if busco_path == 'busco_docker':
					container_db_dir = busco_db_dir
					docker_image = f"ezlabgva/busco:{busco_version}_{container_version}"
					# writing bash script for executing busco installed with docker
					bash_script = f"""#!/bin/bash
							# BUSCO Docker wrapper script with proper cleanup

							# Container name based on process ID and timestamp for uniqueness
							CONTAINER_NAME="busco_$$_$(date +%s)_$(shuf -i 1000-9999 -n 1)"

							# Cleanup function
							cleanup() {{
							echo "Cleaning up Docker container: $CONTAINER_NAME" >&2
							docker stop "$CONTAINER_NAME" 2>/dev/null || true
							docker rm "$CONTAINER_NAME" 2>/dev/null || true
							exit $1
							}}

							# Set up signal handlers for proper cleanup
							trap 'cleanup 130' INT    # Ctrl+C (SIGINT)
							trap 'cleanup 143' TERM   # Termination (SIGTERM)
							trap 'cleanup 1' EXIT     # Any exit
							trap 'cleanup 1' ERR      # Any error

							# Check if Docker is available
							if ! command -v docker &>/dev/null; then
							echo "Error: Docker is not installed or not in PATH." >&2
							exit 1
							fi

							# Run BUSCO in Docker with automatic cleanup
							docker run --rm \\
							--name "$CONTAINER_NAME" \\
							-u $(id -u) \\
							-v "{host_path}:{container_path}" \\
							-v "{host_cache_dir}:{container_cache_dir}" \\
							-v "{busco_db_dir}:{container_db_dir}" \\
							-e XDG_CONFIG_HOME="{container_cache_dir}" \\
							-w "{container_workdir}" \\
							{docker_image} \\
							busco "$@"

							# Capture exit code and exit cleanly
							EXIT_CODE=$?
							exit $EXIT_CODE
							"""
					output_filename = os.path.join(busco_dir_final, "run_busco_docker.sh")
					with open(output_filename, "w") as f:
						f.write(bash_script)
					os.chmod(output_filename, 0o755)
					busco_path_generated = output_filename
					# Pre-download for Docker
					if not os.path.isdir(busco_db_dir):
						pre_download_databases(busco_path_generated, org_type, busco_dir, busco_db_dir, logger)
				else:
					if not os.path.isdir(busco_db_dir):
						# Pre-download for normal BUSCO
						pre_download_databases(busco_path, org_type, busco_dir, busco_db_dir, logger)
				logger.info("starting actual BUSCO runs")
				ploidy_results = []
				busco_single_copy_lists = {}
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir,str(orgname)+'.pkl')
					if not os.path.isfile(busco_pickle_file):
						if orgname != ref_name:
							if busco_path != 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path, busco_dir, busco_dir_final,str(busco_threads_per_organism), org_type, host_cache_dir,busco_db_dir, busco_pickle_file,logger,buscolineage)
							elif busco_path == 'busco_docker':
								run_busco(orgname, str(fasta_file), busco_path_generated, busco_dir,busco_dir_final, str(busco_threads_per_organism), org_type,container_cache_dir, container_db_dir, busco_pickle_file,logger,buscolineage)
				busco_dic_file_path = Path(busco_dic_dir)
				busco_dic_file_list = list(busco_dic_file_path.iterdir())
				with open(busco_qc_result, 'w') as out:
					out.write('Organism' + '\t' + 'Pseudo ploidy number' + '\t' + 'BUSCO Completeness (%)' + '\t' + 'BUSCO Duplication (%)' + '\n')
					for busco_dic_file in busco_dic_file_list:
						with open(busco_dic_file,'rb') as f:
							busco_dic = pickle.load(f)
							if 'feff' in busco_dic:
								out.write(str(busco_dic['feff']) + '\n')
				logger.info(f"BUSCO QC check completed.")
			else:
				logger.info("BUSCO not found. No QC check possible. Using default method to find the threshold to segregate singletons and duplicated genes in the query organism(s)")
				for fasta_file in peptides_file_list:
					name = get_basename(fasta_file)
					orgname = str(name[0]).split('_no_alt_trans')[0]
					busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
					busco_single_copy_lists = {}
					with open (busco_pickle_file,'wb') as f:
						pickle.dump(busco_single_copy_lists,f)
		else:
			for fasta_file in peptides_file_list:
				name = get_basename(fasta_file)
				orgname = str(name[0]).split('_no_alt_trans')[0]
				busco_pickle_file = os.path.join(busco_dic_dir, str(orgname) + '.pkl')
				busco_single_copy_lists = {}
				with open(busco_pickle_file, 'wb') as f:
					pickle.dump(busco_single_copy_lists, f)
		busco_dic_file_path = Path(busco_dic_dir)
		busco_dic_file_list = list(busco_dic_file_path.iterdir())
		parallelize_singleton_duplicate_segregation(self_hits_dir,self_best_hits_dir,singleton_dir,tmp_dir,self_similarity_cutoff,normalized_bit_score,cores,logger,dpi_no,norm_bit_score_plots_dir,busco_dic_file_list,busco_dic_dir)

		if ref_name != 'NA':
			#debug code
			peptides_folder_path = Path(no_alt_trans_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
			# --- Forward blast of queries vs reference ---#

			database_folder_path = Path(database_dir)
			database_file_list = os.listdir(database_folder_path)
			db_names = set()
			for files in database_file_list:
				if "_db" in files:
					orgname = files.split("_db")[0]
					db_names.add(orgname)

		# --- Forward blast of queries vs reference ---#
			for each in db_names:
				if ref_name == str(each):
					fwd_blast_db = os.path.join(database_folder_path, str(each)+'_db')
			self_dummy = False
			singletons_dummy = None
			for every in peptides_file_list:
				blast_pep_org_name = (os.path.basename(str(every))).split('_no_alt_trans')[0]
				if blast_pep_org_name == ref_name:
					pep_ref = every
			parallel_process_organism_fwd_blast(fwd_blast_db, peptides_file_list, ref_name, tmp_dir, blast_dir,fwd_best_hits_dir, tool, blastp, diamond, mmseqs2, eval, self_dummy,singletons_dummy, normalized_bit_score, cores, logger,pep_ref, number_of_hits, mafft, occupancy, fasttree,fwd_similarity_cutoff, score_ratio_cutoff, self_similarity_cutoff,db_names,database_folder_path, orthologs_dir, synteny_cutoff, flank_number, best_side, possible_orthologs,do_tree, max_bit_score_dic_dir, gff_dic_dir)

	# making a list of files in the no alternate transcripts peptide output, forward blast best hits, self blast second best hits and self blast hit folders for finding the gene duplicates
	if org_type == 'eukaryote':
		peptides_folder_path = Path(no_alt_trans_dir)
		peptides_file_list = list(peptides_folder_path.iterdir())
		trans_peptides_folder_path = Path(pep_dir)
		trans_peptides_file_list = list(trans_peptides_folder_path.iterdir())
	elif org_type == 'prokaryote':
		peptides_folder_path = Path(pep_dir)
		peptides_file_list = list(peptides_folder_path.iterdir())
	fwd_folder_path = Path(fwd_best_hits_dir)
	fwd_file_list = list(fwd_folder_path.iterdir())
	self_best_folder_path = Path(self_best_hits_dir)
	self_file_list = list(self_best_folder_path.iterdir())
	self_folder_path = Path(self_hits_dir)
	self_list = list(self_folder_path.iterdir())
	ortholog_folder_path=Path(orthologs_dir)
	orthologs_file_list = list(ortholog_folder_path.iterdir())

	# code lines to get the final gene duplicates outputs
	# code to get the positional dictionary of the reference organism
	logger.info('Starting to find gene duplicates in query organism(s).')
	gff_posdic_folder_path = Path(gff_dic_dir)
	gff_posdic_file_list = list(gff_posdic_folder_path.iterdir())
	if ref_name != 'NA':
		for gffdic_file in gff_posdic_file_list:
			gffdic_org_name = get_basename(str(gffdic_file))[0]
			if gffdic_org_name == ref_name:
				with open(gffdic_file, 'rb') as f:
					gff_pos_details = pickle.load(f)
					pos_dic_ref, pos_nos_ref, coding_non_coding_ref = gff_pos_details[gffdic_org_name]
		tmp_singleton_files = glob.glob(os.path.join(singleton_dir, "*.tsv"))
		# code to get the positional dictionary of the query organism and final grouping and writing of gene duplicates
		parallelize_duplicates_processing(gff3_input_file, peptides_file_list, self_list, self_file_list, fwd_file_list,orthologs_file_list, tmp_singleton_files, proximity, no_alt_trans_dir, evo_analysis,pos_dic_ref, synteny_cutoff, flank_number, best_side,mafft, kaks_dir, singleton_dir_final, classify, cores, ref_name, tandems_dir,proximals_dir, dispersed_dir, mixed_dir, unclassified_dir_final, output_dir,logger,process_pseudos, gff_posdic_file_list,frequency_plots_dir,int(dpi_no))
	else:
		tmp_singleton_files = glob.glob(os.path.join(singleton_dir, "*.tsv"))
		pos_dic_ref = {}
		# code to get the positional dictionary of the query organism and final grouping and writing of gene duplicates
		parallelize_duplicates_processing(gff3_input_file, peptides_file_list, self_list, self_file_list, fwd_file_list,orthologs_file_list, tmp_singleton_files, proximity, no_alt_trans_dir, evo_analysis,pos_dic_ref, synteny_cutoff, flank_number, best_side,mafft, kaks_dir, singleton_dir_final, classify, cores, ref_name, tandems_dir,proximals_dir, dispersed_dir, mixed_dir, unclassified_dir_final, output_dir,logger,process_pseudos, gff_posdic_file_list, frequency_plots_dir, int(dpi_no))

	#if the user has specified a single or list of desired genes in the presence of a reference organism
	if len(desired_genes)!=0 and ref_name != 'NA':
		desired_gene_names=[]
		for every in gff3_input_file:
			gff_name = get_basename(str(every))[0]
			if org_type == 'eukaryote':
				for pep in trans_peptides_file_list:
					pepname = get_basename(str(pep))[0]
					if gff_name == pepname and gff_name == ref_name:
						if 'all' in gff_config_parameters:
							child_attribute = gff_config_parameters['all']['child_attribute']
							child_parent_linker = gff_config_parameters['all']['child_parent_linker']
							parent_attribute = gff_config_parameters['all']['parent_attribute']
						elif 'default' in gff_config_parameters:
							child_attribute = gff_config_parameters['default']['child_attribute']
							child_parent_linker = gff_config_parameters['default']['child_parent_linker']
							parent_attribute = gff_config_parameters['default']['parent_attribute']
						else:
							child_attribute = gff_config_parameters[ref_name]['child_attribute']
							child_parent_linker = gff_config_parameters[ref_name]['child_parent_linker']
							parent_attribute = gff_config_parameters[ref_name]['parent_attribute']
						transcripts_per_gene, gene_names = tpmclean(str(every), str(pep),child_parent_linker,child_attribute)
						# Step 1: Replace transcript names with gene names
						# Create a mapping of transcript to gene for quick lookup
						transcript_to_gene = {transcript: gene for gene, transcripts in transcripts_per_gene.items() for transcript in transcripts}
						for transcript in desired_genes:
							gene_name = transcript_to_gene.get(transcript,transcript)  # Default to transcript name if not found
							desired_gene_names.append(gene_name)
						desired_genes_set = set(desired_gene_names)
			elif org_type == 'prokaryote':
				for pep in peptides_file_list:
					pepname = (get_basename(str(pep))[0]).strip().split('no_alt_trans')[0]
					if gff_name == pepname and gff_name == ref_name:
						if 'all' in gff_config_parameters:
							child_attribute = gff_config_parameters['all']['child_attribute']
							child_parent_linker = gff_config_parameters['all']['child_parent_linker']
							parent_attribute = gff_config_parameters['all']['parent_attribute']
						elif 'default' in gff_config_parameters:
							child_attribute = gff_config_parameters['default']['child_attribute']
							child_parent_linker = gff_config_parameters['default']['child_parent_linker']
							parent_attribute = gff_config_parameters['default']['parent_attribute']
						else:
							child_attribute = gff_config_parameters[ref_name]['child_attribute']
							child_parent_linker = gff_config_parameters[ref_name]['child_parent_linker']
							parent_attribute = gff_config_parameters[ref_name]['parent_attribute']
						transcripts_per_gene, gene_names = tpmclean(str(every), str(pep),child_parent_linker,child_attribute)
						# Step 1: Replace transcript names with gene names
						# Create a mapping of transcript to gene for quick lookup
						transcript_to_gene = {transcript: gene for gene, transcripts in transcripts_per_gene.items() for transcript in transcripts}
						for transcript in desired_genes:
							gene_name = transcript_to_gene.get(transcript,transcript)  # Default to transcript name if not found
							desired_gene_names.append(gene_name)
						desired_genes_set = set(desired_gene_names)
		desired_genes_output_file = os.path.join(output_dir, "Specified_genes_results.tsv")
		merged_rows_tandems = lists_for_specified_results(tandems_dir, desired_genes_set, ref_name)
		merged_rows_proximals = lists_for_specified_results(proximals_dir, desired_genes_set, ref_name)
		merged_rows_proximals_new = merged_rows_proximals[1:]
		merged_rows_dispersed = lists_for_specified_results(dispersed_dir, desired_genes_set, ref_name)
		merged_rows_dispersed_new = merged_rows_dispersed[1:]
		if classify == 'strict':
			merged_rows_mixed = lists_for_specified_results(mixed_dir, desired_genes_set, ref_name)
			merged_rows_mixed_new = merged_rows_mixed[1:]
		else:
			merged_rows_mixed_new = []
		merged_rows_singletons = lists_for_specified_results(singleton_dir_final, desired_genes_set, ref_name)
		merged_rows_singletons_new = merged_rows_singletons[1:]
		merged_rows_total = merged_rows_tandems + merged_rows_proximals_new + merged_rows_dispersed_new + merged_rows_mixed_new + merged_rows_singletons_new
		# code to check if there are any desired genes that do not have a hit at all
		found_ref_genes = set()
		for row in merged_rows_total[1:]:  # Skip header
			ref_genes = row[0].replace("'", "").split(',')  # Remove quotes and split
			found_ref_genes.update(g.strip() for g in ref_genes)

		# Loop through desired_genes and find the ones that do not match
		for gene in desired_genes_set:
			if gene not in found_ref_genes:
				# Construct a new row: gene in refrence, '--' for each organism column, 'not found' as classification
				num_orgs = len(merged_rows_total[0]) - 2  # Total columns minus reference and classification
				new_row = [f"{gene}"] + ['--'] * num_orgs + ['No hits in queries']
				merged_rows_total.append(new_row)
		# writing to specified_gene_results output file
		with open(desired_genes_output_file, 'w', newline='') as f:
			for row in merged_rows_total:
				line = '\t'.join(row)  # Tab-separated values
				f.write(line + '\n')
		#code lines to condense the above output file like a comparative genomics file with the number of copies of the reference gens in each of the query organisms
		# Step 1: Load the input tab-separated file
		# Load DataFrame
		df = pd.read_csv(desired_genes_output_file,sep="\t")
		if 'classification' in df.columns:
			df = df.drop(columns=['classification'])
		# Step 3: Define reference and organism columns
		ref_col = df.columns[0]
		other_cols = df.columns[1:]

		# Step 4: Build gene map from reference to organism genes
		gene_map = defaultdict(lambda: defaultdict(list))

		for _, row in df.iterrows():
			ref_genes = [g.strip().strip("'") for g in str(row[ref_col]).split(',')]

			for org in other_cols:
				if row[org] != '--':
					org_genes = [g.strip().strip("'") for g in str(row[org]).split(',')]

					# Case 1: Equal number of ref and org genes -> use positional mapping
					if len(ref_genes) == len(org_genes):
						for i, ref_gene in enumerate(ref_genes):
							gene_map[ref_gene][org].append(org_genes[i])

					# Case 2: Single ref gene with multiple org genes -> assign all to ref gene
					elif len(ref_genes) == 1:
						ref_gene = ref_genes[0]
						for org_gene in org_genes:
							gene_map[ref_gene][org].append(org_gene)

					# Case 3: Multiple ref genes with single org gene -> assign to all ref genes
					elif len(org_genes) == 1:
						org_gene = org_genes[0]
						for ref_gene in ref_genes:
							gene_map[ref_gene][org].append(org_gene)

					# Case 4: Check if org_genes is a multiple of ref_genes (combined groups)
					elif len(org_genes) % len(ref_genes) == 0:
						# Multiple groups combined - assign genes in chunks
						chunk_size = len(ref_genes)
						for chunk_start in range(0, len(org_genes), chunk_size):
							chunk_genes = org_genes[chunk_start:chunk_start + chunk_size]
							for i, ref_gene in enumerate(ref_genes):
								gene_map[ref_gene][org].append(chunk_genes[i])

					# Case 5: Mismatched counts - fall back to positional mapping
					else:
						for i, ref_gene in enumerate(ref_genes):
							if i < len(org_genes):
								gene_map[ref_gene][org].append(org_genes[i])

		# Step 5: Load orthologs scores for OGCS calculation
		orthologs_scores = {}
		for org in other_cols:
			orthologs_file = os.path.join(orthologs_dir, f"{org}_orthologs.tsv")
			if os.path.exists(orthologs_file):
				scores_dict = {}
				try:
					with open(orthologs_file, 'r') as f:
						for line_num, line in enumerate(f, 1):
							# Skip header line
							if line_num == 1:
								continue

							parts = line.strip().split('\t')
							if len(parts) >= 3:
								gene_name = parts[0].strip()
								try:
									score = float(parts[2].strip())
									scores_dict[gene_name] = score
								except ValueError:
									logger.debug(f"Warning: Invalid score in {org}_orthologs.tsv line {line_num}: {parts[2]}")
									continue
					orthologs_scores[org] = scores_dict
					logger.info(f"Loaded {len(scores_dict)} orthology gene scores for {org}")
				except Exception as e:
					logger.error(f"Warning: Could not read orthologs file for {org}: {e}")
					orthologs_scores[org] = {}
			else:
				logger.error(f"Warning: Orthologs file not found for {org}: {orthologs_file}")
				orthologs_scores[org] = {}

		# Step 6: Prepare final output rows with OGCS calculation
		final_rows = []
		for gene in sorted(gene_map.keys()):
			if gene in desired_genes_set:
				row = {ref_col: gene}
				for org in other_cols:
					genes = gene_map[gene].get(org, [])
					if genes:
						gene_list = sorted(set(genes))

						# Calculate OGCS if orthologs scores are available
						if org in orthologs_scores and orthologs_scores[org]:
							scores = []
							missing_genes = []
							for g in gene_list:
								if g in orthologs_scores[org]:
									scores.append(orthologs_scores[org][g])
								else:
									missing_genes.append(g)

							if scores:  # Only calculate OGCS if we have scores for at least some genes
								ogcs = sum(scores) / len(gene_list)  # Divide by total genes in cell
								if missing_genes:
									# Note: OGCS calculated using available scores divided by total gene count
									row[org] = f"{', '.join(gene_list)} ({len(gene_list)}) {{OGCS={ogcs:.3f}*}}"
								else:
									row[org] = f"{', '.join(gene_list)} ({len(gene_list)}) {{OGCS={ogcs:.3f}}}"
							else:
								# No scores found for any genes in this cell
								row[org] = f"{', '.join(gene_list)} ({len(gene_list)}) {{OGCS=N/A}}"
						else:
							# No orthologs file or empty orthologs file
							row[org] = f"{', '.join(gene_list)} ({len(gene_list)})"
					else:
						row[org] = "-- (0)"
				final_rows.append(row)

		# Step 7: Add allowed genes that had no matches at all
		for gene in sorted(desired_genes_set):
			if gene not in gene_map:
				row = {ref_col: gene}
				for org in other_cols:
					row[org] = "-- (0)"
				final_rows.append(row)

		# Step 8: Create final DataFrame and write to file
		final_df = pd.DataFrame(final_rows)
		final_df = final_df.sort_values(by=ref_col)
		desired_genes_condensed_output_file = os.path.join(output_dir, "Copy_number_table.tsv")

		# Step 9: Write to output file as tab-separated
		final_df.to_csv(desired_genes_condensed_output_file, sep="\t", index=False)

		#Step 10: Make a separate excel table file with zero copy number cells coloured in red
		def highlight_absent(val):
			return 'background-color: lightcoral' if val == "-- (0)" else ""
		styled = final_df.style.map(highlight_absent)
		# Step 3: Write Excel file
		excel_file = os.path.join(output_dir, "Copy_number_table.xlsx")
		styled.to_excel(excel_file, engine="openpyxl", index=False)

		if specific_analysis=='no':
			pass
		elif specific_analysis=='yes':
			is_empty = len(os.listdir(validated_tpm_files)) == 0
			if is_empty == False:
				# cleaning the TPM files
				if org_type == 'eukaryote':
					validated_tpm_input_file = sorted(glob.glob(os.path.join(validated_tpm_files, "*.txt")) + glob.glob(os.path.join(validated_tpm_files, "*.txt.gz")) + glob.glob(os.path.join(validated_tpm_files, "*.tpms.txt")) + glob.glob(os.path.join(validated_tpm_files, "*.tpms.txt.gz")))
					logger.info('TPM files have been cleaned')
				elif org_type == 'prokaryote':
					validated_tpm_input_file = exp_input_file
				master_folder = os.path.join(output_dir,'Specific_duplicates_analysis')  # Creating the path for the folder for storing the data of particular organism
				os.makedirs(master_folder,exist_ok=True)  # Creates an organism specific folder in the specified path location
				parallelize_rthresh_pseudogenes(master_folder, validated_tpm_input_file, rand_genes, dpi_no, cores,logger,avg_exp_cutoff,rthreshold_pseudogenes_dic_dir)
				rthreshold_pseudogenes_dic_dir_path = Path(rthreshold_pseudogenes_dic_dir)
				rthreshold_pseudogenes_file_list = list(rthreshold_pseudogenes_dic_dir_path.iterdir())
				specific_two_step_parallelization_master(validated_tpm_input_file, master_folder, desired_genes_output_file,dpi_no, rand_genes,rthreshold_pseudogenes_file_list,logger,cores)
	#global definition of specific_dup_list for specific gene expression analysis in the absence of a reference organism
	specific_dup_list = []
	# if the user has specified a single or list of desired genes
	if len(ref_free_desired_genes) != 0 and ref_name == 'NA':
		desired_gene_names = []
		for every in gff3_input_file:
			gff_name = get_basename(str(every))[0]
			if org_type == 'eukaryote':
				for pep in trans_peptides_file_list:
					pepname = get_basename(str(pep))[0]
					if gff_name == pepname:
						if 'all' in gff_config_parameters:
							child_attribute = gff_config_parameters['all']['child_attribute']
							child_parent_linker = gff_config_parameters['all']['child_parent_linker']
							parent_attribute = gff_config_parameters['all']['parent_attribute']
						elif 'default' in gff_config_parameters:
							child_attribute = gff_config_parameters['default']['child_attribute']
							child_parent_linker = gff_config_parameters['default']['child_parent_linker']
							parent_attribute = gff_config_parameters['default']['parent_attribute']
						else:
							child_attribute = gff_config_parameters[ref_name]['child_attribute']
							child_parent_linker = gff_config_parameters[ref_name]['child_parent_linker']
							parent_attribute = gff_config_parameters[ref_name]['parent_attribute']
						transcripts_per_gene, gene_names = tpmclean(str(every), str(pep), child_parent_linker,child_attribute)
						# Step 1: Replace transcript names with gene names
						# Create a mapping of transcript to gene for quick lookup
						transcript_to_gene = {transcript: gene for gene, transcripts in transcripts_per_gene.items() for transcript in transcripts}
						for transcript in ref_free_desired_genes:
							gene_name = transcript_to_gene.get(transcript,transcript)  # Default to transcript name if not found
							desired_gene_names.append(gene_name)
						ref_free_desired_genes_set = set(desired_gene_names)
			elif org_type == 'prokaryote':
				for pep in peptides_file_list:
					pepname = (get_basename(str(pep))[0]).strip().split('no_alt_trans')[0]
					if gff_name == pepname:
						if 'all' in gff_config_parameters:
							child_attribute = gff_config_parameters['all']['child_attribute']
							child_parent_linker = gff_config_parameters['all']['child_parent_linker']
							parent_attribute = gff_config_parameters['all']['parent_attribute']
						elif 'default' in gff_config_parameters:
							child_attribute = gff_config_parameters['default']['child_attribute']
							child_parent_linker = gff_config_parameters['default']['child_parent_linker']
							parent_attribute = gff_config_parameters['default']['parent_attribute']
						else:
							child_attribute = gff_config_parameters[ref_name]['child_attribute']
							child_parent_linker = gff_config_parameters[ref_name]['child_parent_linker']
							parent_attribute = gff_config_parameters[ref_name]['parent_attribute']
						transcripts_per_gene, gene_names = tpmclean(str(every), str(pep), child_parent_linker,child_attribute)
						# Step 1: Replace transcript names with gene names
						# Create a mapping of transcript to gene for quick lookup
						transcript_to_gene = {transcript: gene for gene, transcripts in transcripts_per_gene.items() for transcript in transcripts}
						for transcript in ref_free_desired_genes:
							gene_name = transcript_to_gene.get(transcript,transcript)  # Default to transcript name if not found
							desired_gene_names.append(gene_name)
						ref_free_desired_genes_set = set(desired_gene_names)
		specific_dup_list = list(ref_free_desired_genes_set)

		is_empty = len(os.listdir(validated_tpm_files)) == 0
		if is_empty == False:
			# cleaning the TPM files
			if org_type == 'eukaryote':
				validated_tpm_input_file = sorted(glob.glob(os.path.join(validated_tpm_files, "*.txt")) + glob.glob(os.path.join(validated_tpm_files, "*.txt.gz")) + glob.glob(os.path.join(validated_tpm_files, "*.tpms.txt")) + glob.glob(os.path.join(validated_tpm_files, "*.tpms.txt.gz")))
				logger.info('TPM files have been cleaned')
			elif org_type == 'prokaryote':
				validated_tpm_input_file = exp_input_file
			ref_free_master_folder = os.path.join(output_dir,'Specific_ref_free_duplicates_analysis')  # Creating the path for the folder for storing the data of particular organism
			os.makedirs(ref_free_master_folder,exist_ok=True)  # Creates an organism specific folder in the specified path location
			parallelize_rthresh_pseudogenes(ref_free_master_folder, validated_tpm_input_file, rand_genes, dpi_no, cores,logger, avg_exp_cutoff, rthreshold_pseudogenes_dic_dir)
			rthreshold_pseudogenes_dic_dir_path = Path(rthreshold_pseudogenes_dic_dir)
			rthreshold_pseudogenes_file_list = list(rthreshold_pseudogenes_dic_dir_path.iterdir())
			singleton_files=[]
			for root, dirs, files in os.walk(singleton_dir_final):
				for file in files:
					singleton_files.append(os.path.join(root,file))
			tandems_analysis = []
			for root, dirs, files in os.walk(tandems_dir):
				for file in files:
					tandems_analysis.append(os.path.join(root, file))
			duplicates_analysis = tandems_analysis
			two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no, 'Tandem',rthreshold_pseudogenes_file_list, cores, logger,specific_dup_list,ref_free_master_folder,singleton_files,ref_name)
			proximals_analysis = []
			for root, dirs, files in os.walk(proximals_dir):
				for file in files:
					proximals_analysis.append(os.path.join(root, file))
			duplicates_analysis = proximals_analysis
			two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no, 'Proximal',rthreshold_pseudogenes_file_list, cores, logger,specific_dup_list,ref_free_master_folder,singleton_files,ref_name)
			if classify == 'strict':
				mixed_analysis = []
				for root, dirs, files in os.walk(mixed_dir):
					for file in files:
						mixed_analysis.append(os.path.join(root, file))
				duplicates_analysis = mixed_analysis
				two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no,'Mixed', rthreshold_pseudogenes_file_list, cores, logger,specific_dup_list,ref_free_master_folder,singleton_files,ref_name)
			elif classify == 'overlap':
				pass
			dispersed_analysis = []
			for root, dirs, files in os.walk(dispersed_dir):
				for file in files:
					dispersed_analysis.append(os.path.join(root, file))
			duplicates_analysis = dispersed_analysis
			two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no,'Dispersed', rthreshold_pseudogenes_file_list, cores, logger,specific_dup_list,ref_free_master_folder,singleton_files,ref_name)

	if evo_analysis == 'yes':
		if org_type == 'eukaryote':
			peptides_folder_path = Path(no_alt_trans_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
			cds_folder_path = Path(no_alt_trans_cds_dir)
			cds_file_list = list(cds_folder_path.iterdir())
		elif org_type == 'prokaryote':
			peptides_folder_path = Path(pep_dir)
			peptides_file_list = list(peptides_folder_path.iterdir())
			cds_folder_path = Path(clean_cds_dir)
			cds_file_list = list(cds_folder_path.iterdir())
		tandems_folder_path = Path(tandems_dir)
		tandems_file_list = list(tandems_folder_path.iterdir())
		proximals_folder_path = Path(proximals_dir)
		proximals_file_list = list(proximals_folder_path.iterdir())
		dispersed_folder_path = Path(dispersed_dir)
		dispersed_file_list = list(dispersed_folder_path.iterdir())
		mixed_folder_path = Path(mixed_dir)
		mixed_file_list = list(mixed_folder_path.iterdir())
		time_kaks_start = time.perf_counter()
		parallelize_kaks_calculations(cds_file_list, peptides_file_list, tandems_file_list, proximals_file_list, dispersed_file_list, mixed_file_list, mafft,cores,classify,logger,kaks_dir,ref_name,method,org_type,dpi_no)
		time_kaks_end = time.perf_counter()
		print(f"time taken for ka/ks analyses: {time_kaks_end - time_kaks_start}")

	#Performing statistical analyses and plotting pairwise gene expression plots of the gene duplicates
	if analysis == 'no':
		pass
	elif analysis == 'yes':
		logger.info('Starting gene duplicates expression analysis')
		is_empty = len(os.listdir(validated_tpm_files)) == 0
		if is_empty == False:
			#cleaning the TPM files
			if org_type == 'eukaryote':
				validated_tpm_input_file = sorted(glob.glob(os.path.join(validated_tpm_files, "*.txt")) + glob.glob(os.path.join(validated_tpm_files, "*.txt.gz"))+ glob.glob(os.path.join(validated_tpm_files, "*.tpms.txt")) + glob.glob(os.path.join(validated_tpm_files, "*.tpms.txt.gz")))
			elif org_type == 'prokaryote':
				validated_tpm_input_file = exp_input_file
			parallelize_rthresh_pseudogenes (analyses,validated_tpm_input_file,rand_genes,dpi_no,cores,logger,avg_exp_cutoff,rthreshold_pseudogenes_dic_dir)
			rthreshold_pseudogenes_dic_dir_path = Path(rthreshold_pseudogenes_dic_dir)
			rthreshold_pseudogenes_file_list = list(rthreshold_pseudogenes_dic_dir_path.iterdir())
			singleton_files=[]
			tandems_analysis = []
			for root, dirs, files in os.walk(tandems_dir):
				for file in files:
					tandems_analysis.append(os.path.join(root, file))
			duplicates_analysis = tandems_analysis
			two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no, 'Tandem',rthreshold_pseudogenes_file_list,cores,logger,specific_dup_list,rthreshold_pseudogenes_dic_dir,singleton_files,ref_name)
			proximals_analysis = []
			for root, dirs, files in os.walk(proximals_dir):
				for file in files:
					proximals_analysis.append(os.path.join(root, file))
			duplicates_analysis = proximals_analysis
			two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no, 'Proximal',rthreshold_pseudogenes_file_list,cores,logger,specific_dup_list,rthreshold_pseudogenes_dic_dir,singleton_files,ref_name)
			if classify == 'strict':
				mixed_analysis = []
				for root, dirs, files in os.walk(mixed_dir):
					for file in files:
						mixed_analysis.append(os.path.join(root, file))
				duplicates_analysis = mixed_analysis
				two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no,'Mixed', rthreshold_pseudogenes_file_list,cores,logger,specific_dup_list,rthreshold_pseudogenes_dic_dir,singleton_files,ref_name)
			elif classify == 'overlap':
				pass
			if dispersed_analysis_choice == 'no':
				pass
			elif dispersed_analysis_choice == 'yes':
				dispersed_analysis = []
				for root, dirs, files in os.walk(dispersed_dir):
					for file in files:
						dispersed_analysis.append(os.path.join(root, file))
				duplicates_analysis = dispersed_analysis
				two_step_parallelization_master(duplicates_analysis, validated_tpm_input_file, analyses, rand_genes, dpi_no, 'Dispersed',rthreshold_pseudogenes_file_list,cores,logger,specific_dup_list,rthreshold_pseudogenes_dic_dir,singleton_files,ref_name)

	t_end = time.perf_counter()
	logger.info(f"Time taken for the Dupylicate analysis: {t_end - t_start} seconds" + '\n')
	logger.info("Your Dupylicate analysis is complete!" + '\n')
	# Check if Tmp folder exists first and then delete
	if clean_remnants == 'yes':
		if os.path.exists(tmp_dir):
			shutil.rmtree(tmp_dir)
		if analysis == 'no':
			if os.path.exists(analyses):
				shutil.rmtree(analyses)
		if evo_analysis == 'no':
			if os.path.exists(kaks_dir):
				shutil.rmtree(kaks_dir)
		if not os.listdir(unclassified_dir_final):
			os.rmdir(unclassified_dir_final)
	#code to remove empty subfolders without plots after the duplicates or specific duplicates analysis
	if analysis == 'yes':
		for root, dirs, files in os.walk(analyses, topdown=False):
			for d in dirs:
				full_path = os.path.join(root, d)
				# Check if folder is empty
				if not os.listdir(full_path):
					os.rmdir(full_path)
	if specific_analysis == 'yes':
		pass
		for root, dirs, files in os.walk(master_folder, topdown=False):
			for d in dirs:
				full_path = os.path.join(root, d)
				#Check if folder is empty
				if not os.listdir(full_path):
					os.rmdir(full_path)

#script usage possibilities in reference-free mode
if '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv and '--specific_genes' not in sys.argv and '--specific_duplicates_analysis' not in sys.argv:
	main( sys.argv )
elif '--gff' in sys.argv and '--fasta' in sys.argv and '--to_annotate' in sys.argv and '--gemoma' in sys.argv and '--out' in sys.argv and '--specific_genes' not in sys.argv and '--specific_duplicates_analysis' not in sys.argv:
	main(sys.argv)
elif '--gff' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--specific_genes' not in sys.argv and '--specific_duplicates_analysis' not in sys.argv:
	main( sys.argv )
elif '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' not in sys.argv and '--specific_duplicates_analysis' not in sys.argv:
	main(sys.argv)
elif '--gff' in sys.argv and '--fasta' in sys.argv and '--to_annotate' in sys.argv and '--gemoma' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' not in sys.argv and '--specific_duplicates_analysis' not in sys.argv:
	main(sys.argv)
elif '--gff' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' not in sys.argv and '--specific_duplicates_analysis' not in sys.argv:
	main(sys.argv)
elif '--gff' in sys.argv and '--pep' in sys.argv and '--out' in sys.argv and '--ka_ks' not in sys.argv and '--specific_genes' not in sys.argv and '--specific_duplicates_analysis' not in sys.argv:
	main( sys.argv )

#script usage possibilities in the presence of a reference
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv and '--specific_genes' in sys.argv:
	main( sys.argv )
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--to_annotate' in sys.argv and '--gemoma' in sys.argv and '--out' in sys.argv and '--specific_genes' in sys.argv:
	main( sys.argv )
elif '--ref' in sys.argv and '--gff' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--specific_genes' in sys.argv:
	main( sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' in sys.argv:
	main(sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--to_annotate' in sys.argv and '--gemoma' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' in sys.argv:
	main(sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' in sys.argv:
	main(sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--pep' in sys.argv and '--out' in sys.argv and '--ka_ks' not in sys.argv and '--specific_genes' in sys.argv:
	main( sys.argv )

#script usage possibilities in the presence of a reference; since specific_duplicates_analysis can be carried out only when specific_genes are provided, this set of possibilities is included
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv and '--specific_genes' in sys.argv and '--specific_duplicates_analysis' in sys.argv:
	main( sys.argv )
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--to_annotate' in sys.argv and '--gemoma' in sys.argv and '--out' in sys.argv and '--specific_genes' in sys.argv and '--specific_duplicates_analysis' in sys.argv:
	main( sys.argv )
elif '--ref' in sys.argv and '--gff' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--specific_genes' in sys.argv and '--specific_duplicates_analysis' in sys.argv:
	main( sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' in sys.argv and '--specific_duplicates_analysis' in sys.argv:
	main(sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--fasta' in sys.argv and '--to_annotate' in sys.argv and '--gemoma' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' in sys.argv and '--specific_duplicates_analysis' in sys.argv:
	main(sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv and '--ka_ks' in sys.argv and '--specific_genes' in sys.argv and '--specific_duplicates_analysis' in sys.argv:
	main(sys.argv)
elif '--ref' in sys.argv and '--gff' in sys.argv and '--pep' in sys.argv and '--out' in sys.argv and '--ka_ks' not in sys.argv and '--specific_genes' in sys.argv and '--specific_duplicates_analysis' in sys.argv:
	main( sys.argv )

#if none of the above possibilities are satisfied - exit the script and display script usage 
else:
	sys.exit( __usage__ )
