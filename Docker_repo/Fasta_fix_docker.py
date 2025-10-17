#v0.1
### Shakunthala Natrajan ###
### Boas Pucker ###

__usage__ ="""
			python3 Fasta_fix.py

			MANDATORY:
			--in <full path to a folder containing FASTA files or a single FASTA file>
			--out <full path to output directory>
			--config <full path to config file including the config file name>

			OPTIONAL:
			--cores <number of cores for running the script>

			Config file preparation instructions:
			1). The config file is a simple .txt file

			2). Specify the organism name i.e. the base name without the extension of your file in the first column

			3). If the same set of string manipulation operations are to be performed for all the files in a folder,
				 just specify the word all in the first column and follow it up in the next columns with the desired
				 operations - This single line is enough if the same set of operations are to be performed on all the files

			4). Specify the various string manipulation operations to be performed on the FASTA header of this particular
				organism's file in the subsequent columns separated by white space or tab

			5) IMPORTANT: The operations will be performed in the same order as you specify them in the config wise i.e.
				column wise order of the different operations you specify will be followed by the script. 
				Hence operation ORDER is IMPORTANT

			6) Possible operations and the manner in which they are to be specified are as follows
				i extract: <to extract text of a particular attribue in the header>
				example- extract:ID=:; <This tells extract the text following the attribute ID= in the header
									   until the delimiter ; is encountered><If you specify no delimiter, a
									   default list of delimiters will be searched automatically by the script>

				ii. split: <specify a character at which splitting must be done>
				example- split:_ <This will split the text at undersore>

				iii. take: <specify a single position or multiple positions separated by commas to be taken after splitting>
				example- take:1,3 <This will take the first and third elements after splitting>

				iv. join: <join a list of string elements using a specific character>
				example- join:- <This will join a list of strings into a single string separated by - >

				v. remove: <to remove a particular pattern from the text>
				example- remove:v2 <This will remove v2 from the text>

				vi. add_prefix: <to add a prefix to the text>
				   add_suffix: <to add a suffix to the text>
				   example: add_suffix:ath <This will add the text ath to the text>

				vii. replace: <to replace a specfic character with another character or pattern; Specify the pattern
							  to be replaced first and the pattern to be added in its place - separate these two 
							  with a comma>
				example- replace:phytozome,phyt <This will replace phytozome with phyt>

				viii. uppercase <to capitalize the text fully>
					 lowercase <to write the full text in lower case>
			"""


##### imports #####
import re,sys,os
os.environ["PATH"] = "/tools/dupy_env/bin:" + os.environ.get("PATH", "")
import glob
import gzip
import concurrent.futures
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
##### end of imports #####

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

#function to clean FASTA headers and write the headers and sequences to a new FASTA file
def process_fasta_headers(file, operations, output_dir, filename, extension):
	sequences = {}
	seq_lines = []  # Initialize before the loop
	current_header = None  # Track current header

	if file[-2:].lower() != 'gz':  # uncompressed FASTA file
		with open(file, 'r') as f:
			for line in f:
				line = line.strip()  # Remove whitespace
				if line.startswith(">"):  # Better check for header
					# Store previous sequence before processing new header
					if current_header:
						sequences[current_header] = ''.join(seq_lines)

					# Process new header
					header = line[1:]  # Remove ">"

					# Apply all operations
					for operation in operations:
						if operation.startswith("extract:"):
							extract_params = operation[8:]
							if ':' not in extract_params:
								search_term = extract_params
								pattern = f"{re.escape(search_term)}([^\\s|,:;()\\[\\]{{}}/<>=@#~!&%^*+?]+)"
							else:
								parts = extract_params.split(':', 1)
								search_term = parts[0]
								stop_delimiter = parts[1]
								if not stop_delimiter:
									pattern = f"{re.escape(search_term)}(.+)"
								else:
									escaped_delimiter = re.escape(stop_delimiter)
									pattern = f"{re.escape(search_term)}([^{escaped_delimiter}]+)"
							match = re.search(pattern, header)
							header = match.group(1) if match else header
						elif operation.startswith("split:"):
							split_params = operation[6:]
							header = header.split(split_params)
						elif operation.startswith("take:"):
							indices = operation[5:]
							if isinstance(header, list):
								if ',' in indices:
									idx_list = [(int(i) - 1) for i in indices.split(',')]
									header = [header[i] for i in idx_list if 0 <= i < len(header)]
								else:
									idx = (int(indices) - 1)
									header = header[idx] if 0 <= idx < len(header) else header
						elif operation.startswith("join:"):
							separator = operation[5:]
							if isinstance(header, list):
								header = separator.join(header)
						elif operation.startswith("remove:"):
							pattern = operation[7:]
							header = re.sub(pattern, "", str(header))
						elif operation.startswith("add_prefix:"):
							pattern = operation[11:]
							header = str(pattern) + str(header)
						elif operation.startswith("add_suffix:"):
							pattern = operation[11:]
							header = str(header) + str(pattern)
						elif operation.startswith("replace:"):
							pattern = operation[8:].split(',')
							to_replace = pattern[0]
							replace_with = pattern[1]
							header = re.sub(to_replace, replace_with, str(header))
						elif operation.startswith("uppercase"):
							header = str(header).upper()
						elif operation.startswith("lowercase"):
							header = str(header).lower()

					# Set current header and reset sequence lines
					# Handle case where header is still a list
					if isinstance(header, list):
						current_header = ''.join(header).strip()
					else:
						current_header = str(header).strip()
					seq_lines = []
				else:
					# Accumulate sequence lines
					seq_lines.append(line)

			# Store the last sequence after loop ends
			if current_header:
				sequences[current_header] = ''.join(seq_lines)

	else:  # compressed FASTA file
		with gzip.open(file, 'rt') as f:
			for line in f:
				line = line.strip()
				if line.startswith(">"):
					# Store previous sequence
					if current_header:
						sequences[current_header] = ''.join(seq_lines)

					header = line[1:]

					# Apply all operations (same as above)
					for operation in operations:
						if operation.startswith("extract:"):
							extract_params = operation[8:]
							if ':' not in extract_params:
								search_term = extract_params
								pattern = f"{re.escape(search_term)}([^\\s|,:;()\\[\\]{{}}/<>=@#~!&%^*+?]+)"
							else:
								parts = extract_params.split(':', 1)
								search_term = parts[0]
								stop_delimiter = parts[1]
								if not stop_delimiter:
									pattern = f"{re.escape(search_term)}(.+)"
								else:
									escaped_delimiter = re.escape(stop_delimiter)
									pattern = f"{re.escape(search_term)}([^{escaped_delimiter}]+)"
							match = re.search(pattern, header)
							header = match.group(1) if match else header
						elif operation.startswith("split:"):
							split_params = operation[6:]
							header = header.split(split_params)
						elif operation.startswith("take:"):
							indices = operation[5:]
							if isinstance(header, list):
								if ',' in indices:
									idx_list = [(int(i) - 1) for i in indices.split(',')]
									header = [header[i] for i in idx_list if 0 <= i < len(header)]
								else:
									idx = (int(indices) - 1)
									header = header[idx] if 0 <= idx < len(header) else header
						elif operation.startswith("join:"):
							separator = operation[5:]
							if isinstance(header, list):
								header = separator.join(header)
						elif operation.startswith("remove:"):
							pattern = operation[7:]
							header = re.sub(pattern, "", str(header))
						elif operation.startswith("add_prefix:"):
							pattern = operation[11:]
							header = str(pattern) + str(header)
						elif operation.startswith("add_suffix:"):
							pattern = operation[11:]
							header = str(header) + str(pattern)
						elif operation.startswith("replace:"):
							pattern = operation[8:].split(',')
							to_replace = pattern[0]
							replace_with = pattern[1]
							header = re.sub(to_replace, replace_with, str(header))
						elif operation.startswith("uppercase"):
							header = str(header).upper()
						elif operation.startswith("lowercase"):
							header = str(header).lower()

					# Set current header and reset sequence lines
					# Handle case where header is still a list
					if isinstance(header, list):
						current_header = ''.join(header).strip()
					else:
						current_header = str(header).strip()
					seq_lines = []
				else:
					seq_lines.append(line)

			# Store the last sequence
			if current_header:
				sequences[current_header] = ''.join(seq_lines)

	# Write output
	if sequences:
		output_file = os.path.join(output_dir, str(filename) + str(extension))
		if output_file[-2:].lower() != 'gz':
			with open(output_file, 'w') as out:
				for key, value in sequences.items():
					out.write('>' + str(key) + '\n')
					out.write(str(value) + '\n')
		else:
			with gzip.open(output_file, 'wt') as out:
				for key, value in sequences.items():
					out.write('>' + str(key) + '\n')
					out.write(str(value) + '\n')

#master parallelization function
def parallelize_header_fix(fasta_input_file,cores,config_dic,output_dir):
	futures = []
	num_files = len(fasta_input_file)
	files_to_parallelize = min(len(fasta_input_file), cores)
	with ProcessPoolExecutor(max_workers=files_to_parallelize) as executor:
		if len(config_dic.keys())==1 and 'all' in config_dic.keys(): # if all the fasta files in a folder have the same string manipulation specify the column wise operations by specifying all in the first column
			for file in fasta_input_file:
				filename = get_basename(str(file))[0]
				extension = os.path.basename(file).replace(filename, "")
				futures.append(executor.submit(process_fasta_headers, file, config_dic['all'], output_dir, filename, extension))
		else:
			for file in fasta_input_file:
				filename = get_basename(str(file))[0]
				extension = os.path.basename(file).replace(filename, "")
				if filename in config_dic:
					futures.append(executor.submit(process_fasta_headers,file,config_dic[filename],output_dir,filename,extension))

def main(arguments):
	if '--in' in arguments:
		fasta_folder = arguments[arguments.index('--in') + 1]  # full path to GFF3 file or folder containing GFF3 annotation files
		if fasta_folder[-1] == '/':
			pattern = fasta_folder + "*{,.cds,.pep}*.{fa,fasta,fna}{,.gz}"
			# Use glob.glob with the pattern
			fasta_input_file = sorted(glob.glob(pattern, recursive=False))
		else:
			fasta_input_file = [fasta_folder]
	output_dir = arguments[arguments.index('--out') + 1]  # full path to output folder
	if output_dir[-1] != "/":
		output_dir += "/"
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	if '--config' in arguments:
		config_file = arguments[arguments.index('--config')+1] # full path to config file
	if '--cores' in arguments:
		cores = arguments[arguments.index('--cores')+1]
	else:
		cores = 4
	config_dic = {}
	with open (config_file,'r')as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			parts = line.split()
			if len(parts)>=2:
				config_dic[parts[0]] = parts[1:]
			else:
				sys.stdout.write("Incorrect or insufficient config parameters. Exiting.")
				sys.exit()
	if config_dic:
		parallelize_header_fix(fasta_input_file, cores, config_dic,output_dir)

if '--in' in sys.argv and '--out' in sys.argv and '--config' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)
