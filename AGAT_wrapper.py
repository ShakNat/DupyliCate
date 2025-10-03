### --- AGAT_wrapper.py --- ###
### --- v0.1 --- ###
### Shakunthala Natrajan ###
### Boas Pucker ###
### --- imports --- ###
import subprocess
import os
import shutil
import re, sys, gzip, glob
import concurrent.futures
import multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
try:
	from tqdm import tqdm
	tqdm_available = True
except ImportError:
	tqdm_available = False
### --- end of imports --- ###
__usage__ = """
			python3 AGAT_wrapper.py --gff_in <input_gff_folder> --gff_out <output_gff_folder>
			--gff_in <full path to GFF3 file or folder containing GFF3 annotation files>
			--gff_out <full path to output folder or output file>
			
			Optional:
			--agat <full path to the agat script agat_convert_sp_gxf2gxf.pl including the script name>
			"""

def agat_convert (agat,input_gff,output_gff,output_dir):
	cmd = f'{agat} -g {input_gff} -o {output_gff}'
	result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
	if result.returncode != 0:
		raise RuntimeError(f"AGAT conversion failed: {result.stderr}")
	else:
		out_base = os.path.basename(output_gff)
		final_out_base = str(out_base).replace('agat_','')
		final_out_base = str(final_out_base).replace('.gz', '')
		final_output_gff = os.path.join(output_dir, final_out_base)
		with open(output_gff, 'r') as f:
			with open(final_output_gff, 'w') as out:
				line = f.readline()
				while line:
					if line[0] != '#':
						if 'agat-' in line:
							modified_line = str(line).replace('agat-', '')
							out.write(modified_line)
						else:
							out.write(line)
					else:
						out.write(line)
					line = f.readline()
		os.remove(output_gff)

def parallelize_agat_processing(agat,gff_file_list,output_dir, cores):
	futures = []
	num_files = len(gff_file_list)
	files_to_parallelize = min(len(gff_file_list),cores)
	with ProcessPoolExecutor(max_workers=files_to_parallelize) as executor:
		for file in gff_file_list:
			base_file = os.path.basename(file)
			new_base_file = 'agat_'+base_file
			output_gff = os.path.join(output_dir,new_base_file)
			input_gff = file
			futures.append(executor.submit(agat_convert,agat,input_gff,output_gff,output_dir))
		# Wrap the iterator conditionally
		iterator = as_completed(futures)
		if tqdm_available:
			iterator = tqdm(iterator, total=len(futures), desc="Organisms processed by AGAT")
		# Add tqdm to monitor futures
		for future in iterator:
			try:
				result = future.result()
				# Wait and capture exceptions
			except Exception as e:
				sys.stdout.write(f"Worker crashed with error: {e}")
				sys.stdout.write("The error is as follows:")
				raise
		sys.stdout.write(f"Processing of GFF files complete.")

def main(arguments):
	if '--agat' in arguments:
		agat = arguments[arguments.index('--agat')+1] # full path tp the agat_convert_sp_gxf2gxf.pl script
	else:
		agat = 'agat_convert_sp_gxf2gxf.pl'

	if '--gff_in' in arguments:
		gff3_folder = arguments[arguments.index('--gff_in') + 1]  # full path to GFF3 file or folder containing GFF3 annotation files
		if gff3_folder[-1] == "/":
			gff3_input_file = sorted(glob.glob(gff3_folder + "*.gff") + glob.glob(gff3_folder + "*.gff3") + glob.glob(gff3_folder + "*.gff.gz") + glob.glob(gff3_folder + "*.gff3.gz"))
		else:
			gff3_input_file = [gff3_folder]

	if '--cores' in arguments:
		cores = arguments[arguments.index('--cores')+1]
	else:
		cores = 4

	output_dir = arguments[arguments.index('--gff_out') + 1]  # full path to output folder
	if output_dir[-1] != "/":
		output_dir += "/"
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	parallelize_agat_processing(agat, gff3_input_file, output_dir, cores)

if '--gff_in' in sys.argv and '--gff_out' in sys.argv:
	main( sys.argv )
else:
	sys.exit(__usage__)



