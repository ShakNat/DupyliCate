# DupyliCate
Light weight python tool for mining, classification and analysis of gene duplications.

<img width="3933" height="2583" alt="Dupylicate_github_schematic" src="https://github.com/user-attachments/assets/cbae189d-19d2-45a3-b75f-fb0dbc5c4a25" />

## Workflow
<img width="4638" height="5484" alt="Dupylicate_workflow_diagram" src="https://github.com/user-attachments/assets/a818ed9f-530f-4b54-bf4a-8cc0fa5f2ce7" />

## Gene duplicates types output by DupyliCate
<img width="3120" height="1851" alt="Gene_duplicates_types" src="https://github.com/user-attachments/assets/ca081760-2e30-4eea-bfb3-94877b443d6f" />

```
Usage:

  python3 Dupylicate.py

MANDATORY:

	--gff <full path to folder containing GFF files>

	--fasta <full path to folder containing WGS FASTA files> | --cds <full path to folder containing CDS FASTA files> | --pep <full path to folder containing PEP files>

	--out <full path to output folder>

OPTIONAL:

	--ref <name of the organism to be taken as reference> Default is NA and the script runs in reference-free mode

	--prokaryote <use the flag if the organisms for analysis are prokaryotes>

	--pseudos_gff <yes | no for optional inclusion or exclusion of pseudogenes with coding features in the GFF file> Default is no

	--mode <overlap | strict modes; in the overlap mode genes repeat among the different duplicates classification; 
		in the strict mode there is no gene repetition and you have a new classification group called mixed 
		duplicates containing the related connected components from the other three duplicate classes> DEFAULT is overlap

	--to_annotate <full path to file containing names of queries and the corresponding reference organism for GeMoMa annotation separated by comma - Query,Reference -> one pair per line>

	--seq_aligner <choose one among blast | diamond | mmseqs2 > DEFAULT is blast

	--blast <full path to BLAST if not already in yopur PATH environment variable>

	--diamond <full path to diamond if not already in yopur PATH environment variable>

	--mmseqs <full path to mmseqs2 if not already in yopur PATH environment variable>

	--mafft <full path to MAFFT if not already in yopur PATH environment variable>

	--evalue <evalue for alignment> DEFAULT is 1e-5

	--gemoma <full path to GeMoMa if not already in yopur PATH environment variable>

	--qc < yes | no for quality check with BUSCO> Default is no

	--busco  <full path to BUSCO or choose 'busco_docker' for docker-based BUSCO installtion> Default is 'busco'

	--busco_version <BUSCO version in the format vx.x.x - needed only if you have docker-based BUSCO installation> Deafult is v5.8.2

	-- container_version <docker container version of BUSCO> Default is cv1

	--docker_host_path < full host folder path - needed for docker-based BUSCO installation>

	--docker_container_path < full mount path in the docker container - needed for docker-based BUSCO installation>

	--score <auto for automatic threshold finding and <float number between 0 and 1> for manual threshold finding for segregating singletons and duplicates> DEFAULT is auto

	--self_simcut <float similarity percentage to remove self alignment hits with low similarity percentage> Default is 50.0

	--hits <number of top hits to be considered for finding suitable orthologs in the reference organism> Default is 10

	--occupancy <float number cutoff for MAFFT aligned file trimming> Default is 0.1

	--scoreratio < float ratio of forward alignment bit score and self alignment bit score of query to assess if the forward hit is valid> Default is 0.3

	--fwd_simcut < float similarity percentage to remove forward alignment hits with low similarity percentage> Default is 40.0

	--cores <number of cores needed to run Dupylicate analysis> DEFAULT is 4

	--proximity <integer value for the number of intervening genes to detect proximal duplications> DEFAULT is 10

	--synteny_score <float value which is used as a cut-off or threshold for synteny analysis> DEFAULT is 0.5

	--flank <integer value specifying the number of flanking genes to be considered to determine the synteny window size in synteny analysis>

	--side <integer value for synteny support from either side of a flanking region of a synteny window>

	--ka_ks <yes | no to calculate ka, ks values> DEFAULT is 'no'

	--ka_ks_method < MYN | NG methods for Ka/Ks ratio caclulation> Default is NG

	--duplicates_analysis <yes | no for further statistical analysis of identified gene duplicates> DEFAULT is no

	--specific_duplicates_analysis <yes | no for further statistical analysis of specified ref genes' gene duplicates> DEFAULT is no

	--dpi <resolution value desired for plots low | moderate | high | very high > DEFAULT is moderate

	--analyse_disperse_duplicates <yes | no for statistical analysis of dispersed gene duplicates> DEFAULT is no

	--exp <full path to the folder with expression counts files>

	--avg_exp_cut < float value cutoff for average expression value across samples for classifying a gene as pseudogene based on gene expression> Default is 1

	--genes <integer value referring to the number of genes to be taken for random pairing to determine threshold 
		value for assessing functional divergence of gene duplicates based on gene expression values> DEFAULT is 10000

	--specific_genes <user given list of specific genes from the reference with one transcript name per line, for analysing gene duplication outputs>

	--clean up <yes | no cleans up intermediate files/ folders> DEFAULT is YES

ALLOWED FILE EXTENSIONS:

	<'.gff', '.gff.gz', '.gff3', '.gff3.gz', '.fa', '.fa.gz', '.fna', '.fna.gz','.fasta', '.fasta.gz',
	'.genome.fasta', '.genome.fa', '.genome.fasta.gz', '.genome.fa.gz', '.genome.fna',
	'.genome.fna.gz','.cds.fa', '.cds.fasta','.cds.fa.gz', '.cds.fasta.gz', '.cds.fna',
	'.cds.fna.gz','.pep.fa','.pep.fa.gz','.pep.fasta','.pep.fasta.gz', '.pep.fna',
	'.pep.fna.gz','.tsv','.txt','.txt.gz','.tpms.txt','.tpms.txt.gz'>

```

```
Usage:

	python3 AGAT_wrapper.py

MANDATORY:

	--gff_in <full path to GFF3 file or folder containing GFF3 annotation files>

	--gff_out <full path to output folder or output file>

OPTIONAL:

	--agat <full path to the agat_convert_sp_gxf2gxf.pl script including the script name>

```
	
