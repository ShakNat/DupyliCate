# DupyliCate

<img width="3933" height="2583" alt="Dupylicate_github_schematic" src="https://github.com/user-attachments/assets/cbae189d-19d2-45a3-b75f-fb0dbc5c4a25" />


<div align="justify">
	
DupyliCate is a light weight python tool for mining, classification and analysis of gene duplications. It can help find and classify gene duplicates in a number of organisms concurrently and is designed to have high throughput. It can be used with a reference organism for comparative analysis or in a reference-free manner for intra-species gene duplicates identification. The tool moreover can be used with a reference organism to obtain a comparative gene table of specific genes in the reference, whose copy number variation that the user wants to know in other sample organisms, thus facilitating in-depth comparative genomic analyses.


There are two main modes of DupyliCate - 


**'Overlap'**: Gene duplicates repeat across groups and this mode also produces an output file called 'Duplicates_relationships' that connects the different genes repeating across the duplicate classes - tandems, proximals, dispersed

**'Strict'** : Gene duplicates do not repeat across groups and this mode produces a fourth class of duplicates type called 'Mixed_duplicates' where the gene duplicates from different duplicate classes - tandems, proximals, dispersed are merged together to retain their individual classification an also their relationship across the groups

DupyliCate also facilitates optional gene expression analysis of gene duplicates using count table data. Further, it provides the user, option to calculate Ka/Ks values of the duplicates. 
</div>

## Workflow
<img width="4638" height="5484" alt="Dupylicate_workflow_diagram" src="https://github.com/user-attachments/assets/a818ed9f-530f-4b54-bf4a-8cc0fa5f2ce7" />

<div align="justify">
	
(1) DupyliCate needs the structural annotation file(s) (GFF3) along with one of assembly or coding sequence or peptide sequence FASTA file(s). In case, structural annotation is not available for a sample organism in which duplicates need to be analysed, there is an option to provide the annotation of a related organism's annotation as reference. This will then be used to produce the required annotation and carry out further analysis. (**Step 1**)
</div>

<div align="justify">
	
(2) The input files are first checked and validated before the actual run starts. If, the check fails, the script exits and the errors will be recorded and displayed organism-wise in the path Tmp/Errors_warnings. In case of GFF errors, there is a helper script provided along with the main script that uses AGAT to process and correct the GFF files. The corrected and validated files then enter the main analysis and processed to give PEP files without alternate transcripts. (**Steps 2,3**)
</div>

<div align="justify">
	
(3) Since the output depends on the quality of the input files and is also influenced by the ploidy of the organisms being analysed, a BUSCO-based QC step is included. This provides detailed QC reports containing the BUSCO completeness, duplication and BUSCO-derived pseudo-ploidy number. (**Step 4**)
</div>

<div align="justify">
	
(4) Before moving on to the duplication analysis, it is important to segregate singletons and duplicates correctly. For this, two cut-offs - one based on normalized bit score and self-similarity are offered. A more detailed information on the cut-offs and parameters can be obtained here. By default, a BUSCO-based auto threshold method is chosen, where BUSCO single copy genes are used to identify the normalized bit score threshold to segregate singletons and duplicates. If BUSCO is not available or if the organism has a very low number of BUSCO single copy genes, then the default fallback is to go for a default self-similarity threshold instead of normalized bit score threshold. There is an option to manually set the normalized bit score threshold and self-similarity threshold as well. (**Step 5**)
</div>

<div align="justify">
	
(5) Self alignment of sample organisms is performed and if a reference organism is provided, forward alignment of each sample is performed against this reference organism. This is followed by a comprehensive ortholog detection step, where a synergy of local, global alignments and phylogeny is utilized to obtain orthologs of genes in the samples against the reference genes. (**Steps 5,6**)
</div>

<div align="justify">
	
(6) Next, the duplicates clustering and grouping step is performed sample organism-wise. In the presence of a reference organism, synteny analysis is also carried out for small scale gene duplicate groups of tandem and proximal duplicates, to add a confidence layer of genomic positional context. Also, if a reference is involved, detailed gene duplicate group nature details like gene group expansion, de novo duplication etc., are provided for small scale gene duplicate groups of tandem and proximal duplicates in the output. (**Steps 7,8**)
</div>

<div align="justify">
	
(7)Finally, as a clustering approach is used to obtain gene duplicate groups/ arrays, for all the duplicate groups in the output, there is an internal scoring scheme used to classify the group as low, moderate and high confidence group that can help in interpreting the results accordingly. Along with the different duplicate group files, the singletons in a sample organism are also provided organism-wise. It is important to note that, in the reference-free mode, ortholog detection, synteny analysis, gene group nature analysis steps are absent. (**Steps 7,8**)
</div>

<div align="justify">
	
(9) Ka/Ks computation can also be done for all gene duplicates on an individual gene basis using either Nei Gojobori or Modified Ynag Nielsen methods. Please note that assembly FASTA file(s) or coding sequence file(s) are required along with the annotation file for Ka/Ks analysis. (**Step 9**)
</div>

<div align="justify">
	
(10) If expression data is available for the sample organism(s), expression analysis of gene duplicates can also be peformed with the script. This step gives out comprehensive information about the correlation among genes in duplicate groups, generates pairwise gene expression plots in a matrix figure and also helps determine divergent expression among duplicates (**Step 10**)
</div>


## Gene duplicates types output by DupyliCate
<img width="3120" height="1851" alt="Gene_duplicates_types" src="https://github.com/user-attachments/assets/ca081760-2e30-4eea-bfb3-94877b443d6f" />

## Installation and dependencies

### (1) Manual installation

Clone this repository

**Mandatory dependencies:**

- **Tools**: BLAST, DIAMOND, MMSeqs2

- **Python libraries**: pandas (v2.3.1 or greater), numpy (v2.3.2 or greater), seaborn (v0.13.2 or greater), matplotlib (v3.10.5 or greater), scipy (v1.16.1 or greater)

**Optional dependencies:**

- **Tools**: GeMoMa, BUSCO, AGAT, MAFFT, FastTree
  
- **Python libraries**: dendropy (v5.0.8), tqdm (v4.67.1)

### (2) Docker installation




## Running the script

```
Usage:

  python3 DupyliCate.py

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
## Description of output files

- **md5sum.tsv**: Contains the md5sums of every input file used in a particular run of the tool

- **run_parameters.json**: json file that keeps a record of all the parameters used for a particular run of the tool

- **Duplication landscape plots**: Histogram plots of normalized bit scores of a gene's second best hit and self hit in self alignment of the sample; Depicts the genome level gene duplication status.

- **BUSCO QC Summary file**: TSV file with information about BUSCO completeness, duplication, Feff (pseudo-ploidy number - approximate indication of ploidy)

- **Tandem duplicates, Proximal duplicates, Dispersed duplicates**: Folders containing organism-wise tandem, proximal and dispersed gene duplicates output TSV files

- **Mixed duplicates**: Mixed duplicates folder containing organism-wise mixed duplicates output TSV files (output in strict mode)

- **Duplicates_relationships**: Duplicates relationships folder containing organism-wise duplicates relationships across duplicate gruops (output in overlap mode)

- **Singletons**: Singletons folder containing organism-wise singleton gene output TSV files

- **Comparative_duplicates_table**: Comparative genomics table depicting the specific genes and their copy numbers in sample organisms; These specific genes are orthologs of user specified genes in the reference organism whose copy number 										variation the user wants to know

- **Ka_Ks_analysis**: Folder containing organism-wise Ka/Ks computation TSV files of gene duplicates on an individual gene basis along with statistical significance and nature of selection pressure

- **Duplicates_analysis**: Folder containing organism wise gene expression output folders; Each organism folder shows Plots folder, Stats folder, the kernel density estimation plot of correlation coefficients used for determining the 						           divergent expression threshold, and TXT file of perceived pseudogenes (genes that have low expression across the different RNASeq samples used for generating the counts table file) in that organism

- **Specific_duplicates_analysis**: Folder similar to the Duplicates_analysis folder, except that contains the respective output folders, and files for specific gene duplicates

