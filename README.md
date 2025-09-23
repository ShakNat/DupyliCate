<img width="100" height="100" alt="draft_dupylicate_logo" src="https://github.com/user-attachments/assets/b46fb8c0-c6b7-440f-9323-5a8ad76c79c0" />

# DupyliCate 

<img width="3933" height="2583" alt="Dupylicate_github_schematic" src="https://github.com/user-attachments/assets/cbae189d-19d2-45a3-b75f-fb0dbc5c4a25" />


<div align="justify">
	
DupyliCate is a python tool for mining, classification and analysis of gene duplications. It can help find and classify gene duplicates in a number of organisms concurrently and is designed to have high throughput. It can be used in a reference-free manner for intra-species gene duplicates identification, and classification or with a reference organism for comparative genomic analysis. The gene duplicates will however be identified and classified using intra-species local alignment in both the cases. The only difference is that, in the presence of a reference, the orthologs for the sample organism genes in the reference organism genome will be assigned, and further analysis pertaining synteny and gene copy number variation will be carried out. The tool moreover can be used with a reference organism to obtain a comparative gene table of specific genes in the reference, whose copy number variation that the user wants to know in other sample organisms, thus facilitating in-depth comparative genomic analyses.


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
	
(7) Finally, as a clustering approach is used to obtain gene duplicate groups/ arrays, for all the duplicate groups in the output, there is an internal scoring scheme used to classify the group as low, moderate and high confidence group that can help in interpreting the results accordingly. Along with the different duplicate group files, the singletons in a sample organism are also provided organism-wise. It is important to note that, in the reference-free mode, ortholog detection, synteny analysis, gene group nature analysis steps are absent. (**Steps 7,8**)
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

- **Tools**: BLAST, DIAMOND, MMSeqs2 (latest versions preferred)

- **Python libraries**: pandas (v2.3.1 or greater), numpy (v2.3.2 or greater), seaborn (v0.13.2 or greater), matplotlib (v3.10.5 or greater), scipy (v1.16.1 or greater)

**Optional dependencies:**

- **Tools**: GeMoMa, BUSCO, AGAT, MAFFT, FastTree (latest versions preferred)
  
- **Python libraries**: dendropy (v5.0.8), tqdm (v4.67.1)

### (2) Docker installation




**Note:** If you are using a docker image of the tool, you need not specify the dependencies' full paths while running the tool; 
          The dependencies are built-in in the docker image and that makes it simple to use the tool across systems without the 
		  need for local installations.


## Running the script

```
Usage:

  python3 dupylicate.py --gff <GFF_FILE_OR_DIR>
                        [--fasta <ASSEMBLY_FASTA_FILE_OR_DIR> | --cds <CDS_FASTA_FILE_OR_DIR> | --pep <PEP_FASTA_FILE_OR_DIR>]
                        --out <DIR>


MANDATORY:

	--gff 							STR		full path to folder containing GFF files

	Choose one of FASTA/PEP/CDS

	[--fasta						STR		full path to folder containing WGS FASTA files

    --cds							STR		full path to folder containing CDS FASTA files

    --pep							STR		full path to folder containing PEP files]

    Output directory

	--out							STR		full path to output folder

OPTIONAL:

	--ref							STR			<Name of the organism to be taken as reference <Default is NA and the script runs in reference-free mode>

	--prokaryote					STR			<Use the flag if the organisms for analysis are prokaryotes>

	--pseudos_gff					STR			yes | no	<Optional inclusion or exclusion of pseudogenes with coding features in the GFF file> [DEFAULT is no]

	--gff_config					STR			<Full path to TXT file containing the gff parameters> 

	--mode							STR			overlap | strict [DEFAULT is overlap] 

	--to_annotate					STR			<Full path to TXT file> <containing names of queries and the corresponding reference organism for GeMoMa annotation separated by comma - Query,Reference -> one pair per line

	--seq_aligner					STR			blast | diamond | mmseqs2 [DEFAULT is diamond]

	--blast							STR			<Full path to BLAST if not already in yopur PATH environment variable>

	--diamond						STR			<Full path to diamond if not already in your PATH environment variable>

	--mmseqs						STR			<Full path to mmseqs2 if not already in your PATH environment variable>

	--mafft							STR			<Full path to MAFFT if not already in your PATH environment variable>

	--evalue						FLOAT		<evalue for alignment> [DEFAULT is 1e-5]

	--gemoma						STR			<Full path to GeMoMa if not already in your PATH environment variable>

	--qc							STR			yes | no	<Quality check with BUSCO> [DEFAULT is no]

	--busco							STR			Full path to BUSCO | busco_docker [DEFAULT is busco]

	--busco_version					STR			<BUSCO version in the format vx.x.x> [DEFAULT is v5.8.2]  <needed only if you have docker-based BUSCO installation> 

	--container_version				STR			<Docker container version of BUSCO> [DEFAULT is cv1]

	--docker_host_path				STR			<Full host folder path <needed for docker-based BUSCO installation>

	--docker_container_path			STR			<Full mount path in the docker container> <needed for docker-based BUSCO installation>

	--score							STR | FLOAT	auto | float number between 0 and 1	<DEFAULT is auto <for automatic threshold finding and for manual threshold finding for segregating singletons and duplicates> 

	--self_simcut					FLOAT		<Similarity percentage to remove self alignment hits with low similarity percentage> [DEFAULT is 50.0]

    --hits							INT			<Number of top hits to be considered for finding suitable orthologs in the reference organism> [DEFAULT is 10]

	--ortho_candidates				INT			<User specified integer value for listing the potential ortholog candidates> [DEFAULT is 3]

	--occupancy						FLOAT		<Occupancy cutoff for MAFFT aligned file trimming> [DEFAULT is 0.1]

    --scoreratio					FLOAT		<Ratio of forward alignment bit score and self alignment bit score of query to assess if the forward hit is valid> [DEFAULT is 0.3]

	--fwd_simcut					FLOAT		<Similarity percentage to remove forward alignment hits with low similarity percentage> [DEFAULT is 40.0]

    --cores							INT			<Number of cores needed to run Dupylicate analysis> [DEFAULT is 4]

	--proximity						INT			<Value for the number of intervening genes to detect proximal duplications> [DEFAULT is 10]

    --synteny_score					FLOAT		<Value which is used as a cut-off or threshold for synteny analysis> [DEFAULT is 0.5]

	--flank							INT			<Value specifying the number of flanking genes to be considered to determine the synteny window size in synteny analysis> [DEFAULT is 5]

	--side							INT			<Value for synteny support from either side of a flanking region of a synteny window> [DEFAULT is 1]

	--ka_ks							STR			yes | no	<To calculate ka, ks values> [DEFAULT is no]

	--ka_ks_method					STR			MYN | NG	<Methods for Ka/Ks ratio caclulation> [Default is NG]

	--duplicates_analysis			STR			yes | no	<For further statistical analysis of identified gene duplicates> [DEFAULT is no]

	--specific_duplicates_analysis	STR			yes | no	<For further statistical analysis of specified ref genes' gene duplicates> [DEFAULT is no]

	--dpi							STR			low | moderate | high | very high	<Resolution level desired for plots> [DEFAULT is moderate]

	--analyse_disperse_duplicates	STR			yes | no	<For statistical analysis of dispersed gene duplicates> [DEFAULT is no]

	--exp							STR			<Full path to the folder with expression counts files>

	--avg_exp_cut					FLOAT		<Value cutoff for average expression value across samples> <for classifying a gene as pseudogene based on gene expression ; Default is 1

	--genes							INT			<Value referring to the number of genes to be taken for random pairing> [DEFAULT is 10000] 

	--specific_genes				STR			<Full path to TXT file> <user given list of specific genes from the reference with one transcript name per line, for analysing gene duplication outputs>

	--clean up						STR			yes | no	<Cleans up intermediate files/ folders> [DEFAULT is YES]

ALLOWED FILE EXTENSIONS:

	<'.gff', '.gff.gz', '.gff3', '.gff3.gz', '.fa', '.fa.gz', '.fna', '.fna.gz','.fasta', '.fasta.gz',
	'.genome.fasta', '.genome.fa', '.genome.fasta.gz', '.genome.fa.gz', '.genome.fna',
	'.genome.fna.gz','.cds.fa', '.cds.fasta','.cds.fa.gz', '.cds.fasta.gz', '.cds.fna',
	'.cds.fna.gz','.pep.fa','.pep.fa.gz','.pep.fasta','.pep.fasta.gz', '.pep.fna',
	'.pep.fna.gz','.tsv','.txt','.txt.gz','.tpms.txt','.tpms.txt.gz'>

```

### GFF fields config file preparation instructions for Dupylicate.py:
	
 1) This is a simple .txt file containing the gff parameters in different columns

 2) The columns can be separated by tabs or spaces

 3) There are four columns mandatorily needed in the config file in the following order:

	(i) base file name - same as the base name you use for the gff file | all in case all the gff have the same gff pattern
		
  	(ii) child_attribute: attribute field of the mRNA or transcript feature in the file like ID

    <div align="justify">
	(iii) child_parent_linker: attribute field of the mRNA or transcript, CDS, exon features that link them with their 
			  respective parent feature like Parent - Note: base assumption by the tool is that all child levels
			  have the same child-parent linker attribute fields. For eg., if Parent is the child-parent linker in the mRNA feature line,
			  then Parent will be the child-parent linker for all other child-level feature lines in the GFF
    </div>
							
	(iv) parent_attribute: attribute field of the gene feature like ID

    (v) By default in the script, child_attribute is ID, child_parent_linker is Parent, and parent_attribute is ID

**Sample config file and GFF file example:**
```
If the config looks like this -

all	Name	Parent	ID

And the corresponding GFF file looks as below - 

##gff-version 3

##annot-version Araport11

##species Arabidopsis thaliana columbia

Chr1    	phytozomev12    	gene    	3631    	5899    	.       	+       	.       	ID=AT1G01010.Araport11.447;Name=AT1G01010

Chr1    	phytozomev12    	mRNA    	3631    	5899    	.       	+       	.       	ID=AT1G01010.1.Araport11.447;Name=AT1G01010.1;pacid=37401853;longest=1;geneName=NAC001;Parent=AT1G01010.Araport11.447
```

Understanding the GFF config file:

- The first column of the file says all. This means all the files in the analysis will have the same GFF file format/ pattern

- The second column that is the child attribute column is Name. So the text following the Name field in the last column of
  the mRNA feature will be picked which is AT1G01010.1 above

- The third column that is the child-parent linker column is Parent. It is the field Parent in mRNA that links it to its parent gene,
  and in the above example its corresponding text picked will be AT1G01010.Araport11.447

- The fourth column is the parent attribute that is mentioned as ID. In the above example it is the ID field in the gene feature line which is AT1G01010.Araport11.447

<div align="justify">
	
- **IMPORTANT:** It is important to note that the child-parent linker and the parent must be chosen in such a way that they point to the same text. For example,
  both the child-parent linker and the parent attribute in the above example point to AT1G01010.Araport11.447; This is important to ensure that the transcripts
  correctly map to the parent gene especially in the alternate transcript removal step
</div>

### Preparing list of reference genes for specific analysis

1) If you have a known list of genes in the reference organism whose copy number variation you want to analyze in the sample,
   the --specific_genes flag can be used. 

2) This needs the full path to a simple .txt file

3) This .txt file should have one transcript name per line

### Preparing the txt file for GeMoMa annotation

1) In case, some of your input files lack structural annotation, the --to_annotate flag in the script can be used

2) This needs the full path to a simple .txt file

3) The name of the organism to be annotated followed by comma and the name of the reference organism to
   be used for annotation must be mentioned in a line

4) If there are multiple organisms for annotation, specify each of them along with their respective reference
   in a single line in the format as specified above

   eg.
   
    Vamurensis,Vvinifera
   
   	Vrotundifolia,Vvinifera

## Helper scripts

```
Usage:

	python3 AGAT_wrapper.py

MANDATORY:

	--gff_in <full path to GFF3 file or folder containing GFF3 annotation files>

	--gff_out <full path to output folder or output file>

OPTIONAL:

	--agat <full path to the agat_convert_sp_gxf2gxf.pl script including the script name>

```

```
Usage:

	python3 Fasta_fix.py

MANDATORY:

	--in <full path to a folder containing FASTA files or a single FASTA file>

	--out <full path to output directory>

	--config <full path to config file including the config file name>

OPTIONAL:

	--cores <number of cores for running the script>
```

### Config file preparation instructions for Fasta_fix.py:

1. The config file is a simple .txt file

2. Specify the organism name i.e. the base name without the extension of your file in the first column

	<div align="justify">
3. If the same set of string manipulation operations are to be performed for all the files in a folder,
	just specify the word all in the first column and follow it up in the next columns with the desired
	operations - This single line is enough if the same set of operations are to be performed on all the files
	</div>

	<div align="justify">
4. Specify the various string manipulation operations to be performed on the FASTA header of this particular
	organism's file in the subsequent columns separated by white space or tab
	</div>
	
5. **IMPORTANT**: The operations will be performed in the same order as you specify them in the config wise i.e.
    column wise order of the different operations you specify will be followed by the script; 
    Hence operation ORDER is IMPORTANT

6. Possible operations and the manner in which they are to be specified are as follows
	<div align="justify">
   i. extract: <to extract text of a particular attribue in the header>
      example- extract:ID=:; <This tells extract the text following the attribute ID= in the header
								until the delimiter ; is encountered><If you specify no delimiter, a
								default list of delimiters will be searched automatically by the script>
	</div>

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
  
    &nbsp;&nbsp;&nbsp;&nbsp;lowercase <to write the full text in lower case>

## Description of output files

- **md5sum.tsv**: Contains the md5sums of every input file used in a particular run of the tool

- **run_parameters.json**: json file that keeps a record of all the parameters used for a particular run of the tool

- **Duplication landscape plots**: Histogram plots of normalized bit scores of a gene's second best hit and self hit in self alignment of the sample; Depicts the genome level gene duplication status.

- **BUSCO QC Summary file**: TSV file with information about BUSCO completeness, duplication, Feff (pseudo-ploidy number - approximate indication of ploidy)

- **Tandem duplicates, Proximal duplicates, Dispersed duplicates**: Folders containing organism-wise tandem, proximal and dispersed gene duplicates output TSV files

- **Mixed duplicates**: Mixed duplicates folder containing organism-wise mixed duplicates output TSV files (output in strict mode)

- **Duplicates_relationships**: Duplicates relationships folder containing organism-wise duplicates relationships across duplicate gruops (output in overlap mode)

- The different gene duplicates classification output files have a number of columns. In all the duplicates results files, both in the presence or absence of a reference organism, the following columns are present:
  
  (i)   Pairwise gene distance (bp) - Genetic distance in bp between two neighbouring genes in an identified duplicates group

  <div align="justify">
  (ii)  Actual intervening gene number: The tool considers only protein coding genes for identifying and classifying gene duplicates. But it also important to
        know the actual number of protein coding and non-coding genes in-between gene duplicates, and the actual intervening gene number gives this information
  </div>
  
  (iii) Apparent intervening gene number: The number of protein coding genes between neighbouring genes in an identified duplicates group
  
  (iv) Group confidence: The gene duplicates are clustered into groups. This column gives a confidence score for the reliability of the identified groups -
  
  &nbsp;confidence_score <= 0.3 -> Low confidence gene duplicates group
  
  &nbsp;confidence_score <= 0.5 -> Moderate confidence gene duplicates group
  
  &nbsp;confidence_score > 0.5 -> High confidence gene duplicates group
  
  (v) Nature of duplicate group: This column is present only in the small scale duplicate groups - tandems and proximals output files in the presence of a reference organism;
  
  &nbsp;Small scale duplicates are mainly responsible for the evolutionary innovations seen in stress response mechanisms and biosynthetic  pathways;
  
  &nbsp;Understanding their nature like whether they lead to gene expansion, conservation or de novo duplication can lead to crucial biological insights;
  
  &nbsp;The nature of duplicate group provides this information as detailed in the image below
  

<img width="2985" height="1743" alt="Gene_duplicates_group_nature" src="https://github.com/user-attachments/assets/8dbc31da-955f-40c1-ab9f-78c8ff902ae7" />


- **Singletons**: Singletons folder containing organism-wise singleton gene output TSV files

-  In the presence of a reference organism, the small scale gene duplicates (tandems, proximals) files and the singletons files have a column called Synteny information
  
    (vi) Synteny information: This column evaluates the synteny between the identified gene duplicates and their corresponding orthoolog genes in the reference organism

<div align="justify">
	
- **Copy_number_table.tsv**: This is found only in the presence of a reference organism and if the user has given specific reference genes or analysis using the --specific_genes flag;
  
  &nbsp;Comparative genomics table depicting the specific user-given reference genes and their copy numbers in sample organisms;
  
  &nbsp;These specific genes are orthologs of user specified genes in the reference organism whose copy number variation the user wants to know;
  
  &nbsp;Every cell in each organism-wise column has the orthologs corresponding to specific reference organism genes in that sample organism;
  
  &nbsp;Each such identified orthologous group has a score called Ortholog group confidence score (OGCS) appended to it to show the reliability of the identified orthologous group.
  
  &nbsp;The OGCS scoring system is as follows:
  
  &nbsp;OGCS <= 0.5 -> Low confidence orthologous group
  
  &nbsp;0.5 < OGCS <= 0.8 -> Moderate confidence orthologous group
  
  &nbsp;OGCS > 0.8 -> High confidence orthologous group

</div>
  
  The ortholog confidence scoring is more stringent than the gene duplicates confidence scoring system because of the species specific variations and distance considerations between the reference and sample organisms

<div align="justify">
	
- **Copy_number_table.xlsx**: This output file is also found only in the presence of a reference organism and if the user has given specific reference genes or analysis using the --specific_genes flag;
  
  &nbsp;It is the same file as Copy_number_table.tsv, except that, it is red colour coded in the cells that do not have orthologs for specific reference genes in the sample organism(s);
  
  &nbsp;Being an xlsx file, this colour difference can be used as a quick presence-absence variation assessment of the user-specified genes for comparative genomics between the reference and the sample organism(s)
</div>

- **Summary.tsv**: Summary file consolidating the number of different gene duplicates in each analysed sample organism; It also lists the number of gene duplicate groups classified as low/ moderate/ high confidence groups 

- **Ka_Ks_analysis**: Folder containing organism-wise Ka/Ks computation TSV files of gene duplicates on an individual gene basis along with statistical significance and nature of selection pressure

- **Duplicates_analysis**: Folder containing organism wise gene expression output folders; Each organism folder shows Plots folder, Stats folder, the kernel density estimation plot of correlation coefficients used for determining the
  divergent expression threshold, and TXT file of perceived pseudogenes (genes that have low expression across the different RNASeq samples used for generating the counts table file) in that organism;
  
  &nbsp;It is important to note that some gene duplicate groups can have statistics folder present, but might not have corresponding plots folder due to the possibility of all or most of them being perceived pseudogenes,
  disabling their correlation analysis and gene expression plotting

- **Specific_duplicates_analysis**: Folder similar to the Duplicates_analysis folder, except that contains the respective output folders, and files for specific gene duplicates

## Common error prone steps and some recommendations

- Please remember to add / at the end of the folder path when giving a folder of input files to DupyliCate

- Please ensure that the files have the allowed extensions as listed above in the main script usage

- The base names of the files excluding the extensions must be the same across all input files, the config files, and the annotation txt file for GeMoMa annotation

- Gene duplications classification needs the GFF file(s) as the inputs. Since GFF files have wide varying formats, the analysis can be interrupted at the validation step due to GFF file formatting issues
  
- Next, the FASTA headers need to match specific GFF attribute fields. Otherwise, the analysis would stop after the validation and PEP file generation step

- Please look into the GFF config file preparation instructions to avoid such errors during analysis

<div align="justify">
- The script is designed in such a way that it can start from the point, an analysis was stopped or interrupted. But if the output files of the previous steps are truncated or empty,
  this will not be captured and may cause errors downstream. It is safe to remove the output files of the step that was interrupted while retaining the files generated in the earlier steps to facilitate a seamless run despite interruptions.
</div>

## Third party tool references

- Altschul, S.F., Gish, W., Miller, W., Myers, E.W., Lipman, D.J. (1990) “Basic local alignment search tool.” J. Mol. Biol. 215:403-410. PubMed

- Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x

- Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi: 10.1038/nbt.3988 (2017).

- Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. doi:10.1371/journal.pone.0009490

- Kazutaka Katoh, Daron M. Standley, MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability, Molecular Biology and Evolution, Volume 30, Issue 4, April 2013, Pages 772–780, https://doi.org/10.1093/molbev/mst010

- J. Keilwagen, F. Hartung, M. Paulini, S. O. Twardziok, and J. Grau Combining RNA-seq data and homology-based gene prediction for plants, animals and fungi. BMC Bioinformatics, 2018. doi: 10.1186/s12859-018-2203-5

- Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  
(Version v1.5.1). Zenodo. https://www.doi.org/10.5281/zenodo.3552717

- Fredrik Tegenfeldt, Dmitry Kuznetsov, Mosè Manni, Matthew Berkeley, Evgeny M Zdobnov, Evgenia V Kriventseva, OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes, Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D516–D522, https://doi.org/10.1093/nar/gkae987

- Jeet Sukumaran, Mark T. Holder, DendroPy: a Python library for phylogenetic computing, Bioinformatics, Volume 26, Issue 12, June 2010, Pages 1569–1571, https://doi.org/10.1093/bioinformatics/btq228


## When using DupyliCate in your research please cite




  

