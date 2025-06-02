# DupyliCate
Light weight python tool for mining, classification and analysis of gene duplications.


![twinfinder drawio](https://github.com/user-attachments/assets/57f00b9a-e2e9-4d48-9266-d453e4342818)


```
Usage:

  python3 Dupylicate.py

Mandatory:

  --gff <full path to folder containing GFF files>
  --fasta <full path to folder containing WGS FASTA files> | --cds <full path to folder containing CDS FASTA files> | --pep <full path to folder containing PEP files>
  --ref <name of the organism to be taken as reference>

Optional:

  --score <auto for automatic threshold finding and <float number between 0 and 1> for manual threshold finding for segregating singletons and duplicates> DEFAULT is auto
  --mode <overlap | strict modes; in the overlap mode genes repeat among the different duplicates classification; 
	 in the strict mode there is no gene repetition and you have a new classification group called mixed 
	 duplicates containing the related connected components from the other three duplicate classes> DEFAULT is overlap
  --to_annotate <full path to file containing names of queries to be annotated separated by commas>
  --seq_aligner <choose one among blast | diamond | mmseqs2 > DEFAULT is blast
  --blast <full path to BLAST if not already in yopur PATH environment variable>
  --diamond <full path to diamond if not already in yopur PATH environment variable>
  --mmseqs <full path to mmseqs2 if not already in yopur PATH environment variable>
  --mafft <full path to MAFFT if not already in yopur PATH environment variable>
  --evalue <evalue for alignment> DEFAULT is 1e-5
  --gemoma <full path to GeMoMa if not already in yopur PATH environment variable>
  --cores <number of cores needed to run Dupylicate analysis> DEFAULT is 4
  --proximity <integer value for the number of intervening genes to detect proximal duplications> DEFAULT is 10
  --synteny_score <float value which is used as a cut-off or threshold for synteny analysis> DEFAULT is 0.5
  --flank <integer value specifying the number of flanking genes to be considered to determine the synteny window size in synteny analysis>
  --side <integer value for synteny support from either side of a flanking region of a synteny window>
  --ka_ks <yes | no to calculate ka, ks values> DEFAULT is 'no'
  --duplicates_analysis <yes | no for further statistical analysis of identified gene duplicates> DEFAULT is no
  --specific_duplicates_analysis <yes | no for further statistical analysis of specified ref genes' gene duplicates> DEFAULT is no
  --dpi <resolution value desired for plots low | moderate | high | very high > DEFAULT is moderate
  --analyse_disperse_duplicates <yes | no for statistical analysis of dispersed gene duplicates> DEFAULT is no
  --exp <full path to the folder with expression counts files>
  --genes <integer value referring to the number of genes to be taken for random pairing to determine threshold 
	  value for assessing functional divergence of gene duplicates based on gene expression values> DEFAULT is 10000
  --specific_genes <user given list of specific genes from the reference with one transcript name per line, for analysing gene duplication outputs>
  --clean up <yes | no cleans up intermediate files/ folders> DEFAULT is YES

Allowed file extensions:

			<'.gff', '.gff.gz', '.gff3', '.gff3.gz', '.fa', '.fa.gz', '.fna', '.fna.gz','.fasta', '.fasta.gz',
			'.genome.fasta', '.genome.fa', '.genome.fasta.gz', '.genome.fa.gz', '.genome.fna',
			'.genome.fna.gz','.cds.fa', '.cds.fasta','.cds.fa.gz', '.cds.fasta.gz', '.cds.fna',
			'.cds.fna.gz','.pep.fa','.pep.fa.gz','.pep.fasta','.pep.fasta.gz', '.pep.fna',
			'.pep.fna.gz','.tsv','.txt','.txt.gz','.tpms.txt','.tpms.txt.gz'>

```
