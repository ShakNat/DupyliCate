**Steps to test DupyliCate**

git clone https://github.com/ShakNat/DupyliCate

cd test_data

There are different possible test cases with the test files

```
1) Mining gene duplicates in multiple organisms with assembly FASTA/ PEP FASTA/ CDS FASTA input**

python3 dupylicate.py [ --fasta <FASTA_DIR> | --cds <CDS_DIR> | --pep <PEP_DIR>]
                      [--gff <GFF_DIR> ]
                      [--out <OUTPUT_DIR>]


2) Mining gene duplicates in sample organism(s) with assembly FASTA/ PEP FASTA/ CDS FASTA input and comparative genomics with a reference organism

python3 dupylicate.py [ --fasta <FASTA_DIR> | --cds <CDS_DIR> | --pep <PEP_DIR>]
                      [--gff <GFF_DIR> ]
                      [--ref <REF_NAME>]
                      [--out <OUTPUT_DIR>]

3) Mining gene duplicates in sample organism(s) with assembly FASTA/ PEP FASTA/ CDS FASTA input and expression analysis of gene duplicates

python3 dupylicate.py [ --fasta <FASTA_FILE> | --cds <CDS_FILE> | --pep <PEP_FILE>]
                      [--gff <GFF_FILE>]
                      [--exp <EXP_FILE>]
                      [--duplicates_analysis <yes>]
                      [--out <OUTPUT_DIR>]
```
**Further details**

- The test_data folder has subfolders for assembly FASTA, PEP FASTA and CDS FASTA inputs of A. thaliana Col-0 and Nd-1 accessions. The contigs and the genes correspond to the first 0.1 Mb of chromosome 1 in both the accessions

- Based on the test file runs, with the default fall-back threshold settings, a pair of tandem duplicate genes were identified in both of the above datasets

- Since, it is a very small dataset, using BUSCO-based thresholding (detailed in the main README) will fail and fall-back to default threshold settings; If needed, you could tweak this by changing the self normalized bit score and self similarity cutoff values

- When testing with the input directories, add a / at the end of the directories' paths

- The test file runs will work with the default GFF parameters config settings in the script

- The test expression dataset is available only for the Col-0 accession; Hence when testing for gene expression analysis of gene duplicates, specify the Col0 file paths in the respective flags for DupyliCate
