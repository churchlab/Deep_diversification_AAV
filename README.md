# Bioinformatics pipeline overview 
This document summarizes the two main components of the bioinformatics analysis that was used to generate and parsed data for the paper [Deep diversification of an AAV capsid protein by machine learning](https://www.nature.com/articles/s41587-020-00793-4#Abs1). For machine learning models [see this](https://github.com/google-research/google-research/tree/master/aav). 


A processed version of the data is available in the data folder (this should look similar to what the processing pipeline outputs). For additional annotation (e.g. model scores), and training data browse through [these datasets](https://github.com/alibashir/aav). For raw sequencing data [see NCBI](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA673640/). Additional meta-data and artifacts to reproduce the results can be found in this [Dropbox link](https://www.dropbox.com/sh/lmd8wmgibz24r2h/AABnfUHq6k-H_sT7I0UGR_D5a?dl=0) (too big to host on github, NCBI did not support these directory structure). 

+   [Synthesis pipeline](#synthesis-pipeline)
	+   *Step 1*: Assembles the nucleotide sequence for the corresponding protein sequence variants such that it can be generated and processe with the desired cloning strategy. 
	+   *Step 2*: Tests the dataframe produced by Step 1 to ensure that the library has the intended composition, while additionally testing that the correct RE sites are in each sequence. This step produces the files that are sent to Agilent for synthesis. 
	+   *Step 3*: Simulates the cloning process in silico, to ensure that the library can be successfully produced with the set of primers, plasmid backbone, and other molecular parameters.  
+   [Parsing pipeline](#parsing-pipeline)
	+   *Step 1*: Merge fastq files using PEAR.
	+   *Step 2*: Count the number of variants across sequencing files.
	+   *Step 3*: Compute selection scores based on the raw count files. 



Details below. 
# Synthesis Pipeline

Takes the AA sequences designed by ML and produces nucleotide sequences to be printed for synthesis such that it is compatible with our cloning strategy. 

### Requirements
``` Python
Pandas
Numpy
BioPython
PyDNA
editdistance
```

## Description

### Step 1
Assembles the nucleotide sequence for the corresponding protein sequence variants such that it can be generated and processe with the desired cloning strategy. 

Note: We used barcodes in our original design but actually never used them as identifiers for variants.

#### Input files: 

Barcode designs:

+   `barcodes16-1.txt`  from John A. Hawkins et al. PNAS 2018 https://www.pnas.org/content/115/27/E6217 (not used for analysis)
or if barcodes already chosen:

+   `c1barcodes16-1_app_BsrBI.txt` these are a selected group of barcodes compatible with our cloning strategy.

Designed Variants:

+   `chip1_GAS_nredundant.csv` the ML designed variants

+   `backfill_random_doubles.csv` random doubles to backfill the chip if there is room

+   `singles.csv`  set of all single mutations to the WT 

Primer files:
+   `skpp15-forward.fasta` forward primers

+   `skpp15-reverse.fasta` reverse primers

#### Output files:

+   `chip_df.csv` contains the library sequences

+   [Optional] `c1barcodes16-1_app_BsrBI.txt` as selected barcodes



### Step 2
Tests the dataframe produced by Step 1 to ensure that the library has the intended composition, while additionally testing that the correct RE sites are in each sequence. This step produces the files that are sent to Agilent for synthesis. 

#### Input files: 
+   `chip_df.csv`  contains the library sequences

#### Output files:

+   `chip_for_agilent.txt`  this is what is sent to Agilent



### Step 3
Simulates the cloning process in silico, to ensure that the library can be successfully produced with the set of primers, plasmid backbone, and other molecular parameters.   

#### Input files: 

Primer files

+   `skpp15-forward.fasta`  forward primer

+   `skpp15-reverse.fasta`  reverse primer 


+   `chip_df.csv` contains the library sequences

# Parsing Pipeline

Takes the fastq nucleotide sequences from experimental sequencing runs and maps them back to original AA sequences and computes selection scores (We performed two sequencing runs, hence step 1 and 2 should be run on both sets before combining them on step 3)
### Requirements
```
PEAR
Pandas
Biopython
```

## Description


### Step 1
Merge fastq files using PEAR.

#### Input files: 
+   `fastq files in experimental run folder` contains all the fastq files 
+   `manifest file for samples`  contains the mapping between file names and the relevant samples

#### Output files:

+   `merged files in Parsed_data/merged`  merged fastq files




### Step 2
Count the number of variants across sequencing files.

#### Input files: 
+   `merged files in Parsed_data/merged`  merged fastq files
+   `designed_variants.csv`   set of designed AAs and corresponding coding nucleotides 

#### Output files:

+   ` files in Parsed_data/library`  merged fastq files
+ 	`raw_counts_raw_counts_NextSeq_run<run_num>.csv` raw counts 

### Step 3
Compute selection scores based on the count files. 

#### Input files: 
+ 	`raw_counts_raw_counts_NextSeq_run1.csv` raw counts from run1 sequencing
+ 	`raw_counts_raw_counts_NextSeq_run2.csv` raw counts from run2 sequencing (3x)
+   `chip_df.csv`   *[this is the output of the synthesis pipeline]* set of designed AAs and corresponding coding nucleotides 

#### Output files:

+   `library_w_selection_scores.csv`  computed selection scores for the libraries together.


