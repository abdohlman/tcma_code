# Code for generating the TCMA database

### Description

The pipeline used to generate the TCGA database is broken into five steps:

  1. Metagenomic screening (`pathseq.sh`)
  2. Tabulation of results (`collate_pathseq_output.py`)
  3. Analysis of taxa prevalence (`compute_prevalence.py`)
  4. Removal of putative contaminants (`decontaminate.py`)
  5. Averaging and normalization (`collapse_levels.py`)

### Python requirements

* pandas
* numpy
* skbio
* scipy
* statsmodels


## Scripts

### pathseq.sh

This script downloads TCGA sequencing data and runs PathSeq to screen them for microbial sequences. The BAM files to be analyzed are given by a [manifest](https://docs.gdc.cancer.gov/Encyclopedia/pages/Manifest_File/), which contains the list of files corresponding to a given TCGA project (eg. COAD) and experimental assay (eg. WGS) and must be of the format "`gdc_manifest.$assay.$project.txt`". Results for each file are saved to `./results/$assay/$project`.

The script relies on the [Slurm workload manager](https://slurm.schedmd.com/documentation.html) to executes an `sbatch` array, starting a job for each BAM file in the manifest. Additionally, because raw sequencing data is considered "controlled access" a GDC user token is neccessary. One can apply for access to raw TCGA sequencing data from [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v1.p1). The tools [GATK](https://gatk.broadinstitute.org/hc/en-us) 4.0< and samtools are also required.

Usage: 
~~~
sbatch --array 0-$(expr $(wc -l < $manifest) - 2)%25 pathseq.sh $token_file $manifest $min_clipped_read_length
~~~

**Parameters**
  * `$manifest` - GDC manifest file which must be of the format "`gdc_manifest.$assay.$project.txt`"
  * `$token_file` - text file containing GDC user token for controlled access data. This expires one month after download.
  * `$min_clipped_read_length` - this corresponds to a PathSeq parameter which defines the minimum trimmed read length required for alignment.


### collate_pathseq_output.py

This script parses the PathSeq output files and collates unambiguously-aligned reads for each taxa into an *n* x *m* table, with *n* taxa (NCBI taxonomy ID) and *m* sequencing runs ([UUIDs](https://docs.gdc.cancer.gov/Encyclopedia/pages/UUID/)). This table is saved to `./results/`.

Usage:
~~~
python collate_pathseq_output.py $project $assay $min_clipped_read_length
~~~

**Parameters**
  * `$project` - TCGA sequencing project (eg. COAD)
  * `$assay` - TCGA sequencing strategy (eg. WGS)
  * `$min_clipped_read_length` - this corresponds to a PathSeq parameter which defines the minimum trimmed read length required for alignment.


### compute_prevalence.py

This script compares the prevalence of taxa across various sample types. Metagenomic results from the TCGA brain cancer projects (GBM, LGG) must be present. For a given TCGA sequencing `$project`, it makes three comparisons:

  1. Tissue from `$project` vs. blood from `$project`
  2. Tissue from `$project` vs. tissue from brain cancer projects
  3. Blood from  `$project` vs. blood from brain cancer projects

Results for each comparison are saved to `./prevalence`.

Usage: 
~~~
python compute_prevalence.py $project $assay
~~~

**Parameters**
  * `$project` - TCGA sequencing project (eg. COAD)
  * `$assay` - TCGA sequencing strategy (eg. WGS)
  

### decontaminate.py


Usage: 
~~~
python decontaminate.py $project $assay
~~~

**Parameters**
  * `$project` - TCGA sequencing project (eg. COAD)
  * `$assay` - TCGA sequencing strategy (eg. WGS)


### collapse_levels.py

Usage: 
~~~
python collapse_levels.py $project $assay $statistic
~~~
