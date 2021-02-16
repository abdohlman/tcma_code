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

### 1. pathseq.sh
~~~
sbatch --array 0-$(expr $(wc -l < $manifest) - 2)%25 pathseq.sh $token_file $manifest $min_clipped_read_length
~~~

This script downloads TCGA sequencing data and runs PathSeq to screen them for microbial sequences. The BAM files to be analyzed are given by a [manifest](https://docs.gdc.cancer.gov/Encyclopedia/pages/Manifest_File/), which contains the list of files corresponding to a given TCGA project (eg. COAD) and experimental assay (eg. WGS) and must be of the format "`gdc_manifest.$assay.$project.txt`". Results for each file are saved to `./results/$assay/$project`.

The script relies on the [Slurm workload manager](https://slurm.schedmd.com/documentation.html) to executes an `sbatch` array, starting a job for each BAM file in the manifest. Additionally, because raw sequencing data is considered "controlled access" a GDC user token is neccessary. One can apply for access to raw TCGA sequencing data from [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v1.p1). The tools [GATK](https://gatk.broadinstitute.org/hc/en-us) 4.0< and samtools are also required.

*Parameters*
  * `$manifest` - GDC manifest file which must be of the format "`gdc_manifest.$assay.$project.txt`"
  * `$token_file` - text file containing GDC user token for controlled access data. This expires one month after download.
  * `$min_clipped_read_length` - this corresponds to a PathSeq parameter which defines the minimum trimmed read length required for alignment.


### 2. collate_pathseq_output.py
~~~
python collate_pathseq_output.py $project $assay $min_clipped_read_length
~~~
This script parses the PathSeq output files and collates unambiguously-aligned reads for each taxa into an *n* x *m* table, with *n* taxa (NCBI taxonomy ID) and *m* sequencing runs ([UUIDs](https://docs.gdc.cancer.gov/Encyclopedia/pages/UUID/)). This table is saved to `./results/`.

*Parameters*
  * `$project` - TCGA sequencing project (eg. COAD)
  * `$assay` - TCGA sequencing strategy (eg. WGS)
  * `$min_clipped_read_length` - this corresponds to a PathSeq parameter which defines the minimum trimmed read length required for alignment.


### 3. compute_prevalence.py
~~~
python compute_prevalence.py $project $assay
~~~
This script compares the prevalence of taxa across various sample types. Metagenomic results from the TCGA brain cancer projects (GBM, LGG) must be present. For a given TCGA sequencing `$project`, it makes three comparisons:

  1. Tissue from `$project` vs. blood from `$project`
  2. Tissue from `$project` vs. tissue from brain cancer projects
  3. Blood from  `$project` vs. blood from brain cancer projects

For each comparison, a FDR-normalized fisher exact test is used to estimate the likelihood a taxon is equiprevalent (found at similar rates in different sample types). Results are saved to `./prevalence`.

*Parameters*
  * `$project` - TCGA sequencing project (eg. COAD)
  * `$assay` - TCGA sequencing strategy (eg. WGS)
  

### 4. decontaminate.py
~~~
python decontaminate.py $project $assay
~~~
This script uses the prevalence results from `compute_prevalence.py` to decompose an abundance matrix into "tissue-resident" and "contaminant" fractions. The decomposition of observed metagenomic data (*K*) into tissue-resident (*T*) and contaminant (*C*) components for a given taxon in a given sample can be described using a mixture with two components of the form *K* = *mT* + *nC*, where *m* and *n* represent the estimated fractions of tissue-resident or contaminant sequencing reads belonging to a given taxon, respectively (such that *m* + *n* = 1). Species are binarily classified as tissue-resident or contamination according to their relative prevalence in tissue and blood (by default, classification of species is done using WGS prevalence, regardless of of the value of `$assay`). For each sequencing run, the proportions *m* and *n* are then assigned using the relative fractions of unambiguously aligned sequencing reads from species classified as tissue-resident or contaminant within the corresponding clade. Proportions for taxa with insufficient reads assigned to them are imputed using plate-wise, center-wise median, and dataset-wise medians, depending on how many samples are available.

The tissue-resident and contaminant components are saved to `./results/$project/$assay` and are denoted with file extensions `decontam.txt` and `contam.txt`, respectively. Imputed and pre-imputed tissue-resident and contamination proportions are saved to `./mixtures`, and denoted with denoted file extensions `imputed.txt` and `raw.txt`, respectively.

*Parameters*
  * `$project` - TCGA sequencing project (eg. COAD)
  * `$assay` - TCGA sequencing strategy (eg. WGS)


### 5. collapse_levels.py
~~~
python collapse_levels.py $project $assay $statistic
~~~
This script performs various normalization and averaging functions, preparing the TCMA database for downstream analysis using TCGA endpoints.  Microbial compositions are collected and decontaminated by sequencing run. Since samples may be sequenced multiple times and multiple samples may be donated from a given patient, microbial compositions must be grouped and averaged (1) sample-wise and (2) case-wise. This allows the comparison of TCMA to clinically relevant metadata collected at the sample-level (eg. anatomic site) or patient-level (eg. clinical stage) and well as molecular expression assays (eg. RNA-seq, RPPA). 

First, counts are normalized using reads-per-million (RPM) by default. Subsequently, matched samples and cases are collapsed using the mean relative abundance (for case-level averaging, tissue and blood samples are averaged spearately). Finally, a centered-log-ratio (CLR) transform is used to remove the positivity constraint and map the microbial compositions to a normalize distribution. Output from each step is saved to  `./results/$project/$assay`, formatted with the file extension `${sample_type}.${level}.${normalization}.txt` (eg. ...`.tissue.case.clr.txt`).

*Parameters*
  * `$project` - TCGA sequencing project (eg. COAD)
  * `$assay` - TCGA sequencing strategy (eg. WGS)
