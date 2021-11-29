# Code for generating the TCMA database

### Description

The pipeline used to generate the TCGA database is broken into five steps:

  1. Metagenomic screening (`pathseq.sh`)
  2. Tabulation of results (`collate_pathseq_output.py`)
  3. Analysis of taxa prevalence (`compute_prevalence.py`)
  4. Removal of putative contaminants (`decontaminate.py`)
  5. Averaging and normalization (`collapse_levels.py`)
  6. Combine projects (`combine_projects.py`)
  7. package_phyloseq.R (`package_phyloseq.R`)

### Python requirements

* pandas
* numpy
* skbio
* scipy
* statsmodels

### R requirements

* optparse
* phyloseq

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
python collate_pathseq_output.py -p $project -a $assay -m $min_clipped_read_length
~~~
This script parses the PathSeq output files and collates unambiguously-aligned reads for each taxa into an *n* x *m* table, with *n* taxa (NCBI taxonomy ID) and *m* sequencing runs ([UUIDs](https://docs.gdc.cancer.gov/Encyclopedia/pages/UUID/)). This table is saved to `./results/$project/$assay/$domain.$statistic.txt`. By default, unambiguously aligned reads from pathseq are taken for bacteria only.

*Required Parameters*
  * `-p`, `--project` - TCGA sequencing project (eg. `COAD`)
  * `-a`, `--assay` - TCGA sequencing strategy (eg. `WGS`)
  * `--min-clipped-read-length` - this corresponds to a PathSeq parameter which defines the minimum trimmed read length required for alignment (eg. `50`).

*Optional Parameters*
  * `-s`, `--statistic` - Pathseq statistic to analyze (default is `unambiguous`)
  * `-d`, `--domain` - Microbial domain to analyze (default is `bacteria`)


### 3. compute_prevalence.py
~~~
python compute_prevalence.py -p $project -a $assay [ -d $domain -s $statistic ]
~~~
This script compares the prevalence of taxa across various sample types. Metagenomic results from the TCGA brain cancer projects (GBM, LGG) must be present. For a given TCGA sequencing `$project`, it makes three comparisons:

  1. Tissue from `$project` vs. blood from `$project`
  2. Tissue from `$project` vs. tissue from brain cancer projects
  3. Blood from  `$project` vs. blood from brain cancer projects

For each comparison, a FDR-normalized fisher exact test is used to estimate the likelihood a taxon is equiprevalent (found at similar rates in different sample types). Results are saved to `./prevalence`.

*Required Parameters*
  * `-p`, `--project` - TCGA sequencing project (eg. `COAD`)
  * `-a`, `--assay` - TCGA sequencing strategy (eg. `WGS`)

*Optional Parameters*
  * `-s`, `--statistic` - Pathseq statistic to analyze (default is `unambiguous`)
  * `-d`, `--domain` - Microbial domain to analyze (default is `bacteria`)
  * `--tissue-comparison` – Tissue type to use as comparision. (options are `brain` (default), `ovary`.
		Requires data for GBM/LGG and OV TCGA projects, respectively.


### 4. decontaminate.py
~~~
python decontaminate.py -p $project -a $assay
~~~
This script uses the prevalence results from `compute_prevalence.py` to decompose an abundance matrix into "tissue-resident" and "contaminant" fractions. The decomposition of observed metagenomic data (*K*) into tissue-resident (*T*) and contaminant (*C*) components for a given taxon in a given sample can be described using a mixture with two components of the form *K* = *mT* + *nC*, where *m* and *n* represent the estimated fractions of tissue-resident or contaminant sequencing reads belonging to a given taxon, respectively (such that *m* + *n* = 1). Species are classified as tissue-resident or contamination (ie. *m*, *n* in {0, 1}) according to their relative prevalence in tissue and blood (by default, classification of species is done using WGS prevalence data for the same project). For higher-level taxa, the proportions *m* and *n* are calculated recursively from the species- to phylum- level, using the relative fractions of reads classified as tissue-resident or contaminant from taxonomic level below (ie. *m*, *n* in [0, 1]). That is, genus-level mixtures are estimated from the proportion of tissue-resident vs. contaminant species-level reads, family-level mixtures are then estimated from the proportion of tissue-resident vs. contaminant genus-level reads, and so on. The tissue-resident and contaminant components are saved to `./results/$project/$assay` and are denoted with file extensions `decontam.txt` and `contam.txt`, respectively. Tissue-resident and contamination proportions are saved to `./mixtures`.

*Required Parameters*
  * `-p`, `--project` - TCGA sequencing project (eg. `COAD`)
  * `-a`, `--assay` - TCGA sequencing strategy (eg. `WGS`)

*Optional Parameters*
  * `-s`, `--statistic` - Pathseq statistic to analyze (default is `unambiguous`)
  * `-d`, `--domain` - Microbial domain to analyze (default is `bacteria`)
  * `--q-value` - Maximum FDR q-value at which to classify taxa as tissue-resident (default is `0.05`)
  * `--max-prevalence-blood` – Maximum blood prevalence at which to classify taxa as tissue-resident (default is `0.2` or 20%)
  * `--rank-for-classification` – Taxonomic rank to be used for classification. Abundance for higher-order taxa will be adjusted recursively. No features below this rank will be preserved (default is `species`)
  * `--stat-for-classification` – Pathseq output statistic to be be used for classification (default is `unambiguous`).
  * `--assay-for-classification` – Assay to be be used for classification (default is `WGS`).
  * `--custom-tissue-resident-list` – Use a custom list of tissue-resident taxa. Taxa not in this list will be classified as contamination. Provide a file path to a text file containing the list of taxa separated by newlines.


### 5. collapse_levels.py
~~~
python collapse_levels.py -p $project -a $assay -s $statistic [ -d $domain ]
~~~
This script performs various normalization and averaging functions, preparing the TCMA database for downstream analysis using TCGA endpoints.  Microbial compositions are collected and decontaminated by sequencing run. Since samples may be sequenced multiple times and multiple samples may be donated from a given patient, microbial compositions must be grouped and averaged (1) sample-wise and (2) case-wise. This allows the comparison of TCMA to clinically relevant metadata collected at the sample-level (eg. anatomic site) or patient-level (eg. clinical stage) and well as molecular expression assays (eg. RNA-seq, RPPA). 

First, counts are normalized using reads-per-million (RPM) by default. Subsequently, matched samples and cases are collapsed using the mean relative abundance (for case-level averaging, tissue and blood samples are averaged speparately). Finally, a centered-log-ratio (CLR) transform is used to remove the positivity constraint and map the microbial compositions to a normalize distribution. Output from each step is saved to  `./results/$project/$assay`, formatted with the file extension `${sample_type}.${level}.${normalization}.txt` (eg. ...`.tissue.case.clr.txt`).

*Required Parameters*
  * `-p`, `--project` - TCGA sequencing project (eg. `COAD`)
  * `-a`, `--assay` - TCGA sequencing strategy (eg. `WGS`)

*Optional Parameters*
  * `-s`, `--statistic` - Pathseq statistic to analyze (default is `unambiguous.decontam`)
  * `-d`, `--domain` - Microbial domain to analyze (default is `bacteria`)
  * `--no-rpm` – Do not calculate RPM prior to normalizations. 
  * `--min-clr-prevalence` – Remove OTUs with lesser than this fraction of prevalence prior to calculating CLR. The CLR transform performs best with lower sparsity (default is `0.25` or 25%).


### 6. combine_projects.py
~~~
python combine_projects.py -n $combined_project_name -p $project_list -a $assay
~~~

This script concatenates together the data, metadata (in `./metadata/`), and read statistics from projects in `$project_list` and creates a new project name `$combined_project_name` in `./results`. To create the TCMA database, we set `combined_project_name=TCMA` and `project_list=COAD,READ,HNSC,STAD,ESCA`.

*Required Parameters*
  * `-a`, `--assay` - TCGA sequencing strategy (eg. WGS)
  * `-n`, `--combined-project-name` - Name for newly combined project set (eg. `TCMA`)
  * `-l`, `--projects-list` - Comma-separated list of projects to combine (eg. `COAD,READ,HNSC,STAD,ESCA`)

*Optional Parameters*
  * `-s`, `--statistic` - Pathseq statistic to analyze (default is `unambiguous.decontam`)
  * `-d`, `--domain` - Microbial domain to analyze (default is `bacteria`)


### 6. package_phyloseq.R
~~~
Rscript convert_to_phyloseq.R -p $project  -a $assay -l $level -n $normalization [ -s $statistic -d $domain -t $sample_type ]
~~~

This script combines collapsed data (in `./collapsed`), metadata (in `./metadata/`), and taxonomy (in `./taxonomy`) to create a phyloseq `.Rds` object, saved to `./collapsed`. 

*Required Parameters*
  * `-p`, `--project` - TCGA sequencing project (eg. `COAD`)
  * `-a`, `--assay` - TCGA sequencing strategy (eg. `WGS`)
  * `-l`, `--level` – Collapsed level to use (eg. `file`, `sample`, or `case`)
  * `-n`, `--normalization` – Data normalization to use (eg. `rpm`, `reads`, `rpm.relabund`, `reads.relabund`, `rpm.clr`, `reads.clr`)

*Optional Parameters*
  * `-s`, `--statistic` - Pathseq statistic to analyze (default is `unambiguous.decontam`)
  * `-d`, `--domain` - Microbial domain to analyze (default is `bacteria`)
  * `-t`, `--sample-type` – Sample type to use (eg. `tissue` (default) or `blood`)





