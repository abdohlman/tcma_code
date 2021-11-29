#!/usr/bin/env Rscript

# USAGE: Rscript convert_to_phyloseq.R -p $project  -a $assay -l sample -n rpm

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 11-29-21                      #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################

library("optparse")
suppressMessages(library(phyloseq))

args = commandArgs(trailingOnly=TRUE)

option_list = list(

  make_option(c("-p", "--project"), type="character", default=NULL,
              help="TCGA sequencing project, e.g. COAD)", metavar="character"),

  make_option(c("-a", "--assay"), type="character", default=NULL,
              help="TCGA experimental strategy, e.g. WGS)", metavar="character"),
  
  make_option(c("-l", "--level"), type="character", default=NULL,
              help="Collapsed level to use, e.g. sample", metavar="character"),

  make_option(c("-n", "--normalization"), type="character", default=NULL,
              help="Data normalization, e.g. rpm", metavar="character"),
  
  make_option(c("-t", "--sample-type"), type="character", default="tissue",
              help="Sample type to use. Default: %default]", metavar="character"),
  
  make_option(c("-s", "--statistic"), type="character", default="unambiguous.decontam",
              help="Statistic to acquire from pathseq output. Default: %default]", metavar="character"),
  
  make_option(c("-d", "--domain"), type="character", default="bacteria",
              help="Domain to apply pipeline to, Default: %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$project)) {
  print("Please provide a project")
  print_help(opt_parser)
} else if (is.null(opt$assay)) {
  print("Please provide an assay")
  print_help(opt_parser)
} else if (is.null(opt$level)) {
  print("Please provide a level")
  print_help(opt_parser)
} else if (is.null(opt$normalization)) {
  print("Please provide a normalization")
  print_help(opt_parser)
}

project = opt$project
assay = opt$assay
level = opt$level
norm = opt$normalization
sample_type = opt$`sample-type`
stat = opt$statistic
domain = opt$domain

results_path =  paste('./collapsed',project,assay,'',sep='/')

# NEWICK TREE
tree_file = paste('./taxonomy/newick/taxa',domain,'species.nw',sep='.')
print(tree_file)
phyloTree <- read_tree(tree_file)

# TAXONOMY TABLE
tax_table_file <- paste('./taxonomy/tables/tax_table',domain,'names.txt',sep='.')
tax_df <- read.delim(tax_table_file,row.names=1)
taxonomyTable <- tax_table(as.matrix(tax_df))

# METADATA
meta_file <- paste('./metadata/metadata',project,level,'txt',sep='.')
meta_df <- read.delim(meta_file,row.names=1)
if (level == 'sample') { row.names(meta_df) = meta_df$bcr_sample_barcode }
if (level == 'case') { row.names(meta_df) = meta_df$case.bcr_patient_barcode }
sampleData <- sample_data(meta_df)

# ABUNDANCE DATA
abundance_file = paste0(results_path, paste(domain,stat,sample_type,level,norm,'txt',sep='.') )
data_df <- data.frame(t(read.delim(abundance_file,row.names=1, check.names = FALSE)), check.names = FALSE)
if (norm != 'clr') { data_df[is.na(data_df)] <- 0 }
data_df <- data.frame(t(data_df), check.names = FALSE)
otuData <- otu_table(data_df, taxa_are_rows = TRUE)
print(dim(otuData))

# PUT IT ALL TOGETHER
physeq_file <- paste0(results_path, paste('physeq',domain,stat,sample_type,level,norm,'Rds',sep='.'))
print(paste("Saving to",physeq_file))
physeq <- phyloseq(phyloTree, otuData, sampleData, taxonomyTable)
saveRDS(physeq, file = physeq_file)

print("Done.")

