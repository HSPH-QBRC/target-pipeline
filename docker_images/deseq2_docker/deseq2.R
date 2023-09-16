suppressMessages(suppressWarnings(
    library("DESeq2", character.only=T, warn.conflicts = F, quietly = T)))


###########################################################################

# args from command line:
args<-commandArgs(TRUE)

# the raw counts for the TCGA cohort, in TSV format
RAW_COUNT_MATRIX <- args[1]

# annotation file with an 'expression_state' column featuring
# values of 'low' or 'high'
ANNOTATIONS <- args[2]

# the output file for differential expression results
DGE_OUTPUT <- args[3]

###########################################################################

# Don't mangle the column names quite yet. Keep them as-is and then we can convert them after we
# create a column name mapping
counts <- read.table(RAW_COUNT_MATRIX, 
                    sep='\t', 
                    stringsAsFactors=F, 
                    row.names=1, 
                    header=T, 
                    check.names=T)

# read the annotation file which has the high/low designations
ann_df = read.table(ANNOTATIONS, sep='\t', row.names=1, header=TRUE)

# to align the annotation and counts, we need to have a common set of names.
# Hence, we have to modify the rownames of the annotations so they will match
rownames(ann_df) <- make.names(rownames(ann_df))

# we will need to have the "expression_state" as a factor
# for use with DESeq2
ann_df$expression_state <- factor(ann_df$expression_state, levels=c('low', 'high'))

# Subset the raw count matrix to only include the low
# and high-expressed samples:
counts = counts[,rownames(ann_df)]


dds <- DESeqDataSetFromMatrix(countData = counts,
							  colData = ann_df,
							  design = ~expression_state)
dds <- DESeq(dds)
res <- results(dds, name='expression_state_high_vs_low')
resOrdered <- res[order(res$padj),]

# Write the dge results:
write.table(resOrdered, DGE_OUTPUT, sep='\t', quote=F)
