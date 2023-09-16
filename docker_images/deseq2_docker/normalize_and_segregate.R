suppressMessages(suppressWarnings(
    library("DESeq2", character.only=T, warn.conflicts = F, quietly = T)))

###########################################################################

# args from command line:
args<-commandArgs(TRUE)

# the raw counts for the TCGA cohort, in TSV format
RAW_COUNT_MATRIX <- args[1]

# the full annotation matrix:
FULL_ANN_MATRIX <- args[2]

# the gene we are using to segment the samples
SELECTED_GENE <- args[3]

# the quantile for the lower and upper cutoffs
LOW_Q <- as.numeric(args[4])
HIGH_Q <- as.numeric(args[5])

# An annotation file giving the samples that are high or low
# on the distribution (as well as other clinical metadata)
ANN_OUTPUT <- args[6]

# file name for normalized counts
NORM_COUNTS_OUTPUT <- args[7]

###########################################################################

# Don't mangle the column names quite yet. Keep them as-is and then we can convert them after we
# create a column name mapping
counts <- read.table(RAW_COUNT_MATRIX, 
                    sep='\t', 
                    stringsAsFactors=F, 
                    row.names=1, 
                    header=T, 
                    check.names=F)
orig_colnames = colnames(counts)

# cast the column names to "proper" names for R. This way we don't run into trouble with
# any called functions
new_colnames = make.names(orig_colnames)
colnames(counts) = new_colnames

# Create a map of the column names so we can re-assign them at the end
colname_mapping = setNames(orig_colnames, new_colnames)



# Now run the normalization so we can find the "extremes"
# for the gene of interest:
deseq_sf = estimateSizeFactorsForMatrix(counts)
norm_mtx = sweep(counts,2,deseq_sf,'/')

# Note that so far we have ALL the samples normalized.
# At this point, we subset to keep only our cohort of interest.
# We do this by joining with our annotation matrix
# now load the full annotation matrix :
full_annotations = read.table(FULL_ANN_MATRIX, row.names=1, sep='\t', header=T)
rownames(full_annotations) = make.names(rownames(full_annotations))


# extract the expressions for the gene of interest and find
# the cutoffs for the desired quantiles. 
gene_array <- as.matrix(norm_mtx[SELECTED_GENE, rownames(full_annotations)])

cutoffs <- unname(
                  quantile(
                           gene_array, 
                           na.rm=T, 
                           probs=c(LOW_Q, HIGH_Q)
                  )
           )
low_cutoff <- cutoffs[1]
high_cutoff <- cutoffs[2]


# now get the samples corresponding to the low and high-expression groups
low_samples <- names(gene_array[,gene_array < low_cutoff])
high_samples <- names(gene_array[,gene_array > high_cutoff])



# complete the annotation dataframe so we can have for later
# in case we need it:
# create a "dummy" annotation matrix which we will later fill in
# with high or low expression status
ann_df <- data.frame(expression_state=rep(NA, length(gene_array)))
rownames(ann_df) <- rownames(full_annotations)
ann_df[low_samples, 'expression_state'] <- 'low'
ann_df[high_samples, 'expression_state'] <- 'high'
# drop uninteresting ones in the middle:
ann_df <- na.omit(ann_df)


# merge with the dataframe describing the expression status (high or low)
joined_ann = merge(ann_df, full_annotations, by=0)
rownames(joined_ann) <- colname_mapping[joined_ann[,'Row.names']]
joined_ann <- subset(joined_ann, select=-c(Row.names))

write.table(joined_ann, ANN_OUTPUT, sep='\t', quote=F)

# write normalized counts
colnames(norm_mtx)<-colname_mapping[colnames(norm_mtx)]
write.table(norm_mtx, NORM_COUNTS_OUTPUT, sep='\t', quote=F)