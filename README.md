# TARGET data processing pipeline

This repository contains a basic nextflow pipeline for processing RNA-seq data from TARGET-NBL cohort, although the pipeline should work on any count matrix + annotation file (see below).

To use, you will need to have installed Nextflow software (https://www.nextflow.io/docs/latest/getstarted.html). 
You should also have Docker installed, since the pipeline components require containerized environments. 
By default, the Docker images refer to prebuilt images available from this repository.

Once installed, you will need two input files:
- A tab-delimited (TSV) raw counts file. This should have a first column containing Ensembl gene IDs and a first row including aliquot/sample IDs. All entries must be integers.
- A tab-delimited (TSV) annotation file. This should have a first column containing aliquot/sample IDs that have some intersection with the column headers of your count matrix. As written, the pipeline expects this file to have a column named `inss_stage` from which we extract Stage 4 subjects.
  However, additional key-value pairs may be used with adjustments to the `cohort_selection` process.

These two files are referenced in `pipeline.nf` as `raw_counts_ch` and `full_annotations_ch`, respectively. Change the path as required.

Additional parameters include:
- `params.gene`: the Ensembl ID of the gene of interest. The pipeline stratifies into low and high-expression cohorts based on the normalized expression of this gene.
- `params.low_q` and `params.high_q`: the quantiles ( $q \in (0,1)$ ) at which we stratify low and high-expressing cohorts.
- `params.output_dir`: The name of the results folder.

By default, we create result folders *inside* `params.output_dir` that are time-stamped so that subsequent executions are not overwritten.
