/*

##############################################################################
##############################################################################

 TARGET exploration pipeline:

 This pipeline allows us to parallelize a process over the entire TCGA
 project.

 Note that this "simplified" version does not include the large network tools
 like panda and lioness. Simply does the "basic" DGE processes.

##############################################################################
##############################################################################

*/

import java.text.SimpleDateFormat

/*
    START default params:
*/

// remove genes with cohort-wise means lower than this threshold:
// (where cohort is the subset of samples AFTER subsetting for our
// high and low expressed gene of interest)
params.min_reads = 5

params.gene = "ENSG00000169550"
params.low_q = 0.25
params.high_q = 0.75
params.output_dir = "results"

/*
    END default params:
*/

// Define an output directory based on a timestamp:
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd.HHmmss")
def timestamp =  sdf.format(date)
output_dir = params.output_dir + "/" + timestamp + "/" + params.gene + "_" + params.low_q + "_" + params.high_q


process cohort_selection {

    tag "Subset the data for our cohort of interest"
    publishDir "${output_dir}/annotations", mode:"copy"
    container "hsph-qbrc/target-pipeline/pandas"
    cpus 2
    memory '4 GB'

    input:
        path(full_annotations)

    output:
        path("subset_annotations.tsv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/cohort_selection.py \
            -i ${full_annotations} \
            -o subset_annotations.tsv \
            "inss_stage=Stage 4"
        """
}

process segregate_by_expression {

    tag "Segregate expression on $raw_counts"
    publishDir "${output_dir}/normalized_counts", mode:"copy", pattern: "deseq2_norm_counts.all.tsv"
    publishDir "${output_dir}/annotations", mode:"copy", pattern: "*annotations*"
    container "hsph-qbrc/target-pipeline/deseq2"
    cpus 2
    memory '8 GB'

    input:
        path(raw_counts)
        path(full_annotations)

    output:
        path("final_annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv")
        path("deseq2_norm_counts.all.tsv")

    script:
        """
        /opt/software/miniconda/envs/deseq2/bin/Rscript \
            /opt/software/scripts/normalize_and_segregate.R \
            ${raw_counts} \
            ${full_annotations} \
            ${params.gene} \
            ${params.low_q} \
            ${params.high_q} \
            final_annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv \
            deseq2_norm_counts.all.tsv
        """
}


process run_dge {

    tag "Run differential expression on $raw_counts"
    publishDir "${output_dir}/dge_results", mode:"copy", pattern: "deseq2_results*"
    container "hsph-qbrc/target-pipeline/deseq2"
    cpus 4
    memory '8 GB'

    input:
        path(raw_counts)
        path(annotations)

    output:
        path("deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv")

    script:
        """
        /opt/software/miniconda/envs/deseq2/bin/Rscript /opt/software/scripts/deseq2.R \
            ${raw_counts} \
            ${annotations} \
            deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv
        """
}

process prep_gsea_inputs {

    tag "Prep the GSEA inputs"
    publishDir "${output_dir}/gsea", mode:"copy"
    container "hsph-qbrc/target-pipeline/gsea"
    cpus 2
    memory '4 GB'

    input:
        path(norm_counts)
        path(annotations)
        path(dge_results)

    output:
        path("${gct_file}")
        path("${cls_file}")
        path("${params.gene}_${params.low_q}_${params.high_q}.rnk")

    script:
        def gct_file_template = "%s_%s_%s.high_vs_low.gct"
        def cls_file_template = "%s_%s_%s.high_vs_low.cls"
        gct_file = String.format(gct_file_template, params.gene, params.low_q, params.high_q)
        cls_file = String.format(cls_file_template, params.gene, params.low_q, params.high_q)
        """
        /usr/bin/python3 /opt/software/scripts/prep_files.py \
            -f ${norm_counts} \
            -a ${annotations} \
            -g ${gct_file} \
            -c ${cls_file} \
            -t ${params.min_reads}

        # Create the RNK file and run GSEA preranked
        /usr/bin/python3 /opt/software/scripts/create_rnk_file.py \
            -f ${dge_results} \
            -o ${params.gene}_${params.low_q}_${params.high_q}.rnk
        """

}


process run_gsea {

    tag "Run GSEA on ${label}"
    publishDir "${output_dir}/gsea/${label}/preranked", mode:"copy", pattern: "*preranked_results.zip"
    publishDir "${output_dir}/gsea/${label}/original", mode:"copy", pattern: "*gsea_results.zip"
    container "hsph-qbrc/target-pipeline/gsea"
    cpus 4
    memory '8 GB'

    input:
        path gct_file
        path cls_file
        path rnk_file
        tuple val(key), val(label), val(gmx_url)

    output:
        path("${params.gene}_${params.low_q}_${params.high_q}.${label}.gsea_results.zip")
        path("${params.gene}_${params.low_q}_${params.high_q}.${label}.gsea_preranked_results.zip")

    script:
        """
        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEA \
            -res "${gct_file}" \
            -cls "${cls_file}#high_versus_low" \
            -gmx  ${gmx_url} \
            -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -out /gsea/ \
            -rpt_label "${label}" \
            -zip_report true \
            -collapse Collapse \
            -mode Max_probe \
            -norm meandiv \
            -nperm 1000 \
            -permute phenotype \
            -rnd_seed timestamp \
            -rnd_type no_balance \
            -scoring_scheme weighted \
            -metric Signal2Noise \
            -sort real \
            -order descending \
            -create_gcts false \
            -create_svgs false \
            -include_only_symbols true \
            -make_sets true \
            -median false \
            -num 100 \
            -plot_top_x 20 \
            -save_rnd_lists false \
            -set_max 500 \
            -set_min 15

        /usr/bin/python3 /opt/software/scripts/move_final_files.py \
            -p "/gsea/*/*.zip" \
            -o ${params.gene}_${params.low_q}_${params.high_q}.${label}.gsea_results.zip


        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEAPreranked \
            -gmx ${gmx_url} \
            -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -rnk ${rnk_file} \
            -collapse Collapse \
            -mode Abs_max_of_probes \
            -norm meandiv \
            -nperm 1000 \
            -rnd_seed timestamp \
            -scoring_scheme weighted \
            -rpt_label "preranked" \
            -create_svgs false \
            -include_only_symbols true \
            -make_sets true \
            -plot_top_x 20 \
            -set_max 500 \
            -set_min 15 \
            -zip_report true \
            -out /gsea_preranked/

        /usr/bin/python3 /opt/software/scripts/move_final_files.py \
            -p "/gsea_preranked/*/*.zip" \
            -o ${params.gene}_${params.low_q}_${params.high_q}.${label}.gsea_preranked_results.zip
        """
}


process map_ensg_to_symbol {
    tag "Run ENSG to symbol gene mapping on $exp_mtx and $dge_results"
    publishDir "${output_dir}/dge_results", mode:"copy", pattern: "deseq2_results*"
    publishDir "${output_dir}/normalized_counts", mode:"copy", pattern: "deseq2_norm_counts.symbol_remapped.*"
    container "hsph-qbrc/target-pipeline/pandas"
    cpus 2
    memory '4 GB'

    input:
        path(exp_mtx)
        path(dge_results)

    output:
        path("deseq2_norm_counts.symbol_remapped.all.tsv")
        path("deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.symbol_remapped.tsv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${exp_mtx} \
            -m /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -o deseq2_norm_counts.symbol_remapped.all.tsv

        /usr/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${dge_results} \
            -m /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -o deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.symbol_remapped.tsv
        """
}


workflow {

    raw_count_ch = Channel.fromPath('nbl_counts.target-rnaseq.tsv')
    full_annotations_ch = Channel.fromPath('nbl_ann.target-rnaseq.tsv')

    // subset to the group of interest (e.g. Stage 4-- see process)
    subset_ann_ch = cohort_selection(full_annotations_ch)

    // given the high and low quantile thresholds, segregate into high and low expression groups:
    (curated_ann_ch, norm_counts_ch)=segregate_by_expression(raw_count_ch, subset_ann_ch)

    // run the DGE on those groups
    dge_results_ch = run_dge(raw_count_ch, curated_ann_ch)

    // create the GCT, CLS, and RNK files for GSEA:
    (gct_ch, cls_ch, rnk_ch) = prep_gsea_inputs(norm_counts_ch, curated_ann_ch, dge_results_ch)
    gct_ch = gct_ch.collect()
    cls_ch = cls_ch.collect()
    rnk_ch = rnk_ch.collect()
    
    // these are the key-value pairs defining the labels and paths to GMX files.
    // This lets us run GSEA over different sets:
    label_ch = Channel.of(['hallmark', 'HALLMARK'], ['c2', 'C2'])
    gmx_path_ch = Channel.of(
        ['hallmark', 'ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/h.all.v2023.1.Hs.symbols.gmt'],
        ['c2', 'ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/c2.all.v2023.1.Hs.symbols.gmt']
    )
    gsea_set_ch = label_ch.combine(gmx_path_ch, by:0)

    // finally run gsea:
    run_gsea(gct_ch, cls_ch, rnk_ch, gsea_set_ch)
    
    // for ease, convert to gene symbols:
    (norm_counts_symbol_remapped_ch, dge_results_remapped_ch) = map_ensg_to_symbol(norm_counts_ch, dge_results_ch)

}
