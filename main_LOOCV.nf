#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define parameter combinations
params.positive_cutoffs = [ 0.1, 0.2,0.5,0.75]
params.fisher_cutoffs = [0.1, 0.01, 0.001,0.0001,1e-5,1e-10,1e-15,1e-20]
params.wilcox_cutoffs = [0.1, 0.01, 0.001,0.0001,1e-5,1e-10,1e-15,1e-20]
params.compare = ["c_vs_ps", "c_vs_p", "p_vs_s"]
params.target = ["aa_j", "aa", "aa_v", "aa_j_v"]
params.meta_file = "meta.txt"
params.data_file = "final_samples_for_classification.csv.bz2"
params.script_file = "select_LOOCV.R"


process run_analysis {
    publishDir 'result/', mode: 'copy'
    module 'R'
    executor 'lsf'
    queue 'medium'
    time '24h'
    memory '256 GB'
    cpus 28



    input:
    tuple val(positive_cutoff), val(fisher_cutoff), val(wilcox_cutoff), val(compare), val(target)
    path script
    path data_file
    path meta_file, name: 'meta.txt'

    output:
    path "*.csv"
    path "*.RDS"

    script:
    """
    Rscript ${script} ${positive_cutoff} ${fisher_cutoff} ${wilcox_cutoff} ${target} ${compare}
    """
}

workflow {
    script_ch = Channel.fromPath(params.script_file)
    data_ch = Channel.fromPath(params.data_file)
    meta_ch = Channel.fromPath(params.meta_file)

    param_combinations = Channel.from(params.positive_cutoffs)
        .combine(Channel.from(params.fisher_cutoffs))
        .combine(Channel.from(params.wilcox_cutoffs))
        .combine(Channel.from(params.compare))
        .combine(Channel.from(params.target))

    run_analysis(param_combinations, script_ch.first(), data_ch.first(), meta_ch.first())
}

