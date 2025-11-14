#!/usr/bin/env nextflow

// Enable Nextflow DSL v2
nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// Parameter definitions
// These parameter lists define the grid of values that will be used to
// construct combinations of arguments for the analysis. Each combination
// will be fed into the `run_analysis` process.
// -----------------------------------------------------------------------------
params.positive_cutoffs = [ 0.1, 0.2,0.5,0.75]                     // thresholds for positive calls
params.fisher_cutoffs = [0.1, 0.01, 0.001,0.0001,1e-5,1e-10,1e-15,1e-20] // p-value cutoffs for Fisher test
params.wilcox_cutoffs = [0.1, 0.01, 0.001,0.0001,1e-5,1e-10,1e-15,1e-20] // p-value cutoffs for Wilcoxon test
params.compare = ["c_vs_ps", "c_vs_p", "p_vs_s"]                   // comparison types used by the R script
params.target = ["aa_j", "aa", "aa_v", "aa_j_v"]                   // target feature sets
params.meta_file = "meta.txt"                                      // metadata file path (downloaded from 10.5281/zenodo.13141051)
params.data_file = "final_samples_for_classification.csv.bz2"      // input data file path (downloaded from 10.5281/zenodo.13141051)
params.script_file = "select_LOOCV.R"                              // R script that runs the analysis


// -----------------------------------------------------------------------------
// Process: run_analysis
// This process runs the R script for a single combination of parameters.
// It publishes results into the `result/` directory (copy mode).
// Resource directives (module, executor, queue, time, memory, cpus) are set
// to match the compute environment (LSF cluster with 'medium' queue here).
// -----------------------------------------------------------------------------
process run_analysis {
    // publish output CSV/RDS files to result/ by copying them
    publishDir 'result/', mode: 'copy'

    // load R module in the cluster environment (adjust if your environment differs)
    module 'R'

    // execution engine and job parameters (LSF). (adjust if you do not run it on LSF)
    executor 'lsf'
    queue 'medium'
    time '24h'
    memory '256 GB'
    cpus 28

    // Inputs:
    // - tuple of values for positive_cutoff, fisher_cutoff, wilcox_cutoff, compare, target
    // - script: path to the R script (select_LOOCV.R)
    // - data_file: path to the data file (compressed CSV)
    // - meta_file: metadata file, explicitly named 'meta.txt' in the process (keeps name inside job)
    input:
    tuple val(positive_cutoff), val(fisher_cutoff), val(wilcox_cutoff), val(compare), val(target)
    path script
    path data_file
    path meta_file, name: 'meta.txt'

    // Outputs:
    // - any generated CSV and RDS files produced by the R script
    output:
    path "*.csv"
    path "*.RDS"

    // The script executed inside the job.
    // Order of arguments matches those expected by select_LOOCV.R.
    script:
    """
    Rscript ${script} ${positive_cutoff} ${fisher_cutoff} ${wilcox_cutoff} ${target} ${compare}
    """
}


// -----------------------------------------------------------------------------
// Workflow
// Create channels from the parameter values and input files, then generate
// all combinations of parameter values using chained `combine` calls.
// Finally call the run_analysis process for each combination.
// -----------------------------------------------------------------------------
workflow {
    // Create channels for the provided files (script, data, metadata)
    script_ch = Channel.fromPath(params.script_file)
    data_ch = Channel.fromPath(params.data_file)
    meta_ch = Channel.fromPath(params.meta_file)

    // Build the Cartesian product of parameter lists by chaining combine()
    // Note: The chaining style results in nested tuples; the final tuple
    // destructures in the process input declaration.
    param_combinations = Channel.from(params.positive_cutoffs)
        .combine(Channel.from(params.fisher_cutoffs))
        .combine(Channel.from(params.wilcox_cutoffs))
        .combine(Channel.from(params.compare))
        .combine(Channel.from(params.target))

    // Launch the run_analysis process for each parameter combination.
    // Use .first() on the file channels to pass a single path (the file) into
    // each process invocation rather than streaming multiple copies.
    run_analysis(param_combinations, script_ch.first(), data_ch.first(), meta_ch.first())
}
