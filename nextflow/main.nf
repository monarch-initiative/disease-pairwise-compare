#!/usr/bin/env nextflow

include { get_env } from './modules/get_env.nf'
include { download_and_setup } from './modules/download_and_setup.nf'
include { compute_pairwise_comps } from './modules/pairwise_comparisons.nf'
include { cluster_and_plot } from './modules/cluster_and_plot.nf'


// Phenologs calculation parameters 
params.cpu_cores = 10 // Not actuall used (any more... config takes care of this)
params.sim_metric = "aic"


workflow {

    // Input parameters
    Channel.value(params.sim_metric).set{ sim_metric }

    // Setup data environment
    get_env()
    download_and_setup(get_env.out.env_path)

    // Compute pairwise comparisons
    compute_pairwise_comps(get_env.out.env_path, 
                           download_and_setup.out.project_path)

    // Cluster and plot
    cluster_and_plot(get_env.out.env_path, 
                     download_and_setup.out.project_path,
                     compute_pairwise_comps.out.comp_sig,
                     sim_metric)
}