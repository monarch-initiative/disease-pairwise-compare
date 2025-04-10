#!/usr/bin/env nextflow

process cluster_and_plot {
    tag 'cluster_and_plot'
    publishDir "./", mode: 'copy'

    input:
    path env_dir
    path data_dir
    val comp_sig
    val sim_metric

    output:
    path data_dir, emit: project_path

    script:
    """
    source ${env_dir}/.venv/bin/activate
    python disease-pairwise-compare/python/disease_cluster_and_plot.py -p ./disease-comparisons \
                                                                                  -m ${sim_metric} \
    """
}