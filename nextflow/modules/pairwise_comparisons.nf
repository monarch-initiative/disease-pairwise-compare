#!/usr/bin/env nextflow

process compute_pairwise_comps {
    tag 'compute_pairwise_comps'
    publishDir "./", mode: 'copy'

    input:
    path env_dir
    path data_dir

    output:
    path data_dir, emit: project_path
    val "done", emit: comp_sig

    script:
    """
    source ${env_dir}/.venv/bin/activate
    python disease-pairwise-compare/python/disease_pairwise_similarity.py -p ./disease-comparisons \
                                                                                  -c ${task.cpus} \
    """
}