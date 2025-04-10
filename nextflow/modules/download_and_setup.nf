#!/usr/bin/env nextflow

process download_and_setup {
    tag "get_monarch_kg_data"

    input:
    path env_dir

    output:
    path "pairwise-comparisons-data", emit: project_path

    script:
    """
    source ${env_dir}/.venv/bin/activate
    python pairwise-disease-compare/python/download_and_setup.py -p ./disease-comparisons
    """
}