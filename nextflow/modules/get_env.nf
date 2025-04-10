
process get_env {
    tag "get_env"
    publishDir "./", mode: 'copy'

    output:
    path "pairwise-disease-compare", emit: env_path

    script:
    """
    git clone git@github.com:monarch-initiative/pairwise-disease-compare.git
    cd pairwise-disease-compare
    poetry config virtualenvs.in-project true
    poetry install
    """
}