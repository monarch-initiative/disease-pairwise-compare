
process get_env {
    tag "get_env"
    publishDir "./", mode: 'copy'

    output:
    path "disease-pairwise-compare", emit: env_path

    script:
    """
    git clone git@github.com:monarch-initiative/disease-pairwise-compare.git
    cd disease-pairwise-compare
    poetry config virtualenvs.in-project true
    poetry install
    """
}