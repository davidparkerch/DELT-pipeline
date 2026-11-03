nextflow.enable.types = true

/*
 * Interactive plotly figures (PCA, correlation, distributions, importance) as JSON.
 */
process VISUALISE {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        clean: Path,
        embedding: Path,
        stats: Path,
        id_column: String,
        target_column: String?
    )

    output:
    record(id: id, figures: file("${id}.figures.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('plotly', eval("python -c 'import plotly; print(plotly.__version__)'")) >> 'versions'

    script:
    def target = target_column ? "--target-column '${target_column}'" : ''
    """
    delt_pipeline visualise \\
        --input ${clean} \\
        --id-column '${id_column}' ${target} \\
        --embedding ${embedding} \\
        --stats ${stats} \\
        --output ${id}.figures.json
    """

    stub:
    """
    touch ${id}.figures.json
    """
}
