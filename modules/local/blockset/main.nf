nextflow.enable.types = true

/*
 * Building-block set enrichment: rank a screen's compounds by Poisson lower bound and test, for every
 * building block, whether the compounds containing it rank above the rest (Mann-Whitney from rank
 * sums, Benjamini-Hochberg over all blocks). The gene-set-enrichment idea, applied to synthons.
 */
process BLOCKSET {
    tag id
    label 'process_low'

    input:
    record(
        id: String,
        scores: Path,
        top: Path
    )

    output:
    record(id: id, blockset: file("${id}.blockset.parquet"), blockset_qc: file("${id}.blockset_qc.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('scipy', eval("python -c 'import scipy; print(scipy.__version__)'")) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    """
    export POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline blockset \\
        --scores ${scores} \\
        --top ${top} \\
        --output ${id}.blockset.parquet \\
        --qc ${id}.blockset_qc.json
    """

    stub:
    """
    touch ${id}.blockset.parquet ${id}.blockset_qc.json
    """
}
