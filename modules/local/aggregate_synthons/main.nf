nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Pool compound counts per monosynthon (A, B, C) and disynthon (AB, AC, BC) and score every
 * synthon with the same exposure-mapped Poisson rate ratio as single compounds.
 */
process AGGREGATE_SYNTHONS {
    tag id
    label 'process_medium'

    input:
    record(
        id: String,
        scores: Path,
        counts_qc: Path
    )
    options: Options

    output:
    record(
        id: id,
        synthons: file("${id}.synthons.parquet"),
        synthon_qc: file("${id}.synthon_qc.json"),
    )

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    """
    export POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline synthons \\
        --scores ${scores} \\
        --counts-qc ${counts_qc} \\
        --alpha ${options.alpha} \\
        --level ${options.ci_level} \\
        --output ${id}.synthons.parquet \\
        --qc ${id}.synthon_qc.json
    """

    stub:
    """
    touch ${id}.synthons.parquet ${id}.synthon_qc.json
    """
}
