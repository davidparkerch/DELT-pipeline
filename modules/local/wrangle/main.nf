nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Clean a raw table with polars: normalise names, drop unusable features, impute, sort.
 */
process WRANGLE {
    tag id
    label 'process_low'

    input:
    record(
        id: String,
        csv: Path,
        id_column: String,
        target_column: String?
    )
    options: Options

    output:
    record(id: id, clean: file("${id}.clean.parquet"), qc: file("${id}.qc.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('python', eval("python --version | sed 's/^Python //'")) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    def target = target_column ? "--target-column '${target_column}'" : ''
    """
    export POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline wrangle \\
        --input ${csv} \\
        --id-column '${id_column}' ${target} \\
        --impute ${options.impute} \\
        --max-missing ${options.max_missing} \\
        --output ${id}.clean.parquet \\
        --qc ${id}.qc.json
    """

    stub:
    """
    touch ${id}.clean.parquet ${id}.qc.json
    """
}
