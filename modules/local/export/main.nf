nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Export the cleaned table for downstream tools: Parquet, Arrow, CSV, TSV, JSON Lines,
 * NumPy (npz) and svmlight, as chosen by options.export_formats.
 */
process EXPORT {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        clean: Path,
        id_column: String,
        target_column: String?
    )
    options: Options

    output:
    record(id: id, exports: files("${id}.*"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('numpy', eval("python -c 'import numpy; print(numpy.__version__)'")) >> 'versions'

    script:
    def target = target_column ? "--target-column '${target_column}'" : ''
    """
    delt_pipeline export \\
        --input ${clean} \\
        --id-column '${id_column}' ${target} \\
        --formats '${options.export_formats}' \\
        --prefix ${id}
    """

    stub:
    """
    touch ${id}.parquet
    """
}
