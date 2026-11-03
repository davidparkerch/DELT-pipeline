nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Export the top compounds (the published hit table) for downstream tools, in the formats
 * chosen by options.export_formats (Parquet, Arrow, CSV, TSV, JSON Lines, NumPy, svmlight).
 */
process EXPORT {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        top: Path
    )
    options: Options

    output:
    record(id: id, exports: files("${id}_hits.*"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('numpy', eval("python -c 'import numpy; print(numpy.__version__)'")) >> 'versions'

    script:
    """
    delt_pipeline export \\
        --input ${top} \\
        --id-column molecule_hash \\
        --formats '${options.export_formats}' \\
        --prefix ${id}_hits
    """

    stub:
    """
    touch ${id}_hits.parquet
    """
}
