nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Export a published table (hits, building blocks) for downstream tools, as <stem>.<format>, in the formats
 * chosen by options.export_formats (Parquet, Arrow, CSV, TSV, JSON Lines, NumPy, svmlight).
 */
process EXPORT {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        top: Path,
        id_column: String,
        stem: String
    )
    options: Options

    output:
    record(id: id, exports: files("${stem}.*"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('numpy', eval("python -c 'import numpy; print(numpy.__version__)'")) >> 'versions'

    script:
    """
    delt_pipeline export \\
        --input ${top} \\
        --id-column ${id_column} \\
        --formats '${options.export_formats}' \\
        --prefix ${stem}
    """

    stub:
    """
    touch ${stem}.parquet
    """
}
