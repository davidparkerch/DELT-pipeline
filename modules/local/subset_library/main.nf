nextflow.enable.types = true

/*
 * Stream one full KinDEL screen (Parquet) into the compact sub-library count table. The QC file
 * carries the whole-library depth of every count column, which normalises all later rates.
 */
process SUBSET_LIBRARY {
    tag id
    label 'process_high'

    input:
    record(
        id: String,
        table: Path
    )
    synthons: Path
    sublibrary: Path

    output:
    record(id: id, counts: file("${id}.counts.parquet"), counts_qc: file("${id}.counts_qc.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    """
    export POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline subset \\
        --input ${table} \\
        --synthons ${synthons} \\
        --sublibrary ${sublibrary} \\
        --output ${id}.counts.parquet \\
        --qc ${id}.counts_qc.json
    """

    stub:
    """
    touch ${id}.counts.parquet ${id}.counts_qc.json
    """
}
