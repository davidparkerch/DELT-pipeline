nextflow.enable.types = true

/*
 * The building-block resource: one row per building block (378 A, 1,128 B, 191 C in KinDEL) with RDKit
 * descriptors of the fragment, a Murcko scaffold, the blocks that encode the same fragment, and the
 * pre-selection load of the block's compounds as a synthesis-yield proxy. Runs once over all screens.
 */
process BUILDING_BLOCKS {
    label 'process_low'

    input:
    synthons: Path
    sublibrary: Path
    counts: Bag<Path>

    output:
    record(
        id: 'building_blocks',
        table: file('blocks.parquet'),
        blocks_qc: file('blocks_qc.json'),
    )

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('rdkit', eval("python -c 'import rdkit; print(rdkit.__version__)'")) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    def tables = counts.collect { path -> path.name }.toSorted().join(' ')
    """
    export POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline blocks \\
        --synthons ${synthons} \\
        --sublibrary ${sublibrary} \\
        --counts ${tables} \\
        --output blocks.parquet \\
        --qc blocks_qc.json
    """

    stub:
    """
    touch blocks.parquet blocks_qc.json
    """
}
