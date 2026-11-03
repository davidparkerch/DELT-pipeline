nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Choose the BB-A building blocks of the sub-library: every BB-A used by a held-out (validated)
 * compound, plus random BB-As up to options.n_bb_a. Run once for all screens.
 */
process SELECT_SUBLIBRARY {
    label 'process_single'

    input:
    synthons: Path
    heldout: Bag<Path>
    options: Options

    output:
    file('sublibrary.json')

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    def tables = heldout.collect { path -> path.name }.toSorted().join(' ')
    """
    delt_pipeline select-sublibrary \\
        --synthons ${synthons} \\
        --heldout ${tables} \\
        --n-bb-a ${options.n_bb_a} \\
        --seed ${options.seed} \\
        --output sublibrary.json
    """

    stub:
    """
    touch sublibrary.json
    """
}
