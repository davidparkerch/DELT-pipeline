nextflow.enable.types = true

include { Options ; shellQuote } from '../../../common.nf'

/*
 * A fingerprint model against the same measured KD as the scores: a gradient-boosted Poisson regression of
 * the target counts on Morgan fingerprints, once ignoring the control and once weighted by it, under a random
 * and a disynthon split, next to KinDEL's published machine-learning results. Screens without held-out data
 * get an empty result list.
 */
process MODEL {
    tag id
    label 'process_medium'

    input:
    record(
        id: String,
        target: String,
        scores: Path,
        heldout_ondna: Path?,
        heldout_offdna: Path?
    )
    synthons: Path
    options: Options

    output:
    record(id: id, model: file("${id}.model.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('scikit-learn', eval("python -c 'import sklearn; print(sklearn.__version__)'")) >> 'versions'
    tuple('rdkit', eval("python -c 'import rdkit; print(rdkit.__version__)'")) >> 'versions'

    script:
    def ondna = heldout_ondna ? "--ondna ${heldout_ondna}" : ''
    def offdna = heldout_offdna ? "--offdna ${heldout_offdna}" : ''
    """
    export OMP_NUM_THREADS=${task.cpus} POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline model \\
        --scores ${scores} \\
        --synthons ${synthons} \\
        --target ${shellQuote(target.toLowerCase())} ${ondna} ${offdna} \\
        --seed ${options.seed} \\
        --output ${id}.model.json
    """

    stub:
    """
    touch ${id}.model.json
    """
}
