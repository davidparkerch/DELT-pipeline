nextflow.enable.types = true

include { shellQuote } from '../../../common.nf'

/*
 * Correlate every enrichment score with the measured KD of resynthesised compounds (on- and
 * off-DNA), next to KinDEL's published baselines. Screens without held-out data get an empty list.
 */
process VALIDATE {
    tag id
    label 'process_medium'

    input:
    record(
        id: String,
        target: String,
        scores: Path,
        score_qc: Path,
        heldout_ondna: Path?,
        heldout_offdna: Path?
    )

    output:
    record(id: id, validation: file("${id}.validation.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('scipy', eval("python -c 'import scipy; print(scipy.__version__)'")) >> 'versions'

    script:
    def ondna = heldout_ondna ? "--ondna ${heldout_ondna}" : ''
    def offdna = heldout_offdna ? "--offdna ${heldout_offdna}" : ''
    """
    delt_pipeline validate \\
        --scores ${scores} \\
        --score-qc ${score_qc} \\
        --target ${shellQuote(target.toLowerCase())} ${ondna} ${offdna} \\
        --output ${id}.validation.json
    """

    stub:
    """
    touch ${id}.validation.json
    """
}
