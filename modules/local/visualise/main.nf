nextflow.enable.types = true

include { shellQuote } from '../../../common.nf'

/*
 * Interactive plotly figures as JSON: the building-block cube of the top compounds (coloured by
 * triage tier), replicate agreement per compound and per disynthon, an MA plot and the validation
 * of every score against measured KD.
 */
process VISUALISE {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        target: String,
        scores: Path,
        top: Path,
        synthons: Path,
        validation: Path
    )
    hits: Path

    output:
    record(id: id, figures: file("${id}.figures.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('plotly', eval("python -c 'import plotly; print(plotly.__version__)'")) >> 'versions'

    script:
    """
    delt_pipeline visualise \\
        --scores ${scores} \\
        --top ${top} \\
        --synthons ${synthons} \\
        --hits ${hits} \\
        --id ${shellQuote(id)} \\
        --validation ${validation} \\
        --screen ${shellQuote(target)} \\
        --output ${id}.figures.json
    """

    stub:
    """
    touch ${id}.figures.json
    """
}
