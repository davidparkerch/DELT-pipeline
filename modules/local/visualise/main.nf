nextflow.enable.types = true

include { shellQuote } from '../../../common.nf'

/*
 * Interactive plotly figures as JSON: the building-block cube of the top compounds and the
 * validation of every score against measured KD.
 */
process VISUALISE {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        target: String,
        top: Path,
        validation: Path
    )

    output:
    record(id: id, figures: file("${id}.figures.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('plotly', eval("python -c 'import plotly; print(plotly.__version__)'")) >> 'versions'

    script:
    """
    delt_pipeline visualise \\
        --top ${top} \\
        --validation ${validation} \\
        --screen ${shellQuote(target)} \\
        --output ${id}.figures.json
    """

    stub:
    """
    touch ${id}.figures.json
    """
}
