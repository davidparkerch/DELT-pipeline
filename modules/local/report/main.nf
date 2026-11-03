nextflow.enable.types = true

include { RunInfo ; shellQuote } from '../../../common.nf'

/*
 * One self-contained HTML report (jinja2 + inline plotly.js) covering every dataset.
 *
 * This process consumes versions.yml, which is built from the 'versions' topic, so it must
 * NOT emit to that topic itself: the topic would never close and the run would hang.
 */
process REPORT {
    label 'process_single'

    input:
    files: Bag<Path>
    versions: Path
    title: String
    run: RunInfo

    stage:
    stageAs versions, 'software_versions.yml'

    output:
    file('report.html')

    script:
    """
    delt_pipeline report \\
        --input-dir . \\
        --versions software_versions.yml \\
        --title ${shellQuote(title)} \\
        --pipeline ${shellQuote(run.pipeline)} \\
        --pipeline-version ${shellQuote(run.version)} \\
        --homepage ${shellQuote(run.homepage)} \\
        --output report.html
    """

    stub:
    """
    touch report.html
    """
}
