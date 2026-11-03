nextflow.enable.types = true

include { Options ; shellQuote } from '../../../common.nf'

/*
 * Tier and flag the top compounds of every target screen: replicate and disynthon support,
 * counter-screen and bead (matrix) binders, and hits shared between target screens.
 * Runs once over all screens; files are staged as <id>.{top,scores,synthons}.parquet and
 * <id>.counts_qc.json, and `roles` lists every screen as 'id=role'.
 */
process TRIAGE {
    label 'process_medium'

    input:
    files: Bag<Path>
    roles: Bag<String>
    options: Options

    output:
    record(
        id: 'triage',
        hits: file('hits.parquet'),
        triage_qc: file('triage_qc.json'),
    )

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    def screens = roles
        .toSorted()
        .collect { entry ->
            def parts = entry.tokenize('=')
            "--screen ${shellQuote(parts[0])} ${shellQuote(parts[1])}"
        }
        .join(' ')
    """
    export POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline triage \\
        ${screens} \\
        --input-dir . \\
        --level ${options.ci_level} \\
        --compare-top-n ${options.top_n} \\
        --output hits.parquet \\
        --qc triage_qc.json
    """

    stub:
    """
    touch hits.parquet triage_qc.json
    """
}
