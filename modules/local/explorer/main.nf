nextflow.enable.types = true

include { RunInfo ; shellQuote } from '../../../common.nf'

/*
 * The DELT Explorer: one static HTML page (plotly.js inline, RDKit.js from a CDN) to browse the triaged
 * hits, the resynthesis shortlist and the building blocks. Filters, a scatter plot and the building-block
 * cube with linked marking, substructure and similarity search, and a CSV export, all in the browser.
 * Runs once over staged files: hits.parquet, picks.parquet, blocks.parquet, triage.metadata.jsonld (the data
 * license is read from it) and one <screen>.blockset.parquet per screen.
 */
process EXPLORER {
    label 'process_single'

    input:
    files: Bag<Path>
    title: String
    run: RunInfo

    output:
    file('explorer.html')

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('rdkit', eval("python -c 'import rdkit; print(rdkit.__version__)'")) >> 'versions'

    script:
    def sets = files
        .findAll { path -> path.name.endsWith('.blockset.parquet') }
        .collect { path -> "--blockset ${shellQuote(path.name.replace('.blockset.parquet', ''))} ${shellQuote(path.name)}" }
        .toSorted()
        .join(' ')
    """
    delt_pipeline explore \\
        --hits hits.parquet \\
        --picks picks.parquet \\
        --blocks blocks.parquet \\
        ${sets} \\
        --metadata triage.metadata.jsonld \\
        --title ${shellQuote(title)} \\
        --report-url index.html \\
        --pipeline ${shellQuote(run.pipeline)} \\
        --pipeline-version ${shellQuote(run.version)} \\
        --homepage ${shellQuote(run.homepage)} \\
        --output explorer.html
    """

    stub:
    """
    touch explorer.html
    """
}
