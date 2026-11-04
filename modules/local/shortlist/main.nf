nextflow.enable.types = true

include { ShortlistOptions } from '../../../common.nf'

/*
 * A resynthesis shortlist from the triaged hits: tier-1 hits are clustered by structure (Butina on Morgan
 * fingerprints) and picked across the clusters, best first, with a cap per cluster and per disynthon.
 * The QC file compares the shortlist with the naive top of the ranking.
 */
process SHORTLIST {
    label 'process_low'

    input:
    hits: Path
    picking: ShortlistOptions

    output:
    record(
        id: 'shortlist',
        picks: file('picks.parquet'),
        picks_qc: file('picks_qc.json'),
    )

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('rdkit', eval("python -c 'import rdkit; print(rdkit.__version__)'")) >> 'versions'
    tuple('polars', eval("python -c 'import polars; print(polars.__version__)'")) >> 'versions'

    script:
    """
    export POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline shortlist \\
        --hits ${hits} \\
        --n-picks ${picking.n_picks} \\
        --similarity ${picking.similarity} \\
        --max-per-cluster ${picking.max_per_cluster} \\
        --max-per-disynthon ${picking.max_per_disynthon} \\
        --output picks.parquet \\
        --qc picks_qc.json
    """

    stub:
    """
    touch picks.parquet picks_qc.json
    """
}
