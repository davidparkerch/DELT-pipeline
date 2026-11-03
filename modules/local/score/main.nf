nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Score every compound with five enrichment methods (count difference, enrichment ratio, Poisson
 * lower bound, Faver z-score, empirical-Bayes shrunken ratio) and keep the top compounds.
 */
process SCORE {
    tag id
    label 'process_medium'

    input:
    record(
        id: String,
        counts: Path,
        counts_qc: Path
    )
    sublibrary: Path
    options: Options

    output:
    record(
        id: id,
        scores: file("${id}.scores.parquet"),
        top: file("${id}.top.parquet"),
        score_qc: file("${id}.score_qc.json"),
    )

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('scipy', eval("python -c 'import scipy; print(scipy.__version__)'")) >> 'versions'

    script:
    """
    export POLARS_MAX_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}
    delt_pipeline score \\
        --input ${counts} \\
        --counts-qc ${counts_qc} \\
        --sublibrary ${sublibrary} \\
        --alpha ${options.alpha} \\
        --level ${options.ci_level} \\
        --top-n ${options.top_n} \\
        --output ${id}.scores.parquet \\
        --top ${id}.top.parquet \\
        --qc ${id}.score_qc.json
    """

    stub:
    """
    touch ${id}.scores.parquet ${id}.top.parquet ${id}.score_qc.json
    """
}
