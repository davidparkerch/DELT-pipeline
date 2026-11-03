nextflow.enable.types = true

include { Options } from '../../../common.nf'

/*
 * Statistics (scipy), PCA + k-means (scikit-learn) and a cross-validated baseline model.
 */
process ANALYSE {
    tag id
    label 'process_low'

    input:
    record(
        id: String,
        clean: Path,
        id_column: String,
        target_column: String?
    )
    options: Options

    output:
    record(id: id, embedding: file("${id}.embedding.parquet"), stats: file("${id}.stats.json"))

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'
    tuple('scipy', eval("python -c 'import scipy; print(scipy.__version__)'")) >> 'versions'
    tuple('scikit-learn', eval("python -c 'import sklearn; print(sklearn.__version__)'")) >> 'versions'

    script:
    def target = target_column ? "--target-column '${target_column}'" : ''
    """
    export OMP_NUM_THREADS=${task.cpus} POLARS_MAX_THREADS=${task.cpus}
    delt_pipeline analyse \\
        --input ${clean} \\
        --id-column '${id_column}' ${target} \\
        --n-components ${options.n_components} \\
        --n-clusters ${options.n_clusters} \\
        --cv-folds ${options.cv_folds} \\
        --seed ${options.seed} \\
        --threads ${task.cpus} \\
        --embedding ${id}.embedding.parquet \\
        --stats ${id}.stats.json
    """

    stub:
    """
    touch ${id}.embedding.parquet ${id}.stats.json
    """
}
