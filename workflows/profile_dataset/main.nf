nextflow.enable.types = true

include { AnalysedDataset ; Dataset ; Options ; RunInfo } from '../../common.nf'
include { ANALYSE } from '../../modules/local/analyse/main.nf'
include { CHECKSUM } from '../../modules/local/checksum/main.nf'
include { DESCRIBE } from '../../modules/local/describe/main.nf'
include { EXPORT } from '../../modules/local/export/main.nf'
include { PREPARE_TABLE } from '../../modules/local/prepare_table/main.nf'
include { VISUALISE } from '../../modules/local/visualise/main.nf'
include { WRANGLE } from '../../modules/local/wrangle/main.nf'

/*
 * Profile each dataset independently: prepare, wrangle, analyse, visualise, export,
 * checksum and describe. Records are re-assembled by joining on their 'id' field.
 */
workflow PROFILE_DATASET {
    take:
    datasets: Channel<Dataset>
    options: Options
    run: RunInfo

    main:
    prepared = PREPARE_TABLE(datasets)
    wrangled = WRANGLE(prepared.join(datasets, by: 'id'), options)
    annotated = wrangled.join(datasets, by: 'id')
    analysed = ANALYSE(annotated, options)
    figures = VISUALISE(analysed.join(annotated, by: 'id'))
    exported = EXPORT(annotated, options)
    checksums = CHECKSUM(exported)
    described = DESCRIBE(
        annotated.join(analysed, by: 'id').join(exported, by: 'id').join(checksums, by: 'id'),
        run,
    )
    profiled = wrangled
        .join(analysed, by: 'id')
        .join(figures, by: 'id')
        .join(exported, by: 'id')
        .join(checksums, by: 'id')
        .join(described, by: 'id')

    emit:
    results: Channel<AnalysedDataset> = profiled
}
