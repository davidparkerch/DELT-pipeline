nextflow.enable.types = true

include { Options ; RunInfo ; Screen ; ScreenResult } from '../../common.nf'
include { CHECKSUM } from '../../modules/local/checksum/main.nf'
include { DESCRIBE } from '../../modules/local/describe/main.nf'
include { EXPORT } from '../../modules/local/export/main.nf'
include { SCORE } from '../../modules/local/score/main.nf'
include { SELECT_SUBLIBRARY } from '../../modules/local/select_sublibrary/main.nf'
include { SUBSET_LIBRARY } from '../../modules/local/subset_library/main.nf'
include { VALIDATE } from '../../modules/local/validate/main.nf'
include { VISUALISE } from '../../modules/local/visualise/main.nf'

/*
 * Score DEL selection screens: choose one sub-library for all screens, stream each screen into a
 * count table, score and validate every compound, then publish the top compounds as FAIR datasets.
 * Records are re-assembled by joining on their 'id' field.
 */
workflow DELT_SELECTION {
    take:
    screens: Channel<Screen>
    synthons: Path
    dictionary: Path
    options: Options
    run: RunInfo

    main:
    heldout = screens
        .flatMap { screen -> [screen.heldout_ondna, screen.heldout_offdna] }
        .filter { path -> path != null }
        .collect()
    sublibrary = SELECT_SUBLIBRARY(synthons, heldout, options)
    counts = SUBSET_LIBRARY(screens, synthons, sublibrary)
    scored = SCORE(counts, sublibrary, options)
    validated = VALIDATE(scored.join(screens, by: 'id'))
    figures = VISUALISE(scored.join(validated, by: 'id').join(screens, by: 'id'))
    exported = EXPORT(scored, options)
    checksums = CHECKSUM(exported)
    described = DESCRIBE(
        scored.join(exported, by: 'id').join(checksums, by: 'id').join(screens, by: 'id'),
        dictionary,
        run,
    )
    results = counts
        .map { r -> record(id: r.id, counts_qc: r.counts_qc) }
        .join(scored, by: 'id')
        .join(validated, by: 'id')
        .join(figures, by: 'id')
        .join(exported, by: 'id')
        .join(checksums, by: 'id')
        .join(described, by: 'id')

    emit:
    results: Channel<ScreenResult> = results
    sublibrary: Value<Path> = sublibrary
}
