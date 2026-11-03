nextflow.enable.types = true

include { DatasetInfo ; Options ; RunInfo ; Screen ; ScreenResult ; TriageResult } from '../../common.nf'
include { AGGREGATE_SYNTHONS } from '../../modules/local/aggregate_synthons/main.nf'
include { CHECKSUM } from '../../modules/local/checksum/main.nf'
include { DESCRIBE } from '../../modules/local/describe/main.nf'
include { EXPORT } from '../../modules/local/export/main.nf'
include { SCORE } from '../../modules/local/score/main.nf'
include { SELECT_SUBLIBRARY } from '../../modules/local/select_sublibrary/main.nf'
include { SUBSET_LIBRARY } from '../../modules/local/subset_library/main.nf'
include { TRIAGE } from '../../modules/local/triage/main.nf'
include { VALIDATE } from '../../modules/local/validate/main.nf'
include { VISUALISE } from '../../modules/local/visualise/main.nf'

/*
 * Score DEL selection screens: choose one sub-library for all screens, stream each screen into a
 * count table, score, aggregate and validate every compound, triage the top compounds of the target
 * screens against the counter-screens, then publish the hit tables as FAIR datasets.
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
    aggregated = AGGREGATE_SYNTHONS(scored.join(counts, by: 'id'), options)
    validated = VALIDATE(scored.join(screens, by: 'id'))

    // Triage looks at every screen at once: counter-screens flag artefacts among the targets' hits.
    triage_files = scored
        .join(counts, by: 'id')
        .join(aggregated, by: 'id')
        .flatMap { r -> [r.top, r.scores, r.synthons, r.counts_qc] }
        .collect()
    roles = screens.map { screen -> "${screen.id}=${screen.role}" }.collect()
    triaged = TRIAGE(triage_files, roles, options)

    figures = VISUALISE(
        scored.join(validated, by: 'id').join(aggregated, by: 'id').join(screens, by: 'id'),
        triaged.map { r -> r.hits },
    )

    // Published hit tables: each screen's top compounds, plus the triaged hits of all targets.
    tables = scored
        .map { r -> record(id: r.id, top: r.top, id_column: 'molecule_hash', exports_base: "../../exports/${r.id}/") }
        .mix(triaged.flatMap { r -> [record(id: r.id, top: r.hits, id_column: 'hit_id', exports_base: "../exports/${r.id}/")] })
    info = screens
        .map { s -> screenInfo(s) }
        .mix(screens.collect().flatMap { all -> [triageInfo(all)] })
    exported = EXPORT(tables, options)
    checksums = CHECKSUM(exported)
    described = DESCRIBE(
        tables.join(exported, by: 'id').join(checksums, by: 'id').join(info, by: 'id'),
        dictionary,
        run,
    )
    results = counts
        .map { r -> record(id: r.id, counts_qc: r.counts_qc) }
        .join(scored, by: 'id')
        .join(aggregated, by: 'id')
        .join(validated, by: 'id')
        .join(figures, by: 'id')
        .join(exported, by: 'id')
        .join(checksums, by: 'id')
        .join(described, by: 'id')
    triage = triaged
        .flatMap { r -> [r] }
        .join(exported, by: 'id')
        .join(checksums, by: 'id')
        .join(described, by: 'id')

    emit:
    results: Channel<ScreenResult> = results
    triage: Channel<TriageResult> = triage
    sublibrary: Value<Path> = sublibrary
}

/*
 * FAIR metadata of the triaged hits, derived from the screens they come from: a license or source
 * is inherited only when every screen declares the same one (otherwise the run's default license
 * applies, and no source is given).
 */
def triageInfo(screens: Bag<Screen>) -> DatasetInfo {
    def targets = screens.findAll { s -> s.role == 'target' }.collect { s -> s.id }.toSorted().join(', ')
    def counters = screens.findAll { s -> s.role == 'counter' }.collect { s -> s.id }.toSorted().join(', ')
    def licenses = screens.collect { s -> s.license }.toSet()
    def sources = screens.collect { s -> s.source }.toSet()
    return record(
        id: 'triage',
        description: "Triaged top compounds of the target screens ${targets}" + (counters ? " (counter-screens: ${counters})" : ''),
        source: sources.size() == 1 ? sources.toList()[0] : null,
        license: licenses.size() == 1 ? licenses.toList()[0] : null,
        keywords: 'DNA-encoded library;hit triage',
    )
}

// The descriptive metadata of one screen.
def screenInfo(s: Screen) -> DatasetInfo {
    return record(id: s.id, description: s.description, source: s.source, license: s.license, keywords: s.keywords)
}
