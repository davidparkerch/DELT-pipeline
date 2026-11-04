nextflow.enable.types = true

include { BlockResult ; DatasetInfo ; Options ; RunInfo ; Screen ; ScreenResult ; ShortlistOptions ; ShortlistResult ; TriageResult } from '../../common.nf'
include { AGGREGATE_SYNTHONS } from '../../modules/local/aggregate_synthons/main.nf'
include { BLOCKSET } from '../../modules/local/blockset/main.nf'
include { BUILDING_BLOCKS } from '../../modules/local/building_blocks/main.nf'
include { CHECKSUM } from '../../modules/local/checksum/main.nf'
include { DESCRIBE } from '../../modules/local/describe/main.nf'
include { EXPORT } from '../../modules/local/export/main.nf'
include { SCORE } from '../../modules/local/score/main.nf'
include { SELECT_SUBLIBRARY } from '../../modules/local/select_sublibrary/main.nf'
include { SHORTLIST } from '../../modules/local/shortlist/main.nf'
include { SUBSET_LIBRARY } from '../../modules/local/subset_library/main.nf'
include { TRIAGE } from '../../modules/local/triage/main.nf'
include { VALIDATE } from '../../modules/local/validate/main.nf'
include { VISUALISE } from '../../modules/local/visualise/main.nf'

/*
 * Score DEL selection screens: choose one sub-library for all screens, stream each screen into a
 * count table, score, aggregate and validate every compound, test every building block for enrichment,
 * triage the top compounds of the target screens against the counter-screens, then publish the hit
 * tables and the building-block resource as FAIR datasets.
 * Records are re-assembled by joining on their 'id' field.
 */
workflow DELT_SELECTION {
    take:
    screens: Channel<Screen>
    synthons: Path
    dictionary: Path
    options: Options
    run: RunInfo
    picking: ShortlistOptions

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
    blocksets = BLOCKSET(scored)

    // The building-block resource is built once; the pre-selection load is shared by all screens.
    blocks = BUILDING_BLOCKS(synthons, sublibrary, counts.map { r -> r.counts }.collect())

    // Triage looks at every screen at once: counter-screens flag artefacts among the targets' hits.
    triage_files = scored
        .join(counts, by: 'id')
        .join(aggregated, by: 'id')
        .flatMap { r -> [r.top, r.scores, r.synthons, r.counts_qc] }
        .collect()
    roles = screens.map { screen -> "${screen.id}=${screen.role}" }.collect()
    triaged = TRIAGE(triage_files, roles, options)
    shortlisted = SHORTLIST(triaged.map { r -> r.hits }, picking)

    figures = VISUALISE(
        scored.join(validated, by: 'id').join(aggregated, by: 'id').join(counts, by: 'id').join(screens, by: 'id'),
        triaged.map { r -> r.hits },
    )

    // Published tables: each screen's top compounds, the triaged hits of all targets and the building blocks.
    tables = scored
        .map { r -> record(id: r.id, top: r.top, id_column: 'molecule_hash', stem: "${r.id}_hits", exports_base: "../../exports/${r.id}/") }
        .mix(triaged.flatMap { r -> [record(id: r.id, top: r.hits, id_column: 'hit_id', stem: "${r.id}_hits", exports_base: "../exports/${r.id}/")] })
        .mix(blocks.flatMap { r -> [record(id: r.id, top: r.table, id_column: 'block_id', stem: r.id, exports_base: "../exports/${r.id}/")] })
        .mix(shortlisted.flatMap { r -> [record(id: r.id, top: r.picks, id_column: 'hit_id', stem: r.id, exports_base: "../exports/${r.id}/")] })
    info = screens
        .map { s -> screenInfo(s) }
        .mix(screens.collect().flatMap { all -> [triageInfo(all), blocksInfo(all), shortlistInfo(all)] })
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
        .join(blocksets, by: 'id')
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
    shortlist = shortlisted
        .flatMap { r -> [r] }
        .join(exported, by: 'id')
        .join(checksums, by: 'id')
        .join(described, by: 'id')
    block_results = blocks
        .flatMap { r -> [r] }
        .join(exported, by: 'id')
        .join(checksums, by: 'id')
        .join(described, by: 'id')

    emit:
    results: Channel<ScreenResult> = results
    triage: Channel<TriageResult> = triage
    blocks: Channel<BlockResult> = block_results
    shortlist: Channel<ShortlistResult> = shortlist
    sublibrary: Value<Path> = sublibrary
}

/*
 * FAIR metadata of a dataset derived from several screens: a license or source is inherited only when
 * every screen declares the same one (otherwise the run's default license applies, and no source is given).
 */
def derivedInfo(id: String, description: String, keywords: String, screens: Bag<Screen>) -> DatasetInfo {
    def licenses = screens.collect { s -> s.license }.toSet()
    def sources = screens.collect { s -> s.source }.toSet()
    return record(
        id: id,
        description: description,
        source: sources.size() == 1 ? sources.toList()[0] : null,
        license: licenses.size() == 1 ? licenses.toList()[0] : null,
        keywords: keywords,
    )
}

// The triaged hits, derived from the target screens (counter-screens are named in the description).
def triageInfo(screens: Bag<Screen>) -> DatasetInfo {
    def targets = screens.findAll { s -> s.role == 'target' }.collect { s -> s.id }.toSorted().join(', ')
    def counters = screens.findAll { s -> s.role == 'counter' }.collect { s -> s.id }.toSorted().join(', ')
    def description = "Triaged top compounds of the target screens ${targets}" + (counters ? " (counter-screens: ${counters})" : '')
    return derivedInfo('triage', description, 'DNA-encoded library;hit triage', screens)
}

// The resynthesis shortlist, derived from the triaged hits of the target screens.
def shortlistInfo(screens: Bag<Screen>) -> DatasetInfo {
    def targets = screens.findAll { s -> s.role == 'target' }.collect { s -> s.id }.toSorted().join(', ')
    def description = "Diverse resynthesis shortlist of the trusted hits of the target screens ${targets}: best hit of each structural cluster, capped per disynthon"
    return derivedInfo('shortlist', description, 'DNA-encoded library;hit triage;resynthesis;diversity', screens)
}

// The building-block resource, derived from the library tables the screens share.
def blocksInfo(screens: Bag<Screen>) -> DatasetInfo {
    def ids = screens.collect { s -> s.id }.toSorted().join(', ')
    def description = "Building blocks of the DNA-encoded library with RDKit descriptors, fragments shared between blocks and pre-selection load, from the screens ${ids}"
    return derivedInfo('building_blocks', description, 'DNA-encoded library;building blocks;cheminformatics;RDKit', screens)
}

// The descriptive metadata of one screen.
def screenInfo(s: Screen) -> DatasetInfo {
    return record(id: s.id, description: s.description, source: s.source, license: s.license, keywords: s.keywords)
}
