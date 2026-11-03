#!/usr/bin/env nextflow
/*
 * DELT-pipeline: Score DNA-encoded library (DEL) selections, validate the scores against measured affinities and publish FAIR hit tables.
 * https://github.com/davidparkerch/delt-pipeline
 */
nextflow.enable.types = true

include { Screen ; ScreenResult } from './common.nf'
include { DELT_SELECTION } from './workflows/delt_selection/main.nf'
include { REPORT } from './modules/local/report/main.nf'
include { WRITE_VERSIONS } from './modules/local/write_versions/main.nf'

params {
    // Samplesheet (CSV, JSON or YAML), one row per selection screen; see assets/example/samplesheet.csv
    input: Channel<Screen>
    // Folder with the building-block tables reacted_a.csv, reacted_b.csv and reacted_c.csv
    synthons: Path
    // Building blocks at position A in the sub-library (all held-out ones, plus random ones)
    n_bb_a: Integer = 60
    // Pseudocount of the enrichment ratio
    alpha: Float = 0.5
    // Confidence level of the Poisson rate-ratio interval
    ci_level: Float = 0.95
    // Compounds kept in each screen's published hit table
    top_n: Integer = 2000
    // Comma-separated export formats of the hit tables: parquet, arrow, csv, tsv, jsonl, npz, svmlight
    export_formats: String = 'parquet,tsv'
    // Random seed for every stochastic step
    seed: Integer = 42
    // Title of the HTML report
    report_title: String = 'DELT-pipeline report'
    // License applied to datasets that do not declare one (SPDX id or URL)
    data_license: String = 'CC-BY-4.0'
    // Data creator recorded in the FAIR metadata
    creator_name: String = 'David Parker'
    // ORCID iD of the data creator (optional)
    creator_orcid: String = '0009-0009-2897-5030'
}

workflow {

    main:
    options = record(
        n_bb_a: params.n_bb_a,
        alpha: params.alpha,
        ci_level: params.ci_level,
        top_n: params.top_n,
        seed: params.seed,
        export_formats: params.export_formats,
    )
    run = record(
        pipeline: workflow.manifest.name,
        version: workflow.manifest.version,
        homepage: workflow.manifest.homePage,
        data_license: params.data_license,
        creator_name: params.creator_name,
        creator_orcid: params.creator_orcid,
    )
    dictionary = file("${projectDir}/assets/hits_dictionary.csv")
    selection = DELT_SELECTION(params.input, params.synthons, dictionary, options, run)
    results = selection.results
    versions = WRITE_VERSIONS(channel.topic('versions').collect())
    report_files = results
        .flatMap { screen -> [screen.counts_qc, screen.score_qc, screen.validation, screen.figures, screen.metadata] }
        .collect()
    report = REPORT(report_files, versions, params.report_title, run)

    publish:
    screens = results
    sublibrary = selection.sublibrary
    report = report
    versions = versions
}

output {
    screens: Channel<ScreenResult> {
        path { screen ->
            screen.counts_qc >> "screens/${screen.id}/"
            screen.scores >> "screens/${screen.id}/"
            screen.top >> "screens/${screen.id}/"
            screen.score_qc >> "screens/${screen.id}/"
            screen.validation >> "screens/${screen.id}/"
            screen.figures >> "screens/${screen.id}/"
            screen.metadata >> "screens/${screen.id}/"
            screen.exports >> "exports/${screen.id}/"
            screen.sha256sums >> "exports/${screen.id}/"
            screen.datapackage >> "exports/${screen.id}/"
        }
        index {
            path 'screens.json'
        }
    }
    sublibrary: Path {
        path '.'
    }
    report: Path {
        path '.'
    }
    versions: Path {
        path 'pipeline_info'
    }
}
