#!/usr/bin/env nextflow
/*
 * DELT-pipeline: Process DNA-Encoded Library Technology data and exports FAIR, analysis-ready datasets.
 * https://github.com/davidparkerch/delt-pipeline
 */
nextflow.enable.types = true

include { AnalysedDataset ; Dataset } from './common.nf'
include { PROFILE_DATASET } from './workflows/profile_dataset/main.nf'
include { REPORT } from './modules/local/report/main.nf'
include { WRITE_VERSIONS } from './modules/local/write_versions/main.nf'

params {
    // Samplesheet (CSV, JSON or YAML), one row per feature table; see assets/example/samplesheet.csv
    input: Channel<Dataset>
    // Comma-separated export formats: parquet, arrow, csv, tsv, jsonl, npz, svmlight
    export_formats: String = 'parquet,tsv,npz'
    // How to fill missing feature values: median, mean or drop
    impute: String = 'median'
    // Drop features with a larger fraction of missing values
    max_missing: Float = 0.5
    // Principal components to keep
    n_components: Integer = 2
    // Number of k-means clusters
    n_clusters: Integer = 3
    // Cross-validation folds for the supervised baseline
    cv_folds: Integer = 5
    // Random seed for every stochastic step
    seed: Integer = 42
    // Title of the HTML report
    report_title: String = 'DELT-pipeline report'
    // SPDX license applied to datasets that do not declare one
    data_license: String = 'CC-BY-4.0'
    // Data creator recorded in the FAIR metadata
    creator_name: String = 'David Parker'
    // ORCID iD of the data creator (optional)
    creator_orcid: String = '0009-0009-2897-5030'
}

workflow {

    main:
    options = record(
        export_formats: params.export_formats,
        impute: params.impute,
        max_missing: params.max_missing,
        n_components: params.n_components,
        n_clusters: params.n_clusters,
        cv_folds: params.cv_folds,
        seed: params.seed,
    )
    run = record(
        pipeline: workflow.manifest.name,
        version: workflow.manifest.version,
        homepage: workflow.manifest.homePage,
        data_license: params.data_license,
        creator_name: params.creator_name,
        creator_orcid: params.creator_orcid,
    )
    results = PROFILE_DATASET(params.input, options, run)
    versions = WRITE_VERSIONS(channel.topic('versions').collect())
    report_files = results
        .flatMap { dataset -> [dataset.qc, dataset.stats, dataset.figures, dataset.metadata] }
        .collect()
    report = REPORT(report_files, versions, params.report_title, run)

    publish:
    datasets = results
    report = report
    versions = versions
}

output {
    datasets: Channel<AnalysedDataset> {
        path { dataset ->
            dataset.clean >> "datasets/${dataset.id}/"
            dataset.qc >> "datasets/${dataset.id}/"
            dataset.embedding >> "datasets/${dataset.id}/"
            dataset.stats >> "datasets/${dataset.id}/"
            dataset.figures >> "datasets/${dataset.id}/"
            dataset.metadata >> "datasets/${dataset.id}/"
            dataset.exports >> "exports/${dataset.id}/"
            dataset.sha256sums >> "exports/${dataset.id}/"
            dataset.datapackage >> "exports/${dataset.id}/"
        }
        index {
            path 'datasets.json'
        }
    }
    report: Path {
        path '.'
    }
    versions: Path {
        path 'pipeline_info'
    }
}
