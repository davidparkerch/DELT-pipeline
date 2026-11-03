nextflow.enable.types = true

/*
 * Shared record types and helper functions, included by main.nf, workflows and modules.
 */

// One row of the input samplesheet: a feature table plus its descriptive metadata.
record Dataset {
    id: String
    table: Path
    id_column: String
    target_column: String?
    description: String?
    source: String?
    license: String?
    keywords: String?
    dictionary: Path?
}

// Analysis settings shared by every dataset (built from params in main.nf).
record Options {
    export_formats: String
    impute: String
    max_missing: Float
    n_components: Integer
    n_clusters: Integer
    cv_folds: Integer
    seed: Integer
}

// Who and what produced the outputs, recorded in the FAIR metadata.
record RunInfo {
    pipeline: String
    version: String
    homepage: String
    data_license: String
    creator_name: String
    creator_orcid: String
}

// Everything produced for one dataset.
record AnalysedDataset {
    id: String
    identifier: String
    clean: Path
    qc: Path
    embedding: Path
    stats: Path
    figures: Path
    exports: Set<Path>
    sha256sums: Path
    metadata: Path
    datapackage: Path
}

/*
 * Single-quote free text (titles, descriptions, names) for safe use in a shell command.
 * Prefer this over `stage: env` for user-supplied text: env values are not quoted
 * when tasks run in containers.
 */
def shellQuote(text: String?) -> String {
    return "'" + (text ?: '').replace("'", "'\"'\"'") + "'"
}
