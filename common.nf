nextflow.enable.types = true

/*
 * Shared record types and helper functions, included by main.nf, workflows and modules.
 */

// One row of the input samplesheet: a DEL selection screen (KinDEL layout) plus its metadata.
record Screen {
    id: String
    table: Path
    target: String
    role: String
    heldout_ondna: Path?
    heldout_offdna: Path?
    description: String?
    source: String?
    license: String?
    keywords: String?
}

// Analysis settings shared by every screen (built from params in main.nf).
record Options {
    n_bb_a: Integer
    alpha: Float
    ci_level: Float
    top_n: Integer
    seed: Integer
    export_formats: String
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

// Everything produced for one screen.
record ScreenResult {
    id: String
    identifier: String
    counts_qc: Path
    scores: Path
    top: Path
    score_qc: Path
    validation: Path
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
