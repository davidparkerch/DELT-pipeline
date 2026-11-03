nextflow.enable.types = true

include { RunInfo ; shellQuote } from '../../../common.nf'

/*
 * FAIR metadata for one published hit table: schema.org JSON-LD (identifier, license, provenance,
 * variables, distribution) and a Frictionless data package for the exports. `exports_base` is the
 * relative URL from the published JSON-LD file to the exports folder.
 */
process DESCRIBE {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        top: Path,
        id_column: String,
        exports: Set<Path>,
        sha256sums: Path,
        exports_base: String,
        description: String?,
        source: String?,
        license: String?,
        keywords: String?
    )
    dictionary: Path
    run: RunInfo

    output:
    record(
        id: id,
        identifier: eval('cat identifier.txt'),
        metadata: file("${id}.metadata.jsonld"),
        datapackage: file('datapackage.json'),
    )

    topic:
    tuple('delt_pipeline', eval('delt_pipeline --version')) >> 'versions'

    script:
    def names = exports.collect { path -> path.name }.toSorted().join(' ')
    """
    delt_pipeline describe \\
        --input ${top} \\
        --id-column ${id_column} \\
        --id ${id} \\
        --exports ${names} \\
        --exports-base ${shellQuote(exports_base)} \\
        --sha256sums ${sha256sums} \\
        --dictionary ${dictionary} \\
        --description ${shellQuote(description)} \\
        --source ${shellQuote(source)} \\
        --keywords ${shellQuote(keywords)} \\
        --license ${shellQuote(license ?: run.data_license)} \\
        --creator-name ${shellQuote(run.creator_name)} \\
        --creator-orcid ${shellQuote(run.creator_orcid)} \\
        --pipeline ${shellQuote(run.pipeline)} \\
        --pipeline-version ${shellQuote(run.version)} \\
        --homepage ${shellQuote(run.homepage)} \\
        --jsonld ${id}.metadata.jsonld \\
        --datapackage datapackage.json \\
        --identifier identifier.txt
    """

    stub:
    """
    echo "urn:sha256:stub" > identifier.txt
    touch ${id}.metadata.jsonld datapackage.json
    """
}
