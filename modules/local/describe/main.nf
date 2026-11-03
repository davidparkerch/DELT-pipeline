nextflow.enable.types = true

include { RunInfo ; shellQuote } from '../../../common.nf'

/*
 * FAIR metadata for one dataset: schema.org JSON-LD (identifier, license, provenance,
 * variables, distribution) and a Frictionless data package for the exports.
 */
process DESCRIBE {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        clean: Path,
        stats: Path,
        exports: Set<Path>,
        sha256sums: Path,
        id_column: String,
        target_column: String?,
        description: String?,
        source: String?,
        license: String?,
        keywords: String?,
        dictionary: Path?
    )
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
    def target = target_column ? "--target-column '${target_column}'" : ''
    def names = exports.collect { path -> path.name }.toSorted().join(' ')
    def data_dictionary = dictionary ? "--dictionary ${dictionary}" : ''
    """
    delt_pipeline describe \\
        --input ${clean} \\
        --id-column '${id_column}' ${target} \\
        --id ${id} \\
        --exports ${names} \\
        --exports-base ../../exports/${id}/ \\
        --sha256sums ${sha256sums} \\
        --stats ${stats} ${data_dictionary} \\
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
