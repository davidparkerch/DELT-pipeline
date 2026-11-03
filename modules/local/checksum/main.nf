nextflow.enable.types = true

/*
 * SHA-256 checksums of the exported files (sha256sum format) for integrity and provenance.
 */
process CHECKSUM {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        exports: Set<Path>
    )

    output:
    record(id: id, sha256sums: file('SHA256SUMS'))

    topic:
    tuple('coreutils', eval("sha256sum --version | head -n1 | sed 's/.* //'")) >> 'versions'

    script:
    def names = exports.collect { path -> path.name }.toSorted().join(' ')
    """
    LC_ALL=C sha256sum ${names} > SHA256SUMS
    """

    stub:
    """
    touch SHA256SUMS
    """
}
