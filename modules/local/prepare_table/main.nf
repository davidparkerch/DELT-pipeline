nextflow.enable.types = true

/*
 * Decompress (if gzipped) and normalise line endings of a raw table with standard Linux tools.
 * Binary formats (Parquet, Arrow) pass through unchanged.
 */
process PREPARE_TABLE {
    tag id
    label 'process_single'

    input:
    record(
        id: String,
        table: Path
    )

    output:
    record(id: id, csv: file("${id}.input.*"))

    topic:
    tuple('gzip', eval("gzip --version | head -n1 | sed 's/^gzip //'")) >> 'versions'

    script:
    def name = table.name.replaceFirst(/\.gz$/, '')
    def output = "${id}.input.${name.tokenize('.').last()}"
    """
    case "${table}" in
        *.gz) gzip -cd "${table}" ;;
        *) cat "${table}" ;;
    esac | case "${output}" in
        *.csv | *.tsv | *.txt) tr -d '\\r' ;;
        *) cat ;;
    esac > "${output}"
    """

    stub:
    def name = table.name.replaceFirst(/\.gz$/, '')
    """
    touch "${id}.input.${name.tokenize('.').last()}"
    """
}
