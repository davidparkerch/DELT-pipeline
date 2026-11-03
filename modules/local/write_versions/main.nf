nextflow.enable.types = true

/*
 * Collate the (tool, version) pairs sent to the 'versions' topic into one sorted YAML file.
 * Runs inside Nextflow (exec), so it needs no container.
 */
process WRITE_VERSIONS {
    input:
    versions: Bag<Tuple<String, String>>

    output:
    file('versions.yml')

    exec:
    def lines = versions
        .collect { version -> "${version[0]}: '${version[1].trim()}'" }
        .toUnique()
        .toSorted()
    task.workDir.resolve('versions.yml').text = lines.join('\n') + '\n'
}
