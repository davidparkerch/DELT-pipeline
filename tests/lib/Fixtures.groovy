import groovy.json.JsonSlurper

import java.nio.file.Files
import java.nio.file.Path

/**
 * Helpers shared by the nf-test suites.
 */
class Fixtures {

    /**
     * Files nf-prov copies into the output directory for the RO-Crate. Relative paths with these
     * names would resolve to the project's own files in snapshots, so they are listed separately.
     */
    static final List<String> CRATE_COPIES = ['main.nf', 'nextflow.config', 'README.md', 'ro-crate-metadata.json']

    /**
     * Link the project's assets/ into an nf-test launch directory, so the samplesheet's
     * relative paths (assets/example/...) resolve exactly as they do from the project root.
     */
    static void linkAssets(Object launchDir, Object projectDir) {
        Path link = Path.of(launchDir.toString(), 'assets')
        if( !Files.exists(link) )
            Files.createSymbolicLink(link, Path.of(projectDir.toString(), 'assets'))
    }

    /** The pipeline's own version, read from pyproject.toml. */
    static String ownVersion(Object projectDir) {
        def line = Path.of(projectDir.toString(), 'pyproject.toml').readLines().find { it.startsWith('version = ') }
        return line.split('"')[1]
    }

    /**
     * Tool versions as a sorted map. The pipeline's own package and version are masked so that
     * renaming or releasing does not invalidate every snapshot; Python is reduced to X.Y.
     * Accepts topic entries ([tool, version] lists) or a versions.yml path.
     */
    static Map versions(Object entries, Object projectDir) {
        def own = ownVersion(projectDir)
        def pairs = entries instanceof Collection
            ? entries.collect { entry -> [entry[0], entry[1]] }
            : Path.of(entries.toString()).readLines().findAll { it.contains(': ') }.collect { line ->
                def (tool, version) = line.split(': ', 2)
                [tool, version.replaceAll(/^'|'$/, '')]
            }
        return pairs.collectEntries { tool, version ->
            if( version == own )
                return ['<pipeline package>', '<pipeline version>']
            if( tool == 'python' )  // patch releases change with base-image rebuilds
                return [tool, version.tokenize('.').take(2).join('.')]
            return [tool, version]
        }.sort()
    }

    /**
     * Sort nested file collections so that unordered outputs (e.g. the Set<Path> from files())
     * snapshot identically on every run. Map keys are sorted too.
     */
    static Object canonical(Object value) {
        if( value instanceof Map )
            return new TreeMap(value.collectEntries { key, item -> [key.toString(), canonical(item)] })
        if( value instanceof Collection ) {
            def items = value.collect { item -> canonical(item) }
            // Records keep their order (sort them explicitly); files sort by name, not by work-dir path.
            return items.every { item -> item instanceof Map }
                ? items
                : items.sort { item -> item.toString().tokenize('/').last() }
        }
        return value
    }

    /**
     * Compact, stable summary of the HTML report (the file itself embeds plotly.js):
     * title, dataset sections, figure ids and whether it is fully offline.
     */
    static Map report(Object path) {
        def html = Path.of(path.toString()).text
        return [
            title: (html =~ /<title>(.*?)<\/title>/)[0][1],
            datasets: (html =~ /<section class="dataset" id="ds-([^"]+)"/).collect { match -> match[1] },
            figures: (html =~ /<div class="figure" id="(fig-[^"]+)"/).collect { match -> match[1] },
            offline: !html.contains('<script src'),
            bytes_over_1mb: html.length() > 1_000_000,
        ]
    }

    /** Parse a JSON file, dropping top-level keys that change with every release. */
    static Map json(Object path, List<String> drop = ['prov:wasGeneratedBy']) {
        def data = new JsonSlurper().parse(Path.of(path.toString()).toFile()) as Map
        return data.findAll { key, value -> !(key in drop) }
    }
}
