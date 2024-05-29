
params.input = '/home/horneflablinux/Documents/workflows/investigut_pipeline/my_test.fa' // input fasta file
params.output = ''
params.single = false // single mode
params.multi = false // multi mode
params.dmnd_db = '/home/horneflablinux/Documents/workflows/investigut_pipeline/output/5_diamonddb/my_dbbb.dmnd' // DIAMOND database
params.query_cover = 90.0 // DIAMOND query-coverage
params.subject_cover = 90.0 // DIAMOND subject-coverage
params.id = 90 // DIAMOND %id
params.threads = Runtime.runtime.availableProcessors() // maximum number available
params.help = false

if (params.s && params.m) {
    log.error("You cannot select both single and multi mode at the same time.")
    exit 1
}
if (params.low) {
    params.query_cover = 50.0
    params.subject_cover = 50.0
    params.id = 50
}

// Define a help message
def helpMessage = """
Usage: nextflow run investigut.nf --input <fasta_file> [--output <output_folder>] [--single | --mutli] [--dmnd_tsv_output {DMND_TSV_OUTPUT}] [--query-cover {NUM}] [--subject-cover {NUM}] [--id {NUM}] [--low] [--threads <num>]

Options:
  --input <fasta_file>        Input fasta file
  --output <output_folder>    Output folder (default: output/<current_date_time>_investigut_output)
  --single                    Single mode (default)
  --multi                     Multi mode
  --dmnd_tsv_output <file>    Optional DIAMOND .tsv output file to skip alignment
  --query_cover <num>         DIAMOND query-coverage (default 90.0)
  --subject_cover <num>       DIAMOND subject-coverage (default 90.0)
  --id <num>                  DIAMOND %id (default 90)
  --low                       Overrides DIAMOND query/subject coverage and %id, giving all a value of 50
  --threads <num>             Number of threads to use (default: ${params.threads} (all available))
  --help                      Display this help message
"""

// Display the help message if the --help option was provided
if (params.help) {
    println helpMessage
    exit 0
}

process DIAMOND_BLASTP {
    tag "$meta.id"
    label 'process_medium'
    cpus params.threads

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    val blast_columns

    output:
    tuple val(meta), path('*.txt')  , optional: true, emit: txt
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def columns = blast_columns ? "${blast_columns}" : ''
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi


    diamond \\
        blastp \\
        --threads ${task.cpus} \\
        --db ${db} \\
        --query ${fasta_name} \\
        --outfmt 6 ${columns} \\
        ${args} \\
        --out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}

workflow {
    // Define input with metadata
    def input_data = [meta: [id: 'sample1', source: 'experiment_1'], fasta: params.input]
    def db_data = [meta2: [id: 'db1', source: 'reference_db'], db: params.dmnd_db]

    // Run the DIAMOND_BLASTP process
    DIAMOND_BLASTP(input_data, db_data, 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp full_sseq')
}