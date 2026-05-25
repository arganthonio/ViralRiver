nextflow.enable.dsl=2

workflow {

    if( !params.host_ref ) {
        error "Please provide host reference FASTA using --host_ref /path/to/hg38.fa"
    }

    Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq_1),
                file(row.fastq_2)
            )
        }
        .set { reads_ch }

    module1_out = MODULE1(reads_ch)

    candidate_reads_m2 = module1_out.candidate_reads
    candidate_reads_m3 = module1_out.candidate_reads

    MODULE2(candidate_reads_m2)
    MODULE3(candidate_reads_m3)

    MODULE4_DEPTH(module1_out.fastp_json.collect())
}
process MODULE1 {

    tag "$sample"

    publishDir "${params.outdir}/module1", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}/${sample}_candidate_reads_1.fq.gz"),
          path("${sample}/${sample}_candidate_reads_2.fq.gz"),
          emit: candidate_reads

    path("${sample}/fastp_reports/${sample}_fastp.json"),
          emit: fastp_json

    script:
    """
    mkdir -p input ${sample}

    ln -s ${read1} input/${sample}_1.fastq.gz
    ln -s ${read2} input/${sample}_2.fastq.gz

    ViralRiver.module1.sh \\
      -i input \\
      -o . \\
      -d ${params.kraken_db} \\
      -r ${params.host_ref} \\
      -t ${task.cpus} \\
      -x ${params.kraken_taxid}
    """
}

process MODULE2 {

    tag "$sample"

    publishDir "${params.outdir}/module2", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("${sample}/${sample}_viral_counts.tsv"), emit: minimap2_counts
    path("${sample}/${sample}_rescued_viral_reads.fasta"), emit: minimap2_fasta

    script:
    """
    mkdir -p ${sample}

    cp ${read1} ${sample}/${sample}_candidate_reads_1.fq.gz
    cp ${read2} ${sample}/${sample}_candidate_reads_2.fq.gz

    ViralRiver.module2.sh \\
      -i . \\
      -o . \\
      -v ${params.viral_fasta} \\
      -r ${params.host_ref} \\
      -t ${task.cpus} \\
      -q ${params.mapq_minimap2}
    """
}

process MODULE3 {

    tag "$sample"

    publishDir "${params.outdir}/module3", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("${sample}/${sample}_bowtie2_viral_counts.tsv"), emit: bowtie2_counts
    path("${sample}/${sample}_rescued_high_qual_reads.fasta"), emit: bowtie2_fasta

    script:
    """
    mkdir -p ${sample}

    cp ${read1} ${sample}/${sample}_candidate_reads_1.fq.gz
    cp ${read2} ${sample}/${sample}_candidate_reads_2.fq.gz

    ViralRiver.module3.sh \\
      -i . \\
      -o . \\
      -v ${params.viral_fasta} \\
      -t ${task.cpus} \\
      -q ${params.mapq_bowtie2}
    """
}

process MODULE4_DEPTH {

    publishDir "${params.outdir}/module4_depth", mode: 'copy'

    input:
    path fastp_jsons

    output:
    path("viralriver_read_depth.tsv"), emit: read_depth

    script:
    """
    mkdir -p module1_fastp_jsons

    cp ${fastp_jsons} module1_fastp_jsons/

    ViralRiver.module4.depth.sh \\
      -i module1_fastp_jsons \\
      -o viralriver_read_depth.tsv
    """
}
