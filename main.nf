nextflow.enable.dsl=2

workflow {

    // Check for the updated parameter name: host_ref
    if( !params.host_ref ) {
        error "Please provide host reference FASTA using --host_ref /path/to/host.fa"
    }

    // 1. Detect the actual directory where the CSV file is located (e.g., .../ViralRiver-main2)
    def csvFile = file(params.samples)
    def csvDir = csvFile.parent

    Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample,
                file("${csvDir}/${row.fastq_1}"), // Force searching inside the directory where the CSV resides
                file("${csvDir}/${row.fastq_2}")
            )
        }
        .set { reads_ch }

    // Invoke the updated modules
    module1_out = MODULE1(reads_ch)

    module2_out = MODULE2(module1_out.candidate_reads) 
    
    module3_out = MODULE3(module1_out.candidate_reads) 

    MODULE4(module1_out.fastp_json.collect())
	
	all_mqc_files = module1_out.all_outputs
        .mix(module2_out.all_outputs, module3_out.all_outputs)
        .flatten()
        .filter { file -> !file.name.endsWith('.fq.gz') && !file.name.endsWith('.fasta') && !file.name.endsWith('.fa') }
        .collect()

    MODULE5(
        all_mqc_files,
        file("${projectDir}/bin/multiqc_config.yaml")
    )
}

process MODULE1 {

    tag "$sample"

    publishDir "${params.outdir}/module1", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    // Mandatory channels for subsequent Nextflow processes
    tuple val(sample),
          path("${sample}/${sample}_candidate_reads_1.fq.gz"),
          path("${sample}/${sample}_candidate_reads_2.fq.gz"),
          emit: candidate_reads

    path("${sample}/fastp_reports/${sample}_fastp.json"),
          emit: fastp_json

    // Rescue all files and subdirectories for publishDir
    path("${sample}/**"), emit: all_outputs

    script:
    """
    # Create the input directory expected by the script
    mkdir -p input

    # Create correct local symlinks pointing to valid Nextflow files
    ln -s ../${read1} input/${sample}_1.fastq.gz
    ln -s ../${read2} input/${sample}_2.fastq.gz

    # Convert to real absolute paths of the working environment using Bash \$PWD
    INPUT_DIR_ABS="\$PWD/input"
    OUTPUT_DIR_ABS="\$PWD"

    # Using \${file(...)} forces Nextflow to automatically resolve relative parameters 
    # into absolute paths relative to your launch directory before executing Bash.
    ViralRiver.module.1.sh \\
      -i "\${INPUT_DIR_ABS}" \\
      -o "\${OUTPUT_DIR_ABS}" \\
      -d "${file(params.kraken_db)}" \\
      -r "${file(params.host_ref)}" \\
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
    
    // Rescue all files and subdirectories for publishDir
    path("${sample}/**"), emit: all_outputs

    script:
    """
    # Create the sample directory expected by the script
    mkdir -p ${sample}
    
    # Create local symlinks pointing to the candidate reads instead of risking 'cp' path errors
    ln -s ../${read1} ${sample}/${sample}_candidate_reads_1.fq.gz
    ln -s ../${read2} ${sample}/${sample}_candidate_reads_2.fq.gz

    # Convert to real absolute paths of the working environment using Bash \$PWD
    INPUT_DIR_ABS="\$PWD"
    OUTPUT_DIR_ABS="\$PWD"

    # Using \${file(...)} to safely resolve reference paths from the isolated work folder
    ViralRiver.module.2.sh \\
      -i "\${INPUT_DIR_ABS}" \\
      -o "\${OUTPUT_DIR_ABS}" \\
      -v "${file(params.viral_fasta)}" \\
      -r "${file(params.host_ref)}" \\
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
    
    // Rescue all files and subdirectories for publishDir
    path("${sample}/**"), emit: all_outputs

    script:
    """
    # Create the sample directory expected by the script
    mkdir -p ${sample}
    
    # Create local symlinks pointing to the candidate reads instead of risking 'cp' path errors
    ln -s ../${read1} ${sample}/${sample}_candidate_reads_1.fq.gz
    ln -s ../${read2} ${sample}/${sample}_candidate_reads_2.fq.gz

    # Convert to real absolute paths of the working environment using Bash \$PWD
    INPUT_DIR_ABS="\$PWD"
    OUTPUT_DIR_ABS="\$PWD"

    # Using \${file(...)} to safely resolve reference paths from the isolated work folder
    ViralRiver.module.3.sh \\
      -i "\${INPUT_DIR_ABS}" \\
      -o "\${OUTPUT_DIR_ABS}" \\
      -v "${file(params.viral_fasta)}" \\
      -t ${task.cpus} \\
      -q ${params.mapq_bowtie2}
    """
}

process MODULE4 {

    publishDir "${params.outdir}/module4", mode: 'copy'

    input:
    path fastp_jsons

    output:
    path("viralriver_read_depth.tsv"), emit: read_depth
    
    // Rescue all files in the top-level process directory for publishDir
    path("*")

    script:
    """
    # Recreate the nested directory structure expected by the Module 4 script
    mkdir -p mock_results
    
    for json in ${fastp_jsons}; do
        # Extract the sample name (e.g., from SRR16948824_fastp.json we get SRR16948824)
        sample=\$(basename "\$json" _fastp.json)
        
        # Create the specific nested subdirectory structure for this sample
        mkdir -p "mock_results/\${sample}/fastp_reports"
        
        # Place the json report inside its expected subfolder
        cp "\$json" "mock_results/\${sample}/fastp_reports/"
    done

    # Define absolute paths for the execution script
    INPUT_DIR_ABS="\$PWD/mock_results"
    OUTPUT_FILE_ABS="\$PWD/viralriver_read_depth.tsv"

    # Execute the script using the recreated directory structure
    ViralRiver.module.4.sh \\
      -i "\${INPUT_DIR_ABS}" \\
      -o "\${OUTPUT_FILE_ABS}"
    """
}

process MODULE5 {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path 'mqc_inputs/*'
    path 'multiqc_config.yaml'

    output:
    path("viralriver_multiqc_report.html"), emit: html_report
    path("viralriver_multiqc_report_data")
    path("*")

    script:
    """
    python3 -c "
	import glob, os

	def process_files(files, out_filename):
		if not files: return
		all_v = set()
		data = {}
		for f in files:
			if 'summary_' in f: continue
			
			sample = os.path.basename(f).split('_')[0]
			data[sample] = {}
			with open(f) as fh:
				for line in fh:
					line = line.strip()
					if not line or line.startswith('#') or line.startswith('VIRUS_TAXON'):
						continue
					parts = line.split('\\t')
					if len(parts) >= 2:
						try:
							v = parts[0]
							c = int(parts[-1])
							if c > 0:
								data[sample][v] = c
								all_v.add(v)
						except ValueError:
							continue
		if not all_v: return
		sorted_v = sorted(list(all_v))
		
		with open(out_filename, 'w') as out:
			out.write('\\t'.join(['Sample'] + sorted_v) + '\\n')
			for s in sorted(data.keys()):
				row = [s] + [str(data[s].get(v, 0)) for v in sorted_v]
				out.write('\\t'.join(row) + '\\n')
				
	
	bowtie2_files = glob.glob('mqc_inputs/*_bowtie2_viral_counts.tsv')
	process_files(bowtie2_files, 'mqc_inputs/summary_bowtie2_mqc.tsv')

	mm2_files = [f for f in glob.glob('mqc_inputs/*_viral_counts.tsv') if '_bowtie2' not in f and 'summary_' not in f]
	process_files(mm2_files, 'mqc_inputs/summary_minimap2_mqc.tsv')
	"
	
    ViralRiver.module.5.sh \\
      -i mqc_inputs \\
      -o .
    """
}