
rule unzipFastq:
	input:
		flow_cell_reads = lambda wildcards: config[wildcards.flow_cell]['fastq']
	output:
		'output/unzippedReads/{flow_cell}.fastq'
	shell:'''
	gzip -d -c {input.flow_cell_reads} > {output}
	'''

rule runFastQC:
	conda:
		'../envs/fastQC.yml'
	input:
		flow_cell_reads = 'output/unzippedReads/{flow_cell}.fastq'
	output:
		output_dir = directory('output/fastqc/{flow_cell}')
	threads: 6
	shell:'''
	mkdir {output.output_dir}
	fastqc -o {output.output_dir} -t {threads} {input.flow_cell_reads}
	'''

