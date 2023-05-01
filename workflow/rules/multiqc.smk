
rule run_multiqc:
    conda:
        '../envs/multiQC.yml'
    input:
        fastqc=expand('output/fastqc/{flow_cell}', flow_cell=FLOW_CELLS)
    output:
        output_dir=directory('output/multiqc')
    shell:'''
    multiqc . -o {output.output_dir}
    '''