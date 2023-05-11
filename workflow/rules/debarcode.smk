

rule make_primer_fasta:
    conda:
        '../envs/py.yml'
    input:
        sample_df=config['SAMPLE_TABLE']
    output:
        'output/barcode/refs/primers.fa'
    script:'../scripts/primersToFasta.py'


rule bowtie_primers:
    conda:
        '../envs/bowtie.yml'
    input:
        primers='output/barcode/refs/primers.fa',
        bowtie_db='output/bowtie/dbs/{flow_cell}/{flow_cell}_{file_num}'
    output:
        "output/barcode/primerBowtie/sam/{flow_cell}.{file_num}.sam"
    params:
        db_path=lambda wildcards: f"output/bowtie/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}"
    shell:'''
    bowtie2 -x {params.db_path} -f {input.primers} -S {output} -a
    '''


rule assign_primers:
    conda:
        '../envs/py.yml'
    input:
        sam="output/barcode/primerBowtie/sam/{flow_cell}.{file_num}.sam",
        reads='output/reads/seperatedTrim/{flow_cell}/{file_num}.prune.fastq',
        samples=config['SAMPLE_TABLE']
    output:
        'output/barcode/primerAssignment/{flow_cell}/{file_num}.primerAssignment.tsv'
    script:'../scripts/assignPrimers.py'


rule concat_barcodes:
    input:
        expand(
            'output/barcode/primerAssignment/{flow_cell}/{file_num}.primerAssignment.tsv',
            file_num=DIVS, allow_missing=True
        )
    output:
        'output/barcode/allBarcodes/{flow_cell}/{flow_cell}.all.barcodes.tsv'
    shell:'''
    cat {input} > {output}
    '''


rule plot_barcodes:
    conda:
        '../envs/R.yml'
    input:
        barcodes='output/barcode/allBarcodes/{flow_cell}/{flow_cell}.all.barcodes.tsv',
        read_lengths='output/reads/readLengths/{flow_cell}/all.lengths.tsv',
        exp_read_lengths='../resources/productLengths/productLengthsExpected.tsv'
    output:
        heatmaps='output/plots/barcodes/{flow_cell}/barcodes.heatmaps.{flow_cell}.png',
        length_dist='output/plots/barcodes/{flow_cell}/barcodes.readLength.{flow_cell}.png',
        counts='output/plots/barcodes/{flow_cell}/barcodes.count.{flow_cell}.png'
    script:'../scripts/plotBarcodes.R'


# rule blast_all_primers:
#     input:
#         agg_primer_blast
#     output:
#         'blast.primers.{flow_cell}'
#     shell:'''
#     touch {output}
#     '''
