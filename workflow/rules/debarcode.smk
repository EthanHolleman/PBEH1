

rule make_primer_fasta:
    conda:
        '../envs/py.yml'
    input:
        sample_df=config['SAMPLE_TABLE']
    output:
        'output/barcode/refs/primers.fa'
    script:'../scripts/primersToFasta.py'



rule blast_primers:
    conda:
        '../envs/blast.yml'
    input:
        query='output/barcode/refs/primers.fa',
        db='output/blast/dbs/{flow_cell}/{flow_cell}_{file_num}'
    output:
        "output/barcode/primerBLAST/{flow_cell}/{file_num}.blastn.search.tsv"
    params:
        db_path=lambda wildcards: f"output/blast/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}.blastn.db"
    threads: 6
    shell:
        """
    cat {input.query} | blastn -db {params.db_path} -outfmt 6 \
                        -perc_identity 80 -task blastn -num_threads {threads} \
                        > {output}
    """


rule assign_samples:
    conda:
        '../envs/py.yml'
    input:
        blast="output/barcode/primerBLAST/{flow_cell}/{file_num}.blastn.search.tsv",
        sample_table=config['SAMPLE_TABLE']
    output:
        'output/barcode/primerAssignment/{flow_cell}/{file_num}.primerAssignment.tsv'
    script:'../scripts/debarcodeBlast.py'



def agg_barcodes(wildcards):
    checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
    file_num = glob_wildcards(
            Path(checkpoint_output).joinpath('{file_num}/sep.fastq')
        ).file_num
    return expand(
        'output/barcode/primerAssignment/{flow_cell}/{file_num}.primerAssignment.tsv',
        flow_cell=wildcards.flow_cell,
        file_num=file_num
    )



rule concat_barcodes:
    input:
        agg_barcodes
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
