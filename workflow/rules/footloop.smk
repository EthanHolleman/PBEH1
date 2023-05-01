

rule aggregate_genome_files:
    # combine all reference fasta files into an aggregated
    # fasta file for mapping
    conda:
        '../envs/py.yml'
    output:
        genome='output/mapping/genome/aggregatedGenome.fa',
        bed6='output/mapping/genome/aggregatedGenome.bed'
    params:
        genome_dir=config['GENOME_DIR']
    script:
        '../scripts/aggregateRefs.py'


def flow_cell_name_to_footloop_label(flow_cell):
    # Converts a flowcell label into a label that footloop likes
    d1 = flow_cell.split('-')[0][-1]
    d2 = flow_cell.split('-')[-1][-1]
    return f'PCB{d1}{d2}'



rule gzip_split_files:
    input:
        'output/groupedReads/{flow_cell}/{file_num}/{ref_name}.group.fastq'
    output:
        'output/groupedReads/{flow_cell}/{file_num}/{ref_name}.group.fastq'
    resources:
        mem_mb=2048,
        cpu=1,
        disk_mb=5000,
        time_hr=2
    shell:'''
    gzip {input}
    '''


# rule mapping:
#     # map all reads of a flowcell against all references included in the
#     # aggregated genome file
#     conda:
#         '../envs/footloop.yml'
#     input:
#         reads='output/blast/groupedReads/{flow_cell}/{file_num}/{ref_name}.group.fastq.gz',
#         genome=lambda wildcards: 
#             Path(config['GENOME_DIR']).joinpath(wildcards.ref_name),
#         index='output/mapping/genome/aggregatedGenome.bed',
#     output:
#         output_dir=directory('output/mapping/{flow_cell}/{file_num}/{ref_name}.group')
#     params:
#         label=lambda wildcards: 
#             flow_cell_name_to_footloop_label(wildcards.flow_cell)
#     shell:'''
#     submodules/footLoop/footLoop.pl -r {input.reads} \
#         -n {output.output_dir} \
#         -l {params.label} \
#         -g {input.genome} \
#         -i {input.index}  \
#         -Z
#     '''

rule mapping_sucks:
    # map all reads of a flowcell against all references included in the
    # aggregated genome file
    conda:
        '../envs/footloop.yml'
    input:
        reads='output/groupedReads/{flow_cell}/{file_num}',
        index='output/mapping/genome/aggregatedGenome.bed'
    output:
        output_dir=directory('output/mapping/{flow_cell}/{file_num}')
    params:
        label=lambda wildcards: 
            flow_cell_name_to_footloop_label(wildcards.flow_cell),
        genome_dir='../resources/referenceDNA/fasta',
        plasmid_id=config['PLASMID_ID_TABLE']
    script:'../scripts/mapper.py'

    

def agg_mapping(wildcards):
    checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
    file_num = glob_wildcards(
            Path(checkpoint_output).joinpath('{file_num}/sep.fastq')
        ).file_num
    return expand(
        'output/mapping/{flow_cell}/{file_num}',
        flow_cell=wildcards.flow_cell,
        file_num=file_num
    )

rule agg_all_mapping:
    input:
        agg_mapping
    output:
        'test.balls.{flow_cell}.txt'
    shell:''