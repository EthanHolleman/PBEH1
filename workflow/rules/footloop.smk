

rule aggregate_genome_files:
    # combine all reference fasta files into an aggregated
    # fasta file for mapping
    conda:
        '../envs/py.yml'
    output:
        genome='output/mapping/genome/aggregatedGenome.fa',
        bed6='output/mapping/genome/aggregatedGenome.bed',
        shifted_dir=directory('output/mapping/shiftedRefs/')
    params:
        genome_dir=f"{config['GENOME_DIR']}/fasta"
    script:
        '../scripts/shiftFa.py'


def flow_cell_name_to_footloop_label(flow_cell):
    # Converts a flowcell label into a label that footloop likes
    d1 = flow_cell.split('-')[0][-1]
    d2 = flow_cell.split('-')[-1][-1]
    return f'PCB{d1}{d2}'



# rule gzip_split_files:
#     input:
#         'output/groupedReads/{flow_cell}/{file_num}/{ref_name}.group.fastq'
#     output:
#         'output/groupedReads/{flow_cell}/{file_num}/{ref_name}.group.fastq'
#     resources:
#         mem_mb=2048,
#         cpu=1,
#         disk_mb=5000,
#         time_hr=2
#     shell:'''
#     gzip {input}
#     '''


# rule map_reads:
#     # map all reads of a flowcell against all references included in the
#     # aggregated genome file
#     conda:
#         '../envs/footloop.yml'
#     input:
#         reads='output/groupedReads/{flow_cell}/{file_num}',
#         index='output/mapping/genome/aggregatedGenome.bed'
#     output:
#         output_dir=directory('output/mapping/{flow_cell}/{file_num}')
#     params:
#         label=lambda wildcards: 
#             flow_cell_name_to_footloop_label(wildcards.flow_cell),
#         genome_dir='../resources/referenceDNA/fasta',
#         plasmid_id=config['PLASMID_ID_TABLE']
#     resources:
#         time=120,  # 120 mins
#         mem_mb=16000
#     script:'../scripts/mapper.py'


rule trim_read_ends:
    conda:
        '../envs/blast.yml'
    input:
        reads='output/groupedReads/{flow_cell}/{file_num}/{plasmid}.gb.group.fastq'
    output:
        clipped='output/groupedReadsClipped/{flow_cell}/{file_num}/{plasmid}.gb.group.fastq'
    shell:'''
    seqtk trimfq -b 70 -e 70 {input.reads} > {output.clipped}
    '''


rule map_reads:
    # map all reads of a flowcell against all references included in the
    # aggregated genome file
    conda:
        '../envs/footloop.yml'
    input:
        reads='output/groupedReads/{flow_cell}/{file_num}/{plasmid}.gb.group.fastq',
        index='output/mapping/genome/aggregatedGenome.bed'
    output:
        output_dir=directory('output/mapping/{flow_cell}/{file_num}/{plasmid}')
    params:
        label=lambda wildcards: 
            flow_cell_name_to_footloop_label(wildcards.flow_cell),
        genome=lambda wildcards: f'output/mapping/shiftedRefs/{wildcards.plasmid}.shift.fasta',
        genome_index=lambda wildcards: f'../resources/referenceDNA/fasta/{wildcards.plasmid}.fasta.fai',
        plasmid_id=config['PLASMID_ID_TABLE'],
        rel_path='../../../../../'
    resources:
        time=360,  # 120 mins
        mem_mb=16000
    shell:'''
    mkdir -p {output.output_dir}
    cd {output.output_dir}
    {params.rel_path}submodules/footLoop/footLoop.pl -r {params.rel_path}{input.reads} -n {params.rel_path}{output.output_dir} \
    -l {params.label} -g {params.rel_path}{params.genome} -i {params.rel_path}{input.index} -Z
    '''
    # this command deletes the index file created by the fastaFromBed command
    # in case fasta name has changed. Bedtool will not detect this and then 
    # it will never find the correct fasta file. I learned this the
    # hard way :(


rule call_peak:
    conda:
        '../envs/footloop.yml'
    input:
        'output/mapping/{flow_cell}/{file_num}/{plasmid}'
    output:
        directory('output/peakCall/{flow_cell}/{file_num}/{plasmid}')
    resources:
        time=120,
        mem_mb=16000
    shell:'''
    submodules/footLoop/footPeak.pl -n {input} -o {output}
    '''
    

def agg_mapping(wildcards):
    checkpoint_output = checkpoints.group_plasmids.get(**wildcards).output[0]
    print(wildcards, "wildcards")
    plasmid = glob_wildcards(
        Path(checkpoint_output).joinpath("{plasmid}.gb.group.fastq")
    ).plasmid
    return expand(
        "output/peakCall/{flow_cell}/{file_num}/{plasmid}",
        flow_cell=wildcards.flow_cell,
        file_num=wildcards.file_num,
        plasmid=plasmid
    )

rule agg_all_mapping:
    input:
        agg_mapping
    output:
        'output/aggMapping/{flow_cell}.{file_num}.done'
    shell:'''
    touch {output}
    '''