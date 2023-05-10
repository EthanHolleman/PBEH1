

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
    # aggregated genome file. It is CRITICAL that commands issues by snakemake
    # are executed in the output directory. This is because some program that is
    # run as a part of footloop.pl dumps files into the directory that the command
    # was issued from. If there are multible footloop.pl instances running
    # at the same time (lots of snakemake jobs) then the different jobs can end up
    # writing to the same file causing confusion and resulting in 0% mapping
    # Also need to be aware of the .fai index files that are created by bedToFasta command
    # as part of the footloop.pl script. If fasta header names change these need to be
    # deleted manually otherwise nothing will map. HEED THESE WARNINGS; they were learned
    # the hard way. 
    conda:
        '../envs/footloop.yml'
    input:
        reads='output/groupedReadsClipped/{flow_cell}/{file_num}/{plasmid}.gb.group.fastq',
        index='output/mapping/genome/aggregatedGenome.bed'
    output:
        output_dir=directory('output/mapping/{flow_cell}/{file_num}/{plasmid}'),
        done='output/mapping/{flow_cell}/{file_num}/{plasmid}/updatedMapping.done.touch'
        # extra file used to trigger remapping without deleting all output
    params:
        label=lambda wildcards: 
            flow_cell_name_to_footloop_label(wildcards.flow_cell),
        genome=lambda wildcards: f'output/mapping/shiftedRefs/{wildcards.plasmid}.shift.fasta',
        genome_index=lambda wildcards: f'../resources/referenceDNA/fasta/{wildcards.plasmid}.fasta.fai',
        plasmid_id=config['PLASMID_ID_TABLE'],
        rel_path='../../../../../',
        updatedMappingName='updatedMapping.done.touch',
        min_amplicon_length='40p'  # min percent length of read compared to amplicon length
    resources:
        time=360,  # 120 mins
        mem_mb=16000
    shell:'''
    mkdir -p {output.output_dir}
    cd {output.output_dir}
    {params.rel_path}submodules/footLoop/footLoop.pl -r {params.rel_path}{input.reads} -n {params.rel_path}{output.output_dir} \
    -l {params.label} -g {params.rel_path}{params.genome} -i {params.rel_path}{input.index} \
    -L {params.min_amplicon_length} -Z
    touch {params.updatedMappingName}
    '''


rule call_peak:
    conda:
        '../envs/footloop.yml'
    input:
        mapping_dir='output/mapping/{flow_cell}/{file_num}/{plasmid}',
        update='output/mapping/{flow_cell}/{file_num}/{plasmid}/updatedMapping.done.touch'
    output:
        peak_call = directory('output/peakCall/{flow_cell}/{file_num}/{plasmid}'),
    params:
        label=lambda wildcards: 
            flow_cell_name_to_footloop_label(wildcards.flow_cell),
        min_amplicon_length='40p'  # min percent length of read compared to amplicon length
    resources:
        time=60,
        mem_mb=10000
    shell:'''
    mkdir -p {output.peak_call}
    submodules/footLoop/footPeak.pl -n {input.mapping_dir} -o {output.peak_call}
    '''


rule assign_sample_to_peaks:
    conda:
        '../envs/py.yml'
    input:
        mapping='output/mapping/{flow_cell}/{file_num}/{plasmid}',
        samples='output/barcode/primerAssignment/{flow_cell}/{file_num}.primerAssignment.tsv'
    output:
        'output/peakMerge/{flow_cell}/{file_num}/{plasmid}.merge.tsv'
    script:'../scripts/assignSampleToPeaks.py'



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
