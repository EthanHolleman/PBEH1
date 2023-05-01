from pathlib import Path
import glob



checkpoint break_up_fastq:
    conda:
        '../envs/py.yml'
    input:
        flow_cell_reads = 'output/unzippedReads/{flow_cell}.fastq'
    output:
        directory('output/reads/seperated/{flow_cell}')
    params:
        lines_per_file=40000,  # 4 lines per fastq entry
        flow_cell=lambda wildcards:
            wildcards.flow_cell,
        filename='sep.fastq'
    shell:'''
    mkdir -p {output}/temp
    split -d -l {params.lines_per_file} \
        {input} {output}/temp/
    python scripts/processSplit.py {output}/temp {params.filename}
    '''


rule prune_reads_by_length:
    conda:
        '../envs/blast.yml'
    input:
        "output/reads/seperated/{flow_cell}/{file_num}/sep.fastq"
    output:
        "output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq"
    shell:'''
    seqtk trimfq -b 1000 -e 2700 {input} > {output}
    '''



rule count_reads_and_length:
    conda:
        '../envs/py.yml'
    input:
        "output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq"
    output:
        read_counts='output/reads/readCounts/{flow_cell}/sep/{file_num}.count.tsv',
        read_lengths='output/reads/readLengths/{flow_cell}/sep/{file_num}.lengths.tsv'
    params:
        flow_cell=lambda wildcards: wildcards.flow_cell
    script:'../scripts/countUniqueReads.py'


rule count_reads_no_plasmid_id:
    conda:
        '../envs/py.yml'
    input:
        'output/plasmidID/noID/{flow_cell}/{file_num}.no.ID.fastq'
    output:
        read_counts='output/plasmidID/noID/{flow_cell}/readCounts/{file_num}.no.ID.count.tsv',
        read_lengths='output/plasmidID/noID/{flow_cell}/readLengths/{file_num}.no.ID.lengths.tsv'
    params:
        flow_cell=lambda wildcards: wildcards.flow_cell
    script:
        '../scripts/countUniqueReads.py'

def agg_count_reads(wildcards):
    checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
    print(wildcards, 'wildcards')
    file_num = glob_wildcards(
            Path(checkpoint_output).joinpath('{file_num}/sep.fastq')
        ).file_num
    return expand(
        'output/reads/readCounts/{flow_cell}/sep/{file_num}.count.tsv',
        flow_cell=wildcards.flow_cell,
        file_num=file_num
    )


def agg_count_reads_no_id(wildcards):
    checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
    print(wildcards, 'wildcards')
    file_num = glob_wildcards(
            Path(checkpoint_output).joinpath('{file_num}/sep.fastq')
        ).file_num
    return expand(
        'output/plasmidID/noID/{flow_cell}/readLengths/{file_num}.no.ID.lengths.tsv',
        flow_cell=wildcards.flow_cell,
        file_num=file_num
    )


rule concat_no_plasmid_id_lengths:
    input:
        agg_count_reads_no_id
    output:
        'output/plasmidID/noID/{flow_cell}/allReadLengths/all.no.ID.lengths.tsv'
    shell:'''
    cat {input} > {output}
    '''


rule plot_no_plasmid_id_lengths:
    input:
        read_lengths='output/plasmidID/noID/{flow_cell}/allReadLengths/all.no.ID.lengths.tsv'
    output:
        ''
    script:'../scripts/plotNoIDLengths.R'

        



rule concat_read_counts:
    input:
        agg_count_reads
    output:
        'output/reads/readCounts/{flow_cell}/{flow_cell}.all.count.tsv'
    shell:'''
    cat {input} > {output}
    '''

rule concat_flow_cell_counts:
    input:
        expand(
             'output/reads/readCounts/{flow_cell}/{flow_cell}.all.count.tsv',
             flow_cell=FLOW_CELLS
        )
    output:
        'output/reads/readCounts/allFlowCells.count.tsv'
    shell:'''
    cat {input} > {output}
    '''


rule plot_read_counts:
    conda:
        '../envs/R.yml'
    input:
        'output/reads/readCounts/allFlowCells.count.tsv'
    output:
        'output/plots/readCounts/totalReadCounts.png'
    script:'../scripts/plotReadCounts.R'



rule fastq_to_fasta:
    conda:
        '../envs/blast.yml'
    input:
        "output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq"
    output:
        'output/reads/seperatedFasta/{flow_cell}/{file_num}/sep.fasta'
    shell:'''
    seqtk seq -a {input} > {output}
    '''


rule make_blast_db:
    # build blastn databases for the split reads
    conda:
        "../envs/blast.yml"
    input:
        seqs="output/reads/seperatedFasta/{flow_cell}/{file_num}/sep.fasta",
    output:
        directory("output/blast/dbs/{flow_cell}/{flow_cell}_{file_num}"),
    params:
        db_name=lambda wildcards: f"{wildcards.flow_cell}_{wildcards.file_num}",
        db_path=lambda wildcards: f"output/blast/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}.blastn.db",
    shell:
        """
    mkdir -p {output}
    makeblastdb -in {input.seqs} -parse_seqids \
    -title "{params.db_name}" -dbtype nucl -out {params.db_path}
    """


rule local_blast_plasmids:
    # blast the plasmid identifying sequences against the database
    # of reads
    conda:
        "../envs/blast.yml"
    input:
        db="output/blast/dbs/{flow_cell}/{flow_cell}_{file_num}",
        query=config["PLASMID_ID_SEQS"]
    output:
        "output/plasmidID/search/{flow_cell}/{file_num}/blastn.search.tsv"
    params:
        db_path=lambda wildcards: f"output/blast/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}.blastn.db",
    shell:
        """
    cat {input.query} | blastn -db {params.db_path} -outfmt 6 -perc_identity 85 -task blastn > {output}
    """


rule identify_plasmid:
    conda:
        '../envs/py.yml'
    input:
        blast_tsv='output/plasmidID/search/{flow_cell}/{file_num}/blastn.search.tsv'
    output:
        'output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv'
    params:
        ref_dir=config['GENOME_DIR'],
        id_df=config['PLASMID_ID_TABLE']
    script:'../scripts/assignPlasmidID.py'


rule merge_read_length_plasmid_id:
    conda:
        '../envs/py.yml'
    input:
        plasmid_id= 'output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv',
        lengths='output/reads/readLengths/{flow_cell}/sep/{file_num}.lengths.tsv'
    output:
        'output/plasmidID/{flow_cell}/mergedIDs/{flow_cell}.{file_num}.lengths.ID.tsv'
    script:'../scripts/addLengthToIDs.py'


def agg_read_lengths(wildcards):
    checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
    print(wildcards, 'wildcards')
    file_num = glob_wildcards(
            Path(checkpoint_output).joinpath('{file_num}/sep.fastq')
        ).file_num
    return expand(
        'output/plasmidID/{flow_cell}/mergedIDs/{flow_cell}.{file_num}.lengths.ID.tsv',
        flow_cell=wildcards.flow_cell,
        file_num=file_num
    )


rule concat_read_lengths:
    input:
        agg_read_lengths
    output:
        'output/reads/readLengths/{flow_cell}/all.lengths.tsv'
    shell:'''
    cat {input} > {output}
    '''



def agg_plasmid_ids(wildcards):
    checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
    print(wildcards, 'wildcards')
    file_num = glob_wildcards(
            Path(checkpoint_output).joinpath('{file_num}/sep.fastq')
        ).file_num
    return expand(
        'output/plasmidID/{flow_cell}/mergedIDs/{flow_cell}.{file_num}.lengths.ID.tsv',
        flow_cell=wildcards.flow_cell,
        file_num=file_num
    )


def agg_plasmid_no_ids(wildcards):
    checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
    print(wildcards, 'wildcards')
    file_num = glob_wildcards(
            Path(checkpoint_output).joinpath('{file_num}/sep.fastq')
        ).file_num
    return expand(
        'output/plasmidID/noID/{flow_cell}/{file_num}.no.ID.fastq',
        flow_cell=wildcards.flow_cell,
        file_num=file_num
    )
    

rule group_no_id_read:
    conda:
        '../envs/py.yml'
    input:
        plasmid_id='output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv',
        reads='output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq'
    output:
        'output/plasmidID/noID/{flow_cell}/{file_num}.no.ID.fastq'
    script:'../scripts/findUnIDPlasmids.py'


rule concat_plasmid_ids:
    input:
        agg_plasmid_ids
    output:
        'output/plasmidID/allPlasmidID/{flow_cell}/{flow_cell}.concat.plasmidID.tsv'
    shell:'''
    cat {input} > {output}
    '''


rule agg_plasmid_no_ids:
    input:
        agg_plasmid_no_ids
    output:
        'output/plasmidID/allPlasmidNoIDs/{flow_cell}/{flow_cell}.concat.no.ID.tsv'
    shell:'''
    cat {input} > {output}
    '''



rule plot_plasmid_id_count:
    conda:
        '../envs/R.yml'
    input:
        plasmid_id='output/plasmidID/allPlasmidID/{flow_cell}/{flow_cell}.concat.plasmidID.tsv',
        no_plasmid_id='output/plasmidID/allPlasmidNoIDs/{flow_cell}/{flow_cell}.concat.no.ID.tsv',
        read_lengths='../resources/productLengths/productLengthsExpected.tsv'
    output:
        countPlot='output/plots/plasmidID/{flow_cell}.plasmid.count.png',
        lengthPlot='output/plots/plasmidID/{flow_cell}.plasmid.length.count.png'
    script:'../scripts/plotPlasmidID.R'


rule group_plasmids:
    conda:
        '../envs/py.yml'
    input:
        plasmid_id_df='output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv',
        reads='output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq',
        id_count='output/plots/plasmidID/{flow_cell}.plasmid.count.png',
        read_count='output/plots/readCounts/totalReadCounts.png'
    output:
        directory('output/groupedReads/{flow_cell}/{file_num}/')
    script:'../scripts/groupReadsByPlasmid.py'





















































# checkpoint copy:
#     input:
#         'output/groupedReads/{flow_cell}/{file_num}'
#     output:
#         '.plasmidCopy'
#     shell:'''
#     touch {output}
#     '''

# # checkpoint code to read count and specify all the outputs
# class Checkpoint_MakePattern:
#     def __init__(self, pattern):
#         self.pattern = pattern

#     def get_names(self):
#         with open('names.csv', 'rt') as fp:
#             names = [ x.rstrip() for x in fp ]
#         return names

#     def __call__(self, w):
#         global checkpoints

#         # wait for the results of 'check_csv'; this will trigger an
#         # exception until that rule has been run.
#         checkpoints.check_csv.get(**w)

#         # the magic, such as it is, happens here: we create the
#         # information used to expand the pattern, using arbitrary
#         # Python code.
#         names = self.get_names()

#         pattern = expand(self.pattern, name=names, **w)
#         return pattern


# def aggregate(wildcards):
#     #outputs_i = glob.glob(f"{checkpoints.break_up_fastq.get(flow_cell=wildcards.flow_cell).output}/*/")
#     #print(outputs_i, 'test')
#     c = checkpoints.break_up_fastq.get(**wildcards).output[0]
#     #outputs_i = glob.glob(f"{checkpoints.first.get().output}/*/")

#     # iterate through the split file directory to get the filenames
#     split_path = Path(f'output/reads/seperated/{wildcards.flow_cell}')
#     split_nums = [p.stem for p in split_path.iterdir()]
#     for each_num in split_nums:
#         outputs_group = checkpoints.groupPlasmids.get(file_num=each_num).output[0]
#         print(outputs_group)
#     return ''



#     # outputs_i = glob.glob(f"{c}/*/")
#     # outputs_i = [output.split('/')[-2] for output in outputs_i]
#     # split_files = []
#     # for i in outputs_i:
#     #     outputs_j = glob.glob(f"{checkpoints.second.get(i=i).output}/*/")
#     #     outputs_j = [output.split('/')[-2] for output in outputs_j]
#     #     for j in outputs_j:
#     #         split_files.append(f"copy/{i}/{j}/test2.txt")

#     # return split_files



# # NOTE: The script that groups files should create 1 file for each
# # reference plasmid for each sub fastq file that is produced even if no 
# # reads are assigned to a given plasmid. This was done so we know what the output
# # of this rule will be before it is run. I needed to do this because I was
# # having huge problems trying to get snakemake checkpoints to work when 
# # unknown outputs also spawn other unknown outputs. So the way I 

# def aggregate_grouped_reads(wildcards):
#     checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
#     return expand(
#         'output/blast/groupedReads/{flow_cell}/{file_num}',
#         flow_cell=wildcards.flow_cell,
#         file_num=glob_wildcards(
#             str(Path(checkpoint_output).joinpath(
#                 f'{wildcards.flow_cell}_{{file_num}}.fastq')
#                 )
#         ).file_num
#     )


# checkpoint want:
#     input:
#         aggregate
#     output:
#         'check.{flow_cell}.txt'
#     shell:'''
    
#     '''


# checkpoint tester_tock:
#     input:
#         inner = 'out.{flow_cell}.txt',
#         direct='output/blast/groupedReads/{flow_cell}/{file_num}'
#     output:
#         directory('testdir/{flow_cell}/{file_num}')
#     shell:'cp files'


# def agg_tester(wildcards):
#     checkpoint_output = checkpoints.tester_tock.get(**wildcards).output[0]
#     return '10.txt'


# rule aggregate_and_move_groups:
#     input:
#         agg_tester
#     output:
#         'agg.test.{flow_cell}.txt'
#     shell:'''
#     # copies files that have 3 or more lines to new location
#    # bash -c 'for file in "$@"; do if [ $(wc -l < "$file") -gt 3 ]; \
#     #then mv "$file" {output}; fi; done' bash {input}
#     '''

# def identify_groups(wildcards):
#     checkpoint_output = checkpoints.aggregate_and_move_groups.get(**wildcards).output[0]
#     glob_substrate = str(Path(checkpoint_output).join('{file_num}.{ref_name}.group.fastq'))
#     glob = glob_wildcards(glob_substrate, '{file_num}.{ref_name}.group.fastq')
#     ex = expand(
#         '/output/blast/prunedGroupedReads/{flow_cell}/{file_num}.{ref_name}.group.fastq',
#         flow_cell=wildcards.flow_cell, file_num=glob.file_num, 
#         ref_name=glob.ref_name
#     )
#     print(ex)
#     return ex
