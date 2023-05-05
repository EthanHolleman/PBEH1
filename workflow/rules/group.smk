# from pathlib import Path
# import glob


# # checkpoint break_up_fastq:
# #     conda:
# #         "../envs/py.yml"
# #     input:
# #         flow_cell_reads="output/unzippedReads/{flow_cell}.fastq",
# #     output:
# #         directory("output/reads/seperated/{flow_cell}"),
# #     params:
# #         lines_per_file=40000,  # 4 lines per fastq entry
# #         flow_cell=lambda wildcards: wildcards.flow_cell,
# #         filename="sep.fastq",
# #     shell:
# #         """
# #     mkdir -p {output}/temp
# #     split -d -l {params.lines_per_file} \
# #         {input} {output}/temp/
# #     python scripts/processSplit.py {output}/temp {params.filename}
# #     """

# checkpoint break_up_fastq:
#     conda:
#         "../envs/py.yml"
#     input:
#         flow_cell_reads="output/unzippedReads/{flow_cell}.fastq",
#     output:
#         expand("output/reads/seperated/{flow_cell}/{file_num}.fastq",
#                 range(1, NUM_DIVS+1))
#     params:
#         lines_per_file=40000,  # 4 lines per fastq entry
#         flow_cell=lambda wildcards: wildcards.flow_cell,
#         filename="sep.fastq",
#     script:'../splitFastq.py'


# rule prune_reads_by_length:
#     conda:
#         "../envs/py.yml"
#     input:
#         "output/reads/seperated/{flow_cell}/{file_num}.fastq",
#     output:
#         "output/reads/seperatedTrim/{flow_cell}/{file_num}.prune.fastq",
#     params:
#         lower=1100,
#         upper=2700,
#     script:
#         "../scripts/lengthFilter.py"


# rule count_reads_and_length:
#     conda:
#         "../envs/py.yml"
#     input:
#         "output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq",
#     output:
#         read_counts="output/reads/readCounts/{flow_cell}/sep/{file_num}.count.tsv",
#         read_lengths="output/reads/readLengths/{flow_cell}/sep/{file_num}.lengths.tsv",
#     params:
#         flow_cell=lambda wildcards: wildcards.flow_cell,
#     script:
#         "../scripts/countUniqueReads.py"


# rule count_reads_no_plasmid_id:
#     conda:
#         "../envs/py.yml"
#     input:
#         "output/plasmidID/noID/{flow_cell}/{file_num}.no.ID.fastq",
#     output:
#         read_counts="output/plasmidID/noID/{flow_cell}/readCounts/{file_num}.no.ID.count.tsv",
#         read_lengths="output/plasmidID/noID/{flow_cell}/readLengths/{file_num}.no.ID.lengths.tsv",
#     params:
#         flow_cell=lambda wildcards: wildcards.flow_cell,
#     script:
#         "../scripts/countUniqueReads.py"


# def agg_count_reads(wildcards):
#     checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
#     print(wildcards, "wildcards")
#     file_num = glob_wildcards(
#         Path(checkpoint_output).joinpath("{file_num}/sep.fastq")
#     ).file_num
#     return expand(
#         "output/reads/readCounts/{flow_cell}/sep/{file_num}.count.tsv",
#         flow_cell=wildcards.flow_cell,
#         file_num=file_num,
#     )


# def agg_count_reads_no_id(wildcards):
#     checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
#     print(wildcards, "wildcards")
#     file_num = glob_wildcards(
#         Path(checkpoint_output).joinpath("{file_num}/sep.fastq")
#     ).file_num
#     return expand(
#         "output/plasmidID/noID/{flow_cell}/readLengths/{file_num}.no.ID.lengths.tsv",
#         flow_cell=wildcards.flow_cell,
#         file_num=file_num,
#     )


# rule concat_no_plasmid_id_lengths:
#     input:
#         agg_count_reads_no_id,
#     output:
#         "output/plasmidID/noID/{flow_cell}/allReadLengths/all.no.ID.lengths.tsv",
#     shell:
#         """
#     cat {input} > {output}
#     """


# rule plot_no_plasmid_id_lengths:
#     input:
#         read_lengths="output/plasmidID/noID/{flow_cell}/allReadLengths/all.no.ID.lengths.tsv",
#     output:
#         "",
#     script:
#         "../scripts/plotNoIDLengths.R"


# rule concat_read_counts:
#     input:
#         agg_count_reads,
#     output:
#         "output/reads/readCounts/{flow_cell}/{flow_cell}.all.count.tsv",
#     shell:
#         """
#     cat {input} > {output}
#     """


# rule concat_flow_cell_counts:
#     input:
#         expand(
#             "output/reads/readCounts/{flow_cell}/{flow_cell}.all.count.tsv",
#             flow_cell=FLOW_CELLS,
#         ),
#     output:
#         "output/reads/readCounts/allFlowCells.count.tsv",
#     shell:
#         """
#     cat {input} > {output}
#     """


# rule plot_read_counts:
#     conda:
#         "../envs/R.yml"
#     input:
#         "output/reads/readCounts/allFlowCells.count.tsv",
#     output:
#         "output/plots/readCounts/totalReadCounts.png",
#     script:
#         "../scripts/plotReadCounts.R"


# rule fastq_to_fasta:
#     conda:
#         "../envs/blast.yml"
#     input:
#         "output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq",
#     output:
#         "output/reads/seperatedFasta/{flow_cell}/{file_num}/sep.fasta",
#     shell:
#         """
#     seqtk seq -a {input} > {output}
#     """


# rule make_bowtie2_db:
#     conda:
#         "../envs/bowtie.yml"
#     input:
#         seqs="output/reads/seperatedFasta/{flow_cell}/{file_num}/sep.fasta",
#     output:
#         directory("output/bowtie/dbs/{flow_cell}/{flow_cell}_{file_num}"),
#     params:
#         db_name=lambda wildcards: f"{wildcards.flow_cell}_{wildcards.file_num}",
#         db_path=lambda wildcards: f"output/bowtie/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}",
#     shell:
#         """
#     mkdir -p {output}
#     bowtie2-build {input.seqs} {params.db_path}
#     """


# rule search_bowtie2_db:
#     conda:
#         "../envs/bowtie.yml"
#     input:
#         db="output/bowtie/dbs/{flow_cell}/{flow_cell}_{file_num}",
#         plasmid_id=config["PLASMID_ID_SEQS"],
#     output:
#         "output/plasmidID/bowtieSearch/sam/{flow_cell}/{flow_cell}.{file_num}.sam",
#     params:
#         db_path=lambda wildcards: f"output/bowtie/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}",
#     shell:
#         """
#     bowtie2 -x {params.db_path} -f {input.plasmid_id} -S {output} --local -a
#     """


# rule sam_to_tsv:
#     input:
#         "output/plasmidID/bowtieSearch/sam/{flow_cell}/{flow_cell}.{file_num}.sam",
#     output:
#         "output/plasmidID/bowtieSearch/tsv/{flow_cell}/{flow_cell}.{file_num}.tsv",
#     shell:
#         """
#     awk '{{if($1 ~ /^@/) next; print $1 \"\\t\" $3 \"\\t\" $4 \"\\t\" $10}}' {input} > {output}
#     """


# rule make_blast_db:
#     # build blastn databases for the split reads
#     conda:
#         "../envs/blast.yml"
#     input:
#         seqs="output/reads/seperatedFasta/{flow_cell}/{file_num}/sep.fasta",
#     output:
#         directory("output/blast/dbs/{flow_cell}/{flow_cell}_{file_num}"),
#     params:
#         db_name=lambda wildcards: f"{wildcards.flow_cell}_{wildcards.file_num}",
#         db_path=lambda wildcards: f"output/blast/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}.blastn.db",
#     shell:
#         """
#     mkdir -p {output}
#     makeblastdb -in {input.seqs} -parse_seqids \
#     -title "{params.db_name}" -dbtype nucl -out {params.db_path}
#     """


# rule local_blast_plasmids:
#     # blast the plasmid identifying sequences against the database
#     # of reads
#     conda:
#         "../envs/blast.yml"
#     input:
#         db="output/blast/dbs/{flow_cell}/{flow_cell}_{file_num}",
#         query=config["PLASMID_ID_SEQS"],
#     output:
#         "output/plasmidID/search/{flow_cell}/{file_num}/blastn.search.tsv",
#     params:
#         db_path=lambda wildcards: f"output/blast/dbs/{wildcards.flow_cell}/{wildcards.flow_cell}_{wildcards.file_num}/{wildcards.flow_cell}_{wildcards.file_num}.blastn.db",
#     shell:
#         """
#     cat {input.query} | blastn -db {params.db_path} -outfmt 6 -perc_identity 85 -task blastn > {output}
#     """


# rule identify_plasmid:
#     conda:
#         "../envs/py.yml"
#     input:
#         bowtie_tsv="output/plasmidID/bowtieSearch/tsv/{flow_cell}/{flow_cell}.{file_num}.tsv",
#     output:
#         "output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv",
#     params:
#         ref_dir=config["GENOME_DIR"],
#         id_df=config["PLASMID_ID_TABLE"],
#     script:
#         "../scripts/assignPlasmid.py"


# rule merge_read_length_plasmid_id:
#     conda:
#         "../envs/py.yml"
#     input:
#         plasmid_id="output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv",
#         lengths="output/reads/readLengths/{flow_cell}/sep/{file_num}.lengths.tsv",
#     output:
#         "output/plasmidID/{flow_cell}/mergedIDs/{flow_cell}.{file_num}.lengths.ID.tsv",
#     script:
#         "../scripts/addLengthToIDs.py"


# def agg_read_lengths(wildcards):
#     checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
#     print(wildcards, "wildcards")
#     file_num = glob_wildcards(
#         Path(checkpoint_output).joinpath("{file_num}/sep.fastq")
#     ).file_num
#     return expand(
#         "output/plasmidID/{flow_cell}/mergedIDs/{flow_cell}.{file_num}.lengths.ID.tsv",
#         flow_cell=wildcards.flow_cell,
#         file_num=file_num,
#     )


# rule concat_read_lengths:
#     input:
#         agg_read_lengths,
#     output:
#         "output/reads/readLengths/{flow_cell}/all.lengths.tsv",
#     shell:
#         """
#     cat {input} > {output}
#     """


# def agg_plasmid_ids(wildcards):
#     checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
#     print(wildcards, "wildcards")
#     file_num = glob_wildcards(
#         Path(checkpoint_output).joinpath("{file_num}/sep.fastq")
#     ).file_num
#     return expand(
#         "output/plasmidID/{flow_cell}/mergedIDs/{flow_cell}.{file_num}.lengths.ID.tsv",
#         flow_cell=wildcards.flow_cell,
#         file_num=file_num,
#     )


# def agg_plasmid_no_ids(wildcards):
#     checkpoint_output = checkpoints.break_up_fastq.get(**wildcards).output[0]
#     print(wildcards, "wildcards")
#     file_num = glob_wildcards(
#         Path(checkpoint_output).joinpath("{file_num}/sep.fastq")
#     ).file_num
#     return expand(
#         "output/plasmidID/noID/{flow_cell}/{file_num}.no.ID.fastq",
#         flow_cell=wildcards.flow_cell,
#         file_num=file_num,
#     )


# rule group_no_id_read:
#     conda:
#         "../envs/py.yml"
#     input:
#         plasmid_id="output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv",
#         reads="output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq",
#     output:
#         "output/plasmidID/noID/{flow_cell}/{file_num}.no.ID.fastq",
#     script:
#         "../scripts/findUnIDPlasmids.py"


# rule concat_plasmid_ids:
#     input:
#         agg_plasmid_ids,
#     output:
#         "output/plasmidID/allPlasmidID/{flow_cell}/{flow_cell}.concat.plasmidID.tsv",
#     shell:
#         """
#     cat {input} > {output}
#     """


# rule agg_plasmid_no_ids:
#     input:
#         agg_plasmid_no_ids,
#     output:
#         "output/plasmidID/allPlasmidNoIDs/{flow_cell}/{flow_cell}.concat.no.ID.tsv",
#     shell:
#         """
#     cat {input} > {output}
#     """


# rule plot_plasmid_id_count:
#     conda:
#         "../envs/R.yml"
#     input:
#         plasmid_id="output/plasmidID/allPlasmidID/{flow_cell}/{flow_cell}.concat.plasmidID.tsv",
#         no_plasmid_id="output/plasmidID/allPlasmidNoIDs/{flow_cell}/{flow_cell}.concat.no.ID.tsv",
#         read_lengths="../resources/productLengths/productLengthsExpected.tsv",
#     output:
#         countPlot="output/plots/plasmidID/{flow_cell}.plasmid.count.png",
#         lengthPlot="output/plots/plasmidID/{flow_cell}.plasmid.length.count.png",
#         lengthNoID="output/plots/plasmidID/{flow_cell}.plasmid.length.no.ID.png",
#         pFC9Plot="output/plots/plasmidID/{flow_cell}.pFC9.length.png"
#     script:
#         "../scripts/plotPlasmidID.R"


# checkpoint group_plasmids:
#     conda:
#         "../envs/py.yml"
#     input:
#         plasmid_id_df="output/plasmidID/{flow_cell}/{file_num}/plasmidID.tsv",
#         reads="output/reads/seperatedTrim/{flow_cell}/{file_num}/sep.fastq",
#         id_count="output/plots/plasmidID/{flow_cell}.plasmid.count.png",
#         read_count="output/plots/readCounts/totalReadCounts.png",
#     output:
#         directory("output/groupedReads/{flow_cell}/{file_num}/{plasmid}.fastq"),
#     script:
#         "../scripts/groupReadsByPlasmid.py"
    

# def agg_plasmids(wildcards):
#     checkpoint_output = checkpoints.group_plasmids.get(**wildcards).output[0]
#     print(wildcards, "wildcards")
#     file_num = glob_wildcards(
#         Path(checkpoint_output).joinpath("{file_num}/sep.fastq")
#     ).file_num
#     return expand(
#         "output/plasmidID/noID/{flow_cell}/{file_num}.no.ID.fastq",
#         flow_cell=wildcards.flow_cell,
#         file_num=file_num,
#     )

# rule agg_plas_names:
#     input:
#         agg_plasmids
#     output:
#         'output/testids.{flow_cell}.{plasmid}'
#     shell:'touch'
