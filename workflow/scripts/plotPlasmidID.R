library(ggplot2)
library(ggpubr)
library(RColorBrewer)

plasmid.id <- snakemake@input[['plasmid_id']]
expected.read.lengths <- snakemake@input[['read_lengths']]
plasmid.no.id.len <- snakemake@input[['no_plasmid_lengths']]


output.plot.path.count <- snakemake@output[['countPlot']]
output.plot.path.length <- snakemake@output[['lengthPlot']]
output.plot.path.noID <- snakemake@output[['lengthNoID']]


# Read the concat dataframe
id.df <- as.data.frame(
    read.table(file = plasmid.id, sep = '\t', header = FALSE)
)
exp.len <- as.data.frame(
    read.table(file = expected.read.lengths, sep = ',', header = TRUE)
)
no.id.len <- as.data.frame(
    read.table(file = expected.read.lengths, sep = ',', header = TRUE)
)

colnames(id.df) <- c('ref_name', 'read', 'length')

plasmid.species.count <- ggplot(id.df, aes(x=ref_name, fill=ref_name)) + geom_bar(color='black') +
                         scale_fill_viridis_d() + theme_pubr() +
                         theme(legend.position='none') + 
                         labs(y='Read count', x='', title='Read count by plasmid') +
                         coord_flip()


length.by.id <- ggplot(id.df, aes(x=length, fill=ref_name)) +
                geom_density(color='black') +
                theme_pubr() + theme(legend.position='none') +
                scale_fill_viridis_d() + labs(x='Plasmid', y='Read length') +
                geom_vline(data=exp.len, aes(xintercept=length, color=contruct_type),
                           linetype="dashed") +
                facet_wrap(~ref_name)


length.no.id <- ggplot(no.id.len, aes(x=length), fill='tan') +
                geom_density(color='black') +
                theme_pubr() + theme(legend.position='none') +
                scale_fill_viridis_d() + labs(x='Plasmid', y='Read length') +
                geom_vline(data=exp.len, aes(xintercept=length, color=contruct_type),
                           linetype="dashed")




ggsave(output.plot.path.count , plasmid.species.count, dpi=300, width=10, height=10, unit='in')
ggsave(output.plot.path.length, length.by.id, dpi=1000, width=30, height=30, unit='in')
ggsave(output.plot.path.noID , length.no.id, dpi=300, width=10, height=5, unit='in')
                         

