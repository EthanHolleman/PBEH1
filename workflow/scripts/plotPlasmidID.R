library(ggplot2)
library(ggpubr)
library(RColorBrewer)

plasmid.id <- snakemake@input[['plasmid_id']]
expected.read.lengths <- snakemake@input[['read_lengths']]
plasmid.no.id.len <- snakemake@input[['no_plasmid_lengths']]


output.plot.path.count <- snakemake@output[['countPlot']]
output.plot.path.length <- snakemake@output[['lengthPlot']]
output.plot.path.noID <- snakemake@output[['lengthNoID']]
output.plot.path.pFC9 <- snakemake@output[['pFC9Plot']]


# Read the concat dataframe
id.df <- as.data.frame(
    read.table(file = plasmid.id, sep = '\t', header = FALSE)
)
exp.len <- as.data.frame(
    read.table(file = expected.read.lengths, sep = ',', header = TRUE)
)
# no.id.len <- as.data.frame(
#     read.table(file = expected.read.lengths, sep = ',', header = TRUE)
# )

colnames(id.df) <- c('ref_name', 'read', 'length')

remove_low_count_plasmids <- function(df) {
  # Calculate the total count for each plasmid
  plasmid_counts <- aggregate(df$length, by = list(df$ref_name), FUN = length)
  print(plasmid_counts)
  names(plasmid_counts) <- c("name", "count")
  
  # Calculate the mean and standard deviation of the counts
  mean_count <- mean(plasmid_counts$count)
  sd_count <- sd(plasmid_counts$count)
  
  # Calculate the threshold below which plasmids will be removed
  #threshold <- mean_count - 2 * sd_count
  threshold <- 3000
  
  # Subset the data frame to keep only plasmids with counts above the threshold
  df_filtered <- subset(df, ref_name %in% plasmid_counts$name[plasmid_counts$count > threshold])
  
  # Return the filtered data frame
  return(df_filtered)
}


plasmid.species.count <- ggplot(id.df, aes(x=ref_name, fill=ref_name)) + geom_bar(color='black') +
                         scale_fill_viridis_d() + theme_pubr() +
                         theme(legend.position='none') + 
                         labs(y='Read count', x='', title='Read count by plasmid') +
                         coord_flip()


id.df.len <- remove_low_count_plasmids(id.df)
print('======================================')
print(nrow(id.df))
print(nrow(id.df.len))

length.by.id <- ggplot(id.df.len, aes(x=length, y=ref_name, fill=ref_name)) +
                geom_violin(color='black', alpha=0.2) +
                theme_pubr() + theme(legend.position='none') +
                scale_fill_viridis_d() + labs(x='Plasmid', y='Read length') +
                geom_vline(data=exp.len, aes(xintercept=length, color=contruct_type),
                           linetype="dashed")


# length.no.id <- ggplot(no.id.len, aes(x=length), fill='tan') +
#                 geom_histogram(color='black') +
#                 theme_pubr() +
#                 scale_fill_viridis_d() + labs(x='Plasmid', y='Read length') +
#                 geom_vline(data=exp.len, aes(xintercept=length, color=contruct_type),
#                            linetype="dashed")

pFC9 <- subset(id.df, ref_name=='pFC9.gb')

pFC9.plt <- ggplot(pFC9, aes(x=length, y=ref_name, fill=ref_name)) +
                geom_violin(color='black', alpha=0.2) +
                theme_pubr() + theme(legend.position='none') +
                scale_fill_viridis_d() + labs(x='Plasmid', y='Read length') +
                geom_vline(data=exp.len, aes(xintercept=length, color=contruct_type),
                           linetype="dashed")



ggsave(output.plot.path.count , plasmid.species.count, dpi=300, width=20, height=10, unit='in')
ggsave(output.plot.path.length, length.by.id, dpi=1000, width=8, height=20, unit='in')
# ggsave(output.plot.path.noID , length.no.id, dpi=300, width=10, height=5, unit='in')
ggsave(output.plot.path.pFC9, pFC9.plt, dpi=300, width=8, height=5)
                         

