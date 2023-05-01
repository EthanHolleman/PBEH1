# Plot unique and total read counts from each flow cell

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

read.counts.path <- snakemake@input[[1]]
output.plot.path <- snakemake@output[[1]]

# Read the concat dataframe
read.counts.df <- as.data.frame(
    read.table(file = read.counts.path, sep = '\t', header = FALSE)
)
colnames(read.counts.df) <- c('filepath', 'uniqueReadCount', 'totalReads', 'flow_cell')


# make barplots

unique.reads.plt <- ggplot(read.counts.df, aes(x=flow_cell, y=uniqueReadCount, fill=flow_cell)) +
                   geom_bar(stat='identity') + labs(title='Unique reads by flow cell', x='', y='Read count') +
                   scale_fill_brewer(palette='Dark2') + theme_pubr() +  theme(legend.position = "none")


total.reads.plt <- ggplot(read.counts.df, aes(x=flow_cell, y=totalReads, fill=flow_cell)) +
                   geom_bar(stat='identity') + labs(title='Total reads by flow cell', x='', y='Read count') +
                   scale_fill_brewer(palette='Dark2') + theme_pubr() + theme(legend.position = "none")


merge.plot <- ggarrange(total.reads.plt, unique.reads.plt, nrow=1, ncol=2)

# write plot to output path
ggsave(output.plot.path, merge.plot, dpi=300, width=10, height=10, units='in')

                
